/**
 * bs2_cooccurrence.js
 * --------------------------------------------------------------------------
 * Decision-support logic for applying ACMG/AMP BS2 ("observed in a healthy
 * individual who should manifest disease") using gnomAD genotype counts and
 * gnomAD variant co-occurrence (inferred phasing).
 *
 * Two pathways are handled:
 *   A. observed_biallelic     - a directly observed healthy homozygote (AR)
 *                               or hemizygote (XLR). No phasing required.
 *                               This is the canonical, strongest BS2 case.
 *   B. inferred_compound_het  - the query VUS is predicted IN TRANS with an
 *                               established P/LP partner in a gnomAD individual,
 *                               via the co-occurrence (EM/p_chet) tool. Phase
 *                               is inferred, not observed, so it is treated
 *                               more conservatively than pathway A.
 *
 * PROVENANCE
 *   - ACMG/AMP 2015 (Richards et al.) BS2.
 *   - ClinGen SVI / VCEP practice: gnomAD incidence usable for severe
 *     paediatric AR; reduced-penetrance genes (e.g. ATM) excluded from BS2.
 *   - gnomAD variant co-occurrence: pairs must be in gnomAD EXOMES, same gene,
 *     global AF <= 5%, and coding / near-splice / UTR. Phase is an EM estimate
 *     (p_trans). Singletons are lower-accuracy.
 *
 * SAFETY
 *   This module is DECISION SUPPORT, not an autonomous classifier. Every
 *   non-null verdict carries requiresExpertReview = true. All numeric
 *   thresholds are STARTING DEFAULTS and must be calibrated per gene/disease
 *   (ideally to the relevant ClinGen VCEP specification) before clinical use.
 *
 * ----------------------------------------------------------------------------
 * PORT NOTE (Lab Variant Search): the upstream source is an ES module. This
 * codebase loads every src/*.js as a plain global <script> (no import/export),
 * so the module is wrapped in an IIFE that exposes a single `window.BS2`
 * namespace. The clinical logic and thresholds below are byte-for-byte
 * identical to the supplied bs2_cooccurrence_logic.js — only the ESM `export`
 * keywords were removed and the public API attached to window.BS2 at the end.
 * --------------------------------------------------------------------------
 */
(function () {
  'use strict';

  // ---------------------------------------------------------------------------
  // Strength -> Bayesian (Tavtigian) points. Classic ACMG benign strengths are
  // Supporting (-1), Strong (-4), Stand-alone (-8). "Moderate" (-2) is a
  // points-system extension some VCEPs use; BS2 never reaches Stand-alone here.
  // ---------------------------------------------------------------------------
  const BENIGN_POINTS = Object.freeze({
    'Supporting': -1,
    'Moderate':   -2,
    'Strong':     -4,
  });

  const DEFAULT_OPTIONS = Object.freeze({
    maxAF: 0.05,                 // gnomAD co-occurrence eligibility ceiling
    pTransMin: 0.5,              // min P(trans) for pathway B (prefer gnomAD's
                                 // own categorical "predicted in trans" if given)
    minIndividualsForStrong: 2,  // observed healthy homozygotes for Strong (A)
    allowUnknownPenetrance: false,
    allowSingleton: true,        // permitted but downgraded one notch
    allowHypomorphicOverride: false, // genotype-phenotype-complex genes (GALC/GBA/CFTR...)
  });

  /**
   * @typedef {Object} GeneDiseaseConfig
   * @property {string}  gene
   * @property {string}  disease
   * @property {'AR'|'XLR'|'AD'} inheritance
   * @property {'full'|'high'|'reduced'|'variable'|'unknown'} penetrance
   * @property {'pediatric_severe'|'pediatric'|'adult'|'variable'} onset
   * @property {boolean} hypomorphicComplexity  // known genotype-phenotype / hypomorph traps
   * @property {Object}  [thresholds]           // optional per-gene overrides of DEFAULT_OPTIONS
   */

  /**
   * Gate 1: is this gene/disease even eligible for a BS2 "should be affected"
   * argument from gnomAD (whose individuals are presumed free of severe
   * paediatric disease)? Returns a hard pass/fail with human-readable reasons.
   */
  function checkGeneDiseaseEligibility(geneConfig, options = {}) {
    const opt = { ...DEFAULT_OPTIONS, ...options, ...(geneConfig.thresholds || {}) };
    const reasons = [];

    const penOk = geneConfig.penetrance === 'full' || geneConfig.penetrance === 'high' ||
      (opt.allowUnknownPenetrance && geneConfig.penetrance === 'unknown');
    if (!penOk) {
      reasons.push(`Penetrance "${geneConfig.penetrance}" is not full/high; ` +
        `BS2 assumes a carrier of the genotype would manifest disease.`);
    }

    const onsetOk = geneConfig.onset === 'pediatric_severe' || geneConfig.onset === 'pediatric';
    if (!onsetOk) {
      reasons.push(`Onset "${geneConfig.onset}" is not paediatric; adult/late-onset ` +
        `disease may be undiagnosed in gnomAD adults, breaking the "healthy" assumption.`);
    }

    if (geneConfig.hypomorphicComplexity && !opt.allowHypomorphicOverride) {
      reasons.push(`Gene flagged for genotype-phenotype/hypomorph complexity; ` +
        `tolerated in-trans combinations can occur (e.g. mild-allele combinations). ` +
        `Manual override required.`);
    }

    return { eligible: reasons.length === 0, reasons };
  }

  /**
   * Gate 2: are the two variants eligible for gnomAD co-occurrence phasing?
   * (Pathway B only.) Both must be in gnomAD exomes, AF <= ceiling, and in a
   * coding / near-splice / UTR region.
   */
  function checkCooccurrenceEligibility(query, partner, options = {}) {
    const opt = { ...DEFAULT_OPTIONS, ...options };
    const reasons = [];
    const okRegions = new Set(['coding', 'splice', 'utr']);

    for (const [label, v] of [['query', query], ['partner', partner]]) {
      if (!v) { reasons.push(`${label} variant missing.`); continue; }
      if (v.inGnomadExomes === false) reasons.push(`${label} not present in gnomAD exomes.`);
      if (typeof v.afGlobal === 'number' && v.afGlobal > opt.maxAF) {
        reasons.push(`${label} global AF ${v.afGlobal} exceeds co-occurrence ceiling ${opt.maxAF}.`);
      }
      if (v.region && !okRegions.has(v.region)) {
        reasons.push(`${label} region "${v.region}" outside coding/splice/UTR scope.`);
      }
    }

    return { eligible: reasons.length === 0, reasons };
  }

  /**
   * @typedef {Object} VariantInput
   * @property {string}  id
   * @property {number}  [afGlobal]
   * @property {'coding'|'splice'|'utr'|'other'} [region]
   * @property {boolean} [inGnomadExomes]
   * @property {'P'|'LP'|'VUS'|'LB'|'B'} [classification]   // partner only
   * @property {string}  [reviewStatus]                     // partner only, free text
   */

  /**
   * @typedef {Object} ObservationInput
   * @property {'homozygous'|'hemizygous'|'compound_het_inferred'} mode
   * @property {number}  individualCount   // # gnomAD individuals with this genotype/pair
   * @property {number}  [pTrans]          // required for compound_het_inferred
   * @property {boolean} [isSingleton]
   * @property {boolean} [qcPass]          // passes gnomAD QC / not in segdup-LCR artifact
   * @property {boolean} [ancestryConcordant] // pair seen within one genetic-ancestry group
   */

  /**
   * Main entry point. Evaluate whether a benign observation code applies to the
   * QUERY variant given the gene/disease context, the (optional) pathogenic
   * partner, and the gnomAD observation.
   *
   * @param {{query: VariantInput, partner?: VariantInput|null,
   *          observation: ObservationInput, geneConfig: GeneDiseaseConfig}} input
   * @param {Object} [options] overrides of DEFAULT_OPTIONS
   * @returns {{
   *   applicable: boolean, code: 'BS2'|null,
   *   strength: 'Strong'|'Moderate'|'Supporting'|null,
   *   bayesianPoints: number, pathway: string|null,
   *   confidence: 'high'|'moderate'|'low'|null,
   *   flags: string[], rationale: string[], requiresExpertReview: boolean
   * }}
   */
  function evaluateBenignObservation(input, options = {}) {
    const { query, partner = null, observation, geneConfig } = input;
    const opt = { ...DEFAULT_OPTIONS, ...options, ...(geneConfig.thresholds || {}) };

    const out = {
      applicable: false, code: null, strength: null, bayesianPoints: 0,
      pathway: null, confidence: null, flags: [], rationale: [],
      requiresExpertReview: true,
    };
    const reject = (msg) => { out.rationale.push(msg); return out; };

    // --- Gate 1: gene/disease eligibility -----------------------------------
    const g1 = checkGeneDiseaseEligibility(geneConfig, opt);
    if (!g1.eligible) { g1.reasons.forEach(r => out.rationale.push(r)); return out; }

    // --- Inheritance vs observation-mode consistency ------------------------
    const mode = observation.mode;
    if (mode === 'homozygous'  && geneConfig.inheritance !== 'AR')
      return reject('Homozygous-observation BS2 expects autosomal recessive inheritance.');
    if (mode === 'hemizygous'  && geneConfig.inheritance !== 'XLR')
      return reject('Hemizygous-observation BS2 expects X-linked recessive inheritance.');
    if (mode === 'compound_het_inferred' && geneConfig.inheritance !== 'AR')
      return reject('Inferred compound-het BS2 expects autosomal recessive inheritance.');

    // --- QC ------------------------------------------------------------------
    if (observation.qcPass === false)
      return reject('Observation failed QC (artifact / low-complexity / segdup).');
    if ((observation.individualCount || 0) < 1)
      return reject('No qualifying healthy individual observed.');

    // ========================================================================
    // PATHWAY A: directly observed biallelic genotype (homozygous / hemizygous)
    // ========================================================================
    if (mode === 'homozygous' || mode === 'hemizygous') {
      out.pathway = 'observed_biallelic';
      out.code = 'BS2';

      let strength = (observation.individualCount >= opt.minIndividualsForStrong)
        ? 'Strong' : 'Moderate';
      out.rationale.push(
        `${observation.individualCount} healthy ${mode} individual(s) observed in gnomAD for a ` +
        `${geneConfig.penetrance}-penetrance, ${geneConfig.onset} ${geneConfig.inheritance} ` +
        `disease — such individuals would be expected to manifest disease.`);

      if (observation.isSingleton && opt.allowSingleton) {
        strength = downgrade(strength);
        out.flags.push('Singleton observation: strength downgraded one notch.');
      } else if (observation.isSingleton && !opt.allowSingleton) {
        return reject('Singleton observations disallowed by configuration.');
      }

      out.confidence = strength === 'Strong' ? 'high' : 'moderate';
      return finalize(out, strength);
    }

    // ========================================================================
    // PATHWAY B: inferred compound heterozygote (VUS in trans with P/LP partner)
    // ========================================================================
    if (mode === 'compound_het_inferred') {
      out.pathway = 'inferred_compound_het';

      // Partner must be an established pathogenic anchor.
      if (!partner) return reject('Compound-het pathway requires a pathogenic partner variant.');
      if (!['P', 'LP'].includes(partner.classification))
        return reject(`Partner classification "${partner.classification}" is not P/LP; ` +
          `cannot anchor a compound-het BS2 argument.`);

      // Co-occurrence eligibility.
      const g2 = checkCooccurrenceEligibility(query, partner, opt);
      if (!g2.eligible) { g2.reasons.forEach(r => out.rationale.push(r)); return out; }

      // Phase: prefer gnomAD's categorical call; otherwise threshold p_trans.
      const pTrans = observation.pTrans;
      if (typeof pTrans !== 'number')
        return reject('pTrans not provided; phase cannot be established for this pair.');
      if (pTrans < opt.pTransMin)
        return reject(`Predicted phase not in trans (pTrans=${pTrans} < ${opt.pTransMin}); ` +
          `a cis/ambiguous pair is a carrier and is uninformative for BS2.`);

      out.code = 'BS2';
      // Inferred phase => one notch more conservative than observed biallelic.
      // Single inferred-trans individual -> Supporting; >=2 independent -> Moderate.
      let strength = (observation.individualCount >= 2) ? 'Moderate' : 'Supporting';
      out.rationale.push(
        `Query VUS predicted IN TRANS (pTrans=${pTrans}) with established ` +
        `${partner.classification} variant ${partner.id} in ${observation.individualCount} ` +
        `healthy gnomAD individual(s); the resulting biallelic genotype would be expected ` +
        `to cause this ${geneConfig.onset} ${geneConfig.inheritance} disease.`);
      out.flags.push('Phase is INFERRED (EM/p_chet), not observed — interpret conservatively.');

      if (partner.classification === 'LP') {
        out.flags.push('Partner is LP (not P): anchor confidence reduced.');
        strength = downgrade(strength);
      }
      if (observation.isSingleton) {
        out.flags.push('Singleton pair: gnomAD phasing accuracy is lower for singletons.');
        strength = downgrade(strength);
      }
      if (observation.ancestryConcordant === false) {
        out.flags.push('Pair not confined to one genetic-ancestry group: phasing estimate weaker.');
      }

      out.confidence = strength === 'Moderate' ? 'moderate' : 'low';
      return finalize(out, strength);
    }

    return reject(`Unsupported observation mode "${mode}".`);
  }

  // ---------------------------------------------------------------------------
  // helpers
  // ---------------------------------------------------------------------------
  function downgrade(strength) {
    if (strength === 'Strong') return 'Moderate';
    if (strength === 'Moderate') return 'Supporting';
    return null; // below Supporting -> not usable
  }

  function finalize(out, strength) {
    if (!strength) {
      out.rationale.push('Cumulative downgrades fell below Supporting; code not applied.');
      return out; // applicable stays false
    }
    out.applicable = true;
    out.strength = strength;
    out.bayesianPoints = BENIGN_POINTS[strength];
    return out;
  }

  // ---------------------------------------------------------------------------
  // Public API — attached to the global namespace (this app is non-ESM).
  // ---------------------------------------------------------------------------
  window.BS2 = {
    BENIGN_POINTS,
    DEFAULT_OPTIONS,
    checkGeneDiseaseEligibility,
    checkCooccurrenceEligibility,
    evaluateBenignObservation,
  };
})();

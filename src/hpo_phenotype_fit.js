/*
  hpo_phenotype_fit.js
  ====================
  Ontology-aware HPO phenotype-fit scorer (Phenomizer / LIRICAL style).
  Pure module exported on window.HpoFit (no implicit data.* reads — mirrors
  window.BS2 in bs2_cooccurrence.js). Consumes three local tables built by
  scripts/build_hpo_background.py and fetched from data/hpo/:
      hpo_background.json         : { "HP:0001250": 0.18, ... }   term -> background p(h)
      hpo_parents.json            : { "HP:0001250": ["HP:0012638"], ... } is_a edges
      hpo_disease_phenotypes.json : { "OMIM:614558": { terms:[{id,f}], completeness } }

  Two consumers:
    - The rigorous LR scorer (createScorer().score / HpoFit.scoreDisease) is
      UNIDIRECTIONAL (true-path rule): a patient term satisfies its ancestors.
      This is the LIRICAL semantics the background table was built for.
    - HpoFit.linSim() is a helper for the SEPARATE omim_local.js F-score matcher,
      which is deliberately bidirectional + IC-floored for candidate surfacing.

  LR per disease phenotype (frequency f, background p):
    present  -> f / p        excluded -> (1-f)/(1-p)        f == null -> 1 (no evidence)
  Combined in log space, applied to the prior in log-odds:
    posterior = logistic( logit(prior) + Σ ln(LR_i) )
*/
window.HpoFit = (function () {
  const EPS = 1e-3;
  const FLOOR_P = 0.05;            // IC floor for the matcher: ignore shared ancestors with p > 5%
  const clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x));

  function createScorer({ background = {}, parents = {}, defaultFallbackP = 0.02 } = {}) {
    const ancMemo = new Map();

    function ancestors(term) {
      // inclusive: term + all is_a ancestors
      const cached = ancMemo.get(term);
      if (cached) return cached;
      const acc = new Set([term]);
      const stack = [...(parents[term] || [])];
      while (stack.length) {
        const p = stack.pop();
        if (!acc.has(p)) {
          acc.add(p);
          const up = parents[p];
          if (up) for (const q of up) stack.push(q);
        }
      }
      ancMemo.set(term, acc);
      return acc;
    }

    function inducedUp(terms) {
      const set = new Set();
      for (const t of terms || []) for (const a of ancestors(t)) set.add(a);
      return set;
    }

    function bgP(term, fallbackP = defaultFallbackP) {
      const p = background[term];
      return clamp(p == null ? fallbackP : p, EPS, 1 - EPS);
    }

    function ic(term) { return -Math.log(bgP(term)); }

    function score({
      diseasePhenotypes = [],
      present = [],
      excluded = [],
      prior = 1e-3,
      fallbackP = defaultFallbackP,
    } = {}) {
      const presentUp = inducedUp(present);   // patient terms + ancestors (true-path)
      const excludedSet = new Set(excluded);

      let sumLn = 0;
      let covNum = 0, covDen = 0;             // frequency-weighted coverage (length-normalized cross-check)
      const rows = diseasePhenotypes.map(({ id, f }) => {
        // f == null/undefined: frequency genuinely unknown -> contributes nothing (LR = 1).
        const known = f != null;
        const ff = known ? clamp(f, EPS, 1 - EPS) : null;
        const p = bgP(id, fallbackP);
        let status, lr;
        if (presentUp.has(id)) {
          status = "present";
          if (known) { lr = ff / p; covNum += ff; covDen += ff; }
          else { lr = 1; }                    // present but unknown f: no evidence, no coverage weight
        } else {
          let isExcluded = false;
          for (const a of ancestors(id)) { if (excludedSet.has(a)) { isExcluded = true; break; } }
          if (isExcluded) {
            status = "excluded";
            lr = known ? (1 - ff) / (1 - p) : 1;
          } else {
            status = "unobserved";
            lr = 1;
            if (known) covDen += ff;          // disease term we didn't observe still counts toward coverage denom
          }
        }
        const ln = Math.log(lr);
        if (status !== "unobserved" && lr !== 1) sumLn += ln;
        return { id, f: ff, p, status, lr, ln };
      });

      const pri = clamp(prior, 1e-12, 1 - 1e-12);
      const logPostOdds = Math.log(pri) - Math.log(1 - pri) + sumLn;
      const posterior = 1 / (1 + Math.exp(-logPostOdds));
      const coverage = covDen > 0 ? covNum / covDen : 0;

      return { posterior, totalLR: Math.exp(sumLn), sumLn, coverage, rows };
    }

    return { score, ancestors, inducedUp, bgP, ic };
  }

  // ── loaded-singleton layer (convenience cache, like omimLocalReady) ──────
  const HpoFit = { createScorer, ready: false };
  let _loadPromise = null;
  const _state = { ready: false, scorer: null, diseaseMap: null };

  HpoFit.load = function (baseUrl = "data/hpo") {
    if (_loadPromise) return _loadPromise;
    _loadPromise = (async () => {
      try {
        const getJson = (name) =>
          fetch(`${baseUrl}/${name}`).then((r) => { if (!r.ok) throw new Error(name); return r.json(); });
        const [background, parents, diseaseMap] = await Promise.all([
          getJson("hpo_background.json"),
          getJson("hpo_parents.json"),
          getJson("hpo_disease_phenotypes.json"),
        ]);
        _state.scorer = createScorer({ background, parents });
        _state.diseaseMap = diseaseMap;
        _state.ready = true;
        HpoFit.ready = true;
      } catch (e) {
        // Local tables absent (not built yet) — degrade silently; callers no-op.
        _state.ready = false;
        HpoFit.ready = false;
      }
      return _state;
    })();
    return _loadPromise;
  };

  /**
   * scoreDisease(diseaseKey, present, prior)
   * Inputs: OMIM/ORPHA key (e.g. "OMIM:312750"), patient HP IDs (explicit), prior.
   * Output: { posterior, totalLR, sumLn, coverage, completeness, matched, informative, rows } or null.
   * Failure: returns null if tables aren't loaded or the disease has no annotation.
   */
  HpoFit.scoreDisease = function (diseaseKey, present, prior = 1e-3) {
    if (!_state.ready || !_state.diseaseMap) return null;
    const entry = _state.diseaseMap[diseaseKey];
    if (!entry || !entry.terms?.length) return null;
    const res = _state.scorer.score({ diseasePhenotypes: entry.terms, present, prior });
    const matched = res.rows.filter((r) => r.status === "present").length;
    const informative = res.rows.filter((r) => r.status === "present" && r.f != null).length;
    return { ...res, matched, informative, completeness: entry.completeness ?? null };
  };

  /**
   * linSim(a, b) — IC-floored Lin similarity between two HP terms, for the
   * omim_local.js bidirectional matcher. 0 when the most-informative common
   * ancestor is too shallow (p > FLOOR_P) — prevents matching on "Phenotypic
   * abnormality" etc. Returns 1.0 for an exact term match below the floor.
   */
  HpoFit.linSim = function (a, b) {
    if (!_state.ready) return 0;
    const A = _state.scorer.ancestors(a);
    const B = _state.scorer.ancestors(b);
    const small = A.size <= B.size ? A : B;
    const big = small === A ? B : A;
    let micaP = 2; // sentinel > 1
    for (const t of small) {
      if (big.has(t)) { const p = _state.scorer.bgP(t); if (p < micaP) micaP = p; }
    }
    if (micaP > FLOOR_P) return 0;                 // IC floor
    const icM = -Math.log(micaP);
    const denom = _state.scorer.ic(a) + _state.scorer.ic(b);
    if (denom <= 0) return 0;
    return clamp((2 * icM) / denom, 0, 1);
  };

  return HpoFit;
})();

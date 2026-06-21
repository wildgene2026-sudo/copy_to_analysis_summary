// ── MAIN APPLICATION LOGIC & STATE ────────────────────────────────

const phenotypeCache = {};

let data = {
  coords: { gene: '', hgvs: '-', chrom: '', pos38: '', pos19: '', ref: '', alt: '', strand: 1, hg38String: '', hg19String: '', userCoord: null, hg19FromVV: false, primaryFromVV: false },
  ptPhenoTexts: [], ptPhenoGroups: [], selectedHpoTerms: {},
  scores: { revel: null, spliceAI: null, spliceAI_AG: null, spliceAI_AL: null, spliceAI_DG: null, spliceAI_DL: null, spliceAI_DP_AG: null, spliceAI_DP_AL: null, spliceAI_DP_DG: null, spliceAI_DP_DL: null, spliceAiManeMatch: true, alphaMissenseScore: null, alphaMissensePred: null },
  ensembl: { transcript: null, protein: null, vepConsequence: null, vepExon: null, vepTotalExons: null, vepAminoAcids: null, vepProteinStart: null, strand: 1, rsId: null, proteinId: null, allTranscripts: [], mutalyzerHgvs: null, _geneFromVep: null, _transcriptFromVep: null, _hgvscFromVep: null, _hgvspFromVep: null, _hg38FromVep: null },
  gnomad: { popmax: null, detailed: null, viewMode: 'total', geneConstraint: null, regionalConstraint: null, regionalConstraintMsg: null, regionalOE: null, regionalRange: null, regionalPValue: null },
  dosage: { haplo: null, triplo: null, source: null, haploId: null },
  geneValidity: { curations: [], gene: null, fetched: false, focusMondo: null },
  geneReviews: { gene: null, chapters: [], fetched: false, focusNbk: null },
  curationSupport: { mechanism: null, phenoMatches: [], chapterTitle: null, chapterUrl: null },
  clinvar: { sig: null, stars: 0, reviewStatus: null, subs: 0, uid: null, accession: null, interpretationId: undefined, evidenceFlags: null, ps1Suggested: null, pm5Suggested: null, evidenceRepoUuid: null, historicalProteinNotations: [], detail: null },
  litvar: { combinedTotal: 0, uniqueNotations: [], _lastKey: null },
  startLoss: { evaluated: false, nextMetCodon: null, nextMetDisplay: null, hasUpstreamPathogenic: false, pvs1Code: null, gChrom: null, gStart: null, gEnd: null, clinVarUrl: null, html: null, error: null },
  splicePS1: { evaluated: false, spliceType: null, anchorBase: null, motifRange: null, querySpliceAI: null, candidates: [], ps1Code: null, clinVarUrl: null, html: null, error: null },
  // Raw ClinGen Allele Registry values (validation source — never feeds back to other modules)
  clingenAR: { gene: null, transcript: null, hgvsC: null, hgvsP: null, hg38String: null, hg19String: null, caId: null },
  // Raw MyVariant.info values (validation source)
  myvariant: { hg38String: null, hg19String: null, rsId: null },
  // Per-field two-source validation (see CLAUDE.md). Each field tracks raw source values,
  // whether 2+ sources agree (validated), whether validation is locked (no further updates),
  // which APIs responded (available), and any active conflict for ⚠ display.
  validation: {
    gene:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
    transcript:     { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
    hgvsC:          { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
    hgvsP:          { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
    hg38:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
    hg19:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
    exonIntron:     { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null, singleSource: true },
    rsId:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
    caId:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null, singleSource: true },
    vepConsequence: { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null, singleSource: true }
  },
  // BS2 decision-support (gnomAD homozygote/hemizygote + ClinVar-P/LP co-occurrence screen). Owns
  // this namespace exclusively; never auto-feeds the ACMG grid (display-only). See bs2_cooccurrence.js.
  bs2: { evaluated: false, geneConfig: null, pathwayA: null, error: null,
         screen: { status: 'idle', screened: 0, total: 0, targetInV2: null, hits: [] } },
  litvarCount: 0,
  selectedZygosity: 'Unknown',
  associatedConditions: [],
  uniprotDomains: null,
  uniprotAccession: null,
  omimSynopsisResults: null,  // per-disease OMIM synopsis match results
  // HPO phenotype-fit (ontology-aware likelihood ratio vs background frequency).
  // Owned exclusively by scorePhenotypeFit(); populated after OMIM synopsis match.
  phenotypeFit: { ready: false, results: [] }
};

window.currentVariantState = data;

let currentSearchGen = 0;
window.currentSearchGen = () => currentSearchGen;

let suggestedCodes = new Set();
let selectedCodes = new Set();
let debounceTimer = null;
let selectionOrder = [];

// ── EVENT HANDLERS ────────────────────────────────────────────────

function handleInput() {
  clearTimeout(debounceTimer);
  if (window.globalSearchController) {
    const abortErr = new Error('AbortError');
    abortErr.name = 'AbortError';
    window.globalSearchController.abort(abortErr);
  }
  window.globalSearchController = new AbortController();

  disableButtons();
  clearDash();
  document.getElementById('statusMsg').innerText = 'Analysing…';
  const gen = ++currentSearchGen;
  debounceTimer = setTimeout(() => decideAndFetch(gen), 800);
}

function handleGeneInput() {
  const g = document.getElementById('geneInput').value.trim().toUpperCase();
  if (!g) return;

  const raw = document.getElementById('rawInput').value.trim();
  if (raw) {
    // A variant is present — the gene field is only the VCEP hint and must not
    // disturb the in-flight variant analysis. Update gene + VCEP only.
    data.coords.gene = g;
    checkVCEP(g);
    return;
  }

  // Gene-only mode. A changed gene must wipe the previous gene's cards so no stale,
  // misleading data lingers: abort in-flight fetches, then clear the whole dashboard
  // before re-fetching (mirrors handleInput's behaviour for the variant field).
  clearTimeout(debounceTimer);
  if (window.globalSearchController) {
    const abortErr = new Error('AbortError');
    abortErr.name = 'AbortError';
    window.globalSearchController.abort(abortErr);
  }
  window.globalSearchController = new AbortController();
  disableButtons();
  clearDash();
  data.coords.gene = g;   // clearDash() resets coords.gene — restore it for checkVCEP
  checkVCEP(g);
  document.getElementById('statusMsg').innerText = 'Analysing…';
  const gen = ++currentSearchGen;
  debounceTimer = setTimeout(() => decideAndFetch(gen), 800);
}

// ── Phenotype edited AFTER a search — re-run only the phenotype-dependent analysis ──
// Bug fix: when the variant/gene is searched first and a phenotype is added afterward, the
// phenotype must still feed the HPO/OMIM match, the NCBI literature search, the GeneReviews
// focus, and the disease summary — WITHOUT re-validating/re-fetching the unchanged variant.
// (The variant & gene inputs trigger decideAndFetch; #hpoInput previously triggered nothing.)
// If no search has run yet (no gene), decideAndFetch will pick up the phenotype on first search.
let phenoRerunTimer = null;
function schedulePhenotypeRerun() {
  clearTimeout(phenoRerunTimer);
  phenoRerunTimer = setTimeout(runPhenotypeRerun, 500);
}
async function runPhenotypeRerun() {
  const gene = data.coords?.gene;
  if (!gene || gene === '-') return;              // nothing searched yet — decideAndFetch handles it
  const gen = window.currentSearchGen();          // reuse current generation; do NOT abort variant fetches
  const phenoRaw = document.getElementById('hpoInput')?.value.trim() || '';
  try {
    if (phenoRaw) { window._phenoParsePromise = parsePatientPhenotypes(phenoRaw); await window._phenoParsePromise; }
    else { data.ptPhenoTexts = []; data.ptPhenoGroups = []; window._phenoParsePromise = Promise.resolve(); }
  } catch (e) { if (e.name === 'AbortError') return; }
  if (gen !== window.currentSearchGen()) return;  // a new variant search superseded this phenotype edit

  if (data.ptPhenoTexts.length > 0) {
    // HPO gene match + OMIM disease match + phenotype-aware PubMed/notations (runPubMedOnly is fired inside)
    runPhenoMatch(gene).then(() => {
      if (gen !== window.currentSearchGen()) return;
      focusGeneReviewsByPhenotype(data.omimSynopsisResults, gen);
      focusGeneValidityByPhenotype(data.omimSynopsisResults, gen);
      refreshCurationSupport(gen);   // phenotype changed → re-match frequency table
      if (typeof populateVariantDiseaseRow === 'function' && document.getElementById('sumVariantDiseaseList'))
        populateVariantDiseaseRow(gene);
    }).catch(() => {});
  } else {
    // phenotype fully cleared → revert to gene-only literature and un-focus GeneReviews (show all)
    runPubMedOnly(gene);
    if (data.geneReviews) { data.geneReviews.focusNbk = null; renderGeneReviewsCard(); }
    if (data.geneValidity) { data.geneValidity.focusMondo = null; renderGeneValidityCard(); }
    refreshCurationSupport(gen);   // phenotype cleared → re-render concordance without phenotype matches
    if (typeof populateVariantDiseaseRow === 'function' && document.getElementById('sumVariantDiseaseList'))
      populateVariantDiseaseRow(gene);
  }
}

function handleCheckboxChange(k, on) {
  if (on) { if (!selectionOrder.includes(k)) selectionOrder.push(k); }
  else { selectionOrder = selectionOrder.filter(i => i !== k); }
  updateRanks();
}

function toggleChip(id, isPath) {
  if (selectedCodes.has(id)) { selectedCodes.delete(id); }
  else { selectedCodes.add(id); }
  renderChips(); updateBadge();
}

function clearACMG() {
  selectedCodes.clear(); renderChips(); updateBadge();
}

function disableButtons() {
  // Database launchers remain enabled for workflow efficiency (V81 change)
  document.querySelectorAll('.db-wrapper').forEach(w => w?.classList.remove('disabled'));
  const btnAll = document.getElementById('btnLaunchAll');
  if (btnAll) btnAll.disabled = false;
  const btnSel = document.getElementById('btnLaunchSel');
  if (btnSel) btnSel.disabled = false;
  ['vv', 'vep', 'mv', 'cv', 'sai'].forEach(d => setDot(d, ''));
  setAPIStatus('ucsc', '');
}

function clearDash() {
  ['dGene', 'dHgvs', 'dProtein', 'dHg38', 'dHg19', 'eManeTranscript', 'eExon', 'eRsId', 'eCAId', 'eConsequence',
    'mvPopmax', 'mvAllPop', 'mvEasPop', 'mvSasPop', 'mvNfePop', 'mvAfrPop', 'mvAmrPop', 'mvRevel', 'mvAlphaMissense', 'mvSpliceAI',
    'mvSpliceAI_AG', 'mvSpliceAI_AL', 'mvSpliceAI_DG', 'mvSpliceAI_DL',
    'mvSuggest', 'gnomadGeneConstraint', 'eRegionalOE', 'eRegionalRange', 'eHaploScore', 'eTriploScore', 'eGeneReviewsBody',
    'cvSig', 'cvStatus', 'cvStars', 'cvSubs', 'cvId', 'cvAcc',
    'pScore', 'pTerms', 'pLin', 'pLitVarMatches', 'pPubTiabRev', 'pPubTiabNon',
    'eInheritance'].forEach(id => {
      const el = document.getElementById(id); if (el) { el.innerText = '—'; el.innerHTML = '—'; if (el.classList.contains('dv')) el.className = 'dv'; }
    });
  const pbd = document.getElementById('pPhenoBreakdown'); if (pbd) pbd.innerHTML = '';
  const omimRow = document.getElementById('pOmimSynopsisRow'); if (omimRow) omimRow.style.display = 'none';
  const sps1 = document.getElementById('splicePS1Card'); if (sps1) sps1.style.display = 'none';
  const sps1c = document.getElementById('splicePS1Content'); if (sps1c) sps1c.innerHTML = '—';
  const txPicker = document.getElementById('rowTranscriptPicker'); if (txPicker) txPicker.style.display = 'none';
  const txDrop = document.getElementById('transcriptDropdown'); if (txDrop) txDrop.innerHTML = '';
  const txDetail = document.getElementById('transcriptDropdownDetail'); if (txDetail) txDetail.innerText = '';
  const mutalyzerRow = document.getElementById('rowMutalyzer'); if (mutalyzerRow) mutalyzerRow.style.display = 'none';
  const mutalyzerEl = document.getElementById('eMutalyzerHgvs'); if (mutalyzerEl) mutalyzerEl.innerText = '—';
  const multiTxRow = document.getElementById('rowPubMedMultiTx'); if (multiTxRow) multiTxRow.style.display = 'none';
  const multiTxEl = document.getElementById('pMultiTxTotal'); if (multiTxEl) multiTxEl.innerHTML = '—';
  const multiTxBd = document.getElementById('pMultiTxBreakdown'); if (multiTxBd) { multiTxBd.style.display = 'none'; multiTxBd.innerHTML = ''; }
  data.hpoMatchedTerms = new Set();
  data.hpoMatchResult = null;
  data.ptPhenoMeSH = {};

  const header = document.getElementById('conditions-toggle-header');
  if (header) header.style.display = 'none';
  const container = document.getElementById('associated-conditions-container');
  if (container) container.innerHTML = '';

  Object.assign(data, {
    coords:  { gene: '', hgvs: '-', chrom: '', pos38: '', pos19: '', ref: '', alt: '', strand: 1, hg38String: '', hg19String: '', hg19FromVV: false, primaryFromVV: false, userCoord: null },
    scores:  { revel: null, spliceAI: null, spliceAI_AG: null, spliceAI_AL: null, spliceAI_DG: null, spliceAI_DL: null, spliceAI_DP_AG: null, spliceAI_DP_AL: null, spliceAI_DP_DG: null, spliceAI_DP_DL: null, spliceAiManeMatch: true, alphaMissenseScore: null, alphaMissensePred: null },
    ensembl: { transcript: null, protein: null, vepConsequence: null, vepExon: null, vepTotalExons: null, vepAminoAcids: null, vepProteinStart: null, strand: 1, rsId: null, proteinId: null, allTranscripts: [], mutalyzerHgvs: null, _geneFromVep: null, _transcriptFromVep: null, _hgvscFromVep: null, _hgvspFromVep: null, _hg38FromVep: null },
    gnomad:  { popmax: null, detailed: null, viewMode: 'total', geneConstraint: null, regionalConstraint: null, regionalConstraintMsg: null, regionalOE: null, regionalRange: null, regionalPValue: null },
    dosage:  { haplo: null, triplo: null, source: null, haploId: null },
    geneValidity: { curations: [], gene: null, fetched: false, focusMondo: null },
    geneReviews: { gene: null, chapters: [], fetched: false, focusNbk: null },
    curationSupport: { mechanism: null, phenoMatches: [], chapterTitle: null, chapterUrl: null },
    clinvar: { sig: null, stars: 0, reviewStatus: null, subs: 0, uid: null, accession: null, interpretationId: undefined, evidenceFlags: null, ps1Suggested: null, pm5Suggested: null, evidenceRepoUuid: null, historicalProteinNotations: [], detail: null },
    litvar:  { combinedTotal: 0, uniqueNotations: [], _lastKey: null },
    startLoss: { evaluated: false, nextMetCodon: null, nextMetDisplay: null, hasUpstreamPathogenic: false, pvs1Code: null, gChrom: null, gStart: null, gEnd: null, clinVarUrl: null, html: null, error: null },
    splicePS1: { evaluated: false, spliceType: null, anchorBase: null, motifRange: null, querySpliceAI: null, candidates: [], ps1Code: null, clinVarUrl: null, html: null, error: null },
    clingenAR: { gene: null, transcript: null, hgvsC: null, hgvsP: null, hg38String: null, hg19String: null, caId: null },
    myvariant: { hg38String: null, hg19String: null, rsId: null },
    validation: {
      gene:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
      transcript:     { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
      hgvsC:          { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
      hgvsP:          { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
      hg38:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
      hg19:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
      exonIntron:     { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null, singleSource: true },
      rsId:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null },
      caId:           { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null, singleSource: true },
      vepConsequence: { value: null, sources: {}, validated: false, locked: false, available: [], conflict: null, singleSource: true }
    },
    ptPhenoTexts: [], ptPhenoGroups: [], selectedHpoTerms: {},
    litvarPMIDs: null, uniprotDomains: null, uniprotAccession: null,
    associatedConditions: [],
    caId: null, vcepNarrative: null, civicEvidenceCount: 0, maveData: null,
    vcepGuideline: null,
    hpoMatchedTerms: new Set(), hpoMatchResult: null, ptPhenoMeSH: {},
    omimSynopsisResults: null,
    bs2: { evaluated: false, geneConfig: null, pathwayA: null, error: null,
           screen: { status: 'idle', screened: 0, total: 0, targetInV2: null, hits: [] } }
  });

  // Reset BS2 decision-support panel (it owns #bs2Card; render only its own state, Rule 9)
  const bs2Card = document.getElementById('bs2Card'); if (bs2Card) bs2Card.style.display = 'none';

  // Reset Population Frequency UI
  const popContainer = document.getElementById('populationBreakdown');
  if (popContainer) popContainer.innerHTML = '';
  document.querySelectorAll('.toggle-btn').forEach(btn => {
    btn.classList.toggle('active', btn.innerText.toLowerCase() === 'total');
  });
  ['exomePopmax', 'genomePopmax', 'mvPopmax'].forEach(id => {
    const el = document.getElementById(id); if (el) el.innerText = '—';
  });
  ['exomeFilterBadge', 'genomeFilterBadge'].forEach(id => {
    const el = document.getElementById(id); if (el) el.innerHTML = '';
  });
  ['maveDBCardContainer', 'domainAnnotationContainer', 'codonViewerContainer'].forEach(id => {
    const el = document.getElementById(id); if (el) el.style.display = 'none';
  });
  const vdc = document.getElementById('variantDistributionContainer');
  if (vdc) { vdc.innerHTML = ''; vdc.style.display = 'none'; }
  suggestedCodes.clear(); renderChips(); updateBadge();
  renderVCEPBanner();
}

async function decideAndFetch(gen) {
  const raw = document.getElementById('rawInput').value.trim();
  const phenoRaw = document.getElementById('hpoInput').value.trim();
  const geneRaw = document.getElementById('geneInput').value.trim().toUpperCase();
  if (geneRaw) data.coords.gene = geneRaw;

  // Start phenotype resolution immediately but don't block the variant pipeline.
  // runPhenoMatch will await this promise before using data.ptPhenoGroups.
  // For gene-only mode (no raw variant), we still await it directly since there
  // is nothing else to overlap with.
  if (phenoRaw) {
    window._phenoParsePromise = parsePatientPhenotypes(phenoRaw);
  } else {
    data.ptPhenoTexts = []; data.ptPhenoGroups = [];
    window._phenoParsePromise = Promise.resolve();
  }

  if (!raw) {
    try {
      await window._phenoParsePromise; // gene-only: phenotype is the critical path
    } catch (e) {
      if (e.name === 'AbortError') return; // superseded by a newer search — abort cleanly
    }
    if (gen !== currentSearchGen) return; // a newer gene/variant search superseded this one
    if (geneRaw) runGeneOnlyMode(gen);
    else document.getElementById('statusMsg').innerText = 'Enter a variant or gene to begin analysis.';
    return;
  }

  // [FIX] 1. Check for Local VCF String Format FIRST using the unmodified 'raw' input
  if (raw.includes(';') && raw.includes('c.')) {
    if (parseLocalVCFString(raw)) {
      document.getElementById('statusMsg').innerText = 'Local Variant Parsed Successfully.';
      enableWrappers(['franklin', 'gnomadv4', 'spliceai', 'liftover', 'clingen', 'clinvar', 'omim', 'gr', 'scholar', 'decipher', 'hgmd', 'gtex', 'mastermind', 'gnomadv2']);
      document.getElementById('btnLaunchAll').disabled = false;
      triggerDownstreamAPIs(gen);
      return;
    }
  }

  // 2. V81 Robust Sanitization: Run ONLY if the local parser didn't catch it
  let clean = raw.replace(/^["'“”'']+|["'“”'']+$/g, '')
    .replace(/\s*\([^)]+\):/g, ':') // Robust removal of (GENE) before colon
    .split(/\s*\(?p\./)[0]           // Splits at protein consequence
    .trim();

  // Step 1: Enforce VCF Hyphenation (Convert "chr7 140753336 A T" -> "7-140753336-A-T")
  if (/\s+/.test(clean) && !clean.includes(':')) {
    // Strip 'chr' prefix first to standardize
    let noChr = clean.replace(/^chr/i, '');
    const parts = noChr.split(/\s+/);

    // Validate it's a 4-part VCF coordinate
    if (parts.length === 4 && /^(\d+|X|Y|MT)$/i.test(parts[0]) && /^\d+$/.test(parts[1])) {
      clean = parts.join('-');
    }
  }

  // 3. API Routing Logic
  const apiMode = document.getElementById('apiSelector')?.value || 'VV';

  if (apiMode === 'VV') {
    try {
      await fetchVV(clean, gen);
    } catch (err) {
      if (err.name === 'AbortError') return; // Silent abort
      console.warn("VariantValidator API Failed:", err.message);
      document.getElementById('statusMsg').innerHTML = `⚠️ <b style="color:var(--amber)">VV API Unavailable. Falling back to Ensembl VEP...</b>`;
      try {
        await runEnsemblPrimary(clean, gen);
      } catch (e) {
        if (e.name === 'AbortError') return;
        console.warn("Ensembl VEP Failed:", e.message);
        try {
          await fetchClinGenFallback(clean, gen);
        } catch (e3) {
          if (e3.name === 'AbortError') return;
          document.getElementById('statusMsg').innerHTML = `❌ <b style="color:var(--red)">All annotation APIs failed: ${e3.message}</b>`;
          return;
        }
      }
      triggerDownstreamAPIs(gen);
      return;
    }
    triggerDownstreamAPIs(gen);
    return;
  } else if (apiMode === 'ENSEMBL') {
    // Forced Ensembl Mode
    try {
      await runEnsemblPrimary(clean, gen);
    } catch (err) {
      if (err.name === 'AbortError') return;
      throw err;
    }
  } else if (apiMode === 'CLINGEN') {
    // Forced ClinGen Mode: coords + CAid from AR, then VEP for gene/transcript/consequence
    try {
      await fetchClinGenFallback(clean, gen);
    } catch (err) {
      if (err.name === 'AbortError') return;
      throw err;
    }
    // VEP enrichment using resolved hg38 coords (best-effort — preserves caId)
    if (data.coords.hg38String) {
      const savedCaId = data.caId;
      try {
        await runEnsemblPrimary(data.coords.hg38String, gen);
        document.getElementById('statusMsg').innerHTML =
          `✅ <b style="color:var(--teal)">ClinGen + Ensembl VEP — coords from Allele Registry, annotation from VEP</b>`;
      } catch (_) {
        document.getElementById('statusMsg').innerHTML =
          `⚠️ <b style="color:var(--amber)">ClinGen coords only — Ensembl VEP unavailable for full annotation</b>`;
      }
      data.caId = savedCaId; // restore CAid overwritten by VEP
      if (data.caId) document.getElementById('eCAId').innerText = data.caId;
    }
  }

  // 4. Fire Downstream
  triggerDownstreamAPIs(gen);
}

function evaluateACMG() {
  // Preserve start-loss codes if they were already added by evaluateStartLoss
  const slCodes = [...suggestedCodes].filter(c => c.startsWith('PVS1_'));
  suggestedCodes.clear();
  slCodes.forEach(c => suggestedCodes.add(c));
  if (data.clinvar.ps1Suggested) suggestedCodes.add(data.clinvar.ps1Suggested);
  if (data.clinvar.pm5Suggested) suggestedCodes.add(data.clinvar.pm5Suggested);
  if (data.vepConsequence && NULL_CONS.some(c => data.vepConsequence.includes(c))) {
    // ClinGen SVI: Skip standard PVS1 for start-loss (it has its own automated check)
    if (!data.vepConsequence.includes('start_lost')) {
      suggestedCodes.add('PVS1');
    }
  }
  if (data.gnomad.popmax !== null && data.gnomad.popmax > 0.05) suggestedCodes.add('BA1');
  if (data.gnomad.popmax === null || data.gnomad.popmax < 0.0001) suggestedCodes.add('PM2');
  let usedInSilico = false;
  if (data.scores.revel !== null) {
    usedInSilico = true; // REVEL available — AlphaMissense must not act as fallback even if no code fires
    const r = data.scores.revel;
    if (r > 0 && r <= 0.003) { suggestedCodes.add('BP4_VS'); }
    else if (r > 0.003 && r <= 0.016) { suggestedCodes.add('BP4_S'); }
    else if (r > 0.016 && r <= 0.183) { suggestedCodes.add('BP4_M'); }
    else if (r > 0.183 && r <= 0.290) { suggestedCodes.add('BP4'); }
    else if (r >= 0.932) { suggestedCodes.add('PP3_S'); }
    else if (r >= 0.773) { suggestedCodes.add('PP3_M'); }
    else if (r >= 0.644) { suggestedCodes.add('PP3'); }
  }
  if (data.scores.alphaMissenseScore !== null && !usedInSilico) {
    const am = data.scores.alphaMissenseScore;
    if (am >= 0.984) { suggestedCodes.add('PP3_S'); usedInSilico = true; }
    else if (am >= 0.869) { suggestedCodes.add('PP3_M'); usedInSilico = true; }
    else if (am >= 0.761) { suggestedCodes.add('PP3'); usedInSilico = true; }
    else if (am <= 0.073) { suggestedCodes.add('BP4_VS'); usedInSilico = true; }
    else if (am <= 0.147) { suggestedCodes.add('BP4_S'); usedInSilico = true; }
    else if (am <= 0.331) { suggestedCodes.add('BP4_M'); usedInSilico = true; }
  }
  if (data.scores.spliceAI !== null && data.scores.spliceAI > 0.20 && !usedInSilico) suggestedCodes.add('PP3');

  // MaveDB Functional Assay Evaluation (PS3 / BS3)
  if (data.maveData && data.maveData.interpretation) {
    const interp = data.maveData.interpretation.toLowerCase();
    if (interp.includes('damaging') || interp.includes('pathogenic') || interp.includes('loss of function')) {
      suggestedCodes.add('PS3');
    } else if (interp.includes('normal') || interp.includes('benign') || interp.includes('wild type')) {
      suggestedCodes.add('BS3');
    }
  }

  // PM1: variant in critical/high-tier functional domain (UniProt)
  if (data.vepProteinStart && data.uniprotDomains?.length) {
    const pos = parseInt(data.vepProteinStart, 10);
    const pmHit = data.uniprotDomains.find(
      d => pos >= d.start && pos <= d.end && (d.tier === 'critical' || d.tier === 'high')
    );
    if (pmHit) suggestedCodes.add('PM1');
  }

  // BP7: synonymous with no predicted splice impact
  if (data.vepConsequence?.includes('synonymous_variant')) {
    if (data.scores.spliceAI === null || data.scores.spliceAI < 0.1) suggestedCodes.add('BP7');
  }

  // PP2: missense in gene with low benign missense rate (mis_z > 3.09)
  if (data.vepConsequence?.includes('missense_variant')) {
    if (data.geneConstraint?.mis_z != null && data.geneConstraint.mis_z > 3.09) suggestedCodes.add('PP2');
  }

  // PM4: in-frame indel or stop-loss
  if (data.vepConsequence?.includes('inframe_insertion') || data.vepConsequence?.includes('inframe_deletion') || data.vepConsequence?.includes('stop_lost')) {
    suggestedCodes.add('PM4');
  }

  if (selectedCodes.size === 0) suggestedCodes.forEach(c => selectedCodes.add(c));
  const codes = [...suggestedCodes];
  const suggestEl = document.getElementById('mvSuggest');
  if (data.startLoss?.html) {
    suggestEl.innerHTML = data.startLoss.html;
  } else {
    suggestEl.textContent = codes.length ? codes.join(', ') : 'None';
  }
  renderChips(); updateBadge();
}

// ── BS2 DECISION SUPPORT (gnomAD homozygote/hemizygote + co-occurrence) ─────────
// Adapter layer between the app's global `data.*` state and the PURE window.BS2 module
// (src/bs2_cooccurrence.js). Per CLAUDE.md Rule 2, the module never reads `data.*`; this
// adapter GATHERS inputs from `data.*` and passes them explicitly. Owns `data.bs2.*` only.
// Display-only: it NEVER writes suggestedCodes/selectedCodes (Rule 10) — the curator toggles
// the BS2 chip manually if they concur with the verdict.

// Genes with well-known genotype–phenotype / hypomorph complexity, where tolerated in-trans
// mild-allele combinations occur — BS2 needs a manual override for these (module Gate 1).
const BS2_HYPOMORPHIC_GENES = new Set([
  'GALC', 'GBA', 'GBA1', 'CFTR', 'HFE', 'GAA', 'HEXA', 'SMPD1', 'PAH', 'DHCR7',
  'ATP7B', 'BTD', 'ACADM', 'NPC1', 'ABCA4'
]);

const bs2Esc = s => String(s ?? '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

// Map a free-text inheritance/MOI string → the module's enum. null if not recessive/dominant.
function bs2MapInheritance(str) {
  if (!str) return null;
  const s = String(str).toUpperCase();
  if (/X-?LINKED RECESSIVE|XLR\b/.test(s)) return 'XLR';
  if (/X-?LINKED(?! DOMINANT)|\bXL\b/.test(s) && !/DOMINANT/.test(s)) return 'XLR';
  if (/AUTOSOMAL RECESSIVE|RECESSIVE|BIALLELIC|\bAR\b/.test(s)) return 'AR';
  if (/AUTOSOMAL DOMINANT|DOMINANT|MONOALLELIC|\bAD\b/.test(s)) return 'AD';
  return null;
}

// Map a free-text age-of-onset clause → the module's onset enum.
function bs2MapOnset(str) {
  if (!str) return null;
  const s = String(str).toLowerCase();
  if (/neonat|congenital|prenat|antenat|at birth|in utero|newborn|infan/.test(s)) return 'pediatric_severe';
  if (/child|juvenile|pediatr|paediatr|adolescen/.test(s)) return 'pediatric';
  if (/adult|late|older|aging|elderly/.test(s)) return 'adult';
  return null;
}

// Classify the QUERY variant region for the module's co-occurrence eligibility gate (Pathway B).
// Pure: classify a consequence string into the module's co-occurrence region scope.
function bs2RegionOf(consequence) {
  const c = (consequence || '').toLowerCase();
  if (!c) return undefined;
  if (c.includes('splice')) return 'splice';
  if (c.includes('utr')) return 'utr';
  if (/missense|synonymous|stop|frameshift|inframe|start_lost|protein_altering|coding|initiator/.test(c)) return 'coding';
  return 'other';
}
function bs2RegionFromConsequence() { return bs2RegionOf(data.ensembl?.vepConsequence); }

/**
 * deriveBS2GeneConfig() — PURE-ish: reads gene/disease context from `data.*` and builds the
 * GeneDiseaseConfig the module expects. Prefers ClinGen gene-validity curations (authoritative
 * MOI + onset), falling back to associated-conditions inheritance. Returns null if no usable
 * recessive/dominant inheritance is known yet.
 * Inputs read: data.coords.gene, data.coords.chrom, data.geneValidity.curations[], data.associatedConditions[]
 */
function deriveBS2GeneConfig() {
  const gene = (data.coords?.gene || '').toUpperCase();
  if (!gene) return null;

  const curations = data.geneValidity?.curations || [];
  const conditions = data.associatedConditions || [];

  // Collect candidate inheritances (curations first, then conditions).
  const inhCandidates = [];
  curations.forEach(c => { const m = bs2MapInheritance(c.moi); if (m) inhCandidates.push({ inh: m, disease: c.disease, onset: c.onset }); });
  conditions.forEach(c => { const m = bs2MapInheritance(c.inheritance); if (m) inhCandidates.push({ inh: m, disease: c.diseaseName, onset: null }); });

  if (!inhCandidates.length) return null;

  // BS2 biallelic pathways need recessive inheritance. Prefer XLR for X-chromosome variants,
  // else AR; fall back to whatever single mode is present (incl. AD → card stays hidden).
  const onX = String(data.coords?.chrom || '').toUpperCase() === 'X';
  const has = inh => inhCandidates.find(c => c.inh === inh);
  let chosen = null;
  if (onX && has('XLR')) chosen = has('XLR');
  else if (has('AR')) chosen = has('AR');
  else if (has('XLR')) chosen = has('XLR');
  else chosen = inhCandidates[0];

  // Earliest / most severe onset across ALL curations (auto-assume severe-AR intent).
  const rank = { pediatric_severe: 0, pediatric: 1, adult: 2 };
  let bestOnset = null, bestRank = Infinity;
  curations.forEach(c => {
    const o = bs2MapOnset(c.onset);
    if (o && (rank[o] ?? 9) < bestRank) { bestRank = rank[o] ?? 9; bestOnset = o; }
  });
  const onset = bestOnset || 'variable';

  // Auto-assume full penetrance ONLY for recessive + paediatric onset; otherwise unknown
  // (which Gate 1 rejects, surfacing the assumption to the curator instead of firing silently).
  const isPaed = onset === 'pediatric_severe' || onset === 'pediatric';
  const isRecessive = chosen.inh === 'AR' || chosen.inh === 'XLR';
  const penetrance = (isRecessive && isPaed) ? 'full' : 'unknown';

  return {
    gene,
    disease: chosen.disease || '(disease not resolved)',
    inheritance: chosen.inh,
    penetrance,
    onset,
    hypomorphicComplexity: BS2_HYPOMORPHIC_GENES.has(gene),
  };
}

/**
 * evaluateBS2(gen) — Pathway A adapter. Builds the observed-biallelic observation from gnomAD
 * homozygote (AR) / hemizygote (XLR) counts and evaluates it. Generation-guarded; idempotent.
 * Inputs read: data.gnomad.detailed.total (hom/hemi/filters/ac), data.gnomad.popmax, data.coords.
 * Output: data.bs2.geneConfig, data.bs2.pathwayA; UI #bs2Card. Failure: stays silent, card hidden.
 */
function evaluateBS2(gen) {
  if (typeof gen === 'number' && gen !== window.currentSearchGen()) return;
  if (!window.BS2 || !data.gnomad?.detailed?.total) return;

  const config = deriveBS2GeneConfig();
  data.bs2.geneConfig = config;
  data.bs2.evaluated = true;

  // AD-only / unknown inheritance → biallelic BS2 not applicable; hide the card.
  if (!config || config.inheritance === 'AD') { data.bs2.pathwayA = null; renderBS2Panel(); return; }

  const total = data.gnomad.detailed.total;
  const mode = config.inheritance === 'XLR' ? 'hemizygous' : 'homozygous';
  const count = mode === 'hemizygous' ? (total.hemi || 0) : (total.hom || 0);
  const filters = total.filters;
  const qcPass = !(Array.isArray(filters) && filters.length > 0);

  // Allele-singleton semantics: a homozygote carries ≥2 alt copies and is NEVER an allele
  // singleton (its small-N weakness is already captured by the Strong-vs-Moderate count
  // threshold). A single hemizygote IS one alt copy → a genuine singleton (lower confidence).
  const isSingleton = mode === 'hemizygous' && count === 1;

  const query = {
    id: data.coords?.hg38String || config.gene,
    afGlobal: typeof data.gnomad.popmax === 'number' ? data.gnomad.popmax : undefined,
    region: bs2RegionFromConsequence(),
    inGnomadExomes: (data.gnomad.detailed.exome?.ac || 0) > 0,
  };
  const observation = { mode, individualCount: count, isSingleton, qcPass };

  try {
    data.bs2.pathwayA = window.BS2.evaluateBenignObservation({ query, observation, geneConfig: config });
  } catch (e) {
    data.bs2.pathwayA = null;
    data.bs2.error = 'Pathway A evaluation failed: ' + (e?.message || e);
  }
  renderBS2Panel();
}

// Guarded entry point fired from multiple data-resolution points (gnomAD ready, conditions ready,
// onset enrichment done). Idempotent — recomputes from whatever is currently in `data`.
function maybeEvaluateBS2(gen) {
  try {
    if (typeof gen === 'number' && gen !== window.currentSearchGen()) return;
    if (!window.BS2 || !data.gnomad?.detailed?.total) return;
    evaluateBS2(gen);
    maybeAutoScreenBS2(gen);   // Pathway B: auto-detect ClinVar-P/LP co-occurrence (AR + variant present)
  } catch { /* non-blocking, Rule 3 */ }
}

// Tracks which search generation has already auto-screened, so the heavy ClinVar-P/LP co-occurrence
// screen fires AT MOST ONCE per search (maybeEvaluateBS2 is called from several resolution points).
let _bs2AutoScreenGen = -1;

/**
 * maybeAutoScreenBS2(gen) — automatically launches the Pathway B co-occurrence screen once the gene
 * is known to be autosomal-recessive and a query variant is present. Cost is bounded: the screen's
 * own first step is a single gnomAD-v2 presence check on the target, so a rare/novel VUS (not in v2)
 * short-circuits after ~1 call without fetching/screening the gene's ClinVar P/LP list.
 */
function maybeAutoScreenBS2(gen) {
  if (typeof gen !== 'number') gen = window.currentSearchGen();
  if (gen === _bs2AutoScreenGen) return;                          // already auto-screened this search
  const cfg = data.bs2.geneConfig;
  if (!cfg || cfg.inheritance !== 'AR') return;                   // co-occurrence BS2 is AR-only
  if (!(data.coords?.hg19String || '').replace(/^chr/i, '')) return; // need a GRCh37 target id
  if (data.bs2.screen && data.bs2.screen.status !== 'idle') return;  // a screen is already running/done
  _bs2AutoScreenGen = gen;
  screenBS2Cooccurrence(gen);   // async, self-rendering, generation-guarded internally
}

// Hard cap on co-occurrence calls per screen (bounds gnomAD load; candidates are pre-sorted
// highest-confidence first, so the cap keeps the most informative partners).
const BS2_SCREEN_CAP = 150;

// Minimal concurrency pool — runs `worker` over `items` at most `concurrency` at a time,
// invoking onProgress after each. Never rejects (per-item errors are the worker's concern).
async function bs2Pool(items, worker, concurrency, onProgress) {
  let i = 0;
  async function runner() {
    while (i < items.length) {
      const idx = i++;
      try { await worker(items[idx], idx); } catch { /* per-item, swallowed */ }
      if (onProgress) onProgress();
    }
  }
  await Promise.all(Array.from({ length: Math.min(concurrency, items.length) }, runner));
}

/**
 * screenBS2Cooccurrence(gen) — Pathway B: AUTOMATIC screen. Pulls the gene's ClinVar P/LP variants
 * that are present in gnomAD v2, then tests each for co-occurrence (inferred compound-het) with the
 * TARGET variant in healthy gnomAD individuals, evaluating BS2 for every co-occurring pair.
 * Curator-triggered (one click) to bound gnomAD API load; the screening itself is fully automatic.
 * Inputs read: data.bs2.geneConfig, data.coords.hg19String, data.gnomad.popmax, data.ensembl.vepConsequence.
 * Output: data.bs2.screen.{status, screened, total, targetInV2, hits[]}; UI #bs2PathwayB. Errors → data.bs2.error.
 */
async function screenBS2Cooccurrence(gen) {
  const btn = document.getElementById('bs2ScreenBtn');
  data.bs2.error = null;
  data.bs2.screen = { status: 'running', screened: 0, total: 0, targetInV2: null, hits: [] };

  const config = data.bs2.geneConfig || deriveBS2GeneConfig();
  if (!config || config.inheritance !== 'AR') {
    data.bs2.error = 'Co-occurrence screening applies to autosomal-recessive genes only.';
    data.bs2.screen.status = 'error'; renderBS2Panel(); return;
  }
  const targetV2 = (data.coords?.hg19String || '').replace(/^chr/i, '');
  if (!/^[\dXYMT]+-\d+-[ACGT]+-[ACGT]+$/i.test(targetV2)) {
    data.bs2.error = 'Target variant has no usable GRCh37 coordinate (needed for gnomAD v2 co-occurrence).';
    data.bs2.screen.status = 'error'; renderBS2Panel(); return;
  }

  if (btn) btn.disabled = true;
  renderBS2Panel();
  try {
    // Gate 1: the target itself must be in gnomAD v2 exomes, else no pair can co-occur.
    const inV2 = await fetchVariantInGnomadV2Exomes(targetV2);
    if (gen != null && gen !== window.currentSearchGen()) return;
    data.bs2.screen.targetInV2 = inV2;
    if (!inV2) { data.bs2.screen.status = 'done'; renderBS2Panel(); return; }

    // Gate 2: gather the gene's ClinVar P/LP partners that are in gnomAD v2 (one call).
    let candidates = await fetchGeneClinVarPLP(config.gene, targetV2);
    if (gen != null && gen !== window.currentSearchGen()) return;
    candidates = candidates
      .filter(c => ['coding', 'splice', 'utr'].includes(bs2RegionOf(c.consequence)))
      .slice(0, BS2_SCREEN_CAP);
    data.bs2.screen.total = candidates.length;
    renderBS2Panel();

    const queryAfGlobal = typeof data.gnomad.popmax === 'number' ? data.gnomad.popmax : undefined;
    const queryRegion = bs2RegionFromConsequence();
    let lastRender = 0;

    // Screen each candidate for co-occurrence; a co-occurring pair (≥1 double-het) is evaluated.
    await bs2Pool(candidates, async (c) => {
      if (gen != null && gen !== window.currentSearchGen()) return;
      let co = null;
      try { co = await fetchGnomadCooccurrence(targetV2, c.variantId, gen); } catch { co = null; }
      if (!co || co.individualCount < 1) return;   // never seen together → not informative

      const partner = { id: c.variantId, classification: c.classification, region: bs2RegionOf(c.consequence), inGnomadExomes: true };
      const query = { id: targetV2, afGlobal: queryAfGlobal, region: queryRegion, inGnomadExomes: true };
      const observation = { mode: 'compound_het_inferred', individualCount: co.individualCount, pTrans: co.pTrans, isSingleton: co.isSingleton, qcPass: true };
      const verdict = window.BS2.evaluateBenignObservation({ query, partner, observation, geneConfig: config });
      data.bs2.screen.hits.push({ partnerId: c.variantId, classification: c.classification, goldStars: c.goldStars, pTrans: co.pTrans, individualCount: co.individualCount, isSingleton: co.isSingleton, verdict });
    }, 4, () => {
      data.bs2.screen.screened++;
      const now = Date.now();
      if (now - lastRender > 250) { lastRender = now; renderBS2Panel(); }  // throttled progress
    });
    if (gen != null && gen !== window.currentSearchGen()) return;

    // Rank: applicable first, then strongest (most-negative points), then highest pTrans.
    data.bs2.screen.hits.sort((a, b) =>
      (Number(b.verdict.applicable) - Number(a.verdict.applicable)) ||
      ((a.verdict.bayesianPoints || 0) - (b.verdict.bayesianPoints || 0)) ||
      ((b.pTrans || 0) - (a.pTrans || 0)));
    data.bs2.screen.status = 'done';
  } catch (e) {
    data.bs2.error = 'Co-occurrence screen failed: ' + (e?.message || e);
    data.bs2.screen.status = 'error';
  } finally {
    if (btn) btn.disabled = false;
  }
  renderBS2Panel();
}

// Render one verdict object (from window.BS2.evaluateBenignObservation) as an HTML block.
function bs2RenderVerdict(v) {
  if (!v) return '<span style="color:var(--dim); font-size:0.72rem;">—</span>';
  const head = v.applicable
    ? `<span style="color:var(--teal); font-weight:800;">BS2 ${bs2Esc(v.strength)}</span>` +
      `<span style="color:var(--dim);"> · ${v.bayesianPoints} pts · confidence ${bs2Esc(v.confidence)}</span>`
    : `<span style="color:var(--dim); font-weight:700;">Not applicable</span>`;
  const rationale = (v.rationale || []).map(r =>
    `<li style="margin:2px 0;">${bs2Esc(r)}</li>`).join('');
  const flags = (v.flags || []).map(f =>
    `<li style="margin:2px 0; color:var(--amber);">⚠ ${bs2Esc(f)}</li>`).join('');
  return (
    `<div style="font-size:0.72rem;">${head}</div>` +
    (rationale ? `<ul style="margin:4px 0 0 16px; padding:0; font-size:0.68rem; color:var(--dim);">${rationale}</ul>` : '') +
    (flags ? `<ul style="margin:4px 0 0 16px; padding:0; font-size:0.68rem;">${flags}</ul>` : '') +
    (v.applicable ? `<div style="font-size:0.62rem; color:var(--amber); margin-top:4px; font-weight:700;">Requires expert review — not auto-applied.</div>` : '')
  );
}

/**
 * renderBS2Panel() — self-contained (CLAUDE.md Rule 9): reads ONLY data.bs2.* and writes ONLY
 * into #bs2Card and its children. Hidden when no recessive config is known (AD/unknown).
 */
function renderBS2Panel() {
  const card = document.getElementById('bs2Card');
  if (!card) return;
  const cfg = data.bs2.geneConfig;

  if (!cfg || cfg.inheritance === 'AD') { card.style.display = 'none'; return; }
  card.style.display = '';

  const ctx = document.getElementById('bs2Context');
  if (ctx) {
    const penLabel = cfg.penetrance === 'full' ? 'full (assumed)' : bs2Esc(cfg.penetrance);
    ctx.innerHTML =
      `<b style="color:#fff;">${bs2Esc(cfg.gene)}</b> · ${bs2Esc(cfg.disease)}<br>` +
      `Inheritance: <b>${bs2Esc(cfg.inheritance)}</b> · Onset: <b>${bs2Esc(cfg.onset)}</b> · ` +
      `Penetrance: <b>${penLabel}</b>` +
      (cfg.hypomorphicComplexity ? ` · <span style="color:var(--amber);">hypomorph-complex gene</span>` : '');
  }

  const aEl = document.getElementById('bs2PathwayA');
  if (aEl) {
    const total = data.gnomad?.detailed?.total || {};
    const cnt = cfg.inheritance === 'XLR' ? (total.hemi || 0) : (total.hom || 0);
    const lbl = cfg.inheritance === 'XLR' ? 'hemizygotes' : 'homozygotes';
    aEl.innerHTML =
      `<div class="dk" style="font-size:0.72rem; margin-bottom:4px;">Pathway A — observed healthy ${lbl} ` +
      `<span style="color:var(--dim);">(gnomAD: ${cnt})</span></div>` +
      bs2RenderVerdict(data.bs2.pathwayA);
  }

  const bEl = document.getElementById('bs2PathwayB');
  if (bEl) {
    const sc = data.bs2.screen || { status: 'idle', screened: 0, total: 0, hits: [] };
    if (data.bs2.error) {
      bEl.innerHTML = `<span style="color:var(--red); font-size:0.68rem;">${bs2Esc(data.bs2.error)}</span>`;
    } else if (sc.status === 'idle') {
      bEl.innerHTML = `<span style="color:var(--dim); font-size:0.66rem;">Screen this gene's ClinVar P/LP variants for co-occurrence with the target in healthy gnomAD individuals (inferred compound-het).</span>`;
    } else if (sc.status === 'running') {
      const denom = sc.total ? sc.total : '…';
      bEl.innerHTML = `<span style="color:var(--dim); font-size:0.66rem;">Screening ClinVar P/LP partners… ${sc.screened}/${denom} · ${sc.hits.length} co-occurring hit(s)</span>`;
    } else if (sc.targetInV2 === false) {
      bEl.innerHTML = `<span style="color:var(--dim); font-size:0.66rem;">Target variant is not in gnomAD v2 exomes — co-occurrence (v2-only) cannot be computed for it.</span>`;
    } else if (!sc.hits.length) {
      bEl.innerHTML = `<span style="color:var(--dim); font-size:0.66rem;">No co-occurring P/LP partner found in gnomAD (screened ${sc.screened} candidate${sc.screened === 1 ? '' : 's'}).</span>`;
    } else {
      bEl.innerHTML =
        `<div style="font-size:0.66rem; color:var(--dim); margin-bottom:4px;">${sc.hits.length} co-occurring P/LP partner(s) — screened ${sc.screened}:</div>` +
        sc.hits.map(bs2RenderHit).join('');
    }
  }
}

// Render one co-occurrence screen hit (partner variant + its inferred-compound-het BS2 verdict).
function bs2RenderHit(h) {
  const link = `https://gnomad.broadinstitute.org/variant/${encodeURIComponent(h.partnerId)}?dataset=gnomad_r2_1`;
  const stars = h.goldStars ? ' · ' + '★'.repeat(h.goldStars) : '';
  return `<div style="border-top:1px solid rgba(255,255,255,0.08); margin-top:6px; padding-top:6px;">` +
    `<div style="font-size:0.68rem;">` +
    `<a href="${link}" target="_blank" rel="noopener" style="color:var(--teal); text-decoration:none; font-weight:600;">${bs2Esc(h.partnerId)}</a>` +
    `<span style="color:var(--dim);"> · ${bs2Esc(h.classification)}${stars} · pTrans=${bs2Esc(h.pTrans)} · ${bs2Esc(h.individualCount)} indiv.${h.isSingleton ? ' · singleton' : ''}</span></div>` +
    bs2RenderVerdict(h.verdict) + `</div>`;
}

function getUrls() {
  const c = data.coords.userCoord || {};
  const fID = data.coords.userCoord ? `chr${c.chrom}-${c.pos}-${c.ref}-${c.alt}-hg38` : `${data.coords.hg38String}-hg38`;
  const gV4 = data.coords.userCoord ? `${c.chrom}-${c.pos}-${c.ref}-${c.alt}` : `${data.coords.chrom}-${data.coords.pos38}-${data.coords.ref}-${data.coords.alt}`;
  const broad = data.coords.userCoord ? `chr${c.chrom}-${c.pos}%20${c.ref}%3E${c.alt}` : `chr${data.coords.chrom}-${data.coords.pos38}%20${data.coords.ref}%3E${data.coords.alt}`;
  const hgvsTerm = (data.ensembl.transcript && data.coords.hgvs && data.coords.hgvs !== '-') ? `${data.ensembl.transcript}:${data.coords.hgvs}` : data.coords.hgvs;
  const cvTerm = (data.coords.gene && hgvsTerm && hgvsTerm !== '-')
    ? `(${data.coords.gene}[Gene]) AND "${hgvsTerm}"[VARNAME]`
    : data.coords.gene;

  const gV2str = data.coords.hg19String ? data.coords.hg19String.replace('chr', '') : '';
  const selectedTerms = Object.entries(data.selectedHpoTerms || {})
    .filter(([, isSelected]) => isSelected)
    .map(([term]) => term);
  const schQ = `intitle:${data.coords.gene}` + (selectedTerms.length ? `("${selectedTerms.join('" OR "')}")` : '');
  const _mmIsSplice = /\d[+\-]\d/.test(data.coords.hgvs || '');
  const _mmAA1L = toOneLetterAA(data.ensembl.protein);
  const mmHgvs = (!_mmIsSplice && data.coords.gene && _mmAA1L)
    ? encodeURIComponent(`${data.coords.gene}:p.${_mmAA1L}`)
    : (data.ensembl.transcript && data.coords.hgvs && data.coords.hgvs !== '-')
      ? encodeURIComponent(`${data.ensembl.transcript}:${data.coords.hgvs}`)
      : encodeURIComponent(data.coords.gene);

  return {
    franklin: `https://franklin.genoox.com/clinical-db/variant/snp/${fID}`,
    clinvar: `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(cvTerm)}`,
    gnomadv4: `https://gnomad.broadinstitute.org/variant/${gV4}?dataset=gnomad_r4`,
    spliceai: `https://spliceailookup.broadinstitute.org/#variant=${broad}&hg=38&bc=basic&distance=500&mask=0&ra=0`,
    clingen: `https://search.clinicalgenome.org/kb/genes?page=1&size=25&search=${data.coords.gene}`,
    omim: `https://www.omim.org/search?search=${data.coords.gene}`,
    gr: `https://www.ncbi.nlm.nih.gov/books/NBK1116/?term=${data.coords.gene}`,
    scholar: `https://scholar.google.com/scholar?q=${encodeURIComponent(schQ)}`,
    decipher: `https://www.deciphergenomics.org/search?q=${data.coords.gene}`,
    hgmd: `https://my.qiagendigitalinsights.com/bbp/view/hgmd/pro/gene.php?gene=${data.coords.gene}`,
    gtex: `https://gtexportal.org/home/gene/${data.coords.gene}`,
    gnomadv2: `https://gnomad.broadinstitute.org/variant/${gV2str}?dataset=gnomad_r2_1`,
    mastermind: `https://mastermind.genomenon.com/articles?mutation=${mmHgvs}`,
    liftover: `https://liftover.broadinstitute.org/#input=${encodeURIComponent(data.coords.hg38String)}&hg=hg38-to-hg19`,
  };
}

function openTab(url) { window.open(url, '_blank'); }
function launchOne(k) { const urls = getUrls(); if (urls[k]) openTab(urls[k]); }
function openSeq(keys) {
  let i = 0; const urls = getUrls(); function next() {
    if (i >= keys.length) return; if (urls[keys[i]]) openTab(urls[keys[i]]);
    i++; setTimeout(next, 220);
  } next();
}
function launchAll() { const order = ['clingen', 'omim', 'franklin', 'gnomadv4', 'gnomadv2', 'clinvar', 'spliceai', 'decipher', 'hgmd', 'mastermind', 'gtex', 'gr', 'scholar', 'liftover']; openSeq([...order].reverse()); }
function launchSelected() { if (!selectionOrder.length) { alert('Check boxes first.'); return; } openSeq([...selectionOrder].reverse()); }

/**
 * promptZygosity() — show the zygosity modal; resolves to the chosen value (or null if cancelled).
 * Shared by copyReport() and saveReportToFolder().
 */
function promptZygosity() {
  return new Promise(resolve => {
    const modal = document.getElementById('zygosityModal');
    modal.style.display = 'flex';
    window.resolveZygosity = (val) => {
      modal.style.display = 'none';
      data.selectedZygosity = val || 'Unknown';
      resolve(val);
    };
  });
}

/**
 * renderLockedSectionsImage(lockedHtml)
 * Inputs:  lockedHtml (string) — an inline-styled HTML fragment (ClinVar+gnomAD, or OMIM).
 * Output:  returns a PNG data-URI rasterised from that fragment via html2canvas, with a
 *          baked-in "retrieved <date time>" footer so the snapshot is tamper-evident.
 * Effects: briefly mounts an off-screen container under document.body, then removes it (finally).
 * Failures: throws if html2canvas is unavailable or rendering fails — caller falls back to inline HTML.
 */
async function renderLockedSectionsImage(lockedHtml) {
  if (typeof html2canvas !== 'function') throw new Error('html2canvas not loaded');
  const host = document.createElement('div');
  host.style.cssText = 'position:absolute; left:-10000px; top:0; width:720px; ' +
    'background:#ffffff; padding:12px; font-family:Calibri, sans-serif; font-size:11pt; color:#000;';
  const stamp = new Date().toLocaleString();
  host.innerHTML = lockedHtml +
    `<div style="margin-top:8px; font-size:8pt; color:#666;">Locked snapshot — retrieved ${stamp}</div>`;
  document.body.appendChild(host);
  try {
    const canvas = await html2canvas(host, {
      scale: 2, backgroundColor: '#ffffff', windowWidth: 760, logging: false,
    });
    return canvas.toDataURL('image/png');
  } finally {
    host.remove();
  }
}

/**
 * buildGnomadDetailHtml(det, hg38, linkHtml)
 * Inputs:  det      — data.gnomad.detailed (read-only; passed explicitly per CLAUDE.md Rule 2).
 *          hg38     — display coordinate string (e.g. "chr17-7669642-G-A").
 *          linkHtml — pre-built "View on gnomAD ↗" anchor placed next to the coordinate ('' if none).
 * Output:  inline-styled HTML — a borderless summary table (Exomes/Genomes/Total) plus a Genetic
 *          Ancestry Group Frequencies table — styled to match the gnomAD v4 variant page itself
 *          (clean light rules, compact font, Pass / No-variant badges) rather than the heavier
 *          ClinVar/ClinGen grid. Pure function: no DOM reads, no state mutation.
 * Used by: buildReportHtml() when the variant is PRESENT in gnomAD (overall total AC > 0).
 */
function buildGnomadDetailHtml(det, hg38, linkHtml = '') {
  // gnomAD-native styling: compact font, no cell grid — just a header rule + subtle row rules.
  const FS   = 'font-size:8.5pt;';
  const TBL  = `border-collapse:collapse; ${FS} margin:5px 0 4px 0;`;
  const TBLW = `${TBL} width:100%;`;   // full-width (matches the ClinVar/ClinGen report tables)
  const SLBL = 'padding:2px 12px 2px 0; text-align:right; white-space:nowrap;';
  const SVAL = 'padding:2px 16px; text-align:right; white-space:nowrap;';
  const SHV  = 'padding:2px 16px; text-align:right; font-weight:bold; border-bottom:1px solid #cfcfcf;';
  const SHB  = 'border-bottom:1px solid #cfcfcf;';
  const AHL  = 'padding:3px 10px 3px 6px; text-align:left;  font-weight:bold; background:#f2f2f2; border-bottom:1px solid #cfcfcf;';
  const AHR  = 'padding:3px 14px; text-align:right; font-weight:bold; background:#f2f2f2; border-bottom:1px solid #cfcfcf;';
  const ADL  = 'padding:2px 10px 2px 6px; text-align:left;  border-bottom:1px solid #ededed; white-space:nowrap;';
  const ADR  = 'padding:2px 14px; text-align:right; border-bottom:1px solid #ededed; white-space:nowrap;';

  const num   = n => (n == null ? '—' : String(n));
  // Summary allele frequency: gnomAD shows a plain "0" when AC=0 (vs "0.000" in the ancestry table).
  const sAF   = (ac, an) => (an > 0) ? (ac === 0 ? '0' : (ac / an).toPrecision(4)) : '—';
  const sFAF  = (popmax, an) => (popmax != null) ? (popmax === 0 ? '0' : Number(popmax).toPrecision(4)) : (an > 0 ? '0' : '—');
  const aAF   = (ac, an) => (an > 0) ? (ac === 0 ? '0.000' : (ac / an).toPrecision(4)) : '—';
  const badge = (txt, bg, fg) => `<span style="background:${bg}; color:${fg}; padding:1px 6px; border-radius:3px; font-weight:bold;">${txt}</span>`;
  const filtCell = a => !a.present
    ? badge('No variant', '#f4d4d4', '#a11414')
    : ((!a.filters || a.filters.length === 0) ? badge('Pass', '#cfe8d6', '#1a7431') : badge(a.filters.join(', '), '#fbe6c5', '#8a5a00'));

  // Aggregate one source. `present` = the variant was actually called in this dataset (processSrc ran
  // and set filters to an array); when false the source object is just the reset placeholder.
  const agg = s => {
    const present = !!s && Array.isArray(s.filters);
    const o = (s && s.populations && s.populations.overall && s.populations.overall.all) || {};
    return {
      present,
      ac: present ? (s.ac ?? o.ac ?? 0) : 0,
      an: present ? (s.an ?? o.an ?? 0) : 0,
      hom: present ? (s.hom ?? o.hom ?? 0) : 0,
      popmax: present ? (s.popmax ?? null) : null,
      filters: present ? s.filters : null
    };
  };
  const ex = agg(det.exome), ge = agg(det.genome), to = agg(det.total);
  // gnomAD reports a source's coverage AN even when the variant isn't called there ("No variant").
  // Derive the absent source's AC/AN/hom as total − the present source so the column matches the site.
  const derive = other => ({ ac: Math.max(0, (to.ac || 0) - (other.ac || 0)), an: Math.max(0, (to.an || 0) - (other.an || 0)), hom: Math.max(0, (to.hom || 0) - (other.hom || 0)) });
  const exN = ex.present ? ex : (ge.present ? derive(ge) : ex);
  const geN = ge.present ? ge : (ex.present ? derive(ex) : ge);

  const sRow = (label, ev, gv, tv) =>
    `<tr><td style="${SLBL}">${label}</td><td style="${SVAL}">${ev}</td><td style="${SVAL}">${gv}</td><td style="${SVAL}">${tv}</td></tr>`;
  const summary =
    `<table style="${TBL}"><tbody>` +
    `<tr><td style="${SHB}"></td><td style="${SHV}">Exomes</td><td style="${SHV}">Genomes</td><td style="${SHV}">Total</td></tr>` +
    sRow('Filters', filtCell(ex), filtCell(ge), '') +
    sRow('Allele Count', num(exN.ac), num(geN.ac), num(to.ac)) +
    sRow('Allele Number', num(exN.an), num(geN.an), num(to.an)) +
    sRow('Allele Frequency', sAF(exN.ac, exN.an), sAF(geN.ac, geN.an), sAF(to.ac, to.an)) +
    sRow('Grpmax Filtering AF (95% confidence)', sFAF(ex.popmax, exN.an), sFAF(ge.popmax, geN.an), sFAF(to.popmax, to.an)) +
    sRow('Number of homozygotes', num(exN.hom), num(geN.hom), num(to.hom)) +
    `</tbody></table>`;

  // Genetic Ancestry Group Frequencies — combined (total) source, sorted by AF desc (gnomAD default).
  const LABELS = {
    nfe: 'European (non-Finnish)', amr: 'Admixed American', remaining: 'Remaining',
    afr: 'African/African American', asj: 'Ashkenazi Jewish', eas: 'East Asian',
    fin: 'European (Finnish)', mid: 'Middle Eastern', ami: 'Amish', sas: 'South Asian'
  };
  const pops = det.total?.populations || {};
  const groups = Object.keys(pops).filter(k => k && k !== 'overall').map(k => {
    const a = pops[k].all || { ac: 0, an: 0, hom: 0 };
    return { label: LABELS[k] || k.toUpperCase(), ac: a.ac || 0, an: a.an || 0, hom: a.hom || 0 };
  }).sort((a, b) => (b.an ? b.ac / b.an : 0) - (a.an ? a.ac / a.an : 0));

  const aRow = (label, ac, an, hom, bold) =>
    `<tr${bold ? ' style="font-weight:bold;"' : ''}><td style="${ADL}">${label}</td><td style="${ADR}">${num(ac)}</td>` +
    `<td style="${ADR}">${num(an)}</td><td style="${ADR}">${num(hom)}</td><td style="${ADR}">${aAF(ac, an)}</td></tr>`;

  const ov = pops.overall || {};
  const sexRow = (label, s) => (s && s.an > 0) ? aRow(label, s.ac, s.an, s.hom, false) : '';
  const ancestry =
    `<table style="${TBLW}"><tbody>` +
    `<tr><td style="${AHL}">Genetic Ancestry Group</td><td style="${AHR}">Allele Count</td>` +
    `<td style="${AHR}">Allele Number</td><td style="${AHR}">Number of Homozygotes</td><td style="${AHR}">Allele Frequency</td></tr>` +
    groups.map(r => aRow(r.label, r.ac, r.an, r.hom, false)).join('') +
    sexRow('XX', ov.XX) +
    sexRow('XY', ov.XY) +
    aRow('Total', to.ac, to.an, to.hom, true) +
    `</tbody></table>`;

  return `<b>gnomAD v4.1.1:</b> ${hg38 || ''}${linkHtml || ''}` +
    summary +
    `<p style="margin:7px 0 1px 0; ${FS}"><b>Genetic Ancestry Group Frequencies</b></p>` +
    ancestry;
}

/**
 * buildReportHtml(zygosity, opts)
 * Inputs:  zygosity (string) + current data.* state (reads only — does not mutate primary coords).
 *          opts.lockSectionsAsImage (bool) — when true, the externally-fetched ClinVar+gnomAD and
 *          OMIM sections are rasterised to flat, timestamped PNGs (pastes into Word as uneditable
 *          pictures for traceability); everything else stays editable text. Default false.
 * Output:  returns the Word-pasteable clinical report HTML string.
 * Effects: lazy-fetches ClinGen validity / GeneReviews / ClinVar detail if not yet loaded (non-blocking).
 * Used by: copyReport() (→ clipboard, image-locked) and saveReportToFolder() (→ disk, plain HTML).
 */
async function buildReportHtml(zygosity, { lockSectionsAsImage = false } = {}) {
  const g = id => (document.getElementById(id)?.innerText || '-').replace(/[✅✔☑️]/g, '').trim();

  // ClinVar: "★☆☆☆ (1/4)" format matching report template
  const starCount = data.clinvar.stars ?? 0;
  const starStr = starCount > 0
    ? ('★'.repeat(starCount) + '☆'.repeat(4 - starCount) + ` (${starCount}/4)`)
    : 'Not rated';

  // gnomAD: "Absent (hg38 coords)" when AC=0; detailed breakdown when present.
  // The "View on gnomAD" link sits next to the coordinate in the heading (not after the tables).
  const _gnomadUrl = (typeof getUrls === 'function') ? getUrls().gnomadv4 : null;
  const gnomadLink = (_gnomadUrl && !/undefined/.test(_gnomadUrl))
    ? ` <a href="${_gnomadUrl}" style="font-size:9pt; color:#0563C1; text-decoration:underline;">View on gnomAD ↗</a>`
    : '';
  let gnomadFreqLine = `<b>gnomAD v4.1.1:</b> —${gnomadLink}`;
  const gnomadDet = data.gnomad.detailed;
  if (gnomadDet?.total?.populations) {
    const all = gnomadDet.total.populations.overall?.all || { ac: 0, an: 0, hom: 0 };
    const hg38 = data.coords.hg38String || g('dHg38');
    if (all.ac === 0) {
      gnomadFreqLine = `<b>gnomAD v4.1.1:</b> ${hg38} — Absent${gnomadLink}`;
    } else {
      // Variant present → paste the full gnomAD breakdown (summary + ancestry-group table).
      gnomadFreqLine = buildGnomadDetailHtml(gnomadDet, hg38, gnomadLink);
    }
  }

  let proteinStr = g('dProtein').split(':').pop().trim();
  if (proteinStr.startsWith('p.') && !proteinStr.includes('(')) {
    proteinStr = 'p.(' + proteinStr.substring(2) + ')';
  }

  const revel = (data.scores.revel !== null && !isNaN(data.scores.revel)) ? parseFloat(data.scores.revel).toFixed(3) : '-';
  const spliceAI = (data.scores.spliceAI !== null && !isNaN(data.scores.spliceAI)) ? parseFloat(data.scores.spliceAI).toFixed(3) : '-';

  // Ensure the externally-fetched sections are loaded before building (lazy — covers a copy issued
  // before the background fetches finished). These groups touch independent namespaces
  // (data.geneValidity / data.geneReviews / data.clinvar / OMIM local cache), so they run
  // concurrently rather than serially (CLAUDE.md Rule 4 — async independence): on a cold copy this
  // turns the wait from the SUM of the network round-trips into the MAX. The validity→evidence and
  // geneReviews→focus orderings are preserved within each group, and every call keeps its own
  // try/catch so one module's failure can't cascade into the others (Rule 3 — error containment).
  await Promise.all([
    // Group A: ClinGen Gene-Disease Validity → GCEP evidence-summary narratives (heavy HTML pages).
    (async () => {
      if (data.coords.gene && (!data.geneValidity?.fetched || data.geneValidity.gene !== data.coords.gene)) {
        try { await fetchGeneValidity(data.coords.gene); } catch { /* non-blocking */ }
      }
      try { await enrichGeneValidityEvidence(data.geneValidity?.curations); } catch { /* non-blocking */ }
      if (data.geneValidity?.fetched && data.ptPhenoTexts?.length && data.omimSynopsisResults) {
        try { await focusGeneValidityByPhenotype(data.omimSynopsisResults); } catch { /* non-blocking */ }
      }
    })(),
    // Group B: GeneReviews disease context → phenotype focus (focus is sync; depends on the fetch).
    (async () => {
      if (data.coords.gene && (!data.geneReviews?.fetched || data.geneReviews.gene !== data.coords.gene)) {
        try { await fetchGeneReviews(data.coords.gene); } catch { /* non-blocking */ }
      }
      if (data.geneReviews?.fetched && data.ptPhenoTexts?.length && data.omimSynopsisResults) {
        try { focusGeneReviewsByPhenotype(data.omimSynopsisResults); } catch { /* non-blocking */ }
      }
    })(),
    // Group C: ClinVar VCV full-record detail (independent).
    (async () => {
      try { await ensureClinVarDetail(); } catch { /* non-blocking */ }
    })(),
    // Group D: OMIM local cache (large genemap2 file, parsed async at startup). getOmimGeneTable()
    // returns [] until it finishes, so a copy issued mid-load would otherwise bake an empty
    // "No OMIM phenotype entry found" into the locked image. omimLocalReady() resolves immediately
    // once loaded (and resolves even on parse failure), so this adds no latency in the common case.
    (async () => {
      try { if (typeof omimLocalReady === 'function') await omimLocalReady(); } catch { /* non-blocking */ }
    })(),
  ]);

  // Curation-support concordance for the report (depends on GeneReviews + variant + phenotype above).
  try { if (typeof refreshCurationSupport === 'function') refreshCurationSupport(); } catch { /* non-blocking */ }

  // Age of onset for definitive curations — OMIM synopsis primary, GeneReviews fallback (both ready above).
  try {
    if (typeof enrichGeneValidityOnset === 'function')
      await enrichGeneValidityOnset(data.geneValidity?.curations, data.geneValidity?.gene, data.geneReviews?.chapters);
  } catch { /* non-blocking */ }

  const esc = s => String(s ?? '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

  const stamp = new Date().toLocaleString();
  const tsLine = `<span style="font-size:8pt; color:#666; font-style:italic;">Retrieved: ${stamp}</span><br>`;

  // ── Table builders (inline styles → render as real tables when pasted into Word) ──
  const TBL = 'border-collapse:collapse; font-size:10pt; margin:4px 0 10px 0; width:100%;';
  const TH  = 'border:1px solid #888; padding:3px 7px; background:#f0f0f0; text-align:left; font-weight:bold;';
  const TD  = 'border:1px solid #888; padding:3px 7px; vertical-align:top;';

  // ClinGen Gene-Disease Validity table
  const cur = data.geneValidity?.curations || [];
  let geneValidityHtml;
  if (cur.length) {
    // Each curation is one main row + an optional full-width evidence-summary sub-row (like ClinVar submissions).
    const rows = cur.map(c => {
      const mainRow =
        `<tr>` +
        `<td style="${TD}font-weight:bold;">${esc(c.classification)}</td>` +
        `<td style="${TD}">${esc(c.disease)}</td>` +
        `<td style="${TD}font-size:9pt;">${esc(typeof expandClingenMoi === 'function' ? expandClingenMoi(c.moi) : c.moi)}</td>` +
        `<td style="${TD}font-size:9pt;">${c.onset ? esc(c.onset) + (c.onsetSource ? ` <span style="font-size:8.5pt;color:#666;">(${esc(c.onsetSource)})</span>` : '') : '—'}</td>` +
        `<td style="${TD}font-size:9pt;">${esc(c.gcep)}</td>` +
        `<td style="${TD}font-size:9pt;">${esc(c.date)}</td>` +
        `</tr>`;
      if (!c.evidenceSummary) return mainRow;
      const m = typeof extractDiseaseMechanism === 'function' ? extractDiseaseMechanism(c.evidenceSummary) : null;
      const mechCount = m ? ['LoF', 'GoF', 'DN'].filter(k => m[k]).length : 0;
      const multiNote = mechCount >= 2
        ? `<br><i style="color:#555;">Multiple disease mechanisms — refer to ClinGen for full details.</i>`
        : '';
      const evidenceRow =
        `<tr><td style="${TD}border-top:none; font-size:9.5pt; background:#f9f9f9;" colspan="6">` +
        `<b>Evidence Summary:</b> ${esc(c.evidenceSummary)}${multiNote}</td></tr>`;
      return mainRow + evidenceRow;
    }).join('');
    geneValidityHtml =
      `<table style="${TBL}"><thead><tr>` +
      `<th style="${TH}">Classification</th><th style="${TH}">Disease</th>` +
      `<th style="${TH}">MOI</th><th style="${TH}">Age of onset</th>` +
      `<th style="${TH}">Expert Panel (GCEP)</th>` +
      `<th style="${TH}">Date</th></tr></thead><tbody>${rows}</tbody></table>`;
  } else {
    geneValidityHtml = `<i>No ClinGen Gene-Disease Validity curation found.</i><br>`;
  }

  // Offline flag for the pasted report: when the Evidence Summary text came from the Go helper's
  // on-disk copy (couldn't reach Google to refresh), note it inline so the reader knows it may not
  // reflect the latest ClinGen curation. Plain text only — survives the clipboard/Word paste.
  const gvProv = typeof clingenEvidenceProvenance === 'function' ? clingenEvidenceProvenance() : null;
  const gvOfflineNote = (gvProv && gvProv.source === 'cache')
    ? ` <i style="color:#b26a00;">(evidence: offline copy${gvProv.builtAt ? `, as of ${esc(gvProv.builtAt)}` : ''})</i>`
    : '';

  // ── GeneReviews — disease mechanism / penetrance / genotype-phenotype / clinical description ──
  // Rendered in the copied (Word) report alongside ClinGen validity. Because the report is a static
  // paste, the dense Clinical Description is truncated with a [full chapter] link (not an interactive
  // expander). One block per disease chapter (a gene can map to several, e.g. BRCA2 → HBOC + Fanconi).
  const grAll = data.geneReviews?.chapters || [];
  const grFocusNbk = data.geneReviews?.focusNbk || null;
  const grFocused = grFocusNbk ? grAll.find(c => c.nbk === grFocusNbk) : null;
  const grChapters = grFocused ? [grFocused] : grAll;   // phenotype focus → only the best-matching disease
  const grOthers = grFocused ? grAll.filter(c => c.nbk !== grFocusNbk) : [];
  let geneReviewsHtml;
  if (grChapters.length) {
    const grField = (label, val, n, url) => {
      const full = String(val || '').trim();
      if (!full) return '';
      const cut = full.length > n;
      const shown = cut ? full.slice(0, n).replace(/\s+\S*$/, '') : full;
      return `<p style="margin:2px 0 6px 0;"><b>${label}:</b> ${esc(shown)}` +
        (cut ? ` … <a href="${esc(url)}">[full chapter]</a>` : '') + `</p>`;
    };
    // Frequency-of-features table → a real Word table; features matching the patient phenotype are bolded + flagged.
    const ptTerms = data.ptPhenoTexts || [];
    const grFreqTable = (freqs) => {
      if (!Array.isArray(freqs) || !freqs.length) return '';
      const matched = matchPhenotypeFrequencies({ phenotypeFrequencies: freqs, ptPhenoTexts: ptTerms });
      const hitFeatures = new Set(matched.filter(m => m.matched).map(m => m.feature));
      const rows = freqs.slice(0, 30).map(f => {
        const hit = hitFeatures.has(f.feature);
        const cell = hit ? `${TD}background:#fff2cc;` : TD;
        const feat = hit ? `<b>${esc(f.feature)}</b> ✓` : esc(f.feature);
        return `<tr><td style="${cell}">${feat}</td><td style="${cell}">${esc(f.freq)}</td><td style="${cell}">${esc(f.comment || '')}</td></tr>`;
      }).join('');
      return `<p style="margin:4px 0 2px 0;"><b>Frequency of Select Features</b>${ptTerms.length ? ` <span style="font-size:9pt; color:#666;">(✓ = matches patient phenotype)</span>` : ''}:</p>` +
        `<table style="${TBL}"><thead><tr><th style="${TH}">Feature</th><th style="${TH}">Frequency</th><th style="${TH}">Comment</th></tr></thead><tbody>${rows}</tbody></table>`;
    };
    geneReviewsHtml = grChapters.map(c => {
      const head = `<p style="margin:6px 0 2px 0;"><b><a href="${esc(c.url)}">${esc(c.title || c.ch_id)}</a></b>` +
        (c.nbk ? ` (${esc(c.nbk)})` : '') + `:</p>`;
      const body =
        grField('Mechanism of disease', c.mechanism, 600, c.url) +
        grField('Mode of Inheritance', c.inheritance, 800, c.url) +
        grField('Penetrance', c.penetrance, 1200, c.url) +
        grField('Genotype-Phenotype', c.genotypePhenotype, 1200, c.url) +
        grField('Clinical Description', c.clinicalDescription, 700, c.url) +
        grField('Management', c.management, 1000, c.url) +
        grFreqTable(c.phenotypeFrequencies);
      return head + (body || `<p style="margin:2px 0 6px 0;"><i>See chapter.</i></p>`);
    }).join('') +
      (grOthers.length ? `<p style="margin:4px 0 8px 0; font-size:9pt; color:#666;">Other GeneReviews chapters for this gene: ` +
        grOthers.map(c => `<a href="${esc(c.url)}">${esc(c.title || c.nbk)}</a>`).join(', ') + `</p>` : '') +
      `<p style="margin:4px 0 8px 0; font-size:9pt; color:#666;"><i>Source: GeneReviews®, University of Washington — noncommercial research / clinical use; see chapter links for full text.</i></p>`;
  } else if (data.geneReviews?.unavailable) {
    geneReviewsHtml = `<i>GeneReviews unavailable — local helper not reachable.</i><br>`;
  } else {
    geneReviewsHtml = `<i>No GeneReviews chapter found.</i><br>`;
  }

  // ── ClinVar full record detail (from VCV XML → data.clinvar.detail) ──
  // Layout goals: one compact "Record overview" table (variant identity + timeline merged),
  // a Conditions table, and — critically — submissions rendered as ONE scannable table
  // (submitter per row) instead of N stacked key/value tables. Long submitter comments drop
  // into a full-width sub-row beneath each entry so the comparison columns stay aligned.
  const cd = data.clinvar.detail;
  const CAP = 'margin:8px 0 2px 0; font-size:10pt;';   // consistent sub-table caption
  const SUB = 'font-size:8.5pt; color:#666;';          // secondary inline detail (SCV, dates)
  let clinvarDetailHtml = '';
  if (cd) {
    const fmtLoc = (l) => l && l.pos ? `${l.chr}:${l.pos} ${l.ref || ''}&gt;${l.alt || ''}`.trim() : '—';
    const lenStr = cd.length ? `${cd.length} bp` : '';

    // Record overview — Accession + Variation ID side by side in a single row (Accession first).
    clinvarDetailHtml +=
      `<p style="${CAP}"><b>Record overview</b></p>` +
      `<table style="${TBL}"><tbody>` +
      `<tr>` +
      `<td style="${TD}width:18%;"><b>Accession</b></td>` +
      `<td style="${TD}width:32%;">${esc(cd.accession ? cd.accession + '.' + cd.version : data.clinvar.accession || '—')}</td>` +
      `<td style="${TD}width:18%;"><b>Variation ID</b></td>` +
      `<td style="${TD}">${esc(cd.variationId || '—')}</td>` +
      `</tr>` +
      `</tbody></table>`;

    // Conditions (RCV) — one row per associated condition.
    if (cd.conditions?.length) {
      const rows = cd.conditions.map(c =>
        `<tr><td style="${TD}">${esc(c.names.join('; '))}</td>` +
        `<td style="${TD}">${esc(c.classification)}${c.submissionCount ? ` <span style="${SUB}">(${esc(c.submissionCount)} sub)</span>` : ''}</td>` +
        `<td style="${TD}font-size:9pt;">${esc(c.reviewStatus || '—')}</td>` +
        `<td style="${TD}">${esc(c.lastEvaluated || '—')}</td>` +
        `<td style="${TD}font-size:9pt;">${esc(c.rcv)}</td></tr>`
      ).join('');
      clinvarDetailHtml +=
        `<p style="${CAP}"><b>Conditions (${cd.conditions.length})</b></p>` +
        `<table style="${TBL}"><thead><tr>` +
        `<th style="${TH}">Condition</th><th style="${TH}">Classification</th>` +
        `<th style="${TH}">Review status</th><th style="${TH}">Last evaluated</th>` +
        `<th style="${TH}">RCV record</th></tr></thead><tbody>${rows}</tbody></table>`;
    }

    // Submitted interpretations (SCV) — ONE table, one submitter per row. Observation metadata
    // (method · origin · affected · assertion method) and the free-text comment fold into a
    // full-width sub-row beneath each submitter so the four comparison columns stay aligned.
    if (cd.submissions?.length) {
      const aggConditions = cd.conditions?.[0]?.names?.join('; ') || '';
      const rows = cd.submissions.map(s => {
        const cond = (s.conditions && s.conditions.length) ? s.conditions.join('; ') : aggConditions;
        const obs = [
          s.method,
          s.origin && `origin: ${s.origin}`,
          s.affectedStatus && `affected: ${s.affectedStatus}`,
          s.assertionMethod && `criteria: ${s.assertionMethod}`
        ].filter(Boolean).join(' · ');
        // Visual priority: Classification and Condition are the focus; Submitter and Review status are secondary.
        const mainRow =
          `<tr>` +
          `<td style="${TD}font-weight:bold;">${esc(s.classification || '—')}${s.lastEvaluated ? `<br><span style="${SUB}">${esc(s.lastEvaluated)}</span>` : ''}</td>` +
          `<td style="${TD}">${esc(cond || '—')}</td>` +
          `<td style="${TD}font-size:9pt;color:#555;">${esc(s.submitter || '—')}${s.scv ? `<br><span style="${SUB}">${esc(s.scv)}</span>` : ''}</td>` +
          `<td style="${TD}font-size:9pt;color:#555;">${esc(s.reviewStatus || '—')}</td>` +
          `</tr>`;
        const detail = [];
        if (obs) detail.push(`<i style="color:#666;">${esc(obs)}</i>`);
        if (s.comment) detail.push(`<span style="color:#1a1a1a;">${esc(s.comment)}</span>`);
        const detailRow = detail.length
          ? `<tr><td style="${TD}border-top:none; font-size:10pt; background:#f9f9f9;" colspan="4">${detail.join('<br>')}</td></tr>`
          : '';
        return mainRow + detailRow;
      }).join('');
      clinvarDetailHtml +=
        `<p style="${CAP}"><b>Submitted interpretations (${cd.submissions.length})</b></p>` +
        `<table style="${TBL}"><thead><tr>` +
        `<th style="${TH}">Classification</th><th style="${TH}">Condition</th>` +
        `<th style="${TH}">Submitter</th><th style="${TH}">Review status</th>` +
        `</tr></thead><tbody>${rows}</tbody></table>`;
    }
  }

  // ClinGen Dosage Sensitivity — decode HI/TS scores to descriptive labels.
  // UCSC's position-based track misses some regions (e.g. PAR1 on chrX/Y);
  // fall back to the local clingenDosageDict when the live value is absent.
  const dosageLabel = (score) => {
    const s = String(score ?? '').trim();
    const map = {
      '0': '0 (No Evidence)', '1': '1 (Little Evidence)',
      '2': '2 (Some Evidence)', '3': '3 (Sufficient Evidence)',
      '30': '30 (Gene Associated with Autosomal Recessive Phenotype)',
      '40': '40 (Dosage Sensitivity Unlikely)'
    };
    return map[s] || (s || 'Not Curated');
  };
  const isMissingDosage = (v) => !v || v === 'Not Curated' || v === '—';
  let reportHaplo = data.dosage.haplo;
  let reportTriplo = data.dosage.triplo;
  if ((isMissingDosage(reportHaplo) || isMissingDosage(reportTriplo)) && data.coords.gene) {
    const localDosage = (typeof getClinGenDosage === 'function') ? getClinGenDosage(data.coords.gene) : null;
    if (localDosage && !isMissingDosage(localDosage.hi)) reportHaplo  = localDosage.hi;
    if (localDosage && !isMissingDosage(localDosage.ts)) reportTriplo = localDosage.ts;
  }
  const dosageEval = cd?.dosage?.haplo?.lastEvaluated || cd?.dosage?.triplo?.lastEvaluated || null;
  const dosageHtml =
    `<table style="${TBL}"><tbody>` +
    `<tr><td style="${TD}width:50%;"><b>Haploinsufficiency</b></td><td style="${TD}">${esc(dosageLabel(reportHaplo))}</td></tr>` +
    `<tr><td style="${TD}"><b>Triplosensitivity</b></td><td style="${TD}">${esc(dosageLabel(reportTriplo))}</td></tr>` +
    (dosageEval ? `<tr><td style="${TD}"><b>Last evaluated</b></td><td style="${TD}">${esc(dosageEval)}</td></tr>` : '') +
    `</tbody></table>`;

  // OMIM phenotype table (local genemap2 data)
  const omimRows = (typeof getOmimGeneTable === 'function') ? getOmimGeneTable(data.coords.gene) : [];
  let omimHtml;
  if (omimRows.length) {
    const rows = omimRows.map(r =>
      `<tr><td style="${TD}">${esc(r.location)}</td>` +
      `<td style="${TD}">${esc(r.phenotype)}</td>` +
      `<td style="${TD}">${esc(r.mimNumber)}</td>` +
      `<td style="${TD}">${esc(r.inheritance)}</td>` +
      `<td style="${TD}">${esc(r.mappingKey)}</td></tr>`
    ).join('');
    omimHtml =
      `<table style="${TBL}"><thead><tr>` +
      `<th style="${TH}">Location</th><th style="${TH}">Phenotype</th>` +
      `<th style="${TH}">Phenotype MIM number</th><th style="${TH}">Inheritance</th>` +
      `<th style="${TH}">Phenotype mapping key</th></tr></thead><tbody>${rows}</tbody></table>`;
  } else {
    omimHtml = `<i>No OMIM phenotype entry found for ${esc(data.coords.gene || '—')}.</i><br>`;
  }

  // ── Source links — defined early so they can be inlined into section headings ──
  // URLs reuse the canonical patterns from getUrls() (single source of truth).
  // ClinVar upgrades to a direct variation page when the VariationID is known.
  const reportUrls = (typeof getUrls === 'function') ? getUrls() : {};
  const _validUrl = u => u && !/undefined/.test(u);
  const clinvarSrcUrl = data.clinvar.detail?.variationId
    ? `https://www.ncbi.nlm.nih.gov/clinvar/variation/${data.clinvar.detail.variationId}/`
    : reportUrls.clinvar;
  const SRCLNK = 'font-size:9pt; color:#0563C1; text-decoration:underline;';
  const srcLink = (url, label) => _validUrl(url) ? `<a href="${esc(url)}" style="${SRCLNK}">${esc(label)} ↗</a>` : '';
  const clinvarSrc = srcLink(clinvarSrcUrl, 'View on ClinVar');
  const gnomadSrc  = srcLink(reportUrls.gnomadv4, 'View on gnomAD');
  const clingenSrc = srcLink(reportUrls.clingen, 'View on ClinGen');
  const omimSrc    = srcLink(reportUrls.omim, 'View on OMIM');

  // ── Traceability lock for externally-fetched evidence ──
  // ClinVar + gnomAD (contiguous) and OMIM (separate, below the curation tables) are the three
  // externally-fetched sections. When lockSectionsAsImage is set (Copy To Analysis Summary), each
  // is rasterised to a flat, timestamped PNG so it pastes into Word as an uneditable picture.
  // Report order is preserved (OMIM stays below ClinGen/GeneReviews), so two images are produced
  // rather than one. The curation sections and the analyst's own fields remain editable text.
  // Headline: aggregate germline classification (prefer parsed detail, then the live UI element,
  // then the data-layer value) + star review status, with a compact condition/submission tally.
  const cvSigText = g('cvSig');
  const cvAgg = data.clinvar.detail?.classification
    || (cvSigText && cvSigText !== '-' && cvSigText !== '—' ? cvSigText : null)
    || data.clinvar.sig
    || '—';
  const cvTally = data.clinvar.detail
    ? ` <span style="font-size:9pt; color:#555;">· ${data.clinvar.detail.conditions?.length || 0} condition(s), ${data.clinvar.detail.submissions?.length || data.clinvar.detail.numSubmissions || 0} submission(s)</span>`
    : '';

  // Clean versions (no hyperlinks) — used only for image rasterisation; links inside a PNG are dead.
  const clinvarBlockHtmlClean =
    `<b>ClinVar — ${esc(cvAgg)}</b> ${starStr}${cvTally}<br>` +
    `<span style="font-size:9pt; color:#555;">Accession: ${esc(data.clinvar.accession || '-')}</span>` +
    clinvarDetailHtml;
  const omimBlockHtmlClean = `<b>OMIM:</b><br>${omimHtml}`;

  // Link-inclusive versions — links inline on the headline / label line (editable HTML output).
  const clinvarBlockHtml =
    `<b>ClinVar — ${esc(cvAgg)}</b> ${starStr}${cvTally}${clinvarSrc ? ` ${clinvarSrc}` : ''}<br>` +
    `<span style="font-size:9pt; color:#555;">Accession: ${esc(data.clinvar.accession || '-')}</span>` +
    clinvarDetailHtml;
  const omimBlockHtml = `<b>OMIM:</b>${omimSrc ? ` ${omimSrc}` : ''}<br>${omimHtml}`;

  // gnomadFreqLine already carries its own "View on gnomAD" link next to the coordinate.
  let clinvarGnomadOut = clinvarBlockHtml + `<br>` + gnomadFreqLine + `<br><br>`;
  let omimOut = omimBlockHtml;
  if (lockSectionsAsImage && typeof html2canvas === 'function') {
    try {
      const imgStyle = 'display:block; max-width:100%; margin:4px 0;';
      // Sequential (not Promise.all): deterministic rendering matters more than the ~20 ms saved,
      // and html2canvas does not guarantee concurrency safety. The cost is small and one-off per copy.
      const cvImg = await renderLockedSectionsImage(clinvarBlockHtmlClean + `<br>${gnomadFreqLine}`);
      const omImg = await renderLockedSectionsImage(omimBlockHtmlClean);
      // After each image, emit source links as clickable text (they cannot survive inside a PNG).
      const cvLinks = [clinvarSrc, gnomadSrc].filter(Boolean).join(' &nbsp;·&nbsp; ');
      clinvarGnomadOut = `<img src="${cvImg}" style="${imgStyle}" alt="ClinVar &amp; gnomAD (locked snapshot)"><br>` +
        (cvLinks ? `<div style="margin:2px 0 6px 0;">${cvLinks}</div>` : '');
      omimOut = `<img src="${omImg}" style="${imgStyle}" alt="OMIM (locked snapshot)">` +
        (omimSrc ? `<div style="margin:2px 0 6px 0;">${omimSrc}</div>` : '');
    } catch (e) {
      console.warn('Locked-section image render failed; falling back to editable HTML:', e);
      // keep the inline-HTML fallbacks already assigned above
    }
  }

  // ── Curation Concordance — variant↔mechanism + patient phenotype↔frequency (decision support) ──
  const cs = data.curationSupport || {};
  let curationSupportHtml = '';
  {
    const mech = cs.mechanism;
    const matched = (cs.phenoMatches || []).filter(m => m.matched);
    const parts = [];
    if (matched.length) {
      parts.push(`<p style="margin:2px 0 2px 0;"><b>Patient phenotype frequency (GeneReviews):</b></p><ul style="margin:2px 0 6px 18px;">` +
        matched.map(m => {
          const featTag = m.feature && m.feature.toLowerCase() !== m.ptTerm.toLowerCase() ? ` <span style="color:#666;">(${esc(m.feature)})</span>` : '';
          const dzTag = m.disease && m.disease !== cs.chapterTitle ? ` <span style="color:#666;">[${esc(m.disease)}]</span>` : '';
          return `<li>${esc(m.ptTerm)} — <b>${esc(m.freq || '—')}</b> of affected${featTag}${dzTag}</li>`;
        }).join('') +
        `</ul>`);
    }
    curationSupportHtml = parts.join('');
  }

  // ── Bioinformatics predictions — variant-type-aware in-silico line ──
  // Mirrors the GG report prompt's variant-type requirements (buildVariantBioinformaticsBlock):
  // REVEL for missense, SpliceAI for splice, an NMD note for nonsense/frameshift. Renders ONE
  // concise line to match the Analysis Summary's fill-in style (see reference PDF). Reads the
  // CANONICAL data.ensembl.vepConsequence (top-level data.vepConsequence is never assigned) plus
  // data.scores — no mutation. `revel` / `spliceAI` are the pre-formatted strings from above.
  const bioPredHtml = (() => {
    const vc = [data.ensembl?.vepConsequence, data.vepConsequence].find(v => typeof v === 'string') || '';
    const isMissense = vc.includes('missense_variant');
    const isSplice   = ['splice_donor', 'splice_acceptor', 'splice_region'].some(t => vc.includes(t));
    const isLoF      = ['stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost'].some(t => vc.includes(t));
    const isInframe  = vc.includes('inframe_insertion') || vc.includes('inframe_deletion');

    if (isMissense) return `<b>REVEL score:</b> ${revel}<br>`;
    if (isSplice)   return `<b>SpliceAI:</b> ${spliceAI}<br>`;
    if (isLoF) {
      const exonNum    = data.ensembl?.vepExon ?? data.vepExon;
      const totalExons = data.ensembl?.vepTotalExons ?? data.vepTotalExons;
      if (exonNum != null && totalExons != null) {
        const en = parseInt(String(exonNum).split('/')[0]);
        const te = parseInt(String(totalExons));
        if (en && te && en === te) return `Predicted NMD escape (variant in last exon ${en}/${te})<br>`;
        if (en && te)              return `Predicted NMD (PTC in exon ${en}/${te})<br>`;
      }
      return `Predicted nonsense-mediated decay (NMD)<br>`;
    }
    if (isInframe) return `<i>REVEL / SpliceAI not applicable (in-frame indel)</i><br>`;
    // Undetermined type → keep both scores so nothing is hidden (prior behaviour).
    if (!vc) return `<b>REVEL score:</b> ${revel}<br><b>SpliceAI:</b> ${spliceAI}<br>`;
    // Synonymous / other → SpliceAI only if elevated, else not applicable (mirrors the prompt fallback).
    if (data.scores.spliceAI != null && parseFloat(data.scores.spliceAI) >= 0.1)
      return `<b>SpliceAI:</b> ${spliceAI} <i>(elevated — possible cryptic splice effect)</i><br>`;
    return `<i>In-silico predictors not applicable for this variant type</i><br>`;
  })();

  // ── Source links to the live databases ────────────────────────────────────
  // One hyperlink per external source so the pasted report can be traced back to the real data.
  // ClinVar/gnomAD/OMIM blocks are rasterised to PNGs in locked mode (links inside an image aren't
  // clickable), so these are emitted as editable text OUTSIDE those blocks — they survive every copy
  // path. URLs reuse the canonical patterns from getUrls() (single source of truth); the ClinVar one
  const htmlContent = `<div style="font-family: Calibri, sans-serif; font-size: 11pt;">` +
    `<b>${zygosity} ${data.ensembl.transcript || '-'}(${data.coords.gene || '-'}):${data.coords.hgvs || '-'} ${proteinStr} ${g('eExon')}</b><br>` +
    `<b>Genomic coordinates:</b> [GRCh38] ${g('dHg38')} | [hg19] ${g('dHg19')}<br>` +
    `<b>DP2:</b> <br>` +
    `<b>IGV check:</b> <br>` +
    `<b>Internal database:</b> <br>` +
    clinvarGnomadOut + tsLine +
    `<b>ClinGen — Gene-Disease Validity:</b>${clingenSrc ? ` ${clingenSrc}` : ''}${gvOfflineNote}<br>${geneValidityHtml}` + tsLine +
    // GeneReviews section HIDDEN from the Analysis Summary per request (2026-06-08).
    // No GeneReviews code deleted — the fetch + geneReviewsHtml builder above stay fully intact
    // (and the live GeneReviews card is unaffected). Re-enable by uncommenting the next line:
    // `<b>GeneReviews — Disease Mechanism / Penetrance / Expressivity:</b><br>${geneReviewsHtml}` + tsLine +
    (curationSupportHtml ? `<b>Curation Concordance:</b><br>${curationSupportHtml}` : '') +
    `<b>ClinGen — Dosage Sensitivity:</b>${clingenSrc ? ` ${clingenSrc}` : ''}<br>${dosageHtml}` + tsLine +
    omimOut + tsLine +
    `<b>HGMD:</b> <br>` +
    `<br>` +
    `<b>Bioinformatics predictions:</b> <br>` +
    bioPredHtml +
    `<b>Literature search:</b> <br>` +
    `<b>Conclusion and classification:</b> </div>`;

  return htmlContent;
}

/**
 * copyReport() — prompt zygosity, build the report, write it to the clipboard as rich HTML.
 */
async function copyReport() {
  const zygosity = await promptZygosity();
  if (!zygosity) return;

  const btn = document.getElementById('btnCopyReport');
  const originalLabel = btn?.innerHTML;
  const resetLabel = () => { if (btn) btn.innerHTML = originalLabel || '<span class="btn-icon">📋</span> Copy To Analysis Summary'; };
  if (btn) btn.innerHTML = '<span class="btn-icon">⏳</span> Preparing…';

  // Bind the clipboard write to the user gesture *now* by handing ClipboardItem a Promise, instead
  // of awaiting the (slow) build first and writing afterwards. For a big gene like NF1 the fetch +
  // html2canvas render can exceed the browser's ~5 s transient-activation window; the old code
  // reached navigator.clipboard.write() after that window closed, so the write was rejected and the
  // OS clipboard kept the PREVIOUS variant — pasting the wrong variant (e.g. TP53) into Word. The
  // browser holds this write pending until htmlBlobPromise resolves, so activation is never lost.
  let buildErr = null;
  const htmlBlobPromise = buildReportHtml(zygosity)
    .then(html => new Blob([html], { type: 'text/html' }))
    .catch(e => { buildErr = e; throw e; });

  try {
    await navigator.clipboard.write([new ClipboardItem({ 'text/html': htmlBlobPromise })]);
    if (btn) {
      btn.innerHTML = '<span class="btn-icon">✅</span> Copied!';
      setTimeout(resetLabel, 2500);
    }
  } catch (err) {
    resetLabel();
    if (buildErr) {
      alert('Failed to build report: ' + buildErr.message);
    } else {
      // The write itself failed (lost focus / activation). Critically, the clipboard may still hold
      // the PREVIOUS variant — warn loudly so a stale paste can't be mistaken for this one.
      alert('Copy FAILED — the clipboard may still contain the PREVIOUS variant. Do NOT paste.\n\n'
        + 'Click anywhere in this page to focus it, then press Copy To Analysis Summary again.');
    }
  }
}

/**
 * saveReportToFolder() — prompt zygosity, build the report, and save it as an .html file.
 * Under the Go binary (/api/* present): native folder picker → server writes the file to disk.
 * Otherwise (vite / file://): falls back to a normal browser download. Works in both.
 */
let _lastSaveDir = '';
async function saveReportToFolder() {
  const zygosity = await promptZygosity();
  if (!zygosity) return;

  const btn = document.getElementById('btnSaveReport');
  const originalLabel = btn?.innerHTML;
  if (btn) btn.innerHTML = '<span class="btn-icon">⏳</span> Preparing…';

  let htmlContent;
  try {
    htmlContent = await buildReportHtml(zygosity);
  } catch (e) {
    if (btn) btn.innerHTML = originalLabel;
    alert('Failed to build report: ' + e.message);
    return;
  }

  // Filename from locked primary state: GENE_chr-pos-ref-alt_DATE.html (sanitised).
  const date = new Date().toISOString().slice(0, 10);
  const c = data.coords.userCoord;
  const variantId = data.coords.hg38String
    || (c ? `${c.chrom}-${c.pos}-${c.ref}-${c.alt}` : 'variant');
  const filename = `${data.coords.gene || 'gene'}_${variantId}_${date}.html`
    .replace(/[^A-Za-z0-9._-]/g, '_');

  try {
    // Ask the local Go helper for a folder (native OS picker). 204 = cancelled; error/refused → fallback.
    let dir = _lastSaveDir;
    if (!dir) {
      const pickRes = await fetch('/api/pick-directory');
      if (pickRes.status === 204) { if (btn) btn.innerHTML = originalLabel; return; }
      if (!pickRes.ok) throw new Error('picker unavailable');
      dir = (await pickRes.json()).path || '';
      if (!dir) { if (btn) btn.innerHTML = originalLabel; return; }
      _lastSaveDir = dir;
    }

    const saveRes = await fetch('/api/save-file', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ directory: dir, filename, content: htmlContent }),
    });
    if (!saveRes.ok) throw new Error((await saveRes.text()) || `HTTP ${saveRes.status}`);

    if (btn) {
      btn.innerHTML = '<span class="btn-icon">✅</span> Saved!';
      setTimeout(() => { btn.innerHTML = originalLabel; }, 2500);
    }
  } catch (e) {
    // Not running under the Go binary (or the write failed) → browser download of the same HTML.
    console.warn('Save-to-folder API unavailable; downloading instead:', e);
    const blob = new Blob([htmlContent], { type: 'text/html' });
    const url = URL.createObjectURL(blob);
    const a = Object.assign(document.createElement('a'), { href: url, download: filename });
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    if (btn) {
      btn.innerHTML = '<span class="btn-icon">⬇</span> Downloaded';
      setTimeout(() => { btn.innerHTML = originalLabel; }, 2500);
    }
  }
}

/**
 * Gathers current variant state, wraps it in the strict AI interpretation prompt,
 * and copies it directly to the user's clipboard as plain text.
 */
async function copyAIPrompt() {
  // 1. Prompt for zygosity using the existing modal
  const zygosity = await new Promise(resolve => {
    const modal = document.getElementById('zygosityModal');
    if (!modal) {
      resolve('Unknown');
      return;
    }
    modal.style.display = 'flex';
    window.resolveZygosity = (val) => {
      modal.style.display = 'none';
      data.selectedZygosity = val || 'Unknown';
      resolve(val);
    };
  });

  if (!zygosity) return; // Exit if user cancels modal

  // 2. Helper functions to format the data (mirrors copyReport logic)
  const g = id => (document.getElementById(id)?.innerText || '-').replace(/[✅✔☑️]/g, '').trim();
  const { label } = classifyACMG();

  // Format gnomAD data cleanly for the AI
  let gnomadFreqLine = 'gnomAD: -';
  if (data.gnomad.detailed?.total?.populations) {
    const eas = data.gnomad.detailed.total.populations.eas?.all || { ac: 0, an: 0, hom: 0 };
    const all = data.gnomad.detailed.total.populations.overall?.all || { ac: 0, an: 0, hom: 0 };
    gnomadFreqLine = `gnomAD v4.1.1: ${eas.ac}/${eas.an} (${eas.hom} homozygous) in EAS, ${all.ac}/${all.an} (${all.hom} homozygous) in All populations. Popmax: ${g('mvPopmax')}`;
  }

  let proteinStr = g('dProtein').split(':').pop().trim();
  if (proteinStr.startsWith('p.') && !proteinStr.includes('(')) {
    proteinStr = 'p.(' + proteinStr.substring(2) + ')';
  }

  // 3. Compile the raw data block
  const rawData = `
Variant: ${zygosity} ${data.ensembl.transcript || '-'}(${data.coords.gene || '-'}):${data.coords.hgvs || '-'} ${proteinStr} ${g('eExon')}
Genomic Coordinates: GRCh38 ${g('dHg38')} | hg19 ${g('dHg19')}
ClinVar Status: ${g('cvSig')} ${g('cvStars')} (Accession: ${data.clinvar.accession || '-'})
ClinVar Alt Missense at Codon: ${g('cvAltCodon')}
Population Frequency: ${gnomadFreqLine}
Bioinformatics Predictions: REVEL: ${(data.scores.revel !== null && !isNaN(data.scores.revel)) ? parseFloat(data.scores.revel).toFixed(3) : '-'} | AlphaMissense: ${(data.scores.alphaMissenseScore !== null && data.scores.alphaMissenseScore !== undefined) ? parseFloat(data.scores.alphaMissenseScore).toFixed(3) : '-'} | SpliceAI: ${(data.scores.spliceAI !== null && !isNaN(data.scores.spliceAI)) ? parseFloat(data.scores.spliceAI).toFixed(3) : '-'}
Gene Constraint: pLI=${data.geneConstraint?.pli != null ? parseFloat(data.geneConstraint.pli).toFixed(2) : '-'} | LOEUF=${data.geneConstraint?.loeuf != null ? parseFloat(data.geneConstraint.loeuf).toFixed(3) : '-'} | mis_z=${data.geneConstraint?.mis_z != null ? parseFloat(data.geneConstraint.mis_z).toFixed(2) : '-'}
ClinGen Dosage: Haploinsufficiency=${data.dosage.haplo ?? '-'} | Triplosensitivity=${data.dosage.triplo ?? '-'}
Associated Conditions: ${data.associatedConditions?.length ? data.associatedConditions.map(c => c.name || c).join('; ') : '-'}
ACMG Codes Applied by Curator: ${[...selectedCodes].join(', ')}
Target Classification: ${label}
  `.trim();

  // 4. Assemble prompt from editable template
  const finalPrompt = getPromptTemplate('iem').replace('{{RAW_DATA}}', rawData);

  // 5. Copy to Clipboard (using writeText for plain text)
  navigator.clipboard.writeText(finalPrompt).then(() => {
    const btn = document.getElementById('btnCopyAIPrompt');
    if (btn) {
      const originalHTML = btn.innerHTML;
      btn.style.background = 'var(--teal)'; // Assumes var(--teal) exists in your CSS
      btn.style.borderColor = 'var(--teal)';
      btn.innerHTML = '<span class="btn-icon" style="filter: none">✅</span> IEM Prompt Copied!';

      // Reset button state after 2.5 seconds
      setTimeout(() => {
        btn.style.background = '';
        btn.style.borderColor = '';
        btn.innerHTML = originalHTML;
      }, 2500);
    }
  }).catch((err) => {
    console.error('Clipboard copy failed:', err);
    alert('Copy failed — please allow clipboard access in your browser.');
  });
}

// ── PROMPT TEMPLATE MANAGEMENT ────────────────────────────────────

const PROMPT_DEFAULTS = {

  gg: `You are an expert Clinical Variant Scientist and Bioinformatician. Your task is to write a clinical variant interpretation summary based on the raw data provided below.

### STRICT INTERPRETATION RULES:
Haploinsufficiency (HI) Scoring: You must accurately interpret the provided ClinGen HI score based on the following framework:
- Score 3: Sufficient evidence for haploinsufficiency.
- Score 2: Some/emerging evidence for haploinsufficiency.
- Score 1: Little evidence for haploinsufficiency.
- Score 0: No evidence available supporting haploinsufficiency.
- Score 30: Gene is associated with an autosomal recessive phenotype. Loss of a single copy is typically tolerated. CRITICAL: Do NOT state there is "sufficient evidence for haploinsufficiency" if the score is 30.
- Score 40: Dosage sensitivity is unlikely.

### INSTRUCTIONS:
Write a cohesive, professional paragraph suitable for a clinical genetics report. DO NOT use a rigid fill-in-the-blank template. Vary the sentence structure naturally while maintaining scientific accuracy and objective clinical tone.

Ensure the overall structure follows this logical flow:
1. Nomenclature & Context: State the gene, c./p. nomenclature, variant type (missense, nonsense, splice, etc.), and exon location. Mention the predicted mechanism if applicable (e.g., NMD for nonsense/frameshift, or exon skipping for splice variants).
2. Population Data: State the gnomAD frequency (mention if it is absent, or provide the allele counts if present).
3. Clinical Databases & Gene Evidence: Detail the ClinVar classification, accession number, and submission count. Mention the ClinGen haploinsufficiency score accurately based on the strict rules above (e.g., if the score is 30, state the gene is associated with an autosomal recessive phenotype rather than haploinsufficiency).
4. Related Variants: Mention other variants at the same codon or splice site if provided.
5. In Silico Predictions: Provide REVEL scores for missense variants, or SpliceAI scores for splice variants, along with their predicted impact.
6. Literature: Include a placeholder exactly like "[Insert relevant PMIDs and clinical history here]" so the curator can manually add literature details.
7. Conclusion: End the paragraph exactly with: "[Insert final classification here]".

### REFERENCE EXAMPLES (Mirror this GG Report style and structural flexibility):

Example 1 (Splice variant):
NF1 c.3315-2A>G affects the canonical splice acceptor site of intron 25. It is predicted to cause skipping of exon 26 and subsequent frameshift, leading to protein truncation or nonsense-mediated mRNA decay. NF1 is curated as a gene with sufficient evidence for haploinsufficiency in ClinGen (Curation ID: CCID: 007547). The variant is absent in control populations (gnomAD v4.1.1 and v2.1.1). It has been deposited as pathogenic / likely pathogenic in ClinVar (Accession: VCV000484086.18, 4 submissions). Other nucleotide changes affecting the same acceptor site with similar prediction scores from SpliceAI (c.3315-2A>C, c.3315-1G>T, c.3315-1G>A and c.3315-1G>C) have been reported as pathogenic / likely pathogenic in ClinVar (Accessions: VCV000917573.8, VCV002137978.4, VCV000657669.9, VCV000565393.7).

Example 2 (Nonsense variant):
ELP1 c.2824C>T p.(Arg942*) is a nonsense variant located in exon 26 (out of 37 exons) of the gene. It creates a premature stop codon and is predicted to cause loss of normal protein function by nonsense-mediated mRNA decay. This variant is present at very low frequency in population databases (gnomAD v3.1.2 (non-cancer): 1 in 147,786 alleles; gnomAD v4.1.1: 12 in 1,612,926 alleles). It has been deposited as pathogenic/likely pathogenic in ClinVar by five submitters (Accession: VCV000848810.12) and detected in a patient with medulloblastoma (PMID: 39184053).

Example 3 (Missense variant):
NR5A1 c.982G>A p.(Gly328Arg) is a missense variant located in exon 5 (out of 7 exons) of the gene, residing in the ligand-binding domain (LBD) of the protein. The variant is present at very low frequency in control populations (gnomAD v.4.1.1: total 2 in 1,600,104 alleles; Absent in East Asian). It has been deposited as pathogenic in ClinVar by a single submitter (VCV000958009.9), and reported in individuals with 46,XY disorders of sex development (DSD), with de novo occurrence in one of the reported cases (PMIDs: 30425642, 32738419). Other single nucleotide substitutions affecting the same codon resulting in the same amino acid change (c.982G>C p.(Gly328Arg): PMID: 31745530) or different amino acid changes (c.983G>T p.(Gly328Val): PMIDs: 22474171, 29935645 and ClinVar accession: VCV001699288.3; c.982G>T p.(Gly328Trp): PMID: 29935645) have also been reported in other individuals with DSD. Although the c.982G>A p.(Gly328Arg) variant is predicted to be deleterious by in silico analysis (REVEL score: 0.961), one study reported no significant change in protein expression level or transcriptional activity by TESCO luciferase assay in vitro (PMID: 32738419).

Example 4 (In-frame deletion):
ATRX c.1501_1503del is an in-frame deletion in the exon 9 of 35 of the ATRX gene. The small deletion is predicted to cause a loss of a single amino acid residue without disrupting the reading frame of the protein. This variant is absent from population controls (gnomAD v4.1.1) and clinical databases (ClinVar and HGMD Professional 2026.1) at the time of reporting. This variant was detected in the patient's mother in the heterozygous state.

Example 5 (Missense VUS):
GNAS c.569A>G p.(Tyr190Cys) is a missense variant located in exon 7. The variant is absent in population controls (gnomAD v4.1.1 and v2.1.1). It has been deposited in ClinVar as Conflicting classifications of pathogenicity, with one submission as likely pathogenic and one submission as uncertain significance (Accession: VCV002577811.4). In silico prediction using REVEL suggested that this variant is deleterious at a strong level (REVEL score 0.973). The variant has been reported in GNAS patient, but no clinical details was provided (PMID: 31886927).

### RAW VARIANT DATA TO INTERPRET:
{{RAW_DATA}}

Generate the final interpretation summary now based on the above raw data. Output ONLY the summary text, with no conversational filler or markdown formatting.`,

  iem: `You are an expert Clinical Variant Scientist and Bioinformatician. Your task is to write concise, professional variant interpretation summaries based on raw data extracted from a variant curator tool.

### STRICT NARRATIVE TEMPLATE:
You MUST format your response using the EXACT sentence structures below. Do not use alternative phrasing. Fill in the bracketed variables using the provided raw data. Omit a sentence only if the data for it is completely missing.

1. **Introduction:** "[Transcript]([Gene]):[c.Nomenclature] [p.Nomenclature] is a [variant type] variant located in exon [Exon Number] of the [Gene] gene."
2. **Population Frequency (If present):** "This variant is reported at low frequency in population databases (gnomAD v4.1.0 All populations: [Total AC] in [Total AN] alleles, [Total AF]; [Popmax Group]: [Popmax AC] in [Popmax AN] alleles, [Popmax AF])."
   *Alternative (If absent):* "This variant is absent from population databases (gnomAD v4.1.0)."
3. **ClinVar Status:** "This variant has been classified as [ClinVar Status] for [Disease, if provided] in ClinVar ([ClinVar Accession], [Submission Count] submissions)."
4. **HGMD Status:** "In HGMD Professional 2025.4, this variant has been reported previously to be disease-causing..." (Summarize related PMIDs and phenotypes here if provided).
5. **Alternative Missense (If provided):** "Other missense changes affecting the [p.Residue] residue has been reported to be [Status] in ClinVar..."
6. **In Silico Predictions:** "In silico prediction using REVEL suggested that this variant is [deleterious/benign] at [strong/moderate/supporting] level (REVEL score [Score])." (Adjust for SpliceAI if applicable).
7. **Conclusion:** "For these reasons, this variant is classified as [Target Classification]."

### EXAMPLES OF DESIRED OUTPUT:

**Example 1 (Present in Pop databases):**
Input Data: Heterozygous NM_002225.5(IVD):c.631A>G p.(Thr211Ala) Exon 6. gnomAD v4.1.0: 6/1,614,176 (0.000003717), EAS 6/44,878 (0.0001337). ClinVar: Pathogenic/Likely Pathogenic (Isovaleryl-CoA dehydrogenase deficiency, 5 subs). HGMD: CM2017928. PMIDs: 32505769, 33210480, 35095998. REVEL: 0.956. Target Classification: Pathogenic.
Summary: NM_002225.5(IVD):c.631A>G p.(Thr211Ala) is a missense variant located in exon 6 of the IVD gene. This variant is reported at low frequency in population databases (gnomAD v4.1.0 All populations: 6 in 1,614,176 alleles, 0.000003717; East Asian: 6 in 44,878 alleles, 0.0001337). This variant has been classified as pathogenic/likely pathogenic for Isovaleryl-CoA dehydrogenase deficiency in ClinVar (VCV002676194.8, five submissions). In HGMD Professional 2025.4, this variant has been reported previously to be disease-causing and, in the compound heterozygous state in at least three patients with isovaleric acidaemia (HGMD accession CM2017928; PMID: 32505769, 33210480 and 35095998). In silico prediction using REVEL suggested that this variant is deleterious at strong level (REVEL score 0.956). For these reasons, this variant is classified as pathogenic.

**Example 2 (Absent from Pop databases):**
Input Data: Heterozygous NM_000531.6(OTC):c.77G>T p.(Arg26Leu) Exon 2. gnomAD v4.1.0: 0/0. ClinVar: Likely Pathogenic (1 star). Alt missense at same codon: Trp, Pro, Gln reported. REVEL: 0.386. Target Classification: VUS.
Summary: NM_000531.6(OTC):c.77G>T p.(Arg26Leu) is a missense variant located in exon 2 of the OTC gene. This variant is absent from population databases (gnomAD v4.1.0). This variant has been classified as likely pathogenic in ClinVar (1 star). Other missense changes affecting the p.Arg26 residue have been reported in ClinVar (resulting in substitutions to tryptophan, proline, and glutamine). In silico prediction using REVEL suggested that this variant is benign (REVEL score 0.386). For these reasons, this variant is classified as a Variant of Uncertain Significance.

### TASK:
Generate a variant interpretation summary for the following raw input data. Do not invent missing data. Stick strictly to the narrative template provided.

{{RAW_DATA}}`,

  update: `You are an expert Clinical Variant Scientist and a meticulous Medical Editor. Your task is to update an older variant interpretation paragraph using the most recently pulled database information.

### INSTRUCTIONS AND STRICT RULES:
1. SURGICAL UPDATES ONLY: Cross-reference the "Old Report" with the "New Variant Data". ONLY update numbers, facts, classifications, or database versions (e.g., gnomAD versions, allele counts, ClinVar status/submissions) that are explicitly mentioned in the Old Report but have changed in the New Data.
2. HIGHLIGHT UPDATES: In the updated paragraph, please **bold** or [bracket] the specific text or values that you have changed so the curator can easily spot them.
3. DO NOT ADD NEW CONCEPTS: If the New Data contains information that is NOT discussed in the Old Report, DO NOT add it. Do not add new sentences.
4. GRAMMAR & CLARITY: Correct any grammatical errors or spelling mistakes in the Old Report.
5. PRESERVE THE ORIGINAL VOICE: Keep the original sentence structure and clinical history/PMID summaries exactly as they are.
6. SUMMARY TABLE: After the updated paragraph, provide a Markdown Table titled "Summary of Changes" with columns: [Category, Old Value, New Value, Rationale].

### NEW VARIANT DATA (Use this as your source of truth):
{{RAW_DATA}}

### OLD REPORT (Update this text):
"""
{{OLD_REPORT}}
"""

Generate the updated interpretation paragraph first, followed by the "Summary of Changes" table. Output only the requested content without conversational filler.`

};

function getPromptTemplate(type) {
  return localStorage.getItem(`aiPromptTemplate_${type}`) || PROMPT_DEFAULTS[type];
}
function savePromptTemplate(type, text) {
  localStorage.setItem(`aiPromptTemplate_${type}`, text);
}
function resetPromptTemplate(type) {
  localStorage.removeItem(`aiPromptTemplate_${type}`);
}

const _PROMPT_EDITOR_HINTS = {
  gg:     'Keep {{RAW_DATA}} where the auto-generated variant data should be inserted.',
  iem:    'Keep {{RAW_DATA}} where the auto-generated variant data should be inserted.',
  update: 'Keep {{RAW_DATA}} (new variant data) and {{OLD_REPORT}} (the pasted old report) in the template.',
};
const _PROMPT_EDITOR_TITLES = {
  gg:     'Edit GG Report Prompt Template',
  iem:    'Edit IEM Report Prompt Template',
  update: 'Edit Update Report Prompt Template',
};
let _promptEditorCurrentType = null;

function openPromptEditor(type) {
  _promptEditorCurrentType = type;
  const modal    = document.getElementById('promptEditorModal');
  const textarea = document.getElementById('promptEditorTextarea');
  const title    = document.getElementById('promptEditorTitle');
  const hint     = document.getElementById('promptEditorHint');
  if (!modal || !textarea) return;
  title.textContent = _PROMPT_EDITOR_TITLES[type] || 'Edit Prompt Template';
  hint.textContent  = _PROMPT_EDITOR_HINTS[type]  || 'Keep {{RAW_DATA}} in the template.';
  textarea.value    = getPromptTemplate(type);
  // Mark edit button as custom if overridden
  _updateEditBtnIndicator(type);
  modal.style.display = 'flex';
  setTimeout(() => textarea.focus(), 50);
}

function closePromptEditor() {
  const modal = document.getElementById('promptEditorModal');
  if (modal) modal.style.display = 'none';
  _promptEditorCurrentType = null;
}

function saveCurrentTemplate() {
  if (!_promptEditorCurrentType) return;
  const textarea = document.getElementById('promptEditorTextarea');
  if (!textarea) return;
  const text = textarea.value.trim();
  if (!text) { alert('Template cannot be empty.'); return; }
  if (!text.includes('{{RAW_DATA}}')) {
    if (!confirm('Warning: {{RAW_DATA}} placeholder not found — variant data will NOT be inserted. Save anyway?')) return;
  }
  if (_promptEditorCurrentType === 'update' && !text.includes('{{OLD_REPORT}}')) {
    if (!confirm('Warning: {{OLD_REPORT}} placeholder not found — the old report text will NOT be inserted. Save anyway?')) return;
  }
  savePromptTemplate(_promptEditorCurrentType, text);
  _updateEditBtnIndicator(_promptEditorCurrentType);
  closePromptEditor();
}

function resetCurrentTemplate() {
  if (!_promptEditorCurrentType) return;
  if (!confirm('Reset to the built-in default template? Your customisations will be lost.')) return;
  resetPromptTemplate(_promptEditorCurrentType);
  const textarea = document.getElementById('promptEditorTextarea');
  if (textarea) textarea.value = PROMPT_DEFAULTS[_promptEditorCurrentType];
  _updateEditBtnIndicator(_promptEditorCurrentType);
}

/** Tints the ✏️ button teal when a custom template is saved, grey when default. */
function _updateEditBtnIndicator(type) {
  const ids = { gg: 'btnEditPromptGG', iem: 'btnEditPromptIEM', update: 'btnEditPromptUpdate' };
  const btn = document.getElementById(ids[type]);
  if (!btn) return;
  const isCustom = !!localStorage.getItem(`aiPromptTemplate_${type}`);
  btn.style.color       = isCustom ? 'var(--teal)' : 'var(--dim)';
  btn.style.borderColor = isCustom ? 'var(--teal)' : 'var(--border)';
  btn.title = isCustom
    ? 'Custom prompt template active — click to edit'
    : ((_PROMPT_EDITOR_TITLES[type] || 'Edit Prompt Template').replace('Edit ', 'Edit default '));
}

function buildVariantBioinformaticsBlock(vc) {
  const s = data.scores;
  const gc = data.geneConstraint || {};
  const fmt = (v, d = 3) => (v != null && !isNaN(v)) ? parseFloat(v).toFixed(d) : 'N/A';

  const isMissense = vc?.includes('missense_variant');
  const isSplice   = ['splice_donor', 'splice_acceptor', 'splice_region'].some(t => vc?.includes(t));
  const isLoF      = ['stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost'].some(t => vc?.includes(t));
  const isInframe  = vc?.includes('inframe_insertion') || vc?.includes('inframe_deletion');

  if (isMissense) {
    const amPred = s.alphaMissensePred ? ` (${s.alphaMissensePred})` : '';
    const spliceNote = (s.spliceAI != null && parseFloat(s.spliceAI) >= 0.1)
      ? ` | SpliceAI: ${fmt(s.spliceAI)} (possible splice-region impact)` : '';
    const domainHit = (data.uniprotDomains?.length && data.vepProteinStart)
      ? data.uniprotDomains.find(d => data.vepProteinStart >= d.start && data.vepProteinStart <= d.end)
      : null;
    const domainLine = domainHit ? `\n  UniProt domain: ${domainHit.description} (aa ${domainHit.start}–${domainHit.end})` : '';
    return [
      `REVEL: ${fmt(s.revel)} | AlphaMissense: ${fmt(s.alphaMissenseScore)}${amPred}${spliceNote}`,
      `Gene constraint (missense): mis_z=${fmt(gc.mis_z, 2)}, pLI=${fmt(gc.pli, 2)}, LOEUF=${fmt(gc.loeuf, 3)}${domainLine}`
    ].join('\n  ');
  }

  if (isSplice) {
    const THRESH = 0.2;
    const deltaMap = [
      ['DS_AG', s.spliceAI_AG, s.spliceAI_DP_AG],
      ['DS_AL', s.spliceAI_AL, s.spliceAI_DP_AL],
      ['DS_DG', s.spliceAI_DG, s.spliceAI_DP_DG],
      ['DS_DL', s.spliceAI_DL, s.spliceAI_DP_DL],
    ];
    const significant = deltaMap
      .filter(([, v]) => v != null && parseFloat(v) >= THRESH)
      .map(([name, v, pos]) => `${name}=${fmt(v)}${pos != null ? `@${pos}` : ''}`)
      .join(', ') || 'No delta score ≥0.2';
    const mechanism = getSpliceConsequence() || 'Predicted splice effect — see SpliceAI deltas';
    return [
      `SpliceAI: max=${fmt(s.spliceAI)} — Significant (≥0.2): ${significant}`,
      `Predicted mechanism: ${mechanism}`,
      `Gene constraint (LoF): pLI=${fmt(gc.pli, 2)}, LOEUF=${fmt(gc.loeuf, 3)}`
    ].join('\n  ');
  }

  if (isLoF) {
    const exonNum    = data.vepExon ?? data.ensembl?.vepExon;
    const totalExons = data.vepTotalExons ?? data.ensembl?.vepTotalExons;
    let nmd = 'NMD: Cannot determine (exon data unavailable)';
    if (exonNum != null && totalExons != null) {
      const en = parseInt(String(exonNum).split('/')[0]);
      const te = parseInt(String(totalExons));
      nmd = en === te
        ? 'NMD: Escape predicted (variant in last exon)'
        : `NMD: Predicted (PTC in exon ${en}/${te} — upstream of last exon-exon junction)`;
    }
    const hi   = data.dosage.haplo ?? '-';
    const hiId = data.dosage.haploId ? ` (Curation ID: ${data.dosage.haploId})` : '';
    const hiLine = `ClinGen Haploinsufficiency: Score ${hi}${hiId}`;
    const grChapters = data.geneReviews?.chapters || [];
    const focusNbk   = data.geneReviews?.focusNbk;
    const grChapter  = (focusNbk ? grChapters.find(c => c.nbk === focusNbk) : null) ?? grChapters[0];
    let grLofLine = '';
    if (grChapter?.mechanism) {
      const m = grChapter.mechanism.toLowerCase();
      if (/loss.of.function|haploinsufficiency|\blof\b|null variant|truncat|nonsense.mediated/.test(m)) {
        const ref = grChapter.url || `NBK${grChapter.nbk}`;
        grLofLine = `\n  GeneReviews: LoF is disease mechanism — ${grChapter.title} (${ref})`;
      }
    }
    return [
      nmd,
      hiLine + grLofLine,
      `Gene constraint (LoF): pLI=${fmt(gc.pli, 2)}, LOEUF=${fmt(gc.loeuf, 3)}`
    ].join('\n  ');
  }

  if (isInframe) {
    return `Gene constraint: pLI=${fmt(gc.pli, 2)}, LOEUF=${fmt(gc.loeuf, 3)}, mis_z=${fmt(gc.mis_z, 2)}`;
  }

  // Synonymous or unknown — SpliceAI only if elevated
  if (s.spliceAI != null && parseFloat(s.spliceAI) >= 0.1) {
    return `SpliceAI: ${fmt(s.spliceAI)} (elevated — possible cryptic splice effect)`;
  }
  return 'N/A (not applicable for this variant type)';
}

/**
 * buildPhenotypeFitBlock()
 * Formats patient HPO terms + the top phenotype-fit disease LRs for the GG prompt.
 * Reads data.ptPhenoTexts + data.phenotypeFit (prompt-string assembly only). Returns ''
 * when no phenotypes are entered. LR-led (prior-independent); flags sparse-frequency rows.
 * First place patient phenotypes enter the AI prompt — grounds disease-context reasoning.
 */
function buildPhenotypeFitBlock() {
  const terms = (data.ptPhenoTexts || []).filter(Boolean);
  if (!terms.length) return '';
  let block = `Patient Phenotypes (HPO): ${terms.join(', ')}`;
  const fit = data.phenotypeFit;
  if (fit && fit.ready && Array.isArray(fit.results) && fit.results.length) {
    const lines = fit.results.slice(0, 3).map(r => {
      const sparse = (r.completeness != null && r.completeness < 0.5)
        ? ' [sparse frequency data — weigh coverage over LR]' : '';
      return `    - ${r.name} (MIM ${r.mimId}): log10LR ${r.log10LR.toFixed(2)}, coverage ${Math.round((r.coverage || 0) * 100)}%${sparse}`;
    });
    block += `\nPhenotype-Fit (ontology-aware HPO likelihood ratio vs population background; higher log10LR = better fit, prior-independent):\n${lines.join('\n')}`;
  }
  return block;
}

/**
 * Gathers variant state and formats a few-shot prompt for an AI
 * to generate a GG-style variant interpretation report.
 */
async function copyAIPromptGG() {
  // 1. Prompt for zygosity using the existing modal
  const zygosity = await new Promise(resolve => {
    const modal = document.getElementById('zygosityModal');
    if (!modal) {
      resolve('Unknown');
      return;
    }
    modal.style.display = 'flex';
    window.resolveZygosity = (val) => {
      modal.style.display = 'none';
      data.selectedZygosity = val || 'Unknown';
      resolve(val);
    };
  });

  if (!zygosity) return; // Exit if user cancels modal

  // 2. Helper functions to format the data
  const g = id => document.getElementById(id)?.innerText || '-';
  const { label } = classifyACMG();

  // Format gnomAD data cleanly
  let gnomadFreqLine = 'Absent in population databases (gnomAD v4.1.1)';
  if (data.gnomad.detailed?.total?.populations) {
    const eas = data.gnomad.detailed.total.populations.eas?.all || { ac: 0, an: 0 };
    const all = data.gnomad.detailed.total.populations.overall?.all || { ac: 0, an: 0 };
    if (all.ac > 0) {
      gnomadFreqLine = `gnomAD v4.1.1: Total ${all.ac} in ${all.an} alleles; East Asian: ${eas.ac > 0 ? `${eas.ac} in ${eas.an} alleles` : 'Absent'}`;
    }
  }

  // Format Protein nomenclature strictly
  let proteinStr = g('dProtein').split(':').pop().replace(/✅/g, '').trim();
  if (proteinStr.startsWith('p.') && !proteinStr.includes('(')) {
    proteinStr = 'p.(' + proteinStr.substring(2) + ')';
  }

  // 3. Compile the raw data block for the AI
  // Canonical consequence is data.ensembl.vepConsequence; top-level data.vepConsequence is never
  // assigned (would leave vc undefined → buildVariantBioinformaticsBlock falls through to its N/A
  // default and "Variant Type:" reads Unknown). Mirrors the fallback pattern used elsewhere.
  const vc = data.ensembl?.vepConsequence || data.vepConsequence;
  const isMissenseVC  = vc?.includes('missense_variant');
  const isSpliceVC    = ['splice_donor', 'splice_acceptor', 'splice_region'].some(t => vc?.includes(t));
  const altCodonLine  = isMissenseVC
    ? `ClinVar Alt Missense at Codon: ${g('cvAltCodon')}`
    : isSpliceVC
      ? `ClinVar Alt Variants at Splice Site: ${g('cvAltCodon')}`
      : '';

  const phenoFitBlock = buildPhenotypeFitBlock();

  const rawData = `
Gene: ${data.coords.gene || '-'}
Transcript: ${data.ensembl.transcript || '-'}
HGVS.c: ${data.coords.hgvs || '-'}
HGVS.p: ${proteinStr}
Exon: ${g('eExon').replace(/✅ |✔ |☑️ /g, '')} out of ${data.vepTotalExons || '?'} exons
Zygosity: ${zygosity}
Variant Type: ${vc || 'Unknown'}
ClinGen Dosage: Haploinsufficiency: ${data.dosage.haplo || '-'} | Triplosensitivity: ${data.dosage.triplo || '-'}
Population Frequency: ${gnomadFreqLine}
ClinVar Status: ${g('cvSig')} (${g('cvStars')} stars; Accession: ${data.clinvar.accession || '-'}, Submissions: ${data.clinvar.subs || '0'})
${phenoFitBlock ? phenoFitBlock + '\n' : ''}${altCodonLine ? altCodonLine + '\n' : ''}Bioinformatics:
  ${buildVariantBioinformaticsBlock(vc)}
  `.trim();

  // 4. Assemble prompt from editable template
  const finalPrompt = getPromptTemplate('gg').replace('{{RAW_DATA}}', rawData);

  // 5. Copy to Clipboard
  navigator.clipboard.writeText(finalPrompt).then(() => {
    const btn = document.getElementById('btnCopyAIPromptGG');
    if (btn) {
      const originalHTML = btn.innerHTML;
      btn.style.background = 'var(--teal)';
      btn.style.borderColor = 'var(--teal)';
      btn.style.color = '#000';
      btn.innerHTML = '<span class="btn-icon" style="filter: none">✅</span> GG AI Prompt Copied!';

      setTimeout(() => {
        btn.style.background = '';
        btn.style.borderColor = '';
        btn.style.color = '';
        btn.innerHTML = originalHTML;
      }, 2500);
    }
  }).catch((err) => {
    console.error('Clipboard copy failed:', err);
    alert('Copy failed — please allow clipboard access in your browser.');
  });
}

/**
 * Opens the Update Report modal and clears previous text.
 */
function openUpdateReportModal() {
  const modal = document.getElementById('updateReportModal');
  const textarea = document.getElementById('oldReportTextarea');
  if (textarea) textarea.value = ''; // Clear old paste
  if (modal) modal.style.display = 'flex';
}

/**
 * Closes the Update Report modal.
 */
function closeUpdateReportModal() {
  const modal = document.getElementById('updateReportModal');
  if (modal) modal.style.display = 'none';
}

/**
 * Grabs the pasted text, extracts current variant data, and builds a strict 
 * editor prompt for the LLM.
 */
async function generateUpdateAIPrompt() {
  const textarea = document.getElementById('oldReportTextarea');
  const oldReportText = textarea ? textarea.value.trim() : '';

  if (!oldReportText) {
    alert("Please paste an old report to update.");
    return;
  }

  // 1. Prompt for zygosity (required for the new data block)
  const zygosity = await new Promise(resolve => {
    // Hide update modal temporarily while getting zygosity
    closeUpdateReportModal();

    const zModal = document.getElementById('zygosityModal');
    if (!zModal) return resolve('Unknown');
    zModal.style.display = 'flex';
    window.resolveZygosity = (val) => {
      zModal.style.display = 'none';
      data.selectedZygosity = val || 'Unknown';
      resolve(val);
    };
  });

  if (!zygosity) return; // Exit if cancelled

  // 2. Format Current Data (same as Analysis Summary extraction)
  const g = id => document.getElementById(id)?.innerText || '-';

  // Format gnomAD data cleanly (v4.1.1)
  let gnomadFreqLine = 'Absent in population databases (gnomAD v4.1.1)';
  if (data.gnomad.detailed?.total?.populations) {
    const eas = data.gnomad.detailed.total.populations.eas?.all || { ac: 0, an: 0 };
    const all = data.gnomad.detailed.total.populations.overall?.all || { ac: 0, an: 0 };
    if (all.ac > 0) {
      gnomadFreqLine = `gnomAD v4.1.1: Total ${all.ac} in ${all.an} alleles; East Asian: ${eas.ac > 0 ? `${eas.ac} in ${eas.an} alleles` : 'Absent'}`;
    }
  }

  let proteinStr = g('dProtein').split(':').pop().replace(/✅/g, '').trim();
  if (proteinStr.startsWith('p.') && !proteinStr.includes('(')) {
    proteinStr = 'p.(' + proteinStr.substring(2) + ')';
  }

  // 3. Compile the NEW basic raw data block
  const updateVc = data.vepConsequence;
  const rawData = `
Gene: ${data.coords.gene || '-'}
Transcript: ${data.ensembl.transcript || '-'}
HGVS.c: ${data.coords.hgvs || '-'}
HGVS.p: ${proteinStr}
Variant Type: ${updateVc || 'Unknown'}
Population Frequency: ${gnomadFreqLine}
ClinVar Status: ${g('cvSig')} (${g('cvStars')} stars; Accession: ${data.clinvar.accession || '-'}, Submissions: ${data.clinvar.subs || '0'})
Bioinformatics:
  ${buildVariantBioinformaticsBlock(updateVc)}
  `.trim();

  // 4. Assemble prompt from editable template
  const finalPrompt = getPromptTemplate('update')
    .replace('{{RAW_DATA}}', rawData)
    .replace('{{OLD_REPORT}}', oldReportText);

  // 5. Copy to Clipboard
  navigator.clipboard.writeText(finalPrompt).then(() => {
    const btn = document.getElementById('btnUpdateReport');
    if (btn) {
      const originalHTML = btn.innerHTML;
      btn.style.background = 'var(--teal)';
      btn.style.borderColor = 'var(--teal)';
      btn.style.color = '#000';
      btn.innerHTML = '<span class="btn-icon" style="filter: none">✅</span> Update Prompt Copied!';

      setTimeout(() => {
        btn.style.background = '';
        btn.style.borderColor = '';
        btn.style.color = '';
        btn.innerHTML = originalHTML;
      }, 2500);
    }
  }).catch((err) => {
    console.error('Clipboard copy failed:', err);
    alert('Copy failed — please allow clipboard access in your browser.');
  });
}

// ── INITIALIZATION ────────────────────────────────────────────────

function toggleTheme() {
  const isLavender = document.body.dataset.theme === 'lavender';
  const next = isLavender ? '' : 'lavender';
  document.body.dataset.theme = next;
  const icon  = document.getElementById('themeToggleIcon');
  const label = document.getElementById('themeToggleLabel');
  if (icon)  icon.textContent  = next === 'lavender' ? '🌙' : '🌸';
  if (label) label.textContent = next === 'lavender' ? 'Dark'  : 'Lavender';
  localStorage.setItem('preferredTheme', next);
}

function initTheme() {
  const saved = localStorage.getItem('preferredTheme');
  if (saved === 'lavender') {
    document.body.dataset.theme = 'lavender';
    const icon  = document.getElementById('themeToggleIcon');
    const label = document.getElementById('themeToggleLabel');
    if (icon)  icon.textContent  = '🌙';
    if (label) label.textContent = 'Dark';
  }
}

document.addEventListener('DOMContentLoaded', () => {
  buildACMGMatrix(); buildDBGrid(); fetchEnsemblVersion(); initializeAPIBadges();
  prefetchClinGenSVIS();
  warmupColdStartServices();   // pre-spin SpliceAI/VariantValidator cold starts off the first query
  loadHistory();
  initTheme();
  // Restore custom-template tint on edit buttons
  ['gg', 'iem', 'update'].forEach(_updateEditBtnIndicator);

  // Restore API Mode Preference
  const savedMode = localStorage.getItem('primaryApiMode');
  if (savedMode) {
    const sel = document.getElementById('apiSelector');
    if (sel) sel.value = savedMode;
  }

  const urlParams = new URLSearchParams(window.location.search);
  const q = urlParams.get('query'); if (q) { document.getElementById('rawInput').value = q; handleInput(); }
  const g = urlParams.get('gene'); if (g) { document.getElementById('geneInput').value = g; handleGeneInput(); }
  // Enter key in variant input bar silenced — no longer triggers database search
});

function saveHistory() {
  localStorage.setItem('clinicalChatHistory', JSON.stringify(data.clinicalChatHistory));
}

function loadHistory() {
  const saved = localStorage.getItem('clinicalChatHistory');
  if (saved) {
    data.clinicalChatHistory = JSON.parse(saved);
    renderChatMessages('clinical');
  }
}

function getSpliceConsequence() {
  const THRESHOLD = 0.20;
  // Extract strand from global data
  const { spliceAI_AG: ag, spliceAI_AL: al, spliceAI_DG: dg, spliceAI_DL: dl,
    spliceAI_DP_AG: pag, spliceAI_DP_AL: pal, spliceAI_DP_DG: pdg, spliceAI_DP_DL: pdl,
    spliceAiManeMatch: isMane, vepExon, vepTotalExons, strand } = data;

  const scores = [ag || 0, al || 0, dg || 0, dl || 0];
  if (Math.max(...scores) < THRESHOLD) return null;

  let res = isMane ? "" : "(Non-MANE) ";
  let impact = "";
  let dist = 0;

  // GEOMETRIC FIX: SpliceAI coordinates are forward-strand anchored.
  const isMinusStrand = strand === -1;

  // 1. Cryptic Donor
  if (dl > THRESHOLD && dg > THRESHOLD) {
    dist = Math.abs(pdg - pdl);

    // Default Plus Strand: pdg > pdl means the new donor is downstream (Intron Retention).
    let isIntronRet = pdg > pdl;
    if (isMinusStrand) isIntronRet = !isIntronRet; // Invert for Minus Strand

    impact = `Cryptic Donor: ${isIntronRet ? 'Intron Ret.' : 'Exon Del.'} (${dist}bp) - ${dist % 3 === 0 ? 'Predicted In-frame' : 'Predicted Frameshift'}`;
  }
  // 2. Cryptic Acceptor
  else if (al > THRESHOLD && ag > THRESHOLD) {
    dist = Math.abs(pag - pal);

    // Default Plus Strand: pag < pal means the new acceptor is upstream (Intron Retention).
    let isIntronRet = pag < pal;
    if (isMinusStrand) isIntronRet = !isIntronRet; // Invert for Minus Strand

    impact = `Cryptic Acceptor: ${isIntronRet ? 'Intron Ret.' : 'Exon Del.'} (${dist}bp) - ${dist % 3 === 0 ? 'Predicted In-frame' : 'Predicted Frameshift'}`;
  }
  // 3. Exon Skipping
  else if (dl > THRESHOLD || al > THRESHOLD) {
    impact = `Exon Skipping (${dl > al ? 'Donor' : 'Acceptor'} Loss)`;
    // Warning added because we don't know the exact length of the skipped exon purely from VEP
    impact += " (Frame impact requires manual review)";
  }

  if (!impact) return null;
  res += impact;

  // 4. NMD 50-Nucleotide Rule Heuristic
  if (impact.includes("Frameshift") && vepExon && vepTotalExons) {
    const currentExon = parseInt(vepExon, 10);
    const totalExons = parseInt(vepTotalExons, 10);

    if (currentExon < totalExons - 1) {
      res += " - Predicted NMD";
    } else if (currentExon === totalExons - 1) {
      res += " - Possible NMD (Check 50-nt rule manually)";
    } else if (currentExon === totalExons) {
      res += " - Escapes NMD (Last Exon)";
    }
  }

  return res;
}

// Global initialization
// loadClinGenTSV removed in favor of static dictionary

/**
 * Parses a combined VCF/HGVS format for zero-latency UI population.
 * Example: chr7 140753336 A T; NM_004333.6(BRAF):c.1799T>A p.(Val600Glu)
 */
function parseLocalVCFString(inputStr) {
  try {
    const parts = inputStr.split(';');
    if (parts.length !== 2) return false;

    const coords = parts[0].trim().split(/\s+/);
    if (coords.length < 4) return false;

    data.coords.chrom = coords[0].replace(/chr/i, '');
    data.coords.pos38 = coords[1];
    data.coords.ref = coords[2].toUpperCase();
    data.coords.alt = coords[3].toUpperCase();
    data.coords.hg38String = `chr${data.coords.chrom}-${data.coords.pos38}-${data.coords.ref}-${data.coords.alt}`;

    const hgvsPart = parts[1].trim();
    const hgvsRegex = /([^()]+)\(([^)]+)\):(c\.[^\s]+)\s*(p\..+)?/;
    const match = hgvsPart.match(hgvsRegex);

    if (!match) return false;

    data.ensembl.transcript = match[1];
    data.coords.gene = match[2].toUpperCase();
    data.coords.hgvs = match[3];
    data.ensembl.protein = match[4] || '-';

    if (data.ensembl.protein !== '-') {
      const pMatch = data.ensembl.protein.match(/p\.\(([a-zA-Z]{3})(\d+)/);
      if (pMatch) {
        data.vepProteinStart = pMatch[2];
        const ref3 = pMatch[1];
        data.vepAminoAcids = AA_MAP_3TO1[ref3] || null;
      }
    }

    document.getElementById('dGene').innerText = data.coords.gene;
    document.getElementById('geneInput').value = data.coords.gene;
    document.getElementById('eManeTranscript').innerText = data.ensembl.transcript;
    document.getElementById('dHgvs').innerText = data.coords.hgvs;
    document.getElementById('dProtein').innerText = data.ensembl.protein;
    document.getElementById('dHg38').innerText = data.coords.hg38String;

    return true;
  } catch (error) {
    console.error("Local Parsing Failed:", error);
    return false;
  }
}

/**
 * Consolidates secondary background API calls to ensure consistent UI population.
 */
// Called from both fetchClinVar (after cvStars set) and fetchCanonicalAllele.then (after caId set).
// Fires only when both conditions are met; deduplicates with undefined/null sentinel.
function tryFetchClinGenERepoLink() {
  console.log('[eRepo] tryFetch — cvStars:', data.clinvar.stars, 'gene:', data.coords.gene, 'hgvs:', data.coords.hgvs, 'cvInterpretationId:', data.clinvar.interpretationId);
  if (data.clinvar.interpretationId !== undefined) { console.log('[eRepo] skip: already done'); return; }
  if (data.clinvar.stars < 3 || !data.coords.gene || !data.coords.hgvs || data.coords.hgvs === '-') { console.log('[eRepo] skip: condition not met'); return; }
  data.clinvar.interpretationId = null;
  console.log('[eRepo] Querying with gene:', data.coords.gene, 'hgvs:', data.coords.hgvs);
  fetchClinGenInterpretationId(data.coords.gene, data.coords.hgvs).then(id => {
    console.log('[eRepo] Got UUID:', id);
    if (!id) return;
    data.clinvar.interpretationId = id;
    const eRepoUrl = `https://erepo.clinicalgenome.org/evrepo/ui/interpretation/${id}`;
    const caEl = document.getElementById('cvAcc');
    if (caEl) caEl.innerHTML = `<a href="${eRepoUrl}" target="_blank" style="color:inherit;text-decoration:underline;" title="ClinGen Expert Panel Review">${data.clinvar.accession}</a>`;
    const evEl = document.getElementById('cvEvidenceLink');
    console.log('[eRepo] cvEvidenceLink element:', !!evEl);
    if (evEl) { evEl.href = eRepoUrl; evEl.title = 'ClinGen Expert Panel Review'; }
  });
}

// Re-render the gene-disease condition cards with the authoritative ClinGen validity
// overlay. Safe to call whenever EITHER fetchAssociatedConditions OR fetchGeneValidity
// resolves — whichever lands last produces the reconciled view. No-op until conditions exist.
function reconcileConditionCards(gen) {
  if (gen != null && gen !== window.currentSearchGen()) return;
  if (!data.associatedConditions?.length || !data.coords.gene) return;
  renderConditionCards(data.associatedConditions, data.coords.gene, data.geneValidity?.curations || []);
}

// Curation-support orchestrator (CLAUDE.md Rule 10 — cross-module coordination lives only at the top
// level). Reads the variant consequence, the focused GeneReviews chapter (mechanism + frequency
// table) and the patient phenotype, passes them EXPLICITLY into the pure evaluators (Rule 2), stores
// the result in its own namespace and renders. Idempotent — safe to call whenever inputs change.
function refreshCurationSupport(gen) {
  if (gen != null && gen !== window.currentSearchGen()) return;
  const gr = data.geneReviews;
  const chapter = gr?.chapters?.length
    ? ((gr.focusNbk && gr.chapters.find(c => c.nbk === gr.focusNbk)) || gr.chapters[0])
    : null;
  const mechanism = evaluateMechanismConcordance({
    vepConsequence: data.ensembl?.vepConsequence || data.vepConsequence || null,
    mechanism: chapter?.mechanism || '',
    haplo: data.dosage?.haplo ?? null,
  });
  // Phenotype↔frequency: match across ALL of the gene's chapters (not just the focused one), so a
  // patient feature is found wherever GeneReviews documents it (e.g. focus picks a chapter with a
  // sparse table while the main chapter's table carries the feature). Focused chapter's rows go
  // first, so on a tie the focused disease's frequency wins; each row is tagged with its disease.
  const orderedChapters = chapter
    ? [chapter, ...(gr?.chapters || []).filter(c => c !== chapter)]
    : (gr?.chapters || []);
  const mergedFreqs = orderedChapters.flatMap(c =>
    (c.phenotypeFrequencies || []).map(f => ({ ...f, disease: c.title || null })));
  const phenoMatches = matchPhenotypeFrequencies({
    phenotypeFrequencies: mergedFreqs,
    ptPhenoTexts: data.ptPhenoTexts || [],
    ptPhenoGroups: data.ptPhenoGroups || [],
  });
  data.curationSupport = {
    mechanism, phenoMatches,
    chapterTitle: chapter?.title || null,
    chapterUrl: chapter?.url || null,
  };
  renderCurationSupportCard();
}

function runGeneOnlyMode(gen) {
  document.getElementById('statusMsg').innerHTML =
    `<span style="color:var(--teal);">Gene-only mode — gene pathogenicity, variant distribution, protein viewer &amp; phenotype analysis</span>`;
  enableWrappers(['clingen', 'clinvar', 'omim', 'gr', 'scholar', 'decipher', 'hgmd', 'gtex']);
  checkVCEP(data.coords.gene);

  // Always trigger (regardless of phenotype)
  fetchClinVarGeneDistribution(data.coords.gene);
  fetchGnomADConstraint(gen);
  fetchDosageSensitivity(null, null, data.coords.gene);
  const gvPromise = fetchGeneValidity(data.coords.gene).then(() => {
    renderGeneValidityCard();             // dedicated live card (own namespace)
    reconcileConditionCards(gen);         // overlay authoritative classification onto condition cards
    return enrichGeneValidityEvidence(data.geneValidity?.curations);
  }).then(() => renderGeneValidityCard()).catch(() => {});
  const grPromise = fetchGeneReviews(data.coords.gene, gen).catch(() => {});
  grPromise.then(() => refreshCurationSupport(gen)).catch(() => {});   // phenotype↔frequency (no variant in gene-only)
  // Age of onset for definitive curations (OMIM synopsis primary, GeneReviews fallback once chapters arrive).
  Promise.all([gvPromise, grPromise]).then(() =>
    enrichGeneValidityOnset(data.geneValidity?.curations, data.coords.gene, data.geneReviews?.chapters)
  ).then(() => { if (gen === window.currentSearchGen()) { renderGeneValidityCard(); maybeEvaluateBS2(gen); } }).catch(() => {});

  // Protein Variant Viewer (Phase A + B)
  const domainPromise = fetchUniprotDomains(data.coords.gene, gen).then(domains => {
    if (gen !== window.currentSearchGen()) return [];
    data.uniprotDomains = domains;
    return domains;
  }).catch(err => {
    console.warn('[Domain Annotation] Error:', err);
    data.uniprotDomains = [];
    return [];
  });

  Promise.all([
    fetchCodonData(data.coords.gene),
    domainPromise,
  ]).then(([codonData, uniprotDomains]) => {
    if (gen !== window.currentSearchGen()) return;
    const date = codonData?._meta?.date || null;
    const capturedGene = data.coords.gene;
    renderCodonViewer(capturedGene, codonData, uniprotDomains || [], date, { queryCodon: null, exons: [] });
    // Phase B: fetch exon boundaries via gene symbol (no transcript in gene-only mode)
    fetchExonCodingPositions(null, capturedGene).then(exons => {
      if (gen !== window.currentSearchGen() || !exons.length) return;
      renderCodonViewer(capturedGene, codonData, uniprotDomains || [], date, { queryCodon: null, exons });
    }).catch(() => {});
  }).catch(err => console.warn('[Codon Viewer] Error:', err));

  // Literature Review — only if phenotype present
  if (data.ptPhenoTexts.length > 0) {
    const pmPromise = runPhenoMatch(data.coords.gene);
    // When chapters + phenotype match are both ready, focus GeneReviews + ClinGen on the best-matching disease.
    Promise.all([grPromise, pmPromise]).then(() => {
      focusGeneReviewsByPhenotype(data.omimSynopsisResults, gen);
      refreshCurationSupport(gen);                                       // phenotype↔frequency match
    }).catch(() => {});
    Promise.all([gvPromise, pmPromise]).then(() => focusGeneValidityByPhenotype(data.omimSynopsisResults, gen)).catch(() => {});
  } else runPubMedOnly(data.coords.gene);

  fetchAssociatedConditions(data.coords.gene, gen).then(res => {
    if (gen !== window.currentSearchGen()) return;
    data.associatedConditions = res.conditions;
    maybeEvaluateBS2(gen);   // inheritance now known — refine the BS2 decision-support panel
    renderConditionCards(res.conditions, data.coords.gene, data.geneValidity?.curations || []);
    renderLitSearchLinks();
    _scheduleLitNotations(gen);
    // Conditions can resolve after populateSummaryRows already ran with an empty
    // list (it races runPhenoMatch). Re-populate the Variant+Disease row now.
    if (data.ptPhenoTexts.length > 0 && document.getElementById('sumVariantDiseaseList')) {
      populateSummaryRows(data.coords.gene);
    }
    const inhEl = document.getElementById('eInheritance');
    if (inhEl) inhEl.innerText = res.aggregatedInheritance?.trim() || '—';
  });
}

function triggerDownstreamAPIs(gen) {
  if (data.coords.gene) {
    // ── 1. PubMed / phenotype search — fired FIRST for fastest user-visible result ──
    let pmPromise = null;
    if (data.ptPhenoTexts.length > 0) pmPromise = runPhenoMatch(data.coords.gene);
    else runPubMedOnly(data.coords.gene);

    // ── 2. Variant-level API calls ────────────────────────────────────────────
    fetchBroadSpliceAI(gen);
    fetchGnomAD(gen);
    fetchGnomADConstraint(gen);
    fetchUCSC(gen);
    fetchDosageSensitivity(data.coords.chrom, data.coords.pos38, data.coords.gene);
    fetchClinVar(gen);
    checkVCEP(data.coords.gene);
    fetchLitVar(gen);
    const gvPromise = fetchGeneValidity(data.coords.gene).then(() => {
    renderGeneValidityCard();             // dedicated live card (own namespace)
    reconcileConditionCards(gen);         // overlay authoritative classification onto condition cards
    return enrichGeneValidityEvidence(data.geneValidity?.curations);
  }).then(() => renderGeneValidityCard()).catch(() => {});
    const grPromise = fetchGeneReviews(data.coords.gene, gen).catch(() => {});
    grPromise.then(() => refreshCurationSupport(gen)).catch(() => {});   // mechanism↔variant concordance
    // Age of onset for definitive curations (OMIM synopsis primary, GeneReviews fallback once chapters arrive).
    Promise.all([gvPromise, grPromise]).then(() =>
      enrichGeneValidityOnset(data.geneValidity?.curations, data.coords.gene, data.geneReviews?.chapters)
    ).then(() => { if (gen === window.currentSearchGen()) { renderGeneValidityCard(); maybeEvaluateBS2(gen); } }).catch(() => {});
    if (pmPromise) {
      Promise.all([grPromise, pmPromise]).then(() => {
        focusGeneReviewsByPhenotype(data.omimSynopsisResults, gen);
        refreshCurationSupport(gen);                                     // + phenotype↔frequency match
      }).catch(() => {});
      Promise.all([gvPromise, pmPromise]).then(() => focusGeneValidityByPhenotype(data.omimSynopsisResults, gen)).catch(() => {});
    }

    // Standard VEP fallback: if primary wasn't Ensembl, or if we need consequences
    if (data.coords.hgvs !== '-' && data.ensembl.transcript) {
      fetchEnsembl(data.ensembl.transcript + ':' + data.coords.hgvs, gen).then(() => {
        if (gen !== window.currentSearchGen()) return;
        // VEP may have resolved protein position — update viewer marker if it was missing.
        // Always patch _viewerArgs if pos is newly available; Phase B may have wiped it.
        const pos = Number(data.ensembl?.vepProteinStart) || null;
        if (pos && typeof _viewerArgs !== 'undefined' && _viewerArgs && _viewerArgs.queryCodon !== pos) {
          _viewerArgs.queryCodon = pos;
          _renderViewerInner();
        }
        // VEP has now set data.ensembl.vepConsequence — re-run mechanism↔variant concordance, which
        // may have run earlier (on grPromise) before the consequence was known (race). Idempotent.
        refreshCurationSupport(gen);
        // If Phase B's exon lookup previously failed (NM_→gene-expand timed out), retry
        // now that VEP has run. Pass the transcript (ENST or NM_) so Part A can resolve
        // NM_→ENST via the fast xref endpoint instead of gene-expand.
        const txAfterVep = data.ensembl?.transcript;
        if (txAfterVep && txAfterVep !== '-' && (!_viewerArgs?.exons?.length)) {
          fetchExonCodingPositions(txAfterVep, data.coords.gene).then(exons => {
            if (gen !== window.currentSearchGen() || !exons.length) return;
            if (typeof _viewerArgs !== 'undefined' && _viewerArgs) {
              _viewerArgs.exons = exons;
              _renderViewerInner();
            }
          }).catch(() => {});
        }
      }).catch(() => {});
    }

    // ── 3. Gene-level enrichment ──────────────────────────────────────────────
    fetchAssociatedConditions(data.coords.gene, gen).then(res => {
      if (gen !== window.currentSearchGen()) return;
      data.associatedConditions = res.conditions;
      maybeEvaluateBS2(gen);   // inheritance now known — refine the BS2 decision-support panel
      renderConditionCards(res.conditions, data.coords.gene, data.geneValidity?.curations || []);
      renderLitSearchLinks();
      // Re-trigger literature search now that disease names are available
      _scheduleLitNotations(gen);
      // Conditions can resolve after populateSummaryRows already ran with an empty
      // list (it races runPhenoMatch). Re-run only the Variant+Disease row —
      // not the full populateSummaryRows which would re-fire 3 extra NCBI calls.
      if (data.ptPhenoTexts.length > 0 && document.getElementById('sumVariantDiseaseList')) {
        populateVariantDiseaseRow(data.coords.gene);
      }
      const inhEl = document.getElementById('eInheritance');
      if (inhEl) {
        const inhText = res.aggregatedInheritance;
        inhEl.innerText = inhText && inhText.trim() ? inhText : '—';
      }
    });

    fetchClinVarGeneDistribution(data.coords.gene);
    fetchCivicData(gen);

    // ── 4. Protein variant viewer — two-phase progressive render ─────────────
    //
    // Phase A (fast): render as soon as the local codon file + UniProt domains
    //   are ready (~200–400 ms). This unblocks the viewer immediately.
    //
    // Phase B (progressive): re-render with Pfam structural domains + exon
    //   annotations once the slower Ensembl calls finish (~1–3 s). renderCodonViewer
    //   is idempotent — calling it again simply replaces the SVG in-place.
    // ─────────────────────────────────────────────────────────────────────────
    const domainPromise = fetchUniprotDomains(data.coords.gene, gen).then(domains => {
      if (gen !== window.currentSearchGen()) return [];
      data.uniprotDomains = domains;
      evaluateACMG(); // re-evaluate so PM1 is never missed
      return domains;
    }).catch(err => {
      console.warn('[Domain Annotation] Error:', err);
      data.uniprotDomains = [];
      return [];
    });

    // Phase A — render viewer immediately with fast data
    Promise.all([
      fetchCodonData(data.coords.gene),
      domainPromise,
    ]).then(([codonData, uniprotDomains]) => {
      if (gen !== window.currentSearchGen()) return;
      const date   = codonData?._meta?.date || null;
      const qCodon = Number(data.vepProteinStart) || Number(data.ensembl?.vepProteinStart) || null;
      renderCodonViewer(data.coords.gene, codonData, uniprotDomains || [], date,
        { queryCodon: qCodon, exons: [] });

      // Phase B — upgrade with Pfam + exons once Ensembl data arrives
      Promise.all([
        fetchExonCodingPositions(data.ensembl.transcript, data.coords.gene).catch(() => []),
        fetchEnsemblProteinDomains(data.ensembl.transcript, data.coords.gene, gen).catch(() => []),
      ]).then(([exons, pfamDomains]) => {
        if (gen !== window.currentSearchGen()) return;
        const allDomains = [...pfamDomains, ...(uniprotDomains || [])];
        // Only re-render if Ensembl actually returned something worth updating
        if (pfamDomains.length > 0 || exons.length > 0) {
          // Re-read protein position: VEP may have resolved it after Phase A ran (VV-primary path).
          // Also preserve any position already written into _viewerArgs by the VEP fallback handler.
          const finalQCodon = _viewerArgs?.queryCodon
            || Number(data.vepProteinStart)
            || Number(data.ensembl?.vepProteinStart)
            || qCodon;
          // If this Phase B exon lookup timed out (returned []) but a concurrent retry via
          // fetchEnsembl().then() already succeeded, preserve those exons rather than wiping them.
          const finalExons = exons.length > 0 ? exons : (_viewerArgs?.exons || []);
          renderCodonViewer(data.coords.gene, codonData, allDomains, date,
            { queryCodon: finalQCodon, exons: finalExons });
        }
      });
    });

    // ── 5. VCEP + MaveDB (depend on canonical allele resolution) ─────────────
    if (data.ensembl.transcript && data.coords.hgvs && data.coords.hgvs !== '-') {
      renderVCEPCard(); // Show initial 'Checking' state
      fetchCanonicalAllele(data.ensembl.transcript + ':' + data.coords.hgvs, gen).then(() => {
        if (gen !== window.currentSearchGen()) return;
        Promise.allSettled([
          fetchVCEPNarrative(gen),
          fetchMaveDBData(gen) // renderMaveDBCard called at end of fetchMaveDBData
        ]).then(() => evaluateACMG()); // Re-evaluate ACMG after MaveDB data is potentially fetched
      });
    }
  }
}

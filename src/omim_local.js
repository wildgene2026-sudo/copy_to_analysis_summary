// ── OMIM LOCAL DATA MODULE ───────────────────────────────────────────
// Loads OMIM synopsis cache + genemap2 at startup.
// Provides per-disease HPO matching against patient phenotype.
// Used as independent corroboration of the JAX HPO gene-level match.

const omimLocalCache = {
  synopses: new Map(),        // mimId (string) → Set<'HP:XXXXXXX'>
  synopsisTexts: new Map(),   // mimId (string) → lowercase plain text (for text-fallback matching)
  onset: new Map(),           // mimId (string) → age-of-onset clause(s) (cleaned, original casing)
  geneDiseaseMims: new Map(), // GENE_SYMBOL_UPPER → [{mimId, name, type, evidenceCode, moi}]
  geneCyto: new Map(),        // GENE_SYMBOL_UPPER → cytogenetic location string (e.g. "7q21.3")
  initialized: false,
  loading: false
};

// Promise that resolves once initOmimLocal() completes. Callers can
// await omimLocalReady() instead of polling omimLocalCache.initialized.
let _omimReadyResolve;
const _omimReadyPromise = new Promise(res => { _omimReadyResolve = res; });
function omimLocalReady() { return _omimReadyPromise; }

// ── SYNOPSIS HP EXTRACTION ───────────────────────────────────────────

/**
 * Extracts HP term IDs from an OMIM synopsis oldFormat string.
 * HP IDs appear inline in {…} cross-reference blocks, e.g.:
 *   "Short stature {UMLS C0349588 HP:0004322} {HPO HP:0004322}"
 * Lines starting with negation words are skipped to avoid false positives.
 */
function extractSynopsisHpTerms(text) {
  if (!text) return new Set();
  const hpSet = new Set();
  const NEGATION_PREFIX = /^(no |not |without |absent |absence of |lack of |negative for |excluded? )/i;

  for (const line of text.split('\n')) {
    const trimmed = line.trim();
    if (!trimmed || NEGATION_PREFIX.test(trimmed)) continue;

    // Fresh regex per line to avoid shared lastIndex state
    const blockRe = /\{([^}]+)\}/g;
    let m;
    while ((m = blockRe.exec(trimmed)) !== null) {
      const hpMatches = m[1].match(/HP:\d{7}/g);
      if (hpMatches) hpMatches.forEach(hp => hpSet.add(hp));
    }
  }
  return hpSet;
}

// ── AGE-OF-ONSET EXTRACTION ──────────────────────────────────────────

/**
 * Extracts the disease age-of-onset clause(s) from an OMIM synopsis source string.
 * OMIM records age of onset in the "Miscellaneous" body-system field (and, for older
 * entries, in the flat oldFormat blob), one item per ";"/newline-separated segment, e.g.:
 *   "Adult onset {UMLS C1853562 HP:0003581} {HPO HP:0003581}"
 *   "Onset first to seventh decade with 30 to 40 year mode"
 *   "Variable age at onset (childhood to adult)"
 * Returns a cleaned string (cross-reference {…} braces stripped) of up to two onset
 * clauses, the earliest "onset" mention first (the real disease onset usually leads the
 * segment, whereas incidental mentions like "early-onset form" sit deeper). null if none.
 */
function extractOmimOnset(text) {
  if (!text) return null;
  const clauses = [];
  for (const seg of String(text).split(/;|\n/)) {
    const s = seg.trim();
    if (!s) continue;
    const m = /onset/i.exec(s);
    if (!m) continue;
    const clean = s.replace(/\{[^}]*\}/g, '').replace(/\s+/g, ' ').trim();
    if (clean) clauses.push({ text: clean, pos: m.index });
  }
  if (!clauses.length) return null;
  clauses.sort((a, b) => a.pos - b.pos);
  const seen = new Set();
  const out = [];
  for (const c of clauses) {
    const k = c.text.toLowerCase();
    if (seen.has(k)) continue;
    seen.add(k);
    out.push(c.text);
    if (out.length >= 2) break;
  }
  return out.join('; ');
}

// ── GENEMAP2 PARSING ─────────────────────────────────────────────────

/**
 * Classifies a phenotype entry by its OMIM prefix marker:
 *   {}  → susceptibility/multifactorial
 *   []  → non-disease variation (excluded entirely)
 *   ?   → tentative/uncertain
 *   none→ confirmed disease
 */
function classifyPhenotypeEntry(raw) {
  const t = raw.trim();
  if (t.startsWith('[')) return 'non-disease';
  if (t.startsWith('{')) return 'susceptibility';
  if (t.startsWith('?')) return 'uncertain';
  return 'confirmed';
}

/**
 * Parses the genemap2 "Phenotypes" column (field 13).
 * Format: "Name, MIM (code), MOI; Name2, MIM2 (code2), MOI2"
 * Returns array of disease objects, excluding non-disease [] entries.
 */
function parseGenemap2Phenotypes(phenoStr) {
  if (!phenoStr) return [];
  const diseases = [];

  for (const entry of phenoStr.split(';')) {
    const raw = entry.trim();
    if (!raw) continue;

    const type = classifyPhenotypeEntry(raw);
    if (type === 'non-disease') continue;

    // Extract 6-digit MIM number and evidence code (1-4)
    const mimMatch = raw.match(/(\d{6})\s*\(([1-4])\)/);
    if (!mimMatch) continue;

    const mimId = mimMatch[1];
    const evidenceCode = parseInt(mimMatch[2]);

    // Disease name: everything before the MIM number, strip prefix markers
    let name = raw.substring(0, raw.indexOf(mimId)).trim();
    name = name.replace(/^[{[?]/, '').replace(/[}\]]$/, '').replace(/,\s*$/, '').trim();

    // MOI: text after "(code)", before next semicolon, strip leading comma
    const afterCode = raw.slice(raw.indexOf(`(${evidenceCode})`) + 3).trim();
    const moi = afterCode.replace(/^,\s*/, '').split(';')[0].trim() || null;

    diseases.push({ mimId, name, type, evidenceCode, moi });
  }
  return diseases;
}

// ── INITIALISATION ───────────────────────────────────────────────────

async function decompressGzip(buffer) {
  const ds = new DecompressionStream('gzip');
  const writer = ds.writable.getWriter();
  const reader = ds.readable.getReader();
  writer.write(new Uint8Array(buffer));
  writer.close();
  const chunks = [];
  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
  }
  let offset = 0, total = chunks.reduce((s, c) => s + c.length, 0);
  const out = new Uint8Array(total);
  chunks.forEach(c => { out.set(c, offset); offset += c.length; });
  return new TextDecoder().decode(out);
}

async function initOmimLocal() {
  if (omimLocalCache.initialized || omimLocalCache.loading) return;
  omimLocalCache.loading = true;

  try {
    // 1. Load clinical synopsis cache (gzipped JSON, ~4.7 MB)
    const cacheRes = await fetch('data/omim/omim_api_fallback_cache.json.gz');
    if (cacheRes.ok) {
      const buf = await cacheRes.arrayBuffer();
      // Dev servers (Vite) often send Content-Encoding: gzip, so the browser
      // already decompressed the body — manual DecompressionStream would then
      // fail on plain JSON. Detect gzip magic bytes (0x1f 0x8b) to decide.
      const bytes = new Uint8Array(buf);
      const isStillGzipped = bytes.length > 2 && bytes[0] === 0x1f && bytes[1] === 0x8b;
      const jsonText = isStillGzipped
        ? await decompressGzip(buf)
        : new TextDecoder().decode(bytes);
      const json = JSON.parse(jsonText);
      const synopsisData = json.clinical_synopsis?.clinical_synopses || {};
      // Non-clinical metadata fields to skip when joining structured entries
      const SKIP_FIELDS = new Set(['mimNumber', 'preferredTitle', 'prefix', 'molecularBasis']);

      let loaded = 0;
      for (const [mimId, entry] of Object.entries(synopsisData)) {
        // Two cache formats coexist:
        // 1. oldFormat — flat text blob (older entries, e.g. MIM 100050)
        // 2. Structured — body-system fields: growthHeight, skeletalLimbs, etc. (e.g. MIM 127300)
        let text = (entry.oldFormat || '').trim();
        if (!text) {
          // Join all non-metadata string field values into one searchable blob
          text = Object.entries(entry)
            .filter(([k, v]) => !SKIP_FIELDS.has(k) && typeof v === 'string' && v.trim())
            .map(([, v]) => v)
            .join('\n');
        }
        const hpSet = extractSynopsisHpTerms(text);
        // Always store text for fallback matching, even if no HP tags found
        omimLocalCache.synopsisTexts.set(mimId, text.toLowerCase());
        if (hpSet.size > 0) { omimLocalCache.synopses.set(mimId, hpSet); loaded++; }
        // Age of onset lives in the structured "miscellaneous" field (authoritative,
        // disease-level), else the flat oldFormat blob for legacy entries.
        const onset = extractOmimOnset(entry.miscellaneous || entry.oldFormat || '');
        if (onset) omimLocalCache.onset.set(mimId, onset);
      }
      console.log(`[OMIM Local] ${loaded} disease HP sets loaded from synopsis cache`);
    }

    // 2. Load genemap2.txt for gene→disease MIM mapping
    const gmRes = await fetch('data/omim/genemap2.txt');
    if (gmRes.ok) {
      const text = await gmRes.text();
      let loaded = 0;
      for (const line of text.split('\n')) {
        if (line.startsWith('#') || !line.trim()) continue;
        const f = line.split('\t');
        // genemap2 columns (0-indexed):
        // 0:chr 1:start 2:end 3:cyto 4:computedCyto 5:MIM 6:geneSymbols
        // 7:geneName 8:approvedSymbol 9:entrezId 10:ensemblId 11:comments
        // 12:phenotypes 13:mouseGene
        const approvedSymbol = f[8]?.trim();
        const cyto           = f[3]?.trim();   // cytogenetic location (e.g. "7q21.3")
        const phenotypesStr  = f[12]?.trim();
        if (!approvedSymbol || !phenotypesStr) continue;
        const diseases = parseGenemap2Phenotypes(phenotypesStr);
        if (diseases.length) {
          const key = approvedSymbol.toUpperCase();
          omimLocalCache.geneDiseaseMims.set(key, diseases);
          if (cyto) omimLocalCache.geneCyto.set(key, cyto);
          loaded++;
        }
      }
      console.log(`[OMIM Local] ${loaded} gene-disease mappings loaded from genemap2`);
    }

    omimLocalCache.initialized = true;
    console.log('[OMIM Local] Ready.');
    _omimReadyResolve(); // unblock any callers awaiting omimLocalReady()
  } catch (e) {
    console.warn('[OMIM Local] Init failed:', e.message);
    _omimReadyResolve(); // resolve even on failure so callers are never permanently blocked
  } finally {
    omimLocalCache.loading = false;
  }
}

// ── MATCHING LOGIC ───────────────────────────────────────────────────

/**
 * Matches patient HPO groups against a disease's synopsis HP set + text.
 * Three-tier matching:
 *   Tier 1a — exact HP ID: patient's resolved HP IDs intersect synopsis HP IDs (weight 1).
 *   Tier 1b — IC-floored ontology overlap (only when `sim` is provided): bidirectional —
 *            a patient term and a synopsis term that share an information-rich common
 *            ancestor (background p ≤ 5%) match at a Lin-similarity-discounted weight in
 *            (0,1). Closes the parent↔child gap (e.g. patient "Ataxia" vs synopsis
 *            "Cerebellar ataxia") WITHOUT pretending the specificity is confirmed, and the
 *            IC floor stops spurious matches on shallow ancestors ("Abnormality of the
 *            nervous system"). `sim` = window.HpoFit once its tables have loaded.
 *   Tier 2 — Text fallback: patient term plain-text substring in synopsis text (weight 1).
 *
 * Returns { matched, matchScore, recall, precision, fScore, matchedIndices, matchTypes, matchWeights }
 *   matched      = integer count of matched patient terms (for "X/total" display)
 *   matchScore   = Σ weights (fractional; drives recall/precision/fScore)
 *   matchTypes[i]= 'hp' | 'ic' | 'text' for each matched index
 */
function matchPatientToDisease(ptPhenoGroups, ptPhenoTexts, mimId, diseaseHpSet, sim = null) {
  const hasHpSet  = diseaseHpSet?.size > 0;
  const synText   = omimLocalCache.synopsisTexts.get(mimId) || '';
  const hasText   = synText.length > 0;
  const empty = { matched: 0, matchScore: 0, recall: 0, precision: 0, fScore: 0, matchedIndices: [], matchTypes: [], matchWeights: [] };

  if (!hasHpSet && !hasText) return { ...empty };
  if (!ptPhenoGroups?.length) return { ...empty };

  const useIc = !!(sim && sim.ready && hasHpSet && typeof sim.linSim === 'function');
  const diseaseIds = useIc ? [...diseaseHpSet] : null;

  let matched = 0;       // integer count of matched patient terms
  let matchScore = 0;    // Σ weights (fractional)
  const matchedIndices = [];
  const matchTypes     = [];
  const matchWeights   = [];

  for (let i = 0; i < ptPhenoGroups.length; i++) {
    const ids = ptPhenoGroups[i] || [];

    // Tier 1a: exact HP ID match
    if (hasHpSet && ids.some(id => diseaseHpSet.has(id))) {
      matched++; matchScore += 1;
      matchedIndices.push(i); matchTypes.push('hp'); matchWeights.push(1);
      continue;
    }

    // Tier 1b: IC-floored bidirectional ontology overlap (discounted partial match)
    if (useIc && ids.length) {
      let best = 0;
      for (const pid of ids) {
        for (const did of diseaseIds) {
          const w = sim.linSim(pid, did);
          if (w > best) best = w;
          if (best >= 0.999) break;
        }
        if (best >= 0.999) break;
      }
      if (best > 0) {
        matched++; matchScore += best;
        matchedIndices.push(i); matchTypes.push('ic'); matchWeights.push(best);
        continue;
      }
    }

    // Tier 2: text fallback — patient term substring in synopsis text
    if (hasText && ptPhenoTexts?.[i]) {
      const needle = ptPhenoTexts[i].toLowerCase().trim();
      // Require at least 4 chars to avoid noise from very short terms
      if (needle.length >= 4 && synText.includes(needle)) {
        matched++; matchScore += 1;
        matchedIndices.push(i); matchTypes.push('text'); matchWeights.push(1);
      }
    }
  }

  // recall/precision use the (fractional) matchScore so a discounted ontology hit
  // contributes proportionally; `matched` stays an integer for display.
  const recall    = matchScore / ptPhenoGroups.length;
  const precision = hasHpSet ? matchScore / diseaseHpSet.size : (matchScore > 0 ? 0.5 : 0);
  const fScore    = (recall + precision > 0)
    ? (2 * recall * precision) / (recall + precision)
    : 0;
  return { matched, matchScore, recall, precision, fScore, matchedIndices, matchTypes, matchWeights };
}

/**
 * Main entry point.
 * Returns array of disease match results sorted by fScore desc.
 * Each result has: mimId, name, type, evidenceCode, moi, matched, recall,
 *   precision, fScore, hasSynopsis, isSusceptibility, isUncertain,
 *   matchedTerms: [{text, discriminating}]
 */
function matchOmimSynopsis(geneSymbol, ptPhenoGroups, ptPhenoTexts, sim = (typeof window !== 'undefined' ? window.HpoFit : null)) {
  if (!omimLocalCache.initialized || !geneSymbol || !ptPhenoGroups?.length) return [];

  const diseases = omimLocalCache.geneDiseaseMims.get(geneSymbol.toUpperCase()) || [];

  // Only process confirmed (code 3+) and susceptibility diseases
  // Skip non-disease [] and low-evidence (1,2) entries
  const candidates = diseases.filter(d => d.type !== 'non-disease' && d.evidenceCode >= 3);
  if (!candidates.length) return [];

  // Step 1: match each disease independently
  const results = candidates.map(disease => {
    const synopsisHpSet = omimLocalCache.synopses.get(disease.mimId);
    const hasText       = omimLocalCache.synopsisTexts.has(disease.mimId);
    const { matched, recall, precision, fScore, matchedIndices, matchTypes } =
      matchPatientToDisease(ptPhenoGroups, ptPhenoTexts, disease.mimId, synopsisHpSet, sim);
    return {
      ...disease,
      synopsisHpCount: synopsisHpSet?.size || 0,
      matched,
      recall,
      precision,
      fScore,
      matchedIndices,
      matchTypes,
      matchedTerms: [],       // filled in Step 2
      hasSynopsis: !!synopsisHpSet || hasText,
      isSusceptibility: disease.type === 'susceptibility',
      isUncertain: disease.type === 'uncertain'
    };
  });

  // Step 2: compute discriminating terms
  // A term is "discriminating" if it matches only 1 confirmed disease for this gene.
  // A term is "shared" if it matches all confirmed diseases.
  const confirmedCount = results.filter(r => !r.isSusceptibility).length;
  const termMatchCount = new Map(); // ptPhenoTexts[i] → count of diseases that match it

  for (const r of results) {
    if (r.isSusceptibility) continue;
    for (const idx of r.matchedIndices) {
      const t = ptPhenoTexts[idx];
      termMatchCount.set(t, (termMatchCount.get(t) || 0) + 1);
    }
  }

  for (const r of results) {
    r.matchedTerms = r.matchedIndices.map((idx, pos) => {
      const text  = ptPhenoTexts[idx];
      const count = termMatchCount.get(text) || 1;
      return {
        text,
        matchType:      r.matchTypes?.[pos] || 'hp',  // 'hp' or 'text'
        discriminating: count === 1,
        shared:         count === confirmedCount && confirmedCount > 1
      };
    });
  }

  // Step 3: sort — fScore desc, then evidenceCode desc, then hasSynopsis
  results.sort((a, b) => {
    if (b.fScore !== a.fScore) return b.fScore - a.fScore;
    if (b.evidenceCode !== a.evidenceCode) return b.evidenceCode - a.evidenceCode;
    return (b.hasSynopsis ? 1 : 0) - (a.hasSynopsis ? 1 : 0);
  });

  return results;
}

/**
 * Produces the top-level summary for the Phenotype & Gene Match section.
 * Only considers confirmed (non-susceptibility) diseases.
 */
function getOmimSynopsisSummary(results) {
  const confirmed = (results || []).filter(r => !r.isSusceptibility && r.fScore > 0 && r.matched > 0);
  if (!confirmed.length) return null;

  const topScore   = confirmed[0].fScore;
  const topGroup   = confirmed.filter(r => r.fScore === topScore);
  const pct        = Math.round(topScore * 100);

  if (topGroup.length === 1) {
    const d = topGroup[0];
    const termList = d.matchedTerms.map(t => (t.discriminating ? '★ ' : '') + t.text).join(' · ');
    return {
      type: 'single',
      label: `Best synopsis match: ${d.name}`,
      detail: `F:${pct}% · ${d.matched}/${d.matched} terms · ${termList}`,
      diseases: topGroup
    };
  }
  if (topGroup.length <= 3) {
    return {
      type: 'multiple',
      label: `Top synopsis matches: ${topGroup.map(d => d.name).join(', ')}`,
      detail: `F:${pct}% each — additional phenotype terms needed to discriminate`,
      diseases: topGroup
    };
  }
  return {
    type: 'ambiguous',
    label: `No discriminating synopsis match — ${topGroup.length} diseases score equally`,
    detail: `F:${pct}% · phenotype consistent with ${topGroup[0].name.split(' ')[0]}… gene spectrum`,
    diseases: topGroup
  };
}

/**
 * getOmimGeneTable(geneSymbol)
 * Inputs: omimLocalCache.geneDiseaseMims + geneCyto (local genemap2 data)
 * Outputs: array of rows for the OMIM phenotype table, mirroring the OMIM
 *   gene page layout: { location, phenotype, mimNumber, inheritance, mappingKey }
 * Failures: returns [] if cache not ready or gene not found.
 */
function getOmimGeneTable(geneSymbol) {
  if (!omimLocalCache.initialized || !geneSymbol) return [];
  const key = geneSymbol.toUpperCase();
  const diseases = omimLocalCache.geneDiseaseMims.get(key) || [];
  if (!diseases.length) return [];
  const location = omimLocalCache.geneCyto.get(key) || '';

  // Abbreviate full MOI text to OMIM-style codes (AD, AR, XL, etc.)
  const abbrevMoi = (moi) => {
    if (!moi) return '';
    const m = moi.toUpperCase();
    if (m.includes('X-LINKED DOMINANT')) return 'XLD';
    if (m.includes('X-LINKED RECESSIVE')) return 'XLR';
    if (m.includes('X-LINKED')) return 'XL';
    if (m.includes('Y-LINKED')) return 'YL';
    if (m.includes('SEMIDOMINANT') || m.includes('SEMI-DOMINANT')) return 'SD';
    // Pseudoautosomal MUST be tested before autosomal: "PSEUDOAUTOSOMAL DOMINANT"
    // contains the substring "AUTOSOMAL DOMINANT" (e.g. SHOX → PD/PR, not AD/AR).
    if (m.includes('PSEUDOAUTOSOMAL DOMINANT')) return 'PD';
    if (m.includes('PSEUDOAUTOSOMAL RECESSIVE')) return 'PR';
    if (m.includes('AUTOSOMAL DOMINANT')) return 'AD';
    if (m.includes('AUTOSOMAL RECESSIVE')) return 'AR';
    if (m.includes('MITOCHONDRIAL')) return 'Mit';
    if (m.includes('DIGENIC')) return 'DD';
    if (m.includes('SOMATIC')) return 'SMu';
    if (m.includes('MULTIFACTORIAL')) return 'Mu';
    return moi;
  };

  // OMIM phenotype name formatting: braces for susceptibility, '?' for uncertain.
  // Strip any stray markers left on the parsed name before re-wrapping.
  const formatName = (d) => {
    const clean = (d.name || '').replace(/^[{[?]+/, '').replace(/[}\]]+$/, '').trim();
    if (d.type === 'susceptibility') return `{${clean}}`;
    if (d.type === 'uncertain') return `?${clean}`;
    return clean;
  };

  return diseases.map(d => ({
    location,
    phenotype: formatName(d),
    mimNumber: d.mimId,
    inheritance: abbrevMoi(d.moi),
    mappingKey: d.evidenceCode   // OMIM phenotype mapping key (1-4)
  }));
}

/**
 * getOmimOnset(mimId)
 * Inputs: omimLocalCache.onset (local synopsis-derived age-of-onset map).
 * Outputs: cleaned age-of-onset clause string for the disease, or null if none recorded.
 * Failures: returns null if cache not ready or the MIM has no onset clause.
 */
function getOmimOnset(mimId) {
  if (!omimLocalCache.initialized || !mimId) return null;
  return omimLocalCache.onset.get(String(mimId)) || null;
}

// ── AUTO-INIT ON PAGE LOAD ───────────────────────────────────────────
(function () {
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initOmimLocal);
  } else {
    initOmimLocal();
  }
})();

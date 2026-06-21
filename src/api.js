// ── DATA FETCHING & API LOGIC ──────────────────────────────────────

/**
 * Robust fetch helper with exponential backoff and timeout.
 * @param {string} url - The endpoint to fetch.
 * @param {object} options - Fetch options (method, headers, etc).
 * @param {number} maxRetries - Maximum number of retry attempts.
 * @returns {Promise<Response>}
 */

// Override native fetch to automatically bind the global abort signal
const originalFetch = window.fetch;

const VEP_SEVERITY = {
  'transcript_ablation': 1, 'splice_acceptor_variant': 2, 'splice_donor_variant': 3,
  'stop_gained': 4, 'frameshift_variant': 5, 'stop_lost': 6, 'start_lost': 7,
  'transcript_amplification': 8, 'inframe_insertion': 9, 'inframe_deletion': 10,
  'missense_variant': 11, 'protein_altering_variant': 12, 'splice_region_variant': 13,
  'incomplete_terminal_codon_variant': 14, 'start_retained_variant': 15,
  'stop_retained_variant': 16, 'synonymous_variant': 17, 'coding_sequence_variant': 18,
  'mature_miRNA_variant': 19, '5_prime_UTR_variant': 20, '3_prime_UTR_variant': 21,
  'non_coding_transcript_exon_variant': 22, 'intron_variant': 23,
  'NMD_transcript_variant': 24, 'non_coding_transcript_variant': 25,
  'upstream_gene_variant': 26, 'downstream_gene_variant': 27, 'TFBS_ablation': 28,
  'TFBS_amplification': 29, 'TF_binding_site_variant': 30,
  'regulatory_region_ablation': 31, 'regulatory_region_amplification': 32,
  'feature_elongation': 33, 'regulatory_region_variant': 34,
  'feature_truncation': 35, 'intergenic_variant': 36
};

function getSeverity(terms) {
  if (!terms || !terms.length) return 99;
  return Math.min(...terms.map(t => VEP_SEVERITY[t] || 99));
}

/**
 * Safely parses a value into a float, returning null if invalid.
 */
const parseSafe = (val) => {
  const f = parseFloat(val);
  return isNaN(f) ? null : f;
};

/**
 * Extracts the best REVEL score for a specific MANE transcript.
 * Robustly handles dedicated 'revel' tracks and dbNSFP's nested arrays/strings.
 */
function getBestRevel(h, maneTx = '') {
  const candidates = [];
  const sources = [
    { scores: h.revel?.revel_score, ids: h.revel?.refseq_id },
    { scores: h.dbnsfp?.revel_score, ids: h.dbnsfp?.refseq_id },
    { scores: h.dbnsfp?.revel?.score, ids: h.dbnsfp?.refseq_id || h.dbnsfp?.refseq_id_rev }
  ];

  const cleanMane = maneTx ? String(maneTx).split('.')[0] : '';

  for (const src of sources) {
    if (!src.scores) continue;
    const res = getBestDbnsfpVal(src.scores, src.ids, maneTx);
    if (res !== null) candidates.push(res);
  }

  return candidates.length ? Math.max(...candidates) : null;
}

/**
 * Extracts the best AlphaMissense score for a specific MANE transcript.
 * Built to safely handle dbNSFP's unpredictable arrays, sub-arrays, and flat strings.
 */
function getBestAlphaMissense(h, maneTx = '') {
  if (!h || !h.dbnsfp) return null;
  const db = h.dbnsfp;
  const scoreData = db.alphamissense_score ?? db.alphamissense?.score;
  const refseqData = db.refseq_id;
  const cleanMane = maneTx ? String(maneTx).split('.')[0] : '';
  return getBestDbnsfpVal(scoreData, refseqData, cleanMane);
}

/**
 * Extracts the AlphaMissense prediction label (categorical) for a specific transcript.
 * Prioritizes exact transcript version matches.
 */
function getBestAlphaMissensePred(h, maneTx = '') {
  if (!h || !h.dbnsfp) return null;
  const db = h.dbnsfp;
  const predData = db.alphamissense_pred ?? db.alphamissense?.pred;
  const refseqData = db.refseq_id;
  const fullMane = maneTx ? String(maneTx).trim() : '';
  const cleanMane = fullMane.split('.')[0];

  if (!predData) return null;

  if (Array.isArray(predData)) {
    const idArr = Array.isArray(refseqData) ? refseqData : [refseqData];
    let exactMatch = null;
    let prefixMatch = null;

    for (let i = 0; i < predData.length; i++) {
      if (idArr[i]) {
        const idsAtIndex = Array.isArray(idArr[i]) ? idArr[i] : String(idArr[i]).split(',');
        if (fullMane && idsAtIndex.some(id => String(id).trim() === fullMane)) {
          exactMatch = predData[i]; break;
        }
        if (cleanMane && !prefixMatch && idsAtIndex.some(id => String(id).split('.')[0] === cleanMane)) {
          prefixMatch = predData[i];
        }
      }
    }
    return exactMatch || prefixMatch || predData[0];
  }
  return predData;
}

/**
 * Generic dbNSFP value extractor that handles nested arrays, comma-strings, and MANE matching.
 */
function getBestDbnsfpVal(valData, refseqData, maneTx) {
  if (valData === undefined || valData === null) return null;

  const fullMane = maneTx ? String(maneTx).trim() : '';
  const cleanMane = fullMane.split('.')[0];

  let bestExactScore = null;
  let bestPrefixScore = null;
  let highestOverallScore = null;

  if (Array.isArray(valData)) {
    const idArr = Array.isArray(refseqData) ? refseqData : [refseqData];
    for (let i = 0; i < valData.length; i++) {
      const val = parseSafe(valData[i]);
      if (val === null) continue;
      if (highestOverallScore === null || val > highestOverallScore) highestOverallScore = val;

      if (idArr[i]) {
        const idsAtIndex = Array.isArray(idArr[i]) ? idArr[i] : String(idArr[i]).split(',');

        // Priority 1: Exact Version Match (e.g. NM_020442.6)
        if (fullMane && idsAtIndex.some(id => String(id).trim() === fullMane)) {
          if (bestExactScore === null || val > bestExactScore) bestExactScore = val;
        }
        // Priority 2: Prefix Match (e.g. NM_020442)
        if (cleanMane && idsAtIndex.some(id => String(id).split('.')[0] === cleanMane)) {
          if (bestPrefixScore === null || val > bestPrefixScore) bestPrefixScore = val;
        }
      }
    }
  } else {
    const val = parseSafe(valData);
    if (val !== null) {
      highestOverallScore = val;
      if (refseqData) {
        const idArr = Array.isArray(refseqData) ? refseqData.flat(Infinity) : String(refseqData).split(',');
        if (fullMane && idArr.some(id => String(id).trim() === fullMane)) bestExactScore = val;
        if (cleanMane && idArr.some(id => String(id).split('.')[0] === cleanMane)) bestPrefixScore = val;
      }
      // If no transcript IDs are available, return the score as a fallback since the coordinate matched
      if (bestExactScore === null && bestPrefixScore === null) {
        highestOverallScore = val;
      }
    }
  }

  return bestExactScore || bestPrefixScore || highestOverallScore;
}

function normalizeVcfAllele(a) {
  if (a == null || a === '' || a === '-') return '';
  return String(a);
}

/**
 * Converts HGVS-style indel alleles (one allele is '-'/empty, as returned by
 * ClinGen Allele Registry and Ensembl VEP) to proper VCF-style coordinates by
 * fetching the anchor base from the reference genome.
 *
 * VCF format requires an anchor base (last unchanged base before the indel):
 *   Insertion: ref = anchor,          alt = anchor + inserted_seq, pos = anchor_pos
 *   Deletion:  ref = anchor + del_seq, alt = anchor,               pos = anchor_pos
 *
 * Anchor position (1-based):
 *   - ClinGen AR:   anchorPos = c3x.start  (0-based AR start equals 1-based VCF anchor)
 *   - Ensembl VEP:  anchorPos = main.start - 1
 *
 * Strategy: try Ensembl first (8s timeout); fall back to NCBI efetch on any
 * failure. If both fail, returns null and the caller shows a UI warning.
 *
 * Why originalFetch + explicit AbortController (not fetchWithRetry):
 *   fetchWithRetry injects globalSearchController.signal, which is aborted on
 *   every new user input. A mid-flight anchor lookup must survive a new search
 *   being typed, otherwise the hung-connection pool bug reappears.
 *
 * @param {string}        chr       Chromosome without 'chr' prefix
 * @param {number|string} anchorPos 1-based position of the anchor base
 * @param {string}        hgvsRef   '' for insertions; deleted bases for deletions
 * @param {string}        hgvsAlt   inserted bases for insertions; '' for deletions
 * @returns {Promise<{pos:string, ref:string, alt:string}|null>}
 */
async function hgvsIndelToVcf(chr, anchorPos, hgvsRef, hgvsAlt) {
  const cleanChr = String(chr).replace(/^chr/i, '');
  const aPos = parseInt(anchorPos);
  if (!aPos || aPos < 1) return null;

  // Build VCF coords from a resolved anchor base.
  const buildResult = anchor => {
    if (!anchor) return null;
    if (!hgvsRef) return { pos: String(aPos), ref: anchor, alt: anchor + hgvsAlt.toUpperCase() };
    return { pos: String(aPos), ref: anchor + hgvsRef.toUpperCase(), alt: anchor };
  };

  // ── Primary: Ensembl sequence endpoint (8 s timeout) ─────────────────
  try {
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 8000);
    try {
      const res = await originalFetch(
        `https://rest.ensembl.org/sequence/region/human/${cleanChr}:${aPos}..${aPos}:1?content-type=application/json`,
        { signal: controller.signal }
      );
      clearTimeout(timeoutId);
      if (res.ok) {
        const json = await res.json();
        const result = buildResult((json.seq || '').toUpperCase().charAt(0));
        if (result) return result;
      }
    } finally {
      clearTimeout(timeoutId);
    }
  } catch (e) {
    console.warn('[hgvsIndelToVcf] Ensembl anchor fetch failed:', e.message);
  }

  // ── Fallback: NCBI efetch (uses existing token-bucket rate limiter) ───
  try {
    const GRCh38_ACCESSION = {
      '1':'NC_000001.11','2':'NC_000002.12','3':'NC_000003.12','4':'NC_000004.12',
      '5':'NC_000005.10','6':'NC_000006.12','7':'NC_000007.14','8':'NC_000008.11',
      '9':'NC_000009.12','10':'NC_000010.11','11':'NC_000011.10','12':'NC_000012.12',
      '13':'NC_000013.11','14':'NC_000014.9','15':'NC_000015.10','16':'NC_000016.10',
      '17':'NC_000017.11','18':'NC_000018.10','19':'NC_000019.10','20':'NC_000020.11',
      '21':'NC_000021.9','22':'NC_000022.11','X':'NC_000023.11','Y':'NC_000024.10','MT':'NC_012920.1'
    };
    const accession = GRCh38_ACCESSION[cleanChr.toUpperCase()];
    if (accession) {
      await ncbiAcquireToken();
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), 10000);
      try {
        const url = `${NCBI_BASE}/efetch.fcgi?db=nuccore&id=${accession}&seq_start=${aPos}&seq_stop=${aPos}&rettype=fasta&retmode=text&api_key=${NCBI_KEY}`;
        const res = await originalFetch(url, { signal: controller.signal });
        clearTimeout(timeoutId);
        if (res.ok) {
          const text = await res.text();
          // FASTA: first line is ">header", second line is the sequence base
          const base = text.split('\n').find(l => l && !l.startsWith('>'))?.trim().toUpperCase().charAt(0);
          const result = buildResult(base);
          if (result) return result;
        }
      } finally {
        clearTimeout(timeoutId);
      }
    }
  } catch (e) {
    console.warn('[hgvsIndelToVcf] NCBI anchor fetch failed:', e.message);
  }

  // Both APIs failed — caller must handle null and warn the user.
  return null;
}

window.fetch = async function (url, options = {}) {
  if (!options.signal && window.globalSearchController?.signal) {
    options.signal = window.globalSearchController.signal;
  }
  return originalFetch.call(window, url, options);
};

// Wraps originalFetch with a hard abort timeout so stalled remote connections
// are released within timeoutMs. Prevents Windows TCP pool exhaustion when an
// API server stops responding (6 hung connections fill the per-host pool).
// Rule: use fetchWithTimeout for ALL remote originalFetch calls; bare
// originalFetch is reserved for cases that manage their own AbortController.
function fetchWithTimeout(url, options = {}, timeoutMs = 10000) {
  const controller = new AbortController();
  const id = setTimeout(() => controller.abort(), timeoutMs);
  return originalFetch(url, { ...options, signal: controller.signal })
    .finally(() => clearTimeout(id));
}

const responseCache = new Map();
const RESPONSE_CACHE_MAX = 200;

// ── NCBI RATE LIMITER (token bucket) ──────────────────────────────────
// NCBI E-utilities allow 10 req/s with an API key. We cap at 9/s with a
// burst capacity of 9 so short bursts pass instantly while sustained load
// is throttled below the limit — eliminating reactive 429 retries.
// NCBI responses occasionally contain raw control characters (0x00-0x1F) inside
// JSON string values — e.g. variant titles, webenv tokens, or gene descriptions.
// RFC 8259 forbids unescaped control characters in JSON strings; browsers throw
// "Bad control character in string literal" on JSON.parse. This helper strips
// the offending bytes before parsing.
// Keeps: 0x09 TAB, 0x0A LF, 0x0D CR — valid JSON structural whitespace.
// Strips: 0x00-0x08, 0x0B VT, 0x0C FF, 0x0E-0x1F, 0x7F DEL.
async function ncbiSafeJson(res) {
  const text = await res.text();
  const clean = text.replace(/[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]/g, '');
  return JSON.parse(clean);
}

const _ncbiBucket = { tokens: 9, max: 9, refillPerMs: 9 / 1000, last: Date.now() };
async function ncbiAcquireToken() {
  while (true) {
    const now = Date.now();
    _ncbiBucket.tokens = Math.min(
      _ncbiBucket.max,
      _ncbiBucket.tokens + (now - _ncbiBucket.last) * _ncbiBucket.refillPerMs
    );
    _ncbiBucket.last = now;
    if (_ncbiBucket.tokens >= 1) { _ncbiBucket.tokens -= 1; return; } // check+decrement is atomic (no await between)
    const waitMs = Math.ceil((1 - _ncbiBucket.tokens) / _ncbiBucket.refillPerMs);
    await new Promise(r => setTimeout(r, waitMs));
  }
}

// ── IN-FLIGHT REQUEST DEDUP ───────────────────────────────────────────
// Concurrent GETs for the same URL share a single network request. Each
// caller receives an independent clone, so body reads never collide.
const _pendingGets = new Map(); // url -> Promise<Response>


async function fetchWithRetry(url, options = {}, maxRetries = 3) {
  if (!options.signal && window.globalSearchController?.signal) {
    options.signal = window.globalSearchController.signal;
  }

  const isGet = options.method !== 'POST';

  // Cache hit — return immediately (POSTs skip cache to avoid query collisions).
  if (isGet && responseCache.has(url)) {
    return new Response(responseCache.get(url), { status: 200, headers: { 'Content-Type': 'application/json' } });
  }

  // In-flight dedup — concurrent identical GETs share one network request.
  // Each caller gets a clone so body reads never collide.
  if (isGet && _pendingGets.has(url)) {
    return (await _pendingGets.get(url)).clone();
  }

  const work = _fetchWithRetryInner(url, options, maxRetries);
  if (isGet) {
    _pendingGets.set(url, work);
    try {
      const resp = await work;
      return resp.clone();
    } finally {
      _pendingGets.delete(url);
    }
  }
  return work;
}

async function _fetchWithRetryInner(url, options, maxRetries) {
  const targetUrl = url;
  const isNcbi = url.includes('eutils.ncbi.nlm.nih.gov');

  let lastError;
  for (let i = 0; i < maxRetries; i++) {
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(new Error('Timeout')), 15000);

    if (options.signal) {
      options.signal.addEventListener('abort', () => {
        controller.abort(options.signal.reason);
      }, { once: true });
    }

    try {
      if (i > 0) {
        const delay = i * 1000;
        document.getElementById('statusMsg').innerText = `VV Busy... Retrying (Attempt ${i + 1}/${maxRetries})...`;
        await new Promise(res => setTimeout(res, delay));
      }

      // Proactive NCBI throttle: stay under the 10 req/s key limit.
      if (isNcbi) await ncbiAcquireToken();

      const fetchOpts = { ...options, signal: controller.signal };
      const response = await fetch(targetUrl, fetchOpts);
      clearTimeout(timeoutId);

      if (response.ok) {
        const cloned = response.clone();
        cloned.text().then(text => {
          responseCache.set(url, text); // Cache the original URL
          if (responseCache.size > RESPONSE_CACHE_MAX) {
            responseCache.delete(responseCache.keys().next().value);
          }
        }).catch(e => console.warn('Cache write failed:', e));
        return response;
      }

      // If we get here, it's a non-2xx response
      lastError = new Error(`HTTP ${response.status}: ${response.statusText}`);
      lastError.status = response.status;

      // Don't retry on certain errors (e.g., 404, 400)
      if (response.status >= 400 && response.status < 500 && response.status !== 429) {
        throw lastError;
      }
    } catch (err) {
      clearTimeout(timeoutId);
      lastError = err;
      if (err.name === 'AbortError') {
        if (options.signal && options.signal.aborted) {
          throw err;
        }
        lastError = new Error('Request timed out (15s)');
      }
      console.warn(`Attempt ${i + 1} failed: ${lastError.message}`);
    }
  }
  throw lastError;
}

async function fetchVV(query, gen) {
  setDot('vv', 'loading');
  document.getElementById('statusMsg').innerText = 'Querying VariantValidator...';

  try {
    // 1. Sanitize: Handle clinical strings (e.g., "Transcript(Gene):c.Coord p.Protein")
    let sanitized = query.trim();

    // a. Remove embedded gene names in parentheses (e.g., "NM_019842.4(KCNQ5):c.2228del")
    sanitized = sanitized.replace(/\(([^)]+)\):/, ':');

    // b. Remove protein consequences (both spaced and concatenated)
    // Supports "c.2228del p.(Ala...) " and "c.2228delp.Ala..."
    if (sanitized.includes('c.')) {
      sanitized = sanitized.split(/\s*\(?p\./)[0].trim();
    }

    // 2. Encode for URL safety
    const safeQuery = encodeURIComponent(sanitized);

    const url = `https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/${safeQuery}/mane_select`;

    // 3. Fetch with Resilient Logic
    const res = await fetchWithRetry(url);
    const json = await res.json();
    if (gen !== window.currentSearchGen()) return;

    let found = false;
    for (let key in json) {
      const rec = json[key];
      if (rec.gene_symbol) {
        const geneStr = (rec.gene_symbol || '').trim().toUpperCase();
        data.coords.gene = geneStr;
        document.getElementById('dGene').innerText = geneStr;

        data.ensembl.transcript = (rec.hgvs_transcript_variant || '').split(':')[0];
        document.getElementById('eManeTranscript').innerText = data.ensembl.transcript || '-';
        data.coords.hgvs = rec.hgvs_transcript_variant ? (rec.hgvs_transcript_variant.split(':')[1] || '') : '-';
        document.getElementById('dHgvs').innerText = data.coords.hgvs;
        if (rec.hgvs_predicted_protein_consequence?.tlr) {
          data.ensembl.protein = rec.hgvs_predicted_protein_consequence.tlr;
          document.getElementById('dProtein').innerText = data.ensembl.protein;
        }
        const g38 = rec.primary_assembly_loci?.grch38 || rec.primary_assembly_loci?.hg38;
        if (g38) {
          if (g38.vcf) {
            const v = g38.vcf;
            data.coords.chrom = v.chr.replace('chr', '');
            data.coords.pos38 = String(v.pos);
            data.coords.ref = normalizeVcfAllele(v.ref);
            data.coords.alt = normalizeVcfAllele(v.alt);
            if (data.coords.ref && data.coords.alt) {
              data.coords.hg38String = `chr${data.coords.chrom}-${data.coords.pos38}-${data.coords.ref}-${data.coords.alt}`;
              document.getElementById('dHg38').innerText = data.coords.hg38String;
            }
          }
        }
        const g19 = rec.primary_assembly_loci?.hg19 || rec.primary_assembly_loci?.grch37;
        if (g19?.vcf) {
          const v = g19.vcf;
          data.coords.pos19 = v.pos;
          const ref19 = normalizeVcfAllele(v.ref);
          const alt19 = normalizeVcfAllele(v.alt);
          if (ref19 && alt19) {
            data.coords.hg19String = `chr${v.chr.replace('chr', '')}-${v.pos}-${ref19}-${alt19}`;
            data.coords.hg19FromVV = true;
            document.getElementById('dHg19').innerText = data.coords.hg19String;
          }
        }
        found = true; break;
      }
    }

    if (found) {
      // Cross-validate VV coords with ClinGen Allele Registry (mutual validation: both are ground truth sources)
      const hgvsForAR = (data.ensembl.transcript && data.coords.hgvs && data.coords.hgvs !== '-')
        ? `${data.ensembl.transcript}:${data.coords.hgvs}` : sanitized;
      await confirmCoordsWithClinGenAR(hgvsForAR)
        .catch(e => console.warn('[VV→ClinGen AR] validation skipped:', e.message));

      setDot('vv', 'ok'); document.getElementById('statusMsg').innerText = 'Ready.';
      runValidationPass();
      enableWrappers(['franklin', 'gnomadv4', 'spliceai', 'liftover', 'clingen', 'clinvar', 'omim', 'gr', 'scholar', 'decipher', 'hgmd', 'gtex', 'mastermind', 'gnomadv2']);
      document.getElementById('btnLaunchAll').disabled = false;
      document.getElementById('btnLaunchSel').disabled = false;
      document.getElementById('geneInput').value = data.coords.gene;
      data.coords.primaryFromVV = true;
      // All secondary background fetches moved to app.js decideAndFetch for strict sequencing
    } else {
      setDot('vv', 'warn');
      document.getElementById('statusMsg').innerText = 'Not found in VariantValidator.';
      throw new Error("Variant not found in VariantValidator.");
    }
  } catch (e) {
    if (e.name === 'AbortError') return;
    console.error('VariantValidator Persistent Failure:', e);
    const statusInfo = e.status ? ` (HTTP ${e.status})` : '';
    setDot('vv', 'error');
    document.getElementById('statusMsg').innerText = `VariantValidator error${statusInfo}: ${e.message}`;
    throw e; // Bubble up
  }
}

async function fetchEnsembl(hgvs, gen) {
  setDot('vep', 'loading');
  try {
    // Ensembl REQUIRES transcript version for RefSeq lookups — do NOT strip it
    const res = await fetchWithRetry(`https://rest.ensembl.org/vep/human/hgvs/${encodeURIComponent(hgvs)}?content-type=application/json&mane=1&numbers=1&hgvs=1&refseq=1&af_gnomad=1&CADD=1&AlphaMissense=1`);
    if (!res.ok) throw new Error('VEP fail');
    const json = await res.json(); if (!json || !json.length) throw new Error('No VEP');
    if (gen !== window.currentSearchGen()) return;
    const txList = json[0].transcript_consequences || [];
    txList.sort((a, b) => getSeverity(a.consequence_terms) - getSeverity(b.consequence_terms));
    const best = selectByTranscript(txList, {
      gene:        data.coords.gene,
      transcript:  data.ensembl.transcript,
      getGene:     t => t.gene_symbol,
      getTxId:     t => t.transcript_id,
      isCanonical: t => !!(t.mane && t.mane.includes('MANE_Select') &&
                          (t.transcript_id.startsWith('NM_') || t.transcript_id.startsWith('XM_'))) ||
                        t.canonical === 1,
    });
    if (best) {
      if (best.exon) {
        document.getElementById('eExon').innerText = 'Exon ' + best.exon;
        const parts = best.exon.split('/');
        if (parts.length === 2) {
          data.ensembl.vepExon = parseInt(parts[0]);
          data.ensembl.vepTotalExons = parseInt(parts[1]);
        }
      } else if (best.intron) {
        document.getElementById('eExon').innerText = 'Intron ' + best.intron;
      } else {
        document.getElementById('eExon').innerText = '-';
      }
      const refId = best.transcript_id.startsWith('NM_') || best.transcript_id.startsWith('XM_') ? best.transcript_id : best.mane_select || best.transcript_id;
      document.getElementById('eManeTranscript').innerText = refId;
      data.ensembl.transcript = refId;
      data.proteinId = best.protein_id || null;
      console.log('[fetchEnsembl] MANE Protein ID:', data.proteinId);
      const cons = (best.consequence_terms || []).join(', ');
      document.getElementById('eConsequence').innerText = cons || '-';
      data.ensembl.vepConsequence = cons;
      if (best.amino_acids) data.ensembl.vepAminoAcids = best.amino_acids;
      if (best.protein_start) data.ensembl.vepProteinStart = best.protein_start;
      if (best.hgvsp) data.ensembl.protein = best.hgvsp;
      // Populate gene from VEP if VV hasn't set it yet (e.g. VV is still retrying)
      if (!data.coords.gene && best.gene_symbol) data.coords.gene = best.gene_symbol.toUpperCase();

      // Raw VEP values for two-source validation engine (do not feed back into other modules)
      if (best.gene_symbol)  data.ensembl._geneFromVep       = best.gene_symbol.toUpperCase();
      if (refId)             data.ensembl._transcriptFromVep = refId;
      if (best.hgvsc)        data.ensembl._hgvscFromVep      = best.hgvsc.split(':').pop();
      if (best.hgvsp)        data.ensembl._hgvspFromVep      = best.hgvsp;

      console.log('[fetchEnsembl] Best transcript consequence:', best);

      // ── AlphaMissense Extraction (Zero-Latency Ensembl Integration) ──
      const getAMData = (t) => t.alphamissense || t.am_pathogenicity; // Handles both nested and flat variants
      const getAMScore = (data) => (typeof data === 'object') ? data.am_pathogenicity : data;
      const getAMClass = (data, t) => (typeof data === 'object') ? data.am_class : (t.am_class || t.alphamissense_class || '');

      let rawData = getAMData(best);
      let amScore = getAMScore(rawData);
      let amCls = getAMClass(rawData, best);

      // Fallback: Search all transcripts if the best one is missing AM data
      if (amScore === undefined || amScore === null) {
        const altWithAM = txList.find(t => getAMData(t) !== undefined);
        if (altWithAM) {
          rawData = getAMData(altWithAM);
          amScore = getAMScore(rawData);
          amCls = getAMClass(rawData, altWithAM);
        }
      }

      if (amScore !== undefined && amScore !== null) {
        data.scores.alphaMissenseScore = parseFloat(amScore);
        // Capitalize the first letter (e.g., "likely_benign" -> "Likely_benign")
        let formattedCls = String(amCls).replace(/_/g, ' ');
        data.scores.alphaMissensePred = formattedCls ? formattedCls.charAt(0).toUpperCase() + formattedCls.slice(1) : null;
      }
      renderMetricPanels();
    }

    // ── Transcript Picker Dropdown ──────────────────────────────────────────
    // Build a dropdown of all transcripts with their HGVS c./p. and exon/intron.
    // Selecting an alternate transcript updates the display fields without altering
    // data.ensembl.transcript (MANE) or triggering re-analysis.
    const allTxData = txList.map(t => ({
        id: t.transcript_id,
        hgvsc: t.hgvsc ? t.hgvsc.split(':').pop() : null,
        hgvsp: t.hgvsp ? t.hgvsp.split(':').pop() : null,
        exon: t.exon || null,
        intron: t.intron || null,
        consequence: (t.consequence_terms || []).join(', '),
        isMane: !!(t.mane && t.mane.includes('MANE_Select'))
      }));
    data.ensembl.allTranscripts = allTxData;

    if (allTxData.length > 1) {
      const maneId = data.ensembl.transcript;
      const select = document.getElementById('transcriptDropdown');
      select.innerHTML = allTxData.map(t => {
        const label = `${t.id}${t.isMane ? ' ★' : ''} — ${t.hgvsc || t.consequence}`;
        return `<option value="${t.id}"${t.id === maneId ? ' selected' : ''}>${label}</option>`;
      }).join('');

      // onchange: only writes to transcriptDropdownDetail (ensembl-owned).
      // Never touches dHgvs / dProtein / eExon / eConsequence — those are owned
      // by the coords/VV module and the validation engine; cross-writing them
      // causes state mismatches when renderValidationMarkers() re-fires.
      select.onchange = () => {
        const tx = allTxData.find(t => t.id === select.value);
        const detail = document.getElementById('transcriptDropdownDetail');
        if (!detail) return;
        if (!tx || tx.isMane) { detail.innerText = ''; return; }
        const loc = tx.exon ? `Exon ${tx.exon}` : tx.intron ? `Intron ${tx.intron}` : '';
        detail.innerText = [loc, tx.hgvsp || '', tx.consequence].filter(Boolean).join(' | ');
      };

      document.getElementById('rowTranscriptPicker').style.display = '';
    }

    // Trigger multi-transcript literature search with current Ensembl data.
    // ClinVar VCV may add historical notations later via a second _scheduleLitNotations call.
    _scheduleLitNotations(gen);

    // ── Mutalyzer Genomic HGVS (splice / intronic variants) ────────────────
    // For variants with intronic offset in c. notation, fetch the NC_-anchored
    // genomic HGVS from Mutalyzer (e.g. NC_000010.11(NM_015634.4):c.606-1G>T).
    const SPLICE_TERMS = new Set([
      'splice_donor_variant', 'splice_acceptor_variant', 'splice_region_variant',
      'intron_variant', 'splice_donor_region_variant', 'splice_polypyrimidine_tract_variant',
      'splice_donor_5th_base_variant'
    ]);
    const bestCons = best ? (best.consequence_terms || []) : [];
    const isIntronic = !!best?.intron || bestCons.some(t => SPLICE_TERMS.has(t));
    // Use best.hgvsc (full NM_...:c.xxx from VEP) — data.coords.hgvs is only the c. portion
    if (isIntronic && best?.hgvsc) {
      fetchMutalyzer(best.hgvsc, gen);
    }

    data.coords.strand = json[0].strand || 1;
    data.ensembl.vepConsequence = data.ensembl.vepConsequence || json[0].most_severe_consequence || null;
    document.getElementById('eConsequence').innerText = data.ensembl.vepConsequence || '-';

    let rsId = '-';
    if (json[0].id && json[0].id.startsWith('rs')) rsId = json[0].id;
    else if (json[0].colocated_variants) {
      const rv = json[0].colocated_variants.find(v => v.id && v.id.startsWith('rs'));
      if (rv) rsId = rv.id;
    }
    document.getElementById('eRsId').innerText = rsId;
    data.ensembl.rsId = rsId !== '-' ? rsId : null;

    // CROSS-CHECK VALIDATION (VV vs Ensembl MANE)
    if (data.coords.primaryFromVV && best && best.mane_select) {
      const eMane = document.getElementById('eManeTranscript');
      // Validate Gene
      if (best.gene_symbol?.toUpperCase() === data.coords.gene?.toUpperCase()) markValidated('dGene');
      // Validate Transcript
      if (best.transcript_id === eMane?.innerText || best.transcript_id === data.ensembl.transcript) markValidated('eManeTranscript');
      // Validate HGVS c.
      if (best.hgvsc === data.coords.hgvs) markValidated('dHgvs');
      // Validate HGVS p.
      const normalizedPv = (document.getElementById('dProtein')?.innerText || data.ensembl.protein || '').replace(/\(.*\)/g, '');
      if (best.hgvsp && best.hgvsp.includes(normalizedPv)) markValidated('dProtein');
    }

    // fetchMyVariant removed; now triggered by fetchBroadSpliceAI via triggerDownstreamAPIs
    fetchMyVariant(gen); // Always called here, after data.ensembl.rsId is guaranteed to be set
    evaluateStartLoss(); // Fire the Start-Loss check only after VEP has resolved the consequence
    evaluateSplicePS1(); // Splice PS1 audit (owns splice variants — hybrid split from fetchClinVarCodon)
    fetchClinVarCodon(gen);
    if (json[0].seq_region_name) {
      const chr = json[0].seq_region_name, pos = json[0].start;
      let ref = '', alt = '';
      if (json[0].allele_string) {
        const p = json[0].allele_string.split('/');
        ref = p[0];
        alt = p[1] || '';
      }

      let ensRef = normalizeVcfAllele(ref.replace('-', ''));
      let ensAlt = normalizeVcfAllele(alt.replace('-', ''));
      const cleanChr = chr.replace('chr', '');

      // [FIX] Smart Reverse Complement (VV Mode):
      // If we are on the minus strand, check if Ensembl's alleles need RC to match VV.
      if (json[0].strand === -1) {
        const comp = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '' };
        const reverseComp = seq => seq.split('').map(n => comp[n] || n).reverse().join('');

        const ensemblHg38_orig = `chr${cleanChr}-${pos}-${ensRef}-${ensAlt}`;
        const pV = s => {
          const p = s.split('-');
          return { c: p[0].replace('chr', ''), p: parseInt(p[1]), r: p[2] || '', a: p[3] || '' };
        };
        const isCoordMatch = (e, v) => {
          try {
            const ev = pV(e), vv = pV(v);
            return ev.c === vv.c && Math.abs(ev.p - vv.p) <= 10;
          } catch (err) { return false; }
        };

        if (data.coords.hg38String && isCoordMatch(ensemblHg38_orig, data.coords.hg38String)) {
          const vv = pV(data.coords.hg38String);
          const rcRef = reverseComp(ensRef), rcAlt = reverseComp(ensAlt);
          if ((rcRef.length - rcAlt.length) === (vv.r.length - vv.a.length)) {
            ensRef = rcRef; ensAlt = rcAlt;
          }
        }
      }

      // [FIX] Robust Cross-validation matching (tolerant of indel representations)
      const ensemblHg38 = `chr${cleanChr}-${pos}-${ensRef}-${ensAlt}`;
      data.ensembl._hg38FromVep = ensemblHg38;
      const isMatch = (e, v) => {
        if (!e || !v) return false;
        if (e === v) return true;
        const pV = s => {
          const p = s.split('-');
          return { c: p[0].replace('chr', ''), p: parseInt(p[1]), r: (p[2] || '').replace('-', ''), a: (p[3] || '').replace('-', '') };
        };
        try {
          const ev = pV(e), vv = pV(v);
          if (ev.c !== vv.c) return false;
          return Math.abs(ev.p - vv.p) <= 10 && (ev.r.length - ev.a.length) === (vv.r.length - vv.a.length);
        } catch (err) { return false; }
      };

      if (data.coords.hg38String && isMatch(ensemblHg38, data.coords.hg38String)) {
        markValidated('dHg38');
      }

      // Pass local ensRef/ensAlt to liftover helper to keep global data object untouched
      await fetchEnsemblMap(cleanChr, pos, ensRef, ensAlt);
    }
    setDot('vep', 'ok');
    runValidationPass();
  } catch (e) {
    if (e.name === 'AbortError') return;
    console.error(e); setDot('vep', 'warn'); document.getElementById('eConsequence').innerText = '(VEP Error)';
    // VEP failed but VV may have set protein/consequence — still attempt codon/start-loss/splice lookup
    if (data.ensembl.protein || data.coords.hgvs !== '-') {
      fetchClinVarCodon(gen);
      evaluateStartLoss();
      evaluateSplicePS1();
    }
  }
}

/**
 * Appends a verification checkmark to a UI element to confirm 
 * that primary data (VV) has been cross-validated by a secondary source (Ensembl).
 */
// markValidated is now a no-op — all field markers (✅/⚠/ℹ️) are driven by the
// two-source validation engine (see validation.js + renderValidationMarkers in ui.js).
// Kept as a stub so legacy call sites compile without errors.
function markValidated(_id) { /* no-op */ }

/**
 * fetchMutalyzer(hgvs, gen)
 * Inputs:  NM_-based HGVS string (data.coords.hgvs), gen guard
 * Outputs: data.ensembl.mutalyzerHgvs, #eMutalyzerHgvs, #rowMutalyzer
 * Failures: Silent — non-blocking, row stays hidden on error.
 *
 * Called for splice/intronic variants to retrieve the NC_-anchored genomic
 * HGVS (e.g. NC_000010.11(NM_015634.4):c.606-1G>T) from Mutalyzer normalize.
 */
async function fetchMutalyzer(hgvs, gen) {
  try {
    const res = await fetch(`https://mutalyzer.nl/api/normalize/${encodeURIComponent(hgvs)}`);
    // Mutalyzer returns 422 for intronic NM_ input but embeds the GRCh38 answer in the body
    if (res.status !== 200 && res.status !== 422) return;
    if (gen !== window.currentSearchGen()) return;
    const json = await res.json();

    let genomicHgvs = null;

    if (res.status === 422) {
      // Intronic variant: GRCh38 option is in custom.errors[].options
      const errors = json.custom?.errors || [];
      for (const err of errors) {
        const opt = (err.options || []).find(o => o.assembly_id === 'GRCh38');
        if (opt?.description) { genomicHgvs = opt.description; break; }
      }
    } else {
      // Successful normalization: prefer NC_-anchored coding description
      const cDescs = json.equivalent_descriptions?.c || [];
      genomicHgvs = (cDescs.find(d => d.description?.startsWith('NC_')) || {}).description;
      // Fallback: pure genomic g. notation
      if (!genomicHgvs) {
        const gDescs = json.equivalent_descriptions?.g || [];
        genomicHgvs = (gDescs[0] || {}).description;
      }
    }

    if (genomicHgvs) {
      data.ensembl.mutalyzerHgvs = genomicHgvs;
      document.getElementById('eMutalyzerHgvs').innerText = genomicHgvs;
      document.getElementById('rowMutalyzer').style.display = '';
    }
  } catch (e) {
    console.warn('[fetchMutalyzer] Failed:', e.message);
  }
}

// ── MULTI-TRANSCRIPT LITERATURE SEARCH ───────────────────────────────────────
// Namespace: data.litvar.* (new module, post 2026-05-21 — strict isolation rules apply)
// Rate limit: all NCBI calls use fetchWithRetry → ncbiAcquireToken (9 req/s token bucket)

/**
 * parseVCVProteinNotations(xmlText)
 * Pure function — no side effects.
 * Extracts all ProteinExpression notations from ClinVar VCV XML, including
 * historical transcript versions. Returns array of unique p. strings.
 */
function parseVCVProteinNotations(xmlText) {
  try {
    const doc = new DOMParser().parseFromString(xmlText, 'application/xml');
    const notations = new Set();
    doc.querySelectorAll('ProteinExpression Expression').forEach(el => {
      const text = el.textContent.trim();
      const p = text.includes(':') ? text.split(':')[1] : text;
      if (p && p.startsWith('p.') && p !== 'p.?' && p !== 'p.=') notations.add(p);
    });
    return [...notations];
  } catch (e) {
    console.warn('[parseVCVProteinNotations] Parse error:', e.message);
    return [];
  }
}

/**
 * expandProteinNotation(p)
 * Pure function. For a given p. notation, returns all equivalent search terms:
 *   Missense:   p.Arg186Cys  → + p.R186C
 *   Stop codon: p.Arg186*    → + p.Arg186Ter, p.R186*, p.R186X
 *   Frameshift: p.Arg186fs*5 → + p.Arg186fs, p.R186fs
 * Uses AA3TO1 from ui.js (loads before api.js in script order).
 */
function expandProteinNotation(p) {
  const terms = new Set([p]);
  // Missense: 3-letter → 1-letter
  const mMis = p.match(/^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$/);
  if (mMis) {
    const r1 = AA3TO1[mMis[1]], a1 = AA3TO1[mMis[3]];
    if (r1 && a1) terms.add(`p.${r1}${mMis[2]}${a1}`);
  }
  // Stop codon: *, Ter, X — all three forms in both 1-letter and 3-letter
  const mStop = p.match(/^p\.([A-Z][a-z]{2}|[A-Z])(\d+)(\*|Ter|X)$/);
  if (mStop) {
    const is3 = mStop[1].length === 3;
    const ref3 = is3 ? mStop[1] : Object.keys(AA3TO1).find(k => AA3TO1[k] === mStop[1]);
    const ref1 = is3 ? (AA3TO1[mStop[1]] || '') : mStop[1];
    if (ref3) { terms.add(`p.${ref3}${mStop[2]}*`); terms.add(`p.${ref3}${mStop[2]}Ter`); }
    if (ref1) { terms.add(`p.${ref1}${mStop[2]}*`); terms.add(`p.${ref1}${mStop[2]}X`); }
  }
  // Frameshift: shorten to gene+pos+fs (papers rarely cite the full fs*N form)
  const mFs = p.match(/^p\.([A-Z][a-z]{2}|[A-Z])(\d+).+fs/);
  if (mFs) {
    const is3 = mFs[1].length === 3;
    const ref3 = is3 ? mFs[1] : Object.keys(AA3TO1).find(k => AA3TO1[k] === mFs[1]);
    const ref1 = is3 ? (AA3TO1[mFs[1]] || '') : mFs[1];
    if (ref3) terms.add(`p.${ref3}${mFs[2]}fs`);
    if (ref1) terms.add(`p.${ref1}${mFs[2]}fs`);
  }
  return [...terms];
}

/**
 * buildPubMedOrTerm(gene, notations)
 * Pure function. Expands all notations and builds a PubMed OR query string.
 * Returns null if no notations available.
 *
 * PubMed tokenizes the dot in "p.Arg273Cys" as a sentence boundary, so the
 * "p." prefix must be stripped for [tiab] queries. The asterisk (*) in stop-
 * codon notation is a PubMed wildcard and must be replaced with "Ter".
 */
function buildPubMedOrTerm(gene, notations, diseases = []) {
  if (!gene || !notations.length) return null;
  const allTerms = [...new Set(notations.flatMap(p => expandProteinNotation(p)))];
  const sanitize = t => t
    .replace(/^p\./, '')          // strip p. prefix (dot breaks PubMed phrase matching)
    .replace(/\*/g, 'Ter');       // * is a PubMed wildcard — use Ter instead
  const orClause = allTerms.map(t => `${sanitize(t)}[tiab]`).join(' OR ');
  let term = `${gene}[tiab] AND (${orClause})`;
  if (diseases.length > 0) {
    const diseaseClause = diseases.map(d => `"${d}"[tiab]`).join(' OR ');
    term += ` AND (${diseaseClause})`;
  }
  return term;
}

/**
 * buildPubMedUrl(gene, notations)
 * Builds a clickable PubMed URL. Caps included notations to avoid exceeding
 * the ~2000-char browser URL limit. Always uses the top N notations (sorted
 * alphabetically for stability, MANE-first via caller's ordering).
 */
function buildPubMedUrl(gene, notations, diseases = []) {
  let subset = notations;
  let url;
  do {
    const term = buildPubMedOrTerm(gene, subset, diseases);
    if (!term) return null;
    url = `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(term)}`;
    if (url.length <= 1800 || subset.length <= 1) break;
    subset = subset.slice(0, Math.floor(subset.length * 0.7));
  } while (true);
  return url;
}

/**
 * _scheduleLitNotations(gen)
 * Sync gate called from both fetchEnsembl and fetchClinVar.
 * Fires fetchLiteratureNotations only when the set of available p. notations
 * has actually changed (key-based dedup prevents duplicate NCBI calls).
 * Inputs passed explicitly to fetchLiteratureNotations — no implicit reads there.
 */
function _scheduleLitNotations(gen) {
  if (gen !== window.currentSearchGen()) return;
  const allP = [
    ...(data.ensembl.allTranscripts || []).map(t => t.hgvsp).filter(Boolean),
    ...(data.clinvar.historicalProteinNotations || [])
  ];
  const unique = [...new Set(allP.filter(p => p && p !== 'p.?' && p !== 'p.='))].sort();
  if (!unique.length) return;
  // Include disease/phenotype terms in the key so re-search fires when they change
  const diseases = _buildDiseaseTerms();
  const key = unique.join('|') + '\x00' + diseases.join('|');
  if (key === data.litvar._lastKey) return;
  data.litvar._lastKey = key;
  fetchLiteratureNotations(
    data.coords.gene,
    data.ensembl.allTranscripts || [],
    data.clinvar.historicalProteinNotations || [],
    diseases,
    gen
  );
}

// Builds disease/phenotype terms for the PubMed query.
// User-typed HPO phenotype labels are prioritized (natural language, matches papers).
// Falls back to ClinGen disease names only when no phenotypes are entered.
function _buildDiseaseTerms() {
  const pheno = (data.ptPhenoTexts || []).slice(0, 3);
  if (pheno.length) return pheno;
  return (data.associatedConditions || []).map(c => c.diseaseName).filter(Boolean).slice(0, 2);
}

/**
 * fetchLiteratureNotations(gene, currentTxs, historical, gen)
 * Inputs:  gene, currentTxs (from ensembl namespace), historical (from clinvar namespace) — all explicit
 * Outputs: data.litvar.combinedTotal, data.litvar.uniqueNotations → #rowPubMedMultiTx, #pMultiTxTotal
 * Failures: Silent — non-blocking. Row stays hidden on error.
 * Rate limit: fetchWithRetry → ncbiAcquireToken handles 9 req/s automatically.
 */
async function fetchLiteratureNotations(gene, currentTxs, historical, diseases, gen) {
  if (!gene) return;
  const allP = [...new Set([
    ...currentTxs.map(t => t.hgvsp).filter(p => p && p !== 'p.?' && p !== 'p.='),
    ...historical.filter(p => p && p !== 'p.?' && p !== 'p.=')
  ])];
  if (!allP.length) return;
  data.litvar.uniqueNotations = allP;

  const orTerm = buildPubMedOrTerm(gene, allP, diseases);
  if (!orTerm) return;

  try {
    const countUrl = `${NCBI_BASE}/esearch.fcgi?db=pubmed&retmode=json&retmax=0&term=${encodeURIComponent(orTerm)}&api_key=${NCBI_KEY}&email=${NCBI_EMAIL}`;
    const r = await fetchWithRetry(countUrl); // auto-throttled via ncbiAcquireToken
    if (gen !== window.currentSearchGen()) return;
    if (!r.ok) return;
    const j = await r.json();
    const total = parseInt(j.esearchresult?.count || 0);
    data.litvar.combinedTotal = total;
    const pmUrl = buildPubMedUrl(gene, allP, diseases);
    renderLitMultiTxUI(total, allP, pmUrl, gene);
  } catch (e) {
    if (e.name === 'AbortError') return;
    console.warn('[fetchLiteratureNotations] Failed:', e.message);
  }
}

/**
 * loadLitBreakdown(gene)
 * Lazy per-notation PubMed counts, fired on user expand click.
 * Sequential loop — each call goes through ncbiAcquireToken, safe under 9 req/s.
 * Inputs: gene (from DOM context, safe for event handler), data.litvar.uniqueNotations (read)
 * Outputs: #pMultiTxBreakdown only — no data.* writes.
 */
async function loadLitBreakdown(gene) {
  const breakdownEl = document.getElementById('pMultiTxBreakdown');
  if (!breakdownEl) return;
  const notations = data.litvar.uniqueNotations || [];
  if (!notations.length) return;

  const btn = document.querySelector('#pMultiTxTotal button');
  if (btn) btn.remove();

  breakdownEl.style.display = '';
  breakdownEl.innerHTML = '<div style="font-size:0.65rem; color:var(--dim); padding:2px 0;">Loading per-notation counts…</div>';

  const gen = window.currentSearchGen();
  const rows = [];

  for (const p of notations) {
    if (gen !== window.currentSearchGen()) return;
    const expanded = expandProteinNotation(p);
    const sanitize = t => t.replace(/^p\./, '').replace(/\*/g, 'Ter');
    const orClause = expanded.map(t => `${sanitize(t)}[tiab]`).join(' OR ');
    const term = `${gene}[tiab] AND (${orClause})`;
    const pmUrl = `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(term)}`;
    try {
      const r = await fetchWithRetry(`${NCBI_BASE}/esearch.fcgi?db=pubmed&retmode=json&retmax=0&term=${encodeURIComponent(term)}&api_key=${NCBI_KEY}&email=${NCBI_EMAIL}`);
      if (!r.ok) continue;
      const j = await r.json();
      rows.push({ p, count: parseInt(j.esearchresult?.count || 0), pmUrl });
    } catch (e) {
      if (e.name === 'AbortError') return;
    }
  }

  rows.sort((a, b) => b.count - a.count);
  const withHits = rows.filter(r => r.count > 0);
  const zeroCount = rows.length - withHits.length;

  if (!withHits.length) {
    breakdownEl.innerHTML = '<div style="font-size:0.65rem; color:var(--dim); padding:2px 0;">No individual notation results.</div>';
    return;
  }

  breakdownEl.innerHTML = withHits.slice(0, 5).map(r =>
    `<div style="display:flex; justify-content:space-between; align-items:center; font-size:0.65rem; padding:1px 0;">
       <span style="font-family:'JetBrains Mono',monospace; color:var(--dim);">${r.p}</span>
       <a href="${r.pmUrl}" target="_blank" style="color:var(--teal); text-decoration:none; font-weight:600;">${r.count} ↗</a>
     </div>`
  ).join('');

  if (zeroCount > 0)
    breakdownEl.innerHTML += `<div style="font-size:0.6rem; color:var(--dim); opacity:0.6; margin-top:2px;">${zeroCount} notation${zeroCount > 1 ? 's' : ''} with 0 hits hidden</div>`;
}

/**
 * Captures raw ClinGen Allele Registry values into the isolated data.clingenAR.*
 * namespace. Used by the two-source validation engine (see validation.js).
 * Does NOT write to data.coords.* or data.ensembl.* — read-only data extraction.
 */
function captureClinGenARFields(ar, gene, transcript) {
  if (!ar || !data.clingenAR) return;

  // Genomic coordinates
  const g38 = ar.genomicAlleles?.find(g => g.referenceGenome === 'GRCh38');
  const c38 = g38?.coordinates?.[0];
  if (c38 && g38?.chromosome) {
    const ref = normalizeVcfAllele(c38.referenceAllele);
    const alt = normalizeVcfAllele(c38.allele);
    if (ref && alt) {
      data.clingenAR.hg38String = `chr${String(g38.chromosome).replace('chr', '')}-${c38.start + 1}-${ref}-${alt}`;
    }
  }
  const g37 = ar.genomicAlleles?.find(g => g.referenceGenome === 'GRCh37');
  const c37 = g37?.coordinates?.[0];
  if (c37 && g37?.chromosome) {
    const ref = normalizeVcfAllele(c37.referenceAllele);
    const alt = normalizeVcfAllele(c37.allele);
    if (ref && alt) {
      data.clingenAR.hg19String = `chr${String(g37.chromosome).replace('chr', '')}-${c37.start + 1}-${ref}-${alt}`;
    }
  }

  // Transcript / gene / HGVS — prefer MANE Select
  const txAlleles = ar.transcriptAlleles || [];
  const isManeSelect = t => {
    const m = t.MANE ?? t.mane ?? '';
    return typeof m === 'string' ? m.toUpperCase().includes('SELECT') : !!m;
  };
  const maneTx = selectByTranscript(txAlleles, {
    gene,
    transcript,
    getGene:     t => t.geneSymbol || t.gene?.label || t.gene?.symbol,
    getTxId:     t => t.hgvsMatchStrings?.[0]?.split(':')[0],
    isCanonical: isManeSelect,
  });

  if (maneTx) {
    data.clingenAR.gene = maneTx.geneSymbol || maneTx.gene?.label || maneTx.gene?.symbol || null;
    const hgvsStr = maneTx.hgvsMatchStrings?.[0] || '';
    if (hgvsStr.includes(':')) {
      const colonIdx = hgvsStr.indexOf(':');
      data.clingenAR.transcript = hgvsStr.slice(0, colonIdx);
      data.clingenAR.hgvsC = hgvsStr.slice(colonIdx + 1);
    }
    const pEffect = maneTx.proteinEffect?.hgvs || maneTx.proteinEffect?.hgvsString;
    if (pEffect) data.clingenAR.hgvsP = pEffect.includes(':') ? pEffect.split(':').pop() : pEffect;
  }

  // CA ID
  data.clingenAR.caId = ar.caid || ar['@id']?.split('/').pop() || null;
}

/**
 * Cross-validates Ensembl-resolved hg38/hg19 coords against the ClinGen Allele Registry.
 * AR always returns forward-strand VCF alleles (strand-reliable ground truth).
 * On match → adds ✅. On mismatch → corrects data and shows ⚠.
 * Also populates CAid if Ensembl didn't resolve it.
 */
async function confirmCoordsWithClinGenAR(hgvsInput, gen) {
  if (!hgvsInput || hgvsInput === '-') return;
  try {
    const clean = hgvsInput.replace(/\([^)]*\)/g, '').trim();
    const res = await fetchWithTimeout(`https://reg.genome.network/allele?hgvs=${encodeURIComponent(clean)}`);
    if (!res.ok) return;
    const ar = await res.json();

    // Capture raw AR values for two-source validation (separate from coords correction logic)
    captureClinGenARFields(ar, data.coords.gene, data.ensembl.transcript);

    const pV = s => {
      const p = s.split('-');
      return { c: p[0]?.replace('chr', ''), pos: parseInt(p[1]), r: (p[2] || '').replace('-', ''), a: (p[3] || '').replace('-', '') };
    };
    const coordsAgree = (arStr, localStr) => {
      if (!arStr || !localStr || localStr === '—') return false;
      if (arStr === localStr) return true;
      try {
        const av = pV(arStr), lv = pV(localStr);
        return av.c === lv.c && Math.abs(av.pos - lv.pos) <= 2 && (av.r.length - av.a.length) === (lv.r.length - lv.a.length);
      } catch { return false; }
    };

    // AR returns insertions/deletions HGVS-style (ref or alt = '-'), which is
    // not equivalent to local VCF form (anchored with previous base). When AR
    // is HGVS-style, skip the VCF string comparison entirely — verify semantically
    // (same chrom + same indel size) and trust the already-correct local coords.
    const isHgvsStyle = (c) => c?.referenceAllele === '-' || c?.allele === '-'
                            || c?.allele == null || c?.allele === '';
    const sizeDelta = (c) => (c?.referenceAllele || '').replace('-', '').length
                           - (c?.allele || '').replace('-', '').length;
    const localSizeDelta = (s) => {
      if (!s || s === '—') return null;
      const p = s.split('-');
      return (p[2] || '').length - (p[3] || '').length;
    };

    // GRCh38
    const g38 = ar.genomicAlleles?.find(g => g.referenceGenome === 'GRCh38');
    const c38 = g38?.coordinates?.[0];
    if (c38 && g38?.chromosome) {
      const arChrom = String(g38.chromosome).replace('chr', '');
      const arPos   = String(c38.start + 1);
      const arRef38 = normalizeVcfAllele(c38.referenceAllele);
      const arAlt38 = normalizeVcfAllele(c38.allele);
      const arHg38  = `chr${arChrom}-${arPos}-${arRef38}-${arAlt38}`;

      if (isHgvsStyle(c38)) {
        // AR is HGVS-style; verify semantically against local VCF
        const localDelta = localSizeDelta(data.coords.hg38String);
        const localChrom = data.coords.hg38String?.split('-')[0]?.replace('chr', '');
        if (localChrom === arChrom && localDelta === sizeDelta(c38)) {
          markValidated('dHg38');
        }
        // else: silent — local VCF is authoritative
      } else if (coordsAgree(arHg38, data.coords.hg38String)) {
        markValidated('dHg38');
      } else if (!data.coords.hg38String || data.coords.hg38String.includes('--')) {
        // hg38String was empty or invalid — populate from AR (only safe when AR is VCF-style)
        data.coords.chrom = arChrom; data.coords.pos38 = arPos;
        data.coords.ref = arRef38; data.coords.alt = arAlt38;
        data.coords.hg38String = arHg38;
        const el38 = document.getElementById('dHg38');
        if (el38) el38.innerText = arHg38;
      } else {
        console.warn(`[ClinGen AR] hg38 correction: ${data.coords.hg38String} → ${arHg38}`);
        data.coords.chrom = arChrom; data.coords.pos38 = arPos;
        data.coords.ref = arRef38; data.coords.alt = arAlt38;
        data.coords.hg38String = arHg38;
        const el38 = document.getElementById('dHg38');
        if (el38) el38.innerHTML = `${arHg38} <span style="color:var(--amber);font-size:0.75rem;vertical-align:middle;" title="Corrected by ClinGen Allele Registry">⚠</span>`;
      }
    }

    // GRCh37
    const g37 = ar.genomicAlleles?.find(g => g.referenceGenome === 'GRCh37');
    const c37 = g37?.coordinates?.[0];
    if (c37 && g37?.chromosome) {
      const arChrom37 = String(g37.chromosome).replace('chr', '');
      const arRef37 = normalizeVcfAllele(c37.referenceAllele);
      const arAlt37 = normalizeVcfAllele(c37.allele);
      const arHg19 = `chr${arChrom37}-${String(c37.start + 1)}-${arRef37}-${arAlt37}`;
      const current19 = document.getElementById('dHg19')?.innerText;

      if (isHgvsStyle(c37)) {
        const localDelta = localSizeDelta(data.coords.hg19String);
        const localChrom = data.coords.hg19String?.split('-')[0]?.replace('chr', '');
        if (localChrom === arChrom37 && localDelta === sizeDelta(c37)) {
          markValidated('dHg19');
        }
      } else if (coordsAgree(arHg19, data.coords.hg19String)) {
        markValidated('dHg19');
      } else if (!data.coords.hg19String || data.coords.hg19String === '—' || !current19 || current19 === '—') {
        // Ensembl liftover produced nothing — use AR as primary source (VCF-style only)
        data.coords.pos19 = String(c37.start + 1);
        data.coords.ref = arRef37; data.coords.alt = arAlt37;
        data.coords.hg19String = arHg19;
        document.getElementById('dHg19').innerText = arHg19;
        markValidated('dHg19');
      } else {
        console.warn(`[ClinGen AR] hg19 correction: ${data.coords.hg19String} → ${arHg19}`);
        data.coords.pos19 = String(c37.start + 1);
        data.coords.ref = arRef37; data.coords.alt = arAlt37;
        data.coords.hg19String = arHg19;
        const el19 = document.getElementById('dHg19');
        if (el19) el19.innerHTML = `${arHg19} <span style="color:var(--amber);font-size:0.75rem;vertical-align:middle;" title="Corrected by ClinGen Allele Registry">⚠</span>`;
      }
    }

    // Populate CAid if Ensembl didn't resolve it
    const caid = ar.caid || ar['@id']?.split('/').pop() || null;
    if (caid && !data.caId) {
      data.caId = caid;
      const caEl = document.getElementById('eCAId');
      if (caEl) caEl.innerText = caid;
    }
    runValidationPass();
  } catch (e) {
    if (e.name === 'AbortError') return;
    console.warn('[ClinGen AR confirm] Non-blocking error:', e.message);
  }
}

async function fetchEnsemblMap(chr, pos, ref, alt, gen) {
  try {
    const end = parseInt(pos) + Math.max(ref.length, 1) - 1;
    const res = await fetchWithRetry(`https://rest.ensembl.org/map/human/GRCh38/${chr}:${pos}..${end}:1/GRCh37?content-type=application/json`);
    if (!res.ok) return;
    const json = await res.json();
    if (json.mappings?.length) {
      const m = json.mappings[0].mapped;
      // Safety guard: don't write malformed '..--ACGT' strings if ref/alt are empty
      // (can happen if called with HGVS-style alleles before VCF normalisation runs)
      if (!ref || !alt) return;
      const hg19Str = `chr${chr.replace('chr', '')}-${m.start}-${ref}-${alt}`;
      const currentVal = document.getElementById('dHg19').innerText;
      // CLINICAL GUARD: Never overwrite VV data if it's already clinical-grade.
      if (!data.coords.hg19FromVV && (!currentVal || currentVal === '—' || currentVal === '')) {
        document.getElementById('dHg19').innerText = hg19Str;
        data.coords.hg19String = hg19Str;
        data.coords.pos19 = m.start;
      }
      // [FIX] Cross-validate hg19 coordinates with VV output
      if (data.coords.hg19FromVV && data.coords.hg19String && hg19Str === data.coords.hg19String) {
        markValidated('dHg19');
      } else if (data.coords.hg19FromVV && data.coords.hg19String) {
        // Indel-tolerant matching for hg19
        const pV = s => { const p = s.split('-'); return { c: p[0].replace('chr', ''), p: parseInt(p[1]), r: (p[2] || '').replace('-', ''), a: (p[3] || '').replace('-', '') }; };
        try {
          const ev = pV(hg19Str), vv = pV(data.coords.hg19String);
          if (ev.c === vv.c && Math.abs(ev.p - vv.p) <= 10 && (ev.r.length - ev.a.length) === (vv.r.length - vv.a.length)) {
            markValidated('dHg19');
          }
        } catch (e) { }
      }
    }
  } catch (e) { }
}

async function fetchGnomAD(gen) {
  if (!data.coords.chrom || !data.coords.pos38 || !data.coords.ref || !data.coords.alt) return;
  setAPIStatus('gnomad', 'loading');
  const query = `query { variant(variantId: "${data.coords.chrom}-${data.coords.pos38}-${data.coords.ref}-${data.coords.alt}", dataset: gnomad_r4) { 
    genome { 
      filters ac an homozygote_count hemizygote_count 
      faf95 { popmax } 
      populations { id ac an homozygote_count hemizygote_count } 
    } 
    exome { 
      filters ac an homozygote_count hemizygote_count 
      faf95 { popmax } 
      populations { id ac an homozygote_count hemizygote_count } 
    }
    joint { 
      filters ac an homozygote_count hemizygote_count 
      faf95 { popmax } 
      populations { id ac an homozygote_count hemizygote_count } 
    }
  } }`;
  try {
    const res = await fetchWithRetry('https://gnomad.broadinstitute.org/api/', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ query }) });
    if (!res.ok) throw new Error('gnomAD fail');
    const j = await res.json();
    if (gen !== window.currentSearchGen()) return;

    // Reset structured data
    data.gnomad.detailed = {
      exome: { filters: null, popmax: null, populations: {} },
      genome: { filters: null, popmax: null, populations: {} },
      total: { filters: null, popmax: null, populations: {} }
    };

    if (j.data && j.data.variant) {
      const v = j.data.variant;

      const processSrc = (src, type) => {
        if (!src) return;
        data.gnomad.detailed[type].filters = src.filters || [];
        data.gnomad.detailed[type].popmax = src.faf95?.popmax ?? null;
        data.gnomad.detailed[type].ac = src.ac || 0;
        data.gnomad.detailed[type].an = src.an || 0;
        data.gnomad.detailed[type].hom = src.homozygote_count || 0;
        data.gnomad.detailed[type].hemi = src.hemizygote_count || 0;

        // Initialize 'overall' group with global totals for this source (exome or genome)
        data.gnomad.detailed[type].populations['overall'] = {
          all: {
            ac: src.ac || 0,
            an: src.an || 0,
            hom: src.homozygote_count || 0,
            hemi: src.hemizygote_count || 0
          },
          XX: { ac: 0, an: 0, hom: 0, hemi: 0 },
          XY: { ac: 0, an: 0, hom: 0, hemi: 0 }
        };

        if (src.populations) {
          src.populations.forEach(p => {
            // Exclude only the synthetic 'global' aggregate (remaining IS stored so it feeds the true total)
            if (['global'].includes(p.id)) return;
            // Skip empty/missing ids — gnomAD occasionally returns an unlabelled entry carrying the
            // global totals, which would otherwise create a phantom "" ancestry group and inflate the
            // exome/genome overall.all recompute below.
            if (!p.id) return;

            // Determine ancestral group and sex
            let group = p.id;
            let sex = 'all';
            if (p.id.includes('_')) {
              const parts = p.id.split('_');
              group = parts[0];
              sex = parts[1]; // XX or XY
            }
            if (['XX', 'XY'].includes(group)) {
              sex = group;
              group = 'overall';
            }

            if (!data.gnomad.detailed[type].populations[group]) {
              data.gnomad.detailed[type].populations[group] = {
                all: { ac: 0, an: 0, hom: 0, hemi: 0 },
                XX: { ac: 0, an: 0, hom: 0, hemi: 0 },
                XY: { ac: 0, an: 0, hom: 0, hemi: 0 }
              };
            }
            const target = data.gnomad.detailed[type].populations[group][sex];
            target.ac = p.ac || 0;
            target.an = p.an || 0;
            target.hom = p.homozygote_count || 0;
            target.hemi = p.hemizygote_count || 0;
          });

          // ── POST-PROCESSING: Recalculate 'overall.all' from the sum of ALL individual
          // ethnic populations (including hidden groups) for manual aggregations only.
          // Rule: Skip this override if processing authoritative 'total' joint data.
          if (type !== 'total') {
            const ethnicKeys = Object.keys(data.gnomad.detailed[type].populations).filter(k => k !== 'overall');
            if (ethnicKeys.length > 0) {
              const overall = data.gnomad.detailed[type].populations['overall'];
              let sumAC = 0, sumAN = 0, sumHom = 0, sumHemi = 0;
              ethnicKeys.forEach(k => {
                const p = data.gnomad.detailed[type].populations[k];
                sumAC += p.all.ac || 0;
                sumAN += p.all.an || 0;
                sumHom += p.all.hom || 0;
                sumHemi += p.all.hemi || 0;
              });
              if (sumAN > (overall.all.an || 0)) {
                overall.all = { ac: sumAC, an: sumAN, hom: sumHom, hemi: sumHemi };
              }
              ['XX', 'XY'].forEach(sex => {
                if ((overall[sex]?.an || 0) === 0) {
                  let sexAC = 0, sexAN = 0;
                  ethnicKeys.forEach(k => {
                    sexAC += data.gnomad.detailed[type].populations[k][sex]?.ac || 0;
                    sexAN += data.gnomad.detailed[type].populations[k][sex]?.an || 0;
                  });
                  if (sexAN > 0) overall[sex] = { ac: sexAC, an: sexAN, hom: 0, hemi: 0 };
                }
              });
            }
          }
        }
      };

      processSrc(v.exome, 'exome');
      processSrc(v.genome, 'genome');

      if (v.joint) {
        processSrc(v.joint, 'total');
      } else {
        // Fallback block if joint is unavailable: manual exome + genome summation loop.
        const pops = new Set([
          ...Object.keys(data.gnomad.detailed.exome.populations),
          ...Object.keys(data.gnomad.detailed.genome.populations)
        ]);

        pops.forEach(group => {
          data.gnomad.detailed.total.populations[group] = {
            all: { ac: 0, an: 0, hom: 0, hemi: 0 },
            XX: { ac: 0, an: 0, hom: 0, hemi: 0 },
            XY: { ac: 0, an: 0, hom: 0, hemi: 0 }
          };
          ['all', 'XX', 'XY'].forEach(sex => {
            const e = data.gnomad.detailed.exome.populations[group]?.[sex] || { ac: 0, an: 0, hom: 0, hemi: 0 };
            const g = data.gnomad.detailed.genome.populations[group]?.[sex] || { ac: 0, an: 0, hom: 0, hemi: 0 };
            const t = data.gnomad.detailed.total.populations[group][sex];
            t.ac = e.ac + g.ac;
            t.an = e.an + g.an;
            t.hom = e.hom + g.hom;
            t.hemi = e.hemi + g.hemi;
          });
        });
      }

      // ACMG Safety Fallback for Popmax
      const jMax = data.gnomad.detailed.total.popmax;
      const eMax = data.gnomad.detailed.exome.popmax;
      const gMax = data.gnomad.detailed.genome.popmax;

      if (jMax !== null && jMax !== undefined) {
        data.gnomad.popmax = jMax;
      } else if (eMax === null && gMax === null) {
        data.gnomad.popmax = null;
      } else {
        data.gnomad.popmax = Math.max(eMax || 0, gMax || 0);
      }
    } else {
      data.gnomad.popmax = null;
    }

    setAPIStatus('gnomad', 'ok');
    renderMetricPanels();
    evaluateACMG();
    // BS2 Pathway A depends on the homozygote/hemizygote counts just populated above.
    if (typeof maybeEvaluateBS2 === 'function') maybeEvaluateBS2(gen);
  } catch (e) {
    if (e.name === 'AbortError') return;
    console.error(e);
    setAPIStatus('gnomad', 'error');
  }
}

/**
 * fetchGnomadCooccurrence(v1Id, v2Id, gen)
 * Inputs: two GRCh37 variant ids (chrom-pos-ref-alt, no 'chr'). Passed explicitly (Rule 2).
 * Outputs: returns { pTrans, individualCount, isSingleton, raw } or null — does NOT mutate data.*.
 * Side effects: none (caller owns data.bs2.*).
 * Failure: throws on network/GraphQL error (caller catches and surfaces in data.bs2.error, Rule 3).
 *
 * gnomAD variant co-occurrence is available ONLY for gnomAD v2 exomes (dataset gnomad_r2_1).
 * genotype_counts is a 9-element array ordered by (variant1 genotype × variant2 genotype):
 *   [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]  (A/B = ref, a/b = alt)
 * The inferred-compound-het count is the DOUBLE-HETEROZYGOUS cell AaBb (index 4) — the genotype
 * the p_compound_heterozygous (phase) estimate actually applies to. Cells where an individual is
 * homozygous-alt for one variant are a different (non-compound-het) genotype and are excluded.
 */
async function fetchGnomadCooccurrence(v1Id, v2Id, gen) {
  if (!v1Id || !v2Id) return null;
  const query = `query {
    variant_cooccurrence(variants: ["${v1Id}", "${v2Id}"], dataset: gnomad_r2_1) {
      variant_ids
      genotype_counts
      p_compound_heterozygous
    }
  }`;
  // gnomAD endpoint is rate-limited; mirror fetchGnomAD's fetchWithRetry usage (Rule 11).
  const res = await fetchWithRetry('https://gnomad.broadinstitute.org/api/', {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ query })
  });
  if (!res.ok) throw new Error('gnomAD co-occurrence HTTP ' + res.status);
  const j = await res.json();
  if (typeof gen === 'number' && gen !== window.currentSearchGen()) return null;
  if (j.errors?.length) throw new Error(j.errors[0]?.message || 'gnomAD co-occurrence error');

  const co = j.data?.variant_cooccurrence;
  if (!co) return null;

  const gc = co.genotype_counts || [];
  // Double-heterozygous individuals (AaBb, index 4) — the inferred-compound-het candidates.
  const individualCount = gc[4] || 0;
  const pTrans = (typeof co.p_compound_heterozygous === 'number') ? co.p_compound_heterozygous : undefined;

  return { pTrans, individualCount, isSingleton: individualCount === 1, raw: co };
}

/**
 * fetchVariantInGnomadV2Exomes(v2Id)
 * Inputs: a GRCh37 variant id (chrom-pos-ref-alt, no 'chr'). Passed explicitly (Rule 2).
 * Output: boolean — true if the variant is present in gnomAD v2 exomes (a precondition for
 *         co-occurrence, which is v2-exomes-only). A "variant not found" GraphQL error → false.
 * Used by the BS2 co-occurrence screen to gate work: if the TARGET isn't in v2, no pair can co-occur.
 */
async function fetchVariantInGnomadV2Exomes(v2Id) {
  if (!v2Id) return false;
  const query = `query { variant(variantId: "${v2Id}", dataset: gnomad_r2_1) { exome { ac } } }`;
  const res = await fetchWithRetry('https://gnomad.broadinstitute.org/api/', {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ query })
  });
  if (!res.ok) throw new Error('gnomAD v2 variant HTTP ' + res.status);
  const j = await res.json();
  if (j.errors?.length) return false;            // "Variant not found in gnomAD" → not in v2
  return (j.data?.variant?.exome?.ac || 0) > 0;
}

/**
 * fetchGeneClinVarPLP(gene, excludeV2Id)
 * Inputs: gene symbol; the target's own GRCh37 id to exclude. Passed explicitly (Rule 2). No data.* reads.
 * Output: array of candidate partner variants present in gnomAD v2, classified P/LP, sorted
 *         highest-confidence first: [{ variantId (GRCh37), classification:'P'|'LP', goldStars, consequence }].
 * Source: gnomAD gene(GRCh37).clinvar_variants — one call returns ClinVar variants overlapping the gene
 *         WITH their GRCh37 id, clinical significance, ClinVar gold stars, gnomAD presence, and consequence.
 * Failure: throws on network/GraphQL error (caller catches and surfaces it; Rule 3).
 *
 * Only `in_gnomad` variants are kept — co-occurrence is computable only when both variants are in gnomAD.
 * Conflicting / benign / uncertain significances are excluded; only Pathogenic / Likely pathogenic anchor BS2.
 */
async function fetchGeneClinVarPLP(gene, excludeV2Id) {
  if (!gene || gene === '-') return [];
  const query = `query { gene(gene_symbol: "${gene}", reference_genome: GRCh37) {
    clinvar_variants { variant_id clinical_significance gold_stars in_gnomad major_consequence }
  } }`;
  const res = await fetchWithRetry('https://gnomad.broadinstitute.org/api/', {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ query })
  });
  if (!res.ok) throw new Error('gnomAD gene ClinVar HTTP ' + res.status);
  const j = await res.json();
  if (j.errors?.length) throw new Error(j.errors[0]?.message || 'gnomAD gene ClinVar error');

  const list = j.data?.gene?.clinvar_variants || [];
  const out = [];
  for (const v of list) {
    if (!v.in_gnomad) continue;                       // must be in gnomAD v2 for co-occurrence
    if (excludeV2Id && v.variant_id === excludeV2Id) continue;  // not in trans with itself
    const sig = (v.clinical_significance || '').toLowerCase();
    if (sig.includes('conflict') || sig.includes('benign') || sig.includes('uncertain')) continue;
    let cls = null;
    if (sig.includes('pathogenic/likely pathogenic')) cls = 'P';
    else if (sig.includes('likely pathogenic')) cls = 'LP';
    else if (sig.includes('pathogenic')) cls = 'P';   // "Pathogenic" (alone or in combos)
    if (!cls) continue;
    out.push({ variantId: v.variant_id, classification: cls, goldStars: v.gold_stars || 0, consequence: v.major_consequence || '' });
  }
  // Highest-confidence partners first: P before LP, then more ClinVar gold stars.
  out.sort((a, b) => ((b.classification === 'P') - (a.classification === 'P')) || (b.goldStars - a.goldStars));
  return out;
}

/**
 * Performs a zero-latency synchronous lookup for ClinGen scores.
 * Returns defaults if the gene is not curated in the local dictionary.
 */
function getClinGenDosage(gene) {
  const dict = window.clingenDosageDict;
  if (!dict || !gene) return { hi: "—", ts: "—" };
  const match = dict[gene.trim().toUpperCase()];
  if (match) return { hi: match.hi, ts: match.ts };
  return { hi: "Not Curated", ts: "Not Curated" };
}

/**
 * fetchDosageSensitivity(chrom, pos38, gene)
 * Inputs: coords.chrom, coords.pos38 (1-based VCF), coords.gene
 * Outputs: data.dosage.haplo, data.dosage.triplo, data.dosage.source ('ucsc' | 'local')
 * Side effects: re-renders metric panels (eHaploScore, eTriploScore, eDosageSource)
 * Failure: falls back to local dict silently; logs warning
 */
async function fetchDosageSensitivity(chrom, pos38, gene) {
  if (!chrom || !pos38) return;
  const chromStr = `chr${chrom}`;
  const start = parseInt(pos38, 10) - 1;
  const end = start + 1;
  const base = `https://api.genome.ucsc.edu/getData/track?genome=hg38&chrom=${chromStr}&start=${start}&end=${end}`;

  setAPIStatus('ucsc', 'loading');
  try {
    const [hRes, tRes] = await Promise.all([
      fetchWithTimeout(`${base}&track=clinGenHaplo`),
      fetchWithTimeout(`${base}&track=clinGenTriplo`)
    ]);
    if (!hRes.ok || !tRes.ok) throw new Error(`UCSC HTTP ${hRes.status}/${tRes.status}`);
    const [hJson, tJson] = await Promise.all([hRes.json(), tRes.json()]);

    const haploRecs = hJson.clinGenHaplo ?? [];
    const triploRecs = tJson.clinGenTriplo ?? [];
    data.dosage.haplo = haploRecs.length ? (String(haploRecs[0].haploScore ?? '') || 'Not Curated') : 'Not Curated';
    data.dosage.triplo = triploRecs.length ? (String(triploRecs[0].triploScore ?? '') || 'Not Curated') : 'Not Curated';
    data.dosage.haploId = haploRecs[0]?.name || haploRecs[0]?.id || null;
    data.dosage.source = 'ucsc';
    setAPIStatus('ucsc', 'ok');
  } catch (e) {
    console.warn('[fetchDosageSensitivity] UCSC failed, using local:', e.message);
    if (gene) {
      const local = getClinGenDosage(gene);
      data.dosage.haplo = local.hi;
      data.dosage.triplo = local.ts;
    }
    data.dosage.source = 'local';
    setAPIStatus('ucsc', 'warn');
  }
  renderMetricPanels();
}

async function fetchGnomADConstraint(gen) {
  if (!data.coords.gene) return;
  const gConstraintEl = document.getElementById('gnomadGeneConstraint');
  if (gConstraintEl) gConstraintEl.innerText = 'Loading...';
  const query = `{
    hg38: gene(gene_symbol: "${data.coords.gene}", reference_genome: GRCh38) {
      gnomad_constraint { pLI oe_lof_upper mis_z }
    }
    hg19: gene(gene_symbol: "${data.coords.gene}", reference_genome: GRCh37) {
      gnomad_v2_regional_missense_constraint {
        passed_qc
        regions { start stop obs_exp p_value }
      }
    }
  }`;
  try {
    const res = await fetchWithRetry('https://gnomad.broadinstitute.org/api/', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ query }) });
    if (!res.ok) throw new Error('gnomAD constraint fail');
    const j = await res.json();
    if (j.errors) console.warn('gnomAD GraphQL Errors:', j.errors);
    const hg38 = j.data?.hg38?.gnomad_constraint;
    if (hg38) {
      data.gnomad.geneConstraint = {
        pli: hg38.pLI !== null ? hg38.pLI.toFixed(2) : '—',
        loeuf: hg38.oe_lof_upper !== null ? hg38.oe_lof_upper.toFixed(2) : '—',
        mis_z: hg38.mis_z !== null && hg38.mis_z !== undefined ? hg38.mis_z.toFixed(2) : '—'
      };
    }
    const hg19 = j.data?.hg19?.gnomad_v2_regional_missense_constraint;
    if (hg19 && hg19.passed_qc && hg19.regions) {
      data.rawGnomadRegionalData = hg19.regions;
      if (data.coords.pos19) {
        const numericPos = parseInt(data.coords.pos19, 10);
        const match = hg19.regions.find(r => {
          const minPos = Math.min(r.start, r.stop);
          const maxPos = Math.max(r.start, r.stop);
          return numericPos >= minPos && numericPos <= maxPos;
        });

        if (match) {
          data.gnomad.regionalOE = match.obs_exp;
          data.gnomad.regionalPValue = match.p_value ?? null;
          const minPos = Math.min(match.start, match.stop);
          const maxPos = Math.max(match.start, match.stop);
          data.gnomad.regionalRange = `${minPos}-${maxPos}`;
        } else {
          data.gnomad.regionalOE = null;
          data.gnomad.regionalPValue = null;
          data.gnomad.regionalRange = "Not in constrained interval";
        }
      }
    } else {
      data.rawGnomadRegionalData = null;
      data.gnomad.regionalOE = null;
      data.gnomad.regionalPValue = null;
      data.gnomad.regionalRange = "Not in constrained interval";
      data.gnomad.regionalConstraintMsg = (hg19 && hg19.passed_qc === false) ? 'Did not pass QC' : 'N/A';
    }
    renderMetricPanels();
  } catch (e) {
    if (e.name === 'AbortError') return;
    console.error(e);
  }
}

async function fetchUCSC(gen) {
  if (!data.coords.chrom || !data.coords.pos38 || !data.coords.ref || !data.coords.alt) return;
  if (data.coords.ref.length > 1 || data.coords.alt.length > 1) { document.getElementById('mvRevel').innerText = 'N/A (Indel)'; return; }
  const track = `revel${data.coords.alt}`;
  const start = parseInt(data.coords.pos38) - 1;
  const url = `https://api.genome.ucsc.edu/getData/track?genome=hg38&track=${track}&chrom=chr${data.coords.chrom}&start=${start}&end=${data.coords.pos38}`;
  try {
    const res = await fetch(url); if (!res.ok) throw new Error('UCSC fail');
    const d = await res.json();
    if (d[track] && d[track].length > 0) {
      const ucscVal = parseFloat(d[track][0].value);
      if (!data.scores.revel || ucscVal > data.scores.revel) {
        data.scores.revel = ucscVal;
      }
    }
    renderMetricPanels();
  } catch (e) { console.error(e); }
  evaluateACMG();
}

// Generalized multi-transcript/multi-gene result picker.
//
// APIs that query by genomic position (VEP, SpliceAI, MyVariant.info) return results
// for EVERY gene and transcript overlapping that locus. The first entry is NOT
// guaranteed to belong to the gene of interest. This function applies a consistent
// 5-level priority to select the most relevant entry regardless of API shape.
//
// Priority (most to least specific):
//   1. Exact transcript ID match  — unambiguous even across overlapping genes
//   2. MANE Select of target gene — canonical + gene name both agree
//   3. Any entry of target gene   — right gene, non-canonical transcript
//   4. Any canonical/MANE entry  — canonical but possibly different gene
//   5. list[0]                   — last resort
//
// Callers supply field extractors matching their API's response shape:
//   getGene(r)     → gene symbol string (or null if the field doesn't exist)
//   getTxId(r)     → transcript ID string (or null)
//   isCanonical(r) → boolean — true for MANE Select / canonical entries
//
// Version suffixes (.8, .3, etc.) are stripped before transcript ID comparison.
//
// @param {Array}    list         Array of API result objects
// @param {object}   opts
// @param {string}   opts.gene         Gene symbol to prefer
// @param {string}   opts.transcript   Transcript ID to prefer
// @param {Function} opts.getGene      Extractor: r → gene symbol
// @param {Function} opts.getTxId      Extractor: r → transcript ID
// @param {Function} opts.isCanonical  Extractor: r → boolean
function selectByTranscript(list, { gene, transcript, getGene, getTxId, isCanonical } = {}) {
  if (!list?.length) return null;
  const geneUp = (gene || '').toUpperCase();
  const txBase = (transcript || '').split('.')[0].toUpperCase();

  const byTx        = txBase  && list.find(r => getTxId?.(r)?.split('.')[0]?.toUpperCase() === txBase);
  const byMsGene    = geneUp  && list.find(r => isCanonical?.(r) && (getGene?.(r) || '').toUpperCase() === geneUp);
  const byGene      = geneUp  && list.find(r => (getGene?.(r) || '').toUpperCase() === geneUp);
  const byCanonical =            list.find(r => isCanonical?.(r));

  return byTx || byMsGene || byGene || byCanonical || list[0];
}

// SpliceAI-specific wrapper — Broad API uses its own field names (t_id, g_name, t_priority).
// @param {Array}  scoresList   json.scores from the Broad API
// @param {string} gene         gene symbol to prefer
// @param {string} transcriptId transcript ID to prefer
function pickSpliceAIScore(scoresList, gene, transcriptId) {
  return selectByTranscript(scoresList, {
    gene,
    transcript:  transcriptId,
    getGene:     s => s.g_name,
    getTxId:     s => s.t_id,
    isCanonical: s => s.t_priority === 'MS',
  });
}

/**
 * markSpliceAIUnavailable(gen)  — terminal failure state for the SpliceAI panel (CLAUDE.md Rule 9:
 * the SpliceAI grid owns its own UI). Replaces the transient "Calculating…" message so a failed,
 * timed-out, or rate-limited (429) fetch never leaves the panel frozen. No-ops when (a) the search is
 * stale — an AbortError from a newer search must not clobber it — or (b) a score has since arrived
 * (e.g. from the MyVariant precomputed fallback), in which case renderMetricPanels owns the message.
 */
function markSpliceAIUnavailable(gen) {
  if (gen !== window.currentSearchGen()) return;   // stale search → leave the newer search's UI alone
  if (data.scores.spliceAI !== null) return;       // a score arrived → renderMetricPanels owns the text
  const consEl = document.getElementById('spliceAiConsequence');
  if (consEl) consEl.innerText = 'Unavailable — re-run or use ↗ lookup';
  setDot('sai', 'warn');
}

async function fetchBroadSpliceAI(gen) {
  const consEl = document.getElementById('spliceAiConsequence');
  if (!data.coords.chrom || !data.coords.pos38 || !data.coords.ref || !data.coords.alt) {
    if (consEl) consEl.innerText = '';
    return;
  }
  setDot('sai', 'loading');

  // Format Variant String
  const variantString = `chr${data.coords.chrom}-${data.coords.pos38}-${data.coords.ref}-${data.coords.alt}`;
  const url = `https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&distance=500&mask=0&variant=${variantString}`;

  // UI status update for cold starts
  if (consEl) consEl.innerText = "Calculating (Model waking up, may take up to 90s)...";

  const controller = new AbortController();
  const timeoutId = setTimeout(() => controller.abort(), 90000);

  try {
    const res = await fetch(url, { signal: controller.signal });
    clearTimeout(timeoutId);

    if (res.status === 429) {
      console.warn("Broad SpliceAI Rate Limited");
      await fetchMyVariant(gen);          // precomputed fallback (SNVs only)
      markSpliceAIUnavailable(gen);        // if the fallback had no splice score, resolve "Calculating…"
      return;
    }

    if (!res.ok) throw new Error(`Broad API HTTP ${res.status}`);

    const json = await res.json();
    if (gen !== window.currentSearchGen()) return;
    if (gen !== window.currentSearchGen()) return;
    if (!json || !json.scores || json.scores.length === 0) throw new Error("No scores returned from Broad API");

    const scores = pickSpliceAIScore(json.scores, data.coords.gene, data.ensembl.transcript);

    data.scores.spliceAI_AG = parseFloat(scores.DS_AG) || 0;
    data.scores.spliceAI_AL = parseFloat(scores.DS_AL) || 0;
    data.scores.spliceAI_DG = parseFloat(scores.DS_DG) || 0;
    data.scores.spliceAI_DL = parseFloat(scores.DS_DL) || 0;

    data.scores.spliceAI_DP_AG = parseInt(scores.DP_AG) || null;
    data.scores.spliceAI_DP_AL = parseInt(scores.DP_AL) || null;
    data.scores.spliceAI_DP_DG = parseInt(scores.DP_DG) || null;
    data.scores.spliceAI_DP_DL = parseInt(scores.DP_DL) || null;

    data.scores.spliceAiManeMatch = true;

    const spScores = [data.scores.spliceAI_AG, data.scores.spliceAI_AL, data.scores.spliceAI_DG, data.scores.spliceAI_DL].filter(v => !isNaN(v));
    data.scores.spliceAI = spScores.length ? Math.max(...spScores) : null;

    // We still need to call fetchMyVariant to get ClinGen/rsId backfilling if needed, but the prompt says 
    // "falling back to myvariant.info only if it fails or hits a rate limit." So we skip it here.
    // If you need the other MyVariant data even on Broad success, you would call it here in the background.

    renderMetricPanels();
    evaluateACMG();
    setDot('sai', 'ok');
  } catch (e) {
    clearTimeout(timeoutId);
    // AbortError here is the 90 s timeout OR a newer search superseding this one. Don't early-return:
    // markSpliceAIUnavailable's gen guard already skips the stale-search case, while a genuine timeout
    // on the current search must resolve the "Calculating…" message instead of freezing it.
    if (e.name !== 'AbortError') console.error("Broad SpliceAI Error:", e);
    markSpliceAIUnavailable(gen);
  }
}

async function fetchMyVariant(gen) {
  const mvEl = document.getElementById('mvSpliceAI');
  const consEl = document.getElementById('spliceAiConsequence');

  let mvQuery = null;
  if (data.ensembl.rsId) {
    mvQuery = data.ensembl.rsId;
  } else if (data.coords.pos38 && data.coords.chrom && data.coords.ref && data.coords.alt) {
    mvQuery = `chr${data.coords.chrom}:g.${data.coords.pos38}${data.coords.ref}>${data.coords.alt}`;
  }

  if (!mvQuery) {
    if (mvEl) mvEl.innerText = 'Unavailable';
    if (consEl) consEl.innerText = '';
    return;
  }

  // Always use the /query endpoint — it handles both rsIDs and genomic HGVS without 404s.
  // The /variant/ endpoint requires the exact canonical _id which can differ from the HGVS coords.
  const FIELDS = 'dbnsfp.spliceai,dbnsfp.spliceai_pred,dbnsfp.refseq_id,dbnsfp.spliceai_coords,dbnsfp.rsid,dbnsfp.alphamissense,dbnsfp.alphamissense_pred,dbnsfp.alphamissense_score,dbnsfp.revel_score,dbnsfp.revel.score,dbsnp.rsid';
  const url = `https://myvariant.info/v1/query?q=${encodeURIComponent(mvQuery)}&fields=${FIELDS}&size=10`;

  try {
    console.log('[fetchMyVariant] Querying:', url);
    const res = await fetchWithRetry(url);
    if (!res.ok) {
      console.error('[fetchMyVariant] HTTP', res.status, 'for URL:', url);
      throw new Error(`MV fail (HTTP ${res.status})`);
    }
    const mv = await res.json();
    if (gen !== window.currentSearchGen()) return;
    // /query endpoint wraps results in { hits: [...] }
    const mvArr = mv.hits ? mv.hits : (Array.isArray(mv) ? mv : [mv]);
    let rec;
    const targetAlt = data.coords.alt ? data.coords.alt.toUpperCase() : null;
    if (targetAlt) {
      // _id format: "chr1:g.102878042C>A" — match the alt allele exactly
      rec = mvArr.find(h => {
        const id = (h._id || '').toUpperCase();
        const gtIdx = id.lastIndexOf('>');
        return gtIdx !== -1 && id.slice(gtIdx + 1) === targetAlt && h.dbnsfp;
      });
    }
    // Fallback: transcript-aware selection from records with dbnsfp data
    if (!rec) {
      const pool = mvArr.filter(h => h.dbnsfp);
      rec = selectByTranscript(pool.length ? pool : mvArr, {
        gene:        data.coords.gene,
        transcript:  data.ensembl.transcript,
        getGene:     h => (Array.isArray(h.dbnsfp?.genename) ? h.dbnsfp.genename[0] : h.dbnsfp?.genename),
        getTxId:     h => (Array.isArray(h.dbnsfp?.refseq_id) ? h.dbnsfp.refseq_id[0] : h.dbnsfp?.refseq_id),
        isCanonical: h => !!h.dbnsfp?.mane,
      });
    }

    if (!rec || !rec.dbnsfp) {
      if (mvEl) mvEl.innerText = 'Unavailable';
      if (consEl) consEl.innerText = '';
      console.warn('[fetchMyVariant] No dbnsfp found. Full response:', mv);
      return;
    }

    // Capture raw MyVariant values for two-source validation (independent of mutation logic below)
    if (data.myvariant) {
      // _id format: "chr1:g.102878042C>A" → "chr1-102878042-C-A"
      const idMatch = (rec._id || '').match(/chr([^:]+):g\.(\d+)([ACGT-]+)>([ACGT-]+)/i);
      if (idMatch) data.myvariant.hg38String = `chr${idMatch[1]}-${idMatch[2]}-${idMatch[3].toUpperCase()}-${idMatch[4].toUpperCase()}`;
      // hg19 from dbsnp.hg19 / dbnsfp.hg19_pos (if available)
      if (rec.dbnsfp?.hg19?.start && rec.dbnsfp?.hg19?.end != null && data.coords?.chrom && data.coords?.ref && data.coords?.alt) {
        data.myvariant.hg19String = `chr${data.coords.chrom}-${rec.dbnsfp.hg19.start}-${data.coords.ref}-${data.coords.alt}`;
      }
      const mvRsId = rec.dbsnp?.rsid || rec.dbnsfp?.rsid || null;
      if (mvRsId) data.myvariant.rsId = mvRsId.startsWith('rs') ? mvRsId : `rs${mvRsId}`;
    }

    // Identity Back-filling Logic
    if (!data.ensembl.rsId || data.ensembl.rsId === '—') {
      let foundId = rec.dbsnp?.rsid || rec.dbnsfp?.rsid;
      if (!foundId && rec._id && rec._id.startsWith('rs')) foundId = rec._id;

      if (foundId) {
        const cleanId = foundId.startsWith('rs') ? foundId : 'rs' + foundId;
        data.ensembl.rsId = cleanId;
        const eRsId = document.getElementById('eRsId');
        if (eRsId) eRsId.innerText = cleanId;

        // Trigger secondary searches
        if (data.ensembl.vepAminoAcids) fetchClinVarCodon(gen);
      }
    }

    const db = rec.dbnsfp;
    const targetMane = document.getElementById('eManeTranscript').innerText.trim().split('.')[0];
    let idx = -1;

    if (db.refseq_id) {
      const ids = Array.isArray(db.refseq_id) ? db.refseq_id : [db.refseq_id];
      idx = ids.findIndex(id => id && id.startsWith(targetMane));
    }

    const getIdxVal = (k, i) => {
      const sp = db.spliceai || db.spliceai_pred;
      const sc = db.spliceai_coords || {};
      let source = (k.startsWith('dp_')) ? sc : sp;
      if (!source) return null;
      let v = source[k];
      if (Array.isArray(v)) {
        if (i !== -1 && v[i] !== undefined) return parseFloat(v[i]);
        return Math.max(...v.map(Number));
      }
      return (v != null && !isNaN(v)) ? parseFloat(v) : null;
    };

    if (idx === -1) data.scores.spliceAiManeMatch = false;
    else data.scores.spliceAiManeMatch = true;

    // Safety: Only backfill SpliceAI if Broad fetch hasn't already populated it with high-confidence data
    if (data.scores.spliceAI === null) {
      data.scores.spliceAI_AG = getIdxVal('ds_ag', idx);
      data.scores.spliceAI_AL = getIdxVal('ds_al', idx);
      data.scores.spliceAI_DG = getIdxVal('ds_dg', idx);
      data.scores.spliceAI_DL = getIdxVal('ds_dl', idx);
      data.scores.spliceAI_DP_AG = getIdxVal('dp_ag', idx);
      data.scores.spliceAI_DP_AL = getIdxVal('dp_al', idx);
      data.scores.spliceAI_DP_DG = getIdxVal('dp_dg', idx);
      data.scores.spliceAI_DP_DL = getIdxVal('dp_dl', idx);

      const spScores = [data.scores.spliceAI_AG, data.scores.spliceAI_AL, data.scores.spliceAI_DG, data.scores.spliceAI_DL].filter(v => v !== null);
      data.scores.spliceAI = spScores.length ? Math.max(...spScores) : null;
    }

    // ── AlphaMissense Extraction (Fallback: Only if Ensembl VEP failed) ──
    if (data.scores.alphaMissenseScore === null) {
      const fullManeTxFallback = document.getElementById('eManeTranscript').innerText.trim();
      data.scores.alphaMissenseScore = getBestAlphaMissense(rec, fullManeTxFallback);
      data.scores.alphaMissensePred = getBestAlphaMissensePred(rec, fullManeTxFallback);
    }

    // ── REVEL Extraction (Strict MANE Match with Source Fallback) ──
    const fullManeTx = document.getElementById('eManeTranscript').innerText.trim();
    const bestScore = getBestRevel(rec, fullManeTx);
    if (bestScore !== null && !isNaN(bestScore)) {
      data.scores.revel = bestScore;
      data.revelIsManeMatch = true;
    }

    renderMetricPanels();
    evaluateACMG();
    setDot('sai', 'ok');
    runValidationPass();
  } catch (e) {
    if (e.name === 'AbortError') return;
    console.error('fetchMyVariant:', e);
    if (mvEl) mvEl.innerText = 'Unavailable';
    if (consEl) consEl.innerText = '';
    setDot('sai', 'warn');
  }
}

async function fetchLitVar(gen) {
  // Anchor priority (transcript-agnostic first):
  //   1. rsID          — transcript-independent, best LitVar match
  //   2. VCF genomic   — transcript-independent, immutable for novel variants
  //   3. GENE + c.HGVS — transcript-specific fallback only
  let query;
  if (data.ensembl.rsId && data.ensembl.rsId !== '—') {
    query = data.ensembl.rsId;
  } else if (data.coords.hg38String) {
    query = data.coords.hg38String; // e.g. chr10-69005731-G-T
  } else if (data.coords.gene && data.coords.hgvs && data.coords.hgvs !== '-') {
    query = `${data.coords.gene} ${data.coords.hgvs}`;
  } else {
    document.getElementById('pLitVarMatches').innerHTML = '—';
    return;
  }

  setAPIStatus('litvar', 'loading');
  const targetUrl = `https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1/entity/litvar/${encodeURIComponent(query)}`;
  // Use CodeTabs proxy, which safely supports file:// (null) origins for local development
  const url = `https://api.codetabs.com/v1/proxy?quest=${encodeURIComponent(targetUrl)}`;

  try {
    const res = await fetch(url);
    if (!res.ok) throw new Error(`LitVar fetch failed with status ${res.status}`);

    const text = await res.text();
    let json = null;
    if (text && text.trim()) {
      try { json = JSON.parse(text); } catch { json = null; }
    }

    let litvarEntity = null;
    if (Array.isArray(json) && json.length > 0) {
      litvarEntity = json.find(j => j.id === query || j.name === query)
                  || json.find(j => (j.gene || '').toUpperCase() === (data.coords.gene || '').toUpperCase())
                  || json[0];
    } else if (json && json.pmids) {
      litvarEntity = json;
    }

    if (litvarEntity && litvarEntity.pmids && Array.isArray(litvarEntity.pmids)) {
      data.litvarPMIDs = litvarEntity.pmids;
    } else {
      data.litvarPMIDs = [];
    }

    setAPIStatus('litvar', 'ok');
    renderLitVarUI();
  } catch (e) {
    if (e.name === 'AbortError') return;
    setAPIStatus('litvar', 'error');
    console.error('fetchLitVar error:', e);
    document.getElementById('pLitVarMatches').innerHTML = '—';
  }
}

const _CV_FLAG_PATTERNS = {
  functional:  /\bfunctional (assay|stud|analys|data|evidence|experiment|test)|in vitro|biochemical|minigene|splicing assay|reporter (assay|gene)|luciferase|western blot|transfect|protein (function|expression|activity|stability)\b/i,
  deNovo:      /\bde novo\b/i,
  segregation: /\bco.?segregat|segregates? with|segregation in (the |this |affected )?family|affected (relative|sibling|sib|member|individual)s? (also )?carr|inherited from (an? )?(affected|carrier)\b/i,
};

function parseClinVarEvidenceFlags(xmlText) {
  const doc = new DOMParser().parseFromString(xmlText, 'application/xml');
  const assertions = [
    ...doc.querySelectorAll('ClinVarAssertion'),
    ...doc.querySelectorAll('ClinicalAssertion')
  ];

  const result = { functional: [], deNovo: [], segregation: [], scvDetails: [] };

  for (const a of assertions) {
    // Prefer full SubmitterName (new schema) → legacy submitter attr → OrgAbbreviation
    const fullName = (
      a.querySelector('ClinVarAccession[Type="SCV"]')?.getAttribute('SubmitterName') ||
      a.querySelector('ClinVarSubmissionID')?.getAttribute('submitter') ||
      a.querySelector('ClinVarAccession[Type="SCV"]')?.getAttribute('OrgAbbreviation') ||
      'Unknown'
    ).trim();
    // Short label for the summary pills — first comma-delimited segment
    const submitter = fullName.split(',')[0].trim() || fullName;

    const classification = (
      a.querySelector('ClinicalSignificance Description, Classification GermlineClassification, Interpretation Description')?.textContent ||
      ''
    ).trim();

    // Display comment: classification-level only (avoids ACMG code strings from observation nodes)
    const displayComment = (
      a.querySelector('ClinicalSignificance > Comment')?.textContent?.trim() ||
      a.querySelector('Classification > Comment')?.textContent?.trim() ||
      a.querySelector('Interpretation > Comment')?.textContent?.trim() ||
      ''
    );

    // Detection text: all nodes (classification + observation descriptions + observation comments)
    const allText = [
      displayComment,
      ...[...a.querySelectorAll('Attribute[Type="Description"]')].map(n => n.textContent.trim()),
      ...[...a.querySelectorAll('ObservedIn Comment, Observation Comment')].map(n => n.textContent.trim()),
    ].filter(Boolean).join(' ').toLowerCase();

    // De novo: Origin element OR keyword in comment
    const origins = [...a.querySelectorAll('Origin')].map(o => o.textContent.trim().toLowerCase());
    const isDeNovo = origins.some(o => o.includes('de novo')) || _CV_FLAG_PATTERNS.deNovo.test(allText);

    // Functional: research MethodType OR keyword
    const methodTypes = [...a.querySelectorAll('MethodType')].map(m => m.textContent.toLowerCase());
    const isFunctional = methodTypes.some(m => m === 'research' || m.includes('in vitro'))
      || _CV_FLAG_PATTERNS.functional.test(allText);

    // Segregation: FamilyData OR keyword
    const segFamilies = parseInt(a.querySelector('FamilyData')?.getAttribute('NumFamiliesWithSegregationObserved') || '0');
    const isSegregation = segFamilies > 0 || _CV_FLAG_PATTERNS.segregation.test(allText);

    if (isDeNovo)      result.deNovo.push(submitter);
    if (isFunctional)  result.functional.push(submitter);
    if (isSegregation) result.segregation.push(submitter);

    const flags = [
      ...(isDeNovo      ? ['deNovo']      : []),
      ...(isFunctional  ? ['functional']  : []),
      ...(isSegregation ? ['segregation'] : []),
    ];
    // Date: prefer DateLastEvaluated on classification element, fall back to DateUpdated on accession
    const date = (
      a.querySelector('ClinicalSignificance')?.getAttribute('DateLastEvaluated') ||
      a.querySelector('Classification')?.getAttribute('DateLastEvaluated') ||
      a.querySelector('ClinVarAccession[Type="SCV"]')?.getAttribute('DateUpdated') ||
      ''
    );

    // Disease: first preferred name from associated TraitSet
    const disease = (
      a.querySelector('TraitSet[Type="Disease"] Trait Name ElementValue[Type="Preferred"]')?.textContent?.trim() ||
      ''
    );

    // Only include display entry when classification comment is substantive
    const isInformative = displayComment.length > 25 && !/^not provided\.?$/i.test(displayComment.trim());
    if (flags.length > 0 && isInformative) {
      result.scvDetails.push({ submitter: fullName, classification, comment: displayComment, flags, date, disease });
    }
  }

  result.functional  = [...new Set(result.functional)];
  result.deNovo      = [...new Set(result.deNovo)];
  result.segregation = [...new Set(result.segregation)];
  return result;
}

function _highlightEvidenceSentences(comment) {
  // Split into sentences, wrap those matching a flag pattern with a colour span
  const sentences = comment.match(/[^.!?]+[.!?]*\s*/g) || [comment];
  return sentences.map(s => {
    const lower = s.toLowerCase();
    let bg = '', color = '';
    if (_CV_FLAG_PATTERNS.deNovo.test(lower)) {
      bg = 'rgba(255,77,109,.13)'; color = 'var(--red)';
    } else if (_CV_FLAG_PATTERNS.functional.test(lower)) {
      bg = 'rgba(255,179,71,.15)'; color = 'var(--amber)';
    } else if (_CV_FLAG_PATTERNS.segregation.test(lower)) {
      bg = 'rgba(0,212,170,.12)'; color = 'var(--teal)';
    }
    if (!bg) return _escHtml(s);
    return `<span style="background:${bg};border-radius:3px;padding:1px 3px;color:${color};">${_escHtml(s)}</span>`;
  }).join('');
}

function _escHtml(s) {
  return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

function _selectTopSCVs(scvDetails, flagType) {
  return scvDetails
    .filter(d => d.flags.includes(flagType))
    .sort((a, b) => {
      const aD = a.disease && !/^not (provided|specified)$/i.test(a.disease.trim()) ? 1 : 0;
      const bD = b.disease && !/^not (provided|specified)$/i.test(b.disease.trim()) ? 1 : 0;
      if (bD !== aD) return bD - aD;
      return (b.date || '').localeCompare(a.date || '');
    })
    .slice(0, 2);
}

function renderClinVarSubmitterComments(flags, uid) {
  const card = document.getElementById('clinvarEvidenceCard');
  const content = document.getElementById('clinvarEvidenceContent');
  if (!card || !content) return;

  const FLAG_META = {
    deNovo:      { label: 'De novo',         bg: 'var(--red-d)',   color: 'var(--red)'   },
    functional:  { label: 'Functional data',  bg: 'var(--amber-d)', color: 'var(--amber)' },
    segregation: { label: 'Co-segregation',   bg: 'var(--teal-d)',  color: 'var(--teal)'  },
  };

  const allDetails = flags.scvDetails || [];
  const sections = ['deNovo', 'functional', 'segregation']
    .map(flagType => ({ flagType, entries: _selectTopSCVs(allDetails, flagType), total: allDetails.filter(d => d.flags.includes(flagType)).length }))
    .filter(s => s.entries.length > 0);

  if (sections.length === 0) { card.style.display = 'none'; return; }

  const assertionsUrl = uid ? `https://www.ncbi.nlm.nih.gov/clinvar/variation/${uid}/#clinical-assertions` : '#';

  content.innerHTML = sections.map((section, si) => {
    const m = FLAG_META[section.flagType];

    const entriesHtml = section.entries.map((d, i) => {
      const sep = i > 0 ? 'margin-top:12px;padding-top:12px;border-top:1px solid rgba(255,255,255,.05);' : '';
      const dateStr = d.date ? `<span style="color:var(--dim);font-size:.65rem;margin-left:4px;">${d.date}</span>` : '';
      const diseaseStr = d.disease ? `<div style="font-size:.68rem;color:var(--dim);margin-bottom:4px;">${_escHtml(d.disease)}</div>` : '';
      return `<div style="${sep}">
        <div style="display:flex;align-items:baseline;gap:6px;flex-wrap:wrap;margin-bottom:3px;">
          <span style="font-weight:700;font-size:.78rem;color:var(--text);">${_escHtml(d.submitter)}</span>
          <span class="sig-badge ${sigClass(d.classification)}" style="font-size:.62rem;">${_escHtml(d.classification)}</span>
          ${dateStr}
        </div>
        ${diseaseStr}
        <div style="font-size:.76rem;line-height:1.65;color:var(--text);">${_highlightEvidenceSentences(d.comment)}</div>
      </div>`;
    }).join('');

    const moreHtml = section.total > 2
      ? `<div style="font-size:.65rem;color:var(--dim);margin-top:8px;">
           +${section.total - 2} more — <a href="${assertionsUrl}" target="_blank" style="color:var(--teal);">view all on ClinVar ↗</a>
         </div>`
      : '';

    const sectionSep = si > 0 ? 'margin-top:20px;padding-top:18px;border-top:1px solid var(--border);' : '';
    const countLabel = section.total > 2 ? `<span style="font-size:.65rem;color:var(--dim);">Showing 2 of ${section.total}</span>` : '';

    return `<div style="${sectionSep}">
      <div style="display:flex;align-items:center;gap:8px;margin-bottom:10px;">
        <span style="background:${m.bg};color:${m.color};border-radius:3px;padding:2px 9px;font-size:.65rem;font-weight:700;text-transform:uppercase;letter-spacing:.06em;">${m.label}</span>
        <span style="flex:1;height:1px;background:var(--border);"></span>
        ${countLabel}
      </div>
      ${entriesHtml}
      ${moreHtml}
    </div>`;
  }).join('');

  card.style.display = '';
}

function renderClinVarEvidenceFlags(flags, uid) {
  const assertionsUrl = `https://www.ncbi.nlm.nih.gov/clinvar/variation/${uid}/#clinical-assertions`;
  const viewLink = `<a href="${assertionsUrl}" target="_blank"
    style="color:var(--teal);text-decoration:none;font-size:0.68rem;font-weight:700;
           border:1px solid var(--teal);padding:1px 5px;border-radius:3px;margin-left:6px;">View ↗</a>`;

  const renderRow = (rowId, labs) => {
    const dv = document.querySelector(`#${rowId} .dv`);
    if (!dv) return;
    if (labs.length === 0) { dv.innerHTML = '—'; return; }
    const shown = labs.slice(0, 3).join(', ') + (labs.length > 3 ? ` +${labs.length - 3}` : '');
    dv.innerHTML = `<span style="color:var(--green);font-weight:700;">Yes</span>`
      + ` <span style="color:var(--dim);font-size:.68rem;">(${shown})</span>${viewLink}`;
  };

  renderRow('cvFlagDeNovo',      flags.deNovo);
  renderRow('cvFlagFunctional',  flags.functional);
  renderRow('cvFlagSegregation', flags.segregation);

  const hasAny = flags.deNovo.length > 0 || flags.functional.length > 0 || flags.segregation.length > 0;
  const container = document.getElementById('cvEvidenceFlags');
  if (container) container.style.display = hasAny ? '' : 'none';
}

/**
 * parseVCVRecordDetail(xmlText)
 * Inputs: raw ClinVar VCV XML (efetch rettype=vcv).
 * Outputs: a structured object mirroring the ClinVar Variation page —
 *   identifiers, type/length, location (cyto/GRCh38/GRCh37), timeline,
 *   embedded ClinGen dosage, RCV conditions, and SCV submissions (with
 *   comment + observation). Consumed by copyReport's ClinVar tables.
 * Failures: returns null on parse error (non-blocking).
 */
function parseVCVRecordDetail(xmlText) {
  try {
    const doc = new DOMParser().parseFromString(xmlText, 'text/xml');
    if (doc.querySelector('parsererror')) return null;
    const va = doc.querySelector('VariationArchive');
    if (!va) return null;

    const txt = (el) => el ? el.textContent.trim() : null;
    const locOf = (s) => s ? {
      assembly: s.getAttribute('Assembly'),
      chr: s.getAttribute('Chr'),
      pos: s.getAttribute('positionVCF') || s.getAttribute('start'),
      ref: s.getAttribute('referenceAlleleVCF'),
      alt: s.getAttribute('alternateAlleleVCF'),
      length: s.getAttribute('variantLength')
    } : null;

    const detail = {
      variationId: va.getAttribute('VariationID'),
      name: va.getAttribute('VariationName'),
      type: va.getAttribute('VariationType'),
      accession: va.getAttribute('Accession'),
      version: va.getAttribute('Version'),
      recordType: va.getAttribute('RecordType'),
      numSubmissions: va.getAttribute('NumberOfSubmissions'),
      numSubmitters: va.getAttribute('NumberOfSubmitters'),
      dateCreated: va.getAttribute('DateCreated'),
      dateLastUpdated: va.getAttribute('DateLastUpdated'),
      mostRecentSubmission: va.getAttribute('MostRecentSubmission')
    };

    // Variant-level location (direct child of SimpleAllele, carries positionVCF)
    const seqLocs = [...doc.querySelectorAll('ClassifiedRecord > SimpleAllele > Location > SequenceLocation')];
    detail.cytogenetic = txt(doc.querySelector('ClassifiedRecord > SimpleAllele > Location > CytogeneticLocation'));
    detail.grch38 = locOf(seqLocs.find(s => s.getAttribute('Assembly') === 'GRCh38'));
    detail.grch37 = locOf(seqLocs.find(s => s.getAttribute('Assembly') === 'GRCh37'));
    detail.length = detail.grch38?.length || detail.grch37?.length || null;

    // Aggregate germline classification + last evaluated
    const gc = doc.querySelector('ClassifiedRecord > Classifications > GermlineClassification');
    detail.classification = txt(gc?.querySelector('Description'));
    detail.reviewStatus   = txt(gc?.querySelector('ReviewStatus'));
    detail.lastEvaluated  = gc?.getAttribute('DateLastEvaluated') || null;

    // Embedded ClinGen dosage (gene-level)
    const hi = doc.querySelector('Gene > Haploinsufficiency');
    const ts = doc.querySelector('Gene > Triplosensitivity');
    detail.dosage = {
      haplo:  hi ? { text: txt(hi), lastEvaluated: hi.getAttribute('last_evaluated') } : null,
      triplo: ts ? { text: txt(ts), lastEvaluated: ts.getAttribute('last_evaluated') } : null
    };

    // RCV conditions
    detail.conditions = [...doc.querySelectorAll('RCVList > RCVAccession')].map(rcv => {
      const desc = rcv.querySelector('RCVClassifications Description');
      return {
        names: [...rcv.querySelectorAll('ClassifiedCondition')].map(c => c.textContent.trim()),
        classification: txt(desc),
        submissionCount: desc?.getAttribute('SubmissionCount'),
        lastEvaluated: desc?.getAttribute('DateLastEvaluated'),
        reviewStatus: txt(rcv.querySelector('ReviewStatus')),
        rcv: rcv.getAttribute('Accession') + '.' + rcv.getAttribute('Version')
      };
    });

    // SCV submissions (with comment + observation)
    detail.submissions = [...doc.querySelectorAll('ClinicalAssertionList > ClinicalAssertion')].map(ca => {
      const accEl = ca.querySelector('ClinVarAccession');
      const cls   = ca.querySelector('Classification');
      const obs   = ca.querySelector('ObservedInList > ObservedIn');
      const assertionMethod = [...ca.querySelectorAll('AttributeSet > Attribute')]
        .find(a => a.getAttribute('Type') === 'AssertionMethod');
      return {
        submitter: accEl?.getAttribute('SubmitterName'),
        scv: accEl ? accEl.getAttribute('Accession') + '.' + accEl.getAttribute('Version') : null,
        orgCategory: accEl?.getAttribute('OrganizationCategory'),
        dateCreated: accEl?.getAttribute('DateCreated'),
        dateUpdated: accEl?.getAttribute('DateUpdated'),
        classification: txt(cls?.querySelector('GermlineClassification, SomaticClinicalImpact, OncogenicityClassification')),
        lastEvaluated: cls?.getAttribute('DateLastEvaluated'),
        reviewStatus: txt(cls?.querySelector('ReviewStatus')),
        comment: txt(cls?.querySelector('Comment')),
        assertionMethod: txt(assertionMethod),
        conditions: [...ca.querySelectorAll('TraitSet > Trait > Name > ElementValue[Type="Preferred"]')].map(n => n.textContent.trim()),
        origin: txt(obs?.querySelector('Sample > Origin')),
        affectedStatus: txt(obs?.querySelector('Sample > AffectedStatus')),
        method: txt(obs?.querySelector('Method > MethodType'))
      };
    });

    return detail;
  } catch (e) {
    console.warn('[ClinVar] VCV detail parse failed:', e.message);
    return null;
  }
}

/**
 * ensureClinVarDetail()
 * Lazy fallback for copyReport: if the VCV detail wasn't parsed yet (user
 * copied before fetchClinVar finished) but a uid is known, re-fetch + parse.
 */
async function ensureClinVarDetail() {
  if (data.clinvar.detail || !data.clinvar.uid) return data.clinvar.detail;
  try {
    const xRes = await fetchWithRetry(`${NCBI_BASE}/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid=true&id=${data.clinvar.uid}&api_key=${NCBI_KEY}`);
    if (xRes.ok) data.clinvar.detail = parseVCVRecordDetail(await xRes.text());
  } catch (e) { console.warn('[ClinVar] detail refetch failed:', e.message); }
  return data.clinvar.detail;
}

async function fetchClinVar(gen) {
  setDot('cv', 'loading');
  // Reset evidence flags and comments card on each new fetch
  const efContainer = document.getElementById('cvEvidenceFlags');
  if (efContainer) efContainer.style.display = 'none';
  ['cvFlagDeNovo', 'cvFlagFunctional', 'cvFlagSegregation'].forEach(id => {
    const dv = document.querySelector(`#${id} .dv`);
    if (dv) dv.innerHTML = '—';
  });
  const evCard = document.getElementById('clinvarEvidenceCard');
  if (evCard) evCard.style.display = 'none';
  const evContent = document.getElementById('clinvarEvidenceContent');
  if (evContent) evContent.innerHTML = '';
  const hgvsTerm = (data.ensembl.transcript && data.coords.hgvs && data.coords.hgvs !== '-') ? `${data.ensembl.transcript}:${data.coords.hgvs}` : data.coords.hgvs;
  const term = (data.coords.gene && hgvsTerm && hgvsTerm !== '-')
    ? `(${data.coords.gene}[Gene]) AND "${hgvsTerm}"[VARNAME]`
    : data.coords.gene ? `${data.coords.gene}[Gene]` : '';

  if (!term) { setDot('cv', 'warn'); return; }

  try {
    const sUrl = `${NCBI_BASE}/esearch.fcgi?db=clinvar&term=${encodeURIComponent(term)}&retmode=json&api_key=${NCBI_KEY}&email=${NCBI_EMAIL}`;
    const sRes = await fetchWithRetry(sUrl);
    if (!sRes.ok) throw new Error(`Search fail (${sRes.status})`);

    const sJson = await ncbiSafeJson(sRes);
    if (gen !== window.currentSearchGen()) return;
    const uid = sJson.esearchresult?.idlist?.[0];

    if (!uid) {
      setDot('cv', 'warn');
      document.getElementById('cvSig').innerHTML = '<span class="sig-badge">Not found</span>';
      return;
    }
    data.clinvar.uid = uid;

    // Small delay to prevent race conditions with NCBI indexing
    await new Promise(r => setTimeout(r, 400));

    const eUrl = `${NCBI_BASE}/esummary.fcgi?db=clinvar&id=${uid}&retmode=json&api_key=${NCBI_KEY}&email=${NCBI_EMAIL}`;
    let eJson = null;

    // Retry loop for summary as it's often more unstable than search
    for (let attempt = 0; attempt < 3; attempt++) {
      const eRes = await fetchWithRetry(eUrl);
      if (eRes.ok) {
        eJson = await ncbiSafeJson(eRes);
        if (eJson.result?.[uid]) break;
      }
      await new Promise(r => setTimeout(r, 500 * (attempt + 1)));
    }

    if (!eJson || !eJson.result?.[uid]) throw new Error('Summary retrieval failed');
    if (gen !== window.currentSearchGen()) return;
    const rec = eJson.result[uid];

    let sig = rec.germline_classification?.description || rec.clinical_significance?.description || rec.variation_set?.[0]?.variation_xrefs?.[0]?.description || 'Unknown';
    const rs = rec.germline_classification?.review_status || rec.clinical_significance?.review_status || '';
    let subs = rec.germline_classification?.num_submitters || rec.num_submitters;

    if (subs === undefined && rec.supporting_submissions) {
      const scv = rec.supporting_submissions.scv ? rec.supporting_submissions.scv.length : 0;
      const rcv = rec.supporting_submissions.rcv ? rec.supporting_submissions.rcv.length : 0;
      subs = scv || rcv || 0;
    }
    subs = subs || 0;

    const nStars = reviewStars(rs);
    const acc = rec.accession_version || rec.accession || '—';
    data.clinvar.accession = acc;

    document.getElementById('cvAcc').innerHTML = `<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${uid}/" target="_blank" style="color:inherit; text-decoration:underline;">${acc}</a>`;

    if (uid) {
      try {
        const xRes = await fetchWithRetry(`${NCBI_BASE}/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid=true&id=${uid}&api_key=${NCBI_KEY}`);
        if (xRes.ok) {
          const xTxt = await xRes.text();
          if (sig.toLowerCase().includes('conflicting')) {
            const match = xTxt.match(/<(?:Explanation|Comment)[^>]*DataSource="ClinVar"[^>]*Type="public"[^>]*>([^<]+)<\/(?:Explanation|Comment)>/);
            if (match && match[1]) sig = `Conflicting: ${match[1]}`;
          }
          const flags = parseClinVarEvidenceFlags(xTxt);
          data.clinvar.evidenceFlags = flags;
          renderClinVarEvidenceFlags(flags, uid);
          renderClinVarSubmitterComments(flags, uid);
          // Full VCV record detail → data.clinvar.detail (for the analysis summary tables)
          data.clinvar.detail = parseVCVRecordDetail(xTxt);
          // Mine historical p. notations from VCV (different transcript versions)
          // and re-trigger literature search with the expanded notation set.
          const historical = parseVCVProteinNotations(xTxt);
          data.clinvar.historicalProteinNotations = historical;
          _scheduleLitNotations(gen);
        }
      } catch (e) { console.warn("ClinVar VCV Fetch failed:", e); }
    }

    data.clinvar.sig = sig; data.clinvar.stars = nStars; data.clinvar.reviewStatus = rs; data.clinvar.subs = subs;

    // Track classification for API badge
    updateAPIClassification('cv', sig);

    const renderEvidenceLink = (interpretationId) => {
      const href = interpretationId
        ? `https://erepo.clinicalgenome.org/evrepo/ui/interpretation/${interpretationId}`
        : `https://www.ncbi.nlm.nih.gov/clinvar/variation/${uid}/#clinical-assertions`;
      const title = interpretationId ? 'ClinGen Expert Panel Review' : 'Expert panel submission — view evidence';
      return ` <a href="${href}" target="_blank" id="cvEvidenceLink"
                  style="color:var(--teal);text-decoration:none;font-size:0.7rem;font-weight:700;border:1px solid var(--teal);padding:1px 5px;border-radius:3px;margin-left:6px;"
                  title="${title}">Evidence ↗</a>`;
    };

    if (nStars >= 3) {
      // Show Evidence Repository row and try to fetch UUID
      document.getElementById('cvEvidenceRepoRow').style.display = 'flex';
      tryFetchClinGenERepoLink();
    }

    document.getElementById('cvSig').innerHTML = `<span class="sig-badge ${sigClass(sig)}">${sig}</span>`;
    document.getElementById('cvStatus').innerText = rs || '—';
    const evidenceLinkHtml = nStars >= 3 ? renderEvidenceLink(null) : '';
    document.getElementById('cvStars').innerHTML = `<span class="stars">${stars(nStars)}</span> (${nStars}/4)${evidenceLinkHtml}`;
    document.getElementById('cvSubs').innerText = subs;
    setDot('cv', 'ok');

  } catch (e) {
    if (e.name === 'AbortError') return;
    console.error("fetchClinVar Error:", e);
    setDot('cv', 'error');
    document.getElementById('cvSig').innerText = 'ClinVar error: ' + (e.message || 'Unknown');
  }
}

/**
 * Fetches ClinGen Evidence Repository UUID for an expert-panel-curated variant (3+ stars).
 * If CA ID not yet available, fetches it from ClinGen Allele Registry first.
 * Then queries: https://erepo.clinicalgenome.org/evrepo/api/classifications?caId={caId}
 * Extracts uuid from response and renders link below Accession.
 */
async function tryFetchClinGenERepoLink() {
  // Ensure CA ID is resolved — reuse existing fetchCanonicalAllele if needed
  if (!data.caId && data.ensembl.transcript && data.coords.hgvs && data.coords.hgvs !== '-') {
    await fetchCanonicalAllele(`${data.ensembl.transcript}:${data.coords.hgvs}`);
  }

  const caId = data.caId;
  if (!caId) return;

  try {
    const url = `https://erepo.clinicalgenome.org/evrepo/api/classifications?caId=${encodeURIComponent(caId)}`;
    const res = await fetchWithTimeout(url);
    if (!res.ok) return;

    const json = await res.json();
    const items = Array.isArray(json) ? json : (json.results || json.classifications || []);
    const uuid = items[0]?.uuid;
    if (!uuid) return;

    data.clinvar.evidenceRepoUuid = uuid;
    const repoEl = document.getElementById('cvEvidenceRepo');
    if (repoEl) {
      const link = `https://erepo.clinicalgenome.org/evrepo/ui/interpretation/${uuid}`;
      repoEl.innerHTML = `<a href="${link}" target="_blank" style="color:inherit; text-decoration:underline;">View ↗</a>`;
    }
  } catch (e) {
    // ERepo not available for this variant — row stays as "—"
  }
}

// ── CLINGEN SVIS CACHE ─────────────────────────────────────────────
let cspecSVISCache = null;

/**
 * Bootstraps the ClinGen SVIS index in the background.
 * Prevents hitting the API repeatedly for every gene query.
 */
async function prefetchClinGenSVIS() {
  try {
    const res = await fetchWithRetry("https://cspec.genome.network/cspec/api/svis");
    const json = await res.json();
    // Accommodate standard REST arrays or wrapped { data: [] } payloads
    cspecSVISCache = Array.isArray(json) ? json : (json.data || []);
    console.log(`[ClinGen] Successfully cached ${cspecSVISCache.length} VCEP guidelines.`);
  } catch (e) {
    console.warn("[ClinGen] Failed to pre-fetch SVIS index:", e);
    cspecSVISCache = []; // Fallback to prevent infinite retries on failure
  }
}

/**
 * warmupColdStartServices()  — fire-and-forget on app load; touches NO data.* (CLAUDE.md Rules 1/5/6).
 * Pre-spins the VariantValidator backend (~30 s cold start) so that cold start doesn't land on the
 * user's first real query. The response is discarded — nothing is parsed or stored.
 *
 * Deliberately does NOT warm Broad SpliceAI: that API is rate-limited to "several requests per user
 * per minute" (per the Broad docs), so a load-time warmup would spend that tiny budget and risk
 * 429-ing the user's real query moments later — the opposite of helpful. SpliceAI warms on first use.
 *
 * Rule 11: fetchWithTimeout, never bare fetch. Rule 3: fully error-contained (silent on failure).
 * Mirrors prefetchClinGenSVIS (also a load-time prefetch, never into data.*).
 */
function warmupColdStartServices() {
  // VariantValidator has no helper proxy; hit it directly. Rule 11: 30 s timeout for VV cold starts.
  fetchWithTimeout(
    'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/NM_000088.3%3Ac.589G%3ET/mane_select',
    {}, 30000
  ).catch(() => { /* warmup is best-effort — silent */ });
}

async function fetchUniprotDomains(gene, gen) {
  if (!gene) return [];
  setAPIStatus('uniprot', 'loading');
  const TIER = {
    'Active site': 'critical', 'Metal binding': 'critical', 'Binding site': 'critical',
    'Zinc finger': 'high', 'Domain': 'high', 'Disulfide bond': 'high', 'Repeat': 'high',
    'Transmembrane': 'moderate', 'Region': 'moderate', 'Site': 'moderate', 'Coiled coil': 'moderate', 'Coil': 'moderate',
    'Motif': 'low', 'Signal peptide': 'low', 'Propeptide': 'low'
  };
  try {
    const fields = 'ft_act_site,ft_binding,ft_zn_fing,ft_domain,ft_disulfid,ft_transmem,ft_region,ft_site,ft_motif,ft_signal,ft_repeat,ft_propep,accession';
    const url = `https://rest.uniprot.org/uniprotkb/search?query=gene_exact:${encodeURIComponent(gene)}+AND+reviewed:true+AND+organism_id:9606&fields=${fields}&format=json&size=1`;
    const res = await fetchWithTimeout(url);
    if (!res.ok) {
      setAPIStatus('uniprot', 'warn');
      return [];
    }
    const json = await res.json();
    const entry = json.results?.[0];
    if (!entry) {
      setAPIStatus('uniprot', 'warn');
      return [];
    }
    data.uniprotAccession = entry.primaryAccession || null;
    const rawFeatures = entry.features || [];
    console.log(`[UniProt] ${gene} (${data.uniprotAccession}): ${rawFeatures.length} features returned`);
    // UniProt v2 API returns location as {start:{value:N}, end:{value:N}} OR {start:N, end:N}
    const getPos = loc => (typeof loc === 'object' && 'value' in loc) ? loc.value : loc;
    setAPIStatus('uniprot', 'ok');
    return rawFeatures
      .filter(f => getPos(f.location?.start) != null && getPos(f.location?.end) != null)
      .map(f => ({
        type: f.type,
        tier: TIER[f.type] || 'low',
        start: getPos(f.location.start),
        end: getPos(f.location.end),
        name: f.description || f.type
      }));
  } catch (e) {
    console.warn('[UniProt] Domain fetch failed:', e.message);
    setAPIStatus('uniprot', 'error');
    return [];
  }
}

/**
 * fetchEnsemblProteinDomains(transcriptOrGene, gene, gen)
 * Inputs:  transcript ID (ENST or NM_) and/or gene symbol; gen counter
 * Outputs: [{type:'Domain', tier:'high', start, end, name, pfamId}] — Pfam structural domains
 * Side effects: none (writes nothing to data.*)
 * Failures: silent — returns [] on any error
 *
 * Uses Ensembl /overlap/translation/{ENSP} endpoint (same source as DECIPHER protein viewer).
 * ENSP is resolved via:
 *   Path 1: ENST → /lookup/id/{ENST}  → Translation.id
 *   Path 2: NM_ or missing → /lookup/symbol/homo_sapiens/{gene}?expand=1 → canonical transcript → Translation.id
 * Pfam domains are normalised to type:'Domain' so codon-viewer's STRIP_TYPES pass them through.
 */
async function fetchEnsemblProteinDomains(transcriptOrGene, gene, gen) {
  if (!transcriptOrGene && !gene) return [];
  try {
    let ensp = null;

    // Path 1: ENST → direct lookup for Translation.id
    const tId = String(transcriptOrGene || '').replace(/\.\d+$/, '');
    if (tId.startsWith('ENST')) {
      const txRes = await fetchWithTimeout(
        `https://rest.ensembl.org/lookup/id/${tId}?content-type=application/json`
      );
      if (txRes.ok) {
        const txJson = await txRes.json();
        ensp = txJson.Translation?.id || null;
      }
    }

    // Path 2: gene symbol → canonical transcript → Translation.id
    if (!ensp && gene && gene !== '-') {
      const gRes = await fetchWithTimeout(
        `https://rest.ensembl.org/lookup/symbol/homo_sapiens/${encodeURIComponent(gene)}?expand=1&content-type=application/json`
      );
      if (gRes.ok) {
        const gJson = await gRes.json();
        const txList = gJson.Transcript || [];
        const canonical = selectByTranscript(txList, {
          gene,
          transcript:  transcriptOrGene,
          getGene:     () => gene,
          getTxId:     t => t.id,
          isCanonical: t => t.is_canonical === 1,
        });
        ensp = canonical?.Translation?.id || null;
        if (ensp) console.log(`[EnsemblDomains] ${gene} → canonical → ${ensp}`);
      }
    }

    if (!ensp) {
      console.warn(`[EnsemblDomains] Could not resolve ENSP for ${transcriptOrGene || gene}`);
      return [];
    }

    // Protein feature overlap — fetch all databases, keep Pfam only
    const featRes = await fetchWithTimeout(
      `https://rest.ensembl.org/overlap/translation/${ensp}?content-type=application/json`
    );
    if (!featRes.ok) return [];
    const features = await featRes.json();

    const pfam = Array.isArray(features) ? features.filter(f => f.type === 'Pfam') : [];
    console.log(`[EnsemblDomains] ${ensp}: ${pfam.length} Pfam domains`);
    return pfam.map(f => ({
      type:   'Domain',   // normalise so codon-viewer STRIP_TYPES picks it up
      tier:   'high',     // Pfam = major structural domain → blue in viewer
      start:  f.start,
      end:    f.end,
      name:   f.description || f.id,
      pfamId: f.id
    }));
  } catch (e) {
    console.warn('[EnsemblDomains] Fetch failed:', e.message);
    return [];
  }
}

/**
 * fetchCodonData(gene)
 * Inputs:  gene symbol string
 * Outputs: parsed data from data/codon/{GENE}.json (HTTP) or data/codon/{GENE}.js (file://), or null
 * Failures: silent — returns null if the file doesn't exist yet
 *
 * Strategy:
 *   1. Return from window._codonData cache if already loaded via <script> tag.
 *   2. Try fetch() on the .json file — works when served via HTTP.
 *   3. Fall back to a dynamic <script> tag loading the companion .js wrapper file.
 *      This bypasses Chrome's file:// CORS restriction that blocks fetch() on local files.
 *      The .js files are pre-generated wrappers: (window._codonData=…)["GENE"]={…};
 */
async function fetchCodonData(gene) {
  if (!gene || gene === '-') return null;

  // Cache hit — already loaded by a previous script-tag injection or earlier fetch
  if (window._codonData?.[gene]) return window._codonData[gene];

  // Fast path: fetch the JSON directly (works under HTTP server)
  try {
    const res = await originalFetch(`data/codon/${encodeURIComponent(gene)}.json`);
    if (res.ok) {
      const json = await res.json();
      (window._codonData = window._codonData || {})[gene] = json;
      return json;
    }
  } catch { /* fall through to script-tag strategy */ }

  // Fallback: dynamic <script> tag — works under file:// where fetch() is blocked by CORS.
  // The .js files set window["__codon_GENE"] = {...}; cache into window._codonData for reuse.
  return new Promise((resolve) => {
    const s = document.createElement('script');
    s.src = `data/codon/${encodeURIComponent(gene)}.js`;
    s.onload = () => {
      const d = window[`__codon_${gene}`] || window._codonData?.[gene] || null;
      if (d) (window._codonData = window._codonData || {})[gene] = d;
      resolve(d);
    };
    s.onerror = () => { console.warn('[CodonData] Not found for gene:', gene); resolve(null); };
    document.head.appendChild(s);
  });
}

/**
 * fetchExonCodingPositions(transcriptOrGene, gene)
 * Inputs:  transcript ID (ENST or NM_) and/or gene symbol
 * Outputs: [{ number, codonStart, codonEnd, cdsStart, cdsEnd, gc1, gc2, strand, chrom }, ...]
 *          ordered by transcript exon rank.
 *          cdsStart/cdsEnd are the c. (CDS) coordinates of the exon's first/last coding base.
 *          gc1/gc2 are the GRCh38 genomic coords of those same first/last coding bases
 *          (gc1↔cdsStart, gc2↔cdsEnd; gc1>gc2 on the minus strand) — used to map a
 *          codon range → genomic window for ClinVar chrpos38 search in the codon viewer.
 * Failures: silent — returns [] on any error
 *
 * Strategy:
 *   1. If transcriptOrGene is ENST* → direct /lookup/id/{ENST}?expand=1
 *   2. Otherwise (NM_ or missing) → /lookup/symbol/homo_sapiens/{gene}?expand=1
 *      then pick the canonical transcript from the gene object.
 * Maps each exon's CDS contribution to codon ranges. Handles +/− strand.
 * UTR-only exons are skipped.
 */
async function fetchExonCodingPositions(transcriptOrGene, gene) {
  if (!transcriptOrGene && !gene) return [];

  let tx = null;  // will hold the transcript object with Exon + Translation

  // ── Path 1: ENST ID → direct transcript lookup ──────────────────────────
  let tId = (transcriptOrGene || '').replace(/\.\d+$/, '');

  // NM_ → resolve to ENST first via a lightweight xref call (avoids the slow gene-expand
  // in Path 2, which fetches all transcripts for the gene and times out under load).
  // NOTE: /xrefs/id/{NM_} returns 400 "not found" — the correct endpoint is
  //       /xrefs/symbol/homo_sapiens/{NM_} which treats the accession as a symbol alias.
  if (tId.startsWith('NM_') || tId.startsWith('NR_') || tId.startsWith('XM_')) {
    try {
      const xrefRes = await fetchWithTimeout(
        `https://rest.ensembl.org/xrefs/symbol/homo_sapiens/${encodeURIComponent(tId)}?content-type=application/json`, {}, 8000
      );
      if (xrefRes.ok) {
        const xrefs = await xrefRes.json();
        const enst = (Array.isArray(xrefs) ? xrefs : []).find(x => x.id?.startsWith('ENST'));
        if (enst?.id) tId = enst.id;
      }
    } catch { /* fall through — will use gene-based lookup in Path 2 */ }
  }

  if (tId.startsWith('ENST')) {
    try {
      const res = await fetchWithTimeout(
        `https://rest.ensembl.org/lookup/id/${tId}?expand=1&content-type=application/json`
      );
      if (res.ok) tx = await res.json();
    } catch { /* fall through to gene-based lookup */ }
  }

  // ── Path 2: Gene symbol → canonical transcript from gene object ─────────
  // This is the heaviest fallback (full gene expand, all transcripts + exons). Use a
  // 20 s timeout — the default 10 s is too tight for slow Ensembl responses under load.
  if (!tx && gene && gene !== '-') {
    try {
      const res = await fetchWithTimeout(
        `https://rest.ensembl.org/lookup/symbol/homo_sapiens/${encodeURIComponent(gene)}?expand=1&content-type=application/json`,
        {}, 20000
      );
      if (!res.ok) return [];
      const geneObj = await res.json();
      const txList  = geneObj.Transcript || [];
      tx = selectByTranscript(txList, {
        gene,
        transcript:  transcriptOrGene,
        getGene:     () => gene,
        getTxId:     t => t.id,
        isCanonical: t => t.is_canonical === 1,
      });
      if (tx) console.log(`[ExonCodingPositions] Gene ${gene} → canonical ${tx.id} (${tx.Exon?.length} exons)`);
    } catch { return []; }
  }

  if (!tx) return [];

  const exons = tx.Exon || [];
  const tr    = tx.Translation;
  if (!tr || !exons.length) return [];

  const cdsHigh = Math.max(tr.start, tr.end);
  const cdsLow  = Math.min(tr.start, tr.end);
  const strand  = tx.strand;
  const chrom   = tx.seq_region_name || null;

  const sorted = [...exons].sort((a, b) =>
    strand > 0 ? a.start - b.start : b.start - a.start
  );

  let cdsAccum = 0;
  const result = [];
  let exonNum  = 0;
  for (const exon of sorted) {
    exonNum++;
    const overlapHigh = Math.min(exon.end, cdsHigh);
    const overlapLow  = Math.max(exon.start, cdsLow);
    const cdsLen      = Math.max(0, overlapHigh - overlapLow + 1);
    if (cdsLen === 0) continue;
    const cdsStart   = cdsAccum + 1;          // c. coord of this exon's first coding base
    const codonStart = Math.floor(cdsAccum / 3) + 1;
    cdsAccum += cdsLen;
    const cdsEnd     = cdsAccum;              // c. coord of this exon's last coding base
    const codonEnd   = Math.ceil(cdsAccum / 3);
    // GRCh38 genomic coords of this exon's first/last CODING base, so the codon
    // viewer can map a codon range → genomic window for ClinVar chrpos38 search.
    const gc1 = strand < 0 ? overlapHigh : overlapLow;  // genomic of cdsStart
    const gc2 = strand < 0 ? overlapLow  : overlapHigh; // genomic of cdsEnd
    result.push({ number: exonNum, codonStart, codonEnd, cdsStart, cdsEnd, gc1, gc2, strand, chrom });
  }
  return result;
}

/**
 * Zero-latency lookup against the cached SVIS index.
 */
async function checkVCEP(geneSymbol) {
  if (!geneSymbol || geneSymbol === '-') {
    data.vcepGuideline = null;
    renderVCEPBanner();
    return;
  }

  // Failsafe: if cache isn't ready yet, await it
  if (cspecSVISCache === null) {
    await prefetchClinGenSVIS();
  }

  const targetGene = geneSymbol.trim().toUpperCase();

  try {
    const foundGuideline = cspecSVISCache.find(item => {
      // Defensive parsing: schemas vary slightly between CSpec API versions
      const affiliationName = String(item.affiliation?.label || item.affiliation?.name || '').toUpperCase();
      const status = String(item.status || item.sviStatus || item.docStatus || '').toLowerCase();

      let matchesGene = false;

      // Match exact gene in affiliation name (like "TP53 VCEP")
      if (affiliationName === targetGene || affiliationName.includes(`${targetGene} `) || affiliationName.includes(`${targetGene}-`)) {
        matchesGene = true;
      }

      // Match exact gene in ruleSets array
      if (!matchesGene && Array.isArray(item.ruleSets)) {
        for (const rs of item.ruleSets) {
          if (Array.isArray(rs.genes)) {
            for (const g of rs.genes) {
              const gLabel = String(g.label || g.name || '').toUpperCase();
              if (gLabel === targetGene) {
                matchesGene = true;
                break;
              }
            }
          }
          if (matchesGene) break;
        }
      }

      // Ensure we only flag published guidelines, ignoring superseded drafts
      const isReleased = status.includes('released') || status.includes('published');

      return matchesGene && isReleased;
    });

    if (foundGuideline) {
      data.vcepGuideline = {
        gene: targetGene,
        version: foundGuideline.version || "Latest",
        condition: foundGuideline.condition || foundGuideline.disease?.name || foundGuideline.disease?.label || "Condition unspecified",
        // Pluck the correct ID for the UI hyperlink mapping
        cspecId: foundGuideline.id || foundGuideline.uuid || foundGuideline.sviId || foundGuideline['@id'] || foundGuideline.url
      };

      // Extract the final identifier segment if it is a URL
      if (data.vcepGuideline.cspecId && String(data.vcepGuideline.cspecId).includes('/')) {
        const parts = String(data.vcepGuideline.cspecId).split('/');
        data.vcepGuideline.cspecId = parts[parts.length - 1];
      }
    } else {
      data.vcepGuideline = null;
    }

    renderVCEPBanner();
  } catch (e) {
    console.error("CSpec Cache Lookup Error:", e);
    data.vcepGuideline = null;
    renderVCEPBanner();
  }
}

// GRCh38 chromosome → RefSeq accession (for VCF → genomic HGVS conversion)
const CHR_TO_NC38 = {
  '1':'NC_000001.11','2':'NC_000002.12','3':'NC_000003.12','4':'NC_000004.12',
  '5':'NC_000005.10','6':'NC_000006.12','7':'NC_000007.14','8':'NC_000008.11',
  '9':'NC_000009.12','10':'NC_000010.11','11':'NC_000011.10','12':'NC_000012.12',
  '13':'NC_000013.11','14':'NC_000014.9','15':'NC_000015.10','16':'NC_000016.10',
  '17':'NC_000017.11','18':'NC_000018.10','19':'NC_000019.10','20':'NC_000020.11',
  '21':'NC_000021.9','22':'NC_000022.11','X':'NC_000023.11','Y':'NC_000024.10','MT':'NC_012920.1'
};

async function fetchClinGenFallback(hgvs, gen) {
  setDot('vv', 'loading');
  setAPIStatus('clingen', 'loading');
  document.getElementById('statusMsg').innerHTML =
    '⚠️ <b style="color:var(--amber)">VV & Ensembl unavailable. Trying ClinGen pipeline…</b>';

  let clean = hgvs.replace(/\([^)]*\)/g, '').replace(/\s*\(?p\..*/i, '').trim();

  // Convert VCF format (1-103017881-G-C or chr1-103017881-G-C) to genomic HGVS
  if (!clean.includes(':')) {
    const parts = clean.replace(/^chr/i, '').split(/[-\s]+/);
    if (parts.length >= 4) {
      const [chr, pos, ref, alt] = parts;
      const acc = CHR_TO_NC38[chr.toUpperCase()] || CHR_TO_NC38[chr];
      if (acc) {
        // SNV: NC_000001.11:g.103017881G>C
        // Deletion: NC_000001.11:g.103017881del  (ref longer than alt)
        // Insertion: NC_000001.11:g.103017881_103017882insATG
        if (ref.length === 1 && alt.length === 1) {
          clean = `${acc}:g.${pos}${ref}>${alt}`;
        } else if (alt === '-' || alt.length < ref.length) {
          const delStart = parseInt(pos) + 1;
          const delEnd = parseInt(pos) + ref.length - 1;
          clean = delStart === delEnd
            ? `${acc}:g.${delStart}del`
            : `${acc}:g.${delStart}_${delEnd}del`;
        } else if (ref === '-' || ref.length < alt.length) {
          const ins = alt.substring(ref.length);
          clean = `${acc}:g.${pos}_${parseInt(pos)+1}ins${ins}`;
        } else {
          clean = `${acc}:g.${pos}${ref}>${alt}`;
        }
      }
    }
  }

  // Step A: ClinGen Allele Registry → CAid + coords + transcript + gene
  const arRes = await fetchWithTimeout(
    `https://reg.genome.network/allele?hgvs=${encodeURIComponent(clean)}`
  );
  if (!arRes.ok) throw new Error(`ClinGen AR ${arRes.status}`);
  const ar = await arRes.json();

  const caid = ar.caid || ar['@id']?.split('/').pop() || null;
  if (!caid) throw new Error('ClinGen AR: no CAid returned');

  // GRCh38 coordinates (required)
  const g38 = ar.genomicAlleles?.find(g => g.referenceGenome === 'GRCh38');
  const c38 = g38?.coordinates?.[0];
  if (!c38 || !g38?.chromosome) throw new Error('ClinGen AR: no GRCh38 coords');

  data.coords.chrom = String(g38.chromosome).replace('chr', '');
  const rawRef38 = c38.referenceAllele || '';
  const rawAlt38 = c38.allele || '';

  if (rawRef38 === '-' || rawAlt38 === '-' || !rawRef38 || !rawAlt38) {
    // [FIX] HGVS-style indel from ClinGen AR (ref or alt is '-'/empty).
    // Anchor position (1-based) = c38.start — see hgvsIndelToVcf for derivation.
    const vcf38 = await hgvsIndelToVcf(
      data.coords.chrom, c38.start,
      rawRef38 === '-' ? '' : rawRef38,
      rawAlt38 === '-' ? '' : rawAlt38
    );
    if (vcf38) {
      data.coords.pos38      = vcf38.pos;
      data.coords.ref        = vcf38.ref;
      data.coords.alt        = vcf38.alt;
      data.coords.hg38String = `chr${data.coords.chrom}-${vcf38.pos}-${vcf38.ref}-${vcf38.alt}`;
      document.getElementById('dHg38').innerText = data.coords.hg38String;
    } else {
      // Both Ensembl and NCBI anchor fetches failed — coords are incomplete.
      data.coords.pos38 = String(c38.start + 1);
      data.coords.ref   = normalizeVcfAllele(rawRef38);
      data.coords.alt   = normalizeVcfAllele(rawAlt38);
      document.getElementById('statusMsg').innerHTML =
        '⚠️ <b style="color:var(--amber)">Indel anchor lookup failed (Ensembl + NCBI). Coordinates may be incomplete — downstream results unreliable.</b>';
    }
  } else {
    data.coords.pos38 = String(c38.start + 1);
    data.coords.ref   = normalizeVcfAllele(rawRef38);
    data.coords.alt   = normalizeVcfAllele(rawAlt38);
    if (data.coords.ref && data.coords.alt) {
      data.coords.hg38String = `chr${data.coords.chrom}-${data.coords.pos38}-${data.coords.ref}-${data.coords.alt}`;
      document.getElementById('dHg38').innerText = data.coords.hg38String;
    }
  }

  // GRCh37 / hg19 coordinates (best-effort)
  const g37 = ar.genomicAlleles?.find(g => g.referenceGenome === 'GRCh37');
  const c37 = g37?.coordinates?.[0];
  if (c37 && g37?.chromosome) {
    const chrom37  = String(g37.chromosome).replace('chr', '');
    const rawRef37 = c37.referenceAllele || '';
    const rawAlt37 = c37.allele || '';
    if (rawRef37 === '-' || rawAlt37 === '-' || !rawRef37 || !rawAlt37) {
      // [FIX] HGVS-style hg19 indel — lift over the VCF-style hg38 coords via Ensembl map
      if (data.coords.hg38String && data.coords.ref && data.coords.alt) {
        await fetchEnsemblMap(data.coords.chrom, data.coords.pos38, data.coords.ref, data.coords.alt);
      }
    } else {
      data.coords.pos19 = String(c37.start + 1);
      const ref37 = normalizeVcfAllele(rawRef37);
      const alt37 = normalizeVcfAllele(rawAlt37);
      if (ref37 && alt37) {
        data.coords.hg19String = `chr${chrom37}-${data.coords.pos19}-${ref37}-${alt37}`;
        document.getElementById('dHg19').innerText = data.coords.hg19String;
      }
    }
  }

  // Transcript — prefer MANE Select; ClinGen AR may return gene under geneSymbol or gene.label
  const txAlleles = ar.transcriptAlleles || [];
  const getGeneFromTx = t => t.geneSymbol || t.gene?.label || t.gene?.symbol || null;
  const isManeSelect = t => {
    const m = t.MANE ?? t.mane ?? '';
    return typeof m === 'string'
      ? m.toUpperCase().includes('SELECT')
      : !!m; // boolean true = MANE Select
  };
  const maneTx = selectByTranscript(txAlleles, {
    gene:        data.coords.gene,
    transcript:  data.ensembl.transcript,
    getGene:     getGeneFromTx,
    getTxId:     t => t.hgvsMatchStrings?.[0]?.split(':')[0],
    isCanonical: isManeSelect,
  });
  if (maneTx) {
    data.coords.gene = getGeneFromTx(maneTx) || data.coords.gene;
    const hgvsStr = maneTx.hgvsMatchStrings?.[0] || '';
    if (hgvsStr.includes(':')) {
      const colonIdx = hgvsStr.indexOf(':');
      const tx    = hgvsStr.slice(0, colonIdx);
      const cHgvs = hgvsStr.slice(colonIdx + 1);
      data.ensembl.transcript = tx;
      data.coords.hgvs = cHgvs;
      document.getElementById('eManeTranscript').innerText = tx;
      document.getElementById('dHgvs').innerText = cHgvs;
    }
    // Protein notation
    const pHgvs = maneTx.proteinAlleles?.[0]?.hgvsMatchStrings?.[0] || '';
    if (pHgvs) {
      data.ensembl.protein = pHgvs.includes(':') ? pHgvs.split(':')[1] : pHgvs;
      document.getElementById('dProtein').innerText = data.ensembl.protein;
      // Extract protein position for domain annotation
      const pMatch = data.ensembl.protein.match(/p\.\(?([A-Z][a-z]{2})(\d+)/);
      if (pMatch) data.ensembl.vepProteinStart = parseInt(pMatch[2], 10);
    }
    // VEP consequence type from HGVS if not yet set
    if (!data.ensembl.vepConsequence && maneTx.hgvsMatchStrings?.[0]) {
      const h = data.coords.hgvs || '';
      if (/del/.test(h)) data.ensembl.vepConsequence = 'deletion';
      else if (/ins/.test(h) || /dup/.test(h)) data.ensembl.vepConsequence = 'insertion';
      else if (/>[A-Z]$/.test(h)) data.ensembl.vepConsequence = 'missense_variant';
    }
  }

  // Gene fallback to geneInput if AR doesn't provide
  data.coords.gene = data.coords.gene || document.getElementById('geneInput').value.trim().toUpperCase();
  if (data.coords.gene) document.getElementById('dGene').innerText = data.coords.gene;

  // CAid
  data.caId = caid;

  // Step B: LDH → REVEL + in silico scores (best-effort, non-blocking)
  try {
    const ldhRes = await fetchWithTimeout(
      `https://ldh.clinicalgenome.org/ldh/Variant/id/${caid}?detail=high`
    );
    if (ldhRes.ok) {
      const ldh = await ldhRes.json();
      const ld = ldh.data?.ld || {};

      // REVEL (MANE Select only)
      const revelEntries = ld.RevelScore?.[0]?.entContent || [];
      const maneRevel = revelEntries.find(e => e.mane === 'MANE Select');
      if (maneRevel?.score != null) {
        data.scores.revel = parseFloat(maneRevel.score);
        document.getElementById('mvRevel').innerText = data.scores.revel.toFixed(3);
      }
    }
  } catch (e) { console.warn('[ClinGen LDH fallback] in silico fetch failed:', e.message); }

  // Enable UI
  enableWrappers(['franklin', 'gnomadv4', 'spliceai', 'liftover', 'clingen', 'clinvar',
                  'omim', 'gr', 'scholar', 'decipher', 'hgmd', 'gtex', 'mastermind', 'gnomadv2']);
  document.getElementById('btnLaunchAll').disabled = false;
  document.getElementById('btnLaunchSel').disabled = false;
  if (data.coords.gene) document.getElementById('geneInput').value = data.coords.gene;
  setDot('vv', 'ok');
  setAPIStatus('clingen', 'ok');
  document.getElementById('statusMsg').innerHTML =
    `✅ <b style="color:var(--teal)">ClinGen pipeline — basic annotation from Allele Registry</b>`;
}

async function resolveHgvsToHg38(hgvs) {
  const clean = hgvs.replace(/\([^)]*\)/g, '').replace(/\s*\(?p\..*/i, '').trim();
  const url = `https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/${encodeURIComponent(clean)}/mane_select`;
  const res = await fetchWithTimeout(url, {}, 30000);
  if (!res.ok) throw new Error(`VV ${res.status}`);
  const json = await res.json();
  for (const key of Object.keys(json)) {
    const rec = json[key];
    const g38 = rec?.primary_assembly_loci?.grch38 || rec?.primary_assembly_loci?.hg38;
    if (g38?.vcf) {
      const v = g38.vcf;
      return { chr: v.chr.replace('chr', ''), pos: parseInt(v.pos, 10), ref: v.ref, alt: v.alt };
    }
  }
  throw new Error('No GRCh38 coords in VV response');
}

async function fetchClinVarCodon(gen) {
  if (!data.coords.gene) {
    document.getElementById('cvAltCodon').innerText = 'N/A';
    document.getElementById('cvCodonText').innerText = 'Codon';
    return;
  }

  let isSplice = false;
  let spliceType = '';
  let anchorBase = 0;

  let cons = data.ensembl.vepConsequence || '';

  // Fallback: detect start-loss from HGVS notation when VEP didn't run.
  // Splice handling moved to evaluateSplicePS1 (hybrid split — see CLAUDE.md).
  if (!cons && data.coords.hgvs) {
    if (/^c\.[123][^0-9]/.test(data.coords.hgvs)) {
      cons = 'start_lost';
    }
  }
  // Additional start-loss check from p. notation (p.Met1, p.M1)
  if (!cons.includes('start_lost') && data.ensembl.protein) {
    if (/p\.\(?(Met1[^0-9]|M1[^0-9])/.test(data.ensembl.protein)) {
      cons = cons ? `${cons},start_lost` : 'start_lost';
    }
  }

  // Hybrid: splice variants are owned by evaluateSplicePS1; skip here entirely.
  // This keeps fetchClinVarCodon focused on missense + start-loss codon lookups.
  const looksLikeSplice = (typeof cons === 'string' && cons.includes('splice'))
                       || (data.coords.hgvs && /c\.\d+[+-]\d+/.test(data.coords.hgvs));
  if (looksLikeSplice) {
    const altBlock = document.getElementById('cvAltBlock');
    const labelTextEl = document.getElementById('cvAltLabelText');
    const codonEl = document.getElementById('cvAltCodon');
    if (altBlock) altBlock.style.display = 'none';
    if (labelTextEl) labelTextEl.innerHTML = '';
    if (codonEl) codonEl.innerHTML = '';
    return;
  }

  if (cons.includes('splice')) {
    isSplice = true;
    if (cons.includes('splice_donor')) spliceType = 'donor';
    else if (cons.includes('splice_acceptor')) spliceType = 'acceptor';
  }

  if (isSplice && data.coords.hgvs) {
    const cMatch = data.coords.hgvs.match(/c\.(\d+)([+-]\d+)?/);
    if (cMatch) {
      anchorBase = parseInt(cMatch[1], 10);
      const offset = cMatch[2] || '';
      if (!spliceType) {
        if (offset.startsWith('+')) spliceType = 'donor';
        else if (offset.startsWith('-')) spliceType = 'acceptor';
        else spliceType = 'donor'; // fallback guess
      }
    } else {
      isSplice = false;
    }
  } else {
    isSplice = false;
  }

  let term = '';
  let searchUrl = '';
  let uiText = '';

  // Variables for Missense fallback
  let pos = '';
  let ref3 = '';
  const hgvsP = (data.ensembl.protein && data.ensembl.protein !== '-') ? data.ensembl.protein : '';

  if (isSplice && anchorBase > 0) {
    const terms = [];
    if (spliceType === 'donor') {
      // Donor motif: last 3 of exon, 3-6 intronic
      for (let i = -2; i <= 0; i++) if (anchorBase + i > 0) terms.push(`"c.${anchorBase + i}"`);
      for (let i = 1; i <= 6; i++) terms.push(`"c.${anchorBase}+${i}"`);
      uiText = 'Alt Splice (Donor)';
    } else {
      // Acceptor motif: intronic -20 to -1, plus first base of exon
      for (let i = -20; i <= -1; i++) terms.push(`"c.${anchorBase}${i}"`);
      terms.push(`"c.${anchorBase}"`);
      uiText = 'Alt Splice (Acceptor)';
    }

    // Attempt ultra-robust genomic range search
    try {
      const b1 = spliceType === 'donor' ? `c.${anchorBase - 2}del` : `c.${anchorBase}-20del`;
      const b2 = spliceType === 'donor' ? `c.${anchorBase}+6del` : `c.${anchorBase}del`;
      const tx = data.ensembl.transcript || '';
      const [l1, l2] = await Promise.all([
        resolveHgvsToHg38(`${tx}:${b1}`).catch(() => null),
        resolveHgvsToHg38(`${tx}:${b2}`).catch(() => null)
      ]);
      if (l1 && l2) {
        const gMin = Math.min(l1.pos, l2.pos);
        const gMax = Math.max(l1.pos, l2.pos);
        // Query by genomic range to capture all variants in the motif
        term = `${data.coords.gene}[gene] AND GRCh38:${data.coords.chrom}:${gMin}-${gMax} AND ("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])`;
      }
    } catch (e) { console.warn('[AltSplice] Range resolution failed:', e); }

    if (!term) {
      const txPrefix = data.ensembl.transcript ? data.ensembl.transcript.split('.')[0] : '';
      term = `${data.coords.gene}[gene] AND (${txPrefix}) AND (${terms.join(' OR ')}) AND ("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])`;
    }
    searchUrl = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(term)}`;

  } else {
    // Standard Codon Logic
    if (!hgvsP && !cons.includes('start_lost')) {
      document.getElementById('cvAltCodon').innerText = 'N/A';
      document.getElementById('cvCodonText').innerText = 'Codon';
      return;
    }
    const pMatch = hgvsP.replace(/[()]/g, '').match(/p\.([A-Z][a-z]{2}|[A-Z])(\d+)/);
    pos = pMatch ? pMatch[2] : '';
    ref3 = pMatch ? pMatch[1] : '';

    // Handle start-loss variants where p. nomenclature might be missing or ambiguous (e.g. p.?)
    if (!pos && cons.includes('start_lost')) {
      pos = '1';
      ref3 = 'Met';
    }

    if (!pos) {
      document.getElementById('cvAltCodon').innerText = 'N/A';
      document.getElementById('cvCodonText').innerText = 'Codon';
      return;
    }
    uiText = `Alternative Missense Variant(s) (p.${ref3}${pos})`;
    // Use ClinVar [VARNAME] + [GENE] field tags (works for both esearch API and browser)
    if (cons.includes('start_lost')) {
      // Start-loss: search for Met1 / M1 variants
      const variantParam = `${ref3}${pos}`;          // e.g. Met1
      const variantShort = `${ref3.charAt(0)}${pos}`; // e.g. M1
      const geneParam = data.coords.gene.toLowerCase();
      term = `("${variantParam}"[VARNAME] OR "${variantShort}"[VARNAME]) AND "${geneParam}"[GENE]`;
      searchUrl = `https://www.ncbi.nlm.nih.gov/clinvar/?variant=${encodeURIComponent(variantParam)}&gene=${encodeURIComponent(geneParam)}&term=${encodeURIComponent(term)}`;
    } else {
      // Standard missense codon search: "Arg1695"[VARNAME] AND "col11a2"[GENE]
      const variantParam = `${ref3}${pos}`;
      const geneParam = data.coords.gene.toLowerCase();
      term = `"${variantParam}"[VARNAME] AND "${geneParam}"[GENE]`;
      searchUrl = `https://www.ncbi.nlm.nih.gov/clinvar/?variant=${encodeURIComponent(variantParam)}&gene=${encodeURIComponent(geneParam)}&term=${encodeURIComponent(term)}`;
    }
  }

  // 4. Update UI Headers
  const altBlock = document.getElementById('cvAltBlock');
  const labelTextEl = document.getElementById('cvAltLabelText');
  const codonWrapper = document.getElementById('cvCodonWrapper');

  if (labelTextEl) {
    const spliceSuffix = spliceType === 'donor' ? '(DONOR)' : '(ACCEPTOR)';
    const fullTitle = isSplice ? `🧬 Alternative Splice Variant(s) ${spliceSuffix}` : `🧬 Alternative Missense Variant(s)`;
    labelTextEl.innerHTML = `<div style="width:100%; font-size:0.85rem; font-weight:800; color:var(--text-bright); border-bottom:1px solid rgba(255,255,255,0.1); padding-bottom:8px; margin-bottom:12px; letter-spacing:0.5px; text-transform:uppercase;">${fullTitle}</div>`;
  }
  if (altBlock) altBlock.style.display = 'block';
  if (codonWrapper) codonWrapper.style.display = 'none';

  const url = `${NCBI_BASE}/esearch.fcgi?db=clinvar&term=${encodeURIComponent(term)}&retmode=json&api_key=${NCBI_KEY}&email=${NCBI_EMAIL}&retmax=500`;

  try {
    const res = await fetchWithRetry(url);
    if (!res.ok) throw new Error('CV Alt search fail ' + res.status);

    const json = await res.json();
  if (gen !== window.currentSearchGen()) return;
    const uids = json.esearchresult?.idlist || [];
    if (!uids.length) {
      document.getElementById('cvAltCodon').innerHTML = `<a href="${searchUrl}" target="_blank" style="color:var(--dim);text-decoration:none;font-size:0.75rem;">None found. Search ↗</a>`;
      return;
    }

    // Limit UIDs to top 100 to prevent URL length issues
    const topUids = uids.slice(0, 100);
    await new Promise(r => setTimeout(r, 400));
    const sumUrl = `${NCBI_BASE}/esummary.fcgi?db=clinvar&id=${topUids.join(',')}&retmode=json&api_key=${NCBI_KEY}&email=${NCBI_EMAIL}`;
    let sumJson = null;

    for (let attempt = 0; attempt < 3; attempt++) {
      const sumRes = await fetchWithRetry(sumUrl);
      if (sumRes.ok) {
        sumJson = await ncbiSafeJson(sumRes);
        if (sumJson.result?.uids) break;
      }
      await new Promise(r => setTimeout(r, 600 * (attempt + 1)));
    }

    if (!sumJson || !sumJson.result?.uids) throw new Error('CV Alt summary failed');

    const altMap = {};
    const ps1Candidates = []; // P/LP variants with same AA change (PS1 evidence)
    (sumJson.result?.uids || []).forEach(uid => {
      const rec = sumJson.result[uid]; if (!rec) return;

      let badgeText = '';
      if (isSplice) {
        if ((rec.title || '').includes(data.coords.hgvs)) return;
        // Strict filter: Only accept titles that look like splicing variants (c. numbering with +/-)
        if (!rec.title.match(/c\.\d+[+-]/)) return;

        const fullMatch = (rec.title || '').match(/c\.[^:\s()]+/);
        if (!fullMatch) return;
        badgeText = fullMatch[0];
      } else {
        // Missense mode: Support both 1-letter and 3-letter matching for robustness
        const missenseRegex = new RegExp("p\\.([A-Z][a-z]{2}|[A-Z])" + pos + "([A-Z][a-z]{2}|[A-Z]|Ter|\\*|\\?)", "i");
        const match = (rec.title || '').match(missenseRegex);
        if (!match) return;

        const altAA = match[2];
        if (hgvsP && hgvsP.toLowerCase().includes(altAA.toLowerCase())) {
          // Same amino acid change — collect as PS1 candidate instead of discarding
          const ps1Sig = rec.germline_classification?.description || rec.clinical_significance?.description || '';
          const ps1SigL = ps1Sig.toLowerCase();
          if (ps1SigL.includes('pathogenic') && !ps1SigL.includes('conflicting')) {
            const ps1Status = rec.germline_classification?.review_status || rec.clinical_significance?.review_status || '';
            ps1Candidates.push({ sig: ps1Sig, stars: reviewStars(ps1Status) });
          }
          return;
        }
        badgeText = (rec.title.match(/p\.[^:\s()]+/)?.[0]) || `p.${ref3}${pos}${altAA}`;
      }

      let sig = rec.germline_classification?.description || rec.clinical_significance?.description || 'Unknown';
      let status = rec.germline_classification?.review_status || rec.clinical_significance?.review_status || '';

      const getStarsHtml = (s) => {
        const n = reviewStars(s);
        const starStr = '★'.repeat(n) + '☆'.repeat(4 - n);
        return `<span style="color:var(--amber); font-size:0.75rem; margin-left:4px; vertical-align:middle;">${starStr}</span> <span style="font-size:0.62rem; opacity:0.8;">(${n}/4)</span>`;
      };
      const starsHtml = getStarsHtml(status);

      if (!altMap[badgeText]) altMap[badgeText] = { sigs: new Set(), hgvsT: null, stars: starsHtml, rsid: null, altAllele: null };
      altMap[badgeText].sigs.add(sig);

      const titleStr = rec.title || '';
      if (!altMap[badgeText].hgvsT) {
        let m = titleStr.match(/(NM_\d+\.\d+)\([^)]+\):(c\.[^\s(,)]+)/);
        if (!m) m = titleStr.match(/(NM_\d+\.\d+):(c\.[^\s(,)]+)/);
        if (!m) m = titleStr.match(/(NM_\d+):(c\.[^\s(,)]+)/);
        if (m) altMap[badgeText].hgvsT = `${m[1]}:${m[2]}`;
      }

      if (altMap[badgeText].hgvsT && !altMap[badgeText].altAllele) {
        const am = altMap[badgeText].hgvsT.match(/[A-Z]>([A-Z])/);
        if (am) altMap[badgeText].altAllele = am[1];
      }

      if (!altMap[badgeText].rsid) {
        const xrefs = rec.variation_set?.[0]?.variation_xrefs || [];
        const dbsnp = xrefs.find(x => x.db_source === 'dbSNP');
        if (dbsnp) altMap[badgeText].rsid = dbsnp.db_id;
      }
    });

    const badgePromises = [];

    // Batch Resolve HGVS via VariantValidator (more accurate for indels than Ensembl)
    const batchResolveHgvs = async (variants) => {
      if (!variants.length) return {};
      const results = {};
      try {
        const url = `https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/${encodeURIComponent(variants.join('|'))}/mane_select`;
        const res = await fetchWithTimeout(url, {}, 30000);
        if (res.ok) {
          const json = await res.json();
          for (const hgvs of variants) {
            const recData = json[hgvs]; // renamed — 'data' would shadow the global variant state
            if (recData && recData.primary_assembly_loci) {
              const v = recData.primary_assembly_loci.grch38?.vcf || recData.primary_assembly_loci.hg38?.vcf;
              if (v) results[hgvs] = { chr: v.chr.replace('chr', ''), pos: v.pos, ref: v.ref, alt: v.alt };
            }
          }
        }
      } catch (e) { console.warn('[AltBatch] VV Batch failed:', e); }
      return results;
    };

    // Helper: Map chromosome to GRCh38 NC accession
    const getNCAcc = (chr) => {
      const c = chr.replace('chr', '');
      const map = {
        '1':'000001.11','2':'000002.12','3':'000003.12','4':'000004.12','5':'000005.10',
        '6':'000006.12','7':'000007.14','8':'000008.11','9':'000009.12','10':'000010.11',
        '11':'000011.10','12':'000012.12','13':'000013.11','14':'000014.9','15':'000015.10',
        '16':'000016.10','17':'000017.11','18':'000018.10','19':'000019.10','20':'000020.11',
        '21':'000021.9','22':'000022.11','X':'000023.11','Y':'000024.11','M':'012920.1'
      };
      return map[c] ? `NC_${map[c]}` : null;
    };

    // Helper: resolve HGVS to ClinGen CAids in batch (Multi-fallback strategy)
    const resolveBatchCaids = async (variantList) => {
      if (!variantList.length) return {};
      const results = {};
      console.log('[ClinGenBatch] Resolving variants:', variantList);
      try {
        const promises = variantList.map(async (v) => {
          if (!v) return;
          const queries = typeof v === 'string' ? [v] : [v.hgvsT, v.gStr, v.pStr].filter(Boolean);
          
          for (let q of queries) {
             // Deep clean: remove gene symbols, parentheses, spaces, AND emojis/checkmarks
             const cleanQ = q.replace(/\(.*\)/g, '')
                             .replace(/[()]/g, '')
                             .replace(/\s+/g, '')
                             .replace(/[^\x20-\x7E]/g, ''); // Strip all non-ASCII (emojis)
             
             const url = `https://reg.genome.network/allele?hgvs=${encodeURIComponent(cleanQ)}`;
             try {
               const res = await fetchWithTimeout(url);
               if (res.ok) {
                 const json = await res.json();
                 let caid = json.caid;
                 if (!caid && json["@id"]) caid = json["@id"].split('/').pop();
                 
                 if (caid) {
                   console.log(`[ClinGenBatch] Found CAid: ${q} -> ${caid}`);
                   results[typeof v === 'string' ? v : v.hgvsT] = caid;
                   break;
                 }
               }
             } catch(err) {}
          }
        });
        await Promise.all(promises);
      } catch (e) { console.warn('[ClinGenBatch] CAid resolution failed:', e); }
      return results;
    };

    // ── REVEL via MyVariant.info (coordinate-based, no CAid needed) ──────────
    // Takes { hgvsT -> { chr, pos, ref, alt } }, returns { hgvsT -> score }.
    // Fires one GET per variant in parallel — reliable for 1 or many variants.
    const fetchRevelByCoords = async (coordMap) => {
      const entries = Object.entries(coordMap)
        .filter(([, c]) => c && String(c.ref).length === 1 && String(c.alt).length === 1);
      if (!entries.length) return {};

      const results = {};
      const maneTx = (document.getElementById('eManeTranscript')?.innerText || '').trim();
      // Fields needed — no dotfield so getBestRevel can traverse the nested object
      const FIELDS = 'dbnsfp.revel_score,dbnsfp.refseq_id,revel.revel_score,revel.refseq_id';

      await Promise.all(entries.map(async ([hgvsT, coords]) => {
        // Build the genomic HGVS query MyVariant understands best
        const q = `chr${coords.chr}:g.${coords.pos}${coords.ref}>${coords.alt}`;
        try {
          const r = await fetchWithTimeout(
            `https://myvariant.info/v1/query?q=${encodeURIComponent(q)}&fields=${FIELDS}&size=5`
          );
          if (!r.ok) {
            console.warn(`[REVEL] HTTP ${r.status} for ${q}`);
            return;
          }
          const j = await r.json();
          const hits = j.hits || [];

          // Alt-allele filter first (prevents wrong-strand hits), then transcript-aware selection.
          const targetAlt = coords.alt.toUpperCase();
          const altPool  = hits.filter(h => (h._id || '').toUpperCase().endsWith('>' + targetAlt) && (h.dbnsfp || h.revel));
          const dataPool = altPool.length ? altPool : hits.filter(h => h.dbnsfp || h.revel);
          const best = selectByTranscript(dataPool.length ? dataPool : hits, {
            gene:        data.coords.gene,
            transcript:  maneTx,
            getGene:     h => (Array.isArray(h.dbnsfp?.genename) ? h.dbnsfp.genename[0] : h.dbnsfp?.genename),
            getTxId:     h => (Array.isArray(h.dbnsfp?.refseq_id) ? h.dbnsfp.refseq_id[0] : h.dbnsfp?.refseq_id),
            isCanonical: h => !!h.dbnsfp?.mane,
          });

          if (best) {
            const score = getBestRevel(best, maneTx);
            if (score !== null) {
              console.log(`[REVEL] ${q} → ${score}`);
              results[hgvsT] = score;
            } else {
              console.warn(`[REVEL] getBestRevel returned null for ${q}. Hit:`, JSON.stringify(best).slice(0, 200));
            }
          } else {
            console.warn(`[REVEL] No usable hit for ${q}`);
          }
        } catch (e) {
          console.warn(`[REVEL] fetch failed for ${q}:`, e.message);
        }
      }));

      return results;
    };

    // Helper: fetch Expert Pathogenicity for a single CAid
    const fetchSingleExpert = async (caid) => {
      if (!caid) return null;
      try {
        const url = `https://ldh.clinicalgenome.org/ldh/Variant/id/${caid}?detail=high`;
        const res = await fetchWithTimeout(url);
        if (res.ok) {
          const resJson = await res.json(); // renamed — 'data' would shadow global variant state
          const ld = resJson.ld || {};
          const pathInfo = ld.PathogenicityClassification?.[0]?.entContent || ld.ClinicalInterpretation?.[0]?.entContent;
          if (pathInfo) {
            return {
              classification: pathInfo.classification || pathInfo.significance,
              vcep: pathInfo.source?.label,
              date: pathInfo.date
            };
          }
        }
      } catch (e) {}
      return null;
    };


    // Helper: fetch SpliceAI scores from Broad API using hg38 VCF coords
    const fetchSpliceAIBroad = async (loc, gene, transcript) => {
      const varStr = `chr${loc.chr}-${loc.pos}-${loc.ref}-${loc.alt}`;
      const url = `https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&distance=500&mask=0&variant=${varStr}`;
      try {
        console.log(`[AltSplice] Broad SpliceAI query: ${varStr}`);
        const controller = new AbortController();
        const timeoutId = setTimeout(() => controller.abort(), 30000);
        const res = await originalFetch(url, { signal: controller.signal });
        clearTimeout(timeoutId);
        if (!res.ok) return null;
        const json = await res.json();
        if (json && json.scores && json.scores.length > 0) {
          const picked = pickSpliceAIScore(json.scores, gene, transcript);
          console.log(`[AltSplice] Broad SpliceAI scores for ${varStr}:`, picked);
          return picked;
        }
      } catch (e) { console.warn(`[AltSplice] Broad SpliceAI failed for ${varStr}:`, e); }
      return null;
    };

    // Helper: format full 2x2 SpliceAI table (Broad API uppercase keys)
    const formatSAITable = (scores) => {
      if (!scores) return '<div style="font-size:0.65rem; color:var(--dim); margin-top:2px;">SpliceAI: N/A</div>';

      const ag = parseFloat(scores.DS_AG) || 0; const dp_ag = parseInt(scores.DP_AG) || 0;
      const al = parseFloat(scores.DS_AL) || 0; const dp_al = parseInt(scores.DP_AL) || 0;
      const dg = parseFloat(scores.DS_DG) || 0; const dp_dg = parseInt(scores.DP_DG) || 0;
      const dl = parseFloat(scores.DS_DL) || 0; const dp_dl = parseInt(scores.DP_DL) || 0;

      const fmt = (v, dp, label) => {
        const color = v > 0.20 ? 'var(--red)' : 'var(--green)';
        const posStr = dp > 0 ? `+${dp}` : `${dp}`;
        return `
            <div style="display: flex; align-items: center; gap: 4px; background: rgba(255,255,255,0.02); padding: 2px 4px; border-radius: 2px; border: 1px solid rgba(255,255,255,0.03);">
               <span style="color: var(--dim); font-size: 0.55rem; width: 14px;">${label}</span>
               <span style="color: ${color}; font-family: 'JetBrains Mono', monospace; font-size: 0.65rem; font-weight: 600;">${v.toFixed(3)}</span>
               <span style="color: var(--dim); font-size: 0.5rem;">(${posStr})</span>
            </div>
          `;
      };

      return `
       <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 4px; margin-top: 4px; width: 100%;">
          ${fmt(ag, dp_ag, 'AG')}
          ${fmt(al, dp_al, 'AL')}
          ${fmt(dg, dp_dg, 'DG')}
          ${fmt(dl, dp_dl, 'DL')}
       </div>`;
    };

    // Build badge entries
    const spliceItems = [];
    const revelItems = [];
    const outBadges = [];

    // Splice: show Pathogenic/Likely Pathogenic only (clinical noise reduction)
    // Missense: show all classifications from ClinVar
    const isPathOrLP = (sigStr) =>
      sigStr.includes('pathogenic') && !sigStr.includes('conflicting');

    for (const badgeText in altMap) {
      const obj = altMap[badgeText];
      const sigs = Array.from(obj.sigs);
      const stars = obj.stars || '';
      let worst = 'Unknown', htmlCls = 'sig-vus';
      const sigStr = sigs.join(' ').toLowerCase();

      if (sigStr.includes('conflicting')) { worst = 'Conflicting'; htmlCls = 'sig-con'; }
      else if (sigStr.includes('likely pathogenic')) { worst = 'Likely Pathogenic'; htmlCls = 'sig-lp'; }
      else if (sigStr.includes('pathogenic')) { worst = 'Pathogenic'; htmlCls = 'sig-p'; }
      else if (sigStr.includes('uncertain')) { worst = 'VUS'; htmlCls = 'sig-vus'; }
      else if (sigStr.includes('likely benign')) { worst = 'Likely Benign'; htmlCls = 'sig-lb'; }
      else if (sigStr.includes('benign')) { worst = 'Benign'; htmlCls = 'sig-b'; }
      else worst = sigs[0] || 'Unknown';

      // Skip non-P/LP results for splice variants
      if (isSplice && !isPathOrLP(sigStr)) continue;

      const idx = outBadges.length;
      if (isSplice && obj.hgvsT) {
        outBadges.push({
          hgvsT: obj.hgvsT,
          title: badgeText,
          plp: (worst === 'Pathogenic' || worst === 'Likely Pathogenic'),
          html: `
          <div style="margin-bottom: 16px; width: 100%; border-bottom: 1px solid rgba(255,255,255,0.05); padding-bottom: 10px;" id="altSplice_${idx}">
            <div style="display: flex; justify-content: space-between; align-items: flex-start; gap: 8px;">
               <div style="flex: 1;">
                 <div style="color: var(--text-bright); font-weight: 800; font-size: 0.85rem; margin-bottom: 3px;">${badgeText}</div>
                 <div style="display: flex; align-items: center; gap: 5px;">
                    <span class="sig-badge ${htmlCls}" style="padding:2px 5px; font-size:0.68rem;" data-tippy-content="${sigs.join(', ')}">${worst}</span>
                    ${stars}
                 </div>
               </div>
               <div id="altSpliceCoords_${idx}" style="font-size: 0.65rem; color: var(--teal); font-family: 'JetBrains Mono', monospace; text-align: right; opacity: 0.8;"></div>
            </div>
            <div id="altSpliceSAI_${idx}" style="font-size: 0.68rem; color: var(--dim); margin-top: 8px;">
               <span style="font-style: italic; opacity: 0.6;">SpliceAI: loading...</span>
            </div>
          </div>`
        });
      } else {
        outBadges.push({
          hgvsT: obj.hgvsT,
          title: badgeText,
          plp: (worst === 'Pathogenic' || worst === 'Likely Pathogenic'),
          html: `
          <div style="margin-bottom: 12px; width: 100%; border-bottom: 1px solid rgba(255,255,255,0.05); padding-bottom: 8px;" id="altMiss_${idx}">
            <div style="display: flex; justify-content: space-between; align-items: flex-start; gap: 8px;">
               <div style="flex: 1;">
                 <div style="display: flex; align-items: center; gap: 8px; margin-bottom: 3px;">
                    <div style="color: var(--text-bright); font-weight: 800; font-size: 0.85rem;">${badgeText}</div>
                    <div id="altMissExpert_${idx}"></div>
                 </div>
                 <div style="display: flex; align-items: center; gap: 5px;">
                    <span class="sig-badge ${htmlCls}" style="padding:2px 5px; font-size:0.68rem;" data-tippy-content="${sigs.join(', ')}">${worst}</span>
                    ${stars}
                 </div>
               </div>
               <div style="text-align: right;">
                 <div id="altMissCoords_${idx}" style="font-size: 0.65rem; color: var(--teal); font-family: 'JetBrains Mono', monospace; opacity: 0.8; margin-bottom: 2px;"></div>
                 <div id="altMissRevel_${idx}" style="font-size: 0.65rem; color: var(--text-dim); font-weight: 600;">REVEL: &hellip;</div>
               </div>
            </div>
          </div>`
        });
      }
    }

    // ── RENDER immediately so element IDs exist before async enrichment fires ──
    const noneMsg = isSplice ? 'No P/LP variants found. Search ↗' : 'None found. Search ↗';
    const flexStyle = `display:flex; flex-direction:column; align-items:flex-start; gap:4px; width:100%;`;
    document.getElementById('cvAltCodon').innerHTML = outBadges.length
      ? `<div style="${flexStyle}">${outBadges.map(b => b.html).join('')}<div style="margin-top:4px;"><a href="${searchUrl}" target="_blank" style="color:var(--amber);text-decoration:none;font-weight:800;font-size:0.8rem;">Search Region ↗</a></div></div>`
      : `<a href="${searchUrl}" target="_blank" style="color:var(--dim);text-decoration:none;font-size:0.75rem;">${noneMsg}</a>`;

    // ── AUTO-SUGGEST PS1 / PM5 ────────────────────────────────────────────────
    if (!isSplice) {
      if (ps1Candidates.length > 0) {
        const ps1MaxStars = Math.max(...ps1Candidates.map(c => c.stars));
        data.clinvar.ps1Suggested = ps1MaxStars >= 2 ? 'PS1' : 'PS1_Moderate';
      } else {
        data.clinvar.ps1Suggested = null;
      }
      data.clinvar.pm5Suggested = outBadges.some(b => b.plp) ? 'PM5' : null;
      evaluateACMG();
    }

    // ── ASYNC ENRICHMENT (ClinGen Allele Registry + LDH) ─────────────────────
    const variantsToResolve = outBadges.map(b => b.hgvsT).filter(Boolean);
    if (variantsToResolve.length > 0) {

      // Two-step per variant: Allele Registry (HGVS → CAid + coords) then
      // LDH (CAid → REVEL MANE Select + expert classification).
      // All variants fire in parallel via Promise.all.
      const fetchClinGenVariantData = async (hgvsT) => {
        const clean = hgvsT.replace(/\([^)]*\)/g, '').replace(/[^\x20-\x7E]/g, '').trim();

        // Step A: Allele Registry — HGVS → CAid + GRCh38 coords (0-based → 1-based)
        let caid = null, coords = null;
        try {
          const arRes = await fetchWithTimeout(
            `https://reg.genome.network/allele?hgvs=${encodeURIComponent(clean)}`
          );
          if (arRes.ok) {
            const ar = await arRes.json();
            caid = ar.caid || ar['@id']?.split('/').pop() || null;
            const gAllele = ar.genomicAlleles?.find(g => g.referenceGenome === 'GRCh38');
            const c = gAllele?.coordinates?.[0];
            if (c && gAllele?.chromosome) {
              coords = {
                chr: String(gAllele.chromosome).replace('chr', ''),
                pos: String(c.start + 1),
                ref: normalizeVcfAllele(c.referenceAllele),
                alt: normalizeVcfAllele(c.allele)
              };
            }
          }
        } catch (e) { console.warn(`[ClinGen AR] Failed for ${hgvsT}:`, e.message); }

        if (!caid) return { caid: null, coords, revel: null, expert: null };

        // Step B: LDH — CAid → REVEL (MANE Select only) + expert classification
        let revel = null, expert = null;
        try {
          const ldhRes = await fetchWithTimeout(
            `https://ldh.clinicalgenome.org/ldh/Variant/id/${caid}?detail=high`
          );
          if (ldhRes.ok) {
            const ldh = await ldhRes.json();
            const ld = ldh.data?.ld || {};

            const revelEntries = ld.RevelScore?.[0]?.entContent || [];
            const maneEntry = revelEntries.find(e => e.mane === 'MANE Select');
            if (maneEntry?.score != null) revel = parseFloat(maneEntry.score);

            const pathInfo = ld.PathogenicityClassification?.[0]?.entContent
                          || ld.ClinicalInterpretation?.[0]?.entContent;
            if (pathInfo?.classification || pathInfo?.significance) {
              expert = {
                classification: pathInfo.classification || pathInfo.significance,
                vcep: pathInfo.source?.label,
                date: pathInfo.date
              };
            }
          }
        } catch (e) { console.warn(`[ClinGen LDH] Failed for ${caid}:`, e.message); }

        return { caid, coords, revel, expert };
      };

      (async () => {
        try {
          // Fire all ClinGen lookups in parallel
          const clingenResults = await Promise.all(
            outBadges.map(b => b.hgvsT ? fetchClinGenVariantData(b.hgvsT) : Promise.resolve({}))
          );

          // Update coords, REVEL, expert badge, and LDH link for each variant
          for (let i = 0; i < outBadges.length; i++) {
            const { coords, revel, caid, expert } = clingenResults[i] || {};

            if (coords) {
              const coordsEl = document.getElementById(
                isSplice ? `altSpliceCoords_${i}` : `altMissCoords_${i}`
              );
              if (coordsEl) coordsEl.innerText = `GRCh38: chr${coords.chr}:${coords.pos} ${coords.ref}>${coords.alt}`;
            }

            if (!isSplice) {
              const revelEl = document.getElementById(`altMissRevel_${i}`);
              if (revelEl) {
                if (revel !== null && revel !== undefined) {
                  const color = revel >= 0.773 ? 'var(--red)'
                              : revel >= 0.644 ? 'var(--amber)'
                              : revel >= 0.183 ? 'var(--dim)'
                              : 'var(--teal)';
                  revelEl.innerHTML = `REVEL: <span style="color:${color}; font-weight:700;">${revel.toFixed(3)}</span>`;
                }
                // else: leave as "…" — MyVariant.info fallback below will fill it in or set N/A
              }

              if (expert?.classification) {
                const expertEl = document.getElementById(`altMissExpert_${i}`);
                if (expertEl) expertEl.innerHTML = `<span style="background:linear-gradient(135deg,var(--amber),#FFD700);color:#000;padding:1px 4px;border-radius:3px;font-size:0.55rem;font-weight:900;text-transform:uppercase;letter-spacing:0.3px;box-shadow:0 0 5px rgba(255,215,0,0.3);" title="Expert Vetted by ${expert.vcep || 'ClinGen'} on ${expert.date || ''}">ClinGen Expert</span>`;
              }
            }
          }

          // REVEL fallback chain for missense variants where LDH returned no score:
          //   1. MyVariant.info (coordinate-based)
          //   2. UCSC revel{Alt} bigWig track (same source as the main variant REVEL)
          if (!isSplice) {
            // Collect variants still needing REVEL
            const missingRevelMap = {};
            for (let i = 0; i < outBadges.length; i++) {
              const { coords, revel } = clingenResults[i] || {};
              if ((revel === null || revel === undefined) && coords) {
                missingRevelMap[outBadges[i].hgvsT] = coords;
              }
            }

            // 1. MyVariant.info
            const myvariantScores = Object.keys(missingRevelMap).length > 0
              ? await fetchRevelByCoords(missingRevelMap)
              : {};

            // 2. UCSC fallback for anything MyVariant couldn't fill
            const ucscScores = {};
            const stillMissing = Object.entries(missingRevelMap)
              .filter(([hgvsT]) => !(hgvsT in myvariantScores));
            await Promise.all(stillMissing.map(async ([hgvsT, coords]) => {
              if (!coords || String(coords.alt).length !== 1) return;
              const track = `revel${coords.alt.toUpperCase()}`;
              const start = parseInt(coords.pos, 10) - 1;
              const end   = parseInt(coords.pos, 10);
              try {
                const r = await fetchWithTimeout(
                  `https://api.genome.ucsc.edu/getData/track?genome=hg38&track=${track}&chrom=chr${coords.chr}&start=${start}&end=${end}`
                );
                if (!r.ok) return;
                const j = await r.json();
                const val = j[track]?.[0]?.value;
                if (val != null) ucscScores[hgvsT] = parseFloat(val);
              } catch { /* non-blocking */ }
            }));

            // Apply scores and render
            const applyRevel = (revelEl, score) => {
              if (score !== null && score !== undefined) {
                const color = score >= 0.773 ? 'var(--red)'
                            : score >= 0.644 ? 'var(--amber)'
                            : score >= 0.183 ? 'var(--dim)'
                            : 'var(--teal)';
                revelEl.innerHTML = `REVEL: <span style="color:${color}; font-weight:700;">${score.toFixed(3)}</span>`;
              } else {
                revelEl.innerText = 'REVEL: N/A';
              }
            };

            for (let i = 0; i < outBadges.length; i++) {
              const hgvsT = outBadges[i].hgvsT;
              const { revel } = clingenResults[i] || {};
              if (revel !== null && revel !== undefined) continue; // already rendered from LDH
              const revelEl = document.getElementById(`altMissRevel_${i}`);
              if (!revelEl) continue;
              if (hgvsT && hgvsT in myvariantScores) {
                applyRevel(revelEl, myvariantScores[hgvsT]);
              } else if (hgvsT && hgvsT in ucscScores) {
                applyRevel(revelEl, ucscScores[hgvsT]);
              } else {
                revelEl.innerText = 'REVEL: N/A';
              }
            }
          }

          // SpliceAI for splice variants (sequential, rate-limited)
          if (isSplice) {
            for (let i = 0; i < outBadges.length; i++) {
              const { coords } = clingenResults[i] || {};
              if (!coords) continue;
              try {
                const scores = await fetchSpliceAIBroad(coords, data.coords.gene, data.ensembl.transcript);
                const el = document.getElementById(`altSpliceSAI_${i}`);
                if (el) el.innerHTML = formatSAITable(scores);
              } catch (e) {
                console.warn(`[AltSplice] SpliceAI failed for badge ${i}:`, e);
              }
              if (i < outBadges.length - 1) await new Promise(r => setTimeout(r, 1500));
            }
          }
        } catch (e) {
          console.error('[AltMiss] Async enrichment error:', e);
        }
      })();
    }

  } catch (e) {
    if (e.name === 'AbortError') return;
    console.error('fetchClinVarCodon:', e);
    document.getElementById('cvAltCodon').innerHTML = `<a href="${searchUrl}" target="_blank" style="color:var(--dim);text-decoration:none;">Search ↗</a>`;
  }
}


async function fetchNcbiCount(term) {
  try { const r = await fetchWithRetry(`${NCBI_BASE}/esearch.fcgi?db=pubmed&retmode=json&term=${term}&api_key=${NCBI_KEY}&email=${NCBI_EMAIL}`); if (r.ok) { const j = await ncbiSafeJson(r); return parseInt(j.esearchresult?.count || 0); } return 0; } catch (e) { if (e.name === 'AbortError') throw e; return 0; }
}

// Resolves MONDO disease ID to OMIM xref via Monarch entity endpoint.
async function resolveMondoToOmim(mondoId) {
  if (!mondoId || !mondoId.startsWith('MONDO:')) return mondoId?.startsWith('OMIM:') ? mondoId : null;
  try {
    const r = await fetchWithTimeout(`https://api.monarchinitiative.org/v3/api/entity/${encodeURIComponent(mondoId)}`);
    if (!r.ok) return null;
    const j = await r.json();
    const omim = (j.xref || []).find(x => x.startsWith('OMIM:'));
    return omim || null;
  } catch (e) { if (e.name === 'AbortError') throw e; return null; }
}

// Fetch HPO terms for a disease from JAX HPO (requires OMIM ID).
async function fetchJaxHpoDisease(omimId) {
  const hpSet = new Set();
  if (!omimId) return hpSet;
  try {
    const r = await fetchWithTimeout(`https://ontology.jax.org/api/network/annotation/${encodeURIComponent(omimId)}`);
    if (r.ok) {
      const j = await r.json();
      Object.values(j.categories || {}).forEach(arr => {
        arr.forEach(p => {
          if (p.id && typeof p.id === 'string' && p.id.startsWith('HP:')) {
            hpSet.add(p.id);
          }
        });
      });
      console.log(`[JAX detail] ${omimId}: found terms:`, Array.from(hpSet));
    }
  } catch (e) { if (e.name === 'AbortError') throw e; }
  return hpSet;
}

// Fetch HPO associations for a disease from Monarch (with object_closure for ancestor matching).
// Populates dSet and depthMap in place — depth 0 = direct annotation.
async function fetchMonarchDiseaseHpo(diseaseId, dSet, depthMap) {
  const PAGE = 500;
  let totalAdded = 0;
  for (let offset = 0; offset <= 1500; offset += PAGE) {
    try {
      const ar = await fetchWithTimeout(`https://api.monarchinitiative.org/v3/api/association/all?category=biolink:DiseaseToPhenotypicFeatureAssociation&subject=${encodeURIComponent(diseaseId)}&limit=${PAGE}&offset=${offset}`);
      if (!ar.ok) break;
      const aj = await ar.json();
      if (!aj.items || aj.items.length === 0) break;
      for (const item of aj.items) {
        if (item.object?.startsWith('HP:')) {
          dSet.add(item.object);
          if (!depthMap.has(item.object)) depthMap.set(item.object, 0);
          totalAdded++;
        }
        if (Array.isArray(item.object_closure)) {
          for (let i = 0; i < item.object_closure.length; i++) {
            const hp = item.object_closure[i];
            if (typeof hp === 'string' && hp.startsWith('HP:')) {
              dSet.add(hp);
              if (!depthMap.has(hp) || depthMap.get(hp) > i) depthMap.set(hp, i);
              totalAdded++;
            }
          }
        }
      }
      if (aj.items.length < PAGE) break;
    } catch (e) { if (e.name === 'AbortError') throw e; break; }
  }
  console.log(`[Monarch HPO detailed] ${diseaseId}: added ${totalAdded} terms, final set size: ${dSet.size}`);
}

// Mirror of gene-phenotype matching: JAX HPO (primary, depth 0) + Monarch closure (secondary).
// Returns { hpoSet: Set<HP:id>, depthMap: Map<HP:id, depth> }
async function fetchDiseaseHPOTerms(diseaseId) {
  const hpoSet = new Set();
  const depthMap = new Map();
  if (!diseaseId) { console.log('[HPO] no diseaseId'); return { hpoSet, depthMap }; }

  console.log(`[HPO] fetchDiseaseHPOTerms(${diseaseId})`);

  // Step 1: JAX HPO disease annotation (primary — requires OMIM)
  const omimId = await resolveMondoToOmim(diseaseId);
  console.log(`[HPO] MONDO→OMIM: ${diseaseId} → ${omimId}`);
  if (omimId) {
    const jaxHpo = await fetchJaxHpoDisease(omimId);
    console.log(`[HPO] JAX HPO (${omimId}): ${jaxHpo.size} terms`);
    jaxHpo.forEach(hp => { hpoSet.add(hp); depthMap.set(hp, 0); });
  }

  // Step 2: Monarch with object_closure (provides ancestor terms for broader matching)
  const beforeMonarch = hpoSet.size;
  await fetchMonarchDiseaseHpo(diseaseId, hpoSet, depthMap);
  console.log(`[HPO] Monarch HPO (${diseaseId}): ${hpoSet.size - beforeMonarch} new terms (total ${hpoSet.size})`);

  // Log final set contents to debug matching
  const hasCardiomyopathy = hpoSet.has('HP:0001638');
  const hasArrhythmia = hpoSet.has('HP:0011675');
  console.log(`[HPO final] ${diseaseId}: has Cardiomyopathy? ${hasCardiomyopathy}, has Arrhythmia? ${hasArrhythmia}`);

  return { hpoSet, depthMap };
}

// Score definitive conditions by patient HPO phenotype overlap.
// Returns { diseases: string[], scores: Map<name, count> }
// Rule: if maxScore >= 2, return ALL conditions tied at maxScore.
//       Otherwise return the single best (or all if no patient phenotypes available).
async function selectDiseasesByPhenoMatch(conditions) {
  if (!conditions || conditions.length === 0) return { diseases: [], scores: new Map() };

  // No patient HPO data — fall back to all passed conditions
  if (!data.ptPhenoGroups || data.ptPhenoGroups.length === 0) {
    return { diseases: conditions.map(c => c.diseaseName), scores: new Map() };
  }

  const ptFlat = data.ptPhenoGroups.flat();
  console.log('[HPO Match] patient texts:', data.ptPhenoTexts);
  console.log('[HPO Match] patient phenotype groups:', JSON.stringify(data.ptPhenoGroups));
  console.log('[HPO Match] flattened patient HP IDs:', ptFlat);
  if (ptFlat.length === 0) console.warn('[HPO Match] WARNING: No patient HP IDs resolved!');

  const scored = await Promise.all(conditions.map(async c => {
    const { hpoSet, depthMap } = await fetchDiseaseHPOTerms(c.diseaseId);
    const matchDetails = [];
    let matchCount = 0;

    // Count patient phenotype groups where at least one HP ID matches within depth ≤ 3
    data.ptPhenoTexts?.forEach((text, idx) => {
      const group = data.ptPhenoGroups[idx] || [];
      console.log(`  [match] ${text}: checking group`, group, 'against hpoSet size', hpoSet.size, 'depthMap size', depthMap.size);
      let bestDepth = 999, bestHp = null;
      for (const id of group) {
        const inSet = hpoSet.has(id);
        const depth = depthMap.get(id);
        console.log(`    checking ${id}: inSet=${inSet}, depth=${depth}`);
        if (inSet) {
          const d = depth ?? 0;
          // For disease HPOs: accept any ancestor depth (patient phenotypes are broad)
          if (d < bestDepth) { bestDepth = d; bestHp = id; if (d === 0) break; }
        }
      }
      if (bestHp) {
        matchCount++;
        const tag = bestDepth === 0 ? 'direct' : `ancestor d${bestDepth}`;
        matchDetails.push(`${text} (${bestHp}, ${tag})`);
        console.log(`    ✓ MATCH: ${bestHp} at depth ${bestDepth}`);
      } else {
        console.log(`    ✗ no match`);
      }
    });
    console.log(`[HPO Match] ${c.diseaseName} (${c.diseaseId}): ${hpoSet.size} HPOs, ${matchCount} matches`, matchDetails);
    return { diseaseName: c.diseaseName, matchCount, matchDetails };
  }));

  const scoreMap = new Map(scored.map(s => [s.diseaseName, s.matchCount]));
  const detailsMap = new Map(scored.map(s => [s.diseaseName, s.matchDetails]));
  const maxScore = Math.max(...scored.map(s => s.matchCount));

  let selected;
  if (maxScore >= 2) {
    // All diseases tied at the best score qualify
    selected = scored.filter(s => s.matchCount === maxScore).map(s => s.diseaseName);
  } else {
    // Single best (even if score is 0 or 1)
    const best = scored.reduce((a, b) => b.matchCount > a.matchCount ? b : a);
    selected = [best.diseaseName];
  }

  return { diseases: selected, scores: scoreMap, details: detailsMap };
}

async function updatePubMedHitCount(query, elementId) {
  if (!query) return;
  try {
    const count = await fetchNcbiCount(query);
    const el = document.getElementById(elementId);
    if (!el) return;
    const color = count === 0 ? 'var(--dim)' : count >= 10 ? 'var(--teal)' : 'var(--amber)';
    const url = `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(query)}`;
    el.innerHTML = `<a href="${url}" target="_blank" style="color:${color}; text-decoration:none; font-weight:700; font-size:0.85rem;">${count.toLocaleString()}</a><span style="color:var(--dim); font-size:0.7rem; margin-left:4px;">articles</span>`;
  } catch (e) {
    console.warn(`[PubMed Hit Count] Failed for ${elementId}:`, e.message);
  }
}

async function fetchEnsemblVersion() {
  try { const r = await fetch('https://rest.ensembl.org/info/software?content-type=application/json'); if (r.ok) { const j = await r.json(); if (j.release) document.getElementById('apiStatus').textContent = `APIs: VV | Ensembl Rel ${j.release} | MyVariant.info | NCBI E-utils (auth) | Monarch HPO`; } } catch (e) { }
}

// --- END OF CLIN GEN DICTIONARY LOGIC ---

/**
 * Parses patient phenotype strings into HPO IDs with optimized search limits.
 * @param {string} textStr - Comma-separated phenotype terms.
 */
/**
 * Resolves a single patient phenotype term to HPO IDs (+ optional MeSH term).
 * JAX→Monarch stay sequential (Monarch is conditional on JAX's result), but the
 * independent MeSH lookup runs concurrently with that chain. Returns { ids, mesh }.
 */
async function resolveOnePhenotype(t) {
  // Cache hit returns ids only (MeSH was never cached — matches prior behaviour).
  if (phenotypeCache[t]) return { ids: phenotypeCache[t], mesh: null };

  const possibleIds = new Set();

  // MeSH translation — independent, fire concurrently with JAX/Monarch.
  const meshPromise = (async () => {
    try {
      const meshRes = await fetchWithRetry(`${NCBI_BASE}/esearch.fcgi?db=mesh&term=${encodeURIComponent(t)}&retmode=json&api_key=${NCBI_KEY}`);
      if (meshRes.ok) {
        const mj = await meshRes.json();
        const trans = mj.esearchresult?.translationset;
        if (trans && trans.length > 0) {
          const match = trans[0].to?.match(/"([^"]+)"\[MeSH/i);
          if (match) return match[1];
        }
      }
    } catch (e) { if (e.name === 'AbortError') throw e; }
    return null;
  })();

  try {
    // 1. JAX HPO Search - limit 20 with exact-match prioritization
    const jaxRes = await fetch(`https://ontology.jax.org/api/hp/search?q=${encodeURIComponent(t)}&page=0&limit=20`);
    if (jaxRes.ok) {
      const json = await jaxRes.json();
      if (json.terms) {
        // Normalise Unicode dashes to hyphen-minus so autocomplete-inserted
        // en-dashes / em-dashes don't break exact-match detection.
        const normDash = s => s.toLowerCase().replace(/[‐-―−]/g, '-');
        const tl = normDash(t);
        const sorted = [...json.terms].sort((a, b) => {
          const aExact = normDash(a.name) === tl ? 1 : 0;
          const bExact = normDash(b.name) === tl ? 1 : 0;
          return bExact - aExact;
        });
        // Exact match: primary name OR any synonym → use only that one term.
        const isExact = term =>
          normDash(term.name) === tl ||
          (term.synonyms || []).some(s => normDash(s) === tl);
        const exactMatch = sorted.find(isExact);
        if (exactMatch) possibleIds.add(exactMatch.id);
      }
    }

    // Monarch removed: JAX exact-match (name + synonyms) is the sole authority.
    // Fuzzy Monarch results caused partial/mistyped terms to get wrong IDs.

    // Colloquial mapping fallbacks
    const tl = t.toLowerCase();
    if (tl.includes('bone deformity')) possibleIds.add('HP:0000924');
    if (tl === 'fractures' || tl === 'fracture') {
      possibleIds.add('HP:0002757');
      possibleIds.add('HP:0002659');
    }
  } catch (e) {
    if (e.name === 'AbortError') throw e; // propagate — don't cache partial results
    console.error("Phenotype parsing error:", e);
  }

  let mesh = null;
  try { mesh = await meshPromise; } catch (e) { if (e.name === 'AbortError') throw e; }

  // Only cache successful (non-empty) resolutions so aborted fetches don't poison the cache
  if (possibleIds.size > 0) phenotypeCache[t] = Array.from(possibleIds);
  return { ids: Array.from(possibleIds), mesh };
}

async function parsePatientPhenotypes(textStr) {
  const terms = textStr.split(',').map(s => s.trim()).filter(s => s && s !== '.' && s !== '-');
  data.ptPhenoMeSH = {};

  // Resolve all terms in parallel — each term's lookups are independent.
  const results = await Promise.all(terms.map(t => resolveOnePhenotype(t)));

  // Set texts and groups atomically so renderHpoPills never sees a mismatched pair.
  data.ptPhenoTexts = terms;
  data.ptPhenoGroups = results.map(r => r.ids);
  results.forEach((r, i) => { if (r.mesh) data.ptPhenoMeSH[terms[i]] = r.mesh; });
  // Re-trigger multi-transcript literature search with updated phenotype terms.
  // _scheduleLitNotations key-change detection ensures it only re-fires if terms changed.
  if (data.coords?.gene) _scheduleLitNotations(window.currentSearchGen());
}

/**
 * Resolves gene symbol to NCBI Entrez Gene ID (human) via NCBI Gene eSearch.
 */
async function resolveEntrezId(gene) {
  try {
    const r = await fetchWithRetry(`${NCBI_BASE}/esearch.fcgi?db=gene&term=${encodeURIComponent(gene)}%5Bgene%5D+AND+"Homo+sapiens"%5Borganism%5D&retmode=json&api_key=${NCBI_KEY}`);
    if (r.ok) {
      const j = await r.json();
      return j.esearchresult?.idlist?.[0] || null;
    }
  } catch (e) { if (e.name === 'AbortError') throw e; }
  return null;
}

/**
 * Fetches HPO phenotype annotations for a gene from JAX HPO (uses Entrez ID).
 * Returns a Set of HP: term IDs directly annotated to this gene.
 */
async function fetchJaxHpoGene(entrezId) {
  const hpSet = new Set();
  try {
    const r = await fetchWithRetry(`https://hpo.jax.org/api/hpo/gene/${entrezId}`);
    if (r.ok) {
      const j = await r.json();
      (j.phenotypes || []).forEach(p => { if (p.ontologyId?.startsWith('HP:')) hpSet.add(p.ontologyId); });
    }
  } catch (e) { if (e.name === 'AbortError') throw e; }
  return hpSet;
}

/**
 * Fetches HPO associations from Monarch using NCBIGene or HGNC subject ID.
 * Populates gSet and depthMap in place.
 */
async function fetchMonarchHpo(subjectId, gSet, depthMap) {
  const PAGE = 500;
  for (let offset = 0; offset <= 1500; offset += PAGE) {
    try {
      const ar = await fetch(`https://api.monarchinitiative.org/v3/api/association/all?category=biolink:GeneToPhenotypicFeatureAssociation&subject=${encodeURIComponent(subjectId)}&limit=${PAGE}&offset=${offset}`);
      if (!ar.ok) break;
      const aj = await ar.json();
      if (!aj.items || aj.items.length === 0) break;
      for (const item of aj.items) {
        if (item.object?.startsWith('HP:')) {
          gSet.add(item.object);
          if (!depthMap.has(item.object)) depthMap.set(item.object, 0);
        }
        if (Array.isArray(item.object_closure)) {
          item.object_closure.forEach((hp, i) => {
            if (hp.startsWith('HP:')) {
              gSet.add(hp);
              if (!depthMap.has(hp) || depthMap.get(hp) > i) depthMap.set(hp, i);
            }
          });
        }
      }
      if (aj.items.length < PAGE) break;
    } catch (e) { if (e.name === 'AbortError') throw e; break; }
  }
}

/**
 * Matches gene-associated phenotypes against patient terms.
 * Sources: JAX HPO gene annotations (primary) + Monarch via NCBIGene ID (secondary).
 * Literature evidence is merged in later by runPubMedOnly.
 */
async function runPhenoMatch(gene) {
  // Re-trigger multi-transcript PubMed search with updated phenotype/disease terms
  _scheduleLitNotations(window.currentSearchGen());
  const pScoreEl = document.getElementById('pScore');
  if (pScoreEl) pScoreEl.innerText = '…';
  const pLinEl = document.getElementById('pLin');
  if (pLinEl) pLinEl.innerText = '…';
  if (!gene || gene === '-') return;
  try {
    // Await phenotype resolution that was fired concurrently with fetchVV.
    // In most cases it has already resolved by the time we reach here.
    if (window._phenoParsePromise) await window._phenoParsePromise;

    // Step 1: Resolve Entrez ID via NCBI Gene (reliable, human-specific)
    const entrezId = await resolveEntrezId(gene);

    const gSet = new Set(), depthMap = new Map();

    if (entrezId) {
      // Steps 2+3: JAX HPO and Monarch both need entrezId but are independent —
      // fire them in parallel and merge results once both resolve.
      const [jaxHpo] = await Promise.all([
        fetchJaxHpoGene(entrezId),
        fetchMonarchHpo(`NCBIGene:${entrezId}`, gSet, depthMap)
      ]);
      jaxHpo.forEach(hp => { gSet.add(hp); depthMap.set(hp, 0); });
    }

    // Step 4: Monarch via HGNC ID fallback (if Ensembl lookup succeeds)
    if (gSet.size === 0) {
      try {
        const er = await fetch(`https://rest.ensembl.org/lookup/symbol/homo_sapiens/${gene}?content-type=application/json`);
        if (er.ok) {
          const ej = await er.json();
          const m = ej.description?.match(/HGNC:(\d+)/);
          if (m) await fetchMonarchHpo(`HGNC:${m[1]}`, gSet, depthMap);
        }
      } catch (e) { if (e.name === 'AbortError') throw e; }
    }

    // Step 5: Monarch gene search fallback — human-only strict filter
    if (gSet.size === 0) {
      try {
        const res = await fetch(`https://api.monarchinitiative.org/v3/api/search?q=${encodeURIComponent(gene)}&category=biolink:Gene&limit=10`);
        if (res.ok) {
          const json = await res.json();
          for (const item of (json.items || [])) {
            const isHuman = item.id?.startsWith('HGNC:') ||
              item.in_taxon?.some(t => t.includes('9606')) ||
              item.in_taxon_label?.toLowerCase().includes('homo sapiens');
            if (isHuman) {
              (item.has_phenotype || []).forEach(hp => {
                if (hp.startsWith('HP:')) { gSet.add(hp); depthMap.set(hp, 0); }
              });
              break;
            }
          }
        }
      } catch (e) { if (e.name === 'AbortError') throw e; }
    }

    // Step 6: Match patient HPO IDs against gene HPO set
    let score = 0, matched = [], directCount = 0, ancestorCount = 0;
    const hpoMatchedTerms = new Set();
    for (let i = 0; i < data.ptPhenoTexts.length; i++) {
      const ids = data.ptPhenoGroups[i] || [];
      let best = 999, hit = false;
      for (const id of ids) {
        if (gSet.has(id)) {
          const d = depthMap.get(id) ?? 0;
          if (d <= 3 && d < best) { best = d; hit = true; if (d === 0) break; }
        }
      }
      if (hit) {
        score++;
        matched.push(`${data.ptPhenoTexts[i]}${best === 0 ? '(D)' : `(C${best})`}`);
        hpoMatchedTerms.add(data.ptPhenoTexts[i]);
        if (best === 0) directCount++; else ancestorCount++;
      }
    }

    // Store for literature merge in runPubMedOnly
    data.hpoMatchedTerms = hpoMatchedTerms;
    data.hpoMatchResult = { score, total: data.ptPhenoTexts.length, matchedArr: matched, directCount, ancestorCount };

    renderPhenoAndLitCards({ score, total: data.ptPhenoTexts.length, matchedArr: matched, directCount, ancestorCount });

    // Step 7: OMIM synopsis corroboration (independent per-disease match).
    // Wait for the local cache if it hasn't finished loading yet (race on cold start).
    // ptPhenoGroups / ptPhenoTexts are passed explicitly — no implicit data.* reads.
    const ptGroups = data.ptPhenoGroups;
    const ptTexts  = data.ptPhenoTexts;
    if (!omimLocalCache.initialized) await omimLocalReady();
    // Load the HPO phenotype-fit tables (background/parents/disease-frequencies).
    // Resolves once; degrades to null if the tables aren't built yet. Used both by
    // the OMIM matcher (IC-floored ontology overlap) and the LR scorer below.
    const hpoFit = (typeof window !== 'undefined' && window.HpoFit)
      ? await window.HpoFit.load().then(() => (window.HpoFit.ready ? window.HpoFit : null)).catch(() => null)
      : null;
    const omimResults = matchOmimSynopsis(gene, ptGroups, ptTexts, hpoFit);
    data.omimSynopsisResults = omimResults;
    renderOmimSynopsisSummary(omimResults);
    applyOmimSynopsisBadges(omimResults, ptTexts.length);

    // Step 8: HPO phenotype-fit LR (ontology-aware, background-frequency-weighted).
    // Owns data.phenotypeFit exclusively; reads only explicit params + the loaded
    // tables. Non-blocking: any failure leaves ready:false and never cascades.
    try {
      data.phenotypeFit = scorePhenotypeFit(gene, ptGroups, hpoFit);
    } catch (e) {
      data.phenotypeFit = { ready: false, results: [] };
    }
    if (typeof renderPhenotypeFitSummary === 'function') renderPhenotypeFitSummary(data.phenotypeFit);

    runPubMedOnly(gene);
  } catch (e) {
    if (e.name !== 'AbortError') console.error("Match error:", e);
  }
}

/**
 * scorePhenotypeFit(gene, ptGroups, hpoFit)
 * Inputs: gene symbol, patient HPO groups (explicit), loaded window.HpoFit (or null).
 *   Reads the gene's confirmed OMIM diseases from omimLocalCache.geneDiseaseMims and
 *   the per-disease frequency tables held by HpoFit. No data.* reads.
 * Outputs: { ready, results:[{ mimId, name, log10LR, totalLR, posterior, coverage,
 *   completeness, matched, informative, isSusceptibility }] } sorted by LR desc.
 *   LR (prior-independent) is the primary rank key; posterior is fixed-prior context only.
 * Failures: returns { ready:false, results:[] } if tables unloaded / no candidates.
 */
function scorePhenotypeFit(gene, ptGroups, hpoFit) {
  if (!hpoFit || !hpoFit.ready || !gene) return { ready: false, results: [] };
  const present = [...new Set((ptGroups || []).flat().filter(Boolean))];
  if (!present.length) return { ready: false, results: [] };

  const candidates = (omimLocalCache.geneDiseaseMims.get(gene.toUpperCase()) || [])
    .filter(d => d.type !== 'non-disease' && d.evidenceCode >= 3);
  if (!candidates.length) return { ready: false, results: [] };

  const results = [];
  for (const d of candidates) {
    const s = hpoFit.scoreDisease(`OMIM:${d.mimId}`, present, 1e-3);
    if (!s || s.matched === 0) continue;   // disease not in the frequency table, or no overlap
    results.push({
      mimId: d.mimId,
      name: d.name,
      isSusceptibility: d.type === 'susceptibility',
      log10LR: s.sumLn / Math.LN10,
      totalLR: s.totalLR,
      posterior: s.posterior,
      coverage: s.coverage,
      completeness: s.completeness,
      matched: s.matched,
      informative: s.informative
    });
  }
  // Rank by likelihood ratio (prior cancels within one gene's candidate set).
  results.sort((a, b) => b.log10LR - a.log10LR);
  return { ready: results.length > 0, results };
}

async function fetchNcbiCountClinVar(term) {
  try {
    const r = await fetchWithRetry(`${NCBI_BASE}/esearch.fcgi?db=clinvar&retmode=json&term=${term}&api_key=${NCBI_KEY}&email=${NCBI_EMAIL}`);
    if (r.ok) { const j = await ncbiSafeJson(r); return parseInt(j.esearchresult?.count || 0); }
    return 0;
  } catch { return 0; }
}

async function runPubMedOnly(gene) {
  if (!gene || gene === '-' || !data.ptPhenoTexts.length) {
    renderPhenoBreakdown([], gene);
    return;
  }
  document.getElementById('pPhenoBreakdown').innerHTML =
    '<div style="font-size:0.7rem;opacity:0.5;margin-top:4px;">Loading per-phenotype counts…</div>';

  // ── Parallel strategy ────────────────────────────────────────────────────────
  // All phenotype terms run concurrently via Promise.all; within each term the
  // three NCBI sub-queries (TIAB, MeSH, ClinVar) also run concurrently.
  // Rate limiting is handled centrally by the NCBI token bucket in fetchWithRetry
  // (≤9 req/s), so no manual per-term stagger is needed here.
  // ────────────────────────────────────────────────────────────────────────────
  const phenoResults = await Promise.all(
    data.ptPhenoTexts.map(async (t, i) => {
      const enc      = encodeURIComponent(t);
      const tiabQ    = `(${gene}%5Btiab%5D)+AND+(%22${enc}%22%5Btiab%5D)`;
      const meshTerm = data.ptPhenoMeSH?.[t] || null;
      const meshQ    = meshTerm
        ? `(${gene}%5Btiab%5D)+AND+(%22${encodeURIComponent(meshTerm)}%22%5BMeSH+Terms%5D)`
        : null;
      const cvQ      = `${gene}%5Bgene%5D+AND+%22${enc}%22%5Bdis%5D` +
                       `+AND+(Pathogenic%5Bclnsig%5D+OR+%22Likely+pathogenic%22%5Bclnsig%5D)`;

      // Fire TIAB + MeSH + ClinVar concurrently — no sequential round-trip waits
      const [tiabCount, meshCount, cvCount] = await Promise.all([
        fetchNcbiCount(tiabQ),
        meshQ ? fetchNcbiCount(meshQ) : Promise.resolve(0),
        fetchNcbiCountClinVar(cvQ),
      ]);

      return { term: t, tiab: tiabCount, tiabQ, mesh: meshCount, meshQ, meshTerm, clinvar: cvCount, cvQ };
    })
  );

  renderPhenoBreakdown(phenoResults, gene);
  mergeLiteratureIntoMatch(phenoResults);
  renderLitSearchLinks();
}

/**
 * Merges per-phenotype literature evidence into the phenotype match score.
 * Adds (L) tags for terms supported by PubMed/ClinVar but not matched via HPO.
 * Re-renders the Phenotype & Gene Match section with the combined result.
 */
function mergeLiteratureIntoMatch(phenoResults) {
  const hpoResult = data.hpoMatchResult || { score: 0, total: data.ptPhenoTexts.length, matchedArr: [], directCount: 0, ancestorCount: 0 };
  const hpoTerms = data.hpoMatchedTerms || new Set();

  let litCount = 0;
  const litMatched = [];
  for (const r of phenoResults) {
    if (!hpoTerms.has(r.term) && (r.tiab > 0 || r.clinvar > 0)) {
      litCount++;
      // Strong lit evidence if ClinVar P/LP > 0, else moderate
      const tag = r.clinvar > 0 ? '(L+)' : '(L)';
      litMatched.push(`${r.term}${tag}`);
    }
  }

  if (litCount === 0 && hpoResult.score === hpoResult.total) return; // nothing to update
  if (litCount === 0) return; // no new literature matches to add

  const combinedScore = hpoResult.score + litCount;
  const combinedMatched = [...hpoResult.matchedArr, ...litMatched];

  renderPhenoAndLitCards({
    score: combinedScore,
    total: hpoResult.total,
    matchedArr: combinedMatched,
    directCount: hpoResult.directCount,
    ancestorCount: hpoResult.ancestorCount,
    litCount
  });
}


/**
 * PDF Text Extraction using pdf.js
 */
async function extractTextFromPDF(file) {
  // Point pdf.js at the locally-vendored worker (offline; no CDN). Set once.
  if (pdfjsLib?.GlobalWorkerOptions && !pdfjsLib.GlobalWorkerOptions.workerSrc) {
    pdfjsLib.GlobalWorkerOptions.workerSrc = 'libs/pdf.worker.min.js';
  }
  const arrayBuffer = await file.arrayBuffer();
  const pdf = await pdfjsLib.getDocument({ data: arrayBuffer }).promise;
  let fullText = "";
  for (let i = 1; i <= pdf.numPages; i++) {
    const page = await pdf.getPage(i);
    const textContent = await page.getTextContent();
    const pageText = textContent.items.map(s => s.str).join(' ');
    fullText += `[Page ${i}] ${pageText}\n\n`;
  }
  return fullText;
}

async function handleFileUploads(event) {
  const files = event.target.files;
  if (!files.length) return;

  const listEl = document.getElementById('uploadedFileList');
  listEl.innerHTML = '<i>Reading PDFs...</i>';
  data.extractedPaperTexts = [];

  for (let file of files) {
    try {
      const text = await extractTextFromPDF(file);
      data.extractedPaperTexts.push({ name: file.name, content: text });
      const chip = document.createElement('span');
      chip.className = 'file-chip';
      chip.textContent = file.name;
      listEl.appendChild(chip);
    } catch (e) {
      console.error("PDF Read Error:", e);
    }
  }
  listEl.querySelector('i')?.remove();
}




/**
 * Gene-Disease Association Fetcher
 * Sources: Monarch Initiative v3 (aggregates OMIM + ClinGen + Orphanet)
 *          + PanelApp (Genomics England)
 *
 * Note: ClinGen direct API is CORS-blocked from file:// origin.
 * ClinGen and OMIM data are accessed via Monarch's aggregation layer,
 * identifiable by primary_knowledge_source containing "clingen" or "omim".
 */
// ── CLINGEN GENE-DISEASE VALIDITY MODULE ─────────────────────────────────
// Session-level cache of the full ClinGen Gene-Disease Validity CSV.
// The download has no CORS headers, so it is routed through the CodeTabs
// proxy already used elsewhere (see fetchLitVar). ~1 MB, fetched once.
let _clingenValidityCsvPromise = null;

/**
 * Parses a single CSV line, honouring double-quoted fields (which may
 * contain commas). ClinGen's export quotes every field.
 */
function _parseCsvLine(line) {
  const out = [];
  let cur = '', inQuotes = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (inQuotes) {
      if (ch === '"') {
        if (line[i + 1] === '"') { cur += '"'; i++; }  // escaped quote
        else inQuotes = false;
      } else cur += ch;
    } else {
      if (ch === '"') inQuotes = true;
      else if (ch === ',') { out.push(cur); cur = ''; }
      else cur += ch;
    }
  }
  out.push(cur);
  return out;
}

function _loadClingenValidityCsv() {
  if (_clingenValidityCsvPromise) return _clingenValidityCsvPromise;
  // Primary: Go helper (/api/clingen-validity) — direct fetch from the user's machine, no CORS issue.
  // Fallback: CodeTabs proxy — used when the helper isn't running (e.g. plain browser open without binary).
  // Rule 11: fetchWithTimeout for fire-and-forget remote calls.
  const helperUrl  = '/api/clingen-validity';
  const proxyUrl   = `https://api.codetabs.com/v1/proxy?quest=${encodeURIComponent('https://search.clinicalgenome.org/kb/gene-validity/download')}`;

  _clingenValidityCsvPromise = fetchWithTimeout(helperUrl, {}, 40000)
    .then(res => {
      // The Go helper returns non-JSON when unavailable (Vite's 404 page or a 502).
      const ct = res.headers.get('content-type') || '';
      if (!res.ok || ct.includes('text/html')) throw new Error('helper unavailable');
      return res.text();
    })
    .catch(() => {
      console.warn('[GeneValidity] Go helper unavailable — falling back to CodeTabs proxy');
      return fetchWithTimeout(proxyUrl, {}, 40000)
        .then(res => { if (!res.ok) throw new Error(`HTTP ${res.status}`); return res.text(); });
    })
    .catch(err => {
      console.warn('[GeneValidity] CSV load failed:', err.message);
      _clingenValidityCsvPromise = null;   // allow a later retry
      return '';
    });
  return _clingenValidityCsvPromise;
}

/**
 * fetchGeneValidity(geneSymbol)
 * Inputs: geneSymbol (explicit parameter — does NOT read data.*)
 * Outputs: data.geneValidity.{curations, gene, fetched}
 *   curations: [{disease, mondo, moi, sop, classification, date, gcep, url}]
 * Side effects: none (no DOM; consumed by copyReport).
 * Failures: Silent (non-blocking). Logs warning. Leaves curations = [].
 */
/**
 * expandClingenMoi(code) — map a ClinGen gene-validity MOI code to a readable label.
 * ClinGen's validity CSV stores MOI as short codes (AD, AR, XL, XLD, XLR, YL, SD = Semidominant,
 * MT = Mitochondrial, UD = Undetermined). Returns the full label; unrecognised input is returned
 * unchanged. PURE (no data.* reads).
 */
function expandClingenMoi(code) {
  const MAP = {
    AD: 'Autosomal Dominant', AR: 'Autosomal Recessive',
    XL: 'X-linked', XLD: 'X-linked Dominant', XLR: 'X-linked Recessive',
    YL: 'Y-linked', SD: 'Semidominant', MT: 'Mitochondrial', UD: 'Undetermined',
  };
  const key = String(code || '').trim().toUpperCase();
  return MAP[key] || code || '';
}

async function fetchGeneValidity(geneSymbol) {
  data.geneValidity = { curations: [], gene: geneSymbol || null, fetched: false, focusMondo: null };
  _clingenEvidenceProvenance = null;   // re-learned from the helper as this gene's evidence is fetched
  renderGeneValidityCard();   // clear stale content immediately on a new search
  if (!geneSymbol) return data.geneValidity;

  const csv = await _loadClingenValidityCsv();
  if (!csv) return data.geneValidity;   // proxy/network failure — stay empty

  const want = geneSymbol.toUpperCase();
  const curations = [];
  for (const line of csv.split('\n')) {
    if (!line.trim() || line.startsWith('"+++') ) continue;
    // Quick reject before the (more expensive) full parse
    if (!line.toUpperCase().includes(`"${want}"`)) continue;
    const f = _parseCsvLine(line);
    if ((f[0] || '').trim().toUpperCase() !== want) continue;
    curations.push({
      disease:        (f[2] || '').trim(),
      mondo:          (f[3] || '').trim(),
      moi:            (f[4] || '').trim(),
      sop:            (f[5] || '').trim(),
      classification: (f[6] || '').trim(),
      url:            (f[7] || '').trim(),
      date:           (f[8] || '').trim().split('T')[0],   // YYYY-MM-DD
      gcep:           (f[9] || '').trim()
    });
  }
  data.geneValidity = { curations, gene: geneSymbol, fetched: true, focusMondo: null };
  console.log(`[GeneValidity] ${geneSymbol}: ${curations.length} ClinGen curation(s)`);
  return data.geneValidity;
}

// The gene-validity Evidence Summary narrative is fetched in two tiers, both via the local Go helper
// (these sources have no CORS header, and the public CodeTabs proxy is blocked on office networks):
//   PRIMARY  — a structured `dc:description` field from GeneGraph's public JSON-LD bundle
//              (/api/clingen-evidence, keyed by the report URL's assertion UUID). Clean text,
//              ~85-90% coverage, with no dependency on the report-page HTML layout.
//   FALLBACK — scrape the ~330 KB per-curation HTML report page (/api/clingen-report → CodeTabs) for
//              curations the bundle doesn't carry, or when the helper is absent.
// Results are cached per report URL for the session.
const _clingenEvidenceCache = {};

// Provenance of the PRIMARY (GeneGraph) Evidence Summary data for THIS session, reported by the Go
// helper: 'live' = freshly downloaded; 'cache' = served from the helper's offline copy because it
// couldn't reach Google (e.g. IT block). Surfaced as an offline flag on the gene-validity card and
// in the Copy-To-Analysis-Summary report. Kept in a session-durable map (keyed by report URL) so a
// same-session re-search — which short-circuits on _clingenEvidenceCache before re-contacting the
// helper — still knows the source. clingenEvidenceProvenance() is the single read accessor (also
// used cross-file by buildReportHtml in app.js).
const _clingenEvidenceSrcCache = {};   // reportUrl -> { source:'live'|'cache', builtAt:'YYYY-MM-DD' }
let _clingenEvidenceProvenance = null; // latest known provenance for the current gene (reset per search)
function clingenEvidenceProvenance() { return _clingenEvidenceProvenance; }

/**
 * fetchGeneValidityEvidenceSummary(reportUrl)
 * Inputs: a ClinGen gene-validity curation report URL (explicit parameter).
 * Transport: PRIMARY /api/clingen-evidence (GeneGraph JSON-LD, structured). FALLBACK /api/clingen-report
 *   then the CodeTabs proxy (HTML scrape) — used when the bundle lacks this curation or the helper is absent.
 * Outputs: the GCEP "Evidence Summary" narrative as plain text; null when the curation genuinely has
 *   none; or undefined when ALL transports fail (so the caller can retry — failures are NOT cached).
 * Failures: Silent — logs a warning.
 */
async function fetchGeneValidityEvidenceSummary(reportUrl) {
  if (!reportUrl) return null;
  if (reportUrl in _clingenEvidenceCache) {
    // Same-session re-fetch: restore provenance so the offline flag survives (the text is cached but
    // we won't re-contact the helper to re-learn live vs cache).
    if (_clingenEvidenceSrcCache[reportUrl]) _clingenEvidenceProvenance = _clingenEvidenceSrcCache[reportUrl];
    return _clingenEvidenceCache[reportUrl];
  }

  // ── PRIMARY: structured Evidence Summary from the GeneGraph JSON-LD bundle (via the Go helper). ──
  // Keyed by the assertion UUID embedded in the report URL (CGGV:assertion_<uuid>-<ts>). A 204/empty
  // means the bundle has no narrative for this curation → fall through to the HTML scrape below.
  const uuidM = /assertion_([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})/i.exec(reportUrl);
  if (uuidM) {
    try {
      const res = await fetchWithTimeout(`/api/clingen-evidence?assertion=${uuidM[1]}`, {}, 30000);
      if (res.ok && res.status !== 204) {
        const j = await res.json();
        const raw = (j && typeof j.description === 'string') ? j.description : '';
        // Record helper-reported provenance (live vs offline copy) for the UI flag — independent of
        // whether THIS curation carries narrative text.
        if (j && (j.source === 'live' || j.source === 'cache')) {
          _clingenEvidenceProvenance = _clingenEvidenceSrcCache[reportUrl] = { source: j.source, builtAt: j.builtAt || '' };
        }
        // Strip GeneGraph's markdown gene-emphasis (*GENE*) and collapse whitespace for parity with the
        // scraped plain text.
        const clean = raw.replace(/\*([^*\n]+)\*/g, '$1').replace(/\s+/g, ' ').trim();
        if (clean) {
          _clingenEvidenceCache[reportUrl] = clean;
          return clean;
        }
      }
    } catch { /* helper/bundle unavailable — fall through to the HTML scrape */ }
  }

  // ── FALLBACK: scrape the per-curation HTML report page. ──────────────────────────────────────────
  // Rule 11: fetchWithTimeout for the fire-and-forget helper call; 30 s for the ~330 KB report page.
  const helperUrl = `/api/clingen-report?url=${encodeURIComponent(reportUrl)}`;
  const proxyUrl  = `https://api.codetabs.com/v1/proxy?quest=${encodeURIComponent(reportUrl)}`;

  let html = null;
  try {
    const res = await fetchWithTimeout(helperUrl, {}, 30000);
    if (res.ok) {
      const text = await res.text();
      // When the helper isn't served (dev without the binary → Vite proxy error or SPA shell) the
      // body won't carry ClinGen markers; treat that as a miss and fall through to CodeTabs.
      if (/Evidence Summary:|clinicalgenome/i.test(text)) html = text;
    }
  } catch { /* helper unavailable — fall through to CodeTabs */ }

  if (html === null) {
    try {
      const res = await fetchWithRetry(proxyUrl);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      html = await res.text();
    } catch (e) {
      console.warn('[GeneValidity] evidence summary fetch failed:', e.message);
      return undefined;   // genuine failure (both transports) — uncached so a later enrich() retries
    }
  }

  let summary = null;
  const i = html.indexOf('Evidence Summary:');
  if (i !== -1) {
    // Grab the <td colspan="…"> cell that follows the label cell
    const m = html.slice(i).match(/<td[^>]*colspan[^>]*>([\s\S]*?)<\/td>/i);
    if (m) {
      const ta = document.createElement('textarea');
      ta.innerHTML = m[1].replace(/<[^>]+>/g, ' ');   // strip tags
      summary = ta.value.replace(/\s+/g, ' ').trim();  // decode entities + collapse ws
    }
  }
  _clingenEvidenceCache[reportUrl] = summary || null;   // success (text) or fetched-but-empty (null)
  return _clingenEvidenceCache[reportUrl];
}

/**
 * enrichGeneValidityEvidence(curations)
 * Fills each curation's `evidenceSummary` field in parallel (in-place).
 * Non-blocking: a failed fetch (undefined) leaves the field undefined so a later call retries;
 * a real result (text, or null when the page has no summary) is stored and not re-fetched.
 */
async function enrichGeneValidityEvidence(curations) {
  if (!curations?.length) return;
  await Promise.all(curations.map(async c => {
    if (c.evidenceSummary === undefined && !c.evidenceSummaryError) {
      const r = await fetchGeneValidityEvidenceSummary(c.url);
      if (r !== undefined) c.evidenceSummary = r;
      else c.evidenceSummaryError = true;   // both transports failed — don't re-try until next search
    }
  }));
}

/**
 * onsetFromGeneReviews(curation, mim, chapters)  — PURE (inputs passed explicitly, Rule 2)
 * Fallback source for age of onset when OMIM has none: locate the GeneReviews chapter for the
 * curation's disease (by shared OMIM id, else fuzzy disease-name match) and pull the first
 * onset-bearing sentence from its clinical description. Returns a trimmed sentence or null.
 */
function onsetFromGeneReviews(curation, mim, chapters) {
  if (!Array.isArray(chapters) || !chapters.length) return null;
  let ch = mim ? chapters.find(c => (c.omim || []).map(String).includes(String(mim))) : null;
  if (!ch) ch = chapters.find(c => sameDiseaseName(curation.disease, c.title));
  if (!ch) return null;
  const text = String(ch.clinicalDescription || ch.genotypePhenotype || '').replace(/\s+/g, ' ').trim();
  if (!text) return null;

  // High precision: GeneReviews prose mixes feature-level timing and recovery durations with disease
  // age of onset, so a bare "onset" mention is too loose (e.g. "onset of aortic dilatation is variable"
  // or "resolution occurs within 15 days"). Accept a sentence only when an onset/presentation CUE
  // co-occurs with an explicit AGE-OF-LIFE token — that pairing is what denotes age of onset. CUE is
  // limited to genuine onset verbs (not "occurs/develops/diagnosed"); AGE is a life-stage or an
  // age-anchored number ("by age 6", "first 2 years of life"), never a bare duration ("within 15 days").
  const CUE = /\b(age (?:at|of) onset|onset|presents?|manifests?|begins?|first (?:symptoms?|signs?)|symptom onset)\b/i;
  const AGE = /\b(neonat\w*|congenital|prenatal\w*|antenatal\w*|in utero|at birth|infan\w*|childhood|adolesc\w*|adult\w*|(?:first|second|third|fourth|fifth|sixth|seventh|eighth|ninth|tenth)\s+decade|decades?|(?:by|at|around|before|after)\s+(?:age\s+)?\d+\s*(?:years?|months?|weeks?)|first\s+[\w-]+\s+(?:years?|months?)\s+of\s+life|\d+\s*(?:years?|months?)\s+of\s+age)\b/i;
  const sentences = text.split(/(?<=[.!?])\s+/).filter(s => CUE.test(s) && AGE.test(s));
  if (!sentences.length) return null;
  const pick = sentences.find(s => /\bonset\b/i.test(s)) || sentences[0];

  // Drop a leading discourse conjunction so the clause reads as a standalone onset statement.
  let out = pick.trim().replace(/^(however|but|in addition|additionally|moreover|furthermore|nevertheless|although|though)\s*,?\s*/i, '');
  out = out.charAt(0).toUpperCase() + out.slice(1);
  return out.length > 220 ? out.slice(0, 217).trim() + '…' : out;
}

/**
 * enrichGeneValidityOnset(curations, geneSymbol, geneReviewsChapters)
 * Fills each DEFINITIVE curation's `onset` + `onsetSource` (in-place). Per user request, age of
 * onset is surfaced only for ClinGen Definitive gene-disease validity entries.
 * Source priority: OMIM clinical synopsis (Miscellaneous field) → GeneReviews clinical description.
 * Inputs: curations + geneSymbol + geneReviewsChapters passed EXPLICITLY (Rule 2). Reads the OMIM
 *   module's own cache via getOmimOnset()/omimLocalCache (the established bridge, cf. resolveCurationToOmim).
 * Outputs: writes c.mimId / c.onset / c.onsetSource on the curation objects only.
 * Failures: non-blocking. c.onset === null means "looked up, none found" (won't refetch);
 *   c.onset === undefined means "not yet resolved" (UI shows "loading…").
 */
async function enrichGeneValidityOnset(curations, geneSymbol, geneReviewsChapters) {
  if (!curations?.length) return;
  if (typeof omimLocalReady === 'function') await omimLocalReady().catch(() => {});
  const geneDiseaseMims = (typeof omimLocalCache !== 'undefined'
    && omimLocalCache?.geneDiseaseMims?.get((geneSymbol || '').toUpperCase())) || [];
  await Promise.all(curations.map(async c => {
    if (!(c.classification || '').toLowerCase().includes('definitive')) return; // definitive only (per request)
    if (c.onset !== undefined && c.onsetSource) return;                          // already resolved with a source
    const mim = await resolveCurationToOmim(c, geneDiseaseMims).catch(() => null);
    if (mim) c.mimId = mim;
    let onset  = mim && typeof getOmimOnset === 'function' ? getOmimOnset(mim) : null;
    let source = onset ? 'OMIM' : null;
    if (!onset) {
      const gr = onsetFromGeneReviews(c, mim, geneReviewsChapters);
      if (gr) { onset = gr; source = 'GeneReviews'; }
    }
    c.onset = onset || null;     // null → looked up, none found (don't show "loading")
    c.onsetSource = source;      // null when nothing found in either source
  }));
}

// ── DISEASE IDENTITY LINKING ─────────────────────────────────────────────────
// Pure helpers (no data.* reads — CLAUDE.md Rule 2) that let the gene-disease
// modules (ClinGen validity, GenCC condition cards, GeneReviews, OMIM synopsis)
// reconcile the SAME disease across sources. Preference: exact ontology id, then
// fuzzy disease-name match. These are the join keys; passing them explicitly keeps
// each module from peeking into another's namespace.

// Significant (≥5-char) words of a disease label, lowercased, punctuation removed.
// Mirrors the matcher used by renderConditionCards (_sigWords) so name reconciliation
// behaves identically across the app.
function diseaseSigWords(name) {
  return String(name || '').toLowerCase().replace(/[^a-z0-9\s]/g, ' ')
    .split(/\s+/).filter(w => w.length >= 5);
}

// True when two disease labels denote the same disease: every significant word of
// the shorter name appears in the longer (handles "Noonan syndrome" ⊂ "Noonan syndrome 1").
function sameDiseaseName(a, b) {
  const wa = diseaseSigWords(a), wb = diseaseSigWords(b);
  if (!wa.length || !wb.length) return false;
  const [shorter, longer] = wa.length <= wb.length ? [wa, wb] : [wb, wa];
  return shorter.every(w => longer.includes(w));
}

// Normalized ontology id (MONDO/OMIM/ORPHA). Returns '' for anything else.
function normDiseaseId(id) {
  const s = String(id || '').trim().toUpperCase().replace(/\s+/g, '');
  return /^(MONDO|OMIM|OMIMPS|ORPHA|ORPHANET):/.test(s) ? s : '';
}

// True when two disease ids are the same recognised ontology id.
function sameDiseaseId(a, b) {
  const na = normDiseaseId(a);
  return !!na && na === normDiseaseId(b);
}

// Generic disease-label words that, on their own, must NOT establish identity. Two labels
// sharing only "syndrome" (e.g. "Marfan syndrome" vs "MASS syndrome") are NOT the same disease.
const _GENERIC_DZ_WORDS = new Set([
  'syndrome', 'disease', 'disorder', 'deficiency', 'susceptibility', 'familial', 'dominant',
  'recessive', 'autosomal', 'linked', 'related', 'congenital', 'infantile', 'juvenile',
  'complex', 'spectrum', 'combined', 'multiple', 'early', 'onset'
]);

// Pick the OMIM disease whose name best matches a curation label. Stricter than sameDiseaseName:
// a match must share ≥1 NON-generic significant word, so a generic-word-only overlap ("syndrome")
// can't win. Scores exact-name equality highest, then non-generic word overlap, then total overlap;
// avoids the "first loose hit wins" bug when a gene has many similarly-named OMIM phenotypes.
function bestOmimNameMatch(name, candidates) {
  const wq = diseaseSigWords(name);
  if (!wq.length || !candidates?.length) return null;
  const qNonGeneric = wq.filter(w => !_GENERIC_DZ_WORDS.has(w));
  let best = null, bestScore = 0;
  for (const d of candidates) {
    if (!d?.mimId) continue;
    const wd = diseaseSigWords(d.name);
    if (!wd.length) continue;
    const shared = wq.filter(w => wd.includes(w));
    const sharedNonGeneric = shared.filter(w => !_GENERIC_DZ_WORDS.has(w));
    if (!sharedNonGeneric.length) continue;                       // generic-only overlap → not a match
    let score = sharedNonGeneric.length * 10 + shared.length;
    if (wq.length === wd.length && wq.every(w => wd.includes(w))) score += 1000;   // exact set equality
    if (qNonGeneric.length && qNonGeneric.every(w => wd.includes(w))) score += 100; // curation ⊆ omim name
    if (score > bestScore) { bestScore = score; best = d; }
  }
  return best;
}

// Resolve a ClinGen validity curation to an OMIM mim number (string, e.g. "162200"),
// bridging into the OMIM-keyed phenotype machinery. Offline-first: match the curation's
// disease name against the gene's local OMIM diseases (passed explicitly), then fall back
// to the curation's MONDO id via Monarch xref (reusing resolveMondoToOmim). Cached by MONDO.
const _curationOmimCache = {};
async function resolveCurationToOmim(curation, geneDiseaseMims) {
  if (!curation) return null;
  const named = bestOmimNameMatch(curation.disease, geneDiseaseMims);
  if (named) return String(named.mimId);
  if (curation.mondo) {
    if (curation.mondo in _curationOmimCache) return _curationOmimCache[curation.mondo];
    const omim = await resolveMondoToOmim(curation.mondo);   // "OMIM:123456" | null
    const mim = omim ? omim.replace(/^OMIM:/i, '') : null;
    _curationOmimCache[curation.mondo] = mim;
    return mim;
  }
  return null;
}

async function fetchAssociatedConditions(geneSymbol) {
  if (!geneSymbol) return { conditions: [], aggregatedInheritance: '' };

  setAPIStatus('panelapp', 'loading');
  setAPIStatus('monarch', 'loading');

  const conditions = [];
  const sourceTracker = new Set();
  const moiTracker = new Set();

  // Predicates that indicate gene-disease causality
  const DISEASE_PREDICATES = [
    'biolink:causes', 'biolink:contributes_to', 'biolink:gene_associated_with_condition',
    'biolink:associated_with', 'biolink:has_phenotype', 'causally_related_to'
  ];

  // Disease ontology prefixes
  const DISEASE_PREFIXES = ['MONDO:', 'OMIM:', 'ORPHA:', 'OMIMPS:', 'Orphanet:'];

  const isDisease = (item) => {
    const obj = item.object || '';
    // Match by ID prefix
    if (DISEASE_PREFIXES.some(pfx => obj.startsWith(pfx))) return true;
    // Match by object_category if available
    const cat = (item.object_category || item.category || '').toLowerCase();
    if (cat.includes('disease') || cat.includes('phenotypic')) return true;
    // Match by predicate containing disease-related terms
    const pred = (item.predicate || '').toLowerCase();
    if (pred.includes('disease') || pred.includes('condition') || pred.includes('cause')) return true;
    return false;
  };

  // HPO MOI map
  const HPO_MOI_MAP = {
    'HP:0000006': 'Autosomal Dominant', 'HP:0000007': 'Autosomal Recessive',
    'HP:0001417': 'X-linked', 'HP:0001419': 'X-linked Recessive',
    'HP:0001423': 'X-linked Dominant', 'HP:0001450': 'Y-linked',
    'HP:0001428': 'Somatic Mutation', 'HP:0003745': 'Sporadic',
    'HP:0001426': 'Multifactorial', 'HP:0032113': 'Semidominant',
  };

  function parseMoiFromName(name) {
    if (!name) return null;
    const n = name.toUpperCase();
    // Semidominant (a.k.a. semi-dominant) — a distinct ClinGen/HPO MOI; check before the
    // dominant patterns so it isn't swallowed (it contains the substring "DOMINANT").
    if (n.includes('SEMIDOMINANT') || n.includes('SEMI-DOMINANT')) return 'Semidominant';
    // Pseudoautosomal before autosomal — "PSEUDOAUTOSOMAL DOMINANT" contains "AUTOSOMAL DOMINANT".
    if (n.includes('PSEUDOAUTOSOMAL DOMINANT')) return 'Pseudoautosomal Dominant';
    if (n.includes('PSEUDOAUTOSOMAL RECESSIVE')) return 'Pseudoautosomal Recessive';
    if (n.includes('AUTOSOMAL DOMINANT') || n.endsWith(', AD')) return 'Autosomal Dominant';
    if (n.includes('AUTOSOMAL RECESSIVE') || n.endsWith(', AR')) return 'Autosomal Recessive';
    if (n.includes('X-LINKED DOMINANT') || n.includes('X-LINKED, DOMINANT')) return 'X-linked Dominant';
    if (n.includes('X-LINKED RECESSIVE') || n.includes('X-LINKED, RECESSIVE')) return 'X-linked Recessive';
    if (n.includes('X-LINKED') || n.endsWith(', XL') || n.endsWith(', XLR')) return 'X-linked';
    if (n.includes('Y-LINKED')) return 'Y-linked';
    if (n.includes('MITOCHONDRIAL')) return 'Mitochondrial';
    if (n.includes('DIGENIC')) return 'Digenic';
    return null;
  }

  function extractMoiFromQualifiers(qualifiers) {
    if (!qualifiers) return null;
    const qList = Array.isArray(qualifiers) ? qualifiers : [qualifiers];
    for (const q of qList) {
      const qStr = typeof q === 'string' ? q : (q?.id || q?.label || '');
      if (HPO_MOI_MAP[qStr]) return HPO_MOI_MAP[qStr];
      const ql = qStr.toLowerCase();
      if (ql.includes('semidominant') || ql.includes('semi-dominant')) return 'Semidominant';
      // Pseudoautosomal before autosomal — "pseudoautosomal dominant" contains "autosomal dominant".
      if (ql.includes('pseudoautosomal dominant')) return 'Pseudoautosomal Dominant';
      if (ql.includes('pseudoautosomal recessive')) return 'Pseudoautosomal Recessive';
      if (ql.includes('autosomal dominant')) return 'Autosomal Dominant';
      if (ql.includes('autosomal recessive')) return 'Autosomal Recessive';
      if (ql.includes('x-linked recessive') || ql.includes('x linked recessive')) return 'X-linked Recessive';
      if (ql.includes('x-linked dominant') || ql.includes('x linked dominant')) return 'X-linked Dominant';
      if (ql.includes('x-linked') || ql.includes('x linked')) return 'X-linked';
      if (ql.includes('mitochondrial')) return 'Mitochondrial';
    }
    return null;
  }

  function standardizeMoi(rawStr) {
    if (!rawStr) return [];
    const s = typeof rawStr === 'string' ? rawStr.toUpperCase() : String(rawStr).toUpperCase();
    const modes = new Set();

    // Semidominant (a.k.a. semi-dominant) — distinct ClinGen/HPO MOI; would otherwise be dropped.
    if (s.includes('SEMIDOMINANT') || s.includes('SEMI-DOMINANT')) modes.add('Semidominant');
    // Pseudoautosomal MUST be tested before autosomal: "PSEUDOAUTOSOMAL DOMINANT" contains
    // the substring "AUTOSOMAL DOMINANT" (e.g. SHOX → Pseudoautosomal Dominant/Recessive).
    if (s.includes('PSEUDOAUTOSOMAL DOMINANT')) modes.add('Pseudoautosomal Dominant');
    else if (s.includes('AUTOSOMAL DOMINANT') || s.includes('MONOALLELIC')) modes.add('Autosomal Dominant');
    if (s.includes('PSEUDOAUTOSOMAL RECESSIVE')) modes.add('Pseudoautosomal Recessive');
    else if (s.includes('AUTOSOMAL RECESSIVE') || s.includes('BIALLELIC')) modes.add('Autosomal Recessive');
    if (s.includes('X-LINKED DOMINANT')) modes.add('X-linked Dominant');
    else if (s.includes('X-LINKED RECESSIVE')) modes.add('X-linked Recessive');
    else if (s.includes('X-LINKED') || s.includes(', XL')) modes.add('X-linked');
    if (s.includes('Y-LINKED')) modes.add('Y-linked');
    if (s.includes('MITOCHONDRIAL')) modes.add('Mitochondrial');
    if (s.includes('DIGENIC')) modes.add('Digenic');

    return Array.from(modes);
  }

  // ── 1. RUN PANELAPP (UK + AU) IN PARALLEL ───────────────────────────────
  const parsePanelAppResults = (results, sourceLabel, baseUrl) => {
    results.forEach(entry => {
      const panel = entry.panel || {};
      const panelName = panel.name || 'Unknown Panel';
      const confidence = String(entry.confidence_level || '');

      let evidenceLevel = 'Limited';
      if (confidence === '3') evidenceLevel = 'Strong';
      else if (confidence === '2') evidenceLevel = 'Moderate';
      else if (confidence === '1') evidenceLevel = 'Limited';

      const rawMoi = entry.mode_of_inheritance || '';
      const standardModes = standardizeMoi(rawMoi);

      let cleanMoi = '—';
      if (standardModes.length > 0) {
        standardModes.forEach(m => moiTracker.add(m));
        cleanMoi = standardModes.join(' / ');
      }
      sourceTracker.add(sourceLabel);

      conditions.push({
        diseaseName: panelName,
        diseaseId: panel.id ? `panel-${panel.id}` : '',
        externalUrl: panel.id ? `${baseUrl}/panels/${panel.id}/` : '',
        evidenceLevel,
        inheritance: cleanMoi,
        sources: `${sourceLabel} (confidence ${confidence}/3)`,
        isClinGen: false,
        _disorders: panel.relevant_disorders || []
      });
    });
  };

  const [paUkResult, paAuResult] = await Promise.allSettled([
    fetch(`https://panelapp.genomicsengland.co.uk/api/v1/genes/?entity_name=${encodeURIComponent(geneSymbol)}&format=json`),
    fetch(`https://panelapp-aus.org/api/v1/genes/?entity_name=${encodeURIComponent(geneSymbol)}&format=json`)
  ]);

  let paAnyOk = false;

  if (paUkResult.status === 'fulfilled' && paUkResult.value.ok) {
    try {
      const paData = await paUkResult.value.json();
      const paResults = paData?.results || [];
      console.log(`[PanelApp UK] Results: ${paResults.length}`);
      parsePanelAppResults(paResults, 'PanelApp UK', 'https://panelapp.genomicsengland.co.uk');
      if (paResults.length > 0) paAnyOk = true;
    } catch (e) { console.warn('[PanelApp UK] Parse failed:', e.message); }
  } else {
    console.warn('[PanelApp UK] Failed:', paUkResult.reason?.message || 'non-ok response');
  }

  if (paAuResult.status === 'fulfilled' && paAuResult.value.ok) {
    try {
      const paData = await paAuResult.value.json();
      const paResults = paData?.results || [];
      console.log(`[PanelApp AU] Results: ${paResults.length}`);
      parsePanelAppResults(paResults, 'PanelApp Australia', 'https://panelapp-aus.org');
      if (paResults.length > 0) paAnyOk = true;
    } catch (e) { console.warn('[PanelApp AU] Parse failed:', e.message); }
  } else {
    console.warn('[PanelApp AU] Failed:', paAuResult.reason?.message || 'non-ok response');
  }

  if (conditions.some(c => c.sources?.startsWith('PanelApp'))) {
    setAPIStatus('panelapp', 'ok');
  } else if (!paAnyOk) {
    setAPIStatus('panelapp', paUkResult.status === 'rejected' && paAuResult.status === 'rejected' ? 'error' : 'warn');
  }

  // ── 1b. LOCAL OMIM genemap2 — authoritative MOI per disease ─────────────
  // genemap2.txt Phenotypes column carries the canonical MOI string (e.g.
  // "Autosomal recessive") directly from OMIM. Use it before the Monarch
  // fallback so every confirmed OMIM disease contributes to moiTracker.
  await omimLocalReady();
  const localDiseases = omimLocalCache.geneDiseaseMims.get(geneSymbol.toUpperCase()) || [];
  for (const d of localDiseases) {
    if (d.evidenceCode >= 3 && d.moi) {
      standardizeMoi(d.moi).forEach(m => moiTracker.add(m));
    }
  }
  if (localDiseases.length > 0) sourceTracker.add('OMIM');

  // Generate a master MOI string from PanelApp + local OMIM (if any) to use as fallback for Monarch
  const masterFallbackMoi = Array.from(moiTracker).join(' / ');

  // ── 2. MONARCH (aggregates OMIM + ClinGen) ───────────────────────────────
  try {
    let hgncId = null;

    // Resolve HGNC ID via Ensembl
    try {
      const ensemblRes = await fetch(
        `https://rest.ensembl.org/lookup/symbol/homo_sapiens/${geneSymbol}?content-type=application/json`
      );
      if (ensemblRes.ok) {
        const ensemblData = await ensemblRes.json();
        const m = ensemblData.description?.match(/HGNC:(\d+)/);
        if (m) hgncId = `HGNC:${m[1]}`;
      }
    } catch (e) {
      console.warn('[Monarch] Ensembl lookup failed, trying Monarch search...');
    }

    // Fallback: Monarch search
    if (!hgncId) {
      const searchRes = await fetch(
        `https://api.monarchinitiative.org/v3/api/search?q=${encodeURIComponent(geneSymbol)}&category=biolink:Gene&limit=1`
      );
      if (searchRes.ok) {
        const searchData = await searchRes.json();
        if (searchData.items?.length > 0) hgncId = searchData.items[0].id;
      }
    }

    if (!hgncId) throw new Error(`Could not resolve HGNC ID for: ${geneSymbol}`);
    console.log(`[Monarch] ${geneSymbol} → ${hgncId}`);

    const monarchRes = await fetch(
      `https://api.monarchinitiative.org/v3/api/association?subject=${encodeURIComponent(hgncId)}&object_category=biolink:Disease&limit=500`
    );
    if (!monarchRes.ok) throw new Error(`HTTP ${monarchRes.status}`);
    const monarchData = await monarchRes.json();

    const allItems = monarchData.items || [];
    console.log(`[Monarch] Total items: ${allItems.length}`);

    // Filter to disease associations using broadened heuristic
    const diseaseItems = allItems.filter(item => isDisease(item));
    console.log(`[Monarch] Disease items: ${diseaseItems.length}`);

    diseaseItems.forEach(item => {
      const sourcesRaw = JSON.stringify({
        p: item.primary_knowledge_source,
        a: item.aggregator_knowledge_source,
        b: item.provided_by
      }).toLowerCase();

      const isClinGen = sourcesRaw.includes('clingen');
      const isOMIM = (item.object || '').startsWith('OMIM:') || sourcesRaw.includes('omim');
      const isOrphanet = (item.object || '').startsWith('ORPHA:') || sourcesRaw.includes('orphanet') || sourcesRaw.includes('orpha');

      // USER REQUEST: Remove Orphanet data
      if (isOrphanet && !isClinGen && !isOMIM) return;

      let sourceName = 'Monarch';
      if (isClinGen) sourceName = 'ClinGen';
      else if (isOMIM) sourceName = 'OMIM';
      else if (isOrphanet) sourceName = 'Orphanet'; // Left here safely in case of mixed origins

      sourceTracker.add(sourceName);

      const predicate = (item.predicate || '').toLowerCase();
      let evidenceLevel = 'Supportive';
      if (predicate.includes('cause')) evidenceLevel = 'Definitive';
      else if (isClinGen) evidenceLevel = 'Strong';
      else if (isOMIM) evidenceLevel = 'Moderate';
      else if (predicate.includes('contributes')) evidenceLevel = 'Limited';

      const diseaseId = item.object || '';
      const diseaseName = item.object_label || item.object_name || diseaseId || 'Unknown Condition';

      let rawMoi = extractMoiFromQualifiers(item.qualifiers)
        || extractMoiFromQualifiers(item.has_evidence)
        || parseMoiFromName(diseaseName)
        || '';

      let standardModes = standardizeMoi(rawMoi);

      // When Monarch carries no MOI for this disease, prefer THIS disease's authoritative OMIM
      // MOI (matched by name) over the gene-wide union: masterFallbackMoi is every source's modes
      // joined, so using it here stamps all of them onto each disease — e.g. a pseudoautosomal SHOX
      // disease would wrongly read "AD / AR / X-linked / PD / PR". Fall back to the union only when
      // no specific OMIM disease matches.
      if (standardModes.length === 0) {
        const omimMatch = localDiseases.find(d => d.evidenceCode >= 3 && d.moi && sameDiseaseName(diseaseName, d.name));
        if (omimMatch) standardModes = standardizeMoi(omimMatch.moi);
        else if (masterFallbackMoi) standardModes = masterFallbackMoi.split(' / ');
      }

      let moi = '—';
      if (standardModes.length > 0) {
        standardModes.forEach(m => moiTracker.add(m));
        moi = standardModes.join(' / ');
      }

      conditions.push({
        diseaseName,
        diseaseId,
        externalUrl: diseaseId ? `https://monarchinitiative.org/disease/${diseaseId}` : '',
        evidenceLevel,
        inheritance: moi,
        sources: sourceName,
        isClinGen
      });
    });

    if (conditions.length > 0) {
      setAPIStatus('monarch', 'ok');
    } else {
      setAPIStatus('monarch', 'warn');
    }

  } catch (err) {
    console.error('[Monarch] Error:', err);
    setAPIStatus('monarch', 'error');
  }

  // ── Result ───────────────────────────────────────────────────────────────
  const aggregatedInheritance = moiTracker.size > 0
    ? Array.from(moiTracker).join(' / ')
    : sourceTracker.size > 0 ? `Via: ${Array.from(sourceTracker).join(', ')}` : 'Unknown';

  console.log(`[Associations] Total: ${conditions.length} | MOI: ${aggregatedInheritance}`);

  return { conditions, aggregatedInheritance };
}



/**
 * Ensembl Primary Anchor: Used when VariantValidator is unavailable or bypassed.
 * Populates core coordinates and gene info from VEP first.
 */
async function runEnsemblPrimary(query) {
  setDot('vv', 'warn');
  const currentMsg = document.getElementById('statusMsg').innerText;
  if (currentMsg.includes('Unavailable')) {
    document.getElementById('statusMsg').innerHTML = '⚠️ <b style="color:var(--amber)">VV Timeout: Displaying Ensembl VEP Backup Data.</b>';
  } else {
    document.getElementById('statusMsg').innerText = 'Querying Ensembl VEP (Backup Mode)...';
  }

  try {
    let ensemblQuery = query;
    if (!ensemblQuery.includes(':')) {
      const parts = ensemblQuery.replace(/chr/i, '').split(/[- \s_]+/);
      if (parts.length >= 4) {
        ensemblQuery = `${parts[0]}:g.${parts[1]}${parts[2]}>${parts[3]}`;
      }
    } else {
      // NEW: Strip version numbers to prevent Ensembl 400 Bad Request errors
      ensemblQuery = ensemblQuery.replace(/\.\d+:/, ':');
    }

    const res = await fetchWithTimeout(`https://rest.ensembl.org/vep/human/hgvs/${encodeURIComponent(ensemblQuery)}?content-type=application/json&mane=1&numbers=1&hgvs=1&refseq=1`, {}, 20000);
    if (!res.ok) throw new Error(`VEP HTTP ${res.status}`);
    const json = await res.json();
    if (!json || !json.length) throw new Error('No VEP results');

    const main = json[0];
    data.coords.strand = main.strand || 1;
    const txList = main.transcript_consequences || [];
    txList.sort((a, b) => getSeverity(a.consequence_terms) - getSeverity(b.consequence_terms));
    let best = txList.find(t => t.mane && t.mane.includes('MANE_Select')) || txList.find(t => t.canonical === 1) || txList[0];

    // Populate core data
    data.coords.chrom = main.seq_region_name;
    const alleles = main.allele_string.split('/');
    let r = alleles[0];
    let a = alleles[1] || '';

    // Reverse complement alleles when Ensembl returns transcript-relative notation.
    // This happens when the input is transcript HGVS (NM_/NR_/ENST) on a minus-strand gene.
    // Genomic coordinate inputs already yield forward-strand alleles, so no RC needed there.
    if (main.strand === -1 && /^(NM_|NR_|NG_|ENST)/i.test(query)) {
      const comp = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-' };
      const reverseComp = seq => seq.split('').map(n => comp[n] || n).reverse().join('');
      r = reverseComp(r);
      a = reverseComp(a);
    }

    // [FIX] HGVS-style indel: Ensembl VEP uses '-' for pure insertions/deletions.
    // VEP main.start = insertion point / first deleted base (1-based). Anchor = main.start - 1.
    if (r === '-' || a === '-') {
      const vcfCoords = await hgvsIndelToVcf(
        data.coords.chrom,
        parseInt(main.start) - 1,   // 1-based anchor = one position before VEP start
        r === '-' ? '' : r,          // empty for insertions
        a === '-' ? '' : a           // empty for deletions
      );
      if (vcfCoords) {
        data.coords.pos38      = vcfCoords.pos;
        data.coords.ref        = vcfCoords.ref;
        data.coords.alt        = vcfCoords.alt;
        data.coords.hg38String = `chr${data.coords.chrom}-${vcfCoords.pos}-${vcfCoords.ref}-${vcfCoords.alt}`;
      } else {
        // Both Ensembl and NCBI anchor fetches failed — coords are incomplete.
        data.coords.pos38 = String(main.start);
        data.coords.ref   = normalizeVcfAllele(r);
        data.coords.alt   = normalizeVcfAllele(a);
        document.getElementById('statusMsg').innerHTML =
          '⚠️ <b style="color:var(--amber)">Indel anchor lookup failed (Ensembl + NCBI). Coordinates may be incomplete — downstream results unreliable.</b>';
      }
    } else {
      data.coords.pos38 = String(main.start);
      data.coords.ref   = normalizeVcfAllele(r);
      data.coords.alt   = normalizeVcfAllele(a);
      if (data.coords.ref && data.coords.alt) {
        data.coords.hg38String = `chr${data.coords.chrom}-${data.coords.pos38}-${data.coords.ref}-${data.coords.alt}`;
      }
    }

    if (best) {
      data.coords.gene = best.gene_symbol?.toUpperCase() || '-';
      data.ensembl.transcript = best.transcript_id || '-';

      // [FIX] Strip the transcript ID if a colon exists, matching VV's behavior
      data.coords.hgvs = best.hgvsc ? (best.hgvsc.includes(':') ? best.hgvsc.split(':')[1] : best.hgvsc) : '-';
      data.ensembl.protein = best.hgvsp ? (best.hgvsp.includes(':') ? best.hgvsp.split(':')[1] : best.hgvsp) : '-';
      if (best.protein_start) data.ensembl.vepProteinStart = best.protein_start;
      if (best.amino_acids) data.ensembl.vepAminoAcids = best.amino_acids;
      if (best.consequence_terms) data.ensembl.vepConsequence = best.consequence_terms.join(', ');
    }

    // Update UI
    document.getElementById('dGene').innerText = data.coords.gene;
    document.getElementById('geneInput').value = data.coords.gene;
    if (data.coords.hg38String) document.getElementById('dHg38').innerText = data.coords.hg38String;
    document.getElementById('eManeTranscript').innerText = data.ensembl.transcript;
    document.getElementById('dHgvs').innerText = data.coords.hgvs;
    document.getElementById('dProtein').innerText = data.ensembl.protein;

    // Await Liftover for hg19 (CRITICAL for downstream gnomAD v2)
    try {
      await fetchEnsemblMap(data.coords.chrom, data.coords.pos38, data.coords.ref, data.coords.alt);
    } catch (liftErr) {
      console.warn("Liftover failed during Ensembl primary:", liftErr);
    }

    // Cross-validate coords with ClinGen Allele Registry (strand-reliable ground truth).
    // Prefers the transcript HGVS Ensembl resolved; falls back to original query.
    const hgvsForAR = (data.ensembl.transcript && data.coords.hgvs && data.coords.hgvs !== '-')
      ? `${data.ensembl.transcript}:${data.coords.hgvs}` : query;
    await confirmCoordsWithClinGenAR(hgvsForAR);

    if (!currentMsg.includes('Unavailable')) {
      document.getElementById('statusMsg').innerText = 'Ready (Ensembl Primary).';
    }
    enableWrappers(['franklin', 'gnomadv4', 'spliceai', 'liftover', 'clingen', 'clinvar', 'omim', 'gr', 'scholar', 'decipher', 'hgmd', 'gtex', 'mastermind', 'gnomadv2']);
    document.getElementById('btnLaunchAll').disabled = false;


  } catch (e) {
    if (e.name === 'AbortError') return;
    console.error("Ensembl Primary Failed:", e);
    document.getElementById('statusMsg').innerText = `Ensembl Error: ${e.message}`;
    throw e;
  }
}


/**
 * CA ID Bridge: Resolves HGVS to Canonical Allele ID.
 * Try transcript HGVS first, then genomic as fallback.
 */
async function fetchCanonicalAllele(hgvs) {
  if (!hgvs) return;
  const attempts = [hgvs];
  // Also try without transcript version number (e.g. NM_000546.6 → NM_000546)
  const unversioned = hgvs.replace(/^(NM_\d+)\.\d+(:)/, '$1$2');
  if (unversioned !== hgvs) attempts.push(unversioned);
  if (data.coords.chrom && data.coords.pos38 && data.coords.ref && data.coords.alt) {
    attempts.push(`chr${data.coords.chrom}:g.${data.coords.pos38}${data.coords.ref}>${data.coords.alt}`);
  }
  for (const attempt of attempts) {
    try {
      const res = await fetch(
        `https://reg.clinicalgenome.org/AlleleRegistry/route/allele?hgvs=${encodeURIComponent(attempt)}`,
        { headers: { 'Accept': 'application/json' } }
      );
      if (res.ok) {
        const json = await res.json();
        if (json['@id']) {
          data.caId = json['@id'].split('=').pop();
          console.log('[ClinGen] Resolved CA ID:', data.caId, 'via:', attempt);
          return;
        }
      }
    } catch (e) {
      // 404 is expected for uncurated variants. Suppress explicit console.warn to reduce clutter.
      if (e.message !== 'Unexpected token \'N\', "Not Found" is not valid JSON') {
        // Only log actual execution errors, not just 404 parsing failures
      }
    }
  }
  data.caId = null;
}

/**
 * Expert VCEP Narratives: Shows a ClinGen audit link using the resolved CA ID.
 */
async function fetchVCEPNarrative() {
  // CA ID lookup is complete. Trigger the UI update.
  renderVCEPCard();
}

/**
 * CIViC Somatic Integration: Shows a direct CIViC search link badge.
 */
async function fetchCivicData() {
  if (!data.coords.gene || !data.ensembl.protein || data.ensembl.protein === '-') return;
  const varName = toOneLetterAA(data.ensembl.protein);
  const el = document.getElementById('civicBadge');
  if (!el) return;

  if (varName) {
    const civicUrl = `https://civicdb.org/variants/${encodeURIComponent(data.coords.gene)}/${encodeURIComponent(varName)}`;
    el.innerHTML = `<a href="${civicUrl}" target="_blank" style="color:white; text-decoration:none;">🧬 Search CIViC ↗</a>`;
    el.style.cssText = 'background:#673ab7; color:white; padding:3px 10px; border-radius:12px; font-size:0.7rem; font-weight:800; margin-left:12px; display:inline-block; vertical-align:middle; box-shadow: 0 0 10px rgba(103, 58, 183, 0.3);';
  } else {
    el.style.display = 'none';
  }
}

/**
 * MaveDB Functional Assay Integration via ClinGen LDH
 */
async function fetchMaveDBData() {
  if (!data.caId) return;

  setAPIStatus('mavedb', 'loading');
  try {
    const url = `https://ldh.clinicalgenome.org/ldh/Variant/caid/${data.caId}`;
    const res = await fetchWithRetry(url, { headers: { 'Accept': 'application/json' } });
    if (!res.ok) {
      setAPIStatus('mavedb', 'warn');
      data.maveData = null;
      renderMaveDBCard();
      return;
    }

    const json = await res.json();
    const payload = json.data || json.payload || json;

    // Defensive parsing for MaveDB provider
    let maveAssertion = null;

    // Standard LDH structure usually has 'contributions' or 'statements'
    const items = payload.contributions || payload.statements || payload.assertions || [];
    if (Array.isArray(items) && items.length > 0) {
      maveAssertion = items.find(item =>
        item?.agent?.name === 'MaveDB' ||
        item?.provider === 'MaveDB' ||
        item?.providedBy === 'MaveDB'
      );
    }

    if (!maveAssertion && Array.isArray(payload)) {
      maveAssertion = payload.find(item =>
        item?.agent?.name === 'MaveDB' ||
        item?.provider === 'MaveDB' ||
        item?.providedBy === 'MaveDB'
      );
    }

    if (maveAssertion) {
      // Extract score and interpretation safely
      const consequence = maveAssertion.realizes?.functional_consequence || maveAssertion.functional_consequence || maveAssertion.interpretation;
      const scoreObj = maveAssertion.realizes?.score || maveAssertion.score;

      const interpretation = consequence?.name || consequence || 'Unknown';
      const score = scoreObj?.value || scoreObj || 'N/A';

      if (interpretation !== 'Unknown' || score !== 'N/A') {
        data.maveData = {
          score: score,
          interpretation: interpretation
        };
        setAPIStatus('mavedb', 'ok');
      } else {
        console.warn('[MaveDB] Found provider block but could not extract score/interpretation structure.', maveAssertion);
        data.maveData = null;
        setAPIStatus('mavedb', 'warn');
      }
    } else {
      data.maveData = null;
      setAPIStatus('mavedb', 'warn');
    }
  } catch (e) {
    console.warn('[MaveDB] Fetch or parse failed:', e.message);
    data.maveData = null;
    setAPIStatus('mavedb', 'error');
  }

  renderMaveDBCard();
}

/**
 * Sanitizes HGVS.p protein string into 1-letter format (e.g. V600E).
 * Strips transcript prefix (e.g. NP_004324.2:) before parsing.
 */
function toOneLetterAA(p) {
  if (!p || p === '-') return '';
  const clean = p.includes(':') ? p.split(':').pop().replace(/✅/g, '').trim() : p.replace(/✅/g, '').trim();
  const m = clean.match(/p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)\)?/);
  if (!m) return '';
  const ref1 = AA_MAP_3TO1[m[1]] || '';
  const alt1 = (m[3] === 'Ter' || m[3] === '*') ? '*' : (AA_MAP_3TO1[m[3]] || '');
  return `${ref1}${m[2]}${alt1}`;
}

/**
 * Automated Start-Loss PVS1 Evaluation (ClinGen SVI Guidelines)
 *
 * Module isolation: writes ONLY to data.startLoss.* namespace.
 * Reads from data.coords.*, data.ensembl.* (legacy reads — see CLAUDE.md Rule #2 exemption).
 *
 * Algorithm (per Tayoun et al. SVI Working Group):
 *  1. Find next in-frame Met in MANE protein sequence (NCBI efetch).
 *  2. Search ClinVar for germline P/LP variants between codon 2 and next Met.
 *  3. If P/LP found → PVS1_Moderate. Otherwise → PVS1_Supporting.
 *  4. Generate ClinVar audit link with genomic range filter (germline P/LP only).
 *
 * @param {Object} [params] Optional explicit input override (preferred per CLAUDE.md Rule #2 for new code).
 *                          Falls back to data.* reads if not provided.
 */
async function evaluateStartLoss(params) {
  const gene       = params?.gene       ?? data.coords.gene;
  const transcript = params?.transcript ?? data.ensembl.transcript;
  const proteinId  = params?.proteinId  ?? data.proteinId;
  const hgvsC      = params?.hgvsC      ?? data.coords.hgvs;
  const hgvsP      = params?.hgvsP      ?? data.ensembl.protein;
  const pos38      = params?.pos38      ?? data.coords.pos38;
  const chrom      = params?.chrom      ?? data.coords.chrom;
  const cons       = params?.cons       ?? data.ensembl.vepConsequence ?? '';

  // Start-loss detection: VEP consequence OR HGVS-based fallback
  let isStartLoss = Array.isArray(cons) ? cons.includes('start_lost')
                                        : (typeof cons === 'string' && cons.includes('start_lost'));
  if (!isStartLoss && hgvsC && /^c\.[123][^0-9]/.test(hgvsC)) isStartLoss = true;
  if (!isStartLoss && hgvsP && /p\.\(?(Met1[^0-9]|M1[^0-9])/.test(hgvsP)) isStartLoss = true;

  if (!isStartLoss) return;
  if (!transcript || !transcript.startsWith('NM_')) return;

  const suggestEl = document.getElementById('mvSuggest');
  if (suggestEl) suggestEl.innerHTML = '<span style="color:var(--dim)">⏱ Evaluating SVI Start-Loss rescue...</span>';

  console.log('[StartLoss] Evaluating SVI Start-Loss rules for:', transcript);

  try {
    // 1. Fetch protein sequence (efetch with proxy fallback)
    const idToFetch = proteinId || transcript;
    const db = proteinId ? 'protein' : 'nuccore';
    const retType = proteinId ? 'fasta' : 'fasta_cds_aa';
    const fastaUrl = `${NCBI_BASE}/efetch.fcgi?db=${db}&id=${idToFetch}&rettype=${retType}&retmode=text`;

    let seqRes;
    try {
      seqRes = await fetchWithRetry(fastaUrl);
    } catch (e) {
      console.warn('[StartLoss] Direct NCBI fetch failed, trying proxy...', e);
      seqRes = await fetchWithRetry(`https://api.codetabs.com/v1/proxy?quest=${encodeURIComponent(fastaUrl)}`);
    }
    if (!seqRes.ok) throw new Error('Failed to retrieve protein sequence from NCBI');

    let sequence = seqRes && (await seqRes.text()).split('\n').slice(1).join('').trim();

    // Fallback to nuccore if proteinId path returned an invalid sequence
    if ((!sequence || !sequence.startsWith('M')) && proteinId) {
      const fallbackUrl = `${NCBI_BASE}/efetch.fcgi?db=nuccore&id=${transcript}&rettype=fasta_cds_aa&retmode=text`;
      const fbRes = await fetchWithRetry(fallbackUrl);
      if (fbRes.ok) sequence = (await fbRes.text()).split('\n').slice(1).join('').trim();
    }
    if (!sequence || !sequence.startsWith('M')) {
      if (suggestEl) suggestEl.innerHTML = '<span style="color:var(--dim)">⚠️ Could not resolve protein sequence for Start-Loss check.</span>';
      return;
    }

    // 2. Find next in-frame Met
    const nextMetIndex = sequence.indexOf('M', 1);
    const nextMetCodon = nextMetIndex === -1 ? null : nextMetIndex + 1;
    const nextMetDisplay = nextMetCodon ? `p.Met${nextMetCodon}` : 'Null';

    // 3. Map genomic range: codon 1 → next-Met codon via dummy variant
    const cdsEnd = (nextMetCodon || sequence.length) * 3;
    const dummyHgvs = `${transcript}:c.${cdsEnd}del`;
    let gStart = null, gEnd = null, gChrom = chrom;
    try {
      const dummyRes = await fetchWithRetry(
        `https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/${dummyHgvs}/mane_select`
      );
      const dummyJson = await dummyRes.json();
      const rec = Object.values(dummyJson)[0];
      const g38 = rec?.primary_assembly_loci?.grch38 || rec?.primary_assembly_loci?.hg38;

      if (g38?.vcf?.pos) {
        const hgvsMatch = (hgvsC || '').match(/c\.(\d+)/);
        const startPosInCds = hgvsMatch ? parseInt(hgvsMatch[1]) : 1;
        const offsetToStart = (startPosInCds - 1) % 3;
        const startCoord = parseInt(pos38) - offsetToStart;
        const endCoord = parseInt(g38.vcf.pos);
        gStart = Math.min(startCoord, endCoord);
        gEnd = Math.max(startCoord, endCoord);
        gChrom = g38.vcf.chr.replace('chr', '');
      }
    } catch (err) {
      console.warn('[StartLoss] Genomic range mapping failed:', err);
    }

    // 4. ClinVar search: germline P/LP variants within the genomic range
    let hasUpstreamPathogenic = false;
    if (gStart && gEnd) {
      const geneLower = gene.toLowerCase();
      // NCBI Entrez requires `N[Chromosome] AND start:end[Base Position for Assembly GRCh38]`
      const apiTerm = `"${geneLower}"[GENE] AND ${gChrom}[Chromosome] AND ${gStart}:${gEnd}[Base Position for Assembly GRCh38] AND ("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties]) AND "germline"[Origin]`;
      const searchUrl = `${NCBI_BASE}/esearch.fcgi?db=clinvar&term=${encodeURIComponent(apiTerm)}&retmode=json&retmax=500`;
      try {
        const cvRes = await fetchWithRetry(searchUrl);
        const cvJson = await ncbiSafeJson(cvRes);
        const uids = cvJson.esearchresult?.idlist || [];
        hasUpstreamPathogenic = uids.length > 0;
        console.log(`[StartLoss] Genomic-range P/LP variants found: ${uids.length}`);
      } catch (err) {
        console.warn('[StartLoss] ClinVar range search failed:', err);
      }
    }

    // 5. Build ClinVar audit URL — uses NCBI Entrez-compatible term so the link
    // returns the SAME variants as our API search (URL location/assembly params
    // alone are not enough; the term carries the filter).
    let clinVarUrl = '';
    if (gStart && gEnd) {
      const geneLower = gene.toLowerCase();
      const browserTerm = `"${geneLower}"[GENE] AND ${gChrom}[Chromosome] AND ${gStart}:${gEnd}[Base Position for Assembly GRCh38] AND ("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties]) AND "germline"[Origin]`;
      clinVarUrl = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(browserTerm)}`;
    } else {
      // Fallback: protein range if genomic mapping unavailable
      const protTerm = `(${gene}[Gene]) AND (${nextMetCodon || sequence.length}[Protein Location End] AND 1[Protein Location Start]) AND ("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties]) AND "germline"[Origin]`;
      clinVarUrl = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(protTerm)}`;
    }

    const pvs1Code  = hasUpstreamPathogenic ? 'PVS1_M' : 'PVS1_S';
    const pvs1Color = hasUpstreamPathogenic ? 'var(--amber)' : 'var(--teal)';
    const auditMsg  = gStart && gEnd
      ? `🔍 Audit ClinVar Range (p.1:${nextMetCodon || sequence.length} / ${gChrom}:${gStart}-${gEnd}) ↗`
      : `🔍 Audit ClinVar Protein Range (p.1:${nextMetCodon || sequence.length}) ↗`;

    // 6. Persist all results in isolated namespace + update UI
    data.startLoss.evaluated = true;
    data.startLoss.nextMetCodon = nextMetCodon;
    data.startLoss.nextMetDisplay = nextMetDisplay;
    data.startLoss.hasUpstreamPathogenic = hasUpstreamPathogenic;
    data.startLoss.pvs1Code = pvs1Code;
    data.startLoss.gChrom = gChrom;
    data.startLoss.gStart = gStart;
    data.startLoss.gEnd = gEnd;
    data.startLoss.clinVarUrl = clinVarUrl;

    suggestedCodes.add('PM2');
    suggestedCodes.add(pvs1Code);

    const startLossHtml = `
        <div style="display:flex; flex-wrap:wrap; align-items:center; gap:12px; font-size:0.75rem; background:rgba(255,255,255,0.03); padding:8px 12px; border-radius:8px; border:1px solid var(--border);">
          <div style="display:flex; align-items:center; gap:6px;">
            <span style="color:var(--dim)">Next In-Frame Met:</span>
            <span style="color:var(--text-bright); font-weight:700;">${nextMetDisplay}</span>
          </div>
          <div style="width:1px; height:12px; background:var(--border);"></div>
          <div style="display:flex; align-items:center; gap:6px;">
            <span style="color:var(--dim)">Suggested SVI Code:</span>
            <span style="background:${pvs1Color}22; color:${pvs1Color}; border:1px solid ${pvs1Color}; padding:1px 6px; border-radius:4px; font-weight:700; font-size:0.7rem;">${pvs1Code}</span>
          </div>
          <div style="width:1px; height:12px; background:var(--border);"></div>
          <a href="${clinVarUrl}" target="_blank" style="color:var(--teal); text-decoration:underline; font-weight:500;">${auditMsg}</a>
        </div>
    `;
    data.startLoss.html = startLossHtml;
    if (suggestEl) suggestEl.innerHTML = startLossHtml;

    evaluateACMG();

  } catch (e) {
    const errMsg = `<span style="color:var(--red)">⚠️ SVI Check Failed: ${e.message}</span>`;
    data.startLoss.error = e.message;
    data.startLoss.html = errMsg;
    if (suggestEl) suggestEl.innerHTML = errMsg;
  }
}

// ── SPLICE PS1 MODULE (ClinGen SVI 2024 Table 2) ───────────────────────────
//
// Owns all splice-variant logic — replaces the splice branch of fetchClinVarCodon.
// Module isolation: writes ONLY to data.splicePS1.* and the #splicePS1Block DOM card.
//
// Algorithm:
//  1. Detect splice variant (donor/acceptor) from VEP consequence or HGVS c.N±M
//  2. Map splice motif (per SVI 2024 Table 2) to a genomic range
//  3. ClinVar API: germline P/LP variants in that range
//  4. Fetch query SpliceAI + per-candidate SpliceAI (sequential, 1.5s gap)
//  5. Strict same-event detection: same delta peak position + same delta sign
//  6. Apply PS1 weight per Table 2 matrix
//  7. Render dedicated UI card

/**
 * resolveExonicSpliceAnchor(queryBase, gene, transcript, totalExons)
 * Resolves the splice motif anchor for an EXONIC variant (one with no intronic c.N±M
 * offset, e.g. a synonymous splice_region_variant at the exon edge). Inputs are passed
 * explicitly (Rule 2); fetches the exon cDNA boundaries via Ensembl.
 *
 * Motif windows (anchor = X):
 *   Donor    c.X−2 … c.X+6   → exonic coverage c.X−2…c.X (X = exon's LAST coding base)
 *   Acceptor c.X−20 … c.X    → exonic coverage c.X only  (X = exon's FIRST coding base)
 * The first transcript exon has no acceptor site; the last has no donor site.
 *
 * Returns:
 *   { spliceType, anchorBase, offsetFromAnchor }  variant lands inside a motif window
 *   { outsideMotif:true }   exon resolved, but variant is outside both exonic windows
 *   null                    boundaries unresolvable (fetch failed, or non-CDS coord like c.-N)
 */
async function resolveExonicSpliceAnchor(queryBase, gene, transcript, totalExons) {
  if (!Number.isFinite(queryBase)) return null;
  const exons = await fetchExonCodingPositions(transcript, gene);
  if (!exons.length) return null;
  const exon = exons.find(e => queryBase >= e.cdsStart && queryBase <= e.cdsEnd);
  if (!exon) return null; // queryBase not within any coding exon (e.g. 5'/3' UTR)

  const lastExonNum = totalExons || exons[exons.length - 1].number;
  const isFirstExon = exon.number === 1;
  const isLastExon  = exon.number === lastExonNum;

  // distance of the variant from each potential splice boundary (0 = the boundary base itself)
  const distToDonor    = isLastExon  ? Infinity : exon.cdsEnd   - queryBase; // toward 3' end of exon
  const distToAcceptor = isFirstExon ? Infinity : queryBase - exon.cdsStart; // toward 5' end of exon

  if (distToDonor >= 0 && distToDonor <= 2) {
    // c.X−2…c.X → anchor is the exon's last coding base; offset is 0/−1/−2 (exonic side)
    return { spliceType: 'donor', anchorBase: exon.cdsEnd, offsetFromAnchor: queryBase - exon.cdsEnd };
  }
  if (distToAcceptor === 0) {
    // Acceptor motif's only exonic base is c.X (the exon's first coding base)
    return { spliceType: 'acceptor', anchorBase: exon.cdsStart, offsetFromAnchor: 0 };
  }
  return { outsideMotif: true };
}

async function evaluateSplicePS1(params) {
  const gene       = params?.gene       ?? data.coords.gene;
  const transcript = params?.transcript ?? data.ensembl.transcript;
  const hgvsC      = params?.hgvsC      ?? data.coords.hgvs;
  const chrom      = params?.chrom      ?? data.coords.chrom;
  const pos38      = params?.pos38      ?? data.coords.pos38;
  const ref        = params?.ref        ?? data.coords.ref;
  const alt        = params?.alt        ?? data.coords.alt;
  const cons       = params?.cons       ?? data.ensembl.vepConsequence ?? '';

  if (!gene || !transcript || !hgvsC) return;

  // Detect splice type
  let spliceType = null;
  let anchorBase = 0;
  let offsetFromAnchor = 0;

  if (typeof cons === 'string' && cons.includes('splice_donor'))    spliceType = 'donor';
  else if (typeof cons === 'string' && cons.includes('splice_acceptor')) spliceType = 'acceptor';

  // Fallback / refine via HGVS c.N+M or c.N-M (simple) or c.N_M+K / c.N+K_M (range indel)
  // Priority 1: range with intronic END  →  c.1314_1315+4del  → anchor=1315, off=+4
  const rangeEndMatch  = hgvsC.match(/c\.\d+_(\d+)([+-])(\d+)/);
  // Priority 2: range with intronic START → c.1315+1_1316del  → anchor=1315, off=+1
  const rangeStartMatch = hgvsC.match(/c\.(\d+)([+-])(\d+)_\d+/);
  // Priority 3: simple intronic position  →  c.1315+1G>T
  const simpleMatch    = hgvsC.match(/c\.(\d+)([+-])(\d+)/);
  const cMatch = rangeEndMatch || rangeStartMatch || simpleMatch;
  if (cMatch) {
    anchorBase = parseInt(cMatch[1], 10);
    const sign = cMatch[2];
    offsetFromAnchor = (sign === '+' ? 1 : -1) * parseInt(cMatch[3], 10);
    if (!spliceType) spliceType = sign === '+' ? 'donor' : 'acceptor';
  } else {
    // Exonic variant (no intronic offset). Per the motif-membership rule, trigger the
    // audit only when the variant actually lands inside a donor/acceptor motif window.
    // The anchor is the exon's true cDNA boundary, so resolve it from the exon structure —
    // the variant's own c.N is NOT the boundary (it can sit 1–2 bp inside the exon).
    const eMatch = hgvsC.match(/^c\.(\d+)/);
    const isSpliceRegion = typeof cons === 'string' && /splice/.test(cons);
    if (!eMatch || (!spliceType && !isSpliceRegion)) return; // not a splice variant we handle
    const queryBase = parseInt(eMatch[1], 10);
    const resolved = await resolveExonicSpliceAnchor(queryBase, gene, transcript, data.ensembl.vepTotalExons);
    if (resolved && !resolved.outsideMotif) {
      spliceType       = resolved.spliceType;
      anchorBase       = resolved.anchorBase;
      offsetFromAnchor = resolved.offsetFromAnchor;
    } else if (resolved?.outsideMotif) {
      return; // exon resolved, but variant is outside both motif windows → do not trigger
    } else if (spliceType) {
      // Boundary lookup failed, but VEP already gave us the side → degrade gracefully to
      // the legacy approximation (variant's own c.N as anchor) so the card still surfaces.
      anchorBase = queryBase;
      offsetFromAnchor = 0;
    } else {
      return; // pure splice_region with no resolvable boundary → cannot place the motif
    }
  }
  if (!spliceType || !anchorBase) return;

  const card = document.getElementById('splicePS1Card');
  const content = document.getElementById('splicePS1Content');
  if (card) card.style.display = 'block';
  if (content) content.innerHTML = '<span style="color:var(--dim)">⏱ Computing splice motif & searching ClinVar...</span>';

  // Defensive: hide & clear the legacy cvAltBlock so its stale "Alternative Splice" content
  // never co-displays alongside the dedicated splice card.
  const altBlock = document.getElementById('cvAltBlock');
  const altLabel = document.getElementById('cvAltLabelText');
  const altCodon = document.getElementById('cvAltCodon');
  if (altBlock) altBlock.style.display = 'none';
  if (altLabel) altLabel.innerHTML = '';
  if (altCodon) altCodon.innerHTML = '';

  const debug = { startedAt: Date.now(), phases: {}, warnings: [] };
  data.splicePS1.debug = debug;

  try {
    // ── 1. Compute genomic motif range arithmetically (no VV calls needed) ──
    // Strand-aware mapping: genomic_pos(c.X + Δ) = anchorGenomic + Δ × strand
    // where anchorGenomic = pos38 − (offsetFromAnchor × strand)
    if (!pos38) throw new Error('Query variant has no hg38 position — cannot compute splice motif range');

    const queryPos = parseInt(pos38, 10);

    // cDNA deltas covering the splice motif (SVI Table 2):
    //   Donor:    c.(X-2) ... c.(X+6)  → 9 bp
    //   Acceptor: c.(X-20) ... c.(X)   → 21 bp
    const motifDeltas = [];
    if (spliceType === 'donor') {
      for (let d = -2; d <= 6; d++) motifDeltas.push(d);
    } else {
      for (let d = -20; d <= 0; d++) motifDeltas.push(d);
    }

    // Strand-aware mapping. If gene strand is unknown (VEP failed), try both
    // and pick the orientation where the query position falls INSIDE the motif.
    const computeRange = (strandSign) => {
      const anchorGenomic = queryPos - (offsetFromAnchor * strandSign);
      const positions = motifDeltas.map(d => anchorGenomic + d * strandSign);
      return { gStart: Math.min(...positions), gEnd: Math.max(...positions), anchorGenomic };
    };
    const candidateStrand = data.coords?.strand === -1 ? -1 : 1;
    let strandSign = candidateStrand;
    let range = computeRange(strandSign);
    // Sanity: query position must lie inside the motif. If not, flip strand.
    if (queryPos < range.gStart || queryPos > range.gEnd) {
      const flipped = computeRange(-strandSign);
      if (queryPos >= flipped.gStart && queryPos <= flipped.gEnd) {
        strandSign = -strandSign;
        range = flipped;
        console.log(`[SplicePS1] Strand inference: flipped to ${strandSign} based on query position`);
      }
    }
    const gChrom = String(chrom).replace('chr', '');
    const { gStart, gEnd } = range;

    const motifCStart = spliceType === 'donor' ? `c.${anchorBase - 2}` : `c.${anchorBase}-20`;
    const motifCEnd   = spliceType === 'donor' ? `c.${anchorBase}+6`  : `c.${anchorBase}`;

    data.splicePS1.spliceType = spliceType;
    data.splicePS1.anchorBase = anchorBase;
    data.splicePS1.motifRange = { cStart: motifCStart, cEnd: motifCEnd, gStart, gEnd, gChrom };
    console.log(`[SplicePS1] ${spliceType} motif: ${motifCStart} to ${motifCEnd} → chr${gChrom}:${gStart}-${gEnd} (strand ${strandSign})`);
    debug.phases.motif = { ms: Date.now() - debug.startedAt, gStart, gEnd, gChrom, strandSign };

    // ── 2. ClinVar search: broad query by gene + position range, filter P/LP client-side ──
    // We do NOT use the `"clinsig pathogenic"[Properties]` filter in Entrez because:
    //   (a) Entrez clinsig filter behavior shifted with the germline_classification rollout (2024+)
    //   (b) Some valid P/LP variants are missed depending on which classification field is set
    // Instead: fetch all variants in motif range, then filter for P/LP germline in esummary post-parse.
    // ClinVar Entrez: search by gene + splice consequence type.
    // Position-range fields ([CHRPOS], [Base Position for Assembly GRCh38]) are unreliable in Entrez.
    // Instead: fetch all splice_donor/acceptor variants for the gene, filter by motif cDNA range client-side.
    const geneLower = gene.toLowerCase();
    const spliceConsTerm = spliceType === 'donor'
      ? `"splice_donor_variant"[Molecular consequence]`
      : `"splice_acceptor_variant"[Molecular consequence]`;
    const apiTerm = `"${geneLower}"[GENE] AND ${spliceConsTerm}`;
    // ClinVar web UI supports [chrpos38] for genomic range; different from the API field.
    const clinVarUrl = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(`"${geneLower}"[gene] AND ${gChrom}[chr] AND ${gStart}:${gEnd}[chrpos38]`)}`;
    data.splicePS1.clinVarUrl = clinVarUrl;
    data.splicePS1.searchQuery = apiTerm;

    const searchUrl = `${NCBI_BASE}/esearch.fcgi?db=clinvar&term=${encodeURIComponent(apiTerm)}&retmode=json&retmax=300`;
    const sRes = await fetchWithRetry(searchUrl);
    const sJson = await ncbiSafeJson(sRes);
    const uids = sJson.esearchresult?.idlist || [];
    debug.phases.esearch = { ms: Date.now() - debug.startedAt, uidCount: uids.length, searchTerm: apiTerm };
    console.log(`[SplicePS1] esearch: ${sJson.esearchresult?.count} total, ${uids.length} returned | term: ${apiTerm}`);

    if (uids.length === 0) {
      data.splicePS1.evaluated = true;
      data.splicePS1.candidates = [];
      data.splicePS1.ps1Code = null;
      data.splicePS1.html = renderSplicePS1Card([], null, spliceType, clinVarUrl, null, gChrom, gStart, gEnd, apiTerm);
      if (content) content.innerHTML = data.splicePS1.html;
      return;
    }

    // ── 3. esummary → extract HGVS, classify, filter for splice + P/LP ──
    await new Promise(r => setTimeout(r, 300));
    const sumUrl = `${NCBI_BASE}/esummary.fcgi?db=clinvar&id=${uids.join(',')}&retmode=json`;
    const sumRes = await fetchWithRetry(sumUrl);
    const sumJson = await ncbiSafeJson(sumRes);

    const splicePattern = /(NM_\d+(?:\.\d+)?)(?:\([^)]+\))?:(c\.\d+(?:[+-]\d+)?(?:_\d+(?:[+-]\d+)?)?(?:[A-Z]>[A-Z]|del[A-Z]*|dup[A-Z]*|ins[A-Z]+|delins[A-Z]+))/;

    // Client-side motif range filter using cDNA positions.
    // Donor motif:   c.(anchorBase-2)  …  c.anchorBase+6
    // Acceptor motif: c.anchorBase-20  …  c.anchorBase
    // Handles simple (c.N+K) and range (c.N_M+K / c.N+K_M) HGVS notations.
    const checkIntronic = (vBase, sign, off) => {
      if (spliceType === 'donor')    return vBase === anchorBase && sign === '+' && off >= 1 && off <= 6;
      if (spliceType === 'acceptor') return vBase === anchorBase && sign === '-' && off >= 1 && off <= 20;
      return false;
    };
    const isInMotifRange = (vName) => {
      // Range with intronic end: c.N_M+K — use M+K as the position
      const rEnd = vName.match(/c\.\d+_(\d+)([+-])(\d+)/);
      if (rEnd) return checkIntronic(parseInt(rEnd[1], 10), rEnd[2], parseInt(rEnd[3], 10));
      // Range with intronic start: c.N+K_M — use N+K as the position
      const rStart = vName.match(/c\.(\d+)([+-])(\d+)_\d+/);
      if (rStart) return checkIntronic(parseInt(rStart[1], 10), rStart[2], parseInt(rStart[3], 10));
      // Simple intronic: c.N+K
      const intronic = vName.match(/c\.(\d+)([+-])(\d+)/);
      if (intronic) return checkIntronic(parseInt(intronic[1], 10), intronic[2], parseInt(intronic[3], 10));
      // Exonic (no offset)
      const exonic = vName.match(/^c\.(\d+)/);
      if (!exonic) return false;
      const vBase = parseInt(exonic[1], 10);
      if (spliceType === 'donor')    return vBase >= anchorBase - 2 && vBase <= anchorBase;
      if (spliceType === 'acceptor') return vBase === anchorBase;
      return false;
    };

    const allRecords = [];
    let plpCount = 0;
    let nonPlpCount = 0;
    for (const uid of uids) {
      const rec = sumJson.result?.[uid];
      if (!rec || !rec.title) continue;
      const m = rec.title.match(splicePattern);
      if (!m) continue;
      const hgvsT = `${m[1]}:${m[2]}`;
      if (!isInMotifRange(m[2])) continue;  // skip variants outside the motif window
      // Skip the query variant itself — comparing a variant to itself is not PS1 evidence
      if (m[2].toUpperCase() === (hgvsC || '').toUpperCase()) continue;
      const sig = rec.germline_classification?.description || rec.clinical_significance?.description || 'Unknown';
      const status = rec.germline_classification?.review_status || rec.clinical_significance?.review_status || '';
      // origin is a numeric bitmask in the ClinVar esummary API:
      //   1=germline  2=somatic  4=de_novo  8=maternal  16=paternal  32=inherited  64=familial  128=unknown
      // Any value without the somatic bit (2) set is treated as germline-compatible.
      const originRaw = rec.origin;
      const originNum = typeof originRaw === 'number' ? originRaw
        : typeof originRaw === 'string' && /^\d+$/.test(originRaw.trim()) ? parseInt(originRaw, 10)
        : null;
      const isPLP = /^pathogenic/i.test(sig) || /likely[\s_-]?pathogenic/i.test(sig);
      const isGermline = originNum === null
        ? (() => { const s = String(originRaw || '').toLowerCase(); return !s || s.includes('germline') || s.includes('unknown'); })()
        : (originNum & 2) === 0;  // somatic bit not set → germline-compatible
      if (isPLP && isGermline) plpCount++; else nonPlpCount++;
      allRecords.push({ uid, hgvsT, vName: m[2], classification: sig, status, accession: rec.accession_version || rec.accession, isPLP, isGermline });
    }
    debug.phases.esummary = { ms: Date.now() - debug.startedAt, totalParsed: allRecords.length, plpCount, nonPlpCount };
    // Debug: show in-motif records with their classifications
    console.log(`[SplicePS1] esummary: ${uids.length} UIDs → ${allRecords.length} in motif (${plpCount} P/LP, ${nonPlpCount} other) | in-motif records:`,
      allRecords.map(r => ({ hgvs: r.vName, sig: r.classification, isPLP: r.isPLP, isGermline: r.isGermline })));

    // Use P/LP germline variants for PS1 inference; fall back to showing all if none
    const candidates = allRecords.filter(r => r.isPLP && r.isGermline);
    data.splicePS1.allRecords = allRecords;

    if (candidates.length === 0) {
      data.splicePS1.evaluated = true;
      data.splicePS1.candidates = [];
      data.splicePS1.html = renderSplicePS1Card([], null, spliceType, clinVarUrl, null, gChrom, gStart, gEnd, apiTerm, allRecords);
      if (content) content.innerHTML = data.splicePS1.html;
      return;
    }

    // ── 4. SpliceAI for the QUERY variant (if not already cached) ──
    // Accept SNVs and short indels (≤50 bp on either allele) — the Broad SpliceAI-38
    // endpoint handles short indels. Skip larger variants to avoid API misuse.
    let querySAI = null;
    if (chrom && pos38 && ref && alt && ref.length <= 50 && alt.length <= 50) {
      querySAI = await fetchSpliceAIForVariant({ chr: chrom, pos: pos38, ref, alt }, data.coords.gene, data.ensembl.transcript);
    }
    data.splicePS1.querySpliceAI = querySAI;
    // Classify query event using SAI-10k-calc-style logic (used for canonical comparison below)
    const queryEvent = classifySpliceAIEvent(querySAI, parseInt(pos38, 10));
    data.splicePS1.querySpliceEvent = queryEvent;
    debug.phases.querySpliceAI = { ms: Date.now() - debug.startedAt, hasScores: !!querySAI, eventType: queryEvent.type };

    // ── 5. Resolve candidate genomic coords + SpliceAI sequentially (1.5s gap) ──
    if (content) content.innerHTML = `<span style="color:var(--dim)">⏱ Fetching SpliceAI for ${candidates.length} candidate(s)...</span>`;

    for (const cand of candidates) {
      try {
        const url = `https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/${encodeURIComponent(cand.hgvsT)}/mane_select`;
        const r = await fetchWithRetry(url);
        const j = await r.json();
        const rec = Object.values(j).find(v => v?.primary_assembly_loci);
        const vcf = rec?.primary_assembly_loci?.grch38?.vcf || rec?.primary_assembly_loci?.hg38?.vcf;
        if (vcf?.pos && vcf?.ref && vcf?.alt) {
          cand.coords = { chr: String(vcf.chr).replace('chr', ''), pos: vcf.pos, ref: String(vcf.ref).toUpperCase(), alt: String(vcf.alt).toUpperCase() };
          cand.spliceAI = await fetchSpliceAIForVariant(cand.coords, data.coords.gene, data.ensembl.transcript);
          // Strict DP equality (kept for backward compatibility — works only when query and
          // candidate sit at the same genomic position).
          cand.sameEvent = isSameSplicingEvent(querySAI, cand.spliceAI, spliceType);
          // SAI-10k-calc-inspired canonical comparison: classify each variant's splice event
          // (cryptic-donor-use, exon-skip, etc.) and require same type + identical absolute
          // coords at every named site. Replaces the naive single-axis absolute-coord check —
          // a variant that creates a cryptic donor (DL+DG hit) is NOT the same biological event
          // as one that causes exon skipping (DL+AL hit), even when their lost-donor positions
          // happen to align.
          cand.spliceEvent = classifySpliceAIEvent(cand.spliceAI, parseInt(cand.coords.pos, 10));
          cand.sameEventCanonical = isSameSpliceEventCanonical(queryEvent, cand.spliceEvent);
          const isAnyMatch = cand.sameEvent || cand.sameEventCanonical;
          cand.ps1Weight = isAnyMatch
            ? computePS1Weight(querySAI, cand.spliceAI, cand.classification, offsetFromAnchor, spliceType)
            : null;
        }
      } catch (err) {
        console.warn(`[SplicePS1] Failed candidate ${cand.hgvsT}:`, err.message);
        debug.warnings.push({ hgvs: cand.hgvsT, msg: err.message });
      }
      await new Promise(r => setTimeout(r, 1500));  // 1.5s gap per user spec
    }
    debug.phases.candidateLoop = { ms: Date.now() - debug.startedAt, withSpliceAI: candidates.filter(c => c.spliceAI).length };

    // Determine strongest PS1 weight across same-event candidates
    const sameEventCands = candidates.filter(c => (c.sameEvent || c.sameEventCanonical) && c.ps1Weight);
    const ps1Code = sameEventCands.reduce((best, c) => {
      const rank = { 'PS1': 4, 'PS1_Moderate': 3, 'PS1_Supporting': 2, 'N/A': 1 };
      return (rank[c.ps1Weight] || 0) > (rank[best] || 0) ? c.ps1Weight : best;
    }, null);

    data.splicePS1.querySpliceAI = querySAI;
    data.splicePS1.candidates = candidates;
    data.splicePS1.ps1Code = ps1Code;
    data.splicePS1.evaluated = true;

    if (ps1Code) {
      suggestedCodes.add(ps1Code);
      if (typeof evaluateACMG === 'function') evaluateACMG();
    }

    // ── 6. Render UI ──
    const queryHgvs = `${transcript}:${hgvsC}`;
    data.splicePS1.html = renderSplicePS1Card(candidates, querySAI, spliceType, clinVarUrl, ps1Code, gChrom, gStart, gEnd, data.splicePS1.searchQuery, data.splicePS1.allRecords, queryHgvs, queryEvent);
    if (content) content.innerHTML = data.splicePS1.html;

  } catch (e) {
    if (e.name === 'AbortError') return;
    console.error('[SplicePS1] error:', e);
    data.splicePS1.error = e.message;
    data.splicePS1.html = `<span style="color:var(--red)">⚠️ Splice PS1 Check Failed: ${e.message}</span>`;
    if (content) content.innerHTML = data.splicePS1.html;
  }
}

// Helper: fetch SpliceAI delta scores from Broad's API. 3-attempt exponential backoff (500/1000/2000 ms) on 429 / 5xx only.
async function fetchSpliceAIForVariant(coords, gene, transcript) {
  if (!coords?.chr || !coords?.pos || !coords?.ref || !coords?.alt) return null;
  const url = `https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&distance=500&mask=0&variant=${coords.chr}-${coords.pos}-${coords.ref}-${coords.alt}`;
  const tag = `${coords.chr}-${coords.pos}-${coords.ref}-${coords.alt}`;
  for (let attempt = 0; attempt < 3; attempt++) {
    try {
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), 30000);
      const res = await originalFetch(url, { signal: controller.signal });
      clearTimeout(timeoutId);
      if (res.status === 429 || res.status >= 500) {
        if (attempt < 2) {
          const delay = 500 * Math.pow(2, attempt);
          console.warn(`[SplicePS1 SpliceAI] ${res.status} for ${tag} — retrying in ${delay}ms (attempt ${attempt + 1}/3)`);
          await new Promise(r => setTimeout(r, delay));
          continue;
        }
        return null;
      }
      if (!res.ok) return null;
      const json = await res.json();
      if (json && json.scores && json.scores.length > 0) return pickSpliceAIScore(json.scores, gene, transcript);
      return null;
    } catch (e) {
      console.warn(`[SplicePS1 SpliceAI] failed for ${tag}:`, e.message);
      return null;
    }
  }
  return null;
}

// Strict same-event detection: same delta peak position + same delta sign.
// For donor variants: compare DS_DL & DP_DL (donor loss) — both must agree.
// For acceptor variants: compare DS_AL & DP_AL (acceptor loss) — both must agree.
// Cryptic gain events: same logic on DS_DG/DP_DG or DS_AG/DP_AG.
function isSameSplicingEvent(q, c, spliceType) {
  if (!q || !c) return false;
  const lossKey  = spliceType === 'donor' ? 'DL' : 'AL';
  const gainKey  = spliceType === 'donor' ? 'DG' : 'AG';
  const STRICT_THRESH = 0.2;  // SpliceAI delta must exceed this for a meaningful event

  const qLoss = parseFloat(q['DS_' + lossKey]) || 0;
  const cLoss = parseFloat(c['DS_' + lossKey]) || 0;
  const qLossPos = parseInt(q['DP_' + lossKey], 10);
  const cLossPos = parseInt(c['DP_' + lossKey], 10);

  if (qLoss >= STRICT_THRESH && cLoss >= STRICT_THRESH && qLossPos === cLossPos) {
    // Both predicted to disrupt the same splice site at the same position
    return true;
  }

  const qGain = parseFloat(q['DS_' + gainKey]) || 0;
  const cGain = parseFloat(c['DS_' + gainKey]) || 0;
  const qGainPos = parseInt(q['DP_' + gainKey], 10);
  const cGainPos = parseInt(c['DP_' + gainKey], 10);

  if (qGain >= STRICT_THRESH && cGain >= STRICT_THRESH && qGainPos === cGainPos) {
    // Both activate the same cryptic site
    return true;
  }

  return false;
}

// ── SAI-10k-calc-inspired splice event classification ──────────────────────
// Adapted from https://github.com/adavi4/SAI-10k-calc to enable biologically
// meaningful comparison between two variants' SpliceAI predictions.
//
// Rationale: comparing absolute genomic coordinates of a SINGLE axis (e.g. DL only)
// is too lenient — a variant that creates a cryptic donor (DL+DG hit) is a different
// molecular event from one that causes exon skipping (DL+AL hit), even when the
// canonical donor loss maps to the same position. The classifier folds the 4-axis
// SpliceAI Δ scores into a canonical event type, then comparison requires same type
// AND identical absolute coordinates at every named site of that event.
//
// Threshold: 0.2 (ClinGen SVI 2023 — matches the existing PS1 invariant; do not change).
function classifySpliceAIEvent(sai, variantPos) {
  if (!sai || !variantPos) return { type: 'no_data', sites: {}, hitAxes: [], absCoords: null };
  const T = 0.2;
  const ds = {
    AG: parseFloat(sai.DS_AG) || 0,
    AL: parseFloat(sai.DS_AL) || 0,
    DG: parseFloat(sai.DS_DG) || 0,
    DL: parseFloat(sai.DS_DL) || 0,
  };
  const abs = {
    AG: variantPos + (parseInt(sai.DP_AG, 10) || 0),
    AL: variantPos + (parseInt(sai.DP_AL, 10) || 0),
    DG: variantPos + (parseInt(sai.DP_DG, 10) || 0),
    DL: variantPos + (parseInt(sai.DP_DL, 10) || 0),
  };
  const hit = { AG: ds.AG >= T, AL: ds.AL >= T, DG: ds.DG >= T, DL: ds.DL >= T };
  const hitAxes = ['DL', 'AL', 'DG', 'AG'].filter(k => hit[k]);

  // Decision tree (priority order). Multi-axis (3-4 hits) compound events take
  // precedence over 2-axis events — per SAI-10k-calc, AL+AG+DL hits produce BOTH
  // cryptic acceptor activation AND exon skipping, which is a biologically
  // distinct event from a pure 2-axis cryptic acceptor use.
  // 4-axis compound
  if (hit.AL && hit.AG && hit.DL && hit.DG) {
    return { type: 'compound_acceptor_donor_switch', sites: { lostAcceptor: abs.AL, newAcceptor: abs.AG, lostDonor: abs.DL, newDonor: abs.DG }, hitAxes, absCoords: abs };
  }
  // 3-axis compound
  if (hit.AL && hit.AG && hit.DL) {
    return { type: 'cryptic_acceptor_with_exon_skip', sites: { lostAcceptor: abs.AL, newAcceptor: abs.AG, lostDonor: abs.DL }, hitAxes, absCoords: abs };
  }
  if (hit.DL && hit.DG && hit.AL) {
    return { type: 'cryptic_donor_with_exon_skip', sites: { lostDonor: abs.DL, newDonor: abs.DG, lostAcceptor: abs.AL }, hitAxes, absCoords: abs };
  }
  if (hit.AL && hit.AG && hit.DG) {
    return { type: 'cryptic_acceptor_with_donor_gain', sites: { lostAcceptor: abs.AL, newAcceptor: abs.AG, newDonor: abs.DG }, hitAxes, absCoords: abs };
  }
  if (hit.DL && hit.DG && hit.AG) {
    return { type: 'cryptic_donor_with_acceptor_gain', sites: { lostDonor: abs.DL, newDonor: abs.DG, newAcceptor: abs.AG }, hitAxes, absCoords: abs };
  }
  // 2-axis
  if (hit.DL && hit.DG) return { type: 'cryptic_donor_use',             sites: { lostDonor: abs.DL, newDonor: abs.DG }, hitAxes, absCoords: abs };
  if (hit.AL && hit.AG) return { type: 'cryptic_acceptor_use',          sites: { lostAcceptor: abs.AL, newAcceptor: abs.AG }, hitAxes, absCoords: abs };
  if (hit.DL && hit.AL) return { type: 'exon_skip_or_intron_retention', sites: { lostDonor: abs.DL, lostAcceptor: abs.AL }, hitAxes, absCoords: abs };
  if (hit.DG && hit.AG) return { type: 'pseudoexon_activation',         sites: { newDonor: abs.DG, newAcceptor: abs.AG }, hitAxes, absCoords: abs };
  // 1-axis
  if (hit.DL)           return { type: 'donor_loss_only',               sites: { lostDonor: abs.DL }, hitAxes, absCoords: abs };
  if (hit.AL)           return { type: 'acceptor_loss_only',            sites: { lostAcceptor: abs.AL }, hitAxes, absCoords: abs };
  if (hit.DG)           return { type: 'donor_gain_only',               sites: { newDonor: abs.DG }, hitAxes, absCoords: abs };
  if (hit.AG)           return { type: 'acceptor_gain_only',            sites: { newAcceptor: abs.AG }, hitAxes, absCoords: abs };
  return { type: 'no_significant_effect', sites: {}, hitAxes: [], absCoords: abs };
}

// Strict comparison aligned with SAI-10k-calc semantics:
//   (1) Same SET of hitting axes — no axis active in one variant but quiet in the other.
//       A candidate with AL+AG+DL (cryptic acceptor + exon skip) is NOT the same event
//       as a query with only AL+AG, even if the acceptor coords coincide.
//   (2) Identical absolute coordinates at every hitting axis (not just the "primary"
//       named sites). Catches subtle off-by-one cryptic-site differences.
function isSameSpliceEventCanonical(qEvent, cEvent) {
  if (!qEvent || !cEvent) return false;
  if (qEvent.type === 'no_data' || qEvent.type === 'no_significant_effect') return false;
  if (cEvent.type === 'no_data' || cEvent.type === 'no_significant_effect') return false;

  // (1) Hit-axis set must be identical
  const qAxes = (qEvent.hitAxes || []).slice().sort().join(',');
  const cAxes = (cEvent.hitAxes || []).slice().sort().join(',');
  if (qAxes !== cAxes) return false;

  // (2) Absolute coords must match at every hitting axis
  const qAbs = qEvent.absCoords || {};
  const cAbs = cEvent.absCoords || {};
  for (const axis of qEvent.hitAxes) {
    if (qAbs[axis] !== cAbs[axis]) return false;
  }
  return true;
}

function spliceEventLabel(type) {
  return ({
    'cryptic_donor_use':                'Cryptic donor activation',
    'cryptic_acceptor_use':             'Cryptic acceptor activation',
    'exon_skip_or_intron_retention':    'Exon skip / intron retention',
    'pseudoexon_activation':            'Pseudoexon activation',
    'donor_loss_only':                  'Donor loss only',
    'acceptor_loss_only':               'Acceptor loss only',
    'donor_gain_only':                  'Cryptic donor gain only',
    'acceptor_gain_only':               'Cryptic acceptor gain only',
    // Compound multi-axis events (per SAI-10k-calc — biologically distinct from 2-axis events)
    'cryptic_acceptor_with_exon_skip':  'Cryptic acceptor activation + exon skip',
    'cryptic_donor_with_exon_skip':     'Cryptic donor activation + exon skip',
    'cryptic_acceptor_with_donor_gain': 'Cryptic acceptor activation + donor gain',
    'cryptic_donor_with_acceptor_gain': 'Cryptic donor activation + acceptor gain',
    'compound_acceptor_donor_switch':   'Compound 4-axis splice switch',
    'no_significant_effect':            'No significant effect',
    'no_data':                          'No SpliceAI data',
  })[type] || type;
}

// PS1 weight per ClinGen SVI 2024 Table 2.
// VUA = variant under assessment (the query); P/LP = comparison variant from ClinVar.
// Strength of query prediction must be >= comparison (prerequisite).
function computePS1Weight(qSAI, cSAI, classification, offsetFromAnchor, spliceType) {
  if (!qSAI || !cSAI) return null;

  // Query prediction strength = max delta score across all 4 axes
  const maxScore = (s) => Math.max(
    parseFloat(s?.DS_AG) || 0,
    parseFloat(s?.DS_AL) || 0,
    parseFloat(s?.DS_DG) || 0,
    parseFloat(s?.DS_DL) || 0
  );
  const qStrength = maxScore(qSAI);
  const cStrength = maxScore(cSAI);
  if (qStrength < cStrength) return null;  // PS1 prerequisite fails

  const isP  = classification && /^pathogenic/i.test(classification);
  const isLP = classification && /likely[\s-]?pathogenic/i.test(classification);
  // The canonical ±1/±2 dinucleotide is INTRONIC: donor +1/+2, acceptor −1/−2. A bare
  // Math.abs() test would wrongly flag an EXONIC donor variant at c.X−1/c.X−2 (offset −1/−2)
  // as a dinucleotide hit, so the test must respect the splice side.
  const inPVS1Dinuc = spliceType === 'acceptor'
    ? (offsetFromAnchor === -1 || offsetFromAnchor === -2)
    : (offsetFromAnchor === 1 || offsetFromAnchor === 2);

  // Table 2 row 4-5: VUA at ±1/±2, comparison at ±1/±2 → PS1_S / N/A
  // Table 2 row 1-3: VUA outside ±1/±2 (or comparison outside) → PS1 / PS1_M / PS1_S
  if (inPVS1Dinuc) {
    return isP ? 'PS1_Supporting' : 'N/A';
  } else {
    if (isP)  return 'PS1_Moderate';
    if (isLP) return 'PS1_Supporting';
  }
  return null;
}

// Render the splice PS1 UI card content — 3-phase workflow view
// candidates: P/LP germline variants used for PS1 inference
// allRecords: every variant found in the motif range (P/LP + VUS + benign)
function renderSplicePS1Card(candidates, querySAI, spliceType, clinVarUrl, ps1Code, gChrom, gStart, gEnd, searchQuery, allRecords, queryHgvs, queryEvent) {
  candidates = candidates || [];
  allRecords = allRecords || [];
  const sameEventCount = candidates.filter(c => c.sameEvent).length;
  const queryPos = parseInt(data.coords.pos38, 10);

  // ── Phase 1: Motif retrieval ──
  const motifCDescription = spliceType === 'donor'
    ? `Donor motif (last 3 exonic + 6 intronic bp)`
    : `Acceptor motif (20 intronic + 1st exonic bp)`;
  const allCount = allRecords.length;
  const plpInAll = allRecords.filter(r => r.isPLP && r.isGermline).length;

  const phase1Html = `
    <div style="font-size:0.78rem; margin-bottom:12px; padding:10px; background:rgba(255,255,255,0.02); border-radius:6px; border-left:3px solid var(--amber);">
      <div style="font-weight:800; color:var(--text-bright); margin-bottom:6px; font-size:0.82rem;">▸ Phase 1: Search ClinVar in ${spliceType.toUpperCase()} motif</div>
      <div style="color:var(--dim); font-size:0.7rem; margin-bottom:4px;">${motifCDescription}</div>
      <div style="color:var(--dim); font-size:0.7rem; margin-bottom:4px;">Genomic range: <strong style="color:var(--text-bright); font-family:JetBrains Mono;">chr${gChrom}:${gStart}–${gEnd}</strong></div>
      <div style="color:var(--dim); font-size:0.7rem;">
        Found <strong style="color:var(--text-bright);">${allCount}</strong> total ClinVar variants in motif;
        <strong style="color:${plpInAll > 0 ? 'var(--amber)' : 'var(--dim)'};">${plpInAll}</strong> are germline P/LP (used for PS1)
      </div>
      <a href="${clinVarUrl}" target="_blank" style="color:var(--teal); text-decoration:underline; font-size:0.7rem; margin-top:4px; display:inline-block;">🔍 Verify in ClinVar ↗</a>
    </div>`;

  // Early empty state: no variants in motif at all
  if (allCount === 0) {
    return `
      ${phase1Html}
      <div style="padding:12px; text-align:center; color:var(--dim); font-size:0.75rem; background:rgba(255,255,255,0.02); border-radius:4px;">
        ℹ️ No ClinVar variants found in this motif range.<br>
        <span style="font-size:0.65rem;">Search query: <code style="background:rgba(0,0,0,0.3); padding:2px 6px; border-radius:3px;">${searchQuery || ''}</code></span>
      </div>`;
  }

  // Empty state when records exist but none are P/LP germline — show what was found anyway
  if (candidates.length === 0 && allCount > 0) {
    const recordRows = allRecords.slice(0, 20).map(r => {
      const sigClass = r.isPLP ? 'sig-p'
                     : /uncertain/i.test(r.classification) ? 'sig-vus'
                     : /benign/i.test(r.classification) ? 'sig-b'
                     : 'sig-vus';
      return `
        <tr style="border-bottom:1px solid rgba(255,255,255,0.05); font-size:0.66rem;">
          <td style="padding:3px 6px; color:var(--text-bright);">${r.hgvsT}</td>
          <td style="padding:3px 6px;"><span class="sig-badge ${sigClass}" style="padding:1px 4px; font-size:0.6rem;">${r.classification}</span></td>
        </tr>`;
    }).join('');
    return `
      ${phase1Html}
      <div style="padding:8px; background:rgba(255,255,255,0.02); border-radius:4px;">
        <div style="font-weight:700; color:var(--text-bright); font-size:0.75rem; margin-bottom:6px;">Variants in motif (no P/LP germline for PS1):</div>
        <table style="width:100%; font-size:0.66rem; border-collapse:collapse;">
          <thead><tr style="font-weight:700; color:var(--dim); border-bottom:1px solid rgba(255,255,255,0.1);">
            <th style="text-align:left; padding:3px 6px;">Variant</th>
            <th style="text-align:left; padding:3px 6px;">Classification</th>
          </tr></thead>
          <tbody>${recordRows}</tbody>
        </table>
      </div>`;
  }

  // ── Phase 2: SpliceAI comparison table ──
  const fmtSAI = (s) => {
    if (!s) return '<span style="color:var(--dim); font-size:0.62rem;">N/A</span>';
    const col  = (v) => parseFloat(v) >= 0.5 ? 'var(--red)' : parseFloat(v) >= 0.2 ? 'var(--amber)' : 'var(--dim)';
    const dp   = (v) => { const n = parseInt(v, 10); return isNaN(n) ? '' : ` <span style="color:var(--dim);">(${n >= 0 ? '+' : ''}${n})</span>`; };
    const sc   = (key) => `<span style="color:${col(s['DS_' + key])};">${key}: ${(parseFloat(s['DS_' + key]) || 0).toFixed(3)}${dp(s['DP_' + key])}</span>`;
    return `<span style="font-size:0.62rem; font-family:JetBrains Mono; display:grid; grid-template-columns:1fr 1fr; gap:1px 10px; line-height:1.6;">${sc('AG')}${sc('AL')}${sc('DG')}${sc('DL')}</span>`;
  };

  const phase2TableRows = candidates.map(c => {
    const sigClass = /^pathogenic/i.test(c.classification) ? 'sig-p'
                   : /likely[\s_-]?pathogenic/i.test(c.classification) ? 'sig-lp'
                   : 'sig-vus';
    const isMatch = c.ps1Weight && c.ps1Weight !== 'N/A' && (c.sameEvent || c.sameEventCanonical);
    const matchLabel = c.sameEvent ? 'SAME EVENT (strict DP)' : 'SAME EVENT (canonical)';
    const ps1Badge = isMatch
      ? `<span style="background:var(--teal)22; color:var(--teal); border:1px solid var(--teal); padding:1px 5px; border-radius:3px; font-weight:700; font-size:0.6rem;" title="${matchLabel}">PS1 → ${c.ps1Weight}</span>`
      : '<span style="color:var(--dim); font-size:0.6rem;">—</span>';
    const rowBg = isMatch ? 'background:rgba(46,229,157,0.05);' : '';
    const evType = c.spliceEvent?.type;
    const evLabel = evType ? spliceEventLabel(evType) : null;
    const evColor = isMatch ? 'var(--teal)' : 'var(--dim)';
    return `
      <tr style="border-bottom:1px solid rgba(255,255,255,0.05); font-size:0.66rem; ${rowBg}">
        <td style="padding:5px 6px; font-family:JetBrains Mono;">
          <span style="color:var(--text-bright); font-weight:600;">${c.hgvsT}</span>
          ${evLabel ? `<br><span style="color:${evColor}; font-size:0.6rem; font-style:italic;">⟶ ${evLabel}</span>` : ''}
        </td>
        <td style="padding:5px 6px;">${fmtSAI(c.spliceAI)}</td>
        <td style="padding:5px 6px; text-align:center;"><span class="sig-badge ${sigClass}" style="padding:1px 4px; font-size:0.6rem;">${c.classification.split(' ').slice(0,2).join(' ')}</span></td>
        <td style="padding:5px 6px; text-align:center;">${ps1Badge}</td>
      </tr>`;
  }).join('');

  const phase2Html = `
    <div style="font-size:0.78rem; margin-bottom:12px; padding:10px; background:rgba(255,255,255,0.02); border-radius:6px; border-left:3px solid var(--amber);">
      <div style="font-weight:800; color:var(--text-bright); margin-bottom:6px; font-size:0.82rem;">▸ Phase 2: SpliceAI Δ Scores (query vs. ${candidates.length} P/LP variant${candidates.length === 1 ? '' : 's'})</div>
      <table style="width:100%; font-size:0.65rem; border-collapse:collapse;">
        <thead>
          <tr style="border-bottom:2px solid rgba(255,255,255,0.1); font-weight:700; color:var(--dim);">
            <th style="text-align:left; padding:5px 6px;">Variant</th>
            <th style="text-align:left; padding:5px 6px;">SpliceAI Δ Scores (≥0.2 highlighted)</th>
            <th style="text-align:center; padding:5px 6px;">ClinVar</th>
            <th style="text-align:center; padding:5px 6px;">PS1 Inference</th>
          </tr>
        </thead>
        <tbody>
          <tr style="border-bottom:2px solid rgba(255,255,255,0.15); background:rgba(245,158,11,0.08); font-size:0.66rem;">
            <td style="padding:5px 6px; font-family:JetBrains Mono;">
              <span style="color:var(--amber); font-weight:800;">QUERY ←</span>
              ${queryHgvs ? `<br><span style="color:var(--text-bright); font-weight:600; font-size:0.62rem;">${queryHgvs}</span>` : ''}
              ${queryEvent?.type && queryEvent.type !== 'no_data' ? `<br><span style="color:var(--amber); font-size:0.6rem; font-style:italic;">⟶ ${spliceEventLabel(queryEvent.type)}</span>` : ''}
            </td>
            <td style="padding:5px 6px;">${fmtSAI(querySAI)}</td>
            <td colspan="2" style="padding:5px 6px; text-align:center; color:var(--dim); font-size:0.6rem;">(patient variant being audited)</td>
          </tr>
          ${phase2TableRows}
        </tbody>
      </table>
    </div>`;

  // ── Phase 3: Full SAI-10k-calc–style axis comparison ──
  // Show every hitting axis (DS ≥ 0.2) for both query and candidate with absolute
  // coordinates. PS1 requires the FULL hit-axis set to match — a candidate with an
  // extra DL hit (= exon skipping) is biologically distinct from a query with only
  // AL+AG hits, even when the acceptor sites coincide.
  const T = 0.2;
  const hasPhase3 = querySAI && (queryEvent?.hitAxes?.length > 0);
  const AXIS_LABELS = { AG: 'Acceptor Gain', AL: 'Acceptor Loss', DG: 'Donor Gain', DL: 'Donor Loss' };

  let phase3Html = '';
  if (hasPhase3) {
    const qHits = queryEvent.hitAxes || [];
    const qAbs  = queryEvent.absCoords || {};
    const qAxisSet = qHits.slice().sort().join(',');

    const phase3Rows = candidates
      .filter(c => c.spliceAI && c.coords && c.spliceEvent?.hitAxes?.length > 0)
      .map(c => {
        const cHits = c.spliceEvent.hitAxes || [];
        const cAbs  = c.spliceEvent.absCoords || {};
        const cAxisSet = cHits.slice().sort().join(',');
        const axisSetMatch = qAxisSet === cAxisSet;

        // Union of all hitting axes from either variant — that's what we display per row
        const allAxes = Array.from(new Set([...qHits, ...cHits]));
        const axisOrder = ['AL', 'AG', 'DL', 'DG'].filter(a => allAxes.includes(a));

        const axisRowsHtml = axisOrder.map(axis => {
          const inQ = qHits.includes(axis);
          const inC = cHits.includes(axis);
          const qCoord = inQ ? qAbs[axis] : null;
          const cCoord = inC ? cAbs[axis] : null;
          const coordMatch = inQ && inC && qCoord === cCoord;

          let cellColor = 'var(--dim)';
          let statusText = '';
          if (!inQ && !inC)      { statusText = '—'; }
          else if (inQ && !inC)  { cellColor = 'var(--red)';   statusText = '✗ Q-only'; }
          else if (!inQ && inC)  { cellColor = 'var(--red)';   statusText = '✗ C-only'; }
          else if (coordMatch)   { cellColor = 'var(--teal)';  statusText = '✓ match'; }
          else                   { cellColor = 'var(--amber)'; statusText = `✗ Δ${cCoord - qCoord}bp`; }

          const qCell = inQ ? `chr${gChrom}:${qCoord}` : '<span style="color:var(--dim);">—</span>';
          const cCell = inC ? `chr${gChrom}:${cCoord}` : '<span style="color:var(--dim);">—</span>';

          return `
            <tr style="font-size:0.6rem; border-bottom:1px solid rgba(255,255,255,0.04);">
              <td style="padding:2px 6px; color:var(--dim);">${AXIS_LABELS[axis]} (${axis})</td>
              <td style="padding:2px 6px; font-family:JetBrains Mono; color:var(--text-bright);">${qCell}</td>
              <td style="padding:2px 6px; font-family:JetBrains Mono; color:var(--text-bright);">${cCell}</td>
              <td style="padding:2px 6px; text-align:center; color:${cellColor}; font-weight:600;">${statusText}</td>
            </tr>`;
        }).join('');

        const verdictColor = axisSetMatch && c.sameEventCanonical ? 'var(--teal)' : 'var(--dim)';
        const verdictText = axisSetMatch && c.sameEventCanonical
          ? '✓ Same event'
          : !axisSetMatch
            ? `✗ Different axis pattern (Q={${qAxisSet}}, C={${cAxisSet}})`
            : '✗ Coord mismatch on hitting axis';

        return `
          <tr style="border-bottom:2px solid rgba(255,255,255,0.08); background:rgba(255,255,255,0.015);">
            <td colspan="4" style="padding:8px 6px 4px 6px;">
              <div style="font-family:JetBrains Mono; color:var(--text-bright); font-weight:700; font-size:0.7rem;">${c.hgvsT}</div>
              <div style="color:${c.sameEventCanonical ? 'var(--teal)' : 'var(--dim)'}; font-size:0.6rem; font-style:italic; margin-top:1px;">→ ${spliceEventLabel(c.spliceEvent.type)}</div>
            </td>
          </tr>
          <tr style="font-size:0.58rem; color:var(--dim); background:rgba(255,255,255,0.015);">
            <td style="padding:2px 6px; font-weight:700;">Axis</td>
            <td style="padding:2px 6px; font-weight:700;">Query coord</td>
            <td style="padding:2px 6px; font-weight:700;">Candidate coord</td>
            <td style="padding:2px 6px; font-weight:700; text-align:center;">Status</td>
          </tr>
          ${axisRowsHtml}
          <tr style="background:rgba(255,255,255,0.015);">
            <td colspan="4" style="padding:4px 6px 8px 6px; text-align:right; color:${verdictColor}; font-size:0.62rem; font-weight:700;">${verdictText}</td>
          </tr>`;
      }).join('');

    const queryEventLabel = spliceEventLabel(queryEvent.type);
    const queryAxesSummary = qHits.map(a => `<strong>${a}</strong> @ chr${gChrom}:${qAbs[a]}`).join(' · ');

    phase3Html = `
      <div style="font-size:0.78rem; margin-bottom:12px; padding:10px; background:rgba(255,255,255,0.02); border-radius:6px; border-left:3px solid var(--teal);">
        <div style="font-weight:800; color:var(--text-bright); margin-bottom:6px; font-size:0.82rem;">▸ Phase 3: SAI-10k-calc–style full event comparison (all hitting axes)</div>
        <div style="color:var(--dim); font-size:0.7rem; margin-bottom:8px;">
          <div><strong>Query event:</strong> <span style="color:var(--amber);">${queryEventLabel}</span></div>
          <div style="margin-top:2px;"><strong>Query hits (DS ≥ ${T}):</strong> <span style="font-family:JetBrains Mono;">${queryAxesSummary || '<em>none</em>'}</span></div>
        </div>
        <table style="width:100%; font-size:0.65rem; border-collapse:collapse;">
          <tbody>${phase3Rows}</tbody>
        </table>
        <div style="font-size:0.62rem; color:var(--dim); margin-top:8px; font-style:italic; line-height:1.5;">
          <strong>PS1 same-event rule (aligned with SAI-10k-calc):</strong> the candidate must hit
          exactly the same set of SpliceAI axes as the query (no extra or missing axis ≥ ${T}) AND
          the absolute genomic coordinate at every hitting axis must match. A candidate with an
          extra DL hit (e.g., exon-skipping signal) is a biologically distinct event from a query
          with only AL+AG hits, so PS1 should not apply.
        </div>
      </div>`;
  }

  // ── PS1 inference summary: WHICH variant supports PS1 ──
  const ps1Supporters = candidates.filter(c => (c.sameEvent || c.sameEventCanonical) && c.ps1Weight && c.ps1Weight !== 'N/A');
  let ps1SummaryHtml;
  if (ps1Code && ps1Supporters.length > 0) {
    const supporterRows = ps1Supporters.map(s => `
      <li style="margin-bottom:3px; font-size:0.72rem;">
        <span style="font-family:JetBrains Mono; color:var(--text-bright); font-weight:600;">${s.hgvsT}</span>
        <span style="color:var(--dim);">— ${s.classification}</span>
        <span style="background:var(--amber)22; color:var(--amber); border:1px solid var(--amber); padding:1px 5px; border-radius:3px; font-weight:700; font-size:0.65rem; margin-left:4px;">${s.ps1Weight}</span>
      </li>`).join('');
    ps1SummaryHtml = `
      <div style="padding:12px; background:rgba(245,158,11,0.08); border-radius:6px; border-left:3px solid var(--amber);">
        <div style="font-weight:800; color:var(--amber); margin-bottom:6px; font-size:0.85rem;">✓ PS1 Inference: <span style="font-size:0.95rem;">${ps1Code}</span></div>
        <div style="color:var(--dim); font-size:0.7rem; margin-bottom:6px;">PS1 can be inferred because the patient variant predicts the <strong>same splicing event</strong> as the following established P/LP variant${ps1Supporters.length > 1 ? 's' : ''}:</div>
        <ul style="margin:0; padding-left:20px; color:var(--text-bright);">${supporterRows}</ul>
      </div>`;
  } else {
    ps1SummaryHtml = `
      <div style="padding:12px; background:rgba(255,255,255,0.02); border-radius:6px; border-left:3px solid var(--dim);">
        <div style="font-weight:800; color:var(--dim); margin-bottom:4px; font-size:0.85rem;">✗ PS1 not applicable</div>
        <div style="color:var(--dim); font-size:0.7rem;">${candidates.length > 0 ? `${candidates.length} P/LP variant${candidates.length === 1 ? ' was' : 's were'} found in the motif, but ${candidates.length === 1 ? 'it does not' : 'none'} predict the same canonical splice event as the patient variant (SAI-10k-calc-style classification: same event type + identical absolute coordinates of affected splice sites, DS≥0.2).` : 'No P/LP variants found in motif to compare against.'}</div>
      </div>`;
  }

  return `
    ${phase1Html}
    ${phase2Html}
    ${phase3Html}
    ${ps1SummaryHtml}
  `;
}

function formatSpliceAIInline(s) {
  if (!s) return '';
  const fmt = (v, dp, key) => {
    const val = parseFloat(v) || 0;
    const color = val >= 0.5 ? 'var(--red)' : val >= 0.2 ? 'var(--amber)' : 'var(--dim)';
    return `<span style="color:${color}; font-weight:600;">${key}=${val.toFixed(2)}@${dp || 0}</span>`;
  };
  return [
    fmt(s.DS_AG, s.DP_AG, 'AG'),
    fmt(s.DS_AL, s.DP_AL, 'AL'),
    fmt(s.DS_DG, s.DP_DG, 'DG'),
    fmt(s.DS_DL, s.DP_DL, 'DL')
  ].join(' &nbsp; ');
}

// ── ClinVar Gene-Level Variant Distribution ────────────────────────────────
let _clinvarIndex = null;
let _clinvarIndexLoading = null;

async function loadClinVarIndex() {
  if (_clinvarIndex) return _clinvarIndex;
  // Prefer the pre-loaded global (works under file:// without a server)
  if (window._clinvarIndex) { _clinvarIndex = window._clinvarIndex; return _clinvarIndex; }
  if (_clinvarIndexLoading) return _clinvarIndexLoading;
  _clinvarIndexLoading = originalFetch('data/clinvar_gene_index.json')
    .then(r => r.ok ? r.json() : null)
    .then(j => { _clinvarIndex = j; _clinvarIndexLoading = null; return j; })
    .catch(() => { _clinvarIndexLoading = null; return null; });
  return _clinvarIndexLoading;
}

async function fetchClinVarGeneDistribution(gene) {
  if (!gene || gene === '-') return;
  const index = await loadClinVarIndex();
  if (!index) return;
  const geneData = index[gene] || index[gene.toUpperCase()] || null;
  const clinvarDate = index._meta?.date || null;
  renderVariantDistribution(gene, geneData, clinvarDate);

  // Update CV badge with gene-level distribution
  if (geneData) {
    updateAPIClassificationFromDistribution('cv', geneData);
  }
}

// ── GeneReviews — disease mechanism / penetrance / genotype-phenotype / clinical description ──
// LIVE: fetched on demand via the local Go helper's /api/genereviews route (genereviews.go),
// which reads NCBI's WEEKLY gene→chapter FTP mapping + the live Bookshelf chapter HTML and parses
// it server-side. The browser can't read Bookshelf directly (no CORS) and public proxies are
// IP-blocked, so the Go binary on the user's own desktop does the fetch. No bundled snapshot —
// if NCBI is unreachable the route returns no chapters and the card shows "—".
/**
 * fetchGeneReviews(gene, gen)
 * Inputs:  gene symbol (+ optional search-generation token for race-guarding). Does NOT read data.*.
 * Outputs: data.geneReviews.{gene, chapters, fetched}
 *   chapters: [{ nbk, ch_id, title, url, mechanism, penetrance, genotypePhenotype, clinicalDescription }]
 * Side effects: writes to #eGeneReviewsBody (Gene Pathogenicity card live UI).
 * Failures: Silent (non-blocking). Unreachable / gene not in mapping → chapters = []. UI shows "—".
 */
async function fetchGeneReviews(gene, gen) {
  data.geneReviews = { gene: gene || null, chapters: [], fetched: false, focusNbk: null, unavailable: false };
  if (!gene || gene === '-') { renderGeneReviewsCard(); return data.geneReviews; }
  let unavailable = false;
  try {
    const res = await fetchWithTimeout(`/api/genereviews?gene=${encodeURIComponent(gene)}`, {}, 30000);
    // /api/genereviews is served ONLY by the Go helper. Under `npm run dev` without the Vite /api
    // proxy (or with the helper down) the request resolves to Vite's index.html (200, text/html) —
    // detect that and surface "unavailable" rather than a misleading "no chapter found".
    const ctype = res.headers.get('content-type') || '';
    if (!res.ok || !ctype.includes('application/json')) {
      unavailable = true;
      throw new Error(res.ok ? `non-JSON response (${ctype || 'unknown content-type'})` : `HTTP ${res.status}`);
    }
    const j = await res.json();
    if (gen != null && gen !== window.currentSearchGen()) return data.geneReviews; // a newer search superseded this
    data.geneReviews = { gene, chapters: Array.isArray(j.chapters) ? j.chapters : [], fetched: true, focusNbk: null, unavailable: false };
    console.log(`[GeneReviews] ${gene}: ${data.geneReviews.chapters.length} chapter(s) [live]`);
  } catch (e) {
    if (gen != null && gen !== window.currentSearchGen()) return data.geneReviews;
    if (e.name === 'TypeError') unavailable = true; // network-level failure → backend not reachable
    // unavailable → fetched:false so a later copy/search retries (recovers if the helper starts later);
    // a valid empty JSON response keeps the normal "no chapter" path (fetched:true above).
    data.geneReviews = { gene, chapters: [], fetched: false, focusNbk: null, unavailable };
    if (unavailable) {
      console.warn('[GeneReviews] local helper not reachable — under `npm run dev` start the Go binary on :8770 and ensure the Vite /api proxy is set (see vite.config.js).');
    } else {
      console.warn('[GeneReviews] live fetch failed:', e.message);
    }
  }
  renderGeneReviewsCard();
  return data.geneReviews;
}

/**
 * focusGeneReviewsByPhenotype(omimResults, gen)
 * When the patient phenotype is entered, narrow the GeneReviews display to the ONE disease
 * chapter that best matches the phenotype, instead of listing all of the gene's chapters.
 * Inputs:  omimResults — ranked OMIM disease matches from matchOmimSynopsis() (passed EXPLICITLY;
 *            this module never reads data.ptPheno* / data.omimSynopsisResults itself — CLAUDE.md Rule 2).
 *          gen — optional search-generation token (race guard).
 * Outputs: sets data.geneReviews.focusNbk (own namespace only) and re-renders.
 * Bridge:  each chapter carries its OMIM ids (from NCBI); a chapter's score = the best fScore among
 *            its OMIM diseases. Picks the top-scoring chapter. Non-destructive — all chapters are kept;
 *            this only selects which to show. No usable match → focus stays null (show all).
 */
function focusGeneReviewsByPhenotype(omimResults, gen) {
  if (gen != null && gen !== window.currentSearchGen()) return;
  const gr = data.geneReviews;
  if (!gr?.fetched || !gr.chapters?.length || !Array.isArray(omimResults) || !omimResults.length) return;

  // mimId → best fScore among diseases the patient phenotype actually matched
  const scoreByMim = new Map();
  for (const d of omimResults) {
    if (!d?.mimId || !(d.matched > 0) || !(d.fScore > 0)) continue;
    const k = String(d.mimId);
    if (!scoreByMim.has(k) || d.fScore > scoreByMim.get(k)) scoreByMim.set(k, d.fScore);
  }
  if (!scoreByMim.size) return; // phenotype matched no disease for this gene → leave all shown

  let best = null, bestScore = 0;
  for (const c of gr.chapters) {
    let s = 0;
    for (const mim of (c.omim || [])) { const v = scoreByMim.get(String(mim)) || 0; if (v > s) s = v; }
    if (s > bestScore) { bestScore = s; best = c; }
  }
  if (best && bestScore > 0) {
    data.geneReviews.focusNbk = best.nbk;
    console.log(`[GeneReviews] phenotype focus → ${best.nbk} ${best.ch_id} (fScore ${bestScore.toFixed(2)})`);
    renderGeneReviewsCard();
  }
}

/**
 * focusGeneValidityByPhenotype(omimResults, gen)
 * Phenotype focus for the ClinGen validity card — the ClinGen analogue of
 * focusGeneReviewsByPhenotype. Highlights the ONE curated disease that best matches the
 * patient phenotype instead of listing all of the gene's curations flat.
 * Inputs:  omimResults — ranked OMIM matches from matchOmimSynopsis() (passed EXPLICITLY; this
 *            module never reads data.ptPheno* / data.omimSynopsisResults itself — CLAUDE.md Rule 2).
 *          gen — optional search-generation token (race guard, re-checked after awaits).
 * Bridge:  each curation's MONDO id → OMIM (resolveCurationToOmim: local OMIM-name match first,
 *            Monarch xref fallback); a curation's score = best fScore among its OMIM diseases.
 * Outputs: sets data.geneValidity.focusMondo (own namespace only) and re-renders. No usable
 *            match → focus stays null (all curations shown).
 */
async function focusGeneValidityByPhenotype(omimResults, gen) {
  if (gen != null && gen !== window.currentSearchGen()) return;
  const gv = data.geneValidity;
  if (!gv?.fetched || !gv.curations?.length || !Array.isArray(omimResults) || !omimResults.length) return;

  // mimId → best fScore among diseases the patient phenotype actually matched
  const scoreByMim = new Map();
  for (const d of omimResults) {
    if (!d?.mimId || !(d.matched > 0) || !(d.fScore > 0)) continue;
    const k = String(d.mimId);
    if (!scoreByMim.has(k) || d.fScore > scoreByMim.get(k)) scoreByMim.set(k, d.fScore);
  }
  if (!scoreByMim.size) return; // phenotype matched no disease for this gene → leave all shown

  await omimLocalReady().catch(() => {});
  const geneDiseaseMims = omimLocalCache?.geneDiseaseMims?.get((gv.gene || '').toUpperCase()) || [];

  let best = null, bestScore = 0;
  for (const c of gv.curations) {
    const mim = await resolveCurationToOmim(c, geneDiseaseMims);
    const s = mim ? (scoreByMim.get(String(mim)) || 0) : 0;
    if (s > bestScore) { bestScore = s; best = c; }
  }
  if (gen != null && gen !== window.currentSearchGen()) return; // a newer search superseded this
  if (best && bestScore > 0 && best.mondo) {
    data.geneValidity.focusMondo = best.mondo;
    console.log(`[GeneValidity] phenotype focus → ${best.mondo} ${best.disease} (fScore ${bestScore.toFixed(2)})`);
    renderGeneValidityCard();
  }
}

/**
 * renderGeneReviewsCard()
 * Renders data.geneReviews into the #eGeneReviewsBody element inside Gene Pathogenicity.
 * One block per chapter. Each block shows mechanism, penetrance, G-P (truncated with link),
 * and a <details> collapsible for Clinical Description.
 * Called by fetchGeneReviews after index lookup completes.
 */
function renderGeneReviewsCard() {
  const el = document.getElementById('eGeneReviewsBody');
  if (!el) return;
  if (data.geneReviews?.unavailable) {
    el.innerHTML = '<span style="color:var(--amber, #e0a030);" title="The GeneReviews route is served by the local Go helper. In dev (npm run dev) run the Go binary on :8770 and keep the Vite /api proxy.">GeneReviews unavailable — local helper not reachable</span>';
    return;
  }
  const allChapters = data.geneReviews?.chapters || [];
  if (!allChapters.length) {
    el.innerHTML = '<span style="color:var(--text-dim);">—</span>';
    return;
  }
  // Phenotype focus: when set, show only the best-matching disease chapter; the rest become links.
  const focusNbk = data.geneReviews?.focusNbk || null;
  const focused  = focusNbk ? allChapters.find(c => c.nbk === focusNbk) : null;
  const chapters = focused ? [focused] : allChapters;
  const others   = focused ? allChapters.filter(c => c.nbk !== focusNbk) : [];

  const esc = s => String(s || '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
  const truncate = (text, n) => {
    const t = String(text || '').trim();
    if (!t) return null;
    return t.length > n ? t.slice(0, n).replace(/\s+\S*$/, '') + '…' : t;
  };

  const html = chapters.map((c, i) => {
    const isFocus = !!focused && c.nbk === focusNbk;
    const badge = isFocus
      ? ` <span style="font-size:0.60rem;color:var(--bg,#0b0f14);background:var(--teal);border-radius:3px;padding:0 4px;font-weight:700;vertical-align:middle;">best phenotype match</span>`
      : '';
    const titleHtml = `<a href="${esc(c.url)}" target="_blank" rel="noopener"
        style="color:var(--teal);font-weight:600;font-size:0.78rem;text-decoration:none;"
        title="Open GeneReviews chapter">${esc(c.title || c.nbk || 'Chapter')}</a>${badge}`;

    const field = (label, val, maxLen) => {
      const t = truncate(val, maxLen);
      if (!t) return '';
      return `<div style="margin-top:3px;">
        <span style="color:var(--dim);font-size:0.70rem;font-weight:600;">${label}:</span>
        <span style="color:var(--text-bright);font-size:0.72rem;"> ${esc(t)}</span>
      </div>`;
    };

    const mechHtml   = field('Mechanism', c.mechanism, 120);
    const penHtml    = field('Penetrance', c.penetrance, 200);
    const gpHtml     = field('Geno-Pheno', c.genotypePhenotype, 200);

    // Clinical Description — collapsible <details> element
    const cdText = truncate(c.clinicalDescription, 500);
    const cdHtml = cdText ? `
      <details style="margin-top:4px;">
        <summary style="cursor:pointer;color:var(--dim);font-size:0.70rem;font-weight:600;list-style:none;user-select:none;">
          ▶ Clinical Description
        </summary>
        <div style="margin-top:3px;color:var(--text-bright);font-size:0.72rem;line-height:1.45;">${esc(cdText)}
          <a href="${esc(c.url)}" target="_blank" rel="noopener"
            style="color:var(--teal);font-size:0.68rem;margin-left:4px;">[full chapter]</a>
        </div>
      </details>` : '';

    const sep = i > 0 ? 'border-top:1px dashed rgba(255,255,255,0.08);margin-top:8px;padding-top:8px;' : '';
    return `<div style="${sep}">${titleHtml}${mechHtml}${penHtml}${gpHtml}${cdHtml}</div>`;
  }).join('');

  const othersHtml = others.length ? `
    <div style="margin-top:8px;border-top:1px dashed rgba(255,255,255,0.08);padding-top:6px;font-size:0.68rem;color:var(--dim);">
      Other GeneReviews chapters for this gene: ` +
      others.map(c => `<a href="${esc(c.url)}" target="_blank" rel="noopener" style="color:var(--teal);text-decoration:none;">${esc(c.title || c.nbk)}</a>`).join(', ') +
    `</div>` : '';

  el.innerHTML = html + othersHtml +
    `<div style="margin-top:6px;font-size:0.65rem;color:var(--text-dim);">
      © University of Washington · noncommercial/clinical use
    </div>`;
}

/**
 * extractDiseaseMechanism(evidenceSummary)  — PURE (no data.* reads, CLAUDE.md Rule 2)
 * Scans a ClinGen GCEP evidence-summary narrative for disease-mechanism statements
 * (loss-of-function, gain-of-function, dominant-negative) and tags each with a confidence
 * inferred from hedging language in the sentence where it appears. ClinGen states confirmed
 * mechanisms declaratively, so an un-hedged mention reads as 'established'; a hedged one
 * ('proposed', 'may', 'unclear', 'limited evidence', …) downgrades to 'proposed'.
 * Returns { LoF, GoF, DN } where each value is 'established' | 'proposed' | null, or null when the
 * text mentions no mechanism at all (caller renders "not stated").
 * Keyword regexes mirror evaluateMechanismConcordance(); kept local to stay pure & decoupled.
 * Heuristic limit: confidence is judged per sentence, so a negated mention sharing a sentence with
 * a positive one is not separately suppressed.
 */
function extractDiseaseMechanism(evidenceSummary) {
  const text = String(evidenceSummary || '').trim();
  if (!text) return null;

  const RX = {
    LoF: /loss[\s-]?of[\s-]?function|haploinsuff|nonsense[\s-]?mediated|null allele|tumou?r suppressor|\bLOF\b/i,
    GoF: /gain[\s-]?of[\s-]?function|activating|constitutive activ|increased (?:activity|function)|\bGOF\b/i,
    DN:  /dominant[\s-]?negative/i,
  };
  const HEDGE = /propos|suggest|hypothesi|postulat|possib|putative|candidate|unclear|uncertain|\bmay\b|\bmight\b|\bcould\b|not (?:yet |fully )?(?:establish|confirm|demonstrat)|remains? to be|insufficient|limited evidence|unknown/i;

  const result = { LoF: null, GoF: null, DN: null };
  const rank = { established: 2, proposed: 1 };
  for (const s of text.split(/(?<=[.!?])\s+/)) {
    const conf = HEDGE.test(s) ? 'proposed' : 'established';
    for (const cls of ['LoF', 'GoF', 'DN']) {
      if (!RX[cls].test(s)) continue;
      if (!result[cls] || rank[conf] > rank[result[cls]]) result[cls] = conf;
    }
  }
  return (result.LoF || result.GoF || result.DN) ? result : null;
}

/**
 * extractCurationVariants(evidenceSummary)  — PURE (no data.* reads, CLAUDE.md Rule 2)
 * Pulls the evidence-base sentence ClinGen GCEP summaries use, e.g.
 *   "At least 12 unique variants (missense, splice-site) ... in 9 publications ... in this curation."
 * The count may be a digit or spelled out ("Five variants"), with adjectives between it and
 * "variants" ("12 unique variants"). Returns { count, types[], pubCount } — count kept as the raw
 * token, pubCount null unless written as a digit — or null when no such sentence is present.
 */
function extractCurationVariants(evidenceSummary) {
  const text = String(evidenceSummary || '').trim();
  if (!text) return null;

  // "<N> [unique/distinct/…] variant(s) (type, type, …)" — the parenthetical is the anchor.
  // N is a digit OR spelled-out word ("Five variants", "Four unique variants"), and adjectives may
  // sit between the count and "variants" ("12 unique variants"). count is kept as the raw token.
  const NUM = '\\b(\\d+|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|nineteen|twenty)';
  const mv = text.match(new RegExp(NUM + '\\s+(?:[a-z-]+\\s+){0,2}variants?\\s*\\(([^)]+)\\)', 'i'));
  if (!mv) return null;
  const count = mv[1];
  const types = [...new Set(mv[2].split(/,|\band\b/i).map(t => t.trim().toLowerCase()).filter(Boolean))];
  if (!types.length) return null;

  const mp = text.match(/(\d+)\s+publications?/i);
  return {
    count,
    types,
    pubCount: mp ? Number(mp[1]) : null,
  };
}

/**
 * renderGeneValidityCard()
 * Renders data.geneValidity into #eGeneValidityBody — a dedicated live card in the Gene
 * Pathogenicity panel. Curations are grouped by disease (the same disease can have >1 curation);
 * every disease is listed. Per curation: classification pill + MOI/GCEP/date, a parsed disease
 * MECHANISM line (extractDiseaseMechanism) tagged established/proposed, and a VARIANTS EVALUATED
 * line (extractCurationVariants). A disease whose curations carry ≥2 distinct mechanism classes
 * gets a "refer to ClinGen for details" hint. When data.geneValidity.focusMondo is set (phenotype
 * focus, Phase 5) the matching disease sorts first with a "best phenotype match" badge.
 * Reads ONLY its own namespace (CLAUDE.md Rule 9).
 */
function renderGeneValidityCard() {
  const el = document.getElementById('eGeneValidityBody');
  if (!el) return;
  const curations = data.geneValidity?.curations || [];
  if (!curations.length) {
    el.innerHTML = data.geneValidity?.fetched
      ? '<span style="color:var(--text-dim);">No ClinGen gene-disease validity curation</span>'
      : '<span style="color:var(--text-dim);">—</span>';
    return;
  }

  const esc = s => String(s || '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  const COLORS = { definitive: '#1b5e20', strong: '#2e7d32', moderate: '#66bb6a', supportive: '#26a69a', limited: '#ff8a80', disputed: '#f44336', refuted: '#b71c1c', animal: '#ff9800' };
  const pill = cls => {
    const key = Object.keys(COLORS).find(k => String(cls || '').toLowerCase().includes(k));
    const bg = COLORS[key] || '#90a4ae';
    return `<span style="background:${bg};color:#fff;border-radius:8px;padding:1px 7px;font-size:0.62rem;font-weight:800;white-space:nowrap;">${esc(cls || '—')}</span>`;
  };
  const MECH_LABEL = { LoF: 'LOF', GoF: 'GOF', DN: 'Dominant-negative' };
  const CONF_COLOR = { established: '#66bb6a', proposed: '#ffb74d' };

  const mechLine = c => {
    if (c.evidenceSummary === undefined && !c.evidenceSummaryError)
      return `<div style="color:var(--text-dim);font-size:0.66rem;margin-top:2px;">Loading mechanism…</div>`;
    const mech = extractDiseaseMechanism(c.evidenceSummary);
    if (!mech)
      return `<div style="color:var(--text-dim);font-size:0.66rem;margin-top:2px;">Mechanism: not stated</div>`;
    const tags = ['LoF', 'GoF', 'DN'].filter(k => mech[k]).map(k =>
      `<b style="color:var(--text-bright);">${MECH_LABEL[k]}</b> <span style="color:${CONF_COLOR[mech[k]]};font-size:0.62rem;">(${mech[k]})</span>`
    ).join(' · ');
    return `<div style="color:var(--dim);font-size:0.68rem;margin-top:2px;">Mechanism: ${tags}</div>`;
  };

  const variantLine = c => {
    const v = extractCurationVariants(c.evidenceSummary);
    if (!v) return '';
    const pubs = v.pubCount ? ` · ${v.pubCount} publication${v.pubCount === 1 ? '' : 's'}` : '';
    return `<div style="color:var(--text-dim);font-size:0.66rem;margin-top:1px;">Variants evaluated: ${v.count} (${v.types.map(esc).join(', ')})${pubs}</div>`;
  };

  // Phenotype focus (Phase 5): the matching curation sorts first so its disease group leads.
  const focusMondo = data.geneValidity?.focusMondo || null;
  const ordered = (focusMondo
    ? [...curations].sort((a, b) => (b.mondo === focusMondo ? 1 : 0) - (a.mondo === focusMondo ? 1 : 0))
    : curations
  ).filter(c => (c.classification || '').toLowerCase().includes('definitive'));

  // Group by disease — the SAME disease can have >1 curation (different MOI / re-curation / SOP).
  // Every disease is listed (incl. ones whose mechanism can't be parsed). Key prefers stable MONDO.
  const groups = [];
  const byKey = new Map();
  for (const c of ordered) {
    const key = c.mondo || (c.disease || '').toLowerCase();
    let g = byKey.get(key);
    if (!g) { g = { key, disease: c.disease, curations: [] }; byKey.set(key, g); groups.push(g); }
    g.curations.push(c);
  }

  const blocks = groups.map((g, gi) => {
    const sep = gi > 0 ? 'border-top:1px dashed rgba(255,255,255,0.08);margin-top:6px;padding-top:6px;' : '';
    const multi = g.curations.length > 1;
    const countBadge = multi
      ? ` <span style="font-size:0.58rem;color:var(--dim);background:rgba(255,255,255,0.08);border-radius:3px;padding:0 4px;font-weight:700;vertical-align:middle;">${g.curations.length} curations</span>`
      : '';
    const focusBadge = (focusMondo && g.key === focusMondo)
      ? ` <span style="font-size:0.58rem;color:var(--bg,#0b0f14);background:var(--teal);border-radius:3px;padding:0 4px;font-weight:700;vertical-align:middle;">best phenotype match</span>`
      : '';
    // Age of onset (disease-level; same MIM across a disease's curations). OMIM synopsis primary,
    // GeneReviews fallback — resolved by enrichGeneValidityOnset. Definitive diseases only.
    const onsetCur  = g.curations.find(c => c.onset);
    const onsetLine = onsetCur
      ? `<div style="color:var(--dim);font-size:0.68rem;margin-top:2px;">Age of onset: <span style="color:var(--text-bright);">${esc(onsetCur.onset)}</span> <span style="color:var(--text-dim);font-size:0.6rem;">(${esc(onsetCur.onsetSource || 'OMIM')})</span></div>`
      : g.curations.some(c => c.onset === undefined)
        ? `<div style="color:var(--text-dim);font-size:0.66rem;margin-top:2px;">Age of onset: loading…</div>`
        : `<div style="color:var(--text-dim);font-size:0.66rem;margin-top:2px;">Age of onset: not stated</div>`;

    // Single curation: link the disease title to its report. Multiple: plain title + per-row links.
    const titleUrl = !multi && g.curations[0].url ? g.curations[0].url : '';
    const title = titleUrl
      ? `<a href="${esc(titleUrl)}" target="_blank" rel="noopener" style="color:var(--text-bright);font-weight:600;text-decoration:none;" title="Open ClinGen evidence report">${esc(g.disease || '—')}</a>`
      : `<span style="color:var(--text-bright);font-weight:600;">${esc(g.disease || '—')}</span>`;

    const rows = g.curations.map((c, ci) => {
      const link = multi && c.url
        ? ` <a href="${esc(c.url)}" target="_blank" rel="noopener" style="color:var(--teal);font-size:0.62rem;text-decoration:none;" title="Open ClinGen evidence report">report ↗</a>`
        : '';
      const rowSep = ci > 0 ? 'margin-top:4px;padding-top:4px;border-top:1px dotted rgba(255,255,255,0.06);' : '';
      return `<div style="${rowSep}${multi ? 'margin-left:6px;' : ''}">
        ${link ? `<div>${link}</div>` : ''}
        ${mechLine(c)}
        ${variantLine(c)}
      </div>`;
    }).join('');

    // Multiple-mechanism hint: distinct mechanism classes for this disease — counted across the
    // disease's curations AND within any single evidence summary (so a lone curation whose summary
    // states e.g. both loss- and gain-of-function also triggers it). ≥2 classes → refer to ClinGen.
    const classes = new Set();
    for (const c of g.curations) {
      const m = extractDiseaseMechanism(c.evidenceSummary);
      if (m) ['LoF', 'GoF', 'DN'].forEach(k => { if (m[k]) classes.add(k); });
    }
    const hint = classes.size >= 2
      ? `<div style="font-size:0.62rem;margin-top:2px;color:var(--text-dim);font-style:italic;">Multiple disease mechanisms reported — refer to ClinGen for more details.</div>`
      : '';

    return `<div style="${sep}">
      <div style="display:flex;align-items:center;gap:6px;flex-wrap:wrap;">${title}${countBadge}${focusBadge}</div>
      ${onsetLine}
      ${rows}
      ${hint}
    </div>`;
  }).join('');

  // Offline flag: the Evidence Summary text came from the Go helper's on-disk copy (couldn't reach
  // Google to refresh). Shown only when stale; live data carries no badge.
  const prov = clingenEvidenceProvenance();
  const offlineBadge = (prov && prov.source === 'cache')
    ? `<div style="margin-bottom:5px;font-size:0.62rem;font-weight:700;color:#0b0f14;background:#ffb74d;border-radius:4px;padding:2px 7px;display:inline-block;">⚠ ClinGen evidence: offline copy${prov.builtAt ? ` · as of ${esc(prov.builtAt)}` : ''}</div>`
    : '';
  el.innerHTML = offlineBadge + blocks + `<div style="margin-top:5px;font-size:0.62rem;color:var(--text-dim);">© ClinGen — Gene-Disease Validity</div>`;
}

// ── CURATION SUPPORT: concordance between the VARIANT/PATIENT and the DISEASE ─────
// Two decision-support checks that pair the query against the gene-disease knowledge:
//   A. variant molecular consequence ↔ disease mechanism (the PVS1 question)
//   B. patient phenotype ↔ disease phenotype frequency (the PP4 question)
// Both evaluators are PURE (inputs passed explicitly, no data.* reads — CLAUDE.md Rule 2).

/**
 * evaluateMechanismConcordance({ vepConsequence, mechanism, haplo })  — PURE
 * Does the query variant's consequence match the disease's known mechanism? A null/LoF variant
 * only earns PVS1 where loss-of-function is an established mechanism; a LoF variant against a
 * gain-of-function-only disease is a PVS1 caveat.
 * Returns { variantClass, mechanismClasses[], verdict, note }. variantClass is null in gene-only
 * mode (no variant). verdict ∈ concordant | plausible | discordant | inconclusive.
 */
function evaluateMechanismConcordance({ vepConsequence, mechanism, haplo } = {}) {
  const cons = String(vepConsequence || '').toLowerCase();
  let variantClass = null;
  if (cons) {
    if (typeof NULL_CONS !== 'undefined' && NULL_CONS.some(c => cons.includes(c))) variantClass = 'LoF';
    else if (cons.includes('missense')) variantClass = 'missense';
    else if (cons.includes('inframe')) variantClass = 'in-frame indel';
    else if (cons.includes('synonymous')) variantClass = 'synonymous';
    else variantClass = 'other';
  }

  const m = String(mechanism || '').toLowerCase();
  const h = String(haplo ?? '').toLowerCase();
  const mech = new Set();
  if (/loss[\s-]?of[\s-]?function|haploinsuff|nonsense[\s-]?mediated|null allele|loss of the|tumou?r suppressor/.test(m)) mech.add('LoF');
  if (Number(haplo) === 3 || /haploinsuff|sufficient evidence/.test(h)) mech.add('LoF'); // ClinGen dosage corroboration
  if (/gain[\s-]?of[\s-]?function|activating|constitutive activ|increased (?:activity|function)/.test(m)) mech.add('GoF');
  if (/dominant[\s-]?negative/.test(m)) mech.add('DN');
  const mechanismClasses = [...mech];

  let verdict = 'inconclusive', note = '';
  if (!variantClass) {
    note = 'No variant in query — mechanism shown for reference.';
  } else if (!mechanismClasses.length) {
    note = mechanism ? 'Disease mechanism not clearly stated in GeneReviews.' : 'No GeneReviews mechanism available.';
  } else if (variantClass === 'LoF') {
    if (mech.has('LoF')) { verdict = 'concordant'; note = 'Null/LoF variant + loss-of-function disease mechanism — supports PVS1 applicability.'; }
    else { verdict = 'discordant'; note = `Disease mechanism is ${mechanismClasses.join('/')}, not loss-of-function — apply PVS1 with caution (LoF may be tolerated here).`; }
  } else if (variantClass === 'missense' || variantClass === 'in-frame indel') {
    if (mech.has('GoF') || mech.has('DN')) { verdict = 'plausible'; note = `${variantClass} is consistent with a ${mechanismClasses.join('/')} mechanism.`; }
    else if (mech.has('LoF')) { verdict = 'plausible'; note = `${variantClass} can act via loss-of-function (e.g. destabilization) — consistent with the disease mechanism.`; }
  } else {
    note = `${variantClass} variant — mechanism concordance not informative.`;
  }
  return { variantClass, mechanismClasses, verdict, note };
}

/**
 * matchPhenotypeFrequencies({ phenotypeFrequencies, ptPhenoTexts, ptPhenoGroups })  — PURE
 * Pairs the patient's phenotype terms against the GeneReviews "Frequency of Select Features" table.
 * A patient term that is a HIGH-frequency feature of the disease strengthens phenotype fit (PP4).
 * ptPhenoGroups (parallel HPO ids) is accepted for future HPO-id matching; matching is text-based now.
 * Returns [{ ptTerm, matched, feature, freq, freqPct, comment }] — one entry per patient term.
 */
function matchPhenotypeFrequencies({ phenotypeFrequencies, ptPhenoTexts } = {}) {
  const freqs = Array.isArray(phenotypeFrequencies) ? phenotypeFrequencies : [];
  const terms = Array.isArray(ptPhenoTexts) ? ptPhenoTexts : [];
  if (!freqs.length || !terms.length) return [];

  const STOP = new Set(['with', 'type', 'this', 'that', 'from', 'have', 'been', 'were', 'also', 'other', 'more', 'than', 'only', 'into', 'onto', 'such', 'disease', 'syndrome', 'disorder', 'abnormal', 'abnormality']);
  const tokens = s => String(s || '')
    .normalize('NFD').replace(/[̀-ͯ]/g, '')          // strip diacritics (café → cafe)
    .toLowerCase().replace(/[^a-z0-9\s]/g, ' ')
    .split(/\s+/).filter(w => w.length >= 4 && !STOP.has(w))
    .map(w => w.length >= 5 && w.endsWith('s') ? w.slice(0, -1) : w);  // light de-plural (macules→macule)
  const overlaps = (a, b) => {
    const ta = tokens(a), tb = tokens(b);
    if (!ta.length || !tb.length) return false;
    const [s, l] = ta.length <= tb.length ? [ta, tb] : [tb, ta];
    return s.every(w => l.includes(w));                         // every significant word of the shorter side present
  };
  const pctOf = freq => { const m = String(freq || '').match(/(\d+(?:\.\d+)?)\s*%/); return m ? parseFloat(m[1]) : null; };

  return terms.map(ptTerm => {
    const hit = freqs.find(f => overlaps(ptTerm, f.feature));
    return hit
      ? { ptTerm, matched: true, feature: hit.feature, freq: hit.freq, freqPct: pctOf(hit.freq), comment: hit.comment || '', disease: hit.disease || null }
      : { ptTerm, matched: false, feature: null, freq: null, freqPct: null, comment: '', disease: null };
  });
}

/**
 * renderCurationSupportCard()
 * Renders data.curationSupport into #eCurationSupportBody; shows/hides #eCurationSupportCard.
 * Reads ONLY its own namespace (Rule 9). Hidden when neither check has anything to say.
 */
function renderCurationSupportCard() {
  const card = document.getElementById('eCurationSupportCard');
  const el = document.getElementById('eCurationSupportBody');
  if (!card || !el) return;
  const cs = data.curationSupport || {};
  const mech = cs.mechanism || null;
  const matches = Array.isArray(cs.phenoMatches) ? cs.phenoMatches : [];
  const matched = matches.filter(m => m.matched);

  const hasMech = !!(mech && mech.variantClass && (mech.mechanismClasses?.length || mech.verdict !== 'inconclusive'));
  const hasPheno = matches.length > 0;
  card.style.display = 'none'; return; // hidden by request

  const esc = s => String(s || '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const V = { concordant: '#2e7d32', plausible: '#26a69a', discordant: '#f44336', inconclusive: '#90a4ae' };
  let html = '';

  if (hasMech) {
    const color = V[mech.verdict] || '#90a4ae';
    const label = mech.verdict.charAt(0).toUpperCase() + mech.verdict.slice(1);
    html += `<div style="margin-bottom:${hasPheno ? '6px' : '0'};">
      <div style="display:flex;align-items:center;gap:6px;flex-wrap:wrap;">
        <span style="background:${color};color:#fff;border-radius:8px;padding:1px 7px;font-size:0.62rem;font-weight:800;white-space:nowrap;">${label}</span>
        <span style="color:var(--dim);font-size:0.68rem;">Variant <b style="color:var(--text-bright);">${esc(mech.variantClass || '—')}</b> vs mechanism <b style="color:var(--text-bright);">${esc(mech.mechanismClasses?.join('/') || '—')}</b></span>
      </div>
      ${mech.note ? `<div style="color:var(--dim);font-size:0.66rem;margin-top:2px;">${esc(mech.note)}</div>` : ''}
    </div>`;
  }

  if (hasPheno) {
    const pctColor = p => p == null ? '#607d8b' : p >= 75 ? '#2e7d32' : p >= 40 ? '#66bb6a' : '#ff8a80';
    const sep = hasMech ? 'border-top:1px dashed rgba(255,255,255,0.08);padding-top:5px;' : '';
    if (matched.length) {
      html += `<div style="${sep}">
        <div style="color:var(--dim);font-size:0.66rem;font-weight:600;margin-bottom:2px;">Patient phenotype frequency (GeneReviews):</div>` +
        matched.map(m => {
          const dz = (m.disease && m.disease !== cs.chapterTitle)
            ? ` <span style="color:var(--dim);font-size:0.58rem;">(${esc(m.disease)})</span>` : '';
          return `<div style="font-size:0.70rem;color:var(--text-bright);margin-top:1px;">
          <span style="background:${pctColor(m.freqPct)};color:#fff;border-radius:6px;padding:0 5px;font-size:0.60rem;font-weight:700;">${esc(m.freq || '—')}</span>
          ${esc(m.ptTerm)} <span style="color:var(--dim);">→ ${esc(m.feature)}</span>${dz}
        </div>`;
        }).join('') + `</div>`;
    } else {
      html += `<div style="${sep}color:var(--dim);font-size:0.66rem;">No patient phenotype matched a listed GeneReviews feature for this gene.</div>`;
    }
  }

  el.innerHTML = html;
}

async function fetchClinGenInterpretationId(gene, hgvs) {
  if (!gene || !hgvs || hgvs === '-') return null;
  try {
    const url = `https://erepo.clinicalgenome.org/evrepo/api/summary/classifications?columns=preferredVarTitle&values=${encodeURIComponent(hgvs)}&matchTypes=contains&matchLimit=20`;
    const res = await fetchWithTimeout(url);
    if (!res.ok) return null;
    const json = await res.json();
    const interps = json.data || [];
    const match = interps.find(r => r.gene === gene);
    return match?.uuid || null;
  } catch (e) {
    console.warn('[ClinGen eRepo] Interpretation ID fetch failed:', e.message);
    return null;
  }
}

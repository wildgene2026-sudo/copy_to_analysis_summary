/**
 * codon-viewer.js — DECIPHER-style protein variant viewer
 *
 * Layout (top → bottom):
 *   1. Exon strip            (alternating colour blocks)
 *   2. Missense P/LP lane   (red  — bars grow UP from domain strip)
 *   3. Missense VUS  lane   (grey — bars grow UP)
 *   4. Missense B/LB lane   (green— bars grow UP)
 *   5. Domain / protein strip (UniProt domains — backbone)
 *   6. LOF track             (stacked, bars grow UP from x-axis toward strip)
 *   7. X-axis with codon ticks
 *
 * Public API:
 *   renderCodonViewer(gene, codonData, domains, clinvarDate, options?)
 *     options.queryCodon — codon to highlight (number)
 *     options.exons      — [{number, codonStart, codonEnd}] from api.js
 *   window.zoomCodonViewer(factor)         // 1.5=in, 1/1.5=out, 0=reset
 *   window.clinvarCodonSearch(gene, codon) // opens ClinVar in new tab
 *   window.toggleCodonViewer()             // collapse / expand
 */

// ── Palette ────────────────────────────────────────────────────────────────
const CV_COLOR = { plp:'#e53935', vus:'#78909c', blb:'#43a047' };
const CV_LABEL = { plp:'P / LP', vus:'VUS', blb:'B / LB' };
const CV_ORDER = ['plp','vus','blb'];

const DOMAIN_BG = {
  critical: { fill:'#ef5350', stroke:'#b71c1c' },
  high:     { fill:'#42a5f5', stroke:'#1565c0' },
  moderate: { fill:'#66bb6a', stroke:'#2e7d32' },
  low:      { fill:'#90a4ae', stroke:'#546e7a' },
};

const QUERY_COLOR  = '#ffeb3b';
const EXON_FILL_A  = 'rgba(99,179,237,0.72)';   // odd  exons — brighter
const EXON_FILL_B  = 'rgba(99,179,237,0.32)';   // even exons

// ── Layout constants ────────────────────────────────────────────────────────
// top → bottom: Exon | Miss PLP | Miss VUS | Miss BLB | Domain | LOF | Axis
const LABEL_W     = 88;
const PAD_R       = 12;
const EXON_H      = 24;
const MISS_LANE_H = 52;
const DOMAIN_H    = 30;
const LOF_H       = 110;
const AXIS_H      = 22;

const Y_EXON     = 0;
const Y_LANE_PLP = EXON_H;                          // 24
const Y_LANE_VUS = Y_LANE_PLP + MISS_LANE_H;        // 76
const Y_LANE_BLB = Y_LANE_VUS + MISS_LANE_H;        // 128
const Y_DOMAIN   = Y_LANE_BLB + MISS_LANE_H;        // 180  ← protein backbone
const Y_LOF      = Y_DOMAIN   + DOMAIN_H;           // 210  ← LOF track top
const Y_AXIS     = Y_LOF      + LOF_H;              // 320
const TOTAL_H    = Y_AXIS     + AXIS_H;             // 342

// ── Module state ────────────────────────────────────────────────────────────
let _viewerArgs      = null;
let _viewerZoom      = 1.0;
let _viewerCollapsed = false;
let _plpFilter       = 'any'; // 'any' | '3plus' | '4only'
const ZOOM_MIN = 1.0, ZOOM_MAX = 12.0;

// Region → ClinVar search (codon range). null until first render computes a default.
let _regionFrom = null, _regionTo = null;
// gene:queryCodon of the last render — region defaults reset only when this changes,
// so a Phase-B re-render (same variant, +exons/domains) preserves the user's selection.
let _lastRenderKey = null;
// Scratch captured each render so the live slider handler can update the overlay
// + readout without a full re-render: { proteinLen, plotW, variantCodons:Set }.
let _cvState = null;
// Protein-VARNAME OR-lists get unwieldy past this many codons; cap the search.
const REGION_MAX = 60;

// ── Helpers ─────────────────────────────────────────────────────────────────
function esc(s) {
  return String(s)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;')
    .replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}
function totalSubs(entry, track) {
  const t = entry?.[track];
  return t ? CV_ORDER.reduce((s,k) => s + (t[k]?.s ?? 0), 0) : 0;
}

// ── ClinVar codon search ────────────────────────────────────────────────────
// Resolve a codon to its p.<aa><pos> VARNAME label.
//   aa is the 3-letter reference amino acid (e.g. 'Ala') stored in codon data
//   as codonData[pos].aa — present after data regeneration with updated pipeline.
//   Codon 1 is the canonical start Met in the reference but ships without an `aa`
//   field, so default it to Met rather than falling back to the bare-number search
//   (which is what kept start-loss / codon-1 LOF bars from finding their variants).
function _codonProteinLabel(codon, aa) {
  const a = aa || (Number(codon) === 1 ? 'Met' : null);
  return a ? `p.${a}${codon}` : null;
}

window.clinvarCodonSearch = function(gene, codon, aa) {
  if (!gene || !codon) return;
  // Format: p.Ala574  →  "p.Ala574"[VARNAME] AND "FLCN"[GENE]
  // Fallback (no label): "574"[VARNAME] AND "FLCN"[GENE]
  const pLabel = _codonProteinLabel(codon, aa);
  const term   = pLabel
    ? `"${pLabel}"[VARNAME] AND "${gene}"[GENE]`
    : `"${codon}"[VARNAME] AND "${gene}"[GENE]`;
  const base   = `https://www.ncbi.nlm.nih.gov/clinvar/`;
  const params = new URLSearchParams({ gene, term });
  if (pLabel) params.set('variant', pLabel);
  window.open(`${base}?${params.toString()}`, '_blank', 'noopener');
};

// ── Codon → genomic mapping (for ClinVar chrpos38 region search) ──────────────
// _cvState.geno (built each render from the transcript exons) carries:
//   { chrom, strand, exons:[{cdsStart,cdsEnd,gc1,gc2}], maxCds }
// gc1/gc2 are the GRCh38 coords of each exon's first/last coding base.
// Map a c. (CDS) base position to its genomic coordinate by locating its exon
// and linearly interpolating gc1→gc2 across cdsStart→cdsEnd.
function _cBaseToGenomic(b, geno) {
  for (const e of geno.exons) {
    if (b >= e.cdsStart && b <= e.cdsEnd) {
      const span = Math.max(e.cdsEnd - e.cdsStart, 1);
      return Math.round(e.gc1 + (e.gc2 - e.gc1) * ((b - e.cdsStart) / span));
    }
  }
  return null;
}
// Codon range [lo,hi] → GRCh38 window {chrom,start,end}. The min/max over the two
// endpoints spans any introns in between — so the window also catches splice/
// intronic variants. padBp widens each edge to catch boundary indels.
function _codonRangeToGenomic(lo, hi, geno, padBp = 0) {
  if (!geno || !geno.exons.length) return null;
  const bLo = Math.max(1, (lo - 1) * 3 + 1);
  const bHi = Math.min(geno.maxCds, hi * 3);
  const g1 = _cBaseToGenomic(bLo, geno), g2 = _cBaseToGenomic(bHi, geno);
  if (g1 == null || g2 == null) return null;
  return { chrom: geno.chrom, start: Math.min(g1, g2) - padBp, end: Math.max(g1, g2) + padBp };
}
function _genomicClinvarUrl(gene, w) {
  const term = `${gene}[gene] AND ${w.chrom}[chr] AND ${w.start}:${w.end}[chrpos38]`;
  return `https://www.ncbi.nlm.nih.gov/clinvar/?${new URLSearchParams({ term }).toString()}`;
}

// ── ClinVar LOF codon search (genomic) ───────────────────────────────────────
// LOF bars count coding indels/frameshifts ClinVar often names by cDNA only
// (e.g. "c.214del" → no p.<aa><pos>), which a protein-name search can't find.
// Search the codon's genomic window instead — finds them regardless of name.
// Falls back to the protein-name search if exon/genomic data isn't loaded yet.
window.clinvarCodonLofSearch = function(gene, codon, aa) {
  if (!gene || !codon) return;
  const geno = _cvState?.geno;
  if (geno) {
    const w = _codonRangeToGenomic(codon, codon, geno, 3);
    if (w) { window.open(_genomicClinvarUrl(gene, w), '_blank', 'noopener'); return; }
  }
  window.clinvarCodonSearch(gene, codon, aa);  // fallback: protein-name
};

// ── ClinVar region search (codon range) ──────────────────────────────────────
// Preferred: map the codon range → GRCh38 window and search by chrpos38, which
// finds EVERY variant type in the region (missense, nonsense, frameshift/indel,
// splice) regardless of how ClinVar names it.
// Fallback (no exon/genomic data yet): OR-list of per-codon protein VARNAME terms
// — same convention as the per-bar links, but misses splice & c.-only indels and
// is capped at REGION_MAX.
window.clinvarRegionSearch = function() {
  if (!_viewerArgs || !_regionFrom || !_regionTo) return;
  const { gene, codonData } = _viewerArgs;
  if (!gene || !codonData) return;
  const lo = Math.min(_regionFrom, _regionTo);
  const hi = Math.max(_regionFrom, _regionTo);

  // ── Preferred: genomic chrpos38 window ──────────────────────────────────
  const geno = _cvState?.geno;
  if (geno) {
    const w = _codonRangeToGenomic(lo, hi, geno, 5);
    if (w) { window.open(_genomicClinvarUrl(gene, w), '_blank', 'noopener'); return; }
  }

  // ── Fallback: protein VARNAME OR-list ───────────────────────────────────
  const withVar = [];
  for (let c = lo; c <= hi; c++) if (_cvState?.variantCodons?.has(c)) withVar.push(c);
  let codons = withVar.length ? withVar : [];
  if (!codons.length) for (let c = lo; c <= hi && codons.length < REGION_MAX; c++) codons.push(c);
  if (codons.length > REGION_MAX) codons = codons.slice(0, REGION_MAX);

  const terms = codons.map(c => {
    const lbl = _codonProteinLabel(c, codonData[c]?.aa);
    return lbl ? `"${lbl}"[VARNAME]` : `"${c}"[VARNAME]`;
  });
  if (!terms.length) return;
  const term   = `(${terms.join(' OR ')}) AND "${gene}"[GENE]`;
  const params = new URLSearchParams({ gene, term });
  window.open(`https://www.ncbi.nlm.nih.gov/clinvar/?${params.toString()}`, '_blank', 'noopener');
};

// Live slider handler — updates readout, slider fill, and SVG band without a
// full re-render (oninput fires rapidly during drag).
window._cvRegionInput = function() {
  const fromEl = document.getElementById('cvRegFrom');
  const toEl   = document.getElementById('cvRegTo');
  if (!fromEl || !toEl || !_cvState) return;
  const lo = Math.min(+fromEl.value, +toEl.value);
  const hi = Math.max(+fromEl.value, +toEl.value);
  _regionFrom = lo; _regionTo = hi;

  let n = 0;
  for (let c = lo; c <= hi; c++) if (_cvState.variantCodons.has(c)) n++;
  const ro = document.getElementById('cvRegionReadout');
  if (ro) {
    // The REGION_MAX cap only applies to the protein-VARNAME fallback; genomic
    // chrpos38 search has no per-codon limit.
    const capNote = (!_cvState.geno && n > REGION_MAX) ? ` · first ${REGION_MAX} searched` : '';
    ro.textContent = `Codons ${lo}–${hi} · ${n} with ClinVar variant${n === 1 ? '' : 's'}${capNote}`;
  }
  const len  = Math.max(_cvState.proteinLen - 1, 1);
  const fill = document.getElementById('cvRegionFill');
  if (fill) {
    const a = (lo - 1) / len * 100, b = (hi - 1) / len * 100;
    fill.style.left = a + '%'; fill.style.width = (b - a) + '%';
  }
  const ov = document.getElementById('cvRegionOverlay');
  if (ov) {
    const xOf = c => LABEL_W + ((c - 1) / len) * _cvState.plotW;
    const x1 = xOf(lo), x2 = xOf(hi);
    ov.setAttribute('x', Math.min(x1, x2).toFixed(1));
    ov.setAttribute('width', Math.max(2, Math.abs(x2 - x1)).toFixed(1));
  }
};

// ── Collapse / expand ────────────────────────────────────────────────────────
window.toggleCodonViewer = function() {
  _viewerCollapsed = !_viewerCollapsed;
  const plot = document.getElementById('codonViewerPlot');
  const row  = document.getElementById('cvRegionRow');
  const btn  = document.getElementById('codonViewerToggleBtn');
  if (plot) plot.style.display = _viewerCollapsed ? 'none' : '';
  if (row)  row.style.display  = _viewerCollapsed ? 'none' : 'flex';
  if (btn)  btn.textContent    = _viewerCollapsed ? '▶ Show' : '▼ Hide';
};

// ── P/LP star filter ─────────────────────────────────────────────────────────
window.setCodonViewerPlpFilter = function(val) {
  _plpFilter = val;
  _renderViewerInner();
};

// ── Zoom ─────────────────────────────────────────────────────────────────────
window.zoomCodonViewer = function(factor) {
  const prev = _viewerZoom;
  _viewerZoom = factor === 0
    ? 1.0
    : Math.max(ZOOM_MIN, Math.min(ZOOM_MAX, _viewerZoom * factor));
  if (Math.abs(prev - _viewerZoom) < 0.001) return;
  _renderViewerInner();
  if (_viewerArgs?.queryCodon) {
    requestAnimationFrame(() => {
      const sc = document.getElementById('codonViewerScroller');
      if (!sc) return;
      const baseW = sc.clientWidth || 960;
      const plotW = baseW * _viewerZoom - LABEL_W - PAD_R;
      const x = LABEL_W + ((_viewerArgs.queryCodon - 1) /
        Math.max(_viewerArgs.proteinLen - 1, 1)) * plotW;
      sc.scrollLeft = Math.max(0, x - baseW / 2);
    });
  }
};

// ── Public entry point ────────────────────────────────────────────────────────
function renderCodonViewer(gene, codonData, domains, clinvarDate, options = {}) {
  const queryCodon = Number(options.queryCodon) || null;
  _viewerArgs = {
    gene, codonData, domains: domains || [], clinvarDate,
    queryCodon,
    exons:      options.exons || [],
    proteinLen: 0,
  };
  _viewerZoom = 1.0;
  // Reset the region default only for a new gene/variant — not for Phase-B upgrades
  // of the same one, which would otherwise wipe a region the user just selected.
  const renderKey = `${gene}:${queryCodon}`;
  if (renderKey !== _lastRenderKey) { _regionFrom = null; _regionTo = null; }
  _lastRenderKey = renderKey;
  _renderViewerInner();
}

// ── Inner render ──────────────────────────────────────────────────────────────
function _renderViewerInner() {
  const container = document.getElementById('codonViewerContainer');
  if (!container || !_viewerArgs) return;

  const { gene, codonData, domains, clinvarDate, queryCodon, exons } = _viewerArgs;

  const codonKeys = Object.keys(codonData || {})
    .filter(k => k !== '_meta' && !isNaN(k)).map(Number);

  if (!codonData || codonKeys.length === 0) { container.style.display='none'; return; }

  // ── P/LP star-filter resolution ──────────────────────────────────────────
  const plpKey = _plpFilter === '4only' ? 'plp_4s'
               : _plpFilter === '3plus' ? 'plp_3s'
               : _plpFilter === '2plus' ? 'plp_2s'
               : _plpFilter === '1plus' ? 'plp_1s'
               : 'plp';
  const hasStarData = codonKeys.some(k =>
    codonData[k]?.missense?.plp_1s !== undefined || codonData[k]?.lof?.plp_1s !== undefined
  );
  const effectivePlpKey = (plpKey !== 'plp' && !hasStarData) ? 'plp' : plpKey;
  // Helper: resolve class key accounting for plp substitution
  const rk = cls => cls === 'plp' ? effectivePlpKey : cls;

  const baseW  = container.clientWidth || 960;
  const outerW = Math.max(baseW, Math.floor(baseW * _viewerZoom));
  const plotW  = outerW - LABEL_W - PAD_R;

  const proteinLen = Math.max(...codonKeys, 100);
  _viewerArgs.proteinLen = proteinLen;

  const xOf   = c  => LABEL_W + ((c - 1) / Math.max(proteinLen - 1, 1)) * plotW;
  const BAR_W = Math.max(2.5, Math.min(12, plotW / proteinLen * 2.2));

  // ── Region → ClinVar search state ─────────────────────────────────────────
  // Codons carrying at least one ClinVar variant (any class, missense or LOF) —
  // drives the readout count and the search OR-list.
  const variantCodons = new Set();
  const anyN = t => t && CV_ORDER.some(c => (t[c]?.n ?? 0) > 0);
  for (const k of codonKeys) {
    if (anyN(codonData[k].missense) || anyN(codonData[k].lof)) variantCodons.add(k);
  }
  // Default region: a window around the query codon, else the protein start.
  if (_regionFrom == null || _regionTo == null) {
    if (queryCodon) {
      _regionFrom = Math.max(1, queryCodon - 10);
      _regionTo   = Math.min(proteinLen, queryCodon + 10);
    } else {
      _regionFrom = 1;
      _regionTo   = Math.min(proteinLen, 50);
    }
  }
  // Clamp to current protein bounds (a re-render may use a different length).
  _regionFrom = Math.max(1, Math.min(_regionFrom, proteinLen));
  _regionTo   = Math.max(_regionFrom, Math.min(_regionTo, proteinLen));

  // Codon → genomic map (GRCh38) from the transcript exons, for chrpos38 search.
  // null until Phase-B exon data with genomic coords arrives; until then LOF bars
  // and the region search fall back to protein-name search.
  const genoExons = (exons || []).filter(e =>
    Number.isFinite(e.gc1) && Number.isFinite(e.gc2) &&
    Number.isFinite(e.cdsStart) && Number.isFinite(e.cdsEnd) && e.chrom != null);
  const geno = genoExons.length ? {
    chrom:  String(genoExons[0].chrom).replace(/^chr/i, ''),
    strand: genoExons[0].strand,
    exons:  genoExons,
    maxCds: Math.max(...genoExons.map(e => e.cdsEnd)),
  } : null;

  _cvState = { proteinLen, plotW, variantCodons, geno };
  let regionN = 0;
  for (let c = _regionFrom; c <= _regionTo; c++) if (variantCodons.has(c)) regionN++;

  // ── Scale maxes ─────────────────────────────────────────────────────────
  let lofMax = 1;
  const missMax = { plp:1, vus:1, blb:1 };
  for (const k of codonKeys) {
    const lofTot = CV_ORDER.reduce((s, cls) => s + (codonData[k]?.lof?.[rk(cls)]?.s ?? 0), 0);
    if (lofTot > lofMax) lofMax = lofTot;
    for (const cls of CV_ORDER) {
      const s = codonData[k].missense?.[rk(cls)]?.s ?? 0;
      if (s > missMax[cls]) missMax[cls] = s;
    }
  }
  // Missense: bars grow UP from lane bottom → height proportional to count
  const yOfMiss = (s, cls) => Math.max(1, (s / missMax[cls]) * (MISS_LANE_H - 6));
  // LOF: bars grow UP from Y_AXIS toward Y_LOF → same scale
  const yOfLof  = s => Math.max(1, (s / lofMax) * (LOF_H - 8));

  const missLaneY = cls => cls==='plp' ? Y_LANE_PLP : cls==='vus' ? Y_LANE_VUS : Y_LANE_BLB;

  // ── Background panels ────────────────────────────────────────────────────
  const trackBg = `
    <rect x="${LABEL_W}" y="${Y_LANE_PLP}" width="${plotW}" height="${MISS_LANE_H}"
      fill="rgba(229,57,53,0.045)" rx="2"/>
    <rect x="${LABEL_W}" y="${Y_LANE_VUS}" width="${plotW}" height="${MISS_LANE_H}"
      fill="rgba(120,144,156,0.045)" rx="2"/>
    <rect x="${LABEL_W}" y="${Y_LANE_BLB}" width="${plotW}" height="${MISS_LANE_H}"
      fill="rgba(67,160,71,0.045)" rx="2"/>
    <rect x="${LABEL_W}" y="${Y_LOF}" width="${plotW}" height="${LOF_H}"
      fill="rgba(255,255,255,0.02)" rx="2"/>`;

  // ── Grid lines ────────────────────────────────────────────────────────────
  const GRID = 4;
  let gridSvg = '';

  // LOF: bars grow UP from Y_AXIS — grid lines rise from Y_AXIS
  for (let i = 1; i <= GRID; i++) {
    const frac     = i / GRID;
    const yLof     = Y_AXIS - LOF_H * frac;   // ← flipped: count from bottom
    const lofLabel = Math.round(lofMax * frac);
    gridSvg += `
      <line x1="${LABEL_W}" y1="${yLof.toFixed(1)}" x2="${outerW - PAD_R}" y2="${yLof.toFixed(1)}"
        stroke="rgba(255,255,255,0.05)" stroke-width="1"/>
      <text x="${LABEL_W - 4}" y="${(yLof - 2).toFixed(1)}" text-anchor="end"
        font-size="7" fill="var(--dim)" font-family="Inter,sans-serif">${lofLabel}</text>`;
  }

  // Missense: bars grow UP from lane bottom → max at top, 0 at bottom
  for (const cls of CV_ORDER) {
    const yLane = missLaneY(cls);
    gridSvg += `
      <line x1="${LABEL_W}" y1="${(yLane + MISS_LANE_H/2).toFixed(1)}"
            x2="${outerW - PAD_R}" y2="${(yLane + MISS_LANE_H/2).toFixed(1)}"
        stroke="rgba(255,255,255,0.03)" stroke-width="1" stroke-dasharray="2,3"/>
      <text x="${LABEL_W - 4}" y="${(yLane + 9).toFixed(1)}" text-anchor="end"
        font-size="7" fill="var(--dim)" font-family="Inter,sans-serif">${missMax[cls]}</text>
      <text x="${LABEL_W - 4}" y="${(yLane + MISS_LANE_H - 2).toFixed(1)}" text-anchor="end"
        font-size="7" fill="var(--dim)" font-family="Inter,sans-serif">0</text>`;
  }

  // ── Bars + hit zones ──────────────────────────────────────────────────────
  let lofBars='', missBars='', hitZones='';
  const hitW = Math.max(BAR_W * 1.8, 8);

  for (const k of codonKeys) {
    const d  = codonData[k];
    const x  = xOf(k);
    const aa = d.aa || null;  // 3-letter ref amino acid, e.g. 'Ala' (present after pipeline regeneration)
    const aaArg = aa ? `,'${esc(aa)}'` : '';
    const codonLabel = aa ? `p.${aa}${k}` : `Codon ${k}`;

    // Hit zone: full column from missense top to x-axis
    hitZones += `<rect x="${(x - hitW/2 + BAR_W/2).toFixed(1)}" y="${Y_LANE_PLP}"
      width="${hitW.toFixed(1)}" height="${(Y_AXIS - Y_LANE_PLP).toFixed(1)}"
      fill="transparent" style="cursor:pointer;"
      onclick="window.clinvarCodonSearch('${esc(gene)}',${k}${aaArg})">
      <title>${codonLabel} — click to search ClinVar</title></rect>`;

    // LOF stacked — bars grow UPWARD from Y_AXIS toward Y_LOF
    let lofBase = Y_AXIS;
    for (const cls of CV_ORDER) {
      const key = rk(cls);
      const s = d.lof?.[key]?.s ?? 0;
      if (!s) continue;
      const h = yOfLof(s);
      lofBase -= h;
      const n = d.lof[key].n;
      lofBars += `<rect x="${x.toFixed(1)}" y="${lofBase.toFixed(1)}"
        width="${BAR_W.toFixed(1)}" height="${h.toFixed(1)}"
        fill="${CV_COLOR[cls]}" rx="1" style="cursor:pointer;"
        onclick="window.clinvarCodonLofSearch('${esc(gene)}',${k}${aaArg})">
        <title>${codonLabel} · LOF · ${CV_LABEL[cls]}: ${n} variants, ${s} submitters — click to search ClinVar by genomic region</title>
      </rect>`;
    }

    // Missense — 3 sub-tracks, bars grow UPWARD from lane bottom
    for (const cls of CV_ORDER) {
      const key = rk(cls);
      const s = d.missense?.[key]?.s ?? 0;
      if (!s) continue;
      const h    = yOfMiss(s, cls);
      const yLane = missLaneY(cls);
      const yBar  = yLane + MISS_LANE_H - h;  // anchor at lane bottom, grow up
      const n = d.missense[key].n;
      missBars += `<rect x="${x.toFixed(1)}" y="${yBar.toFixed(1)}"
        width="${BAR_W.toFixed(1)}" height="${h.toFixed(1)}"
        fill="${CV_COLOR[cls]}" rx="1" style="cursor:pointer;"
        onclick="window.clinvarCodonSearch('${esc(gene)}',${k}${aaArg})">
        <title>${codonLabel} · Missense · ${CV_LABEL[cls]}: ${n} variants, ${s} submitters</title>
      </rect>`;
    }
  }

  // ── Exon strip ────────────────────────────────────────────────────────────
  const exonBgH = EXON_H - 4;
  let exonSvg = `
    <rect x="${LABEL_W}" y="${Y_EXON + 2}" width="${plotW}" height="${exonBgH}"
      rx="3" fill="rgba(40,60,80,0.55)" stroke="rgba(99,179,237,0.35)" stroke-width="1"/>`;

  if (exons && exons.length) {
    for (const ex of exons) {
      const x1   = xOf(ex.codonStart);
      const x2   = xOf(Math.min(ex.codonEnd, proteinLen));
      const w    = Math.max(x2 - x1, 1);
      const fill = ex.number % 2 === 1 ? EXON_FILL_A : EXON_FILL_B;
      exonSvg += `<rect x="${x1.toFixed(1)}" y="${(Y_EXON + 3).toFixed(1)}"
        width="${w.toFixed(1)}" height="${(exonBgH - 2).toFixed(1)}"
        fill="${fill}" stroke="rgba(150,210,255,0.4)" stroke-width="0.5">
        <title>Exon ${ex.number} · aa ${ex.codonStart}–${ex.codonEnd}</title></rect>`;
      if (w > 16) {
        exonSvg += `<text x="${(x1 + w/2).toFixed(1)}" y="${(Y_EXON + EXON_H/2 + 4).toFixed(1)}"
          text-anchor="middle" font-size="8" font-weight="800"
          fill="rgba(255,255,255,0.95)" font-family="Inter,sans-serif"
          style="pointer-events:none;">${ex.number}</text>`;
      }
    }
  } else {
    exonSvg += `<text x="${(LABEL_W + plotW/2).toFixed(1)}"
      y="${(Y_EXON + EXON_H/2 + 4).toFixed(1)}"
      text-anchor="middle" font-size="8" fill="rgba(150,210,255,0.55)"
      font-family="Inter,sans-serif">exon track — waiting for transcript data…</text>`;
  }

  // ── Domain / protein backbone strip ──────────────────────────────────────
  const domY  = Y_DOMAIN + 2;
  const domH  = DOMAIN_H - 4;
  const domCY = Y_DOMAIN + DOMAIN_H / 2;   // vertical centre

  // Backbone line (protein chain)
  let domainSvg = `
    <rect x="${LABEL_W}" y="${domY}" width="${plotW}" height="${domH}"
      rx="3" fill="rgba(30,40,55,0.7)" stroke="rgba(255,255,255,0.18)" stroke-width="1"/>
    <line x1="${LABEL_W}" y1="${domCY.toFixed(1)}" x2="${outerW - PAD_R}" y2="${domCY.toFixed(1)}"
      stroke="rgba(255,255,255,0.25)" stroke-width="2"/>
    <text x="${LABEL_W + 5}" y="${(domCY + 4).toFixed(1)}"
      font-size="8" fill="rgba(255,255,255,0.45)" font-family="Inter,sans-serif">1</text>
    <text x="${(outerW - PAD_R - 4).toFixed(1)}" y="${(domCY + 4).toFixed(1)}"
      text-anchor="end" font-size="8" fill="rgba(255,255,255,0.45)"
      font-family="Inter,sans-serif">${proteinLen} aa</text>`;

  // Show major structural features on the strip.  UniProt uses 'Domain' for some
  // genes but 'Region' / 'Motif' for others (e.g. TP53).  Include both, but
  // exclude noisy interaction/disordered entries that would swamp the strip.
  const STRIP_TYPES = new Set([
    'Domain','Zinc finger','Transmembrane','Signal peptide','Propeptide',
    'Repeat','Region','Motif'
  ]);
  const STRIP_EXCLUDE = /^(Interaction with |Required for interaction|Disordered)/i;
  const domainFiltered = (domains || []).filter(f =>
    STRIP_TYPES.has(f.type) &&
    !STRIP_EXCLUDE.test(f.name) &&
    f.start >= 1 && f.end <= proteinLen + 20 &&
    (f.end - f.start) >= 3
  );
  // Sort so smaller (more specific) domains draw on top of larger ones
  domainFiltered.sort((a, b) => (b.end - b.start) - (a.end - a.start));

  for (const f of domainFiltered) {
    const x1  = xOf(f.start);
    const x2  = xOf(Math.min(f.end, proteinLen));
    const dw  = Math.max(x2 - x1, 6);
    const col = DOMAIN_BG[f.tier] || DOMAIN_BG.low;
    const lbl = f.name.length > 18 ? f.name.slice(0, 16) + '…' : f.name;
    domainSvg += `
      <rect x="${x1.toFixed(1)}" y="${(domY + 1).toFixed(1)}"
        width="${dw.toFixed(1)}" height="${(domH - 2).toFixed(1)}"
        rx="3" fill="${col.fill}" fill-opacity="0.5"
        stroke="${col.stroke}" stroke-opacity="0.9" stroke-width="1.4">
        <title>${esc(f.name)} (${f.type}) · aa ${f.start}–${f.end}</title>
      </rect>`;
    // Only label if the domain rect is wide enough to fit text
    if (dw > 30) {
      domainSvg += `
        <text x="${(x1 + dw/2).toFixed(1)}" y="${(domCY + 4).toFixed(1)}"
          text-anchor="middle" font-size="8.5" font-weight="800"
          fill="#fff" font-family="Inter,sans-serif"
          style="pointer-events:none;">${esc(lbl)}</text>`;
    }
  }

  // ── Variant indicator (query codon) ──────────────────────────────────────
  let indicator = '';
  if (queryCodon && queryCodon >= 1 && queryCodon <= proteinLen) {
    const xi = xOf(queryCodon);
    indicator = `
      <line x1="${xi.toFixed(1)}" y1="${Y_EXON}" x2="${xi.toFixed(1)}" y2="${Y_AXIS}"
        stroke="${QUERY_COLOR}" stroke-width="1.5" stroke-dasharray="4,3" opacity="0.85"
        style="pointer-events:none;"/>
      <polygon points="${xi},${Y_EXON+1} ${xi-5},${Y_EXON-6} ${xi+5},${Y_EXON-6}"
        fill="${QUERY_COLOR}" stroke="rgba(0,0,0,0.4)" stroke-width="0.5"
        style="pointer-events:none;"/>
      <rect x="${xi+6}" y="${Y_EXON-10}" width="52" height="12" rx="2"
        fill="rgba(255,235,59,0.2)" stroke="${QUERY_COLOR}" stroke-width="0.7"
        style="pointer-events:none;"/>
      <text x="${xi+10}" y="${Y_EXON-1}"
        font-size="8.5" fill="${QUERY_COLOR}" font-weight="800"
        font-family="Inter,sans-serif" style="pointer-events:none;">aa ${queryCodon}</text>`;
  }

  // ── Region band (ClinVar region-search selection) ─────────────────────────
  const rbX1 = xOf(_regionFrom), rbX2 = xOf(_regionTo);
  const regionBand = `
    <rect id="cvRegionOverlay" x="${Math.min(rbX1, rbX2).toFixed(1)}" y="${Y_EXON}"
      width="${Math.max(2, Math.abs(rbX2 - rbX1)).toFixed(1)}" height="${Y_AXIS - Y_EXON}"
      fill="rgba(0,212,170,0.09)" stroke="rgba(0,212,170,0.55)" stroke-width="1"
      stroke-dasharray="3,2" style="pointer-events:none;"/>`;

  // ── X-axis ────────────────────────────────────────────────────────────────
  const tickEvery = proteinLen <= 150 ? 25 : proteinLen <= 400 ? 50 : 100;
  let axisSvg = `<line x1="${LABEL_W}" y1="${Y_AXIS}" x2="${outerW-PAD_R}" y2="${Y_AXIS}"
    stroke="rgba(255,255,255,0.2)" stroke-width="1.2"/>`;
  for (let c = tickEvery; c < proteinLen; c += tickEvery) {
    const x = xOf(c);
    axisSvg += `
      <line x1="${x.toFixed(1)}" y1="${Y_AXIS}" x2="${x.toFixed(1)}" y2="${(Y_AXIS+5).toFixed(1)}"
        stroke="rgba(255,255,255,0.25)" stroke-width="1"/>
      <text x="${x.toFixed(1)}" y="${(Y_AXIS+15).toFixed(1)}" text-anchor="middle"
        font-size="9" fill="var(--dim)" font-family="Inter,sans-serif">${c}</text>`;
  }

  // ── Track labels ──────────────────────────────────────────────────────────
  const lSz = (size, weight, color) =>
    `font-size="${size}" font-weight="${weight}" fill="${color}" font-family="Inter,sans-serif"`;
  const labelSvg = `
    <text x="${LABEL_W-8}" y="${(Y_EXON+EXON_H/2+4).toFixed(1)}" text-anchor="end"
      ${lSz(9,700,'var(--text)')}>Exons</text>

    <text x="${LABEL_W-8}" y="${(Y_LANE_PLP+MISS_LANE_H/2-3).toFixed(1)}" text-anchor="end"
      ${lSz(9,700,CV_COLOR.plp)}>P / LP</text>
    <text x="${LABEL_W-8}" y="${(Y_LANE_PLP+MISS_LANE_H/2+8).toFixed(1)}" text-anchor="end"
      ${lSz(7,400,'var(--dim)')}>missense</text>

    <text x="${LABEL_W-8}" y="${(Y_LANE_VUS+MISS_LANE_H/2-3).toFixed(1)}" text-anchor="end"
      ${lSz(9,700,CV_COLOR.vus)}>VUS</text>
    <text x="${LABEL_W-8}" y="${(Y_LANE_VUS+MISS_LANE_H/2+8).toFixed(1)}" text-anchor="end"
      ${lSz(7,400,'var(--dim)')}>missense</text>

    <text x="${LABEL_W-8}" y="${(Y_LANE_BLB+MISS_LANE_H/2-3).toFixed(1)}" text-anchor="end"
      ${lSz(9,700,CV_COLOR.blb)}>B / LB</text>
    <text x="${LABEL_W-8}" y="${(Y_LANE_BLB+MISS_LANE_H/2+8).toFixed(1)}" text-anchor="end"
      ${lSz(7,400,'var(--dim)')}>missense</text>

    <text x="${LABEL_W-8}" y="${(Y_DOMAIN+DOMAIN_H/2+4).toFixed(1)}" text-anchor="end"
      ${lSz(9,700,'var(--text)')}>Protein</text>

    <text x="${LABEL_W-8}" y="${(Y_LOF+14).toFixed(1)}" text-anchor="end"
      ${lSz(10,700,'var(--text)')}>LOF</text>
    <text x="${LABEL_W-8}" y="${(Y_LOF+25).toFixed(1)}" text-anchor="end"
      ${lSz(7,400,'var(--dim)')}>frameshift</text>
    <text x="${LABEL_W-8}" y="${(Y_LOF+34).toFixed(1)}" text-anchor="end"
      ${lSz(7,400,'var(--dim)')}>nonsense</text>
    <text x="${LABEL_W-8}" y="${(Y_LOF+43).toFixed(1)}" text-anchor="end"
      ${lSz(7,400,'var(--dim)')}>splice</text>`;

  // ── Separators ────────────────────────────────────────────────────────────
  const sep = (y, opacity=0.12, w=1, dash='') =>
    `<line x1="${LABEL_W}" y1="${y}" x2="${outerW-PAD_R}" y2="${y}"
      stroke="rgba(255,255,255,${opacity})" stroke-width="${w}"${dash ? ` stroke-dasharray="${dash}"` : ''}/>`;
  const separators = `
    <line x1="${LABEL_W}" y1="${Y_EXON}" x2="${LABEL_W}" y2="${Y_AXIS}"
      stroke="rgba(255,255,255,0.1)" stroke-width="1"/>
    ${sep(Y_LANE_PLP, 0.12)}
    ${sep(Y_LANE_VUS, 0.07, 0.5, '2,3')}
    ${sep(Y_LANE_BLB, 0.07, 0.5, '2,3')}
    ${sep(Y_DOMAIN, 0.25, 1.5)}
    ${sep(Y_LOF, 0.25, 1.5)}`;

  // ── Legend ────────────────────────────────────────────────────────────────
  const legendItems = CV_ORDER.map((cls, i) => `
    <g transform="translate(${i * 90},0)">
      <rect width="11" height="11" rx="2" fill="${CV_COLOR[cls]}"/>
      <text x="15" y="9.5" font-size="10" fill="var(--text)"
        font-family="Inter,sans-serif" font-weight="500">${CV_LABEL[cls]}</text>
    </g>`).join('');

  // ── Summary note ──────────────────────────────────────────────────────────
  const lofTotal  = codonKeys.reduce((s,k) =>
    s + CV_ORDER.reduce((ss, cls) => ss + (codonData[k]?.lof?.[rk(cls)]?.s ?? 0), 0), 0);
  const missTotal = codonKeys.reduce((s,k) =>
    s + CV_ORDER.reduce((ss, cls) => ss + (codonData[k]?.missense?.[rk(cls)]?.s ?? 0), 0), 0);
  const dateLabel = clinvarDate || codonData._meta?.date || '—';
  const queryNote = queryCodon ? ` · query codon ${queryCodon}` : '';
  const filterNote = effectivePlpKey === 'plp_4s' ? ' · P/LP: ★★★★'
                   : effectivePlpKey === 'plp_3s' ? ' · P/LP: ≥3★'
                   : effectivePlpKey === 'plp_2s' ? ' · P/LP: ≥2★'
                   : effectivePlpKey === 'plp_1s' ? ' · P/LP: ≥1★'
                   : '';

  // ── SVG assembly ──────────────────────────────────────────────────────────
  const svg = `<svg width="${outerW}" height="${TOTAL_H+8}" xmlns="http://www.w3.org/2000/svg"
      style="display:block; overflow:visible;">
    ${trackBg}
    ${gridSvg}
    ${separators}
    ${regionBand}
    ${hitZones}
    ${exonSvg}
    ${missBars}
    ${domainSvg}
    ${lofBars}
    ${axisSvg}
    ${labelSvg}
    ${indicator}
  </svg>`;

  // ── Zoom controls ─────────────────────────────────────────────────────────
  const btnS = `padding:2px 9px; font-size:0.7rem; border-radius:4px;
    border:1px solid rgba(255,255,255,0.18); background:rgba(255,255,255,0.04);
    color:var(--text); cursor:pointer; font-family:Inter,sans-serif;`;

  const isCollapsed = _viewerCollapsed;

  container.innerHTML = `
    <div class="data-card" style="padding:14px 18px;">
      <!-- ── Header ── -->
      <div style="display:flex; justify-content:space-between; align-items:flex-start;
                  flex-wrap:wrap; gap:8px;">
        <div style="display:flex; align-items:center; gap:10px; flex:1; min-width:0;">
          <button id="codonViewerToggleBtn"
            onclick="window.toggleCodonViewer()"
            style="${btnS} font-size:0.75rem; min-width:60px;">
            ${isCollapsed ? '▶ Show' : '▼ Hide'}
          </button>
          <div>
            <h3 style="margin:0; padding:0; border:none; font-size:0.82rem; font-weight:700;">
              🧬 Protein Variant Viewer — ${esc(gene)}
            </h3>
            <div style="font-size:0.6rem; color:var(--dim); margin-top:2px;">
              ClinVar ${esc(dateLabel)} · ${proteinLen} aa${queryNote}${filterNote} ·
              LOF ${lofTotal.toLocaleString()} + Missense ${missTotal.toLocaleString()} submitter-counts ·
              <span style="opacity:0.7;">click bar → ClinVar (missense by p.codon · LOF by genomic region) · noncoding &amp; synonymous not plotted</span>
            </div>
          </div>
        </div>
        <div style="display:flex; flex-direction:column; gap:5px; align-items:flex-end; flex-shrink:0;">
          <div style="display:flex; gap:4px; align-items:center;">
            <button title="Zoom out" style="${btnS}" onclick="window.zoomCodonViewer(1/1.5)">−</button>
            <span style="font-size:0.7rem; color:var(--dim); min-width:38px; text-align:center;">
              ${(_viewerZoom*100).toFixed(0)}%</span>
            <button title="Zoom in"  style="${btnS}" onclick="window.zoomCodonViewer(1.5)">+</button>
            <button title="Reset zoom" style="${btnS}; opacity:0.75;"
              onclick="window.zoomCodonViewer(0)">Reset</button>
          </div>
          <div style="display:flex; align-items:center; gap:5px;">
            <span style="font-size:0.65rem; color:var(--dim); white-space:nowrap;">P/LP filter:</span>
            <select id="codonViewerPlpFilter"
              onchange="window.setCodonViewerPlpFilter(this.value)"
              style="${btnS} padding:1px 5px; font-size:0.65rem;">
              <option value="any"   ${_plpFilter === 'any'   ? 'selected' : ''}>Any</option>
              <option value="1plus" ${_plpFilter === '1plus' ? 'selected' : ''}>1+ stars</option>
              <option value="2plus" ${_plpFilter === '2plus' ? 'selected' : ''}>2+ stars</option>
              <option value="3plus" ${_plpFilter === '3plus' ? 'selected' : ''}>3+ stars</option>
              <option value="4only" ${_plpFilter === '4only' ? 'selected' : ''}>4 stars only</option>
            </select>
            ${!hasStarData && plpKey !== 'plp'
              ? '<span style="font-size:0.6rem; color:var(--dim); white-space:nowrap;">showing all</span>'
              : ''}
          </div>
          <svg width="${CV_ORDER.length*90}" height="14" xmlns="http://www.w3.org/2000/svg">
            ${legendItems}
          </svg>
        </div>
      </div>

      <!-- ── Region → ClinVar search ── -->
      <div id="cvRegionRow" style="margin-top:10px; padding:7px 11px; border-radius:6px;
                  background:rgba(0,212,170,0.05); border:1px solid rgba(0,212,170,0.18);
                  align-items:center; gap:12px; flex-wrap:wrap;
                  display:${isCollapsed ? 'none' : 'flex'};">
        <span style="font-size:0.68rem; font-weight:700; color:var(--text); white-space:nowrap;">
          🔎 Region search</span>
        <div class="cv-region" style="flex:1; min-width:170px;">
          <div class="cv-region-track"></div>
          <div class="cv-region-fill" id="cvRegionFill"
            style="left:${((_regionFrom - 1) / Math.max(proteinLen - 1, 1) * 100).toFixed(2)}%;
                   width:${((_regionTo - _regionFrom) / Math.max(proteinLen - 1, 1) * 100).toFixed(2)}%;"></div>
          <input type="range" id="cvRegFrom" min="1" max="${proteinLen}"
            value="${_regionFrom}" oninput="window._cvRegionInput()"
            aria-label="Region start codon">
          <input type="range" id="cvRegTo" min="1" max="${proteinLen}"
            value="${_regionTo}" oninput="window._cvRegionInput()"
            aria-label="Region end codon">
        </div>
        <span id="cvRegionReadout" style="font-size:0.66rem; color:var(--dim);
              white-space:nowrap; min-width:150px; text-align:right;">
          Codons ${_regionFrom}–${_regionTo} · ${regionN} with ClinVar variant${regionN === 1 ? '' : 's'}${(!geno && regionN > REGION_MAX) ? ` · first ${REGION_MAX} searched` : ''}</span>
        <button onclick="window.clinvarRegionSearch()"
          style="${btnS} background:var(--teal-d); border-color:rgba(0,212,170,0.5);
                 font-weight:700; white-space:nowrap;">Search region → ClinVar</button>
      </div>

      <!-- ── Plot (collapsible) ── -->
      <div id="codonViewerPlot"
           style="margin-top:12px; overflow-x:auto; overflow-y:visible; padding-top:12px;
                  ${isCollapsed ? 'display:none;' : ''}">
        ${svg}
      </div>
    </div>`;

  container.style.display = 'block';
}

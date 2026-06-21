// ── UI RENDERING & DOM MANIPULATION ────────────────────────────────

function initTippyDelegate() {
  if (typeof tippy === 'undefined') return;
  tippy.delegate(document.body, {
    target: '[data-tippy-content]',
    delay: [150, 0],
    duration: [150, 100],
    arrow: true,
    theme: 'genelab',
    placement: 'top',
    maxWidth: 300,
  });
}

document.addEventListener('DOMContentLoaded', initTippyDelegate);

function setDot(id, state) {
  const el = document.getElementById('dot_' + id);
  if (el) el.className = 'api-dot ' + state;

  // Also update the new API badge system
  setAPIStatus(id, state);
}

function renderVCEPBanner() {
  let banner = document.getElementById('vcepBanner');

  // Safe injection: Create it if it doesn't exist
  if (!banner) {
    banner = document.createElement('div');
    banner.id = 'vcepBanner';

    // Attempt to insert it above the ACMG matrix or DB Grid
    const targetAnchor = document.getElementById('dbGrid') || document.getElementById('acmgMatrix');
    if (targetAnchor && targetAnchor.parentElement) {
      targetAnchor.parentElement.insertBefore(banner, targetAnchor);
    } else {
      document.body.prepend(banner); // Fallback
    }
  }

  if (data.vcepGuideline) {
    const uiUrl = `https://cspec.genome.network/cspec/ui/svi/doc/${data.vcepGuideline.cspecId}`;
    banner.style.display = 'flex';
    banner.style.cssText = 'background: rgba(0, 200, 150, 0.08); border-left: 4px solid var(--teal); padding: 12px 16px; border-radius: 6px; margin-bottom: 16px; display: flex; justify-content: space-between; align-items: center;';

    banner.innerHTML = `
      <div style="display: flex; flex-direction: column; gap: 4px;">
        <span style="font-weight: 800; color: var(--teal); font-size: 0.85rem;">
          ✨ ClinGen Gene-Specific Guideline Available
        </span>
        <span style="font-size: 0.75rem; color: var(--text-bright);">
          <b>${data.vcepGuideline.gene}</b> (v${data.vcepGuideline.version}) — Review modified ACMG rules before classification.
        </span>
      </div>
      <a href="${uiUrl}" target="_blank" 
         style="background: var(--teal); color: #000; padding: 6px 14px; border-radius: 4px; font-size: 0.75rem; text-decoration: none; font-weight: 800; white-space: nowrap;">
        View Rules ↗
      </a>
    `;
  } else {
    banner.style.display = 'none';
    banner.innerHTML = '';
  }
}

function fmtAF(v) {
  if (v == null) return '—';
  if (v === 0) return '0 (absent)';
  return (v * 100).toFixed(4) + '%';
}

function stars(n) {
  return n ? '★'.repeat(n) + '☆'.repeat(4 - n) : 'Not rated';
}

function reviewStars(reviewStatus) {
  if (!reviewStatus) return 0;
  const statusStr = reviewStatus.toLowerCase().trim();
  return CLINVAR_STAR_MAP[statusStr] ?? 0;
}

function sigClass(s) {
  if (!s) return '';
  const l = s.toLowerCase();
  if (l.includes('pathogenic') && !l.includes('likely')) return 'sig-p';
  if (l.includes('likely pathogenic')) return 'sig-lp';
  if (l.includes('likely benign')) return 'sig-lb';
  if (l.includes('benign') && !l.includes('likely')) return 'sig-b';
  if (l.includes('uncertain')) return 'sig-vus';
  if (l.includes('conflicting')) return 'sig-con';
  return '';
}

function buildACMGMatrix() {
  const rows = [
    { label: 'PVS/PS', pathogen: true, codes: ['PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 'PP3_S'] },
    { label: 'PM', pathogen: true, codes: ['PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6', 'PP3_M', 'PVS1_M'] },
    { label: 'PP', pathogen: true, codes: ['PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'PVS1_S'] },
    { label: 'BA/BS', pathogen: false, codes: ['BA1', 'BP4_VS', 'BS1', 'BS2', 'BS3', 'BS4', 'BP4_S'] },
    { label: 'BP', pathogen: false, codes: ['BP1', 'BP2', 'BP3', 'BP4', 'BP4_M', 'BP5', 'BP6', 'BP7'] },
  ];
  const container = document.getElementById('acmgMatrix');
  container.innerHTML = '';
  rows.forEach(row => {
    const div = document.createElement('div'); div.className = 'acmg-row';
    const lbl = document.createElement('span'); lbl.className = 'acmg-row-label'; lbl.textContent = row.label; div.appendChild(lbl);
    const dvd = document.createElement('div'); dvd.className = 'acmg-divider'; div.appendChild(dvd);
    row.codes.forEach(cid => {
      const def = ACMG_CODES.find(c => c.id === cid);
      const chip = document.createElement('div'); chip.className = 'chip'; chip.id = 'chip_' + cid;
      chip.dataset.cat = def.cat; chip.dataset.pathogen = row.pathogen;
      chip.dataset.tippyContent = def.desc;
      chip.innerHTML = cid;
      chip.onclick = () => toggleChip(cid, row.pathogen);
      div.appendChild(chip);
    });
    container.appendChild(div);
  });
}

function renderChips() {
  ACMG_CODES.forEach(def => {
    const el = document.getElementById('chip_' + def.id); if (!el) return;
    const isPath = ['pvs', 'ps', 'pm', 'pp'].includes(def.cat);
    el.className = 'chip';
    if (selectedCodes.has(def.id)) el.classList.add(isPath ? 'selected-path' : 'selected-ben');
    else if (suggestedCodes.has(def.id)) el.classList.add(isPath ? 'suggested-path' : 'suggested-ben');
  });
}

function classifyACMG() {
  const sel = [...selectedCodes];
  const pvs = sel.filter(c => c === 'PVS1').length;
  const ps = sel.filter(c => ['PS1', 'PS2', 'PS3', 'PS4', 'PP3_S'].includes(c)).length;
  const pm = sel.filter(c => ['PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6', 'PP3_M', 'PVS1_M'].includes(c)).length;
  const pp = sel.filter(c => ['PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'PVS1_S'].includes(c)).length;
  const ba = sel.find(c => ['BA1', 'BP4_VS'].includes(c)) ? 1 : 0;
  const bs = sel.filter(c => ['BS1', 'BS2', 'BS3', 'BS4', 'BP4_S'].includes(c)).length;
  const bp = sel.filter(c => ['BP1', 'BP2', 'BP3', 'BP4', 'BP4_M', 'BP5', 'BP6', 'BP7'].includes(c)).length;
  if (sel.length === 0) return { label: '— Awaiting Data —', cls: '' };
  if (ba >= 1) return { label: 'Benign', cls: 'ben' };
  if (bs >= 2) return { label: 'Benign', cls: 'ben' };
  if (bs >= 1 && bp >= 1) return { label: 'Likely Benign', cls: 'lben' };
  if (bp >= 2) return { label: 'Likely Benign', cls: 'lben' };
  const pathogenic = (pvs >= 1 && ps >= 1) || (pvs >= 1 && pm >= 2) || (pvs >= 1 && pm >= 1 && pp >= 2) || (pvs >= 1 && pp >= 4) || (ps >= 2) || (ps >= 1 && pm >= 3) || (ps >= 1 && pm >= 2 && pp >= 2) || (ps >= 1 && pm >= 1 && pp >= 4);
  if (pathogenic) return { label: 'Pathogenic', cls: 'path' };
  const likelyP = (pvs >= 1 && pm >= 1) || (ps >= 1 && (pm >= 1 || pm >= 2)) || (ps >= 1 && pp >= 2) || (pm >= 3) || (pm >= 2 && pp >= 2) || (pm >= 1 && pp >= 4);
  if (likelyP) return { label: 'Likely Pathogenic', cls: 'lpath' };
  return { label: 'VUS', cls: 'vus' };
}

function updateBadge() {
  const { label, cls } = classifyACMG();
  const badge = document.getElementById('classBadge');
  badge.textContent = label; badge.className = cls;
}

function buildDBGrid() {
  const g = document.getElementById('dbGrid'); g.innerHTML = '';
  DB_DEFS.forEach(db => {
    g.innerHTML += `<div class="db-wrapper" id="wrap_${db.id}">
      <span class="tooltip-text" id="tip_${db.id}"></span>
      <input type="checkbox" class="db-check" value="${db.id}">
      <span class="rank-badge" id="rank_${db.id}"></span>
      <button class="db-btn" onclick="launchOne('${db.id}')">${db.label}</button>
    </div>`;
  });
  document.querySelectorAll('.db-check').forEach(chk => {
    chk.checked = false;
    chk.addEventListener('change', function () { handleCheckboxChange(this.value, this.checked); });
  });
}

function updateRanks() {
  document.querySelectorAll('.rank-badge').forEach(el => { el.innerText = ''; el.classList.remove('rank-visible'); });
  selectionOrder.forEach((k, i) => { const b = document.getElementById('rank_' + k); if (b) { b.innerText = i + 1; b.classList.add('rank-visible'); } });
}

function enableWrappers(ids) {
  ids.forEach(id => { const el = document.getElementById('wrap_' + id); if (el) el.classList.remove('disabled'); });
  updateTooltips();
}


function updateTooltips() {
  const c = data.coords.userCoord || {};
  const fID = data.coords.userCoord ? `chr${c.chrom}-${c.pos}-${c.ref}-${c.alt}-hg38` : `${data.coords.hg38String}-hg38`;
  const gV4 = data.coords.userCoord ? `${c.chrom}-${c.pos}-${c.ref}-${c.alt}` : `${data.coords.chrom}-${data.coords.pos38}-${data.coords.ref}-${data.coords.alt}`;
  const broad = data.coords.userCoord ? `chr${c.chrom}-${c.pos}%20${c.ref}%3E${c.alt}` : `chr${data.coords.chrom}-${data.coords.pos38}%20${data.coords.ref}%3E${data.coords.alt}`;
  const _mmIsSplice = /\d[+\-]\d/.test(data.coords.hgvs || '');
  const _mmAA1L = toOneLetterAA(data.ensembl.protein);
  const mmTip = (!_mmIsSplice && data.coords.gene && _mmAA1L)
    ? `Protein: ${data.coords.gene}:p.${_mmAA1L}`
    : (data.ensembl.transcript && data.coords.hgvs && data.coords.hgvs !== '-')
      ? `HGVS: ${data.ensembl.transcript}:${data.coords.hgvs}`
      : `Gene: ${data.coords.gene}`;
  const tips = { franklin: `ID: ${fID}`, clinvar: `HGVS: ${data.coords.hgvs}`, gnomadv4: `ID: ${gV4}`, spliceai: `Variant: ${broad}`, clingen: `Gene: ${data.coords.gene}`, omim: `Gene: ${data.coords.gene}`, gr: `Gene: ${data.coords.gene}`, scholar: `Gene: ${data.coords.gene}`, decipher: `Gene: ${data.coords.gene}`, hgmd: `Gene: ${data.coords.gene}`, gtex: `Gene: ${data.coords.gene}`, gnomadv2: `ID: ${data.coords.hg19String?.replace('chr', '') || ''}`, mastermind: mmTip, liftover: `Coords: ${data.coords.hg38String}` };
  for (const k in tips) { const el = document.getElementById('tip_' + k); if (el) el.innerText = tips[k]; }
}

// ── HPO PILLS (matched/unmatched term display) ─────────────────────────────
/**
 * Renders matched and unmatched HPO terms below the hpoInput.
 * Matched terms show: checkbox + "Name (HP:0000123)" in teal with ✓.
 * Unmatched terms show: "invalid term" in red/orange background (not selectable).
 * Selected terms are tracked in data.selectedHpoTerms for Scholar search.
 */
function renderHpoPills() {
  const container = document.getElementById('hpoPills');
  if (!container) return;

  const hpoInput = document.getElementById('hpoInput');
  if (!hpoInput || !hpoInput.value.trim()) {
    container.innerHTML = '';
    if (data && typeof data === 'object') data.selectedHpoTerms = {};
    return;
  }

  const inputTerms = hpoInput.value.split(',').map(s => s.trim()).filter(s => s && s !== '.' && s !== '-');
  const matched = []; // { term, ids }
  const unmatched = [];

  for (const term of inputTerms) {
    // phenotypeCache is populated immediately on autocomplete selection and after resolveOnePhenotype.
    // Using it directly avoids any index-alignment race with ptPhenoGroups/ptPhenoTexts.
    const ids = (typeof phenotypeCache !== 'undefined' && phenotypeCache[term]) || [];
    if (ids.length > 0) {
      matched.push({ term, ids });
    } else {
      unmatched.push(term);
    }
  }

  // Initialize selectedHpoTerms if not already present
  if (!data.selectedHpoTerms) data.selectedHpoTerms = {};

  let html = '';

  // Matched terms with checkboxes
  for (const { term, ids } of matched) {
    const idDisplay = ids.join(', ');
    const isChecked = data.selectedHpoTerms[term] ? 'checked' : '';
    const checkId = `hpoCheck_${escapeHtml(term).replace(/[^a-z0-9]/gi, '_')}`;
    html += `<div style="display:inline-block; margin:4px 4px 4px 0; padding:6px 10px; background:rgba(0,255,200,0.1); border:1px solid var(--teal); border-radius:4px; font-size:0.75rem; color:var(--teal); cursor:pointer; user-select:none;" onclick="toggleHpoSelection('${term.replace(/'/g, "\\'")}')">
      <input type="checkbox" id="${checkId}" ${isChecked} style="margin-right:4px; cursor:pointer;" onclick="event.stopPropagation(); toggleHpoSelection('${term.replace(/'/g, "\\'")}')">
      <span style="color:var(--teal); font-weight:500;">${escapeHtml(term)}</span>
      <span style="color:var(--dim); margin-left:4px;">${idDisplay}</span>
      <span style="margin-left:6px; color:var(--teal);">✓</span>
    </div>`;
  }

  // Unmatched terms (not selectable)
  for (const term of unmatched) {
    html += `<div style="display:inline-block; margin:4px 4px 4px 0; padding:6px 10px; background:rgba(255,100,100,0.15); border:1px solid rgba(255,100,100,0.4); border-radius:4px; font-size:0.75rem; color:rgba(255,100,100,0.8); font-style:italic;">
      ${escapeHtml(term)} <span style="color:rgba(255,100,100,0.8);">✗</span>
    </div>`;
  }

  container.innerHTML = html || '';
}

/**
 * Toggle HPO term selection for Scholar search.
 * Called when user clicks a matched term's checkbox.
 */
function toggleHpoSelection(term) {
  if (!data.selectedHpoTerms) data.selectedHpoTerms = {};
  data.selectedHpoTerms[term] = !data.selectedHpoTerms[term];
  renderHpoPills(); // Re-render to update checkbox state
}

// ── TWO-SOURCE VALIDATION MARKERS ──────────────────────────────────────────
// Maps validation field names to their DOM element IDs.
const VALIDATION_FIELD_TO_ELID = {
  gene:           'dGene',
  transcript:     'eManeTranscript',
  hgvsC:          'dHgvs',
  hgvsP:          'dProtein',
  hg38:           'dHg38',
  hg19:           'dHg19',
  exonIntron:     'eExon',
  rsId:           'eRsId',
  caId:           'eCAId',
  vepConsequence: 'eConsequence'
};

function _markerHtml(state, tooltip) {
  if (state === 'validated') {
    return ` <span style="color:var(--teal); font-size:0.75rem; vertical-align:middle;" title="${tooltip}">✅</span>`;
  }
  if (state === 'conflict') {
    return ` <span style="color:var(--amber); font-size:0.75rem; vertical-align:middle;" title="${tooltip}">⚠</span>`;
  }
  if (state === 'single') {
    return ` <span style="color:var(--dim); font-size:0.7rem; vertical-align:middle;" title="${tooltip}">ℹ️</span>`;
  }
  return '';
}

function _stripMarkers(s) {
  return String(s || '').replace(/\s*<span[^>]*>\s*[✅⚠ℹ️]+\s*<\/span>/g, '').trim();
}

/**
 * Reads data.validation.* and updates each annotation field's DOM with the
 * appropriate marker (✅ validated, ⚠ conflict, ℹ️ single-source).
 * Module isolation: only touches DOM IDs in VALIDATION_FIELD_TO_ELID, reads only data.validation.
 */
function renderValidationMarkers() {
  if (!data.validation) return;
  for (const [field, elId] of Object.entries(VALIDATION_FIELD_TO_ELID)) {
    const v = data.validation[field];
    const el = document.getElementById(elId);
    if (!v || !el) continue;

    let state = null;
    let tooltip = '';

    if (v.validated) {
      state = 'validated';
      const agreeingSrcs = Object.entries(v.sources || {})
        .filter(([, val]) => val === v.value)
        .map(([s]) => s);
      tooltip = `Validated by ${agreeingSrcs.join(', ')} = ${v.value}`;
      if (v.conflict && v.conflict.length) {
        tooltip += ' • Conflict: ' + v.conflict
          .map(c => `${c.srcs.join('+')} returned ${c.value}`).join('; ');
      }
    } else if (v.conflict) {
      state = 'conflict';
      tooltip = 'Source disagreement: ' + v.conflict.map(c => `${c.srcs.join('+')}=${c.value}`).join(' vs ');
    } else if (v.singleSource && v.value) {
      state = 'single';
      tooltip = `Single-source field (${v.available[0] || 'API'} only)`;
    }

    // Display value if available, else leave existing text. Append marker.
    if (v.value != null) {
      // Preserve any existing inline structure (e.g., links) — replace marker only
      const current = el.innerHTML;
      const stripped = _stripMarkers(current);
      // If field's current display matches the validated value, only update marker.
      // Otherwise update display to the validated value (canonical form).
      const valueDisplay = String(v.value);
      if (stripped === '—' || stripped === '' || stripped === '-') {
        el.innerHTML = valueDisplay + _markerHtml(state, tooltip);
      } else {
        // Keep existing text (may be richer formatting), just refresh the marker
        el.innerHTML = stripped + _markerHtml(state, tooltip);
      }
    } else if (state === null) {
      // No data — make sure no stale marker hangs around
      const stripped = _stripMarkers(el.innerHTML);
      if (stripped !== el.innerHTML) el.innerHTML = stripped;
    }
  }
}

/**
 * Convenience: collect sources for all fields, validate them, and re-render.
 * Call after any API completes (or fails) to update markers progressively.
 */
function runValidationPass() {
  if (typeof window.Validation === 'undefined') return;
  window.Validation.validateAll(data);
  renderValidationMarkers();
}

/**
 * Renders all metric-related panels (Gene Pathogenicity, Pop Frequency, In Silico)
 * Centralizes logic previously scattered in API callbacks.
 */
function renderMetricPanels() {
  const cgLink = document.getElementById('clingenSearchLink');
  if (cgLink) {
    cgLink.style.display = 'none';
  }

  // --- Panel 1: Gene Pathogenicity ---
  const gCon = document.getElementById('gnomadGeneConstraint');
  if (data.gnomad.geneConstraint) {
    const { pli, loeuf, mis_z } = data.gnomad.geneConstraint;
    const isPliHigh = parseFloat(pli) > 0.9;
    const isLoeufLow = parseFloat(loeuf) < 0.35;
    const isMisZHigh = parseFloat(mis_z) > 3.09;

    // Scale calculations
    const pliVal = parseFloat(pli);
    const pliPct = Math.min(100, Math.max(0, pliVal * 100)); // pLI scale 0 to 1
    const loeufVal = parseFloat(loeuf);
    const loeufPct = Math.min(100, Math.max(0, (loeufVal / 2.0) * 100)); // LOEUF scale 0 to 2
    const misZVal = Math.min(10, parseFloat(mis_z)); // Cap display at 10 for bar
    const misZPct = Math.min(100, Math.max(0, (misZVal / 10.0) * 100)); // MisZ scale 0 to 10

    gCon.innerHTML = `
      <div style="width:100%; display: flex; flex-direction: column; gap: 8px;">
        <!-- pLI Row Group -->
        <div style="display: flex; flex-direction: column;">
          <div style="display: flex; align-items: center; gap: 12px;">
            <div style="min-width: 92px; font-size: 0.82rem; font-weight: 600; color: var(--dim); text-align: right; line-height: 1;">
              pLI: <b class="${isPliHigh ? 'bad-val' : ''}" style="font-size: 0.82rem; margin-left: 2px;">${pli}</b>
            </div>
            <div class="premium-bar-container" data-tippy-content="pLI: ${pli} (Cutoff: 0.9)" style="flex: 1; margin: 0;">
              <div class="premium-bar-tick" style="left: 90%"></div> <!-- 0.9 / 1.0 = 90% -->
              <div class="premium-bar-fill" style="width: ${pliPct}%; background: ${isPliHigh ? '#ff5577' : 'var(--teal)'}"></div>
            </div>
          </div>
          <div style="display: flex; gap: 12px;">
            <div style="min-width: 92px;"></div>
            <div class="premium-scale-indicators" style="flex: 1;">
              <span class="scale-label" style="left: 0%">0</span>
              <span class="scale-label" style="left: 100%">1</span>
            </div>
          </div>
        </div>

        <!-- LOEUF Row Group -->
        <div style="display: flex; flex-direction: column;">
          <div style="display: flex; align-items: center; gap: 12px;">
            <div style="min-width: 92px; font-size: 0.82rem; font-weight: 600; color: var(--dim); text-align: right; line-height: 1;">
              LOEUF: <b class="${isLoeufLow ? 'bad-val' : ''}" style="font-size: 0.82rem; margin-left: 2px;">${loeuf}</b>
            </div>
            <div class="premium-bar-container" data-tippy-content="LOEUF: ${loeuf} (Cutoff: 0.35)" style="flex: 1; margin: 0;">
              <div class="premium-bar-tick" style="left: 17.5%"></div> <!-- 0.35 / 2.0 = 17.5% -->
              <div class="premium-bar-fill" style="width: ${loeufPct}%; background: ${isLoeufLow ? '#ff5577' : 'var(--teal)'}"></div>
            </div>
          </div>
          <div style="display: flex; gap: 12px;">
            <div style="min-width: 92px;"></div>
            <div class="premium-scale-indicators" style="flex: 1;">
              <span class="scale-label" style="left: 0%">0</span>
              <span class="scale-label" style="left: 50%">1</span>
              <span class="scale-label" style="left: 100%">2</span>
            </div>
          </div>
        </div>

        <!-- Missense Z Row Group -->
        <div style="display: flex; flex-direction: column; margin-top: 4px;">
          <div style="display: flex; align-items: center; gap: 12px;">
            <div style="min-width: 92px; font-size: 0.82rem; font-weight: 600; color: var(--dim); text-align: right; line-height: 1;">
              Mis Z: <b class="${isMisZHigh ? 'bad-val' : ''}" style="font-size: 0.82rem; margin-left: 2px;">${mis_z}</b>
            </div>
            <div class="premium-bar-container" data-tippy-content="Missense Z: ${mis_z} (Cutoff: 3.09)" style="flex: 1; margin: 0;">
              <div class="premium-bar-tick" style="left: 30.9%"></div> <!-- 3.09 / 10.0 = 30.9% -->
              <div class="premium-bar-fill" style="width: ${misZPct}%; background: ${isMisZHigh ? '#ff5577' : 'var(--teal)'}"></div>
            </div>
          </div>
          <div style="display: flex; gap: 12px;">
            <div style="min-width: 92px;"></div>
            <div class="premium-scale-indicators" style="flex: 1;">
              <span class="scale-label" style="left: 0%">0</span>
              <span class="scale-label" style="left: 100%">10 <small style="margin-left:2px">${isMisZHigh ? '↑' : ''}</small></span>
            </div>
          </div>
        </div>
      </div>`;
  } else if (data.coords.gene) {
    gCon.innerText = 'Not found';
  }

  // Regional Missense (hg19) is variant-position-specific: it maps a single variant
  // onto a gnomAD constrained interval. In gene-only mode there is no variant, so the
  // row carries no meaning — hide it entirely to avoid a misleading empty "—" row.
  const regionalRow = document.getElementById('rowRegional');
  if (regionalRow) regionalRow.style.display = data.coords.hg38String ? '' : 'none';

  const eRegOE = document.getElementById('eRegionalOE');
  const eRegRange = document.getElementById('eRegionalRange');
  const eRegPVal = document.getElementById('eRegionalPValue');
  const lReg = document.getElementById('labelRegional');
  const fmtPVal = p => {
    if (p === null || p === undefined) return '';
    if (p >= 0.001) return `p = ${p.toFixed(3)}`;
    const s = p.toExponential(2);
    const [m, e] = s.split('e');
    const sup = '⁰¹²³⁴⁵⁶⁷⁸⁹';
    const expStr = e.replace('-', '⁻').replace(/\d/g, d => sup[+d]);
    return `p = ${m}×10${expStr}`;
  };
  if (data.gnomad.regionalOE !== null) {
    const oe = parseFloat(data.gnomad.regionalOE);
    eRegOE.innerText = oe.toFixed(3);
    eRegRange.innerText = `Interval: ${data.gnomad.regionalRange || '—'}`;
    if (eRegPVal) eRegPVal.innerText = fmtPVal(data.gnomad.regionalPValue);
    if (oe < 0.6) {
      eRegOE.className = 'dv bad-val';
      if (lReg && !lReg.innerText.includes('⚠️')) lReg.innerHTML = '⚠️ Regional Missense (hg19)';
    } else if (oe < 0.8) {
      eRegOE.className = 'dv';
      eRegOE.style.color = 'var(--amber)';
      if (lReg) lReg.innerHTML = 'Regional Missense (hg19)';
    } else {
      eRegOE.className = 'dv';
      eRegOE.style.color = 'var(--text-bright)';
      if (lReg) lReg.innerHTML = 'Regional Missense (hg19)';
    }
  } else if (data.gnomad.regionalRange) {
    if (eRegOE) { eRegOE.innerText = '—'; eRegOE.className = 'dv'; }
    if (eRegRange) eRegRange.innerText = `Interval: ${data.gnomad.regionalRange}`;
    if (eRegPVal) eRegPVal.innerText = '';
    if (lReg) lReg.innerHTML = 'Regional Missense (hg19)';
  } else if (data.coords.gene) {
    if (eRegOE) { eRegOE.innerText = '—'; eRegOE.className = 'dv'; }
    if (eRegRange) eRegRange.innerText = 'Interval: —';
    if (eRegPVal) eRegPVal.innerText = '';
    if (lReg) lReg.innerHTML = 'Regional Missense (hg19)';
  }

  // --- ClinGen Dosage Sensitivity ---
  const clingenMap = {
    "3": { label: "3 (Sufficient)", color: "var(--red)", weight: "800", tip: "Sufficient evidence for disease mechanism. PVS1 may be applicable." },
    "2": { label: "2 (Emerging)", color: "var(--amber)", weight: "400", tip: "Emerging evidence. Use caution." },
    "1": { label: "1 (Little Evidence)", color: "var(--text-bright)", weight: "400", tip: "Little evidence." },
    "0": { label: "0 (No Evidence)", color: "var(--text-dim)", weight: "400", tip: "No evidence." },
    "30": { label: "30 (Unlikely)", color: "var(--text-dim)", weight: "400", tip: "Dosage sensitivity unlikely." },
    "4": { label: "4 (Triplo-sensitive Unlikely)", color: "var(--text-dim)", weight: "400", tip: "Triplosensitivity unlikely." }
  };

  const renderDosage = (score, elId) => {
    const el = document.getElementById(elId);
    if (!el) return;
    if (score === null || score === undefined) { el.innerText = '—'; el.style.color = 'var(--text-dim)'; return; }
    const s = String(score).split('/')[0]; // Handle "3/3" or "3"
    const cfg = clingenMap[s];
    if (cfg) {
      el.innerText = cfg.label;
      el.style.color = cfg.color;
      el.style.fontWeight = cfg.weight;
      if (s === "3") el.title = cfg.tip;
    } else {
      el.innerText = score;
      el.style.color = 'var(--text-bright)';
    }
  };

  renderDosage(data.dosage.haplo, 'eHaploScore');
  renderDosage(data.dosage.triplo, 'eTriploScore');

  const dosageSrcEl = document.getElementById('eDosageSource');
  if (dosageSrcEl) {
    if (data.dosage.source === 'ucsc') {
      dosageSrcEl.innerHTML = 'Source: <span style="color:var(--teal)">UCSC Live</span>';
    } else if (data.dosage.source === 'local') {
      dosageSrcEl.innerHTML = 'Source: <span style="color:var(--amber)">Local</span>';
    } else {
      dosageSrcEl.innerText = '';
    }
  }

  // --- Panel 2: Population Frequency ---
  renderPopulationData(data.gnomadViewMode || 'total');

  // --- Panel 3: In Silico Predictions ---
  const rev = document.getElementById('mvRevel');
  if (data.scores.revel !== null) {
    // [FIX] Parse as float to prevent .toFixed crashes during rendering
    const r = parseFloat(data.scores.revel);
    let status = 'Intermediate (No Code)', color = 'var(--dim)', cut = '0.290–0.644';
    if (r >= 0.932) { status = 'PP3_S — Strong Pathogenic'; color = '#ff3366'; cut = '≥ 0.932'; }
    else if (r >= 0.773) { status = 'PP3_M — Moderate Pathogenic'; color = 'var(--red)'; cut = '0.773–0.932'; }
    else if (r >= 0.644) { status = 'PP3 — Supporting Pathogenic'; color = 'var(--amber)'; cut = '0.644–0.773'; }
    else if (r <= 0.003 && r > 0) { status = 'BP4_VS — Very Strong Benign'; color = 'var(--teal)'; cut = '≤ 0.003'; }
    else if (r <= 0.016) { status = 'BP4_S — Strong Benign'; color = 'var(--teal)'; cut = '0.003–0.016'; }
    else if (r <= 0.183) { status = 'BP4_M — Moderate Benign'; color = '#66ccaa'; cut = '0.016–0.183'; }
    else if (r <= 0.290) { status = 'BP4 — Supporting Benign'; color = 'var(--dim)'; cut = '0.183–0.290'; }
    else if (r === 0) { status = 'BP4_S — Strong Benign (≡0)'; color = 'var(--teal)'; cut = '= 0'; }

    rev.innerHTML = `
      <div style="display:flex; flex-direction:column; gap:6px;">
        <div style="display:flex; justify-content:space-between; align-items:baseline;">
          <span style="font-weight:700; color:${color}; font-size:0.95rem;">${r.toFixed(3)}</span>
          <span style="font-size:0.6rem; color:${color}; opacity:0.85; text-align:right; max-width:130px; line-height:1.3;">${status}</span>
        </div>
        <div class="bullet-graph-wrapper" data-tippy-content="REVEL: ${r}">
          <div class="bullet-container">
            <!-- Score Bar -->
            <div class="bullet-bar" style="width: ${Math.min(100, Math.max(0, r * 100))}%; background: ${color}"></div>
          </div>
        </div>
        <div style="font-size:0.58rem; color:var(--dim); text-align:right;">Range: ${cut}</div>
      </div>`;

    rev.className = 'dv';
  } else {
    rev.innerText = 'Not Available';
    rev.className = 'dv';
  }

  // AlphaMissense UI Rendering
  const amEl = document.getElementById('mvAlphaMissense');
  if (amEl) {
    if (data.scores.alphaMissenseScore !== null) {
      const am = parseFloat(data.scores.alphaMissenseScore);
      let status = 'Ambiguous', color = 'var(--dim)', cut = '0.331–0.761';
      
      if (am >= 0.984) { status = 'PP3_S — Strong Pathogenic'; color = '#ff3366'; cut = '≥ 0.984'; }
      else if (am >= 0.869) { status = 'PP3_M — Moderate Pathogenic'; color = 'var(--red)'; cut = '0.869–0.984'; }
      else if (am >= 0.761) { status = 'PP3 — Supporting Pathogenic'; color = 'var(--amber)'; cut = '0.761–0.869'; }
      else if (am <= 0.073) { status = 'BP4_VS — Very Strong Benign'; color = 'var(--teal)'; cut = '≤ 0.073'; }
      else if (am <= 0.147) { status = 'BP4_S — Strong Benign'; color = 'var(--teal)'; cut = '0.073–0.147'; }
      else if (am <= 0.331) { status = 'BP4_M — Moderate Benign'; color = '#66ccaa'; cut = '0.147–0.331'; }


      const predLabel = data.scores.alphaMissensePred ? ` (${data.scores.alphaMissensePred})` : '';

      amEl.innerHTML = `
        <div style="display:flex; flex-direction:column; gap:6px;">
          <div style="display:flex; justify-content:space-between; align-items:baseline;">
            <span style="font-weight:700; color:${color}; font-size:0.95rem;">${am.toFixed(3)}${predLabel}</span>
            <span style="font-size:0.6rem; color:${color}; opacity:0.85; text-align:right; max-width:130px; line-height:1.3;">${status}</span>
          </div>
          <div class="bullet-graph-wrapper" data-tippy-content="AlphaMissense: ${am}">
            <div class="bullet-container">
              <div class="bullet-bar" style="width: ${Math.min(100, Math.max(0, am * 100))}%; background: ${color}"></div>
            </div>
          </div>
          <div style="font-size:0.58rem; color:var(--dim); text-align:right;">Range: ${cut}</div>
        </div>`;
      amEl.className = 'dv';
    } else {
      amEl.innerText = 'Not Available';
      amEl.className = 'dv';
    }
  }

  const spAI = document.getElementById('mvSpliceAI');
  if (data.scores.spliceAI !== null) {
    spAI.innerText = data.scores.spliceAI.toFixed(3);
    spAI.className = 'dv' + (data.scores.spliceAI > 0.20 ? ' bad' : ' good');

    const updateScore = (id, val, pos) => {
      const el = document.getElementById(id);
      if (el) {
        if (val !== null && val !== undefined) {
          const posStr = (pos !== null && pos !== undefined) ? `<span style="font-size:0.7rem; color:var(--dim); font-weight:normal; margin-left:4px;">(${pos > 0 ? '+' + pos : pos})</span>` : '';
          el.innerHTML = `${val.toFixed(3)}${posStr}`;
          el.style.color = val > 0.20 ? 'var(--red)' : 'var(--green)';
        } else {
          el.innerHTML = '—';
          el.style.color = 'inherit';
        }
      }
    };
    updateScore('mvSpliceAI_AG', data.scores.spliceAI_AG, data.scores.spliceAI_DP_AG);
    updateScore('mvSpliceAI_AL', data.scores.spliceAI_AL, data.scores.spliceAI_DP_AL);
    updateScore('mvSpliceAI_DG', data.scores.spliceAI_DG, data.scores.spliceAI_DP_DG);
    updateScore('mvSpliceAI_DL', data.scores.spliceAI_DL, data.scores.spliceAI_DP_DL);

    const consText = getSpliceConsequence();
    const consEl = document.getElementById('spliceAiConsequence');
    if (consEl) consEl.innerText = consText || '';
  }

  // mvSuggest update moved to evaluateACMG for consolidated logic
}


/**
 * Renders a horizontal bullet graph for clinical thresholds.
 */
function renderBulletGraph(type, value, threshold, max, direction = 'greater') {
  const val = parseFloat(value);
  const thresh = parseFloat(threshold);
  const m = parseFloat(max);
  if (isNaN(val)) return '';

  const pct = Math.min(100, Math.max(0, (val / m) * 100));
  const threshPct = (thresh / m) * 100;

  // Directional Breach Zone Logic
  const isBreached = direction === 'greater' ? (val >= thresh) : (val <= thresh);
  const barColor = isBreached ? 'var(--red)' : 'var(--teal)';

  return `
    <div class="bullet-graph-wrapper" data-tippy-content="${type.toUpperCase()}: ${val} (Cutoff: ${thresh})">
      <div class="bullet-container">
        <div class="bullet-marker" style="left: ${threshPct}%"></div>
        <div class="bullet-bar" style="width: ${pct}%; background: ${barColor}"></div>
      </div>
    </div>`;
}

/**
 * Toggles the GenCC expanded details panel for a specific gene.
 * Exposed on window so inline onclick handlers in generated HTML can reach it.
 */
window.toggleGenCCDetails = function (geneSymbol) {
  const detailsDiv = document.getElementById(`gencc-details-${geneSymbol}`);
  const toggleIcon = document.getElementById(`gencc-icon-${geneSymbol}`);
  if (detailsDiv.style.display === 'none') {
    detailsDiv.style.display = 'block';
    if (toggleIcon) toggleIcon.innerText = '⊖';
  } else {
    detailsDiv.style.display = 'none';
    if (toggleIcon) toggleIcon.innerText = '⊕';
  }
};

/**
 * Renders the GenCC-style expandable panel for associated conditions, grouped by validity classification.
 */
function renderConditionCards(conditions, geneSymbol, validityCurations = []) {
  const wrapper = document.getElementById('associated-conditions-wrapper');
  const container = document.getElementById('associated-conditions-container');
  const header = document.getElementById('conditions-toggle-header');
  const countText = document.getElementById('conditions-count-text');

  if (header) header.style.display = 'flex';

  if (!conditions || conditions.length === 0) {
    if (countText) countText.innerHTML = `Gene-Disease Associations for ${geneSymbol} <span style="color:var(--amber); font-size:0.75rem; margin-left:8px;">(No associations found)</span>`;
    return;
  }

  if (countText) countText.innerText = `Gene-Disease Associations for ${geneSymbol}`;

  // --- 1. DATA PROCESSING & GENCC BUCKETING ---
  const uniqueDiseases = new Set();
  const uniqueSubmitters = new Set();

  const evidenceOrder = ['Definitive', 'Strong', 'Moderate', 'Supportive', 'Limited', 'Disputed', 'Refuted', 'Animal', 'No Known'];

  // Canonical name mapping: AU name takes priority when AU and UK refer to the same disease.
  // Matching uses two signals: (1) all significant words of the shorter name appear in the longer,
  // (2) a substantive disorder string from one panel is a substring of one from the other.
  const _sigWords = name => name.toLowerCase().replace(/[^a-z0-9\s]/g, ' ')
    .split(/\s+/).filter(w => w.length >= 5);
  const _isSameDisease = (a, b) => {
    const wa = _sigWords(a.diseaseName), wb = _sigWords(b.diseaseName);
    if (wa.length > 0 && wb.length > 0) {
      const [shorter, longer] = wa.length <= wb.length ? [wa, wb] : [wb, wa];
      if (shorter.every(w => longer.includes(w))) return true;
    }
    const cleanDisorders = arr => (arr || [])
      .filter(s => !/HP:/i.test(s) && !/^\w\d{2,}$/.test(s.trim()))  // drop HP codes and NHS R-codes
      .map(s => s.toLowerCase().replace(/[^a-z0-9\s]/g, ' ').trim())
      .filter(d => d.length >= 8);                                     // require substantive strings
    const aD = cleanDisorders(a._disorders);
    const bD = cleanDisorders(b._disorders);
    return aD.some(da => bD.some(db => da.includes(db) || db.includes(da)));
  };
  const auConds = conditions.filter(c => (c.sources || '').startsWith('PanelApp Australia'));
  const ukConds = conditions.filter(c => (c.sources || '').startsWith('PanelApp UK'));
  const ukToAuName = new Map();
  for (const auC of auConds) {
    for (const ukC of ukConds) {
      if (!ukToAuName.has(ukC.diseaseName) && _isSameDisease(auC, ukC)) {
        ukToAuName.set(ukC.diseaseName, auC.diseaseName);
      }
    }
  }

  // Merge conditions with the same (canonical) disease name, combining sources.
  // PanelApp AU has highest priority: if AU evidence exists, use it; otherwise keep best evidence level.
  const mergedMap = new Map();
  conditions.forEach(item => {
    const key = ukToAuName.get(item.diseaseName) || item.diseaseName;
    const isAU = (item.sources || '').startsWith('PanelApp Australia');
    if (mergedMap.has(key)) {
      const existing = mergedMap.get(key);
      if (item.sources && !existing.allSources.includes(item.sources)) {
        existing.allSources.push(item.sources);
      }
      // Track all panel links; AU links go first
      if (item.diseaseId && !existing._links.some(l => l.id === item.diseaseId)) {
        if (isAU) existing._links.unshift({ id: item.diseaseId, url: item.externalUrl || '' });
        else        existing._links.push({ id: item.diseaseId, url: item.externalUrl || '' });
      }
      // AU always wins on evidence level; otherwise keep the best (lowest index)
      const existingIsAU = existing.allSources.some(s => s.startsWith('PanelApp Australia'));
      if (isAU) {
        existing.evidenceLevel = item.evidenceLevel; // AU overrides
      } else if (!existingIsAU) {
        const existingIdx = evidenceOrder.findIndex(e => existing.evidenceLevel.toLowerCase().includes(e.toLowerCase()));
        const newIdx = evidenceOrder.findIndex(e => (item.evidenceLevel || '').toLowerCase().includes(e.toLowerCase()));
        if (newIdx !== -1 && (existingIdx === -1 || newIdx < existingIdx)) {
          existing.evidenceLevel = item.evidenceLevel;
        }
      }
    } else {
      const entry = { ...item, diseaseName: key, allSources: item.sources ? [item.sources] : [] };
      entry._links = item.diseaseId ? [{ id: item.diseaseId, url: item.externalUrl || '' }] : [];
      mergedMap.set(key, entry);
    }
  });
  const mergedConditions = Array.from(mergedMap.values());

  // ── ClinGen Gene-Disease Validity reconciliation (authoritative override) ──
  // GenCC sources (incl. Monarch's ClinGen signal) carry only a heuristic evidence
  // level — Monarch hard-codes "Strong" for every ClinGen-sourced disease. The REAL
  // classification lives in data.geneValidity (passed in explicitly — Rule 2). Match
  // by ontology id (MONDO) first, then disease name, and override the level so the
  // bucket + pill reflect the authoritative ClinGen call.
  const VALID_LEVELS = ['Definitive', 'Strong', 'Moderate', 'Limited', 'Disputed', 'Refuted', 'Animal', 'No Known'];
  if (Array.isArray(validityCurations) && validityCurations.length) {
    mergedConditions.forEach(cond => {
      const ids = (cond._links || []).map(l => l.id).concat(cond.diseaseId ? [cond.diseaseId] : []);
      const cu = validityCurations.find(c =>
        (c.mondo && ids.some(id => sameDiseaseId(c.mondo, id))) ||
        sameDiseaseName(c.disease, cond.diseaseName)
      );
      if (cu && cu.classification && VALID_LEVELS.some(l => cu.classification.toLowerCase().includes(l.toLowerCase()))) {
        cond.evidenceLevel = cu.classification;       // display-only override (cond is a merged copy)
        cond._clingenValidity = cu;                   // carries gcep / date / moi for the badge tooltip
      }
    });
  }

  // Standard GenCC Categories (Strict, No "Other")
  const buckets = {
    'Definitive': [], 'Strong': [], 'Moderate': [], 'Supportive': [],
    'Limited': [], 'Disputed': [], 'Refuted': [], 'Animal': [], 'No Known': []
  };

  mergedConditions.forEach(item => {
    uniqueDiseases.add(item.diseaseName);
    item.allSources.forEach(s => uniqueSubmitters.add(s));

    const ev = (item.evidenceLevel || '').toLowerCase();
    let assigned = false;
    for (const key of Object.keys(buckets)) {
      if (ev.includes(key.toLowerCase())) {
        buckets[key].push(item);
        assigned = true;
        break;
      }
    }
    if (!assigned) buckets['No Known'].push(item);
  });

  // Within each bucket, AU-sourced diseases first, then others.
  const hasAU = item => (item.allSources || []).some(s => s.startsWith('PanelApp Australia'));
  for (const bucket of Object.values(buckets)) {
    bucket.sort((a, b) => (hasAU(b) ? 1 : 0) - (hasAU(a) ? 1 : 0));
  }

  // --- 2. BUILD UI: SUMMARY HEADER (COLLAPSED VIEW) ---
  // Official GenCC UI Colors
  const pillColors = {
    'Definitive': '#1b5e20', // Dark Green
    'Strong': '#2e7d32', // Green
    'Moderate': '#66bb6a', // Light Green
    'Supportive': '#29b6f6', // Blue
    'Limited': '#ff8a80', // Light Red
    'Disputed': '#f44336', // Red
    'Refuted': '#b71c1c', // Dark Red
    'Animal': '#ff9800', // Orange
    'No Known': '#90a4ae'  // Grey
  };

  let summaryPillsHtml = Object.keys(buckets).map(key => {
    const count = buckets[key].length;
    const isActive = count > 0;
    const bgColor = isActive ? pillColors[key] : 'rgba(255,255,255,0.06)';
    const textColor = isActive ? '#fff' : 'rgba(255,255,255,0.2)';
    const abbr = {
      'Definitive': 'Def', 'Strong': 'Str', 'Moderate': 'Mod', 'Supportive': 'Sup',
      'Limited': 'Lim', 'Disputed': 'Dis', 'Refuted': 'Ref', 'Animal': 'Ani', 'No Known': 'N/K'
    }[key] || key;
    return `
      <div style="display:flex; flex-direction:column; align-items:center; gap:3px; min-width:34px;" data-tippy-content="${key}: ${count}">
        <div style="background:${bgColor}; color:${textColor}; border-radius:10px; padding:2px 6px;
             font-size:0.7rem; font-weight:800; min-width:22px; text-align:center; line-height:1.4;">${count}</div>
        <div style="font-size:0.55rem; color:${isActive ? 'var(--dim)' : 'rgba(255,255,255,0.15)'};
             text-transform:uppercase; letter-spacing:0.02em; white-space:nowrap;">${abbr}</div>
      </div>
    `;
  }).join('');

  let html = `
    <div style="background:var(--bg); border:1px solid var(--border); border-radius:var(--r); margin-bottom:16px;">

      <div style="display:flex; align-items:center; justify-content:space-between; padding:16px 20px; cursor:pointer;" onclick="toggleGenCCDetails('${geneSymbol}')">
        <div style="display:flex; gap:32px; align-items:center;">
          <div style="min-width:80px;">
            <div style="font-weight:800; font-size:1.1rem; color:var(--text-bright);">${geneSymbol}</div>
          </div>
          <div>
            <div style="font-weight:800; font-size:1rem; color:var(--text-bright);">${uniqueDiseases.size}</div>
            <div style="font-size:0.7rem; color:var(--dim); text-transform:uppercase;">Disease Equivalents</div>
          </div>
          <div>
            <div style="font-weight:800; font-size:1rem; color:var(--text-bright);">${uniqueSubmitters.size}</div>
            <div style="font-size:0.7rem; color:var(--dim); text-transform:uppercase;">Submitters</div>
          </div>
        </div>

        <div style="display:flex; align-items:flex-end; gap:8px;">
          <div style="display:flex; gap:4px; margin-right:20px;">
            ${summaryPillsHtml}
          </div>
          <div style="color:var(--dim); font-size:0.85rem; font-weight:600; display:flex; align-items:center; gap:6px;">
            Details <span id="gencc-icon-${geneSymbol}" style="font-size:1.1rem;">⊕</span>
          </div>
        </div>
      </div>

      <div id="gencc-details-${geneSymbol}" style="display:none; border-top:1px solid var(--border); padding:20px; max-height:500px; overflow-y:auto; scrollbar-width:thin;">
  `;

  // --- 3. BUILD UI: EXPANDED DETAILS (GROUPED BY CLASSIFICATION) ---
  Object.keys(buckets).forEach(key => {
    if (buckets[key].length === 0) return;

    html += `<h4 style="margin:0 0 12px 0; font-size:1rem; color:var(--text-bright); padding-bottom:4px; border-bottom:1px solid rgba(255,255,255,0.05);">${key} classifications</h4>`;
    html += `<div style="display:flex; flex-direction:column; gap:8px; margin-bottom:24px; max-height:400px; overflow-y:auto; padding-right:8px; scrollbar-width:thin;">`;

    buckets[key].forEach(item => {
      const allLinks = item._links || (item.diseaseId ? [{ id: item.diseaseId, url: item.externalUrl || '' }] : []);
      const idLinksHtml = allLinks.map(l => l.url
        ? `<a href="${l.url}" target="_blank" style="color:var(--teal);text-decoration:none;font-size:0.68rem;">${l.id} ↗</a>`
        : `<span style="font-size:0.68rem;">${l.id}</span>`).join(' ');
      const bgColor = pillColors[key] || 'var(--dim)';

      // Build source badges — PanelApp AU first, then UK, then others
      const sourceOrder = ['PanelApp Australia', 'PanelApp UK', 'ClinGen', 'OMIM'];
      const sortedSources = (item.allSources || [item.sources || '']).slice().sort((a, b) => {
        const ai = sourceOrder.findIndex(s => a.startsWith(s));
        const bi = sourceOrder.findIndex(s => b.startsWith(s));
        return (ai === -1 ? 99 : ai) - (bi === -1 ? 99 : bi);
      });
      const sourceBadges = sortedSources.map(src => {
        if (item._clingenValidity && src.startsWith('ClinGen')) return ''; // superseded by the dedicated validity badge below
        if (src.startsWith('PanelApp Australia')) return `<span style="background:#0277bd; color:#fff; padding:2px 8px; border-radius:6px; font-size:0.65rem; font-weight:700; white-space:nowrap; letter-spacing:0.03em;">PanelApp AU</span>`;
        if (src.startsWith('PanelApp UK'))        return `<span style="background:#283593; color:#fff; padding:2px 8px; border-radius:6px; font-size:0.65rem; font-weight:700; white-space:nowrap; letter-spacing:0.03em;">PanelApp UK</span>`;
        if (src.startsWith('ClinGen'))             return `<span style="background:#4a148c; color:#fff; padding:2px 8px; border-radius:6px; font-size:0.65rem; font-weight:700; white-space:nowrap; letter-spacing:0.03em;">ClinGen</span>`;
        if (src.startsWith('OMIM'))                return `<span style="background:#4e342e; color:#fff; padding:2px 8px; border-radius:6px; font-size:0.65rem; font-weight:700; white-space:nowrap; letter-spacing:0.03em;">OMIM</span>`;
        return `<span style="background:rgba(255,255,255,0.1); color:var(--dim); padding:2px 8px; border-radius:6px; font-size:0.65rem; font-weight:700; white-space:nowrap;">${src || 'Unknown'}</span>`;
      }).join(' ');

      // Authoritative ClinGen Gene-Disease Validity badge (replaces the generic ClinGen badge).
      const cv = item._clingenValidity;
      const clingenValidityBadge = cv
        ? `<span data-tippy-content="ClinGen Gene-Disease Validity${cv.gcep ? ' · ' + cv.gcep : ''}${cv.date ? ' · ' + cv.date : ''}${cv.moi ? ' · ' + (typeof expandClingenMoi === 'function' ? expandClingenMoi(cv.moi) : cv.moi) : ''}" style="background:#00897b; color:#fff; padding:2px 8px; border-radius:6px; font-size:0.65rem; font-weight:700; white-space:nowrap; letter-spacing:0.03em;">ClinGen Validity ✓</span>`
        : '';

      // Extract 6-digit OMIM MIM ID for synopsis badge placeholder
      const synMimId = allLinks.find(l => l.id?.startsWith('OMIM:'))?.id?.replace('OMIM:', '') ?? null;

      html += `
        <div style="display:flex; flex-direction:column; gap:5px; padding:10px 14px;
             background:rgba(0,0,0,0.15); border-radius:8px; border-left:3px solid ${bgColor};">
          <div style="display:flex; align-items:center; gap:8px; flex-wrap:wrap;">
            <span style="background:${bgColor}; color:#fff; padding:2px 9px; border-radius:8px;
                 font-weight:700; font-size:0.65rem; white-space:nowrap;">${key}</span>
            ${sourceBadges}${clingenValidityBadge}
            <span style="font-weight:700; color:var(--text-bright); font-size:0.82rem;">${item.diseaseName}</span>
          </div>
          <div style="display:flex; gap:14px; flex-wrap:wrap; font-size:0.72rem; color:var(--dim); padding-left:2px;">
            <span>🧬 <b style="color:var(--dim);">${geneSymbol}</b></span>
            <span>⚙️ ${item.inheritance || 'MOI unknown'}</span>
            ${idLinksHtml ? `<span>${idLinksHtml}</span>` : ''}
            ${synMimId ? `<span id="omim-syn-${synMimId}" style="display:none;"></span>` : ''}
          </div>
        </div>
      `;
    });
    html += `</div>`;
  });

  html += `
      </div>
    </div>
  `;

  if (container) {
    container.innerHTML = html;
    // If OMIM synopsis results already arrived (runPhenoMatch faster than fetchAssociatedConditions),
    // apply badges immediately now that the placeholders exist in the DOM.
    if (typeof data !== 'undefined' && data.omimSynopsisResults?.length) {
      applyOmimSynopsisBadges(data.omimSynopsisResults, data.ptPhenoTexts?.length ?? 0);
    }
  }
}

/**
 * Handles the success state for the copy button.
 */
function handleCopySuccess() {
  const btn = document.getElementById('btnCopyReport');
  if (!btn) return;
  const originalContent = btn.innerHTML;

  // Teal background success state
  btn.style.background = 'var(--teal)';
  btn.style.borderColor = 'var(--teal)';
  btn.innerHTML = `<span class="btn-icon" style="color: white; filter: none">✔</span> Copied!`;
  btn.classList.add('copy-success');

  setTimeout(() => {
    btn.style.background = '';
    btn.style.borderColor = '';
    btn.innerHTML = originalContent;
    btn.classList.remove('copy-success');
  }, 1500);
}

/**
 * Updates the LitVar 2.0 DOM container independently.
 */
function renderLitVarUI() {
  const container = document.getElementById('pLitVarMatches');
  if (!container) return;

  if (!data.litvarPMIDs || data.litvarPMIDs.length === 0) {
    container.innerHTML = 'None';
    return;
  }

  const pmids = data.litvarPMIDs;
  const count = pmids.length;
  const displayLimit = 3;
  const topPMIDs = pmids.slice(0, displayLimit);

  const linksHtml = topPMIDs.map(p => `<a href="https://pubmed.ncbi.nlm.nih.gov/${p}" target="_blank" style="color:var(--teal); text-decoration:none;">${p}</a>`).join(', ');

  let html = `<div style="display:flex; flex-direction:column; gap:2px;">`;
  html += `<div style="font-weight:800; font-size:0.85rem; color:var(--text-bright);">${count} Match${count !== 1 ? 'es' : ''}</div>`;
  html += `<div style="font-size:0.7rem; color:var(--dim);">${linksHtml}</div>`;

  if (count > displayLimit) {
    const query = data.rsId && data.rsId !== '—' ? data.rsId : data.coords.hgvs;
    const lvUrl = `https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/#!?query=${encodeURIComponent(query)}`;
    html += `<a href="${lvUrl}" target="_blank" style="font-size:0.65rem; color:var(--teal); text-decoration:none; margin-top:2px;">See all ${count} in LitVar ↗</a>`;
  }

  html += `</div>`;
  container.innerHTML = html;
}

/**
 * renderLitMultiTxUI(total, notations, pmUrl, gene)
 * Inputs:  data.litvar.combinedTotal, data.litvar.uniqueNotations (passed explicitly)
 * Outputs: #rowPubMedMultiTx, #pMultiTxTotal, #pMultiTxBreakdown
 * Side effects: None outside the literature section.
 */
function renderLitMultiTxUI(total, notations, pmUrl, gene) {
  const row = document.getElementById('rowPubMedMultiTx');
  const el  = document.getElementById('pMultiTxTotal');
  if (!row || !el) return;
  row.style.display = '';

  if (total === 0) {
    el.innerHTML = `<span style="color:var(--dim); font-size:0.8rem;">None found</span>`;
    return;
  }

  const caption = notations.length > 1
    ? `${notations.length} notation${notations.length > 1 ? 's' : ''} searched`
    : (notations[0] || '');

  const breakdownBtn = notations.length > 1
    ? `<button onclick="loadLitBreakdown('${gene}')"
         style="font-size:0.65rem; color:var(--blue); background:none; border:none;
                cursor:pointer; padding:0; margin-top:2px; text-align:right;">
         ▸ by notation
       </button>`
    : '';

  el.innerHTML = `
    <div style="display:flex; flex-direction:column; align-items:flex-end; gap:1px;">
      <a href="${pmUrl}" target="_blank"
         style="font-weight:800; font-size:0.85rem; color:var(--teal); text-decoration:none;">
        ${total} Match${total !== 1 ? 'es' : ''} ↗
      </a>
      <span style="font-size:0.63rem; color:var(--dim);">${caption}</span>
      ${breakdownBtn}
    </div>`;
}

const AA3TO1 = {
  Ala:'A', Arg:'R', Asn:'N', Asp:'D', Cys:'C', Gln:'Q', Glu:'E', Gly:'G',
  His:'H', Ile:'I', Leu:'L', Lys:'K', Met:'M', Phe:'F', Pro:'P', Ser:'S',
  Thr:'T', Trp:'W', Tyr:'Y', Val:'V', Ter:'*', Sec:'U', Pyl:'O'
};

function toShortProtein(pLong) {
  // pLong is like p.Arg123Pro or p.Arg123Ter or p.Arg123*
  // Returns p.R123P or null if can't parse
  const m = pLong.match(/^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)$/);
  if (!m) return null;
  const ref1 = AA3TO1[m[1]];
  const alt1 = m[3] === '*' ? '*' : AA3TO1[m[3]];
  if (!ref1 || !alt1) return null;
  return `p.${ref1}${m[2]}${alt1}`;
}

function getTopDiseases(conditions) {
  if (!conditions || conditions.length === 0) return [];
  const priority = ['definitive', 'strong', 'moderate', 'supportive', 'limited'];
  for (const lvl of priority) {
    const matches = conditions.filter(c => (c.evidenceLevel || '').toLowerCase().includes(lvl));
    if (matches.length > 0) return matches.map(c => c.diseaseName).filter(Boolean);
  }
  return conditions[0]?.diseaseName ? [conditions[0].diseaseName] : [];
}

// Renders the Variant + Disease PubMed link into el given the resolved disease list.
// Query: gene[tiab] AND (pLong[tiab] OR pShort[tiab]) AND (diseases)
// scores (Map<name,count>) is shown as a tooltip when phenotype matching was applied.
function renderVariantDiseaseLink(el, gene, pLong, pShort, diseases, scores) {
  if (!diseases || diseases.length === 0) {
    el.innerHTML = `<span style="color:var(--dim); font-size:0.7rem;">Awaiting disease data…</span>`;
    return;
  }
  const proteinParts = [pLong, pShort].filter(Boolean).map(v => `"${v}"[tiab]`).join(' OR ');
  const diseaseParts = diseases.map(d => `"${d}"[tiab]`).join(' OR ');
  const q = proteinParts
    ? `(${gene}[tiab]) AND (${proteinParts}) AND (${diseaseParts})`
    : `(${gene}[tiab]) AND (${diseaseParts})`;
  const url = `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(q)}`;

  const label = diseases.length === 1
    ? diseases[0]
    : `${diseases.length} conditions`;

  const tooltip = scores
    ? diseases.map(d => `${d} (${scores.get(d) ?? 0} HPO match)`).join('; ')
    : diseases.join('; ');

  el.innerHTML = `<a href="${url}" target="_blank" style="color:var(--teal); text-decoration:none; font-size:0.75rem;">Search PubMed ↗</a><span style="font-size:0.62rem; color:var(--dim); margin-left:6px;" title="${tooltip}">${label}</span>`;
}

function renderLitSearchLinks() {
  // Search 1: Gene + All Phenotypes
  const el1 = document.getElementById('pSearchGenePhenotype');
  let q1 = null;
  if (el1) {
    if (data.coords.gene && data.ptPhenoTexts && data.ptPhenoTexts.length > 0) {
      const phenoTerms = data.ptPhenoTexts.map(t => `"${t}"[tiab]`).join(' OR ');
      q1 = `(${data.coords.gene}[tiab]) AND (${phenoTerms}) NOT "Review"[Publication Type]`;
      const url1 = `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(q1)}`;
      el1.innerHTML = `<a href="${url1}" target="_blank" style="color:var(--teal); text-decoration:none; font-size:0.75rem;">Search PubMed ↗</a>`;
    } else {
      el1.innerHTML = `<span style="color:var(--dim); font-size:0.7rem;">—</span>`;
    }
  }

}

/**
 * Renders Literature Review and Phenotype Match results with V81 Premium styling.
 */
function renderPhenoAndLitCards(phenoData = null, litData = null) {
  // 1. Phenotype & Gene Match
  if (phenoData) {
    const { score, total, matchedArr, directCount = 0, ancestorCount = 0 } = phenoData;
    const pct = total > 0 ? (score / total) * 100 : 0;

    // Match Score with Bullet Graph
    const scoreEl = document.getElementById('pScore');
    if (scoreEl) {
      scoreEl.innerHTML = `
        <div style="display:flex; flex-direction:column; gap:6px; width:100%;">
          <div style="display:flex; justify-content:space-between; align-items:center;">
             <span style="font-weight:800; color:var(--teal);">${score}/${total}</span>
             <span style="font-size:0.65rem; color:var(--dim);">${pct.toFixed(0)}% Patient Match</span>
          </div>
          <div class="bullet-container" style="height:6px;">
            <div class="bullet-bar" style="width:${pct}%; background:var(--teal);"></div>
          </div>
        </div>`;
    }

    // Matched Terms as Chips
    const termsEl = document.getElementById('pTerms');
    if (termsEl) {
      if (matchedArr && matchedArr.length > 0) {
        termsEl.innerHTML = matchedArr.map(t => {
          const isDirect = t.includes('(D)');
          const isLitStrong = t.includes('(L+)');
          const isLit = t.includes('(L)');
          const label = t.replace(/\s*\([DCL+\d]+\)\s*$/, '').trim();
          const depth = t.match(/\(C(\d)\)/)?.[1] || '';
          let cls, title;
          if (isDirect)      { cls = 'depth-direct';   title = 'Direct HPO match'; }
          else if (isLitStrong){ cls = 'depth-lit-strong'; title = 'Literature match (PubMed + ClinVar P/LP)'; }
          else if (isLit)    { cls = 'depth-lit';      title = 'Literature match (PubMed)'; }
          else               { cls = 'depth-clinical'; title = `HPO Ancestor match (depth ${depth})`; }
          const suffix = depth ? ` • ${depth}` : '';
          return `<span class="pheno-chip ${cls}" data-tippy-content="${title}">${label}${suffix}</span>`;
        }).join('');
      } else {
        termsEl.innerText = 'None';
      }
    }

    // Match Summary
    const linEl = document.getElementById('pLin');
    if (linEl) {
      const { directCount = 0, ancestorCount = 0, litCount = 0 } = phenoData;
      const parts = [];
      if (directCount > 0)  parts.push(`<span style="color:var(--teal);font-weight:800;">${directCount} HPO-Direct</span>`);
      if (ancestorCount > 0) parts.push(`<span style="color:var(--amber);font-weight:800;">${ancestorCount} HPO-Ancestor</span>`);
      if (litCount > 0)     parts.push(`<span style="color:#a78bfa;font-weight:800;">${litCount} Literature</span>`);
      if (parts.length === 0) parts.push(`<span style="color:var(--dim);">No match found</span>`);
      linEl.innerHTML = parts.join('<span style="opacity:0.4;margin:0 4px;">·</span>');
    }
  }

  // NOTE: The OMIM synopsis row is intentionally NOT cleared here.
  // renderPhenoAndLitCards is re-invoked by mergeLiteratureIntoMatch() after
  // Step 7 has populated the row; clearing it here would wipe the result.
  // Visibility is owned by renderOmimSynopsisSummary() (shows/hides based on
  // match) and reset on a fresh search by clearDash().

  // 2. Literature Review
  if (litData) {
    const { tiab } = litData;
    const base = 'https://pubmed.ncbi.nlm.nih.gov/?term=';

    if (tiab) {
      const el = document.getElementById('pPubTiabRev');
      if (el) {
        el.parentElement.style.backgroundColor = 'transparent';
        const revHtml = tiab.rev > 0
          ? `<a href="${base + tiab.qRev}" target="_blank" class="lit-count-badge lit-badge-rev">${tiab.rev} <small>REV</small></a>`
          : `<span class="lit-count-badge" style="color:var(--dim); opacity:0.4;">0</span>`;
        const nonHtml = tiab.non > 0
          ? `<a href="${base + tiab.qNon}" target="_blank" class="lit-count-badge lit-badge-non">${tiab.non} <small>ORIG</small></a>`
          : `<span class="lit-count-badge" style="color:var(--dim); opacity:0.4;">0</span>`;
        el.innerHTML = `<div style="display:flex; gap:8px;">${revHtml} ${nonHtml}</div>`;
        const tiabNonEl = document.getElementById('pPubTiabNon');
        if (tiabNonEl) tiabNonEl.style.display = 'none';
      }
    }
  }
}

async function populateSummaryRows(gene) {
  const base = 'https://pubmed.ncbi.nlm.nih.gov/?term=';
  const cvBase = 'https://www.ncbi.nlm.nih.gov/clinvar/?term=';

  // ── All phenotypes row ──────────────────────────────────────────────────
  if (gene && data.ptPhenoTexts?.length > 0) {
    const tiabTerms = data.ptPhenoTexts.map(t => `"${t}"[tiab]`).join(' OR ');
    const noReview = ' NOT "Review"[Publication Type]';
    const allTiabQ  = `(${gene}[tiab]) AND (${tiabTerms})${noReview}`;

    const meshValues = Object.values(data.ptPhenoMeSH || {});
    const allMeshQ   = meshValues.length
      ? `(${gene}[tiab]) AND (${meshValues.map(m => `"${m}"[MeSH Terms]`).join(' OR ')})${noReview}` : null;

    const cvDisTerms = data.ptPhenoTexts.map(t => `"${t}"[dis]`).join(' OR ');
    const allCvQ     = `${gene}[gene] AND (${cvDisTerms}) AND (Pathogenic[clnsig] OR "Likely pathogenic"[clnsig])`;

    const [tiabCount, meshCount, cvCount] = await Promise.all([
      fetchNcbiCount(encodeURIComponent(allTiabQ)),
      allMeshQ ? fetchNcbiCount(encodeURIComponent(allMeshQ)) : Promise.resolve(0),
      fetchNcbiCountClinVar(encodeURIComponent(allCvQ))
    ]);

    renderSummaryCell('sumAllPhenoTiab', tiabCount, base + encodeURIComponent(allTiabQ), 'pubmed');
    renderSummaryCell('sumAllPhenoMesh', meshCount, allMeshQ ? base + encodeURIComponent(allMeshQ) : null, 'pubmed');
    renderSummaryCell('sumAllPhenoCv',   cvCount,   cvBase + encodeURIComponent(allCvQ), 'clinvar');
  } else {
    ['sumAllPhenoTiab', 'sumAllPhenoMesh', 'sumAllPhenoCv'].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.innerHTML = `<span style="opacity:0.2;">—</span>`;
    });
  }

  // ── Variant + Disease row ───────────────────────────────────────────────
  let pRaw = data.ensembl.protein && data.ensembl.protein !== '-'
    ? data.ensembl.protein.replace(/^[A-Z0-9_.]+:/, '').replace(/[()]/g, '')
    : null;
  if (pRaw && !pRaw.startsWith('p.')) pRaw = `p.${pRaw}`;
  const pLong = pRaw ? pRaw.replace(/^p\./, '') : null;
  const pShortFull = pRaw ? toShortProtein(pRaw) : null;
  const pShort = pShortFull ? pShortFull.replace(/^p\./, '') : null;

  const definitives = (data.associatedConditions || [])
    .filter(c => (c.evidenceLevel || '').toLowerCase().includes('definitive'));

  let topDiseases, phenoScores = new Map(), phenoDetails = new Map();
  if (definitives.length > 1 && data.ptPhenoGroups?.length) {
    const { diseases, scores, details } = await selectDiseasesByPhenoMatch(definitives);
    topDiseases = diseases.length ? diseases : getTopDiseases(data.associatedConditions);
    phenoScores = scores || new Map();
    phenoDetails = details || new Map();
    console.log('[Summary] HPO match selected:', diseases, '→ topDiseases:', topDiseases);
  } else if (definitives.length > 0) {
    topDiseases = definitives.map(c => c.diseaseName);
    console.log('[Summary] No phenotype groups, using all definitives:', topDiseases);
  } else {
    topDiseases = getTopDiseases(data.associatedConditions);
    console.log('[Summary] No definitives, using getTopDiseases:', topDiseases);
  }

  if (gene && topDiseases.length) {
    // Clean disease names: remove type specifiers like "type 1", "type 1A", or trailing digits
    const cleanDisease = (name) => name
      .replace(/,\s*type\s+[0-9]+[A-Z]*/gi, '')  // Remove ", type 1" or ", type 1A"
      .replace(/,\s*[0-9]+[A-Z]*$/gi, '')         // Remove ", 10" or ", 1A" at end
      .replace(/\s+[0-9]+[A-Z]*$/i, '')           // Remove " 3" or " 1E" at end
      .trim();

    const cleanedDiseases = topDiseases.map(cleanDisease);
    const proteinParts   = [pLong, pShort].filter(Boolean).map(v => `"${v}"`).join(' OR ');
    const diseaseParts = cleanedDiseases.map(d => `"${d}"`).join(' OR ');
    const diseaseMeshParts = topDiseases.map(d => `"${d}"[MeSH Terms]`).join(' OR ');
    const noReview = ' NOT "Review"[Publication Type]';

    const varTiabQ = (proteinParts
      ? `${gene} AND (${proteinParts}) AND (${diseaseParts})`
      : `${gene} AND (${diseaseParts})`) + noReview;
    const varMeshQ = (proteinParts
      ? `${gene} AND (${proteinParts}) AND (${diseaseMeshParts})`
      : `${gene} AND (${diseaseMeshParts})`) + noReview;

    console.log('[Variant+Disease Query]');
    console.log('  gene:', gene);
    console.log('  protein parts:', proteinParts || '(none)');
    console.log('  disease parts:', diseaseParts);
    console.log('  topDiseases:', topDiseases);
    console.log('  TIAB Query:', varTiabQ);
    console.log('  MESH Query:', varMeshQ);
    console.log('  PubMed URL:', `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(varTiabQ)}`);

    const [tiabCount, meshCount] = await Promise.all([
      fetchNcbiCount(encodeURIComponent(varTiabQ)),
      fetchNcbiCount(encodeURIComponent(varMeshQ))
    ]);
    console.log('[Variant+Disease Counts] TIAB:', tiabCount, 'MESH:', meshCount);

    renderSummaryCell('sumVariantTiab', tiabCount, base + encodeURIComponent(varTiabQ), 'pubmed');
    renderSummaryCell('sumVariantMesh', meshCount, base + encodeURIComponent(varMeshQ), 'pubmed');

    // Display ALL definitive diseases with HPO match counts; highlight selected ones
    const diseaseListEl = document.getElementById('sumVariantDiseaseList');
    if (diseaseListEl) {
      const hasPhenoData = data.ptPhenoGroups?.length > 0;
      const selectedSet = new Set(topDiseases);
      const allToShow = definitives.length > 0
        ? definitives.map(c => c.diseaseName).filter(Boolean)
        : topDiseases;

      if (allToShow.length === 0) {
        diseaseListEl.innerHTML = `<span style="opacity:0.3;">—</span>`;
      } else {
        diseaseListEl.innerHTML = allToShow
          .map(d => {
            const score = phenoScores.get(d) ?? 0;
            const matchedTerms = phenoDetails.get(d) || [];
            const matchTooltip = matchedTerms.length > 0 ? matchedTerms.join(' | ') : 'No phenotype matches';
            const isSelected = selectedSet.has(d);
            const nameColor = isSelected ? 'var(--text)' : 'var(--dim)';
            const scoreBadge = hasPhenoData
              ? ` <span style="color:${score > 0 ? 'var(--teal)' : 'var(--dim)'}; font-size:0.65rem; font-weight:${score > 0 ? '700' : '400'};" title="${matchTooltip}">(${score} HPO)</span>`
              : '';
            const queryMark = isSelected
              ? ` <span style="color:var(--teal); font-size:0.6rem;">✓</span>`
              : '';
            return `<div style="color:${nameColor}; padding:1px 0; font-style:normal;" title="${matchTooltip}">${d}${scoreBadge}${queryMark}</div>`;
          })
          .join('');
      }
    }
  } else {
    ['sumVariantTiab', 'sumVariantMesh'].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.innerHTML = `<span style="opacity:0.2;">—</span>`;
    });
    const diseaseListEl = document.getElementById('sumVariantDiseaseList');
    if (diseaseListEl) diseaseListEl.innerHTML = `<span style="opacity:0.3;">—</span>`;
  }
}

/**
 * Re-runs only the Variant + Disease row of populateSummaryRows.
 * Called when associatedConditions arrive after the full populateSummaryRows
 * already ran — avoids re-firing the 3 All-Phenotypes NCBI calls unnecessarily.
 */
async function populateVariantDiseaseRow(gene) {
  if (!document.getElementById('sumVariantDiseaseList')) return;
  // Delegate to the full function but with a guard: if the All Phenotypes cells
  // already have real values, skip re-fetching them by checking if the Variant+
  // Disease list is still "—" (needs update) vs already populated.
  const diseaseListEl = document.getElementById('sumVariantDiseaseList');
  if (!diseaseListEl || !diseaseListEl.innerText.trim().startsWith('—')) return;
  return populateSummaryRows(gene);
}

function renderSummaryCell(id, count, url, type) {
  const el = document.getElementById(id);
  if (!el) return;
  const color = type === 'clinvar' ? 'var(--amber)' : 'var(--teal)';
  el.innerHTML = (count > 0 && url)
    ? `<a href="${url}" target="_blank" style="color:${color};text-decoration:none;font-weight:700;">${count}</a>`
    : `<span style="opacity:0.3;">0</span>`;
}

/**
 * Renders the OMIM synopsis corroboration summary line below the match score.
 * Shows best disease match, tie groups, or ambiguous message.
 * Called from api.js after matchOmimSynopsis() completes.
 */
function renderOmimSynopsisSummary(results) {
  const rowEl = document.getElementById('pOmimSynopsisRow');
  const valEl = document.getElementById('pOmimSynopsis');
  if (!rowEl || !valEl) return;

  const summary = (typeof getOmimSynopsisSummary === 'function')
    ? getOmimSynopsisSummary(results) : null;

  if (!summary) { rowEl.style.display = 'none'; return; }

  const colorMap = {
    single:    'var(--teal)',
    multiple:  'var(--amber)',
    ambiguous: 'rgba(255,255,255,0.35)'
  };
  const color = colorMap[summary.type] || 'var(--dim)';
  const icon  = summary.type === 'single' ? '✦' : summary.type === 'multiple' ? '≈' : '~';

  valEl.innerHTML =
    `<span style="color:${color};font-size:0.72rem;font-weight:700;">${icon} ${summary.label}</span>` +
    `<span style="color:var(--dim);font-size:0.62rem;display:block;margin-top:2px;opacity:0.8;">${summary.detail}</span>`;

  rowEl.style.display = '';
}

/**
 * Renders the HPO phenotype-fit summary row from data.phenotypeFit.
 * LEADS WITH THE LIKELIHOOD RATIO (log₁₀LR) — prior-independent and the correct
 * thing to rank/trust. Posterior is shown only as secondary fixed-prior context.
 * `coverage` (length-normalized) is surfaced as a cross-check, and a sparse-frequency
 * warning appears when the top disease's annotation completeness is low.
 * Reads only data.phenotypeFit; self-contained (no cross-grid DOM).
 */
function renderPhenotypeFitSummary(fit) {
  const rowEl = document.getElementById('pPhenoFitRow');
  const valEl = document.getElementById('pPhenoFit');
  if (!rowEl || !valEl) return;

  const results = (fit && fit.ready && Array.isArray(fit.results)) ? fit.results : [];
  if (!results.length) { rowEl.style.display = 'none'; return; }

  const top = results[0];
  const lr10 = top.log10LR;
  // Positive log₁₀LR favours the disease; negative argues against. Colour by strength.
  const color = lr10 >= 1 ? 'var(--teal)' : lr10 > 0 ? 'var(--amber)' : 'rgba(255,255,255,0.45)';
  const icon  = lr10 >= 1 ? '✦' : lr10 > 0 ? '≈' : '·';
  const lrStr = (lr10 >= 0 ? '+' : '') + lr10.toFixed(2);
  const cov   = Math.round((top.coverage || 0) * 100);
  const post  = Math.round((top.posterior || 0) * 100);
  const sparse = (top.completeness != null && top.completeness < 0.5);
  const moreN = results.length - 1;

  const detailBits = [
    `LR ${top.totalLR >= 100 ? top.totalLR.toExponential(1) : top.totalLR.toFixed(1)}`,
    `coverage ${cov}%`,
    `post ${post}%`
  ];
  if (moreN > 0) detailBits.push(`+${moreN} more`);
  const sparseNote = sparse
    ? `<span style="color:var(--amber);font-size:0.6rem;display:block;margin-top:1px;">⚠ LR from sparse frequency data — rely on coverage</span>`
    : '';

  valEl.innerHTML =
    `<span style="color:${color};font-size:0.72rem;font-weight:700;">${icon} ${top.name} · log₁₀LR ${lrStr}</span>` +
    `<span style="color:var(--dim);font-size:0.62rem;display:block;margin-top:2px;opacity:0.85;">${detailBits.join(' · ')}</span>` +
    sparseNote;

  rowEl.style.display = '';
}

/**
 * Injects synopsis match badges into already-rendered condition cards.
 * Each card has a placeholder span with id="omim-syn-{mimId}".
 * Discriminating terms (unique to this disease) are starred.
 * Called from api.js after matchOmimSynopsis() completes.
 */
function applyOmimSynopsisBadges(results, ptTotal) {
  if (!results?.length) return;
  for (const r of results) {
    const el = document.getElementById(`omim-syn-${r.mimId}`);
    if (!el) continue;

    if (!r.hasSynopsis || r.matched === 0) {
      el.style.display = 'none';
      continue;
    }

    const pct = Math.round(r.fScore * 100);
    const color = pct >= 66 ? 'var(--teal)' : pct >= 40 ? 'var(--amber)' : 'rgba(255,255,255,0.35)';

    const termLabels = r.matchedTerms.map(t => {
      const prefix = t.discriminating ? '★' : '';
      const suffix = t.matchType === 'text'
        ? '<sup style="color:var(--amber);font-size:0.55rem;">T</sup>'
        : t.matchType === 'ic'
          ? '<sup style="color:var(--amber);font-size:0.55rem;" title="ontology-level match — specificity not confirmed">≈</sup>'
          : '';
      return t.discriminating
        ? `<b style="color:var(--teal);">${prefix}${t.text}${suffix}</b>`
        : `${t.text}${suffix}`;
    }).join(' · ');

    const susceptLabel = r.isSusceptibility
      ? `<span style="color:var(--amber);font-size:0.58rem;"> (susceptibility)</span>` : '';
    const uncertainLabel = r.isUncertain
      ? `<span style="color:var(--dim);font-size:0.58rem;"> (?)</span>` : '';

    const tooltip = `OMIM synopsis · F-score ${pct}%\nMatched: ${r.matchedTerms.map(t => (t.discriminating ? '★' : '') + t.text).join(', ')}\nPrecision: ${Math.round(r.precision * 100)}% · Recall: ${Math.round(r.recall * 100)}%`;

    el.style.display = '';
    el.innerHTML =
      `<span style="cursor:default;" data-tippy-content="${tooltip}">` +
        `<span style="color:${color};font-size:0.68rem;font-weight:700;">` +
          `📋 Synopsis ${r.matched}/${ptTotal ?? '?'}` +
        `</span>` +
        susceptLabel + uncertainLabel +
        `<span style="color:var(--dim);font-size:0.62rem;margin-left:4px;">[${termLabels}]</span>` +
      `</span>`;

    if (window.tippy) {
      tippy(el.querySelector('[data-tippy-content]'), { allowHTML: false, placement: 'top' });
    }
  }
}

function renderPhenoBreakdown(results, gene) {
  const container = document.getElementById('pPhenoBreakdown');
  if (!container) return;

  const base = 'https://pubmed.ncbi.nlm.nih.gov/?term=';
  const cvBase = 'https://www.ncbi.nlm.nih.gov/clinvar/?term=';

  let html = `<div style="margin-top:8px; border-top:1px solid rgba(255,255,255,0.07); padding-top:8px;">
    <div style="font-size:0.65rem; color:var(--dim); margin-bottom:6px; letter-spacing:0.05em; text-transform:uppercase;">Per-Phenotype Evidence</div>
    <table style="width:100%; border-collapse:collapse; font-size:0.7rem;">
      <thead>
        <tr style="color:var(--dim); font-size:0.62rem; text-transform:uppercase; letter-spacing:0.04em;">
          <th style="text-align:left; padding:2px 4px; font-weight:500;">Phenotype</th>
          <th style="text-align:center; padding:2px 4px; font-weight:500;">PubMed</th>
          <th style="text-align:center; padding:2px 4px; font-weight:500;">MeSH</th>
          <th style="text-align:center; padding:2px 4px; font-weight:500;">ClinVar P/LP</th>
        </tr>
      </thead>
      <tbody>`;

  for (const r of (results || [])) {
    const hasEvidence = r.tiab > 0 || r.mesh > 0 || r.clinvar > 0;
    const rowColor = hasEvidence ? 'rgba(255,255,255,0.04)' : 'transparent';
    const termLabel = r.term.length > 22 ? r.term.slice(0, 20) + '…' : r.term;
    const indicator = hasEvidence
      ? `<span style="color:var(--teal);">✓</span>`
      : `<span style="color:var(--dim); opacity:0.4;">—</span>`;

    const tiabCell = r.tiab > 0
      ? `<a href="${base + r.tiabQ}" target="_blank" style="color:var(--teal);text-decoration:none;font-weight:700;">${r.tiab}</a>`
      : `<span style="opacity:0.3;">0</span>`;

    const meshCell = r.meshTerm
      ? (r.mesh > 0
          ? `<a href="${base + r.meshQ}" target="_blank" style="color:var(--teal);text-decoration:none;font-weight:700;" title="MeSH: ${r.meshTerm}">${r.mesh}</a>`
          : `<span style="opacity:0.3;">0</span>`)
      : `<span style="opacity:0.2;">—</span>`;

    const cvCell = r.clinvar > 0
      ? `<a href="${cvBase + r.cvQ}" target="_blank" style="color:var(--amber);text-decoration:none;font-weight:700;">${r.clinvar}</a>`
      : `<span style="opacity:0.3;">0</span>`;

    html += `<tr style="background:${rowColor}; border-bottom:1px solid rgba(255,255,255,0.04);">
      <td style="padding:3px 4px; color:var(--text);" title="${r.term}">${indicator} ${termLabel}</td>
      <td style="text-align:center; padding:3px 4px;">${tiabCell}</td>
      <td style="text-align:center; padding:3px 4px;">${meshCell}</td>
      <td style="text-align:center; padding:3px 4px;">${cvCell}</td>
    </tr>`;
  }

  // Summary rows separator
  html += `
      <tr><td colspan="4" style="padding:4px 0; border-top:1px solid rgba(255,255,255,0.12);"></td></tr>
      <tr style="background:rgba(255,255,255,0.03);">
        <td style="padding:3px 4px; color:var(--teal); font-weight:700; font-size:0.68rem;">All phenotypes</td>
        <td id="sumAllPhenoTiab" style="text-align:center; padding:3px 4px; color:var(--dim); font-size:0.68rem;">…</td>
        <td id="sumAllPhenoMesh" style="text-align:center; padding:3px 4px; color:var(--dim); font-size:0.68rem;">…</td>
        <td id="sumAllPhenoCv"   style="text-align:center; padding:3px 4px; color:var(--dim); font-size:0.68rem;">…</td>
      </tr>
      <tr style="border-top:1px solid rgba(255,255,255,0.1);">
        <td colspan="4" id="sumVariantLabel" style="padding:4px; color:var(--teal); font-weight:700; font-size:0.7rem; background:rgba(255,255,255,0.02);">Variant + Disease</td>
      </tr>
      <tr style="background:rgba(255,255,255,0.01);">
        <td style="padding:3px 4px; color:var(--dim); font-size:0.68rem; font-style:italic;" id="sumVariantDiseaseList">—</td>
        <td id="sumVariantTiab" style="text-align:center; padding:3px 4px; color:var(--dim); font-size:0.68rem;">…</td>
        <td id="sumVariantMesh" style="text-align:center; padding:3px 4px; color:var(--dim); font-size:0.68rem;">…</td>
        <td style="text-align:center; padding:3px 4px; opacity:0.25;">—</td>
      </tr>
    </tbody></table></div>`;

  container.innerHTML = html;
  populateSummaryRows(gene);
}

/**
 * Toggles the gnomAD data source view (Total | Exomes | Genomes)
 */
window.setViewMode = (mode) => {
  data.gnomadViewMode = mode;
  document.querySelectorAll('.toggle-btn').forEach(btn => {
    btn.classList.toggle('active', btn.innerText.toLowerCase() === mode);
  });
  renderPopulationData(mode);
};

/**
 * Renders the population data breakdown based on selected view mode.
 */
function renderPopulationData(mode) {
  if (!data.gnomad.detailed) return;
  const d = data.gnomad.detailed;
  const targetSource = d[mode];
  if (!targetSource) return;

  // 1. Update Grpmax Table AFs
  const fmt = (v) => (v != null) ? (v === 0 ? '0' : (v * 100).toFixed(4) + '%') : '—';
  const exomeVal = d.exome.popmax;
  const genomeVal = d.genome.popmax;

  document.getElementById('exomePopmax').innerText = fmt(exomeVal);
  document.getElementById('genomePopmax').innerText = fmt(genomeVal);

  const totalPopmaxEl = document.getElementById('mvPopmax');
  if (totalPopmaxEl) {
    totalPopmaxEl.innerText = fmt(data.gnomad.popmax);
    totalPopmaxEl.className = 'bold-total ' + (data.gnomad.popmax > 0.05 ? 'bad' : (data.gnomad.popmax > 0 && data.gnomad.popmax < 0.0001) ? 'good' : '');
  }

  // 2. Update Filter Badges
  const getBadge = (filters) => {
    if (filters === null) return '<span class="pop-badge none">No variant</span>';
    if (filters.length === 0) return '<span class="pop-badge pass">Pass</span>';
    return `<span class="pop-badge fail" data-tippy-content="${filters.join(', ')}">${filters[0]}</span>`;
  };
  const exFilt = document.getElementById('exomeFilterBadge');
  const genFilt = document.getElementById('genomeFilterBadge');
  if (exFilt) exFilt.innerHTML = getBadge(d.exome.filters);
  if (genFilt) genFilt.innerHTML = getBadge(d.genome.filters);
  const totFilt = document.getElementById('totalFilterBadge');
  if (totFilt) totFilt.innerHTML = getBadge(d.total.filters);

  // 3. Render Ancestry Groups
  const container = document.getElementById('populationBreakdown');
  if (!container) return;
  container.innerHTML = '';

  // Displayed ancestry groups — hidden groups (asj/ami/fin/mid/remaining) are excluded
  // from display but their counts are included in gnomAD's top-level AC/AN which feeds 'overall'.
  const displayGroups = [
    { key: 'eas', label: 'East Asian' },
    { key: 'sas', label: 'South Asian' },
    { key: 'nfe', label: 'European (non-Finnish)' },
    { key: 'afr', label: 'African/African American' },
    { key: 'amr', label: 'Admixed American' },
  ];

  const isSexChr = ['X', 'Y'].includes(String(data.coords.chrom).replace('chr', '').toUpperCase());

  const formatPop = (p) => {
    let s = `${p.ac}/${p.an} & ${p.hom} homo`;
    if (isSexChr) s += ` & ${p.hemi || 0} hemi`;
    return s;
  };

  // Render the 5 main ancestry groups
  displayGroups.forEach(g => {
    const pop = targetSource.populations[g.key];
    if (!pop || pop.all.an === 0) return;

    const details = document.createElement('details');
    details.className = 'pop-group-details';
    details.innerHTML = `
      <summary class="pop-group-summary">
        <span class="dk">${g.label}</span>
        <span class="pop-sub-val bold-total">${formatPop(pop.all)}</span>
      </summary>
      <div class="pop-sub-row">
        <span class="pop-sub-label">XX (Female)</span>
        <span class="pop-sub-val">${formatPop(pop.XX)}</span>
      </div>
      <div class="pop-sub-row">
        <span class="pop-sub-label">XY (Male)</span>
        <span class="pop-sub-val">${formatPop(pop.XY)}</span>
      </div>
    `;
    container.appendChild(details);
  });

  // Render 'Total (All Pop)' row at the bottom with a visual separator.
  // Uses gnomAD's top-level AC/AN (stored in 'overall') which already aggregates
  // ALL populations including hidden ones (Ashkenazi Jewish, Amish, Finnish,
  // Middle Eastern, Remaining) — so this is always the true complete total.
  const overallPop = targetSource.populations['overall'];
  if (overallPop && overallPop.all.an > 0) {
    // Separator line
    const sep = document.createElement('div');
    sep.style.cssText = 'border-top: 1px solid var(--border); margin: 6px 0 4px 0; opacity: 0.5;';
    container.appendChild(sep);

    const totalDetails = document.createElement('details');
    totalDetails.className = 'pop-group-details total-row';
    totalDetails.innerHTML = `
      <summary class="pop-group-summary" style="background: rgba(255,255,255,0.03);">
        <span class="dk" style="font-weight: 800; color: var(--text-bright);">Total (All Pop)</span>
        <span class="pop-sub-val bold-total" style="font-weight: 800; color: var(--text-bright);">${formatPop(overallPop.all)}</span>
      </summary>
      <div class="pop-sub-row">
        <span class="pop-sub-label">XX (Female)</span>
        <span class="pop-sub-val">${formatPop(overallPop.XX)}</span>
      </div>
      <div class="pop-sub-row">
        <span class="pop-sub-label">XY (Male)</span>
        <span class="pop-sub-val">${formatPop(overallPop.XY)}</span>
      </div>
    `;
    container.appendChild(totalDetails);
  }
}

/**
 * Renders the MaveDB Functional Assay Card
 */
function renderMaveDBCard() {
  const container = document.getElementById('maveDBCardContainer');
  if (!container) return;

  if (!data.maveData) {
    container.style.display = 'none';
    return;
  }

  const { score, interpretation } = data.maveData;
  const interpLower = interpretation.toLowerCase();

  let badgeColor = 'var(--dim)';
  let bgHighlight = 'rgba(255,255,255,0.02)';
  let acmgCode = null;

  if (interpLower.includes('damaging') || interpLower.includes('pathogenic') || interpLower.includes('loss of function')) {
    badgeColor = 'var(--red)';
    bgHighlight = 'rgba(255, 68, 68, 0.05)';
    acmgCode = 'PS3';
  } else if (interpLower.includes('normal') || interpLower.includes('benign') || interpLower.includes('wild type')) {
    badgeColor = 'var(--teal)';
    bgHighlight = 'rgba(0, 200, 150, 0.05)';
    acmgCode = 'BS3';
  } else {
    badgeColor = 'var(--amber)';
    bgHighlight = 'rgba(255, 191, 0, 0.05)';
  }

  container.innerHTML = `
    <div class="data-card" style="border-left: 4px solid ${badgeColor}; background: ${bgHighlight}; margin-top:12px; padding: 12px 16px;">
      <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;">
        <h3 style="color:${badgeColor}; margin:0; border:none; padding:0; font-size:0.78rem;">🔬 MaveDB Functional Assay</h3>
        ${acmgCode ? `<span style="background:${badgeColor}; color:#000; font-weight:800; padding:2px 8px; border-radius:4px; font-size:0.65rem;">Auto-Triggered ${acmgCode}</span>` : ''}
      </div>
      <div style="display:flex; justify-content:space-between; margin-bottom: 4px;">
        <span style="font-size:0.7rem; color:var(--dim);">Functional Score</span>
        <span style="font-size:0.7rem; font-weight:800; color:var(--text-bright);">${score}</span>
      </div>
      <div style="display:flex; justify-content:space-between;">
        <span style="font-size:0.7rem; color:var(--dim);">Interpretation</span>
        <span style="font-size:0.7rem; font-weight:800; color:${badgeColor};">${interpretation}</span>
      </div>
    </div>
  `;
  container.style.display = 'block';
}

function renderDomainAnnotation() {
  const container = document.getElementById('domainAnnotationContainer');
  if (!container) return;
  const features = data.uniprotDomains || [];
  const pos = data.vepProteinStart ? parseInt(data.vepProteinStart, 10) : null;
  if (!features.length) { container.style.display = 'none'; return; }

  const TIER_STYLE = {
    critical: { color: 'var(--red)', label: 'Critical' },
    high: { color: 'var(--amber)', label: 'High' },
    moderate: { color: '#29b6f6', label: 'Moderate' },
    low: { color: 'var(--dim)', label: 'Low' },
  };

  const hits = pos ? features.filter(d => pos >= d.start && pos <= d.end) : [];
  const topTier = hits.find(d => d.tier === 'critical') || hits.find(d => d.tier === 'high')
               || hits.find(d => d.tier === 'moderate') || hits[0];

  const borderColor = topTier ? TIER_STYLE[topTier.tier].color : 'var(--dimmer)';
  const accUrl = data.uniprotAccession
    ? `https://www.uniprot.org/uniprotkb/${data.uniprotAccession}/feature-viewer`
    : '';

  let hitHtml = '';
  if (hits.length && pos) {
    hitHtml = hits.map(h => {
      const s = TIER_STYLE[h.tier];
      return `<span style="background:${s.color};color:#000;font-size:0.62rem;font-weight:900;
                    padding:1px 5px;border-radius:3px;text-transform:uppercase;margin-right:4px;">
                ${s.label}</span>
              <b style="color:var(--text-bright);font-size:0.78rem;">${h.name}</b>
              <span style="color:var(--dim);font-size:0.7rem;"> (${h.type}, aa ${h.start}–${h.end})</span>`;
    }).join('<br>');
  } else if (pos) {
    hitHtml = `<span style="color:var(--dim);font-size:0.75rem;">Position ${pos}: no annotated functional domain</span>`;
  }

  const featureRows = features.map(f => {
    const s = TIER_STYLE[f.tier];
    const isHit = pos && pos >= f.start && pos <= f.end;
    return `<div style="display:flex;align-items:center;gap:8px;padding:3px 0;
                  ${isHit ? 'background:rgba(255,255,255,0.05);border-radius:4px;padding:4px 6px;' : ''}">
              <span style="background:${s.color};color:#000;font-size:0.58rem;font-weight:800;
                           padding:1px 4px;border-radius:2px;min-width:52px;text-align:center;">${s.label}</span>
              <span style="color:var(--text);font-size:0.72rem;flex:1;">${f.name}</span>
              <span style="color:var(--dim);font-size:0.65rem;font-family:'JetBrains Mono',monospace;">
                aa ${f.start}–${f.end}</span>
              ${isHit ? '<span style="color:var(--teal);font-size:0.65rem;font-weight:700;">← HERE</span>' : ''}
            </div>`;
  }).join('');

  container.innerHTML = `
    <div class="data-card" style="border-left:4px solid ${borderColor};background:rgba(255,255,255,0.02);padding:12px 16px;margin-top:0;">
      <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px;">
        <h3 style="color:${borderColor};margin:0;border:none;padding:0;font-size:0.78rem;">
          🧩 Protein Domain Annotation${pos ? ` — p.${pos}` : ''}
        </h3>
        ${accUrl ? `<a href="${accUrl}" target="_blank"
           style="color:var(--dim);font-size:0.65rem;text-decoration:none;">
           UniProt ${data.uniprotAccession} ↗</a>` : ''}
      </div>
      ${hitHtml ? `<div style="margin-bottom:10px;padding:8px;background:rgba(0,0,0,0.2);border-radius:6px;line-height:1.8;">
        ${hitHtml}
      </div>` : ''}
      <div style="display:flex;flex-direction:column;gap:1px;max-height:200px;overflow-y:auto;scrollbar-width:thin;">
        ${featureRows}
      </div>
    </div>`;
  container.style.display = 'block';
}

function renderVariantDistribution(gene, geneData, clinvarDate) {
  const container = document.getElementById('variantDistributionContainer');
  if (!container) return;
  container.style.display = 'block';
  if (!geneData) {
    container.innerHTML = `<div style="display:flex;align-items:center;justify-content:center;height:100%;min-height:60px;color:var(--dim);font-size:0.72rem;letter-spacing:0.06em;text-transform:uppercase;opacity:0.25;">📊 ClinVar Variant Distribution</div>`;
    return;
  }

  const ROWS = [
    { key: 'pathogenic',        label: 'Pathogenic',        color: '#e53935' },
    { key: 'likely_pathogenic', label: 'Likely Pathogenic', color: '#ef6c00' },
    { key: 'conflicting',       label: 'Conflicting',       color: '#8e24aa' },
    { key: 'vus',               label: 'VUS',               color: '#757575' },
    { key: 'likely_benign',     label: 'Likely Benign',     color: '#00897b' },
    { key: 'benign',            label: 'Benign',            color: '#1e88e5' },
  ];
  const COLS = [
    { key: 'lof',         label: 'LOF' },
    { key: 'missense',    label: 'Missense + Inframe' },
    { key: 'noncoding',   label: 'Non-coding' },
    { key: 'synonymous',  label: 'Synonymous' },
    { key: 'total',       label: 'Total' },
  ];

  const maxVal = Math.max(1, ...ROWS.map(r => (geneData[r.key]?.total || 0)));

  const headerCells = COLS.map(c =>
    `<th style="padding:6px 10px;font-size:0.65rem;color:var(--dim);font-weight:600;text-align:center;white-space:nowrap;">${c.label}</th>`
  ).join('');

  const bodyRows = ROWS.map(r => {
    const row = geneData[r.key] || {};
    const cells = COLS.map(c => {
      const n = row[c.key] || 0;
      if (n === 0) return `<td style="text-align:center;padding:5px 8px;font-size:0.72rem;color:var(--dim);">—</td>`;
      const intensity = c.key === 'total' ? 0.35 : Math.min(0.6, 0.1 + (n / maxVal) * 0.5);
      const isTotal   = c.key === 'total';
      return `<td style="text-align:center;padding:5px 8px;font-size:${isTotal ? '0.78' : '0.72'}rem;
        font-weight:${isTotal ? '800' : '500'};
        background:${r.color}${Math.round(intensity * 255).toString(16).padStart(2,'0')};
        border-radius:4px;color:#fff;">${n.toLocaleString()}</td>`;
    }).join('');
    return `<tr>
      <td style="padding:5px 10px;font-size:0.7rem;font-weight:600;white-space:nowrap;">
        <span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:${r.color};margin-right:6px;"></span>
        ${r.label}
      </td>
      ${cells}
    </tr>`;
  }).join('');

  const totRow = (() => {
    const t = geneData['total'] || {};
    return `<tr style="border-top:1px solid var(--border);">
      <td style="padding:5px 10px;font-size:0.7rem;font-weight:800;color:var(--dim);">Total</td>
      ${COLS.map(c => `<td style="text-align:center;padding:5px 8px;font-size:0.72rem;font-weight:700;color:var(--text);">${(t[c.key]||0).toLocaleString()}</td>`).join('')}
    </tr>`;
  })();

  const total = geneData['total']?.total || 0;

  container.innerHTML = `
    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;">
      <div>
        <h3 style="margin:0;padding:0;border:none;font-size:0.8rem;">📊 ClinVar Variant Distribution — ${gene}</h3>
        <div style="font-size:0.62rem;color:var(--dim);margin-top:2px;">${total.toLocaleString()} classified variants · ClinVar ${clinvarDate || '—'}</div>
      </div>
    </div>
    <div style="overflow-x:auto;">
      <table style="border-collapse:separate;border-spacing:3px;width:100%;">
        <thead><tr>
          <th style="padding:6px 10px;font-size:0.65rem;color:var(--dim);text-align:left;font-weight:600;">Classification</th>
          ${headerCells}
        </tr></thead>
        <tbody>${bodyRows}${totRow}</tbody>
      </table>
    </div>`;
}

/**
 * Renders the ClinGen VCEP card with permanent audit links.
 */
function renderVCEPCard() {
  const container = document.getElementById('vcepCardContainer');
  if (!container) return;

  if (!data.coords.hgvs || data.coords.hgvs === '-') {
    container.style.display = 'none';
    return;
  }

  // Always update eCAId in the main annotation section
  const caIdEl = document.getElementById('eCAId');
  if (!data.caId) {
    const geneSearchUrl = `https://search.clinicalgenome.org/kb/genes?search=${encodeURIComponent(data.coords.gene || '')}`;
    if (caIdEl) caIdEl.innerText = '—';
    container.innerHTML = `
      <div class="data-card" style="border-left: 4px solid var(--dimmer); margin-top:0; background: rgba(255,255,255,0.02); padding: 10px 15px;">
        <div style="display:flex; justify-content:space-between; align-items:center;">
          <span style="font-size:0.72rem; color:var(--dim);">🔍 ClinGen Registry: No Canonical Allele ID resolved for <b style="color:var(--text)">${data.coords.gene || 'this variant'}</b></span>
          <a href="${geneSearchUrl}" target="_blank" style="font-size:0.7rem; color:var(--teal); white-space:nowrap; margin-left:12px;">Search ClinGen ↗</a>
        </div>
      </div>
    `;
    container.style.display = 'block';
    return;
  }

  // Populate eCAId in main annotation section
  if (caIdEl) {
    const caidUrl = `https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=${data.caId}`;
    caIdEl.innerHTML = `<a href="${caidUrl}" target="_blank" style="color:var(--amber);text-decoration:none;font-weight:700;">${data.caId} ↗</a>`;
  }

  const caidUrl = `https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=${data.caId}`;
  const vcepUrl = `https://search.clinicalgenome.org/kb/genes?search=${encodeURIComponent(data.coords.gene || '')}`;
  const erepoUrl = `https://erepo.clinicalgenome.org/evrepo/?hgncId=&hgvsExpression=${encodeURIComponent((data.ensembl.transcript || '') + ':' + (data.coords.hgvs || ''))}&clinAssertion=&dateRange=Any&gte=&lte=&tableView=false`;

  container.innerHTML = `
    <div class="data-card" style="border-left: 4px solid var(--amber); background: rgba(255, 191, 0, 0.03); margin-top:0; padding: 12px 16px;">
      <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;">
        <h3 style="color:var(--amber); margin:0; border:none; padding:0; font-size:0.78rem;">✨ ClinGen Allele Registry — Permanent Audit Link</h3>
        <a href="${caidUrl}" target="_blank"
           style="background:var(--amber); color:#000; font-weight:800; padding:2px 10px; border-radius:4px; font-size:0.72rem; text-decoration:none; white-space:nowrap;">
          ${data.caId} ↗
        </a>
      </div>
      <div style="display:flex; gap:12px; flex-wrap:wrap; font-size:0.68rem; color:var(--dim);">
        <a href="${vcepUrl}" target="_blank" style="color:var(--teal);">🏛️ VCEP Panels for ${data.coords.gene || ''} ↗</a>
        <a href="${erepoUrl}" target="_blank" style="color:var(--teal);">📋 ClinGen Evidence Repo ↗</a>
      </div>
    </div>
  `;
  container.style.display = 'block';
}

/**
 * Renders the CIViC Somatic Evidence badge in the header metrics.
 */
function renderCivicBadge() {
  // Logic handled directly in fetchCivicData for the search link
}


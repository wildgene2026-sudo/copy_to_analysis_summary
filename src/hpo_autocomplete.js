// HPO/Disease autocomplete dropdown for #hpoInput.
// Queries JAX HPO (phenotypes) + Monarch (diseases) on the comma-separated
// token at the caret. Selecting a phenotype pre-seeds phenotypeCache so the
// existing parsePatientPhenotypes() pipeline doesn't re-lookup the HP ID.

// Global utility for escaping HTML (used by renderHpoPills)
function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, c => (
    { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]
  ));
}

document.addEventListener('DOMContentLoaded', () => {
  const input = document.getElementById('hpoInput');
  const dropdown = document.getElementById('hpoDropdown');
  if (!input || !dropdown) return;

  // Auto-validate and display matched/unmatched HPO terms on input change
  let validateTimer = null;
  input.addEventListener('input', () => {
    // Render immediately using phenotypeCache (populated by dropdown selection instantly).
    // This ensures a just-selected term shows green without waiting for the debounce.
    if (typeof renderHpoPills === 'function') renderHpoPills();
    // Re-run phenotype-dependent analysis (debounced) so a phenotype added AFTER the variant
    // search still feeds the NCBI literature search / HPO-OMIM match / GeneReviews focus.
    if (typeof schedulePhenotypeRerun === 'function') schedulePhenotypeRerun();

    clearTimeout(validateTimer);
    validateTimer = setTimeout(async () => {
      const val = input.value.trim();
      if (val) {
        // Only trigger parsePatientPhenotypes if there's at least one comma or it looks like a list
        if (val.includes(',') || val.length > 30) {
          try {
            await parsePatientPhenotypes(val);
            renderHpoPills();
          } catch (e) {
            if (e.name !== 'AbortError') console.warn('[HPO Parse] Error:', e.message);
          }
        } else {
          renderHpoPills();
        }
      } else {
        renderHpoPills();
      }
    }, 400);
  });

  const DEBOUNCE_MS = 250;
  const MIN_CHARS = 2;
  const MAX_RESULTS_EACH = 10;

  let debounceTimer = null;
  let currentAbort = null;
  let activeIndex = -1;
  let currentItems = []; // {id, name, kind}
  let lastQuery = '';

  function getCurrentToken() {
    const val = input.value;
    const caret = input.selectionStart ?? val.length;
    let start = val.lastIndexOf(',', Math.max(0, caret - 1));
    start = start === -1 ? 0 : start + 1;
    let end = val.indexOf(',', caret);
    if (end === -1) end = val.length;
    const raw = val.slice(start, end);
    const leadingWs = raw.length - raw.trimStart().length;
    const trimmed = raw.trim();
    return {
      start: start + leadingWs,
      end: start + leadingWs + trimmed.length,
      text: trimmed,
    };
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, c => (
      { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]
    ));
  }

  function highlight(name, query) {
    if (!query) return escapeHtml(name);
    const lower = name.toLowerCase();
    const idx = lower.indexOf(query.toLowerCase());
    if (idx === -1) return escapeHtml(name);
    return escapeHtml(name.slice(0, idx))
      + '<mark>' + escapeHtml(name.slice(idx, idx + query.length)) + '</mark>'
      + escapeHtml(name.slice(idx + query.length));
  }

  function closeDropdown() {
    dropdown.classList.remove('open');
    dropdown.innerHTML = '';
    activeIndex = -1;
    currentItems = [];
  }

  function render(query, phenotypes, diseases) {
    currentItems = [
      ...phenotypes.map(p => ({ ...p, kind: 'phenotype' })),
      ...diseases.map(d => ({ ...d, kind: 'disease' })),
    ];

    if (currentItems.length === 0) {
      dropdown.innerHTML = '<div class="dropdown-empty">No matches</div>';
      dropdown.classList.add('open');
      return;
    }

    let html = '';
    if (phenotypes.length) {
      html += `<div class="dropdown-section-header">Phenotypes &mdash; ${phenotypes.length} shown</div>`;
      phenotypes.forEach((p, i) => {
        html += `<div class="dropdown-item" data-idx="${i}">`
          + `<span class="term-name">${highlight(p.name, query)}</span>`
          + `<span class="term-id">${escapeHtml(p.id)}</span>`
          + `</div>`;
      });
    }
    if (diseases.length) {
      const offset = phenotypes.length;
      html += `<div class="dropdown-section-header">Diseases &mdash; ${diseases.length} shown</div>`;
      diseases.forEach((d, i) => {
        html += `<div class="dropdown-item" data-idx="${offset + i}">`
          + `<span class="term-name">${highlight(d.name, query)}</span>`
          + `<span class="term-id">${escapeHtml(d.id)}</span>`
          + `</div>`;
      });
    }
    dropdown.innerHTML = html;
    dropdown.classList.add('open');
    activeIndex = -1;

    dropdown.querySelectorAll('.dropdown-item').forEach(el => {
      el.addEventListener('mousedown', e => {
        e.preventDefault(); // keep input focused so caret stays valid
        selectItem(parseInt(el.dataset.idx, 10));
      });
    });
  }

  function setActive(idx) {
    const items = dropdown.querySelectorAll('.dropdown-item');
    items.forEach(el => el.classList.remove('active'));
    if (idx < 0 || idx >= items.length) { activeIndex = -1; return; }
    activeIndex = idx;
    items[idx].classList.add('active');
    items[idx].scrollIntoView({ block: 'nearest' });
  }

  function selectItem(idx) {
    const item = currentItems[idx];
    if (!item) return;

    // Pre-seed phenotypeCache so parsePatientPhenotypes() skips the HTTP lookup.
    if (item.kind === 'phenotype' && typeof phenotypeCache !== 'undefined') {
      phenotypeCache[item.name] = [item.id];
    }

    const { start, end } = getCurrentToken();
    const val = input.value;
    const before = val.slice(0, start);
    const afterRaw = val.slice(end);
    const afterTrim = afterRaw.trimStart();
    // Always end the inserted term with ", " so the user can keep typing.
    const suffix = afterTrim.length === 0 ? ', ' : ', ';
    const newVal = before + item.name + suffix + afterTrim;
    input.value = newVal;
    const cursor = (before + item.name + suffix).length;
    input.setSelectionRange(cursor, cursor);
    closeDropdown();
    if (typeof renderHpoPills === 'function') renderHpoPills();
    input.focus();
    // Parse the full input value so data.ptPhenoTexts is updated immediately
    // (setting input.value programmatically doesn't fire the 'input' event).
    const fullVal = input.value.trim();
    if (fullVal && typeof parsePatientPhenotypes === 'function') {
      parsePatientPhenotypes(fullVal).then(() => {
        if (typeof renderHpoPills === 'function') renderHpoPills();
      }).catch(() => {});
    }
    // A term was picked from the dropdown — re-run phenotype-dependent analysis (debounced).
    if (typeof schedulePhenotypeRerun === 'function') schedulePhenotypeRerun();
  }

  async function search(query) {
    if (currentAbort) currentAbort.abort();
    const ctrl = new AbortController();
    currentAbort = ctrl;
    lastQuery = query;

    dropdown.innerHTML = '<div class="dropdown-loading">Searching…</div>';
    dropdown.classList.add('open');

    const safeFetch = (url) => fetch(url, { signal: ctrl.signal })
      .then(r => r.ok ? r.json() : null)
      .catch(e => { if (e.name === 'AbortError') throw e; return null; });

    try {
      const [phRes, dxRes] = await Promise.all([
        safeFetch(`https://ontology.jax.org/api/hp/search?q=${encodeURIComponent(query)}&page=0&limit=${MAX_RESULTS_EACH}`),
        safeFetch(`https://api.monarchinitiative.org/v3/api/search?q=${encodeURIComponent(query)}&category=biolink:Disease&limit=${MAX_RESULTS_EACH}`),
      ]);
      if (ctrl.signal.aborted || lastQuery !== query) return;

      const phenotypes = ((phRes && phRes.terms) || [])
        .slice(0, MAX_RESULTS_EACH)
        .map(t => ({ id: t.id, name: t.name }));

      const diseases = ((dxRes && dxRes.items) || [])
        .filter(it => /^(OMIM|ORPHA|MONDO|DOID):/i.test(it.id))
        .slice(0, MAX_RESULTS_EACH)
        .map(it => ({ id: it.id, name: it.name }));

      render(query, phenotypes, diseases);
    } catch (e) {
      if (e.name === 'AbortError') return;
      console.warn('HPO autocomplete error:', e);
      dropdown.innerHTML = '<div class="dropdown-empty">Search error</div>';
      dropdown.classList.add('open');
    }
  }

  function onInput() {
    clearTimeout(debounceTimer);
    const { text } = getCurrentToken();
    if (text.length < MIN_CHARS) {
      if (currentAbort) currentAbort.abort();
      closeDropdown();
      return;
    }
    debounceTimer = setTimeout(() => search(text), DEBOUNCE_MS);
  }

  input.addEventListener('input', onInput);
  input.addEventListener('click', onInput);
  input.addEventListener('keyup', e => {
    if (['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(e.key)) onInput();
  });

  input.addEventListener('keydown', e => {
    if (!dropdown.classList.contains('open')) return;
    const items = dropdown.querySelectorAll('.dropdown-item');
    if (!items.length) return;
    if (e.key === 'ArrowDown') {
      e.preventDefault();
      setActive(activeIndex + 1 >= items.length ? 0 : activeIndex + 1);
    } else if (e.key === 'ArrowUp') {
      e.preventDefault();
      setActive(activeIndex - 1 < 0 ? items.length - 1 : activeIndex - 1);
    } else if (e.key === 'Enter' && activeIndex >= 0) {
      e.preventDefault();
      selectItem(activeIndex);
    } else if (e.key === 'Tab' && activeIndex >= 0) {
      e.preventDefault();
      selectItem(activeIndex);
    } else if (e.key === 'Escape') {
      closeDropdown();
    }
  });

  document.addEventListener('click', e => {
    if (!input.contains(e.target) && !dropdown.contains(e.target)) {
      closeDropdown();
    }
  });
});

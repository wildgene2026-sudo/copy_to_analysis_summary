// Two-source variant-annotation validation engine.
//
// Consumed by runValidationPass() in ui.js, which calls window.Validation.validateAll(data)
// after each API completes. Reads raw per-source values from data.coords / data.ensembl /
// data.clingenAR / data.myvariant; writes per-field state into data.validation[field]:
//   { value, sources, available, validated, locked, conflict, singleSource }
//
// Rules (locked with user 2026-05-28):
//   1. Exact string equality — no normalization.
//   2. VV must be one of the two agreeing sources. If VV is unreachable (didn't report),
//      any two non-VV sources may validate instead.
//   3. Lock on validate: once ≥2 sources agree, freeze the ✅ — later sources cannot
//      flip it. Disagreeing late values are recorded in `conflict` for the tooltip.

window.Validation = (function () {
  // Per-field raw values from each API namespace.
  // Order in each object is informational only; voting logic handles VV specially.
  function collect(data) {
    return {
      gene: {
        VV:      data.coords.gene || null,
        VEP:     data.ensembl._geneFromVep || null,
        ClinGen: data.clingenAR.gene || null
      },
      transcript: {
        VV:      data.ensembl.transcript || null,
        VEP:     data.ensembl._transcriptFromVep || null,
        ClinGen: data.clingenAR.transcript || null
      },
      hgvsC: {
        VV:      (data.coords.hgvs && data.coords.hgvs !== '-') ? data.coords.hgvs : null,
        VEP:     data.ensembl._hgvscFromVep || null,
        ClinGen: data.clingenAR.hgvsC || null
      },
      hgvsP: {
        VV:      data.ensembl.protein || null,
        VEP:     data.ensembl._hgvspFromVep || null,
        ClinGen: data.clingenAR.hgvsP || null
      },
      hg38: {
        VV:        data.coords.hg38String || null,
        VEP:       data.ensembl._hg38FromVep || null,
        ClinGen:   data.clingenAR.hg38String || null,
        MyVariant: data.myvariant.hg38String || null
      },
      hg19: {
        VV:        data.coords.hg19String || null,
        ClinGen:   data.clingenAR.hg19String || null,
        MyVariant: data.myvariant.hg19String || null
      },
      rsId: {
        VEP:       data.ensembl.rsId || null,
        MyVariant: data.myvariant.rsId || null
      },
      // Single-source fields — no cross-validation possible; flagged with ℹ️.
      caId:           { ClinGen: data.clingenAR.caId || null },
      vepConsequence: { VEP: data.ensembl.vepConsequence || null },
      exonIntron:     { VEP: (data.ensembl.vepExon && data.ensembl.vepTotalExons)
                              ? `${data.ensembl.vepExon}/${data.ensembl.vepTotalExons}` : null }
    };
  }

  function validateField(field, sourceMap, slot) {
    // Always refresh available/sources so the UI tooltip is current.
    const present = Object.entries(sourceMap).filter(([, v]) => v != null && v !== '' && v !== '-');
    slot.sources   = Object.fromEntries(present);
    slot.available = present.map(([s]) => s);

    if (slot.singleSource) {
      slot.value = present[0]?.[1] || null;
      slot.validated = false;
      slot.conflict = null;
      return;
    }

    // Group sources by EXACT value (no normalization — rule #1).
    const groups = new Map(); // value -> [sources]
    for (const [src, val] of present) {
      if (!groups.has(val)) groups.set(val, []);
      groups.get(val).push(src);
    }

    // If already locked, just update conflict tracking against the locked value.
    if (slot.locked) {
      const lockedSrcs = groups.get(slot.value) || [];
      const conflicts = [];
      for (const [val, srcs] of groups.entries()) {
        if (val !== slot.value) conflicts.push({ value: val, srcs });
      }
      slot.conflict = conflicts.length ? conflicts : null;
      // Re-affirm winning srcs in `sources` listing
      if (lockedSrcs.length) slot.available = lockedSrcs.concat(conflicts.flatMap(c => c.srcs));
      return;
    }

    if (present.length === 0) {
      slot.value = null; slot.validated = false; slot.conflict = null;
      return;
    }

    // Pick the agreeing group per rule #2:
    //   - If VV reported, the winning group must contain VV AND have ≥2 sources.
    //   - If VV did not report, any group with ≥2 sources wins.
    const vvReported = present.some(([s]) => s === 'VV');
    let winnerVal = null, winnerSrcs = null;

    if (vvReported) {
      for (const [val, srcs] of groups.entries()) {
        if (srcs.includes('VV') && srcs.length >= 2) { winnerVal = val; winnerSrcs = srcs; break; }
      }
    } else {
      for (const [val, srcs] of groups.entries()) {
        if (srcs.length >= 2) { winnerVal = val; winnerSrcs = srcs; break; }
      }
    }

    if (winnerVal != null) {
      slot.value     = winnerVal;
      slot.validated = true;
      slot.locked    = true; // freeze — rule #3
      slot.conflict  = [];
      for (const [val, srcs] of groups.entries()) {
        if (val !== winnerVal) slot.conflict.push({ value: val, srcs });
      }
      if (slot.conflict.length === 0) slot.conflict = null;
      slot.available = winnerSrcs.concat((slot.conflict || []).flatMap(c => c.srcs));
    } else {
      // No agreement yet — surface the VV value if present, else the first available.
      const vvEntry = present.find(([s]) => s === 'VV');
      slot.value     = vvEntry ? vvEntry[1] : present[0][1];
      slot.validated = false;
      slot.conflict  = groups.size > 1
        ? Array.from(groups.entries()).map(([val, srcs]) => ({ value: val, srcs }))
        : null;
    }
  }

  function validateAll(data) {
    if (!data || !data.validation) return;
    const collected = collect(data);
    for (const field of Object.keys(collected)) {
      const slot = data.validation[field];
      if (!slot) continue;
      validateField(field, collected[field], slot);
    }
  }

  return { validateAll };
})();

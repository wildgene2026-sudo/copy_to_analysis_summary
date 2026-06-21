# Variant Curation Engine — Web Edition

A streamlined, browser-only build of the Lab Variant Search tool for ACMG/AMP
variant curation. Enter a variant (and optionally a gene + phenotypes), and the
tool fetches annotation from public genomics APIs in the background to assemble a
copy-ready **Analysis Summary**.

**Live site:** https://wildgene2026-sudo.github.io/lab-variant-search/

## What it does

The interface is intentionally minimal — three inputs and three actions:

- **Variant** (HGVS, e.g. `NM_000059.4:c.5266dupC`, or hg38 coords `12 51662848 C T`)
- **Gene symbol** (for ClinGen VCEP guideline lookup)
- **HPO phenotypes** (autocomplete, for phenotype-fit scoring)

Buttons:
- **📋 Copy To Analysis Summary** — builds a formatted, paste-ready report
- **🤖 Copy AI Prompt for GG Report** — few-shot prompt for AI report drafting
- **🔄 Update / double check reports** — cross-reference an old report with fresh data

## Data sources (all live, public APIs)

VariantValidator · Ensembl VEP · gnomAD · MyVariant.info · SpliceAI ·
ClinGen Allele Registry & VCEP · ClinVar (NCBI E-utils) · JAX HPO ·
Monarch Initiative. ClinGen gene dosage is bundled locally.

## Notes / limitations

- **No backend.** This is a pure static site — no API keys are stored or shipped.
  NCBI E-utils calls run keyless (≈3 req/s), so ClinVar is best-effort and PubMed
  literature search is not used.
- **OMIM data is not included** — it is licensed by Johns Hopkins and cannot be
  redistributed publicly. The OMIM section degrades gracefully.
- **GeneReviews** chapters require a server-side fetch (CORS) and are not shown.
- Report export is a normal browser download (no local folder picker).

## Running locally

It's static — serve the folder with any web server:

```bash
python3 -m http.server 8000
# then open http://localhost:8000
```

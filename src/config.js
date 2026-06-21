// ── CONFIG ──────────────────────────────────────────────────────
// NCBI credentials. On the distributed WINDOWS build, build.sh blanks the two
// literals below (the personal key must not ship to colleagues); the startup
// gate in ncbi_key.js then requires the user to enter their own key, remembered
// in localStorage. The macOS / dev build keeps the baked-in key, so no prompt.
// GitHub Pages build: credentials blanked. ncbi_key.js prompts on first visit; key stored in localStorage.
let NCBI_KEY = '';
let NCBI_EMAIL = '';
if (!NCBI_KEY)   NCBI_KEY   = localStorage.getItem('NCBI_KEY')   || '';
if (!NCBI_EMAIL) NCBI_EMAIL = localStorage.getItem('NCBI_EMAIL') || '';
const NCBI_BASE = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils`;


const POP_NAMES = {
  afr: 'African/African American',
  amr: 'Admixed American',
  asj: 'Ashkenazi Jewish',
  eas: 'East Asian',
  fin: 'Finnish',
  nfe: 'Non-Finnish European',
  sas: 'South Asian',
  mid: 'Middle Eastern',
  ami: 'Amish',
  oth: 'Other',
  remaining: 'Other'
};

const NULL_CONS = [
  'frameshift_variant',
  'stop_gained',
  'splice_acceptor_variant',
  'splice_donor_variant',
  'start_lost',
  'transcript_ablation',
  'stop_lost'
];


const ACMG_CODES = [
  {id:'PVS1',cat:'pvs',desc:'Null variant in gene where LOF is disease mechanism'},
  {id:'PVS1_M',cat:'pm',desc:'ClinGen SVI: Start-loss with next Met, and pathogenic variant in truncated region'},
  {id:'PVS1_S',cat:'pp',desc:'ClinGen SVI: Start-loss with next Met, and no pathogenic variants in truncated region'},
  {id:'PS1',cat:'ps',desc:'Same amino acid change as established pathogenic variant'},
  {id:'PS2',cat:'ps',desc:'De novo (confirmed) in patient with disease, no family history'},
  {id:'PS3',cat:'ps',desc:'Well-established functional studies show damaging effect'},
  {id:'PS4',cat:'ps',desc:'Prevalence significantly higher in affected vs controls'},
  {id:'PM1',cat:'pm',desc:'Located in mutational hotspot or well-established critical domain'},
  {id:'PM2',cat:'pp',desc:'Absent or extremely rare in population databases — applied at PM2_Supporting level'},
  {id:'PM3',cat:'pm',desc:'Detected in trans with pathogenic variant (recessive disease)'},
  {id:'PM4',cat:'pm',desc:'Protein length change due to in-frame indels or stop-loss'},
  {id:'PM5',cat:'pm',desc:'Novel missense at same residue as known pathogenic missense'},
  {id:'PM6',cat:'pm',desc:'Assumed de novo, but not confirmed'},
  {id:'PP1',cat:'pp',desc:'Co-segregation with disease in multiple affected family members'},
  {id:'PP2',cat:'pp',desc:'Missense in gene with low benign missense rate and dominant mechanism'},
  {id:'PP3',cat:'pp',desc:'Computational evidence (REVEL 0.644-0.773 / SpliceAI > 0.20) supports pathogenicity'},
  {id:'PP3_M',cat:'pm',desc:'Computational evidence (REVEL 0.773-0.932) moderately supports pathogenicity'},
  {id:'PP3_S',cat:'ps',desc:'Computational evidence (REVEL >= 0.932) strongly supports pathogenicity'},
  {id:'PP4',cat:'pp',desc:'Patient phenotype/family history highly specific for gene'},
  {id:'PP5',cat:'pp',desc:'Reputable source recently classified as pathogenic'},
  {id:'BA1',cat:'ba',desc:'Allele frequency >5% in gnomAD popmax — stand-alone benign'},
  {id:'BS1',cat:'bs',desc:'Allele frequency greater than expected for disorder'},
  {id:'BS2',cat:'bs',desc:'Observed in healthy adult with full penetrance expected'},
  {id:'BS3',cat:'bs',desc:'Well-established functional studies show no damaging effect'},
  {id:'BS4',cat:'bs',desc:'Lack of segregation in affected members of a family'},
  {id:'BP1',cat:'bp',desc:'Missense variant in gene for which only truncating variants cause disease'},
  {id:'BP2',cat:'bp',desc:'Observed in cis with a pathogenic variant in any inheritance pattern'},
  {id:'BP3',cat:'bp',desc:'In-frame deletions/insertions in repetitive region without known function'},
  {id:'BP4',cat:'bp',desc:'Computational evidence (REVEL 0.183-0.290) suggests benign impact'},
  {id:'BP4_M',cat:'bp',desc:'Computational evidence (REVEL 0.016-0.183) moderately suggests benign impact'},
  {id:'BP4_S',cat:'bs',desc:'Computational evidence (REVEL 0.003-0.016) strongly suggests benign impact'},
  {id:'BP4_VS',cat:'ba',desc:'Computational evidence (REVEL <= 0.003) very strongly suggests benign impact (Standalone)'},
  {id:'BP5',cat:'bp',desc:'Variant found in case with alternate molecular basis for disease'},
  {id:'BP6',cat:'bp',desc:'Reputable source recently classified as benign'},
  {id:'BP7',cat:'bp',desc:'Synonymous variant with no predicted splice impact'},
];

const AA_MAP_1TO3 = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','E':'Glu','Q':'Gln','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val','*':'Ter'};
const AA_MAP_3TO1 = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Glu':'E','Gln':'Q','Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'*','Asx':'B','Glx':'Z','Xaa':'X','STP':'*'};


const DB_DEFS = [
  { id: 'clingen', label: 'ClinGen' },
  { id: 'clinvar', label: 'ClinVar' },
  { id: 'decipher', label: 'Decipher' },
  { id: 'franklin', label: 'Franklin' },
  { id: 'gr', label: 'GeneReviews' },
  { id: 'gnomadv2', label: 'gnomAD v2' },
  { id: 'gnomadv4', label: 'gnomAD v4' },
  { id: 'gtex', label: 'GTEx' },
  { id: 'hgmd', label: 'HGMD' },
  { id: 'liftover', label: 'Liftover' },
  { id: 'mastermind', label: 'Mastermind' },
  { id: 'omim', label: 'OMIM' },
  { id: 'scholar', label: 'Scholar' },
  { id: 'spliceai', label: 'SpliceAI' }
];


const CLINVAR_STAR_MAP = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classifications from unflagged records": 0
};

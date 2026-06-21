// ── API STATUS & VERSION TRACKING ────────────────────────────────

const API_VERSIONS = {
  // Group 1: Variant Annotation
  vv: {
    name: 'VariantValidator',
    abbr: 'VV',
    version: '3.0',
    category: 'Variant Annotation',
    githubUrl: 'https://github.com/openvar/variantValidator'
  },
  vep: {
    name: 'Ensembl VEP',
    abbr: 'VEP',
    version: 'r115',
    category: 'Variant Annotation',
    githubUrl: 'https://github.com/Ensembl/VEP_plugins'
  },
  mv: {
    name: 'MyVariant.info',
    abbr: 'MV',
    version: '4.0',
    category: 'Variant Annotation',
    githubUrl: 'https://github.com/biothings/myvariant.info'
  },
  cv: {
    name: 'ClinVar',
    abbr: 'CV',
    version: 'Mar 2026',
    category: 'Variant Annotation',
    githubUrl: 'https://github.com/ncbi/clinvar'
  },
  sai: {
    name: 'SpliceAI',
    abbr: 'SAI',
    version: '1.3',
    category: 'Variant Annotation',
    githubUrl: 'https://github.com/Illumina/SpliceAI'
  },
  gnomad: {
    name: 'gnomAD',
    abbr: 'gnomAD',
    version: 'v4.1',
    category: 'Variant Annotation',
    githubUrl: 'https://github.com/broadinstitute/gnomad_methods'
  },

  // Group 2: Functional Predictions
  litvar: {
    name: 'LitVar 2.0',
    abbr: 'LitVar',
    version: '2.0',
    category: 'Functional Data',
    githubUrl: 'https://github.com/ncbi/LitVar2'
  },
  mavedb: {
    name: 'MaveDB',
    abbr: 'MaveDB',
    version: 'Live',
    category: 'Functional Data',
    githubUrl: 'https://github.com/ORFeome/MaveDB'
  },
  uniprot: {
    name: 'UniProt',
    abbr: 'UniProt',
    version: '2026-05',
    category: 'Functional Data',
    githubUrl: 'https://github.com/ebi-uniprot/uniprot-rest-api'
  },

  // Group 3: Gene & Phenotype Context
  clingen: {
    name: 'ClinGen',
    abbr: 'ClinGen',
    version: '2026-04-06',
    category: 'Gene Context',
    githubUrl: 'https://github.com/ClinGen/ClinGen_Allele_Registry'
  },
  hpo: {
    name: 'JAX HPO',
    abbr: 'HPO',
    version: 'Live',
    category: 'Gene Context',
    githubUrl: 'https://github.com/TheJacksonLaboratory/hpo_website'
  },
  monarch: {
    name: 'Monarch Initiative',
    abbr: 'Monarch',
    version: 'Live',
    category: 'Gene Context',
    githubUrl: 'https://github.com/monarch-initiative/monarch-ui'
  },
  civic: {
    name: 'CIViC',
    abbr: 'CIViC',
    version: 'Live',
    category: 'Gene Context',
    githubUrl: 'https://github.com/griffithlab/civic-v2'
  },
  panelapp: {
    name: 'PanelApp',
    abbr: 'PanelApp',
    version: 'Live',
    category: 'Gene Context',
    githubUrl: 'https://github.com/genomicsengland/PanelApp'
  },
  ucsc: {
    name: 'UCSC ClinGen Dosage',
    abbr: 'UCSC',
    version: 'hg38',
    category: 'Gene Context',
    githubUrl: 'https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=clinGenHaplo'
  }
};

// Track classification distribution for each API
const apiClassifications = {};

// Initialize all APIs
Object.keys(API_VERSIONS).forEach(apiId => {
  apiClassifications[apiId] = {
    pathogenic: 0,
    likely_pathogenic: 0,
    vus: 0,
    likely_benign: 0,
    benign: 0,
    conflicting: 0,
    status: '',
    color: null
  };
});

/**
 * Calculate dominant color based on variant classification distribution.
 * Returns a color representing the pathogenicity spectrum.
 */
function calculateClassificationColor(classification) {
  const p = classification.pathogenic || 0;
  const lp = classification.likely_pathogenic || 0;
  const vus = classification.vus || 0;
  const lb = classification.likely_benign || 0;
  const b = classification.benign || 0;
  const total = p + lp + vus + lb + b;

  if (total === 0) return '#666'; // Gray for no data

  // Weighted scoring: strong weight to pathogenic end
  const pathogenicityScore = (p * 2 + lp * 1.5 - lb * 1.5 - b * 2) / (total * 2);

  // pathogenicityScore ranges from -1 (all benign) to 1 (all pathogenic)
  if (pathogenicityScore > 0.4) return '#e53935'; // Red (Pathogenic)
  if (pathogenicityScore > 0.1) return '#ef6c00'; // Orange (Likely Pathogenic)
  if (pathogenicityScore > -0.1) return '#757575'; // Gray (VUS)
  if (pathogenicityScore > -0.4) return '#00897b'; // Teal (Likely Benign)
  return '#1e88e5'; // Blue (Benign)
}

/**
 * Updates API classification data based on ClinVar signature.
 * Signature format: "Pathogenic", "Likely pathogenic", "VUS", "Benign", etc.
 */
function updateAPIClassification(apiId, cvSig) {
  if (!apiClassifications[apiId]) return;

  if (!cvSig || cvSig === '-' || cvSig === '—') {
    apiClassifications[apiId].status = 'unknown';
    return;
  }

  const sig = cvSig.toLowerCase();

  // Reset counts
  Object.keys(apiClassifications[apiId]).forEach(k => {
    if (k !== 'status' && k !== 'color') {
      apiClassifications[apiId][k] = 0;
    }
  });

  // Classify based on ClinVar significance
  if (sig.includes('pathogenic') && !sig.includes('likely')) {
    apiClassifications[apiId].pathogenic = 1;
  } else if (sig.includes('likely pathogenic')) {
    apiClassifications[apiId].likely_pathogenic = 1;
  } else if (sig.includes('conflicting')) {
    apiClassifications[apiId].conflicting = 1;
  } else if (sig.includes('benign') && !sig.includes('likely')) {
    apiClassifications[apiId].benign = 1;
  } else if (sig.includes('likely benign')) {
    apiClassifications[apiId].likely_benign = 1;
  } else if (sig.includes('uncertain') || sig.includes('vus')) {
    apiClassifications[apiId].vus = 1;
  } else {
    apiClassifications[apiId].status = 'unknown';
    return;
  }

  // Calculate color and update
  apiClassifications[apiId].color = calculateClassificationColor(apiClassifications[apiId]);
  apiClassifications[apiId].status = 'classified';
  updateAPIBadge(apiId);
}

/**
 * Updates the API badge with current status, version, and classification color.
 */
function updateAPIBadge(apiId) {
  const config = API_VERSIONS[apiId];
  const classification = apiClassifications[apiId];
  const badge = document.getElementById(`api-badge-${apiId}`);

  if (!badge || !config) return;

  const color = classification.color || '#666';
  const status = classification.status || 'pending';

  // Determine display status indicator
  let statusIcon = '⏳';
  let statusColor = '#666';

  if (status === 'error') {
    statusIcon = '✕';
    statusColor = '#ff4d6d';
  } else if (status === 'warn') {
    statusIcon = '⚠';
    statusColor = '#ffb347';
  } else if (status === 'ok') {
    statusIcon = '✓';
    statusColor = '#00d4aa';
  } else if (status === 'loading') {
    statusIcon = '⟳';
    statusColor = '#4d9fff';
  } else if (status === 'classified') {
    statusIcon = '✓';
    statusColor = color;
  }

  badge.innerHTML = `
    <a href="${config.githubUrl}" target="_blank" title="${config.name} v${config.version}" class="api-badge-link">
      <span class="api-badge-icon" style="color: ${statusColor};">${statusIcon}</span>
      <span class="api-badge-name">${config.abbr}</span>
      <span class="api-badge-version">${config.version}</span>
      ${classification.color ? `<span class="api-badge-color-indicator" style="background-color: ${classification.color};"></span>` : ''}
    </a>
  `;
}

/**
 * Set API status (loading, ok, warn, error).
 */
function setAPIStatus(apiId, status) {
  if (!apiClassifications[apiId]) return;
  apiClassifications[apiId].status = status;
  updateAPIBadge(apiId);
}

/**
 * For gene-level analysis: update API classification based on variant distribution.
 */
function updateAPIClassificationFromDistribution(apiId, distributionData) {
  if (!apiClassifications[apiId] || !distributionData) return;

  // Reset counts
  Object.keys(apiClassifications[apiId]).forEach(k => {
    if (k !== 'status' && k !== 'color') {
      apiClassifications[apiId][k] = 0;
    }
  });

  // Sum up classification counts from distribution
  apiClassifications[apiId].pathogenic = distributionData.pathogenic?.total || 0;
  apiClassifications[apiId].likely_pathogenic = distributionData.likely_pathogenic?.total || 0;
  apiClassifications[apiId].conflicting = distributionData.conflicting?.total || 0;
  apiClassifications[apiId].vus = distributionData.vus?.total || 0;
  apiClassifications[apiId].likely_benign = distributionData.likely_benign?.total || 0;
  apiClassifications[apiId].benign = distributionData.benign?.total || 0;

  // Calculate color
  apiClassifications[apiId].color = calculateClassificationColor(apiClassifications[apiId]);
  updateAPIBadge(apiId);
}

/**
 * Initialize all API badges on page load with grouped layout.
 */
function initializeAPIBadges() {
  const container = document.getElementById('api-badges-container');
  if (!container) return;

  // Group APIs by category
  const groupedAPIs = {};
  Object.keys(API_VERSIONS).forEach(apiId => {
    const config = API_VERSIONS[apiId];
    const category = config.category || 'Other';
    if (!groupedAPIs[category]) {
      groupedAPIs[category] = [];
    }
    groupedAPIs[category].push(apiId);
  });

  // Render groups in order
  const categoryOrder = ['Variant Annotation', 'Functional Data', 'Gene Context'];
  container.innerHTML = '';

  categoryOrder.forEach(category => {
    if (!groupedAPIs[category] || groupedAPIs[category].length === 0) return;

    const group = document.createElement('div');
    group.className = 'api-group';

    const label = document.createElement('div');
    label.className = 'api-group-label';
    label.textContent = category;
    group.appendChild(label);

    const badgesDiv = document.createElement('div');
    badgesDiv.className = 'api-group-badges';

    groupedAPIs[category].forEach(apiId => {
      const badge = document.createElement('div');
      badge.id = `api-badge-${apiId}`;
      badgesDiv.appendChild(badge);
      updateAPIBadge(apiId);
    });

    group.appendChild(badgesDiv);
    container.appendChild(group);
  });

  // Render any remaining categories
  categoryOrder.forEach(cat => delete groupedAPIs[cat]);
  Object.keys(groupedAPIs).forEach(category => {
    const group = document.createElement('div');
    group.className = 'api-group';

    const label = document.createElement('div');
    label.className = 'api-group-label';
    label.textContent = category;
    group.appendChild(label);

    const badgesDiv = document.createElement('div');
    badgesDiv.className = 'api-group-badges';

    groupedAPIs[category].forEach(apiId => {
      const badge = document.createElement('div');
      badge.id = `api-badge-${apiId}`;
      badgesDiv.appendChild(badge);
      updateAPIBadge(apiId);
    });

    group.appendChild(badgesDiv);
    container.appendChild(group);
  });
}

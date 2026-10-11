#!/bin/bash
# example_run.sh
# Simple example demonstrating the ReAnnota pipeline
#
# This script shows how to enhance genome annotations by combining:
# - Bakta initial annotation (GBFF format)
# - EggNOG functional annotations (TSV format)
# - InterProScan domain annotations (GFF3 format)
#
# Usage:
#   cd ReAnnota
#   ./examples/example_run.sh

set -e  # Exit immediately if a command exits with a non-zero status
set -u  # Treat unset variables as an error

# ============================================================================
# Configuration
# ============================================================================

# Get the project root directory (parent of examples/)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Input files (from test fixtures)
TEST_FIXTURES="${PROJECT_ROOT}/tests/fixtures"
BAKTA_GBFF="${TEST_FIXTURES}/bakta_test.gbff"
EGGNOG_TSV="${TEST_FIXTURES}/eggnogmapper_test.annotations"
INTERPRO_GFF="${TEST_FIXTURES}/interproscan_test.gff3"
GFF_INPUT="${TEST_FIXTURES}/bakta_test.gff3"
ANTISMASH_OUTPUT="${TEST_FIXTURES}/antismash_test.gbk"
PSEUDOFINDER_INPUT="${TEST_FIXTURES}/pseudofinder_test.gff"
GECCO_INPUT="${TEST_FIXTURES}/gecco_gbks_test.csv"
GECCO_SUMMARY="${TEST_FIXTURES}/gecco_summary_test.clusters.tsv"

# Output directory
OUTPUT_DIR="${PROJECT_ROOT}/examples/output"
mkdir -p "${OUTPUT_DIR}"

# GO Database
GO_OBO="${OUTPUT_DIR}/go-basic.obo"

echo "============================================================"
echo "  ReAnnota Example Pipeline"
echo "============================================================"
echo ""
echo "Project root: ${PROJECT_ROOT}"
echo "Test data:    ${TEST_FIXTURES}"
echo "Output:       ${OUTPUT_DIR}"
echo ""

# ============================================================================
# Step 1: Check Prerequisites
# ============================================================================

echo "=== Step 1: Checking prerequisites ==="
echo ""

# Check if reannota command is available
if ! command -v reannota &> /dev/null; then
    echo "❌ ERROR: 'reannota' command not found!"
    echo ""
    echo "The package is not installed. Please install it first:"
    echo ""
    echo "  cd ${PROJECT_ROOT}"
    echo "  uv venv"
    echo "  source .venv/bin/activate"
    echo "  uv pip install -e ."
    echo ""
    exit 1
fi

echo "✅ ReAnnota command found"

# Check if input files exist
missing_files=0
for file in "${BAKTA_GBFF}" "${EGGNOG_TSV}" "${INTERPRO_GFF}"; do
    if [[ ! -f "${file}" ]]; then
        echo "❌ Missing input file: ${file}"
        missing_files=1
    else
        echo "✅ Found: $(basename "${file}")"
    fi
done

if [[ ${missing_files} -eq 1 ]]; then
    echo ""
    echo "❌ ERROR: Some input files are missing!"
    exit 1
fi

# Check for GO database, download if missing
if [[ ! -f "${GO_OBO}" ]]; then
    echo "📥 Downloading go-basic.obo for the GO DAG network..."
    wget -q -O "${GO_OBO}" http://purl.obolibrary.org/obo/go/go-basic.obo || curl -s -o "${GO_OBO}" http://purl.obolibrary.org/obo/go/go-basic.obo
    echo "✅ Downloaded go-basic.obo"
else
    echo "✅ Found: go-basic.obo"
fi

echo ""
echo "✅ All prerequisites met!"
echo ""

# ============================================================================
# Step 2: Run the Annotation Enhancement Pipeline
# ============================================================================

echo "=== Step 2: Running annotation enhancement ==="
echo ""
echo "This will:"
echo "  1. Parse EggNOG annotations → extract GO, COG, KEGG, PFAM terms"
echo "  2. Parse InterPro annotations → extract InterPro IDs, GO terms"
echo "  3. Merge both into the Bakta GBFF file using maximal GO specificity"
echo "  4. Convert enhanced GBFF to GFF3 format"
echo ""

# Run the main pipeline
reannota annotate \
    --egg-input "${EGGNOG_TSV}" \
    --ipr-input "${INTERPRO_GFF}" \
    --gbff-input "${BAKTA_GBFF}" \
    --gff-input "${GFF_INPUT}" \
    --output "${OUTPUT_DIR}" \
    --antismash-input "${ANTISMASH_OUTPUT}" \
    --pseudofinder-input "${PSEUDOFINDER_INPUT}" \
    --gecco-input "${GECCO_INPUT}" \
    --gecco-summary "${GECCO_SUMMARY}" \
    --go-obo "${GO_OBO}" \
    --prefix test \
    --compare \
    --circos

echo ""
echo "✅ Annotation enhancement complete!"
echo ""

# ============================================================================
# Step 3: Examine Results
# ============================================================================

echo "=== Step 3: Results Summary ==="
echo ""

# List all output files
echo "📁 Output files created:"
ls -lh "${OUTPUT_DIR}/" | tail -n +2 | awk '{printf "  - %-30s %10s\n", $9, $5}'
echo ""

# Analyze the enhanced GFF if it exists
if [[ -f "${OUTPUT_DIR}/enhanced/test_enhanced.gff3" ]]; then
    echo "📊 Enhanced GFF3 Statistics:"
    echo ""
    
    total_features=$(grep -v '^#' "${OUTPUT_DIR}/enhanced/test_enhanced.gff3" | wc -l | tr -d ' ')
    cds_features=$(grep -v '^#' "${OUTPUT_DIR}/enhanced/test_enhanced.gff3" | grep -c $'\tCDS\t' || echo "0")
    go_terms=$(grep -v '^#' "${OUTPUT_DIR}/enhanced/test_enhanced.gff3" | grep -c 'GO:' || echo "0")
    interpro_terms=$(grep -v '^#' "${OUTPUT_DIR}/enhanced/test_enhanced.gff3" | grep -c 'InterPro:' || echo "0")
    kegg_terms=$(grep -v '^#' "${OUTPUT_DIR}/enhanced/test_enhanced.gff3" | grep -c 'KEGG:' || echo "0")
    cog_terms=$(grep -v '^#' "${OUTPUT_DIR}/enhanced/test_enhanced.gff3" | grep -c 'COG:' || echo "0")
    hypotheticals=$(grep -v '^#' "${OUTPUT_DIR}/enhanced/test_enhanced.gff3" | grep -c 'product=hypothetical' || echo "0")
    
    echo "  Total features:        ${total_features}"
    echo "  CDS features:          ${cds_features}"
    echo ""
    echo "  Features with GO:      ${go_terms}"
    echo "  Features with InterPro: ${interpro_terms}"
    echo "  Features with KEGG:    ${kegg_terms}"
    echo "  Features with COG:     ${cog_terms}"
    echo ""
    echo "  Hypothetical proteins: ${hypotheticals}"
    echo ""
fi

# Show top 5 enhanced genes
if [[ -f "${OUTPUT_DIR}/tool_hits/eggnog_hits.tsv" ]]; then
    echo "🧬 Sample EggNOG annotations (first 5 genes):"
    echo ""
    head -6 "${OUTPUT_DIR}/tool_hits/eggnog_hits.tsv" | column -t -s $'\t'
    echo ""
fi

# ============================================================================
# Summary
# ============================================================================

echo "============================================================"
echo "  ✅ Example Pipeline Complete!"
echo "============================================================"
echo ""
echo "📁 Results are in: ${OUTPUT_DIR}/"
echo ""
echo "📄 Key output files:"
echo "  - enhanced/test_enhanced.gbff  → Enhanced GenBank file"
echo "  - enhanced/test_enhanced.gff3  → Enhanced GFF3 file"
echo "  - tool_hits/eggnog_hits.tsv    → Parsed EggNOG annotations"
echo "  - tool_hits/ipr_hits.tsv       → Parsed InterPro annotations"
echo ""
echo "🎯 Next steps:"
echo "  1. Examine the enhanced files in a genome browser"
echo "  2. Compare with original annotation"
echo "  3. Run with your own data:"
echo ""
echo "     reannota annotate \\"
echo "       --egg-input your_eggnog.tsv \\"
echo "       --ipr-input your_interpro.gff3 \\"
echo "       --gbff-input your_bakta.gbff \\"
echo "       --go-obo path/to/go-basic.obo \\"
echo "       --output your_output_dir"
echo ""
echo "============================================================"
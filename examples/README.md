# Cyano-Annotation Examples

This directory contains example scripts demonstrating how to use the cyano-annotation package.

## Quick Start

```bash
# Make sure you've installed the package
cd /path/to/cyano-annotation
source .venv/bin/activate

# Run the example
./examples/example_run.sh
```

## What It Does

The `example_run.sh` script demonstrates the complete annotation enhancement workflow:

1. **Checks prerequisites** - Verifies cyano-annotate is installed and input files exist
2. **Runs the pipeline** - Enhances annotations using EggNOG and InterPro data
3. **Shows results** - Displays statistics about the enhanced annotations

## Input Files

The example uses test data from `tests/fixtures/`:
- `Bakta_anotation.gbff` - Initial genome annotation from Bakta
- `Egg_NOG_anotation.tsv` - Functional annotations from EggNOG-mapper
- `Interpro_anotation.gff3` - Domain annotations from InterProScan
- `Proclorococus_genome.fasta` - Reference genome (Prochlorococcus)

## Output Files

Results are saved to `examples/output/`:
- `enhanced_annotation.gbff` - Enhanced GenBank file with merged annotations
- `enhanced.gff` - Enhanced GFF3 file for genome browsers
- `eggnog_hits.tsv` - Parsed and filtered EggNOG annotations
- `ipr_hits.tsv` - Parsed and filtered InterPro annotations

## Expected Results

For the Prochlorococcus test genome, you should see:
- ~1,800+ CDS features
- Significant increase in GO terms
- Reduction in hypothetical proteins
- Added KEGG, COG, and PFAM annotations

## Using Your Own Data

To run with your own genome:

```bash
cyano-annotate \
  --egg_input path/to/your_eggnog.emapper.annotations \
  --ipr_input path/to/your_interpro.gff3 \
  --gbff_input path/to/your_bakta.gbff \
  --output_file path/to/enhanced_output.gbff
```

## Troubleshooting

### "cyano-annotate: command not found"

Make sure you've installed the package and activated the virtual environment:

```bash
cd /path/to/cyano-annotation
source .venv/bin/activate
uv pip install -e .
```

### "Missing input files"

The test fixtures should be in `tests/fixtures/`. If they're missing:

```bash
# Check if they exist
ls tests/fixtures/

# They should have been copied from data/test_runs/Prochlorococus_run/
```

### Permission denied

Make the script executable:

```bash
chmod +x examples/example_run.sh
```

## Next Steps

- Examine the enhanced files in a genome browser (e.g., IGV, Artemis)
- Compare annotation quality before and after enhancement
- Run benchmarking analysis with `cyano-benchmark` (if available)
- Integrate into your own analysis pipelines


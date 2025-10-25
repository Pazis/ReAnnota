# Test Fixtures

This directory contains test data files used for both **examples** and **automated tests** (pytest).

## Files

### Input Files (for annotation enhancement)

- **`Proclorococus_genome.fasta`** - Prochlorococcus marinus genome (small cyanobacterium, ~1.7 Mb)
- **`Bakta_anotation.gbff`** - Initial genome annotation from Bakta (GenBank format)
- **`Egg_NOG_anotation.tsv`** - Functional annotations from EggNOG-mapper (raw output)
- **`Interpro_anotation.gff3`** - Domain annotations from InterProScan (GFF3 format)

## Data Source

These files were copied from `data/test_runs/Prochlorococus_run/` which contains:
- Original input files (now here in `tests/fixtures/`)
- Expected output files (remain in `data/test_runs/` for comparison)

## Usage

### In Examples

```bash
# The example script uses these files
./examples/example_run.sh
```

### In Tests (Future)

```python
# pytest will use these fixtures
def test_eggnog_parser(test_fixtures_dir):
    eggnog_file = test_fixtures_dir / "Egg_NOG_anotation.tsv"
    result = build_egg_dictionary_clean(str(eggnog_file))
    assert len(result) > 0
```

### In Your Own Tests

```python
from pathlib import Path

# Get path to fixtures
fixtures_dir = Path(__file__).parent / "fixtures"
bakta_file = fixtures_dir / "Bakta_anotation.gbff"
```

## File Formats

### Bakta GBFF
GenBank flat file format with:
- Sequence records
- CDS features with locus_tags
- Initial functional annotations
- Gene predictions

### EggNOG TSV
Tab-separated file with columns:
- query (gene ID)
- seed_ortholog
- evalue
- eggNOG_OGs (orthologous groups)
- COG_category
- Description
- GOs (Gene Ontology terms)
- KEGG_ko (KEGG orthology)
- PFAMs (Pfam domains)
- etc.

### InterPro GFF3
GFF3 format with:
- CDS features
- InterPro domain signatures
- GO terms
- Signature descriptions

## Size Information

```
Proclorococus_genome.fasta  ~1.7 MB (genome sequence)
Bakta_anotation.gbff        ~5-10 MB (annotations)
Egg_NOG_anotation.tsv       ~500 KB (functional data)
Interpro_anotation.gff3     ~2-3 MB (domain data)
```

## Notes

- These are **real biological data**, not synthetic test data
- Suitable for integration testing, not just unit testing
- Keep files small enough for version control (<10 MB each)
- For larger test datasets, store externally and download on-demand

## Reference Output

Expected output from these inputs can be found in:
```
data/test_runs/Prochlorococus_run/
├── final_enchanced.gbff    # Expected enhanced GenBank output
├── final_enchanced.gff3    # Expected enhanced GFF3 output
└── README.md               # Details about this test run
```

These can be used for comparison in tests to verify correct pipeline behavior.


# TRMotifAnnotator

A command-line tool for identifying, annotating, and visualizing motifs in tandem repeats. Detects canonical and variant repeat motifs including substitutions, insertions, and deletions.

---

## Features

- **Motif Detection**: Identifies canonical and non-canonical repeat motifs with single-base resolution
- **Multiple Canonical Motifs**: Support for repeats with multiple motif variants (e.g., `CAG,CAA` for polyglutamine)
- **Degenerate Bases**: IUPAC codes automatically expand (e.g., `GCN` → `GCA,GCC,GCG,GCT`)
- **Rotation Handling**: Configurable treatment of motif rotations (coding vs. non-coding repeats)
- **Population Stratification**: Optional visualization stratified by superpopulation
- **Statistical Output**: Comprehensive TSV with motif counts, percentages, and structures
- **Visualization**: Color-coded sequence composition plots with customizable legends
- **HTML Reports**: Reports with embedded plots and data tables

---

## Installation

**Requirements:**
- Python 3.6+
- NumPy
- Matplotlib

```bash
# Install dependencies
pip install numpy matplotlib

# Clone repository
git clone https://github.com/wf-TRs/TRMotifAnnotator.git
cd TRMotifAnnotator

# Verify installation
python TRMotifAnnotator.py --help
```

---

## Quick Start

```bash
python TRMotifAnnotator.py \
  --input RFC1_sequences.fasta \
  --output RFC1_analysis \
  --canonical-motif AAAAG \
  --max-mers 5 \
  --locus RFC1
```

**Output:**
- `RFC1_analysis.tsv` - Detailed motif annotations
- `RFC1_analysis.png` - Sequence composition visualization

---

## Usage

### Basic Command

```bash
python TRMotifAnnotator.py \
  --input <input.fasta> \
  --output <output_prefix> \
  --canonical-motif <motif> \
  --max-mers <length> \
  --locus <locus_name>
```

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--input` | Input FASTA file | `RFC1.fasta` |
| `--output` | Output prefix (no extension) | `RFC1_results` |
| `--canonical-motif` | Canonical motif(s), comma-separated | `AAAAG` or `CAG,CAA` |
| `--max-mers` | Motif length | `5` |
| `--locus` | Locus name | `RFC1` |

### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--vlines` | Vertical reference lines: `"[(pos, 'color'), ...]"` | None |
| `--num-nc-motifs` | Number of non-canonical motifs in legend | `10` |
| `--non-canonical-motifs` | User-specified non-canonical motifs | Auto-detect |
| `--hide-y-label` | Hide sample IDs on y-axis | Show |
| `--sort-sequences` | Sort order: `input`, `largest`, `smallest` | `input` |
| `--treat-motif-rotations` | Treat rotations as canonical | False |
| `--html` | Generate HTML report | False |
| `--stratify-by-superpop` | Stratify visualization by superpopulation | False |
| `--superpop-file` | TSV mapping if auto-detection fails: `seq_id<TAB>superpop` | Auto-detect |

---

## Examples

### Example 1: Basic Analysis

```bash
python TRMotifAnnotator.py \
  --input RFC1.fasta \
  --output RFC1_basic \
  --canonical-motif AAAAG \
  --max-mers 5 \
  --locus RFC1
```

### Example 2: Multiple Canonical Motifs (Polyglutamine)

```bash
python TRMotifAnnotator.py \
  --input HTT.fasta \
  --output HTT_polyQ \
  --canonical-motif CAG,CAA \
  --max-mers 3 \
  --locus HTT
```

Both `CAG` and `CAA` are treated as canonical and colored identically.

### Example 3: Degenerate Bases (Polyalanine)

```bash
python TRMotifAnnotator.py \
  --input PHOX2B.fasta \
  --output PHOX2B_polyA \
  --canonical-motif GCN \
  --max-mers 3 \
  --locus PHOX2B
```

`GCN` expands to `GCA`, `GCC`, `GCG`, `GCT` (all treated as canonical).

### Example 4: Rotation Handling (Non-Coding)

```bash
python TRMotifAnnotator.py \
  --input RFC1.fasta \
  --output RFC1_noncoding \
  --canonical-motif AAAAG \
  --max-mers 5 \
  --locus RFC1 \
  --treat-motif-rotations
```

All motif rotations (`AAAAG`, `AAAGA`, `AAGAA`, `AGAAA`, `GAAAA`) treated the same. This applies to both canonical and non-canonical motifs.

### Example 5: Allelic Boundaries & HTML Report

```bash
python TRMotifAnnotator.py \
  --input CNBP.fasta \
  --output CNBP_annotated \
  --canonical-motif CCTG \
  --max-mers 4 \
  --locus CNBP \
  --vlines "[(120, 'green'), (250, 'orange'), (400, 'red')]" \
  --num-nc-motifs 15 \
  --html
```

### Example 6: Population Stratification

```bash
python TRMotifAnnotator.py \
  --input cohort.fasta \
  --output cohort_stratified \
  --canonical-motif AAAAG \
  --max-mers 5 \
  --locus RFC1 \
  --stratify-by-superpop
```

**Automatic detection:** TRMotifAnnotator automatically extracts superpopulation codes from FASTA headers using common patterns:
- `[NA18519.p1_AFR]` or `[NA18519.p2_AFR]`
- `HG00321[EUR]` or `HG00321_EUR`
- `NA12878.h1_AFR` or `sample_EAS_CHB`

**Manual mapping (optional):** Provide a TSV file if automatic detection fails:
```bash
--superpop-file population_mapping.tsv
```

**Mapping file format:**
```
NA18519.p1	AFR
NA18519.p2	AFR
HG00321.p1	EUR
```

### Example 7: Large Dataset (Hidden Labels)

```bash
python TRMotifAnnotator.py \
  --input population_cohort.fasta \
  --output cohort_analysis \
  --canonical-motif AAAAG \
  --max-mers 5 \
  --locus RFC1 \
  --sort-sequences smallest \
  --hide-y-label \
  --num-nc-motifs 15
```

---

## Output Files

### TSV File (`<output>.tsv`)

Tab-delimited file with motif annotations:

| Column | Description |
|--------|-------------|
| `Sample ID` | Sequence identifier |
| `Sequence Length` | Total sequence length (bp) |
| `Repeat Copy` | Number of repeat units |
| `Canonical Motif` | Canonical motif(s) analyzed |
| `SortedNon-CanonicalMotifs` | Non-canonical motifs (order of appearance) |
| `Repeat Structure` | Detailed structure (e.g., `(AAAAG)5-AAAGG-(AAAAG)3`) |
| `Canonical Motif Count` | Number of canonical occurrences |
| `Non-Canonical Motif Count` | Number of non-canonical occurrences |
| `%Non-Canonical Base` | Percentage of non-canonical bases |
| `Non-Canonical Same Length` | Same-length variants |
| `Non-Canonical Different Length` | Different-length variants |

### Visualization (`<output>.png`)

High-resolution plot showing:
- **X-axis**: Sequence position
- **Y-axis**: Sample IDs (optional)
- **Colors**: Unique color per motif (canonical = light sky blue)
- **Legend**: Most common motifs (configurable)

### HTML Report (`<output>.html`)

Interactive report with:
- Input parameters and metadata
- Embedded visualization
- Scrollable TSV data table
- Links to output files

---

## Advanced Features

### Motif Rotation Handling

**Coding repeats** (rotations change amino acids):
```bash
--canonical-motif CAG  # Only exact CAG is canonical
```

**Non-coding repeats** (rotations are equivalent):
```bash
--canonical-motif AAAAG --treat-motif-rotations
```

### Custom Non-Canonical Motifs

```bash
--non-canonical-motifs AAAGG,AAAGA,AAGGG
```

### Locus-Specific Processing

Special handling for EIF4A3 (18-23mer motifs):
```bash
--canonical-motif CCTCGCTGCCGCTGCCGA \
--max-mers 18 \
--locus EIF4A3
```

---

## TRGT Integration

### Extracting Sequences from TRGT VCF Files

The included `TRGT-helper.py` script extracts allele sequences from [TRGT](https://github.com/PacificBiosciences/trgt) VCF files and converts them to FASTA format for analysis with TRMotifAnnotator.

**TRGT** (Tandem Repeat Genotyping Tool) genotypes tandem repeats from PacBio HiFi data. The helper script bridges TRGT output to motif-level analysis.

#### Basic Usage

```bash
python TRGT-helper.py <vcf> [--prefix PREFIX]
```

**Arguments:**
- `vcf` - TRGT VCF file (plain or gzipped)
- `--prefix` - Optional output prefix

**Output naming:**
- Single sample VCF: `SAMPLE_TRID.fasta`
- Multi-sample VCF: `TRID.fasta`
- Custom prefix: `PREFIX_TRID.fasta`

#### Examples

**Single Sample:**
```bash
python TRGT-helper.py HG002_trgt.vcf.gz
# Output: HG002_RFC1.fasta, HG002_HTT.fasta, etc.
```

**Multiple Samples:**
```bash
python TRGT-helper.py cohort_trgt.vcf.gz
# Output: RFC1.fasta, HTT.fasta, CNBP.fasta, etc.
```

**Custom Prefix:**
```bash
python TRGT-helper.py cohort_trgt.vcf.gz --prefix EUR_cohort
# Output: EUR_cohort_RFC1.fasta, EUR_cohort_HTT.fasta, etc.
```

### Complete TRGT → TRMotifAnnotator Pipeline

```bash
#!/bin/bash

SAMPLE="HG002"
LOCUS="RFC1"
CANONICAL="AAAAG"

# Step 1: Run TRGT genotyping
trgt \
  --genome reference.fasta \
  --repeats repeats.bed \
  --reads ${SAMPLE}.bam \
  --output-prefix ${SAMPLE}_trgt

# Step 2: Extract allele sequences
python TRGT-helper.py ${SAMPLE}_trgt.vcf.gz
# Output: HG002_RFC1.fasta (and other loci)

# Step 3: Analyze motifs
python TRMotifAnnotator.py \
  --input ${SAMPLE}_${LOCUS}.fasta \
  --output ${SAMPLE}_${LOCUS}_motifs \
  --canonical-motif ${CANONICAL} \
  --max-mers 5 \
  --locus ${LOCUS} \
  --treat-motif-rotations \
  --html

# Results:
# - HG002_RFC1_motifs.tsv
# - HG002_RFC1_motifs.png
# - HG002_RFC1_motifs.html
```

### FASTA Output Format

```
>HG002_RFC1_Allele-1
AAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAG
AAAAGAAAAGAAAAGAAAAG
>HG002_RFC1_Allele-2
AAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAGGAAAAGAAAAGAAAAG
AAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAG
```

**Header format:** `SAMPLE_TRID_Allele-N`

### Helper Script Features

- **Automatic detection**: Extracts all TRIDs from VCF
- **Gzip support**: Handles `.vcf` and `.vcf.gz` files
- **Phasing preserved**: Maintains allele order from GT field
- **Multi-locus**: Processes all tandem repeats in one run
- **No dependencies**: Uses only Python standard library

---

## Citation

If you use TRMotifAnnotator in your research, please cite:

```bibtex
@article{RajanBabu2025,
  title={Population-scale disease-associated tandem repeat analysis reveals locus and ancestry-specific insights},
  author={Rajan-Babu, Indhu-Shree and Chiu, Readman and Weisburd, Ben and Caglayan, Iris and Birol, Inanc and Friedman, Jan M},
  journal={medRxiv},
  year={2025},
  doi={10.1101/2025.10.11.25337795}
}
```

---

## Support

- **GitHub Issues**: [github.com/wf-TRs/TRMotifAnnotator/issues](https://github.com/wf-TRs/TRMotifAnnotator/issues)
- **Contact**: Indhu Shree Rajan Babu - indhu.babu@bcchr.ca

---

## License

MIT License - see [LICENSE](LICENSE) file for details.

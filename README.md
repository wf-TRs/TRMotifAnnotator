# TRMotifAnnotator

## Overview

TRMotifAnnotator is a command-line tool designed for the identification and annotation of tandem repeat (TR) motifs within disease-associated loci. The tool detects both **canonical** and **non-canonical** repeat motifs, including those resulting from **substitutions** or **insertion/deletion (indel) events**. The output includes a **TSV file** with detailed annotations of repeat genotypes and motif structures, along with a **sequence composition plot** highlighting canonical and non-canonical motifs.

## Features

- Identifies and annotates **repeat motifs** in given sequences
- Detects **non-canonical motifs** caused by substitutions or indels
- Outputs a **TSV file** with motif details
- Generates **sequence composition plots**
- Allows **user-defined allelic class thresholds** for visualization

## Installation

### Dependencies

TRMotifAnnotator requires the following:

- Python 3.x
- NumPy
- Matplotlib
- Pandas
- Biopython

Install dependencies using:

```bash
pip install numpy matplotlib pandas biopython
```

### Cloning the Repository

To install TRMotifAnnotator, clone the GitHub repository:

```bash
git clone https://github.com/your-repo/TRMotifAnnotator.git
cd TRSeqExplorer
```

## Usage

Run the script with the following command:

```bash
python TRMotifAnnotator.py --input <sequence.fa> --output <prefix> --canonical-motif <motif> --max-mers <motif_length> --vlines "[(value1, 'color1'), (value2, 'color2')]" --locus <locus-name>
```

### Arguments

- `--input <sequence.fa>`: Input FASTA file containing repeat sequences
- `--output <prefix>`: Prefix for output files
- `--canonical-motif <motif>`: Canonical repeat motif used for comparison
- `--max-mers <motif_length>`: Maximum motif length to consider
- `--vlines "[(value1, 'color1'), (value2, 'color2')]"`: User-defined thresholds for allelic classes, used for visualization
- `--locus <locus-name>`: Name of the locus being analyzed

### Example

```bash
python TRMotifAnnotator.py --input example.fa --output results --canonical-motif CAG --max-mers 3 --vlines "[(30, 'gray'), (50, 'red')]" --locus HTT
```

## Output

- **TSV File (**``**)**: Contains motif annotations and structural details
- **Sequence Composition Plot (**``**)**: Graphical representation of repeat motifs

## Citation

If you use TRMotifAnnotator in your research, please cite:

> Indhu-Shree Rajan Babu, Readman Chiu, Iris Caglayan, Inanc Birol, Jan M. Friedman. In-Depth Characterisation of Disease-Associated Tandem Repeat Loci and Their Local Ancestries and Haplotypes in Diverse Human Populations.

## Contact

For issues or contributions, open a GitHub issue or contact: **Your Name** - [indhu.babu@bcchr.ca](mailto\:indhu.babu@bcchr.ca)

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

- Python 3.6+
- NumPy
- Matplotlib

Install dependencies using:

```bash
pip install numpy matplotlib
```

### Cloning the Repository

To install TRMotifAnnotator, clone the GitHub repository:

```bash
git clone https://github.com/wf-TRs/TRMotifAnnotator.git
cd TRMotifAnnotator
```

## Usage

Run the script with the following command:

```bash
python TRMotifAnnotator.py --input <sequence.fa> --output <prefix> --canonical-motif <motif> --max-mers <motif_length> --vlines "[(value1, 'color1'), (value2, 'color2')]" --locus <locus-name>
```

### Arguments

`--input <sequence.fa>`: Input FASTA file containing repeat sequences  
`--output <prefix>`: Prefix for output files  
`--canonical-motif <motif>`: Canonical repeat motif  
`--max-mers <motif_length>`: Canonical motif length (e.g., 5 if the motif is AAAAG)  
`--vlines "[(value1, 'color1'), (value2, 'color2')]"`: User-defined thresholds for allelic classes  
`--locus <locus-name>`: Name of the locus being analyzed - will be displayed as the plot title  
`--num-nc-motifs <Non-canonical_motif-number>`: Number of non-canonical motifs to be displayed; by default up to 10 motifs will be displayed  


### Example Command

```bash
python TRMotifAnnotator.py --input example.fa --output results --canonical-motif CAG --max-mers 3 --vlines "[(105, 'red')]" --locus ATXN1
```

## Output

- **TSV File:** Contains detailed motif annotations and structural information.

### Example TSV File

| Sample ID        | Sequence Length | Repeat Copy | Canonical Motif | Sorted Non-Canonical Motifs | Repeat Structure                | Canonical Motif Count | Non-Canonical Motif Count | % Non-Canonical Base | Non-Canonical Same Length | Non-Canonical Different Length |
|------------------|-----------------|-------------|-----------------|-----------------------------|--------------------------------|----------------------|---------------------------|----------------------|---------------------------|-------------------------------|
| HG00101_EUR_GBR  | 111             | 37          | CAG             | CAT                         | (CAG)18-CAT-CAG-CAT-(CAG)16   | 35                   | 2                         | 5.41                 | 2                         | 0                             |
| HG00131_EUR_GBR  | 108             | 36          | CAG             | CAT                         | (CAG)17-CAT-CAG-CAT-(CAG)16   | 34                   | 2                         | 5.56                 | 2                         | 0                             |
| HG00261_EUR_GBR  | 108             | 36          | CAG             | CAT                         | (CAG)17-CAT-CAG-CAT-(CAG)16   | 34                   | 2                         | 5.56                 | 2                         | 0                             |
| HG00360_EUR_FIN  | 108             | 36          | CAG             | CAT                         | (CAG)17-CAT-CAG-CAT-(CAG)16   | 34                   | 2                         | 5.56                 | 2                         | 0                             |

The TSV file includes:
- Sequence lengths (in base pairs)  
- Repeat copy numbers  
- Canonical motifs  
- Non-canonical motifs (in the order they occur within the repeat tract)  
- Detailed repeat structure  
- Counts of canonical and non-canonical motifs  
- Percentage of non-canonical bases (proportion of bases from non-canonical motifs)  
- Counts of non-canonical motifs of the same length as the canonical motif or of different length  

- **Sequence Composition Plot:** Graphical representation of repeat motifs of the four sequences above.
![Sequence Composition Plot](examples/ATXN1_sequence_composition.png)

## Citation

If you use TRMotifAnnotator in your research, please cite:

> Indhu-Shree Rajan Babu, Readman Chiu et al. In-Depth Characterisation of Disease-Associated Tandem Repeat Loci and Their Local Ancestries and Haplotypes in Diverse Human Populations.

## Contact

For issues or contributions, open a GitHub issue or contact: **Your Name** - [indhu.babu@bcchr.ca](mailto\:indhu.babu@bcchr.ca)

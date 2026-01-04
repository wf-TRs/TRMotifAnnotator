# FGF14 Example: Effect of Rotation Handling

This example demonstrates the impact of the `--treat-motif-rotations` flag on motif classification for the FGF14 locus.

## Background

**Locus:** FGF14  
**Canonical motif:** GAA  
**Pathogenic threshold:** Varies by study  
**Reference lines:**
- Gray (179 copies): Normal upper limit
- Orange (249 copies): Intermediate upper limit
- Red (335 copies): Pathogenic threshold

## Commands

### With Rotation Treatment

```bash
python TRMotifAnnotator.py \
  --input FGF14_expansions.fasta \
  --output FGF14-treat-rotations \
  --canonical-motif GAA \
  --vlines "[(179, 'gray'), (249,'orange'), (335, 'red')]" \
  --locus FGF14 \
  --sort-sequences smallest \
  --num-nc-motifs 10 \
  --html \
  --treat-motif-rotation
```

**Result:** GAA rotations (AAG, AGA) are treated as the same motif.

### Without Rotation Treatment (Default)

```bash
python TRMotifAnnotator.py \
  --input FGF14_expansions.fasta \
  --output FGF14-no-rotations \
  --canonical-motif GAA \
  --vlines "[(179, 'gray'), (249,'orange'), (335, 'red')]" \
  --locus FGF14 \
  --sort-sequences smallest \
  --num-nc-motifs 10 \
  --html
```

**Result:** Only exact GAA matches are canonical. AAG and AGA rotations are treated as non-canonical variants.

## Output Files

### With Rotations Treated as Canonical
- `FGF14-treat-rotations.png` - Visualization
- `FGF14-treat-rotations.tsv` - Detailed annotations
- `FGF14-treat-rotations.html` - Report displaying plot and TSV data

### Without Rotation Treatment
- `FGF14-no-rotations.png` - Visualization
- `FGF14-no-rotations.tsv` - Detailed annotations
- `FGF14-no-rotations.html` - Report displaying plot and TSV data

## When to Use Each Approach

**Use `--treat-motif-rotations`:**
- Coding sequences where rotations encode the same amino acid (e.g; CAG, CAA)
- Analyses focusing on amino acid composition
- Non-coding sequences

**Don't use `--treat-motif-rotations`:**
- Coding sequences where rotations do not encode the same amino acid

import sys
import os
import ast
import argparse
from typing import Dict, List, Tuple
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend to avoid X11 warnings
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import logging
logging.basicConfig(level=logging.INFO)
from matplotlib.font_manager import FontProperties
from collections import OrderedDict

# Constants
CANONICAL_COLOR = 'lightskyblue'
RANDOM_SEED = 42

# IUPAC degenerate base codes
IUPAC_CODES = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
}

def expand_degenerate_motif(motif: str) -> List[str]:
    """
    Expand degenerate motif to all possible sequences.
    Example: GCN -> [GCA, GCC, GCG, GCT]
    """
    if not motif:
        return []
    
    # Start with a list containing empty string
    expanded = ['']
    
    for base in motif.upper():
        if base not in IUPAC_CODES:
            # Invalid base, treat as literal
            expanded = [seq + base for seq in expanded]
        else:
            # Expand this position
            new_expanded = []
            for seq in expanded:
                for possible_base in IUPAC_CODES[base]:
                    new_expanded.append(seq + possible_base)
            expanded = new_expanded
    
    return expanded

def parse_canonical_motifs(canonical_motif_str: str) -> Tuple[List[str], bool]:
    """
    Parse canonical motif string, handling multiple motifs and degenerate bases.
    Returns: (list of expanded motifs, whether any degenerate bases were used)
    
    Examples:
        "CAG" -> (["CAG"], False)
        "CAG,CAA" -> (["CAG", "CAA"], False)
        "GCN" -> (["GCA", "GCC", "GCG", "GCT"], True)
        "CAG,GCN" -> (["CAG", "GCA", "GCC", "GCG", "GCT"], True)
    """
    motifs = []
    has_degenerate = False
    
    for motif in canonical_motif_str.split(','):
        motif = motif.strip()
        if not motif:
            continue
        
        # Check if this motif has degenerate bases
        has_deg = any(base in motif.upper() for base in IUPAC_CODES.keys() if base not in 'ACGT')
        if has_deg:
            has_degenerate = True
        
        # Expand degenerate bases
        expanded = expand_degenerate_motif(motif)
        motifs.extend(expanded)
    
    # Remove duplicates while preserving order
    seen = set()
    unique_motifs = []
    for m in motifs:
        if m not in seen:
            seen.add(m)
            unique_motifs.append(m)
    
    return unique_motifs, has_degenerate

DISTINCT_COLORS = OrderedDict([
    ("Orange",        "#FF8C00"), ("Magenta",       "#FF00FF"), ("Lime Green",    "#32CD32"),
    ("Gold",          "#FFD700"), ("Vermillion",    "#E34234"), ("Purple",        "#9400D3"),
    ("Teal",          "#008080"), ("Hot Pink",      "#FF69B4"), ("Navy",          "#000080"),
    ("Chocolate",     "#8B4513"), ("Crimson",       "#DC143C"), ("Cyan",          "#00FFFF"),
    ("Olive",         "#6B8E23"), ("Deep Pink",     "#FF1493"), ("Sienna",        "#A0522D"),
    ("Royal Blue",    "#4169E1"), ("Forest Green",  "#228B22"), ("Coral",         "#FF7F50"),
    ("Maroon",        "#800000"), ("Orchid",        "#DA70D6"), ("Slate Gray",    "#708090"),
    ("Tomato",        "#FF6347"), ("Indigo",        "#4B0082"), ("Tan",           "#D2B48C"),
])

def read_fasta(file_path: str) -> Dict[str, str]:
    """Read a FASTA file and return a dictionary of sequences."""
    sequences = {}
    try:
        with open(file_path, 'r') as file:
            seq_id, seq_data = "", ""
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if seq_id and seq_data:
                        sequences[seq_id] = seq_data
                    seq_id, seq_data = line[1:], ""
                else:
                    seq_data += line
            if seq_id and seq_data:
                sequences[seq_id] = seq_data
    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        sys.exit(1)
    return sequences

def extract_superpop_from_id(seq_id: str) -> str:
    """
    Try to extract superpopulation code from sequence ID using common patterns.
    
    Tries multiple patterns in order:
    1. [ID.p1_SUPERPOP] or [ID.p2_SUPERPOP] - e.g., [NA18519.p1_AFR]
    2. ID[SUPERPOP] or ID_POP[SUPERPOP] - e.g., HG00321[EUR]
    3. ID_SUPERPOP_... - e.g., HG00321_EUR_FIN
    
    Returns superpop code or "Unknown" if no pattern matches.
    """
    # Pattern 1: Look for [SAMPLE.p1_SUPERPOP] or [SAMPLE.p2_SUPERPOP]
    if '[' in seq_id and ']' in seq_id:
        try:
            # Extract content between [ and ]
            bracket_content = seq_id.split('[')[-1].split(']')[0].strip()
            
            # Check for underscore in bracket (e.g., NA18519.p1_AFR)
            if '_' in bracket_content:
                superpop = bracket_content.split('_')[-1].strip()
                if superpop and superpop in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
                    return superpop
                elif superpop and 2 <= len(superpop) <= 5 and superpop.isalpha() and superpop.isupper():
                    return superpop
            
            # Fallback: bracket content might be just superpop code (e.g., [EUR])
            elif bracket_content and bracket_content in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
                return bracket_content
            elif bracket_content and 2 <= len(bracket_content) <= 5 and bracket_content.isalpha() and bracket_content.isupper():
                return bracket_content
                
        except (IndexError, ValueError):
            pass
    
    # Pattern 2: Look for _SUPERPOP_ pattern (3-letter code after underscore)
    parts = seq_id.split('_')
    for part in parts:
        if part in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
            return part
        if part and 2 <= len(part) <= 5 and part.isupper() and part.isalpha():
            if len(part) == 3:  # Most superpop codes are 3 letters
                return part
    
    # Pattern 3: Look for .SUPERPOP. pattern
    if '.' in seq_id:
        parts = seq_id.split('.')
        for part in parts:
            if part in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
                return part
    
    return "Unknown"

def read_superpop_mapping(file_path: str) -> Dict[str, str]:
    """
    Read superpopulation mapping from TSV file.
    Format: seq_id<TAB>superpop
    """
    mapping = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    seq_id = parts[0].strip()
                    superpop = parts[1].strip()
                    mapping[seq_id] = superpop
    except FileNotFoundError:
        logging.error(f"Superpopulation mapping file not found: {file_path}")
        sys.exit(1)
    return mapping

def count_non_canonical_motifs(sequence: str, canonical_motifs: List[str], treat_rotations_as_canonical: bool = False) -> List[str]:
    """Determine counts of non-canonical motifs. Canonical_motifs is now a list."""
    non_canonical_motifs = []
    start_index = 0
    
    if treat_rotations_as_canonical:
        # Get all rotations of all canonical motifs
        canonical_set = set()
        for canonical_motif in canonical_motifs:
            for i in range(len(canonical_motif)):
                canonical_set.add(canonical_motif[i:] + canonical_motif[:i])
    else:
        # Just use exact matches
        canonical_set = set(canonical_motifs)
    
    while start_index < len(sequence):
        # Try to find any canonical motif (or rotation)
        found_canonical = False
        for canonical in canonical_set:
            if sequence.startswith(canonical, start_index):
                # Found a canonical motif, skip it
                start_index += len(canonical)
                found_canonical = True
                break
        
        if found_canonical:
            continue
        
        # Not a canonical motif, so it's non-canonical
        # Find where the next canonical motif starts
        next_canonical_index = len(sequence)
        for canonical in canonical_set:
            next_index = sequence.find(canonical, start_index)
            if next_index != -1:
                next_canonical_index = min(next_canonical_index, next_index)
        
        if next_canonical_index > start_index:
            non_canonical_motifs.append(sequence[start_index:next_canonical_index])
            start_index = next_canonical_index
        else:
            break
    
    return [motif for motif in non_canonical_motifs if motif]

def calculate_non_canonical_motifs(sequences: Dict[str, str], canonical_motifs: List[str], treat_rotations_as_canonical: bool = False) -> str:
    """Calculate non-canonical motifs for all sequences."""
    non_can_motifs = []
    for sequence in sequences.values():
        non_can_motifs.extend(count_non_canonical_motifs(sequence, canonical_motifs, treat_rotations_as_canonical))
    return ','.join(set(non_can_motifs))

def format_repeats(sequence: str, canonical_motifs: List[str], submotifs: List[str], locus: str = "") -> str:
    """Format the sequence to show repeated motifs, including submotifs."""
    formatted = []
    i = 0
    length = len(sequence)
    
    # Get max canonical motif length for lookahead
    max_canonical_len = max(len(cm) for cm in canonical_motifs) if canonical_motifs else 0
    
    # Check if all canonical motifs have the same length (reading frame matters)
    # This applies to triplets (3bp), pentamers (5bp), or any fixed-length repeat
    if canonical_motifs:
        canonical_lengths = set(len(cm) for cm in canonical_motifs)
        is_fixed_length_repeat = len(canonical_lengths) == 1
        frame_size = list(canonical_lengths)[0] if is_fixed_length_repeat else 0
    else:
        is_fixed_length_repeat = False
        frame_size = 0
    
    while i < length:
        # Check for any canonical motif first
        canonical_found = False
        for canonical_motif in canonical_motifs:
            if sequence[i:i+len(canonical_motif)] == canonical_motif:
                count = 1
                j = i + len(canonical_motif)
                while j + len(canonical_motif) <= length and sequence[j:j+len(canonical_motif)] == canonical_motif:
                    count += 1
                    j += len(canonical_motif)
                formatted.append(f"({canonical_motif}){count}" if count > 1 else canonical_motif)
                i = j
                canonical_found = True
                break
        
        if canonical_found:
            continue
        
        # BEFORE checking submotifs, check if canonical motif appears nearby
        # This prevents submotifs from breaking up partial repeats (e.g., AAA before AAAAG)
        lookahead_limit = min(max_canonical_len, length - i)
        canonical_nearby = False
        for look_dist in range(1, lookahead_limit + 1):
            for canonical_motif in canonical_motifs:
                if sequence[i+look_dist:i+look_dist+len(canonical_motif)] == canonical_motif:
                    # Found canonical motif nearby! Take chunk before it
                    chunk = sequence[i:i+look_dist]
                    formatted.append(chunk)
                    i += look_dist
                    canonical_nearby = True
                    break
            if canonical_nearby:
                break
        
        if canonical_nearby:
            continue
        
        # For fixed-length repeats, maintain reading frame to avoid frame shifts
        # Works for triplets (3bp), pentamers (5bp), or any fixed-length repeat
        if is_fixed_length_repeat and i + frame_size <= length:
            # Take a frame_size chunk (even if not at perfect frame boundary)
            chunk = sequence[i:i+frame_size]
            # Check for consecutive repeats of this chunk
            count = 1
            j = i + frame_size
            while j + frame_size <= length and sequence[j:j+frame_size] == chunk:
                count += 1
                j += frame_size
            formatted.append(f"({chunk}){count}" if count > 1 else chunk)
            i = j
            continue
        
        # Now check for submotifs (only if not a triplet at frame boundary)
        submotif_found = False
        for submotif in sorted(submotifs, key=len, reverse=True):
            if sequence.startswith(submotif, i):
                formatted.append(submotif)
                i += len(submotif)
                submotif_found = True
                break
        if submotif_found:
            continue
        
        # Collapse homopolymer runs (only if no motif matched and no canonical nearby)
        base = sequence[i]
        if base in "ACGT":
            j = i
            while j < length and sequence[j] == base:
                j += 1
            run_len = j - i
            if run_len >= 2:
                formatted.append(f"({base}){run_len}")
                i = j
                continue
        # Fallback: single base
        formatted.append(sequence[i])
        i += 1
    
    return '-'.join(formatted)

def post_process_repeat_structure(repeat_structure: str, non_canonical_motifs: str) -> str:
    """Post process repeat structures to consolidate non-canonical motifs better and compress homopolymers."""
    elements = repeat_structure.split('-')
    combined = []
    current_motif = ''
    count = 0

    for element in elements:
        # Check if element is a homopolymer (all same base) with length > 1
        if element and len(element) > 1 and len(set(element)) == 1 and element[0] in "ACGT":
            # Convert homopolymer to compressed format (only if more than 1 base)
            element = f"({element[0]}){len(element)}"
        
        if element == current_motif:
            count += 1
        else:
            if count > 0:
                combined.append(f"({current_motif}){count}" if count > 1 else current_motif)
            current_motif = element
            count = 1
    
    if count > 0:
        combined.append(f"({current_motif}){count}" if count > 1 else current_motif)
    
    return '-'.join(combined)

def count_motifs(repeat_structure: str, canonical_motifs: List[str], treat_rotations_as_canonical: bool = False) -> Tuple[int, int, int, int]:
    """Get the counts of canonical and non-canonical motifs, including length comparison.
    
    Counts motifs with special handling for homopolymers:
    - (CGG)13 = 13 canonical copies (expands)
    - (G)2 = 1 non-canonical motif (homopolymer, does NOT expand)
    - G-G = 2 non-canonical motifs (separate elements)
    """
    elements = repeat_structure.split('-')
    canonical_count = 0
    non_canonical_count = 0
    non_canonical_same_length = 0
    non_canonical_different_length = 0
    canonical_length = len(canonical_motifs[0]) if canonical_motifs else 0
    
    # Get all canonical motifs (with rotations if needed)
    if treat_rotations_as_canonical:
        canonical_set = set()
        for canonical_motif in canonical_motifs:
            for i in range(len(canonical_motif)):
                canonical_set.add(canonical_motif[i:] + canonical_motif[:i])
    else:
        canonical_set = set(canonical_motifs)

    for element in elements:
        if '(' in element:
            # Extract motif and count from (MOTIF)N format
            motif = element.strip('()').split(')')[0]
            count = int(element.strip('()').split(')')[1])
        else:
            motif = element
            count = 1
        
        # Check if this is a homopolymer (single base repeated)
        is_homopolymer = len(motif) == 1 and motif in "ACGT"
        
        if motif in canonical_set:
            # Canonical motifs: always expand the count
            canonical_count += count
        else:
            # Non-canonical motifs:
            # - Homopolymers like (G)2: count as 1 motif occurrence
            # - Regular motifs like (CAA)3: expand to 3 occurrences
            if is_homopolymer:
                non_canonical_count += 1
                # Homopolymer length = count (e.g., (G)2 = 2bp)
                if count == canonical_length:
                    non_canonical_same_length += 1
                else:
                    non_canonical_different_length += 1
            else:
                non_canonical_count += count
                if len(motif) == canonical_length:
                    non_canonical_same_length += count
                else:
                    non_canonical_different_length += count

    return canonical_count, non_canonical_count, non_canonical_same_length, non_canonical_different_length

def calculate_non_canonical_base_percentage(repeat_structure, sequence_length, canonical_motifs: List[str]):
    """Calculate non-canonical base proportion in the entire sequence."""
    elements = repeat_structure.split('-')
    non_canonical_length = 0
    canonical_set = set(canonical_motifs)

    for element in elements:
        if element not in canonical_set and '(' not in element:
            non_canonical_length += len(element)
        elif '(' in element:
            motif, count = element.strip('()').split(')')
            if motif not in canonical_set:
                non_canonical_length += len(motif) * int(count)

    return (non_canonical_length / sequence_length) * 100 if sequence_length > 0 else 0.0

def extract_all_non_canonical_motifs(repeat_structure: str,
                                     canonical_motifs: List[str]) -> List[str]:
    elems = repeat_structure.split('-')
    canonical_set = set(canonical_motifs)
    seen = set()
    ordered = []
    for elem in elems:
        if '(' in elem:
            motif, _ = elem.strip('()').split(')')
        else:
            motif = elem
        if motif in canonical_set:
            continue
        if motif not in seen:
            seen.add(motif)
            ordered.append(motif)
    return ordered

def try_parse_eif4a3_chunk(chunk: str, known_patterns: set) -> list:
    """Try to parse a chunk into known EIF4A3 patterns (18-23 mers)."""
    # Try to greedily match known patterns
    result = []
    i = 0
    
    while i < len(chunk):
        matched = False
        # Try patterns from longest to shortest
        for known in sorted(known_patterns, key=len, reverse=True):
            if chunk[i:i+len(known)] == known:
                result.append(known)
                i += len(known)
                matched = True
                break
        
        if not matched:
            # Try generic lengths (20, 21, 22, 23, 19, 18)
            for try_len in [20, 21, 22, 23, 19, 18]:
                if i + try_len <= len(chunk):
                    candidate = chunk[i:i+try_len]
                    result.append(candidate)
                    known_patterns.add(candidate)  # Learn it
                    i += try_len
                    matched = True
                    break
        
        if not matched:
            # Remaining fragment
            if i < len(chunk):
                result.append(chunk[i:])
            break
    
    return result

def process_sequence(seq_id: str, sequence: str, canonical_motifs: List[str], max_mers: int, colors: Dict[str, str], locus: str, treat_rotations_as_canonical: bool = False) -> Tuple[int, int, str, List[str], str]:
    """Process a single sequence and return motif statistics."""
    # Use first canonical motif as the "primary" for display purposes
    primary_canonical = canonical_motifs[0] if canonical_motifs else ""
    
    sequence_length = len(sequence)
    repeat_copy = round(sequence_length / max_mers)


    # Standard extraction for non-EIF4A3 loci
    # Extract regions between EXACT canonical motif occurrences
    # (treat_rotations_as_canonical only affects counting/coloring, not extraction)
    non_canonical_motifs_raw = []
    start_index = 0
    
    # Use ONLY the first canonical motif for extraction (not rotations)
    exact_canonical = canonical_motifs[0]
    
    while start_index < len(sequence):
        # Find next occurrence of EXACT canonical motif
        next_index = sequence.find(exact_canonical, start_index)
        
        if next_index == -1:
            # No more canonical motifs found
            remaining = sequence[start_index:]
            if remaining:
                non_canonical_motifs_raw.append(remaining)
            break
        
        # Extract sequence between current position and next canonical
        between = sequence[start_index:next_index]
        if between:
            non_canonical_motifs_raw.append(between)
        start_index = next_index + len(exact_canonical)

    non_canonical_motifs_str = ','.join(non_canonical_motifs_raw)
    
    # Use colorize_string_tsv to break down motifs for submotif extraction
    # It will handle the max_mers chunking internally
    _, motif_counts = colorize_string_tsv(sequence, canonical_motifs, non_canonical_motifs_str, colors, max_mers, locus)
    submotifs_for_parsing = sorted(motif_counts.keys(), key=len, reverse=True)

    repeat_structure = format_repeats(sequence, canonical_motifs, submotifs_for_parsing, locus)
    post_processed_structure = post_process_repeat_structure(repeat_structure, non_canonical_motifs_str)
    
    # Extract actual non-canonical motifs from the repeat structure
    tsv_submotifs = []
    seen = set()
    canonical_set = set(canonical_motifs)
    
    # If treating rotations as canonical, expand the canonical set to include all rotations
    if treat_rotations_as_canonical:
        expanded_canonical = set()
        for canonical_motif in canonical_motifs:
            for i in range(len(canonical_motif)):
                expanded_canonical.add(canonical_motif[i:] + canonical_motif[:i])
        canonical_set = expanded_canonical
    
    # Parse the repeat structure to extract motifs
    elements = repeat_structure.split('-')
    for element in elements:
        # Handle (MOTIF)N format
        if '(' in element:
            motif = element.split('(')[1].split(')')[0]
        else:
            motif = element
        
        # Only include non-canonical motifs (excluding rotations if applicable)
        if motif and motif not in canonical_set and motif not in seen:
            tsv_submotifs.append(motif)
            seen.add(motif)

    return sequence_length, repeat_copy, primary_canonical, tsv_submotifs, post_processed_structure

def process_sequences(file_path: str, canonical_motifs: List[str], max_mers: int, output_file: str, locus: str, treat_rotations_as_canonical: bool = False) -> str:
    """Process sequences and write results to TSV."""
    sequences = read_fasta(file_path)
    non_can_motifs_str = calculate_non_canonical_motifs(sequences, canonical_motifs, treat_rotations_as_canonical)
    
    # Detect if canonical motifs likely came from degenerate expansion
    # If we have 4+ canonical motifs of the same length, likely from degenerate bases (e.g., GCN -> 4 motifs)
    has_degenerate = len(canonical_motifs) >= 4 and len(set(len(m) for m in canonical_motifs)) == 1
    
    # Use frequency-based color mapping instead of alphabetical
    colors = generate_color_mapping_by_frequency(sequences, canonical_motifs, max_mers, treat_rotations_as_canonical, has_degenerate)
    
    # Canonical motifs display string (for TSV header)
    canonical_display = ",".join(canonical_motifs)
    
    with open(output_file, 'w') as out_file:
        out_file.write("Sample ID\tSequence Length\tRepeat Copy\tCanonical Motif\tSortedNon-CanonicalMotifs\tRepeat Structure\tCanonical Motif Count\tNon-Canonical Motif Count\t%Non-Canonical Base\tNon-Canonical Same Length\tNon-Canonical Different Length\n")
        for seq_id, sequence in sequences.items():
            stats = process_sequence(seq_id, sequence, canonical_motifs, max_mers, colors, locus, treat_rotations_as_canonical)
            submotifs_str = ','.join(stats[3])
            canonical_count, non_canonical_count, non_canonical_same_length, non_canonical_different_length = count_motifs(stats[4], canonical_motifs, treat_rotations_as_canonical)
            non_canonical_base_percentage = calculate_non_canonical_base_percentage(stats[4], stats[0], canonical_motifs)
            # Use canonical_display instead of stats[2] to show all canonical motifs
            out_file.write(f"{seq_id}\t{stats[0]}\t{stats[1]}\t{canonical_display}\t{submotifs_str}\t{stats[4]}\t{canonical_count}\t{non_canonical_count}\t{non_canonical_base_percentage:.2f}\t{non_canonical_same_length}\t{non_canonical_different_length}\n")
    return non_can_motifs_str

def canonical_rotation(motif: str) -> str:
    """Return a canonical representative for all rotations of a motif."""
    if not motif:
        return motif
    doubled = motif + motif
    rotations = [doubled[i:i+len(motif)] for i in range(len(motif))]
    return min(rotations)

def generate_color_mapping(non_can_motifs_str: str, canonical_motifs: List[str], max_mers: int) -> Dict[str, str]:
    """Generate color mapping; rotationally equivalent motifs share a color."""
    motifs_set = set()
    for motif in non_can_motifs_str.split(','):
        if not motif:
            continue
        for i in range(0, len(motif), max_mers):
            sub = motif[i:i+max_mers]
            if sub:
                motifs_set.add(sub)

    # Group submotifs by rotation-canonical representative
    canonical_classes = {}
    for m in motifs_set:
        rep = motif_key(m)
        canonical_classes.setdefault(rep, set()).add(m)

    color_list = list(DISTINCT_COLORS.values())
    color_mapping: Dict[str, str] = {}

    # All canonical motifs get the canonical color, keyed by their rotation-canonical form
    for canonical_motif in canonical_motifs:
        canon_rep = motif_key(canonical_motif)
        color_mapping[canon_rep] = CANONICAL_COLOR

    i = 0
    for rep in sorted(canonical_classes.keys()):
        # Skip if this is any of the canonical motifs
        if any(rep == motif_key(cm) for cm in canonical_motifs):
            continue
        color_mapping[rep] = color_list[i % len(color_list)]
        i += 1

    for rep, col in list(color_mapping.items()):
        base = simplify_motif(rep)   # collapses 'GGGG' -> 'G', 'AAAAAA' -> 'A'
        if base != rep and base not in color_mapping:
            color_mapping[base] = col

    return color_mapping

def generate_color_mapping_by_frequency(sequences: Dict[str, str], canonical_motifs: List[str], max_mers: int, treat_rotations_as_canonical: bool = False, has_degenerate_canonical: bool = False, num_nc_motifs: int = 10) -> Dict[str, str]:
    """Generate color mapping based on motif frequency; most common motifs get first colors.
    
    If has_degenerate_canonical is True, canonical motifs get distinct colors instead of all being CANONICAL_COLOR.
    If multiple canonical motifs are specified (not from degenerate expansion), first gets lightskyblue, others get distinct colors.
    When treat_rotations_as_canonical is False, rotations get different colors.
    Only the top num_nc_motifs non-canonical motifs get unique colors; rest can reuse colors.
    """
    # Count all motif occurrences
    motif_counts: Dict[str, int] = {}
    
    for sequence in sequences.values():
        non_canonical_motifs = count_non_canonical_motifs(sequence, canonical_motifs, treat_rotations_as_canonical)
        for nc_motif in non_canonical_motifs:
            if not nc_motif:
                continue
            # Break into max_mers chunks
            for i in range(0, len(nc_motif), max_mers):
                sub = nc_motif[i:i+max_mers]
                if sub:
                    rep = motif_key_for_coloring(sub, treat_rotations_as_canonical)
                    motif_counts[rep] = motif_counts.get(rep, 0) + 1
    
    # Get canonical reps to exclude (using appropriate key function)
    canonical_reps = set()
    for canonical_motif in canonical_motifs:
        canonical_reps.add(motif_key_for_coloring(canonical_motif, treat_rotations_as_canonical))
    
    # Sort by frequency (most common first)
    sorted_reps = sorted(
        [r for r in motif_counts.keys() if r not in canonical_reps],
        key=lambda r: motif_counts[r],
        reverse=True
    )
    
    # Assign colors by frequency
    color_list = list(DISTINCT_COLORS.values())
    color_mapping: Dict[str, str] = {}
    
    # Handle canonical motif colors
    if has_degenerate_canonical:
        # Each canonical variant from degenerate expansion gets a distinct color from the color palette
        for i, canonical_motif in enumerate(canonical_motifs):
            canon_rep = motif_key_for_coloring(canonical_motif, treat_rotations_as_canonical)
            color_mapping[canon_rep] = color_list[i % len(color_list)]
        # Non-canonical motifs start after canonical colors
        color_offset = len(canonical_motifs)
    elif len(canonical_motifs) > 1:
        # Multiple user-specified canonical motifs: first gets lightskyblue, others get distinct colors
        for i, canonical_motif in enumerate(canonical_motifs):
            canon_rep = motif_key_for_coloring(canonical_motif, treat_rotations_as_canonical)
            if i == 0:
                color_mapping[canon_rep] = CANONICAL_COLOR  # First canonical gets lightskyblue
            else:
                color_mapping[canon_rep] = color_list[(i - 1) % len(color_list)]  # Others get distinct colors
        # Non-canonical motifs start after the canonical colors we've used
        color_offset = len(canonical_motifs) - 1  # -1 because first canonical doesn't use a palette color
    else:
        # Single canonical motif gets the canonical color (traditional behavior)
        for canonical_motif in canonical_motifs:
            canon_rep = motif_key_for_coloring(canonical_motif, treat_rotations_as_canonical)
            color_mapping[canon_rep] = CANONICAL_COLOR
        color_offset = 0
    
    # Only assign unique colors to top displayed motifs to avoid reuse in legend
    # Use ALL available distinct colors before reusing
    available_colors = len(color_list) - color_offset
    
    # Assign unique colors to as many motifs as we have colors for (no reuse)
    for i, rep in enumerate(sorted_reps[:available_colors]):
        color_mapping[rep] = color_list[i + color_offset]
    
    # Remaining motifs beyond our color palette must reuse colors (with modulo)
    for i, rep in enumerate(sorted_reps[available_colors:], start=available_colors):
        color_mapping[rep] = color_list[(i + color_offset) % len(color_list)]
    
    # Add simplified versions
    for rep, col in list(color_mapping.items()):
        base = simplify_motif(rep)
        if base != rep and base not in color_mapping:
            color_mapping[base] = col
    
    return color_mapping

def apply_colors_to_figure(motifs, color_mapping):
    """Apply the color mapping to the motifs in the figure."""
    motif_colors = [color_mapping.get(simplify_motif(motif), "#000000") for motif in motifs]
    return motif_colors

def read_motifs_from_tsv(tsv_file: str) -> Dict[str, List[str]]:
    """Read motifs from TSV file."""
    motifs_dict = {}
    try:
        with open(tsv_file, 'r') as file:
            next(file)  # Skip header
            for line in file:
                parts = line.strip().split('\t')
                motifs_dict[parts[0]] = parts[4].split(',')
    except FileNotFoundError:
        logging.error(f"File not found: {tsv_file}")
        logging.info(f"Current working directory: {os.getcwd()}")
    return motifs_dict

def colorize_string(sequence: str, canonical_motifs: List[str], tsv_motifs: str,
                    colors: Dict[str, str], max_mers: int, locus: str, treat_rotations_as_canonical: bool = False):
    """Colorize a sequence based on motifs, handling rotation coloring based on treat_rotations_as_canonical flag."""
    color_matrix = []
    motif_counts: Dict[str, int] = {}
    start_index = 0
    
    # Build set of canonical motifs for quick lookup
    canonical_set = set(canonical_motifs)
    primary_canonical = canonical_motifs[0] if canonical_motifs else ""
    
    if locus == "EIF4A3" or locus == "PRNP":
        # For EIF4A3 and PRNP, parse using known patterns from TSV (don't chunk into max_mers)
        non_canonical_motifs = sorted(set(filter(None, tsv_motifs.split(','))), key=len, reverse=True)
        known_patterns = set(non_canonical_motifs) if non_canonical_motifs else set()
        
        unmatched_positions = []  # Track where we can't match
        
        while start_index < len(sequence):
            position_start = start_index  # Track starting position
            
            # Check for ANY canonical motif first
            canonical_found = False
            for canonical_motif in canonical_motifs:
                if sequence.startswith(canonical_motif, start_index):
                    # Use the specific canonical motif's color
                    specific_canon_key = motif_key_for_coloring(canonical_motif, treat_rotations_as_canonical)
                    color = colors.get(specific_canon_key, CANONICAL_COLOR)
                    color_matrix.extend([mcolors.to_rgba(color)] * len(canonical_motif))
                    motif_counts[specific_canon_key] = motif_counts.get(specific_canon_key, 0) + 1
                    start_index += len(canonical_motif)
                    canonical_found = True
                    break
            
            if canonical_found:
                continue
            
            # Try to match known non-canonical patterns
            matched = False
            if known_patterns:
                for pattern in sorted(known_patterns, key=len, reverse=True):
                    if sequence.startswith(pattern, start_index):
                        motif_counts[pattern] = motif_counts.get(pattern, 0) + 1
                        color = colors.get(motif_key_for_coloring(pattern, treat_rotations_as_canonical), 'red')
                        color_matrix.extend([mcolors.to_rgba(color)] * len(pattern))
                        start_index += len(pattern)
                        matched = True
                        break
            
            if not matched:
                # Try generic 17-23 base chunks and learn them
                chunk_found = False
                for chunk_len in [20, 21, 22, 23, 19, 18, 17]:
                    if start_index + chunk_len <= len(sequence):
                        chunk = sequence[start_index:start_index+chunk_len]
                        motif_counts[chunk] = motif_counts.get(chunk, 0) + 1
                        known_patterns.add(chunk)
                        color = colors.get(motif_key_for_coloring(chunk, treat_rotations_as_canonical), 'red')
                        color_matrix.extend([mcolors.to_rgba(color)] * len(chunk))
                        start_index += chunk_len
                        chunk_found = True
                        break
                
                if not chunk_found:
                    # Track unmatched positions for debugging if needed
                    remaining = len(sequence) - start_index
                    unmatched_positions.append((start_index, remaining, sequence[start_index:min(start_index+20, len(sequence))]))
                    
                    # Color remaining fragment
                    if remaining > 0:
                        color_matrix.extend([mcolors.to_rgba('lightgray')] * remaining)
                        start_index = len(sequence)
            
    else:
        # Original logic for other loci
        non_canonical_motifs = sorted(set(filter(None, tsv_motifs.split(','))), key=len, reverse=True)

        while start_index < len(sequence):
            # Check for ANY canonical motif
            canonical_found = False
            for canonical_motif in canonical_motifs:
                if sequence.startswith(canonical_motif, start_index):
                    # Use motif_key_for_coloring to respect treat_rotations_as_canonical
                    canon_key = motif_key_for_coloring(canonical_motif, treat_rotations_as_canonical)
                    color = colors.get(canon_key, CANONICAL_COLOR)
                    color_matrix.extend([mcolors.to_rgba(color)] * len(canonical_motif))
                    motif_counts[canon_key] = motif_counts.get(canon_key, 0) + 1
                    start_index += len(canonical_motif)
                    canonical_found = True
                    break
            
            if not canonical_found:
                motif_found = False
                for motif in non_canonical_motifs:
                    if sequence.startswith(motif, start_index):
                        for i in range(0, len(motif), max_mers):
                            sub_motif = motif[i:i+max_mers]
                            if not sub_motif:
                                continue
                            motif_counts[sub_motif] = motif_counts.get(sub_motif, 0) + 1
                            rep = motif_key_for_coloring(sub_motif, treat_rotations_as_canonical)
                            color = colors.get(rep, 'white')
                            color_matrix.extend([mcolors.to_rgba(color)] * len(sub_motif))
                        start_index += len(motif)
                        motif_found = True
                        break
                if not motif_found:
                    # Single unmatched base - color it light gray so it's visible but distinct
                    color_matrix.append(mcolors.to_rgba('lightgray'))
                    start_index += 1

    return color_matrix, motif_counts

def colorize_string_tsv(sequence: str, canonical_motifs: List[str], tsv_motifs: str,
                        colors: Dict[str, str], max_mers: int, locus: str):
    """TSV-only version: keep original behavior so TSV does not change."""
    color_matrix = []
    motif_counts = {}
    start_index = 0
    non_canonical_motifs = sorted(set(tsv_motifs.split(',')), key=len, reverse=True)

    # Use primary canonical for color lookup
    primary_canonical = canonical_motifs[0] if canonical_motifs else ""
    canon_rep = motif_key(primary_canonical)
    canon_color = mcolors.to_rgba(colors.get(canon_rep, CANONICAL_COLOR))

    while start_index < len(sequence):
        # Check for ANY canonical motif
        canonical_found = False
        for canonical_motif in canonical_motifs:
            if sequence.startswith(canonical_motif, start_index):
                color_matrix.extend([canon_color] * len(canonical_motif))
                start_index += len(canonical_motif)
                canonical_found = True
                break
        
        if not canonical_found:
            motif_found = False
            for motif in non_canonical_motifs:
                if sequence.startswith(motif, start_index):
                    for i in range(0, len(motif), max_mers):
                        sub_motif = motif[i:i+max_mers]
                        motif_counts[sub_motif] = motif_counts.get(sub_motif, 0) + 1
                        # Try exact, then simplified, else white
                        color = colors.get(motif_key(sub_motif), 'white')
                        color_matrix.extend([mcolors.to_rgba(color)] * len(sub_motif))
                    start_index += len(motif)
                    motif_found = True
                    break
            if not motif_found:
                color_matrix.append(mcolors.to_rgba('white'))
                start_index += 1

    return color_matrix, motif_counts

def simplify_motif(motif: str) -> str:
    """Simplify repeating motifs to their smallest repeating unit (NO ellipsis)."""
    for k in range(1, len(motif) // 2 + 1):
        if len(motif) % k == 0 and motif == motif[:k] * (len(motif) // k):
            return motif[:k]
    return motif

def label_motif(motif: str, max_len: int = 20) -> str:
    """Short label for legend/UI only."""
    return motif if len(motif) <= max_len else motif[:15] + "..."

def motif_key(m: str) -> str:
    if not m: return m
    if len(set(m)) == 1: return m[0]
    return canonical_rotation(m)

def motif_key_for_coloring(m: str, treat_rotations_as_canonical: bool) -> str:
    """
    Generate motif key for coloring purposes.
    When treat_rotations_as_canonical is False, rotations get different keys (different colors).
    When treat_rotations_as_canonical is True, rotations share the same key (same color).
    """
    if not m: return m
    if len(set(m)) == 1: return m[0]
    if treat_rotations_as_canonical:
        return canonical_rotation(m)
    else:
        # Don't normalize rotations - each rotation gets its own key
        return m

def visualize_motifs(
    sequences: Dict[str, str],
    canonical_motifs: List[str],
    locus: str,
    max_mers: int,
    vlines: List[Tuple[int, str]],
    output_file: str,
    colors: Dict[str, str],
    motifs_dict: Dict[str, List[str]],
    num_nc_motifs: int,
    show_y_label: bool = True,
    sort_sequences: str = "input",
    user_specified_motifs: List[str] = None,
    treat_rotations_as_canonical: bool = False
):
    """Visualize motifs without grouping."""

    # Sequence sorting at the top
    if sort_sequences == 'largest':
        sequences_list = sorted(sequences.items(), key=lambda x: len(x[1]), reverse=True)
    elif sort_sequences == 'smallest':
        sequences_list = sorted(sequences.items(), key=lambda x: len(x[1]))
    else:  # 'input'
        sequences_list = list(sequences.items())

    max_seq_len = max(len(seq) for _, seq in sequences_list)
    max_repeat_count = max_seq_len // max_mers

    middle_repeat = max_repeat_count // 2
    last_repeat = max_repeat_count

    color_matrices = []
    total_motif_counts = {}

    for seq_id, sequence in sequences_list:
        motifs_str = ','.join(motifs_dict.get(seq_id, []))
        color_matrix, seq_motif_counts = colorize_string(sequence, canonical_motifs, motifs_str, colors, max_mers, locus, treat_rotations_as_canonical)
        color_matrices.append(color_matrix)
        for motif, count in seq_motif_counts.items():
            total_motif_counts[motif] = total_motif_counts.get(motif, 0) + count

    plt.rcParams['font.family'] = 'Arial'
    y_labels = [seq_id for seq_id, _ in sequences_list]
    
    max_repeat_count = max(len(seq) // max_mers for _, seq in sequences_list)
    color_matrices = [
        matrix + [(1, 1, 1, 1)] * (max_seq_len - len(matrix))
        for matrix in color_matrices
    ]

    color_stack = np.array(color_matrices)
    n_sequences = len(sequences_list)
    fig_width = 14

    # Dynamic figure height block
    if n_sequences == 1:
        fig_height = 1.4
    elif n_sequences == 2:
        fig_height = 1.8
    elif n_sequences == 3:
        fig_height = 2.2
    elif n_sequences == 4:
        fig_height = 2.6
    elif n_sequences <= 10:
        fig_height = 2.8 + 0.25 * n_sequences
    elif n_sequences <= 20:
        fig_height = 3.0 + 0.25 * n_sequences
    elif n_sequences <= 30:
        fig_height = 3.5 + 0.25 * n_sequences
    elif n_sequences <= 40:
        fig_height = 4 + 0.18 * n_sequences
    else:
        fig_height = min(28, 6.5 + 0.3 * n_sequences)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    ax.imshow(color_stack, aspect="auto", interpolation="nearest", extent=[0, max_seq_len, 0, n_sequences])
    ax.set_xlim(0, max_seq_len)
   
    xticks = [middle_repeat * max_mers, last_repeat * max_mers]
    xtick_labels = [str(middle_repeat), str(last_repeat)]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, fontsize=20, fontfamily='Arial')
    ax.set_xlabel(f'Repeat count', fontsize=22, fontfamily='Arial')

    ax.set_yticks(range(len(y_labels)))
    if show_y_label:
        ax.set_yticklabels(y_labels[::-1], fontsize=19, fontfamily='Arial')
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis='y', left=False, labelleft=False)

    #ax.set_title(f'{locus}', fontsize=20, fontweight='bold', style='italic')
    ax.set_ylim(0, n_sequences)
    ax.invert_yaxis()

    # Optional vertical lines
    if vlines:
        for i, (repeat_idx, color) in enumerate(vlines):
            x_pos = repeat_idx * max_mers
            ax.axvline(x=x_pos, color=color, linestyle='dotted', linewidth=3)
            if i % 2 == 0:
                ax.text(
                    x_pos, 1.0, str(repeat_idx),
                    rotation=90,
                    va='bottom', ha='center',
                    fontsize=18, fontfamily='Arial',
                    transform=ax.get_xaxis_transform()
                )
            else:
                ax.text(
                    x_pos, 0.00, str(repeat_idx),
                    rotation=90,
                    va='top', ha='center',
                    fontsize=18, fontfamily='Arial',
                    transform=ax.get_xaxis_transform()
                )

    # Legend setup
    n = num_nc_motifs
    
    # Collapse counts by motif_key
    collapsed_counts: Dict[str, int] = {}
    for motif, cnt in total_motif_counts.items():
        rep = motif_key(motif)
        collapsed_counts[rep] = collapsed_counts.get(rep, 0) + cnt

    # Get canonical reps (need to use same key function as color mapping)
    canonical_reps = set()
    for canonical_motif in canonical_motifs:
        canonical_reps.add(motif_key_for_coloring(canonical_motif, treat_rotations_as_canonical))
    
    # Detect if using degenerate bases (multiple canonical variants)
    has_degenerate_canon = len(canonical_motifs) >= 4 and len(set(len(m) for m in canonical_motifs)) == 1
    
    # Use first canonical motif as the display label
    primary_canonical = canonical_motifs[0] if canonical_motifs else ""
    
    # Pick a "preferred" label for each rotation class from TSV motifs
    rep_label_freq: Dict[str, Dict[str, int]] = {}
    for seq_id, motif_list in motifs_dict.items():
        for m in motif_list:
            if not m:
                continue
            rep = motif_key(m)
            rep_label_freq.setdefault(rep, {})
            rep_label_freq[rep][m] = rep_label_freq[rep].get(m, 0) + 1

    def preferred_label(rep: str) -> str:
        # If TSV contains an orientation for this rep, use the most common one
        if rep in rep_label_freq and rep_label_freq[rep]:
            return max(rep_label_freq[rep].items(), key=lambda kv: kv[1])[0]
        # Fallback: use the rep itself
        return rep

    # Build legend based on whether degenerate bases were used OR multiple canonical motifs specified
    if has_degenerate_canon or len(canonical_motifs) > 1:
        # Show each canonical variant separately with distinct colors
        legend_count_reps = [motif_key(cm) for cm in canonical_motifs if motif_key(cm) in collapsed_counts]
        legend_reps = [motif_key_for_coloring(cm, treat_rotations_as_canonical) for cm in canonical_motifs if motif_key(cm) in collapsed_counts]
        # Use the ACTUAL canonical motifs from user input, not TSV-derived rotations
        legend_labels = [cm for cm in canonical_motifs if motif_key(cm) in collapsed_counts]
        
        # Add top non-canonical motifs
        if user_specified_motifs:
            user_count_reps = [motif_key(m) for m in user_specified_motifs]
            user_color_reps = [motif_key_for_coloring(m, treat_rotations_as_canonical) for m in user_specified_motifs]
            # Filter based on count reps, but use color reps for legend
            selected_indices = [i for i, r in enumerate(user_count_reps) if r in collapsed_counts and user_color_reps[i] not in canonical_reps][:max(0, n - len(legend_reps))]
            top_nc_reps = [user_color_reps[i] for i in selected_indices]
            top_nc_labels = [preferred_label(user_count_reps[i]) for i in selected_indices]
        else:
            top_count_reps = sorted(
                [r for r in collapsed_counts.keys() if r not in [motif_key(cm) for cm in canonical_motifs]],
                key=lambda r: collapsed_counts[r],
                reverse=True
            )[:max(0, n - len(legend_reps))]
            # Convert count reps to color reps by finding the actual motif
            top_nc_reps = []
            top_nc_labels = []
            for count_rep in top_count_reps:
                # Find first motif with this count_rep
                found = False
                for seq_motifs in motifs_dict.values():
                    for m in seq_motifs:
                        if m and motif_key(m) == count_rep:
                            top_nc_reps.append(motif_key_for_coloring(m, treat_rotations_as_canonical))
                            top_nc_labels.append(preferred_label(count_rep))
                            found = True
                            break
                    if found:
                        break
        
        legend_reps.extend(top_nc_reps)
        legend_labels.extend(top_nc_labels)
        
        # Deduplicate to ensure no duplicate colors in legend
        seen_reps = set()
        unique_legend_reps = []
        unique_legend_labels = []
        for i, rep in enumerate(legend_reps):
            if rep not in seen_reps:
                seen_reps.add(rep)
                unique_legend_reps.append(rep)
                unique_legend_labels.append(legend_labels[i])
        legend_reps = unique_legend_reps
        legend_labels = unique_legend_labels
    else:
        # Traditional behavior: show canonical as one entry
        # IMPORTANT: Use the actual user-specified canonical motif, not any TSV-derived rotation
        if len(canonical_motifs) > 1:
            canonical_display = ','.join(canonical_motifs)  # Show all if multiple (e.g., CAG,CAA)
        else:
            canonical_display = canonical_motifs[0]  # Use FIRST canonical from user input (e.g., CAG)
        
        # Top motifs excluding canonical representatives
        if user_specified_motifs:
            user_count_reps = [motif_key(m) for m in user_specified_motifs]
            user_color_reps = [motif_key_for_coloring(m, treat_rotations_as_canonical) for m in user_specified_motifs]
            # Get indices of motifs that are in collapsed_counts
            selected_indices = [i for i, r in enumerate(user_count_reps) if r in collapsed_counts][:max(0, n - 1)]
            top_reps = [user_color_reps[i] for i in selected_indices]
            top_labels = [preferred_label(user_count_reps[i]) for i in selected_indices]
        else:
            top_count_reps = sorted(
                [r for r in collapsed_counts.keys() if r not in [motif_key(cm) for cm in canonical_motifs]],
                key=lambda r: collapsed_counts[r],
                reverse=True
            )[: max(0, n - 1)]
            # Convert to color reps
            top_reps = []
            top_labels = []
            for count_rep in top_count_reps:
                # Find first actual motif with this count_rep
                found = False
                for seq_motifs in motifs_dict.values():
                    for m in seq_motifs:
                        if m and motif_key(m) == count_rep:
                            top_reps.append(motif_key_for_coloring(m, treat_rotations_as_canonical))
                            top_labels.append(preferred_label(count_rep))
                            found = True
                            break
                    if found:
                        break

        legend_reps = [motif_key_for_coloring(canonical_motifs[0], treat_rotations_as_canonical) if canonical_motifs else None] + top_reps
        legend_labels = [canonical_display] + top_labels


    # Deduplicate legend_reps by key while preserving order
    seen_reps = set()
    unique_legend_reps = []
    unique_legend_labels = []
    for i, rep in enumerate(legend_reps):
        if rep not in seen_reps:
            seen_reps.add(rep)
            unique_legend_reps.append(rep)
            unique_legend_labels.append(legend_labels[i])
    
    legend_labels = unique_legend_labels[:n]
    legend_reps = unique_legend_reps[:n]

    # Deduplicate by actual color values to prevent duplicate colors in legend
    legend_colors = []
    final_legend_labels = []
    seen_colors = set()
    
    for i, rep in enumerate(legend_reps):
        # Check if this rep is in canonical_reps
        if rep in canonical_reps:
            color = colors.get(rep, CANONICAL_COLOR)
            label = legend_labels[i] + '*'
        else:
            color = colors.get(rep, '#000000')
            label = legend_labels[i]
        
        # Only add if this color hasn't been used yet
        if color not in seen_colors:
            seen_colors.add(color)
            legend_colors.append(color)
            final_legend_labels.append(label)

    patches = [
        plt.Line2D([0], [0], marker='o', color='w', label=label,
                   markersize=16, markerfacecolor=color)
        for label, color in zip(final_legend_labels, legend_colors)
    ]
    box = ax.get_position()
    fig.legend(handles=patches, loc='upper center', fontsize=19, ncol=5,
               bbox_to_anchor=(box.x0 + box.width / 2, box.y0 - 0.11), bbox_transform=fig.transFigure, columnspacing=5.0)

    plt.savefig(output_file, dpi=500, bbox_inches='tight')
    plt.close(fig)

def visualize_motifs_stratified(
    sequences: Dict[str, str],
    canonical_motifs: List[str],
    locus: str,
    max_mers: int,
    vlines: List[Tuple[int, str]],
    output_file: str,
    colors: Dict[str, str],
    motifs_dict: Dict[str, List[str]],
    num_nc_motifs: int,
    show_y_label: bool = True,
    sort_sequences: str = "input",
    user_specified_motifs: List[str] = None,
    superpop_mapping: Dict[str, str] = None,
    treat_rotations_as_canonical: bool = False
):
    """Visualize motifs stratified by superpopulation code with shared x-axis."""
    
    # Group sequences by superpopulation code
    populations = {}
    
    for seq_id, sequence in sequences.items():
        # Get superpop from mapping if provided, otherwise auto-detect
        if superpop_mapping and seq_id in superpop_mapping:
            superpop_code = superpop_mapping[seq_id]
        else:
            superpop_code = extract_superpop_from_id(seq_id)
        
        if superpop_code not in populations:
            populations[superpop_code] = []
        populations[superpop_code].append((seq_id, sequence))
    
    # Sort superpops in fixed order if present
    superpop_order = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    sorted_superpops = [
        (sp, populations[sp]) for sp in superpop_order if sp in populations
    ] + [
        (sp, seqs) for sp, seqs in sorted(populations.items()) if sp not in superpop_order
    ]
    
    # Global max length to share x-axis scale
    max_seq_len = max(len(seq) for _, seq in sequences.items())
    max_repeat_count = max_seq_len // max_mers
    middle_repeat = max_repeat_count // 2
    last_repeat = max_repeat_count
    
    # Prepare figure with one row per superpopulation
    num_superpops = len(sorted_superpops)
    fig, axes = plt.subplots(
        num_superpops, 1,
        figsize=(14, 2 * num_superpops + 1.5),  # Add extra space for legend
        sharex=True,
        constrained_layout=False  # Turn off to control spacing manually
    )
    if num_superpops == 1:
        axes = [axes]
    
    # Adjust spacing to make room for legend
    plt.subplots_adjust(bottom=0.15, hspace=0.3)
    
    plt.rcParams['font.family'] = 'Arial'
    
    # Collect motif counts globally for legend
    total_motif_counts: Dict[str, int] = {}
    
    for ax, (superpop_code, seq_list) in zip(axes, sorted_superpops):
        # Apply within-superpop sorting based on sort_sequences flag
        if sort_sequences == 'largest':
            seq_list = sorted(seq_list, key=lambda x: len(x[1]), reverse=True)
        elif sort_sequences == 'smallest':
            seq_list = sorted(seq_list, key=lambda x: len(x[1]))
        # else: keep original order
        
        # Build color matrices for this superpop
        color_matrices = []
        y_labels = []
        for seq_id, sequence in seq_list:
            motifs_str = ','.join(motifs_dict.get(seq_id, []))
            color_matrix, seq_motif_counts = colorize_string(
                sequence, canonical_motifs, motifs_str, colors, max_mers, locus, treat_rotations_as_canonical
            )
            color_matrices.append(color_matrix)
            y_labels.append(seq_id)
            for motif, count in seq_motif_counts.items():
                total_motif_counts[motif] = total_motif_counts.get(motif, 0) + count
        
        # Pad to global max_seq_len
        color_matrices = [
            matrix + [(1, 1, 1, 1)] * (max_seq_len - len(matrix))
            for matrix in color_matrices
        ]
        color_stack = np.array(color_matrices)
        n_sequences = len(seq_list)
        
        ax.imshow(
            color_stack,
            aspect="auto",
            interpolation="nearest",
            extent=[0, max_seq_len, 0, n_sequences]
        )
        ax.set_xlim(0, max_seq_len)
        ax.set_ylim(0, n_sequences)
        ax.invert_yaxis()
        
        # y-axis: labels optional
        ax.set_yticks(range(len(y_labels)))
        if show_y_label:
            ax.set_yticklabels(y_labels[::-1], fontsize=10, fontfamily='Arial')
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis='y', which='both', labelleft=False, left=False)
        
        # Title per superpop
        ax.set_title(
            f'{superpop_code} (n={n_sequences})',
            fontsize=16,
            fontweight='bold',
            fontfamily='Arial'
        )
        
        # Vertical lines in repeat units
        if vlines:
            for repeat_idx, color in vlines:
                x_pos = repeat_idx * max_mers
                if x_pos <= max_seq_len:
                    ax.axvline(x=x_pos, color=color, linestyle='dotted', linewidth=2)
    
    # Shared x-axis ticks (repeat counts)
    xticks = [middle_repeat * max_mers, last_repeat * max_mers]
    xtick_labels = [str(middle_repeat), str(last_repeat)]
    axes[-1].set_xticks(xticks)
    axes[-1].set_xticklabels(xtick_labels, fontsize=14, fontfamily='Arial')
    axes[-1].set_xlabel('Repeat count', fontsize=16, fontfamily='Arial')
    for ax in axes[:-1]:
        ax.tick_params(axis='x', labelbottom=False)
    
    # Build legend using global motif counts
    primary_canonical = canonical_motifs[0]
    collapsed_counts: Dict[str, int] = {}
    for motif, cnt in total_motif_counts.items():
        rep = motif_key(motif)
        collapsed_counts[rep] = collapsed_counts.get(rep, 0) + cnt
    
    canonical_reps = set()
    for cm in canonical_motifs:
        canonical_reps.add(motif_key_for_coloring(cm, treat_rotations_as_canonical))
    
    # Detect if using degenerate bases
    has_degenerate_canon = len(canonical_motifs) >= 4 and len(set(len(m) for m in canonical_motifs)) == 1
    
    # Limit to num_nc_motifs
    n = num_nc_motifs
    
    if has_degenerate_canon or len(canonical_motifs) > 1:
        # Show each canonical variant separately
        legend_count_reps = [motif_key(cm) for cm in canonical_motifs if motif_key(cm) in collapsed_counts]
        legend_reps = [motif_key_for_coloring(cm, treat_rotations_as_canonical) for cm in canonical_motifs if motif_key(cm) in collapsed_counts]
        legend_labels = [cm for cm in canonical_motifs if motif_key(cm) in collapsed_counts]
        
        # Add top non-canonical motifs
        if user_specified_motifs:
            user_count_reps = set(motif_key(m) for m in user_specified_motifs)
            user_specs_list = list(user_specified_motifs)
            top_count_reps = [motif_key(m) for m in user_specs_list if motif_key(m) in collapsed_counts and motif_key_for_coloring(m, treat_rotations_as_canonical) not in canonical_reps][:max(0, n - len(legend_reps))]
            top_nc_reps = [motif_key_for_coloring(m, treat_rotations_as_canonical) for m in user_specs_list if motif_key(m) in top_count_reps]
            top_nc_labels = [m for m in user_specs_list if motif_key(m) in top_count_reps]
        else:
            canonical_count_reps = set(motif_key(cm) for cm in canonical_motifs)
            top_count_reps = sorted(
                [r for r in collapsed_counts.keys() if r not in canonical_count_reps],
                key=lambda r: collapsed_counts[r],
                reverse=True
            )[:max(0, n - len(legend_reps))]
            # Convert to color reps
            top_nc_reps = []
            top_nc_labels = []
            for count_rep in top_count_reps:
                # Find first motif with this count_rep
                found = False
                for seq_motifs in motifs_dict.values():
                    for m in seq_motifs:
                        if m and motif_key(m) == count_rep:
                            top_nc_reps.append(motif_key_for_coloring(m, treat_rotations_as_canonical))
                            top_nc_labels.append(m)
                            found = True
                            break
                    if found:
                        break
        
        legend_reps.extend(top_nc_reps)
        legend_labels.extend(top_nc_labels)
        
        # Deduplicate to ensure no duplicate colors in legend
        seen_reps = set()
        unique_legend_reps = []
        unique_legend_labels = []
        for i, rep in enumerate(legend_reps):
            if rep not in seen_reps:
                seen_reps.add(rep)
                unique_legend_reps.append(rep)
                unique_legend_labels.append(legend_labels[i])
        legend_reps = unique_legend_reps
        legend_labels = unique_legend_labels
    else:
        # Traditional: single canonical entry
        canonical_count_reps = set(motif_key(cm) for cm in canonical_motifs)
        
        if user_specified_motifs:
            user_count_reps = set()
            for m in user_specified_motifs:
                for i in range(len(m)):
                    user_count_reps.add(motif_key(m[i:] + m[:i]))
            
            filtered_counts = {}
            for rep, cnt in collapsed_counts.items():
                if rep in canonical_count_reps or rep in user_count_reps:
                    filtered_counts[rep] = cnt
            collapsed_counts = filtered_counts
        
        canon_count_rep = motif_key(primary_canonical)
        sorted_count_reps = sorted(
            [r for r in collapsed_counts.keys() if r != canon_count_rep],
            key=lambda r: collapsed_counts[r],
            reverse=True
        )
        
        if canon_count_rep in collapsed_counts:
            sorted_count_reps = [canon_count_rep] + sorted_count_reps
        
        legend_labels = []
        legend_reps = []
        
        for count_rep in sorted_count_reps[:n]:
            # For canonical rep, use the user-specified canonical motif
            if count_rep == canon_count_rep:
                if len(canonical_motifs) > 1:
                    label = ','.join(canonical_motifs)  # Show all if multiple
                else:
                    label = canonical_motifs[0]  # Use user-specified (e.g., CAG not AGC)
                # Use the FIRST canonical motif's color key, not arbitrary set element
                color_rep = motif_key_for_coloring(canonical_motifs[0], treat_rotations_as_canonical) if canonical_motifs else count_rep
            else:
                # For non-canonical, use most common rotation from TSV
                label = count_rep
                color_rep = count_rep  # Will be converted below
                for seq_motifs in motifs_dict.values():
                    for m in seq_motifs:
                        if m and motif_key(m) == count_rep:
                            label = m
                            color_rep = motif_key_for_coloring(m, treat_rotations_as_canonical)
                            break
            legend_labels.append(label)
            legend_reps.append(color_rep)
        
        # Deduplicate to ensure no duplicate colors in legend
        seen_reps = set()
        unique_legend_reps = []
        unique_legend_labels = []
        for i, rep in enumerate(legend_reps):
            if rep not in seen_reps:
                seen_reps.add(rep)
                unique_legend_reps.append(rep)
                unique_legend_labels.append(legend_labels[i])
        legend_reps = unique_legend_reps
        legend_labels = unique_legend_labels
    
    legend_colors = []
    final_legend_labels = []
    seen_colors = set()  # Track actual colors used, not just keys
    color_list = list(DISTINCT_COLORS.values())
    next_color_idx = len([v for v in colors.values() if v != CANONICAL_COLOR])
    
    for i, rep in enumerate(legend_reps):
        if rep in canonical_reps:
            color = colors.get(rep, CANONICAL_COLOR)
            label = legend_labels[i] + '*'
        else:
            color = colors.get(rep)
            if color is None:
                # Motif not in original colors dict - assign next available color
                color = color_list[next_color_idx % len(color_list)]
                colors[rep] = color  # Add to colors dict for consistency
                next_color_idx += 1
            label = legend_labels[i]
        
        # Only add if this color hasn't been used yet
        if color not in seen_colors:
            seen_colors.add(color)
            legend_colors.append(color)
            final_legend_labels.append(label)
    
    patches = [
        plt.Line2D([0], [0], marker='o', color='w', label=label,
                   markersize=16, markerfacecolor=color)
        for label, color in zip(final_legend_labels, legend_colors)
    ]
    
    # Place legend at bottom of figure with proper spacing
    fig.legend(handles=patches, loc='lower center', fontsize=19, ncol=5,
               bbox_to_anchor=(0.5, 0.01), columnspacing=5.0)

    plt.savefig(output_file, dpi=500, bbox_inches='tight')
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Process and visualize DNA sequences.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Provide name of the output file")
    parser.add_argument("--canonical-motif", required=True, help="Canonical motif(s). Examples: 'CAG' (single), 'CAG,CAA' (multiple for polyglutamine), 'GCN' (degenerate bases for polyalanine), 'GCN,GCB' (multiple degenerate)")
    parser.add_argument("--locus", required=True, help="Locus name")
    parser.add_argument("--max-mers", type=int, default=None, help="Length of the canonical motif (optional - auto-detected if not provided)")
    parser.add_argument("--vlines", type=str, default=None, help="Thresholds to plot (optional)")
    parser.add_argument("--num-nc-motifs", type=int, default=10, help="Number of top motifs to display on plot (default: 10)")
    parser.add_argument("--non-canonical-motifs", type=str, default=None, help="User-provided non-canonical motifs (comma-separated)")
    parser.add_argument("--hide-y-label", action="store_true", help="Hide the y-axis labels in the plot (default: shown)")
    parser.add_argument("--sort-sequences", choices=['input', 'largest', 'smallest'], default='input', help="Sorting: 'input' (original order), 'largest' (longest at bottom of panel), or 'smallest' (shortest at bottom of panel)")
    parser.add_argument("--treat-motif-rotations", action="store_true", help="Treat rotations of canonical motif as canonical (useful for non-coding repeats). Default: False (rotations are non-canonical, appropriate for coding repeats)")
    parser.add_argument("--html", action="store_true", help="Generate HTML report with embedded visualization and TSV summary")
    parser.add_argument("--stratify-by-superpop", action="store_true", help="Create stratified visualization by superpopulation. Auto-detects from sequence IDs or uses --superpop-file")
    parser.add_argument("--superpop-file", type=str, default=None, help="TSV file mapping sequence IDs to superpopulations (format: seq_id<TAB>superpop). If not provided, auto-detects from FASTA headers")
    args = parser.parse_args()

    output_file_tsv = f"{args.output}.tsv"
    output_file_png = f"{args.output}.png"

    # Parse canonical motifs (handles multiple and degenerate bases)
    canonical_motifs_list, has_degenerate = parse_canonical_motifs(args.canonical_motif)
    logging.info(f"Canonical motifs expanded to: {canonical_motifs_list}")
    if has_degenerate:
        logging.info("Degenerate bases detected - canonical variants will be shown with distinct colors")
    
    # Auto-detect max_mers if not provided
    if args.max_mers is None:
        # Use the length of the first canonical motif
        args.max_mers = len(canonical_motifs_list[0])
        logging.info(f"Auto-detected max_mers: {args.max_mers} (from canonical motif length)")

    non_can_motifs_str = process_sequences(args.input, canonical_motifs_list, args.max_mers, output_file_tsv, args.locus, args.treat_motif_rotations)
    tsv_motifs = read_motifs_from_tsv(output_file_tsv)
    sequences = read_fasta(args.input)

    if args.non_canonical_motifs:
        user_motifs = args.non_canonical_motifs.split(',')
        colors = {motif: list(DISTINCT_COLORS.values())[i % len(DISTINCT_COLORS)] for i, motif in enumerate(user_motifs)}
        # Handle canonical motif colors: first gets lightskyblue, others get distinct colors
        if len(canonical_motifs_list) > 1:
            for i, canon_motif in enumerate(canonical_motifs_list):
                key = motif_key_for_coloring(canon_motif, args.treat_motif_rotations)
                if i == 0:
                    colors[key] = CANONICAL_COLOR
                else:
                    colors[key] = list(DISTINCT_COLORS.values())[(i - 1) % len(DISTINCT_COLORS)]
        else:
            for canon_motif in canonical_motifs_list:
                colors[motif_key_for_coloring(canon_motif, args.treat_motif_rotations)] = CANONICAL_COLOR
    else:
        # Use frequency-based color mapping for visualization
        colors = generate_color_mapping_by_frequency(sequences, canonical_motifs_list, args.max_mers, args.treat_motif_rotations, has_degenerate, args.num_nc_motifs)

    # Ensure all TSV motifs have colors (for all loci)
    # This handles edge cases where TSV has motifs not seen during frequency analysis
    color_list = list(DISTINCT_COLORS.values())
    
    # Find the next available color index (how many distinct color values have been assigned)
    used_colors = set(colors.values()) - {CANONICAL_COLOR}
    color_idx = len(used_colors)
    
    for seq_id, motif_list in tsv_motifs.items():
        for motif in motif_list:
            if motif:
                key = motif_key_for_coloring(motif, args.treat_motif_rotations)
                if key not in colors:
                    # Only use modulo after we've exhausted all 24 colors
                    if color_idx < len(color_list):
                        colors[key] = color_list[color_idx]
                    else:
                        colors[key] = color_list[color_idx % len(color_list)]
                    color_idx += 1

    vlines = ast.literal_eval(args.vlines) if args.vlines else []
    show_y_label = not args.hide_y_label
    
    # Prepare user-specified motifs list (if provided)
    user_motifs_list = args.non_canonical_motifs.split(',') if args.non_canonical_motifs else None
    
    # Load superpopulation mapping if provided
    superpop_mapping = None
    if args.stratify_by_superpop and args.superpop_file:
        superpop_mapping = read_superpop_mapping(args.superpop_file)
        logging.info(f"Loaded superpopulation mapping for {len(superpop_mapping)} sequences")
    
    # Choose visualization function based on stratification flag
    if args.stratify_by_superpop:
        visualize_motifs_stratified(
            sequences, canonical_motifs_list, args.locus, args.max_mers, vlines,
            output_file_png, colors, tsv_motifs, num_nc_motifs=args.num_nc_motifs,
            show_y_label=show_y_label, sort_sequences=args.sort_sequences,
            user_specified_motifs=user_motifs_list, superpop_mapping=superpop_mapping,
            treat_rotations_as_canonical=args.treat_motif_rotations)
    else:
        visualize_motifs(
            sequences, canonical_motifs_list, args.locus, args.max_mers, vlines,
            output_file_png, colors, tsv_motifs, num_nc_motifs=args.num_nc_motifs,
            show_y_label=show_y_label, sort_sequences=args.sort_sequences,
            user_specified_motifs=user_motifs_list,
            treat_rotations_as_canonical=args.treat_motif_rotations)

    # Generate enhanced HTML report (optional)
    if args.html:
        html_path = f"{args.output}.html"
        canonical_display_html = ",".join(canonical_motifs_list)
        
        # Read TSV to get all data
        with open(output_file_tsv, 'r') as f:
            lines = f.readlines()
            num_sequences = len(lines) - 1  # Exclude header
            
            # Parse all sequences for table
            tsv_preview = ""
            for line_idx, line in enumerate(lines):
                cols = line.strip().split('\t')
                if len(cols) > 0:
                    tsv_preview += "<tr>"
                    for col_idx, col in enumerate(cols):
                        tag = "th" if line_idx == 0 else "td"
                        # Add special class for Repeat Structure column (column 5)
                        class_attr = ' class="repeat-structure"' if col_idx == 5 else ''
                        tsv_preview += f"<{tag}{class_attr}>{col}</{tag}>"
                    tsv_preview += "</tr>\n"
        
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset='UTF-8'>
    <title>TRMotifAnnotator Report: {args.locus}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 40px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            font-size: 32px;
            margin-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            font-size: 24px;
            margin-top: 30px;
            margin-bottom: 15px;
            border-bottom: 2px solid #3498db;
            padding-bottom: 5px;
        }}
        .metadata {{
            background-color: #ecf0f1;
            padding: 20px;
            border-radius: 5px;
            margin: 20px 0;
            font-size: 16px;
        }}
        .metadata b {{
            color: #2c3e50;
        }}
        img {{
            max-width: 100%;
            border: 1px solid #bdc3c7;
            border-radius: 4px;
            margin: 20px 0;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
            font-size: 14px;
            font-family: Arial, sans-serif;
        }}
        th, td {{
            border: 1px solid #bdc3c7;
            padding: 10px;
            text-align: left;
            white-space: nowrap;
            font-size: 14px;
        }}
        th {{
            background-color: #3498db;
            color: white;
            font-weight: bold;
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        /* Make Repeat Structure column (6th column) allow wrapping */
        th.repeat-structure, td.repeat-structure {{
            white-space: normal;
            word-break: break-word;
            min-width: 300px;
        }}
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        .summary {{
            font-size: 16px;
            color: #7f8c8d;
            margin: 10px 0;
        }}
        code {{
            background-color: #ecf0f1;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
        }}
        .table-wrapper {{
            max-height: 600px;
            overflow-x: auto;
            overflow-y: auto;
            border: 1px solid #bdc3c7;
            border-radius: 4px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>TRMotifAnnotator Report: {args.locus}</h1>
        
        <div class="metadata">
            <p><b>Input FASTA:</b> <code>{args.input}</code></p>
            <p><b>Canonical Motif(s):</b> <code>{canonical_display_html}</code></p>
            <p><b>Locus:</b> {args.locus}</p>
            <p><b>Max k-mers:</b> {args.max_mers}</p>
            <p><b>Total Sequences:</b> {num_sequences}</p>
            <p><b>Output Files:</b> <code>{os.path.basename(output_file_tsv)}</code>, <code>{os.path.basename(output_file_png)}</code></p>
        </div>
        
        <h2>Visualization</h2>
        <img src="{os.path.basename(output_file_png)}" alt="Motif Visualization">
        
        <h2>Complete TSV Data</h2>
        <p class="summary">Showing all {num_sequences} sequence metadata. Full TSV: <code>{os.path.basename(output_file_tsv)}</code></p>
        <div class="table-wrapper">
            <table>
                {tsv_preview}
            </table>
        </div>
    </div>
</body>
</html>"""
        
        with open(html_path, "w") as f:
            f.write(html_content)
        
        print(f"HTML report written to: {html_path}")

if __name__ == "__main__":
    main()

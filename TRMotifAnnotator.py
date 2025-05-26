import sys
import os
import ast
import argparse
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import logging
logging.basicConfig(level=logging.INFO)
from matplotlib.font_manager import FontProperties
from collections import OrderedDict

# Constants
CANONICAL_COLOR = 'slategray'
RANDOM_SEED = 42

DISTINCT_COLORS = OrderedDict([
    ("Blue",         "#1F77B4"), ("Orange",       "#FF7F0E"), ("Green",        "#2CA02C"),
    ("Red-Orange",   "#D62728"), ("Purple",       "#9467BD"), ("Brown",        "#8C564B"),
    ("Pink",         "#E377C2"), ("Yellow-Green", "#BCBD22"), ("Cyan",         "#17BECF"),
    ("Light Pink",   "#F5A3B8"), ("Pale Green",   "#C7E9A9"), ("Bright Yellow","#F2F200"),
    ("Deep Pink",    "#FF1493"), ("Blue Violet",  "#8A2BE2"), ("Tomato Red",   "#FF6347"),
    ("Goldenrod",    "#DAA520"), ("Amber",        "#F4A300"), ("Dark Turquoise","#00CED1"),
    ("Slate Blue",   "#6A5ACD"), ("Dark Red",     "#8B0000"), ("Medium Orchid","#BA55D3"),
    ("Forest Green", "#347C4C"), ("Light Gray",   "#CFCFCF"), ("Bright Red",   "#D43F00"),
    ("Lavender",     "#CE93D8")
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

def count_non_canonical_motifs(sequence: str, canonical_motif: str) -> List[str]:
    """Determine counts of non-canonical motifs."""
    non_canonical_motifs = []
    start_index = 0
    while start_index < len(sequence):
        next_index = sequence.find(canonical_motif, start_index)
        if next_index == -1:
            non_canonical_motifs.append(sequence[start_index:])
            break
        non_canonical_motifs.append(sequence[start_index:next_index])
        start_index = next_index + len(canonical_motif)
    return [motif for motif in non_canonical_motifs if motif]

def calculate_non_canonical_motifs(sequences: Dict[str, str], canonical_motif: str) -> str:
    """Calculate non-canonical motifs for all sequences."""
    non_can_motifs = []
    for sequence in sequences.values():
        non_can_motifs.extend(count_non_canonical_motifs(sequence, canonical_motif))
    return ','.join(set(non_can_motifs))

def format_repeats(sequence: str, canonical_motif: str, submotifs: List[str]) -> str:
    """Format the sequence to show repeated motifs, including submotifs."""
    formatted = []
    i = 0
    length = len(sequence)
    
    while i < length:
        # Check for canonical motif first
        if sequence[i:i+len(canonical_motif)] == canonical_motif:
            count = 1
            j = i + len(canonical_motif)
            while j + len(canonical_motif) <= length and sequence[j:j+len(canonical_motif)] == canonical_motif:
                count += 1
                j += len(canonical_motif)
            formatted.append(f"({canonical_motif}){count}" if count > 1 else canonical_motif)
            i = j
        else:
            # Check for submotifs
            submotif_found = False
            for submotif in sorted(submotifs, key=len, reverse=True):
                if sequence.startswith(submotif, i):
                    formatted.append(submotif)
                    i += len(submotif)
                    submotif_found = True
                    break
            
            if not submotif_found:
                formatted.append(sequence[i])
                i += 1
    
    return '-'.join(formatted)

def post_process_repeat_structure(repeat_structure: str, non_canonical_motifs: str) -> str:
    """Post process repeat structures to consolidate non-canonical motifs better."""
    elements = repeat_structure.split('-')
    combined = []
    current_motif = ''
    count = 0

    for element in elements:
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

def count_motifs(repeat_structure: str, canonical_motif: str) -> Tuple[int, int, int, int]:
    """Get the counts of canonical and non-canonical motifs, including length comparison."""
    elements = repeat_structure.split('-')
    canonical_count = 0
    non_canonical_count = 0
    non_canonical_same_length = 0
    non_canonical_different_length = 0
    canonical_length = len(canonical_motif)

    for element in elements:
        if '(' in element:
            motif, count = element.strip('()').split(')')
            count = int(count)
        else:
            motif = element
            count = 1

        if motif == canonical_motif:
            canonical_count += count
        else:
            non_canonical_count += count
            if len(motif) == canonical_length:
                non_canonical_same_length += count
            else:
                non_canonical_different_length += count

    return canonical_count, non_canonical_count, non_canonical_same_length, non_canonical_different_length

def calculate_non_canonical_base_percentage(repeat_structure, sequence_length, canonical_motif):
    """Calculate non-canonical base proportion in the entire sequence."""
    elements = repeat_structure.split('-')
    non_canonical_length = 0

    for element in elements:
        if element != canonical_motif and '(' not in element:
            non_canonical_length += len(element)
        elif '(' in element:
            motif, count = element.strip('()').split(')')
            if motif != canonical_motif:
                non_canonical_length += len(motif) * int(count)

    return (non_canonical_length / sequence_length) * 100 if sequence_length > 0 else 0.0

def process_sequence(sequence: str, canonical_motif: str, max_mers: int, colors: Dict[str, str]) -> Tuple[int, int, float, str, float, int, int]:
    """Process a single sequence and return motif statistics."""
    non_canonical_motifs = count_non_canonical_motifs(sequence, canonical_motif)
    sequence_length = len(sequence)
    repeat_copy = sequence_length // max_mers
    non_canonical_motifs_str = ','.join(non_canonical_motifs)
    _, motif_counts = colorize_string(sequence, canonical_motif, non_canonical_motifs_str, colors, max_mers)
    submotifs = list(motif_counts.keys())
    repeat_structure = format_repeats(sequence, canonical_motif, submotifs)
    post_processed_structure = post_process_repeat_structure(repeat_structure, non_canonical_motifs_str)
    return sequence_length, repeat_copy, canonical_motif, submotifs, post_processed_structure

def process_sequences(file_path: str, canonical_motif: str, max_mers: int, output_file: str) -> str:
    """Process sequences and write results to TSV."""
    sequences = read_fasta(file_path)
    non_can_motifs_str = calculate_non_canonical_motifs(sequences, canonical_motif)
    colors = generate_color_mapping(non_can_motifs_str, canonical_motif, max_mers)
    with open(output_file, 'w') as out_file:
        out_file.write("Sample ID\tSequence Length\tRepeat Copy\tCanonical Motif\tSortedNon-CanonicalMotifs\tRepeat Structure\tCanonical Motif Count\tNon-Canonical Motif Count\t%Non-Canonical Base\tNon-Canonical Same Length\tNon-Canonical Different Length\n")
        for seq_id, sequence in sequences.items():
            stats = process_sequence(sequence, canonical_motif, max_mers, colors)
            submotifs_str = ','.join(stats[3])
            canonical_count, non_canonical_count, non_canonical_same_length, non_canonical_different_length = count_motifs(stats[4], canonical_motif)
            non_canonical_base_percentage = calculate_non_canonical_base_percentage(stats[4], stats[0], canonical_motif)
            out_file.write(f"{seq_id}\t{stats[0]}\t{stats[1]}\t{stats[2]}\t{submotifs_str}\t{stats[4]}\t{canonical_count}\t{non_canonical_count}\t{non_canonical_base_percentage:.2f}\t{non_canonical_same_length}\t{non_canonical_different_length}\n")
    return non_can_motifs_str

def generate_color_mapping(non_can_motifs_str: str, canonical_motif: str, max_mers: int) -> Dict[str, str]:
    """Generate color mapping."""
    motifs_set = set()
    for motif in non_can_motifs_str.split(','):
        for i in range(0, len(motif), max_mers):
            motifs_set.add(motif[i:i+max_mers])
    motifs_set.add(canonical_motif)

    color_list = list(DISTINCT_COLORS.values())
    
    color_mapping = {canonical_motif: CANONICAL_COLOR}
    for i, motif in enumerate(sorted(motifs_set - {canonical_motif})):
        color_mapping[motif] = color_list[i % len(color_list)]
    
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

def colorize_string(sequence: str, canonical_motif: str, tsv_motifs: str, colors: Dict[str, str], max_mers: int) -> List[Tuple[float, float, float, float]]:
    """Colorize a sequence based on its motifs, prioritizing longer motifs."""
    color_matrix = []
    motif_counts = {}
    start_index = 0
    non_canonical_motifs = sorted(set(tsv_motifs.split(',')), key=len, reverse=True)
    
    while start_index < len(sequence):
        if sequence.startswith(canonical_motif, start_index):
            color_matrix.extend([mcolors.to_rgba(colors[canonical_motif])] * len(canonical_motif))
            start_index += len(canonical_motif)
        else:
            motif_found = False
            for motif in non_canonical_motifs:
                if sequence.startswith(motif, start_index):
                    # Break down longer motifs into max_mers-sized chunks
                    for i in range(0, len(motif), max_mers):
                        sub_motif = motif[i:i+max_mers]
                        simplified_sub_motif = simplify_motif(sub_motif)
                        motif_counts[simplified_sub_motif] = motif_counts.get(simplified_sub_motif, 0) + 1
                        color = colors.get(sub_motif) or colors.get(simplify_motif(sub_motif), 'white')
                        color_matrix.extend([mcolors.to_rgba(color)] * len(sub_motif))
                    start_index += len(motif)
                    motif_found = True
                    break
            if not motif_found:
                color_matrix.append(mcolors.to_rgba('white'))
                start_index += 1
    
    return color_matrix, motif_counts

def simplify_motif(motif: str) -> str:
    """Simplify repeating motifs to their smallest repeating unit."""
    for k in range(1, len(motif) // 2 + 1):
        if motif == motif[:k] * (len(motif) // k):
            return motif[:k]
    return motif[:15] + "..." if len(motif) > 20 else motif

def visualize_motifs(
    sequences: Dict[str, str],
    canonical_motif: str,
    locus: str,
    max_mers: int,
    vlines: List[Tuple[int, str]],
    output_file: str,
    colors: Dict[str, str],
    motifs_dict: Dict[str, List[str]],
    num_nc_motifs: int
):
    """Visualize motifs without grouping by superpopulation or population."""
    
    sorted_sequences = sorted(sequences.items(), key=lambda x: len(x[1]), reverse=True)
    max_seq_len = max(len(seq) for _, seq in sorted_sequences)

    color_matrices = []
    total_motif_counts = {}

    for seq_id, sequence in sorted_sequences:
        motifs_str = ','.join(motifs_dict.get(seq_id, []))
        color_matrix, seq_motif_counts = colorize_string(sequence, canonical_motif, motifs_str, colors, max_mers)
        color_matrices.append(color_matrix)

        for motif, count in seq_motif_counts.items():
            total_motif_counts[motif] = total_motif_counts.get(motif, 0) + count

    plt.rcParams['font.family'] = 'Arial'
    y_labels = [seq_id for seq_id, _ in sorted_sequences]
    max_color_matrix_length = max(len(matrix) for matrix in color_matrices)
    color_matrices = [
        matrix + [(1, 1, 1, 1)] * (max_color_matrix_length - len(matrix))
        for matrix in color_matrices
    ]
    color_stack = np.array(color_matrices)
    total_height = len(sorted_sequences)

    # Plot
    n_sequences = len(color_stack)
    fig_width = 10
    height_per_seq = 0.3
    fig_height = max(6, min(height_per_seq * n_sequences, 30))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    ax.imshow(color_stack, aspect="auto", interpolation="nearest", extent=[0, max_seq_len, 0, n_sequences])

    # Set x-axis and title
    ax.set_xlim(0, max_seq_len)
    ax.set_xticks([0, max_seq_len // 2, max_seq_len])
    ax.set_xticklabels(['0', str(max_seq_len // 2), str(max_seq_len)], fontsize=14, fontfamily='Arial')
    ax.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
    ax.set_yticks(range(len(y_labels)))
    ax.set_yticklabels(y_labels[::-1], fontsize=14, fontfamily='Arial')
    ax.set_title(f'{locus}', fontsize=17, fontweight='bold', style='italic')
    ax.invert_yaxis()

    if vlines:
        for x_pos, color in vlines:
            ax.axvline(x=x_pos, color=color, linestyle='dotted', linewidth=2)
            ax.text(x_pos, total_height + 0.5, str(x_pos), rotation=90, va='top', ha='center', fontsize=12, fontfamily='Arial')

    # Create a common legend
    n = num_nc_motifs
    top_motifs = sorted(total_motif_counts, key=total_motif_counts.get, reverse=True)[:n-1]
    legend_labels = [canonical_motif] + top_motifs
    legend_labels = list(dict.fromkeys(legend_labels))[:n]
    legend_colors = [colors.get(label, '#000000') for label in legend_labels]

    patches = [
        plt.Line2D([0], [0], marker='o', color='w', label=label,
                   markersize=10, markerfacecolor=color)
        for label, color in zip(legend_labels, legend_colors)
    ]
    fig.legend(handles=patches, loc='upper center', fontsize=14, ncol=7,
               bbox_to_anchor=(0.5, -0.01), bbox_transform=fig.transFigure)

    plt.savefig(output_file, dpi=500, bbox_inches='tight')
    plt.close(fig)

def get_motif(sequence, position, motifs):
    for motif in motifs:
        if sequence.startswith(motif, position):
            return motif
    return "None"

def main():
    parser = argparse.ArgumentParser(description="Process and visualize DNA sequences.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Provide name of the output file")
    parser.add_argument("--canonical-motif", required=True, help="Canonical motif")
    parser.add_argument("--locus", required=True, help="Locus name")
    parser.add_argument("--max-mers", type=int, required=True, help="Length of the canonical motif")
    parser.add_argument("--vlines", type=str, default=None, help="Thresholds to plot (optional)")
    parser.add_argument("--num-nc-motifs", type=int, default=10, help="Number of top motifs to display on plot (default: 10)")
    parser.add_argument("--non-canonical-motifs", type=str, default=None, help="User-provided non-canonical motifs (comma-separated)")
    args = parser.parse_args()

    output_file_tsv = f"{args.output}.tsv"
    output_file_png = f"{args.output}.png"
    output_file_html = f"{args.output}.html"

    non_can_motifs_str = process_sequences(args.input, args.canonical_motif, args.max_mers, output_file_tsv)
    tsv_motifs = read_motifs_from_tsv(output_file_tsv)
    sequences = read_fasta(args.input)

    if args.non_canonical_motifs:
        user_motifs = args.non_canonical_motifs.split(',')
        colors = {motif: DISTINCT_COLORS[i % len(DISTINCT_COLORS)] for i, motif in enumerate(user_motifs)}
        colors[args.canonical_motif] = CANONICAL_COLOR
    else:
        colors = generate_color_mapping(non_can_motifs_str, args.canonical_motif, args.max_mers)

    vlines = ast.literal_eval(args.vlines) if args.vlines else []
    
    visualize_motifs(sequences, args.canonical_motif, args.locus, args.max_mers, vlines, 
                                     output_file_png, colors, tsv_motifs, num_nc_motifs=args.num_nc_motifs)

if __name__ == "__main__":
    main()

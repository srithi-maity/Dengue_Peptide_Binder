#!/usr/bin/env python3
"""
Find HIGHLY CONSERVED regions (≥95% or ≥90%) across ALL sequences
"""

from Bio import SeqIO
from collections import Counter
import numpy as np


def find_highly_conserved_regions(sequences_file, min_conservation=95, min_length=10):
    """
    Find regions with high conservation across sequences.
    """
    print(f" Loading ALL sequences...")

    # Load ALL sequences
    records = list(SeqIO.parse(sequences_file, "fasta"))
    total_sequences = len(records)
    print(f" Loaded ALL {total_sequences} sequences")

    # Get minimum sequence length
    min_len = min(len(rec.seq) for rec in records)
    print(f" Using minimum length: {min_len} aa")

    print(f"\n Scanning {min_len} positions for ≥{min_conservation}% conservation...")

    # First, get conservation at each position
    conservation_data = []

    for pos in range(min_len):
        # Get residues from ALL sequences
        residues = []
        for record in records:
            if pos < len(record.seq):
                residue = str(record.seq[pos]).upper()
                residues.append(residue)

        if residues:
            total = len(residues)
            residue_counts = Counter(residues)

            # Get most common residue and its percentage
            most_common_residue, most_common_count = residue_counts.most_common(1)[0]
            conservation = (most_common_count / total) * 100

            conservation_data.append({
                'position': pos + 1,
                'residue': most_common_residue,
                'conservation': conservation,
                'total': total,
                'most_common_count': most_common_count
            })
        else:
            conservation_data.append({
                'position': pos + 1,
                'residue': 'X',
                'conservation': 0,
                'total': 0,
                'most_common_count': 0
            })

    # Find regions with high conservation
    conserved_regions = []
    current_region = []
    current_conservation_values = []

    for pos_data in conservation_data:
        if pos_data['conservation'] >= min_conservation:
            current_region.append(pos_data)
            current_conservation_values.append(pos_data['conservation'])
        else:
            # End of current region
            if len(current_region) >= min_length:
                avg_conservation = np.mean(current_conservation_values)

                # Get the conserved sequence
                conserved_seq = ''.join([pos['residue'] for pos in current_region])

                conserved_regions.append({
                    'start': current_region[0]['position'],
                    'end': current_region[-1]['position'],
                    'length': len(current_region),
                    'sequence': conserved_seq,
                    'avg_conservation': avg_conservation,
                    'min_conservation': min(current_conservation_values),
                    'max_conservation': max(current_conservation_values)
                })

            current_region = []
            current_conservation_values = []

    # Check last region
    if len(current_region) >= min_length:
        avg_conservation = np.mean(current_conservation_values)
        conserved_seq = ''.join([pos['residue'] for pos in current_region])

        conserved_regions.append({
            'start': current_region[0]['position'],
            'end': current_region[-1]['position'],
            'length': len(current_region),
            'sequence': conserved_seq,
            'avg_conservation': avg_conservation,
            'min_conservation': min(current_conservation_values),
            'max_conservation': max(current_conservation_values)
        })

    return conserved_regions, conservation_data


def main():
    input_file = "DENV2_E_proteins_verified.fasta"

    # Try different conservation thresholds
    thresholds = [100, 99, 98, 95, 90, 85, 80]

    print("=" * 80)
    print("FINDING HIGHLY CONSERVED REGIONS IN DENV2 E PROTEIN")
    print("=" * 80)

    best_results = None
    best_threshold = None

    for threshold in thresholds:
        print(f"\n{'=' * 60}")
        print(f"TESTING: ≥{threshold}% CONSERVATION")
        print(f"{'=' * 60}")

        conserved_regions, conservation_data = find_highly_conserved_regions(
            input_file,
            min_conservation=threshold,
            min_length=10
        )

        if conserved_regions:
            print(f"✅ Found {len(conserved_regions)} region(s) at ≥{threshold}% conservation")

            # Show the longest region
            longest = max(conserved_regions, key=lambda x: x['length'])
            print(f"   Longest region: {longest['length']} aa (positions {longest['start']}-{longest['end']})")
            print(f"   Sequence: {longest['sequence']}")

            if best_results is None or len(conserved_regions) > len(best_results):
                best_results = conserved_regions
                best_threshold = threshold
        else:
            print(f"❌ No regions found at ≥{threshold}% conservation")

    # Save best results
    if best_results:
        print(f"\n{'=' * 80}")
        print(f"🎯 BEST RESULTS: ≥{best_threshold}% CONSERVATION")
        print(f"{'=' * 80}")

        print(f"\nFound {len(best_results)} conserved regions:")

        # Save to file
        with open(f"conserved_regions_{best_threshold}percent.txt", "w") as f:
            f.write(f"DENV2 E PROTEIN - CONSERVED REGIONS (≥{best_threshold}% conservation)\n")
            f.write("=" * 80 + "\n\n")

            for i, region in enumerate(best_results, 1):
                f.write(f"Region {i}:\n")
                f.write(f"  Positions: {region['start']}-{region['end']}\n")
                f.write(f"  Length: {region['length']} amino acids\n")
                f.write(f"  Average conservation: {region['avg_conservation']:.1f}%\n")
                f.write(f"  Conservation range: {region['min_conservation']:.1f}%-{region['max_conservation']:.1f}%\n")
                f.write(f"  Sequence: {region['sequence']}\n\n")

                print(f"\n🔬 Region {i}:")
                print(f"   Positions: {region['start']}-{region['end']}")
                print(f"   Length: {region['length']} aa")
                print(f"   Conservation: {region['avg_conservation']:.1f}%")
                print(f"   Sequence: {region['sequence']}")

        print(f"\n💾 Saved to: conserved_regions_{best_threshold}percent.txt")

        # Also save as FASTA
        with open(f"conserved_regions_{best_threshold}percent.fasta", "w") as f:
            for i, region in enumerate(best_results, 1):
                f.write(f">Region_{i}_Pos{region['start']}-{region['end']}_Cons{region['avg_conservation']:.1f}%\n")
                f.write(f"{region['sequence']}\n")

        print(f"💾 FASTA format: conserved_regions_{best_threshold}percent.fasta")

    print(f"\n{'=' * 80}")
    print("ANALYSIS COMPLETE!")
    print("=" * 80)


if __name__ == "__main__":
    main()
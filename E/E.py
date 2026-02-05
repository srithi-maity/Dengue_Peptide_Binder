import os
from Bio import SeqIO
import pandas as pd


# ==========================
# SIMPLE CONSERVATION CHECKER
# ==========================

def load_protein_sequences(fasta_file):
    """Load protein sequences from FASTA"""
    sequences = []
    seq_ids = []

    print(f"Loading sequences from: {fasta_file}")

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_str = str(record.seq)
        sequences.append(seq_str)
        seq_ids.append(record.id)

    print(f"Loaded {len(sequences)} sequences")

    # Check sequence lengths
    lengths = [len(seq) for seq in sequences]
    print(f"Sequence lengths: {min(lengths)} to {max(lengths)} amino acids")

    return sequences, seq_ids


def find_conserved_regions_simple(sequences, min_length=10):
    """
    Find regions where ALL sequences have EXACTLY the same amino acids
    No gaps, no mismatches allowed
    """
    if not sequences:
        return []

    print(f"\nSearching for conserved regions (≥{min_length} AA)...")

    # Find the shortest sequence length
    min_len = min(len(seq) for seq in sequences)
    print(f"Comparing first {min_len} positions (exists in all sequences)")

    conserved_regions = []

    i = 0
    while i < min_len:
        # Get the amino acid at position i from first sequence
        first_aa = sequences[0][i]

        # Check if ALL sequences have the SAME amino acid at position i
        all_same = True
        for seq in sequences[1:]:
            if seq[i] != first_aa:
                all_same = False
                break

        if all_same:
            # Start of potential conserved region
            start = i
            j = i + 1

            # Find how long this conservation lasts
            while j < min_len:
                next_aa = sequences[0][j]

                # Check if all sequences have same AA at position j
                column_same = True
                for seq in sequences[1:]:
                    if seq[j] != next_aa:
                        column_same = False
                        break

                if not column_same:
                    break
                j += 1

            region_length = j - i
            if region_length >= min_length:
                # Get the conserved sequence
                conserved_seq = sequences[0][i:j]

                conserved_regions.append({
                    'start': i,
                    'end': j - 1,
                    'length': region_length,
                    'sequence': conserved_seq
                })
                print(f"  Found: Positions {i + 1}-{j} ({region_length} AA)")

            i = j  # Skip to end of region
        else:
            i += 1

    print(f"\nTotal conserved regions found: {len(conserved_regions)}")
    return conserved_regions


def verify_conservation(sequences, region):
    """
    Double-check that a region is truly conserved
    """
    start = region['start']
    end = region['end']
    conserved_seq = region['sequence']

    print(f"\nVerifying region at positions {start + 1}-{end + 1}:")
    print(f"Consensus: {conserved_seq}")

    problems = []

    for idx, seq in enumerate(sequences):
        region_seq = seq[start:end + 1]

        if region_seq != conserved_seq:
            # Find the mismatch
            for pos in range(len(conserved_seq)):
                if region_seq[pos] != conserved_seq[pos]:
                    problems.append({
                        'sequence': idx + 1,
                        'position': start + pos + 1,
                        'expected': conserved_seq[pos],
                        'found': region_seq[pos]
                    })
                    break

    if not problems:
        print(f"✓ VERIFIED: All {len(sequences)} sequences match exactly")
        return True
    else:
        print(f"✗ PROBLEMS: {len(problems)} sequences don't match")
        for prob in problems[:5]:  # Show first 5 problems
            print(f"  Sequence {prob['sequence']}: Position {prob['position']}: "
                  f"{prob['expected']} → {prob['found']}")
        return False


def analyze_conservation_patterns(sequences, conserved_regions):
    """
    Analyze the conservation patterns
    """
    if not conserved_regions:
        print("\nNo conserved regions found to analyze")
        return

    print(f"\n{'=' * 60}")
    print("CONSERVATION ANALYSIS")
    print(f"{'=' * 60}")

    total_conserved_aa = sum(region['length'] for region in conserved_regions)
    total_positions = min(len(seq) for seq in sequences)
    conservation_percentage = (total_conserved_aa / total_positions) * 100

    print(f"Total positions analyzed: {total_positions} AA")
    print(f"Total conserved positions: {total_conserved_aa} AA")
    print(f"Overall conservation: {conservation_percentage:.1f}%")

    # Show the conserved regions
    print(f"\nConserved regions found:")
    for idx, region in enumerate(conserved_regions, 1):
        start = region['start'] + 1  # 1-based
        end = region['end'] + 1  # 1-based
        seq = region['sequence']

        print(f"\nRegion {idx}: Positions {start}-{end} ({region['length']} AA)")
        print(f"Sequence: {seq}")

        # Show context (10 AA before and after if available)
        if region['start'] >= 10:
            context_start = region['start'] - 10
        else:
            context_start = 0

        if region['end'] + 10 < len(sequences[0]):
            context_end = region['end'] + 10
        else:
            context_end = len(sequences[0]) - 1

        context_seq = sequences[0][context_start:context_end + 1]

        # Mark the conserved region
        conserved_start_in_context = region['start'] - context_start
        conserved_end_in_context = region['end'] - context_start

        print(f"Context: ", end="")
        for pos, aa in enumerate(context_seq):
            pos_in_context = pos
            if conserved_start_in_context <= pos_in_context <= conserved_end_in_context:
                print(f"[{aa}]", end="")
            else:
                print(f" {aa} ", end="")
        print()


def save_results(sequences, conserved_regions, seq_ids, output_dir, serotype):
    """
    Save the conserved regions to files
    """
    os.makedirs(output_dir, exist_ok=True)

    if not conserved_regions:
        print(f"\nNo conserved regions to save")
        return

    # Save summary CSV
    summary_data = []
    for idx, region in enumerate(conserved_regions, 1):
        summary_data.append({
            'region_id': idx,
            'start': region['start'] + 1,
            'end': region['end'] + 1,
            'length': region['length'],
            'sequence': region['sequence']
        })

    summary_file = os.path.join(output_dir, f"{serotype}_conserved_regions_summary.csv")
    df = pd.DataFrame(summary_data)
    df.to_csv(summary_file, index=False)
    print(f"\nSummary saved to: {summary_file}")

    # Save each region as FASTA
    for idx, region in enumerate(conserved_regions, 1):
        region_file = os.path.join(output_dir, f"region_{idx:03d}.fasta")

        with open(region_file, 'w') as f:
            # Write the conserved sequence
            f.write(f">Consensus_{serotype}_region_{idx}_pos{region['start'] + 1}-{region['end'] + 1}\n")
            f.write(f"{region['sequence']}\n\n")

            # Write all sequences for this region
            for seq_id, seq in zip(seq_ids, sequences):
                region_seq = seq[region['start']:region['end'] + 1]
                f.write(f">{seq_id}_region_{idx}\n")
                f.write(f"{region_seq}\n")

        print(f"  Region {idx} saved to: {region_file}")


# ==========================
# MAIN FUNCTION
# ==========================

def analyze_e_protein_conservation(fasta_file, min_length=10):
    """
    Main function to analyze E protein conservation
    """
    print(f"{'=' * 60}")
    print(f"E PROTEIN CONSERVATION ANALYSIS")
    print(f"{'=' * 60}")

    # Step 1: Load sequences
    sequences, seq_ids = load_protein_sequences(fasta_file)

    if len(sequences) < 2:
        print("Need at least 2 sequences for comparison")
        return

    # Step 2: Find conserved regions
    conserved_regions = find_conserved_regions_simple(sequences, min_length)

    if not conserved_regions:
        print(f"\nNO conserved regions found of length {min_length}+ amino acids")
        print("This means your E protein sequences have variations at every position.")
        return

    # Step 3: Verify the conservation
    print(f"\n{'=' * 60}")
    print("VERIFICATION")
    print(f"{'=' * 60}")

    all_verified = True
    for region in conserved_regions:
        if not verify_conservation(sequences, region):
            all_verified = False

    if not all_verified:
        print("\nWARNING: Some regions failed verification!")
        print("The 'conserved' regions might not be truly conserved.")
        return

    # Step 4: Analyze patterns
    analyze_conservation_patterns(sequences, conserved_regions)

    # Step 5: Save results
    output_dir = "e_protein_conserved_regions"
    serotype = os.path.basename(fasta_file).replace('.fasta', '').replace('.fa', '')
    save_results(sequences, conserved_regions, seq_ids, output_dir, serotype)

    print(f"\n{'=' * 60}")
    print(f"ANALYSIS COMPLETE")
    print(f"{'=' * 60}")


# ==========================
# RUN THE ANALYSIS
# ==========================

if __name__ == "__main__":
    # Ask user for input file
    fasta_file = input("Enter the path to your E protein FASTA file: ").strip()

    if not os.path.exists(fasta_file):
        print(f"File not found: {fasta_file}")
        print("\nPlease provide a valid FASTA file with E protein sequences.")
        exit(1)

    # Ask for minimum length
    min_length = input("Minimum conserved region length (default: 10): ").strip()
    if not min_length:
        min_length = 10
    else:
        min_length = int(min_length)

    # Run analysis
    analyze_e_protein_conservation(fasta_file, min_length)
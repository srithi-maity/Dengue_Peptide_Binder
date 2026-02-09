
import os
from Bio import SeqIO

# MOTIF = "VDRGWGNGCGLFGKG"
MOTIF = "VVVLGSQEGAMHTALTGATEIQ"
OUTPUT_FILE = "motif2_percentages.txt"


def check_serotype(serotype_num, folder_name=None):
    """
    Check a specific DENV serotype.
    Returns (serotype_name, percentage, total_sequences, sequences_with_motif)
    """
    # Try different possible folder structures
    possible_folders = [
        folder_name if folder_name else f"denv{serotype_num}",
        f"DENV{serotype_num}",
        f"serotype{serotype_num}",
        f"DENV{serotype_num}_E"
    ]

    # Try different possible file names
    possible_files = [
        f"DENV{serotype_num}_E_proteins_complete.fasta",
        f"DENV{serotype_num}_E_proteins_verified.fasta",
        f"DENV{serotype_num}_E_proteins.fasta",
        f"denv{serotype_num}_E_complete.fasta",
        f"serotype{serotype_num}_E.fasta",
        "E_proteins.fasta",
        "complete.fasta",
        "verified.fasta"
    ]

    # Search for files
    fasta_file = None
    used_folder = None

    for folder in possible_folders:
        if os.path.exists(folder):
            for file_name in possible_files:
                file_path = os.path.join(folder, file_name)
                if os.path.exists(file_path):
                    fasta_file = file_path
                    used_folder = folder
                    print(f" DENV{serotype_num}: Found {file_name} in {folder}/")
                    break
        if fasta_file:
            break

    if not fasta_file:
        print(f" DENV{serotype_num}: No FASTA file found")
        return f"DENV{serotype_num}", 0.0, 0, 0

    # Count sequences and check motif
    total = 0
    with_motif = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        total += 1
        if MOTIF in str(record.seq):
            with_motif += 1

    percentage = (with_motif / total * 100) if total > 0 else 0.0

    print(f"  DENV{serotype_num}: {with_motif}/{total} sequences have motif ({percentage:.1f}%)")

    return f"DENV{serotype_num}", percentage, total, with_motif


def main():
    print("=" * 60)
    print(f"CHECKING MOTIF: {MOTIF}")
    print("=" * 60)

    # Check all 4 serotypes
    results = []

    for serotype in [1, 2, 3, 4]:
        result = check_serotype(serotype)
        results.append(result)

    print("\n" + "=" * 60)
    print("CREATING SUMMARY FILE")
    print("=" * 60)

    # Create the single txt file
    with open(OUTPUT_FILE, "w") as f:
        f.write("=" * 60 + "\n")
        f.write(f"MOTIF PRESENCE IN DENV SEROTYPES\n")
        f.write(f"Motif: {MOTIF}\n")
        f.write("=" * 60 + "\n\n")

        f.write("SEROTYPE\tTOTAL_SEQ\tWITH_MOTIF\tPERCENTAGE\n")
        f.write("-" * 60 + "\n")

        for serotype_name, percentage, total, with_motif in results:
            f.write(f"{serotype_name}\t{total}\t{with_motif}\t{percentage:.1f}%\n")

        # Add summary
        f.write("\n" + "=" * 60 + "\n")
        f.write("SUMMARY\n")
        f.write("=" * 60 + "\n")

        # Find which serotypes have highest/lowest presence
        valid_results = [(name, perc, total, count) for name, perc, total, count in results if total > 0]

        if valid_results:
            # Sort by percentage
            sorted_by_perc = sorted(valid_results, key=lambda x: x[1], reverse=True)

            highest = sorted_by_perc[0]
            lowest = sorted_by_perc[-1]

            f.write(f"\nHighest motif presence: {highest[0]} ({highest[1]:.1f}%)\n")
            f.write(f"Lowest motif presence: {lowest[0]} ({lowest[1]:.1f}%)\n")

            # Calculate overall average
            total_sequences = sum(r[2] for r in valid_results)
            total_with_motif = sum(r[3] for r in valid_results)
            overall_percentage = (total_with_motif / total_sequences * 100) if total_sequences > 0 else 0

            f.write(f"\nOverall across all serotypes:\n")
            f.write(f"  Total sequences analyzed: {total_sequences}\n")
            f.write(f"  Total with motif: {total_with_motif}\n")
            f.write(f"  Overall percentage: {overall_percentage:.1f}%\n")

    print(f"\n Summary saved to: {OUTPUT_FILE}")

    # Also print the results to console
    print(f"\n{'=' * 60}")
    print("FINAL RESULTS")
    print(f"{'=' * 60}")

    for serotype_name, percentage, total, with_motif in results:
        if total > 0:
            print(f"{serotype_name}: {with_motif}/{total} = {percentage:.1f}%")
        else:
            print(f"{serotype_name}: No sequences found")


if __name__ == "__main__":
    main()
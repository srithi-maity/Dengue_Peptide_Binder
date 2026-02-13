import os
from Bio import SeqIO
from Bio.Seq import Seq

# MOTIF = "VDRGWGNGCGLFGKG"
MOTIF = "VVVLGSQEGAMHTALTGATEIQ"
OUTPUT_FILE = "motif2_counts.txt"

results = []

for serotype in [1, 2, 3, 4]:
    fasta = f"denv{serotype}/DENV{serotype}_E_proteins_verified.fasta"

    if not os.path.exists(fasta):
        print(f" File not found: denv{serotype}/DENV{serotype}_E_proteins_verified.fasta")
        results.append(f"DENV{serotype}\tFILE_NOT_FOUND")
        continue

    total = 0
    with_motif = 0
    total_occurrences = 0
    multiple = []

    for record in SeqIO.parse(fasta, "fasta"):
        total += 1
        seq = str(record.seq)
        count = Seq(seq).count_overlap(MOTIF)

        if count > 0:
            with_motif += 1
            total_occurrences += count
        if count > 1:
            multiple.append(f"{record.id}({count})")

    percentage = (with_motif / total * 100) if total > 0 else 0
    avg = total_occurrences / with_motif if with_motif > 0 else 0
    results.append(
        f"DENV{serotype}\t{total}\t{with_motif}\t{percentage:.1f}%\t{total_occurrences}\t{avg:.1f}\t{', '.join(multiple)}")

# Save results
with open(OUTPUT_FILE, "w") as f:
    f.write("SEROTYPE\tTOTAL\tWITH_MOTIF\t%\tTOTAL_OCC\tAVG\tMULTIPLE\n")
    f.write("-" * 100 + "\n")
    for line in results:
        f.write(line + "\n")

print(f"\n Done! Results saved to {OUTPUT_FILE}")
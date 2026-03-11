import os

folder = "."   # current folder

# valid amino acids accepted by AlphaFold
valid_aa = set("ACDEFGHIKLMNPQRSTVWYX")

print("\nScanning FASTA files for invalid amino acids...\n")

for file in os.listdir(folder):
    if file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".faa"):

        path = os.path.join(folder, file)

        sequence = ""

        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith(">"):   # ignore FASTA header
                    sequence += line.strip()

        print("====================================")
        print("File:", file)
        print("Sequence length:", len(sequence))
        print("Unique characters:", sorted(set(sequence)))

        invalid_positions = []

        for i, aa in enumerate(sequence, start=1):
            if aa not in valid_aa:
                invalid_positions.append((i, aa))

        if invalid_positions:
            print("\n❌ Invalid characters detected:")
            for pos, aa in invalid_positions:
                print(f"Position {pos}: '{aa}'")

        else:
            print("\n✅ Sequence contains only valid amino acids")

print("\nScan finished.")
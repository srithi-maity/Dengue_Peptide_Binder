import os
import subprocess

input_folder = "."
base_output = "alphafold_results"

os.makedirs(base_output, exist_ok=True)

print("Scanning FASTA files...\n")

for file in os.listdir(input_folder):

    if file.endswith(".fasta"):

        fasta_path = os.path.join(input_folder, file)

        name = file.replace(".fasta", "")
        output_folder = os.path.join(base_output, name)

        os.makedirs(output_folder, exist_ok=True)

        print("=================================")
        print("Running AlphaFold for:", name)

        command = ["colabfold_batch", fasta_path, output_folder]

        subprocess.run(command)

        print("Finished:", name)
        print("Results saved in:", output_folder)

print("\nAll predictions completed.")
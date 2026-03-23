# Chimera script for denv_scaffold_11_complex
# Colors: Target (Chain A) - Blue, Binder (Chain B) - Green, Motif (97-111) - Red

from chimera import openModels, runCommand

# Load structure
openModels.open("denv_scaffold_11_complex.pdb")

# Color by chain
runCommand("color blue #0:.A")
runCommand("color green #0:.B")

# Color motif residues
runCommand("color red #0:.A,97-111")

# Show representations
runCommand("repr stick #0:.B")
runCommand("repr stick #0:.A,97-111")

# Highlight motif with spheres
runCommand("repr sphere #0:.A,97-111")
runCommand("setattr m sphereScale 0.3 #0:.A,97-111")

# Label motif
runCommand("labelopt info resName,resNumber")
runCommand("label #0:.A,97-111")

# Save image
runCommand("copy file denv_scaffold_11_complex.png size 1600 1200")

print("=== Visualization Summary ===")
print("Target (Chain A): Blue")
print("Binder (Chain B): Green")
print(f"Motif (97-111): Red")
print("==============================")

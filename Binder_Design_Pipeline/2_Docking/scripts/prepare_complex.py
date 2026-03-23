# import os
#
# RECEPTOR = "../receptors/denv1_clean.pdb"
# LIGAND_DIR = "../ligands"
# OUT_DIR = "../complexes"
#
# os.makedirs(OUT_DIR, exist_ok=True)
#
# ligands = [f for f in os.listdir(LIGAND_DIR) if f.endswith(".pdb")]
#
# for lig in ligands:
#     lig_path = os.path.join(LIGAND_DIR, lig)
#     out_path = os.path.join(OUT_DIR, lig.replace(".pdb", "_complex.pdb"))
#
#     with open(out_path, "w") as out:
#         # Write receptor as chain A
#         for line in open(RECEPTOR):
#             if line.startswith("ATOM"):
#                 out.write(line[:21] + "A" + line[22:])
#
#         # Write ligand as chain B
#         for line in open(lig_path):
#             if line.startswith("ATOM"):
#                 out.write(line[:21] + "B" + line[22:])
#
#     print("Created:", out_path)
#
#
##################################################################################################################################

# import os
#
# RECEPTOR = "../receptors/denv1_clean.pdb"
# LIGAND_DIR = "../ligands"
# OUT_DIR = "../complexes"
#
# os.makedirs(OUT_DIR, exist_ok=True)
#
# ligands = [f for f in os.listdir(LIGAND_DIR) if f.endswith(".pdb")]
#
# # Define your motif residues
# motif_residues = "97-111"
#
# for lig in ligands:
#     lig_path = os.path.join(LIGAND_DIR, lig)
#     complex_name = lig.replace(".pdb", "_complex")
#     out_path = os.path.join(OUT_DIR, complex_name + ".pdb")
#     pml_path = os.path.join(OUT_DIR, complex_name + ".pml")
#     rasmac_path = os.path.join(OUT_DIR, complex_name + ".ras")
#
#     # Write PDB complex
#     with open(out_path, "w") as out:
#         for line in open(RECEPTOR):
#             if line.startswith("ATOM"):
#                 out.write(line[:21] + "A" + line[22:])
#         for line in open(lig_path):
#             if line.startswith("ATOM"):
#                 out.write(line[:21] + "B" + line[22:])
#
#     # Write PyMOL visualization script
#     with open(pml_path, "w") as pml:
#         pml.write(f"""# Load the complex
# load {complex_name}.pdb
#
# # Color settings
# color skyblue, chain A        # Target (receptor)
# color lime, chain B           # Binder (ligand)
#
# # Color the specific motif (residues 97-111) in red
# select motif, chain A and resi {motif_residues}
# color red, motif
# show sticks, motif
# set stick_radius, 0.2, motif
#
# # Highlight the motif with a transparent surface
# show surface, motif
# set transparency, 0.5, motif
#
# # Visual settings for better clarity
# show cartoon, chain A
# show sticks, chain B
# set cartoon_transparency, 0.2, chain A
# bg_color white
#
# # Label the motif
# label motif, "Motif 97-111"
# set label_color, black
# set label_size, 12
#
# # Zoom to the binding interface
# zoom chain B
# zoom chain A, 10
#
# # Save image
# png {complex_name}.png, width=1200, height=900, dpi=300
#
# print("--- Visualization ready ---")
# print("Motif residues 97-111 are colored red")
# print("Target (receptor): blue")
# print("Binder (ligand): green")
# """)
#
#     # Write a simple coloring script for PyMOL (alternative)
#     with open(rasmac_path, "w") as ras:
#         ras.write(f"""# Simple coloring scheme
# load {complex_name}.pdb
# color blue, chain A
# color green, chain B
# color red, chain A and resi {motif_residues}
# show cartoon, chain A
# show sticks, chain B
# show sticks, chain A and resi {motif_residues}
# zoom
# """)
#
#     print(f"Created: {out_path}")
#     print(f"Created: {pml_path}")
#     print(f"Created: {rasmac_path}")
#     print(f"  - Motif residues {motif_residues} will be colored RED")
#     print()
#
# print("=" * 50)
# print("All complexes generated successfully!")
# print("To visualize in PyMOL, run: pymol complex_name.pml")



#######################################################################################################################
import os

RECEPTOR = "../receptors/denv1_clean.pdb"
LIGAND_DIR = "../ligands"
OUT_DIR = "../complexes"

os.makedirs(OUT_DIR, exist_ok=True)

ligands = [f for f in os.listdir(LIGAND_DIR) if f.endswith(".pdb")]

# Your motif residues
motif_residues = "97-111"

for lig in ligands:
    lig_path = os.path.join(LIGAND_DIR, lig)
    complex_name = lig.replace(".pdb", "_complex")
    out_path = os.path.join(OUT_DIR, complex_name + ".pdb")
    chimera_path = os.path.join(OUT_DIR, complex_name + ".py")  # Chimera script
    chimerax_path = os.path.join(OUT_DIR, complex_name + ".cxc")  # ChimeraX script

    # Write PDB complex
    with open(out_path, "w") as out:
        for line in open(RECEPTOR):
            if line.startswith("ATOM"):
                out.write(line[:21] + "A" + line[22:])
        for line in open(lig_path):
            if line.startswith("ATOM"):
                out.write(line[:21] + "B" + line[22:])

    # ChimeraX script (.cxc)
    with open(chimerax_path, "w") as cx:
        cx.write(f"""# ChimeraX script for {complex_name}
# Colors: Target (Chain A) - Blue, Binder (Chain B) - Green, Motif (97-111) - Red

open {complex_name}.pdb

# Color by chain
color chain A #90caf9      # Light blue for target
color chain B #66bb6a      # Light green for binder

# Color motif (residues 97-111) in red
color /A:{motif_residues} red

# Show representations
cartoon style /A helix
show /A cartoon
show /B sticks

# Highlight motif
show /A:{motif_residues} sticks
show /A:{motif_residues} spheres
size /A:{motif_residues} sphereScale 0.3

# Label motif
label /A:{motif_residues} text "Motif 97-111"
label /A:{motif_residues} textColor black
label /A:{motif_residues} bgColor white

# Focus on binding site
view /B
view /A:{motif_residues}

# Save image
save {complex_name}.png width 1600 height 1200 transparentBackground false

# Print info
log info "=== Visualization Summary ==="
log info "Target (Chain A): Blue"
log info "Binder (Chain B): Green"
log info "Motif ({motif_residues}): Red"
log info "=============================="
""")

    # Chimera script (.py) for legacy Chimera
    with open(chimera_path, "w") as ch:
        ch.write(f"""# Chimera script for {complex_name}
# Colors: Target (Chain A) - Blue, Binder (Chain B) - Green, Motif (97-111) - Red

from chimera import openModels, runCommand

# Load structure
openModels.open("{complex_name}.pdb")

# Color by chain
runCommand("color blue #0:.A")
runCommand("color green #0:.B")

# Color motif residues
runCommand("color red #0:.A,{motif_residues}")

# Show representations
runCommand("repr stick #0:.B")
runCommand("repr stick #0:.A,{motif_residues}")

# Highlight motif with spheres
runCommand("repr sphere #0:.A,{motif_residues}")
runCommand("setattr m sphereScale 0.3 #0:.A,{motif_residues}")

# Label motif
runCommand("labelopt info resName,resNumber")
runCommand("label #0:.A,{motif_residues}")

# Save image
runCommand("copy file {complex_name}.png size 1600 1200")

print("=== Visualization Summary ===")
print("Target (Chain A): Blue")
print("Binder (Chain B): Green")
print(f"Motif ({motif_residues}): Red")
print("==============================")
""")

    print(f"Created: {out_path}")
    print(f"  ChimeraX script: {chimerax_path}")
    print(f"  Chimera script: {chimera_path}")

print("\n" + "=" * 50)
print("Visualization in Chimera/ChimeraX:")
print("  🔵 Target (Chain A): Blue")
print("  🟢 Binder (Chain B): Green")
print(f"  🔴 Motif (Residues {motif_residues}): Red")
print("=" * 50)
print("\nTo visualize in ChimeraX:")
print("  open complex_name.cxc")
print("\nTo visualize in Chimera:")
print("  chimera complex_name.py")
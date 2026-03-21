import os
import subprocess
from Bio.PDB import PDBParser, Superimposer
import numpy as np

# ==============================
# PATH SETUP
# ==============================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ALPHAFOLD_DIR = os.path.join(BASE_DIR, "alphafold_str")
RFDIFF_DIR = os.path.join(BASE_DIR, "RFdiffusion")

RF_SCRIPT = os.path.join(RFDIFF_DIR, "scripts", "run_inference.py")

OUTPUT_PREFIX = "outputs/denv_scaffold"

# Cleaned PDBs (IMPORTANT)
SEROTYPES = [
    os.path.join(ALPHAFOLD_DIR, "denv1_clean.pdb"),
    os.path.join(ALPHAFOLD_DIR, "denv2_clean.pdb"),
    os.path.join(ALPHAFOLD_DIR, "denv3_clean.pdb"),
    os.path.join(ALPHAFOLD_DIR, "denv4_clean.pdb"),
]

INPUT_PDB = SEROTYPES[0]

# ==============================
# PARAMETERS
# ==============================

MOTIF_START = 97
MOTIF_END = 111
CHAIN = "A"

NUM_DESIGNS = 200
RMSD_THRESHOLD = 2.0

# ==============================
# STEP 1: RFdiffusion
# ==============================

def run_rfdiffusion():
    print("\nRunning RFdiffusion...")

    contig = f"[{CHAIN}{MOTIF_START}-{MOTIF_END}/0 80-100]"

    cmd = f"""
    python {RF_SCRIPT} \
    inference.input_pdb={INPUT_PDB} \
    inference.output_prefix={OUTPUT_PREFIX} \
    contigmap.contigs='{contig}' \
    inference.num_designs={NUM_DESIGNS} \
    inference.model_directory_path=models \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0
    """

    os.chdir(RFDIFF_DIR)
    os.system(cmd)

# ==============================
# STEP 2: RMSD FUNCTIONS
# ==============================

parser = PDBParser(QUIET=True)

def get_atoms(struct):
    atoms = []
    for model in struct:
        for chain in model:
            if chain.id == CHAIN:
                for res in chain:
                    if res.id[0] == ' ' and MOTIF_START <= res.id[1] <= MOTIF_END:
                        if "CA" in res:
                            atoms.append(res["CA"])
    return atoms

def rmsd(a1, a2):
    if len(a1) != len(a2):
        return 999
    sup = Superimposer()
    sup.set_atoms(a1, a2)
    return sup.rms


# ==============================
# MAIN
# ==============================

if __name__ == "__main__":
    print("\nDENGUE SCAFFOLD DESIGN PIPELINE")
    run_rfdiffusion()
    print("\nDONE")
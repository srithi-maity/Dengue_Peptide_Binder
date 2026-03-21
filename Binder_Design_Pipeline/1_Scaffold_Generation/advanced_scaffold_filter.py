import os
import shutil
from Bio.PDB import PDBParser, Superimposer
import numpy as np

# ==============================
# PATH SETUP
# ==============================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ALPHAFOLD_DIR = os.path.join(BASE_DIR, "alphafold_str")
RFDIFF_DIR = os.path.join(BASE_DIR, "RFdiffusion")

OUTPUT_DIR = os.path.join(RFDIFF_DIR, "outputs")

# ==============================
# INPUT FILES
# ==============================

SEROTYPES = [
    os.path.join(ALPHAFOLD_DIR, "denv1_clean.pdb"),
    os.path.join(ALPHAFOLD_DIR, "denv2_clean.pdb"),
    os.path.join(ALPHAFOLD_DIR, "denv3_clean.pdb"),
    os.path.join(ALPHAFOLD_DIR, "denv4_clean.pdb"),
]

# ==============================
# PARAMETERS
# ==============================

MOTIF_START = 97
MOTIF_END = 111
CHAIN = "A"

MAX_RMSD_THRESHOLD = 0.5
TOP_N = 20

# ==============================
# RMSD FUNCTIONS
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
# FILTER + SCORE
# ==============================

def filter_scaffolds():

    refs = []
    for pdb in SEROTYPES:
        s = parser.get_structure("ref", pdb)
        refs.append(get_atoms(s))

    files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith(".pdb")]

    results = []

    for f in files:
        path = os.path.join(OUTPUT_DIR, f)

        try:
            s = parser.get_structure("scaffold", path)
            atoms = get_atoms(s)

            if len(atoms) == 0:
                continue

            scores = []
            for ref in refs:
                n = min(len(atoms), len(ref))
                scores.append(rmsd(atoms[:n], ref[:n]))

            avg = np.mean(scores)
            max_r = np.max(scores)
            std = np.std(scores)

            if max_r > MAX_RMSD_THRESHOLD:
                continue

            results.append((f, avg, max_r, std, scores))

        except:
            continue

    results.sort(key=lambda x: x[1])
    return results

# ==============================
# SAVE CSV
# ==============================

def save_results(results):

    out1 = os.path.join(BASE_DIR, "all_scaffold_scores.csv")
    out2 = os.path.join(BASE_DIR, "top_scaffolds_detailed.csv")

    with open(out1, "w") as f:
        f.write("File,AvgRMSD,MaxRMSD,StdDev,DENV1,DENV2,DENV3,DENV4\n")
        for r in results:
            f.write(
                f"{r[0]},{r[1]:.4f},{r[2]:.4f},{r[3]:.4f},"
                f"{r[4][0]:.4f},{r[4][1]:.4f},"
                f"{r[4][2]:.4f},{r[4][3]:.4f}\n"
            )

    with open(out2, "w") as f:
        f.write("File,AvgRMSD,MaxRMSD,StdDev,DENV1,DENV2,DENV3,DENV4\n")
        for r in results[:TOP_N]:
            f.write(
                f"{r[0]},{r[1]:.4f},{r[2]:.4f},{r[3]:.4f},"
                f"{r[4][0]:.4f},{r[4][1]:.4f},"
                f"{r[4][2]:.4f},{r[4][3]:.4f}\n"
            )

    print("CSV files saved")

# ==============================
# COPY TOP STRUCTURES
# ==============================

def save_top_structures(results):

    target_dir = os.path.join(BASE_DIR, "top_scaffolds_structures")

    # Create folder
    os.makedirs(target_dir, exist_ok=True)

    print(f"\nCopying top {TOP_N} scaffolds to: {target_dir}")

    for r in results[:TOP_N]:
        filename = r[0]

        src = os.path.join(OUTPUT_DIR, filename)
        dst = os.path.join(target_dir, filename)

        try:
            shutil.copy(src, dst)
        except:
            continue

    print("Top scaffold structures copied successfully")

# ==============================
# MAIN
# ==============================

if __name__ == "__main__":

    print("\nFiltering existing scaffolds...")

    results = filter_scaffolds()

    save_results(results)

    save_top_structures(results)

    print("\nDONE")
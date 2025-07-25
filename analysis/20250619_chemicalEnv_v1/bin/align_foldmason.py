"""
Align with foldMason
salloc --mem 32GB -t 2:00:00 --gres=gpu:a100:1
conda activate mmseqs2
python3 bin/align_foldmason.py 
"""
import os
import subprocess

STRUCTURE_DIR = "/hps/nobackup/agb/research/francesco/20250624_isopep_bond_domains/cropped_structures"
OUTPUT_DIR = "output/MSAs"
TMP_DIR = "tmp"

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(TMP_DIR, exist_ok=True)
    for pfamA_acc in os.listdir(STRUCTURE_DIR):
        print(f"Running domain {pfamA_acc}")
        try:
            result = subprocess.run(
                [
                    "foldmason",
                    "easy-msa",
                    os.path.join(STRUCTURE_DIR, pfamA_acc),
                    os.path.join(OUTPUT_DIR, f"{pfamA_acc}"),
                    "tmp",
                    "--gpu", "1",
                    "--report-mode", "0",
                    "--refine-seed", "42",
                    "--refine-iters", "100"
                ],
                check=True,
                capture_output=True,
                text=True
            )
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            print("Error occurred while running foldmason:")
            print("Return code:", e.returncode)
            print("Standard output:", e.stdout)
            print("Standard error:", e.stderr)

if __name__ == "__main__":
    main()
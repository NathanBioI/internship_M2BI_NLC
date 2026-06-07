#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapping KMA complémentaire sur la base NCBI Borrelia.

Entrée :
  data/fastq/*_R1.fastq.gz + *_R2.fastq.gz
  data/ncbi/borrelia_ncbi_kma_index.*

Sortie :
  results/ncbi_typing/<sample>/<sample>.res, .fsa, ...

Le seuil -ID est fixé à 95 par défaut pour permettre deux lectures ensuite :
  - stricte : 100 % identité, 100 % couverture, profondeur >= 5
  - relâchée : identité >= 95 %, profondeur > 300
"""

import argparse
import re
import subprocess
from pathlib import Path

def sample_name(r1: Path) -> str:
    """847_S54_L001_trimmed_R1.fastq.gz -> 847."""
    return re.sub(r"_R1\.fastq\.gz$", "", r1.name).split("_", 1)[0]

def r2_from_r1(r1: Path) -> Path:
    return Path(str(r1).replace("_R1.fastq.gz", "_R2.fastq.gz"))

def run_kma(r1: Path, r2: Path, db: str, outprefix: Path, identity: float) -> None:
    cmd = [
        "kma",
        "-ipe", str(r1), str(r2),
        "-o", str(outprefix),
        "-t_db", db,
        "-1t1",
        "-ID", str(identity),
        "-and",
        "-ref_fsa",
        "-oa",
        "-ef",
    ]
    subprocess.run(cmd, check=True)

def main() -> None:
    parser = argparse.ArgumentParser(description="Mapping KMA sur la base NCBI Borrelia.")
    parser.add_argument("--fastq-dir", default="data/fastq")
    parser.add_argument("--db", default="data/ncbi/borrelia_ncbi_kma_index")
    parser.add_argument("--out-root", default="results/ncbi_typing")
    parser.add_argument("--identity", type=float, default=95.0)
    args = parser.parse_args()

    fastq_dir = Path(args.fastq_dir)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    n_done = 0

    for r1 in sorted(fastq_dir.glob("*_R1.fastq.gz")):
        sample = sample_name(r1)
        r2 = r2_from_r1(r1)

        sample_dir = out_root / sample
        sample_dir.mkdir(parents=True, exist_ok=True)

        print(f"[RUN] {sample}")
        run_kma(r1, r2, args.db, sample_dir / sample, args.identity)
        n_done += 1

    print(f"\nÉchantillons traités : {n_done}")
    print(f"Résultats dans : {out_root}")

if __name__ == "__main__":
    main()

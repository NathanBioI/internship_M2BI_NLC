#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Concaténation MLST locus par locus.

Entrée  : alignements FASTA MAFFT, un fichier par locus.
Sorties : 1) FASTA concaténé pour IQ-TREE
          2) fichier de partitions au format RAxML/IQ-TREE.

Convention importante :
les headers produits par typing-mapping_with_kma.py sont du type
>sample|locus_allele. Pour l'arbre concaténé, on ne garde que "sample".
"""

import argparse
from pathlib import Path

def read_fasta(path: Path) -> dict[str, str]:
    """Lit un FASTA et retourne {taxon: sequence}."""
    seqs, name, buf = {}, None, []

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                header = line[1:].split()[0]
                name = header.split("|")[0]
                buf = []
            else:
                buf.append(line)

    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs

def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

def locus_from_filename(path: Path) -> str:
    """clpA.aln.fasta -> clpA."""
    return path.name.split(".")[0]

def main() -> None:
    parser = argparse.ArgumentParser(description="Concatène les alignements MLST et écrit les partitions IQ-TREE.")
    parser.add_argument("aln", nargs="+", help="FASTA alignés, un par locus")
    parser.add_argument("-o", "--out", required=True, help="FASTA concaténé")
    parser.add_argument("-p", "--part", required=True, help="Fichier de partitions")
    parser.add_argument("--missing", choices=["N", "-"], default="N", help="Caractère utilisé pour les loci absents")
    args = parser.parse_args()

    per_locus = []
    taxa_order = []

    # Lecture des loci et conservation de l'ordre d'apparition des taxons.
    for file_name in args.aln:
        fasta = Path(file_name)
        seqs = read_fasta(fasta)
        length = len(next(iter(seqs.values()))) if seqs else 0
        per_locus.append((locus_from_filename(fasta), length, seqs))

        for taxon in seqs:
            if taxon not in taxa_order:
                taxa_order.append(taxon)

    # Concaténation. Si un taxon n'a pas un locus, on remplit avec N.
    concatenated = {taxon: "" for taxon in taxa_order}
    partitions = []
    start = 1

    for locus, length, seqs in per_locus:
        end = start + length - 1
        partitions.append((locus, start, end))
        filler = args.missing * length

        for taxon in taxa_order:
            concatenated[taxon] += seqs.get(taxon, filler)

        start = end + 1

    out_fasta = Path(args.out)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w", encoding="utf-8") as handle:
        for taxon, seq in concatenated.items():
            handle.write(f">{taxon}\n{wrap(seq)}\n")

    out_part = Path(args.part)
    out_part.parent.mkdir(parents=True, exist_ok=True)
    with out_part.open("w", encoding="utf-8") as handle:
        for locus, start, end in partitions:
            handle.write(f"DNA,{locus}={start}-{end}\n")

if __name__ == "__main__":
    main()

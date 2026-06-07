#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Construit le fichier minimal de spécificité allèle/espèce utilisé par le typage PubMLST.

Entrée attendue :
  data/pubMLST_profile/BIGSdb_3429203_6349790279_07966.csv

Sortie :
  data/pubMLST_profile/mlst_species_specificity/all_alleles_specificity.tsv

Le fichier produit indique, pour chaque couple locus/allèle, l'espèce majoritaire
dans BIGSdb, le nombre d'isolats, la pureté et le détail par espèce.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

LOCI = ["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"]
DEFAULT_INPUT = Path("data/pubMLST_profile/BIGSdb_3429203_6349790279_07966.csv")
DEFAULT_OUTPUT = Path("data/pubMLST_profile/mlst_species_specificity/all_alleles_specificity.tsv")

FIELDS = [
    "locus", "allele", "major_species", "major_n", "total", "n_species",
    "purity", "purity_pct", "is_multi_species", "species_breakdown",
]

def clean(value) -> str:
    value = "" if value is None else str(value).strip()
    return "" if value.lower() in {"", "nan", "na", "none"} else value

def sniff_delimiter(path: Path) -> str:
    """Détection simple CSV/TSV pour l'export BIGSdb."""
    sample = path.read_text(encoding="utf-8-sig", errors="replace")[:8192]
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        return "\t"

def allele_key(value: str):
    """Tri naturel : 2 avant 10, puis texte si l'allèle n'est pas numérique."""
    try:
        return (0, int(value))
    except ValueError:
        return (1, value)

def read_bigsdb(path: Path):
    """Lit seulement les colonnes species + loci MLST."""
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
        yield from csv.DictReader(handle, delimiter=sniff_delimiter(path))

def build_specificity(rows) -> list[dict[str, str]]:
    counts = {locus: defaultdict(lambda: defaultdict(int)) for locus in LOCI}

    # Comptage : pour chaque locus/allèle, combien d'isolats par espèce.
    for row in rows:
        species = clean(row.get("species"))
        if not species:
            continue
        for locus in LOCI:
            allele = clean(row.get(locus))
            if allele:
                counts[locus][allele][species] += 1

    output = []
    for locus in LOCI:
        for allele in sorted(counts[locus], key=allele_key):
            species_counts = counts[locus][allele]
            ordered = sorted(species_counts.items(), key=lambda x: (-x[1], x[0]))

            major_species, major_n = ordered[0]
            total = sum(species_counts.values())
            n_species = len(species_counts)
            purity = major_n / total

            output.append({
                "locus": locus,
                "allele": allele,
                "major_species": major_species,
                "major_n": str(major_n),
                "total": str(total),
                "n_species": str(n_species),
                "purity": f"{purity:.12g}",
                "purity_pct": f"{100 * purity:.2f}",
                "is_multi_species": "True" if n_species > 1 else "False",
                "species_breakdown": "; ".join(f"{sp} ({n})" for sp, n in ordered),
            })

    return output

def write_tsv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

def main() -> None:
    parser = argparse.ArgumentParser(description="Construit all_alleles_specificity.tsv depuis l'export BIGSdb.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--out", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    rows = build_specificity(read_bigsdb(Path(args.input)))
    write_tsv(rows, Path(args.out))

    if not args.quiet:
        print(f"Fichier écrit : {args.out}")
        print(f"Allèles résumés : {len(rows)}")

if __name__ == "__main__":
    main()

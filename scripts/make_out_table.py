#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Construit le tableau du taux d'infection à partir du TSV final d'identification.

Entrée :
  results/metadata_genaspe_species_final.tsv

Sortie :
  results/table_infection.tsv

"""

import argparse
import csv
from pathlib import Path

SPECIES = [
    "Borrelia afzelii",
    "Borrelia burgdorferi sensu stricto",
    "Borrelia garinii",
    "Borrelia lusitaniae",
    "Borrelia valaisiana",
    "Borrelia miyamotoi",
]

BBSL_SPECIES = [
    "Borrelia afzelii",
    "Borrelia burgdorferi sensu stricto",
    "Borrelia garinii",
    "Borrelia lusitaniae",
    "Borrelia valaisiana",
]

BASE_COLUMNS = ["n°", "n°ADN EPIA", "Site", "date collecte"]
OUT_COLUMNS = BASE_COLUMNS + SPECIES + ["Bbsl", "sp inconnu"]

def present(species_final: str, species: str) -> int:
    return int(species in (species_final or ""))

def main() -> None:
    parser = argparse.ArgumentParser(description="Construit la table de prévalence depuis metadata_genaspe_species_final.tsv.")
    parser.add_argument("--input", default="results/metadata_genaspe_species_final.tsv")
    parser.add_argument("--out", default="results/table_infection.tsv")
    args = parser.parse_args()

    input_tsv = Path(args.input)
    output_tsv = Path(args.out)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    out_rows = []
    counts = {col: 0 for col in SPECIES + ["Bbsl", "sp inconnu"]}

    with input_tsv.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            # On conserve uniquement les tiques positives Borrelia.
            if (row.get("OK Borrelia", "") or "").strip().lower() != "x":
                continue

            species_final = row.get("species_final_retained", "") or "inconnu"

            out = {
                "n°": row.get("n°", ""),
                "n°ADN EPIA": row.get("n°ADN EPIA", ""),
                "Site": row.get("Type échantillon", ""),
                "date collecte": row.get("date collecte", ""),
            }

            for species in SPECIES:
                out[species] = present(species_final, species)

            out["Bbsl"] = int(any(out[species] for species in BBSL_SPECIES))
            out["sp inconnu"] = int(species_final.strip().lower() == "inconnu")

            for col in counts:
                counts[col] += int(out[col])

            out_rows.append(out)

    with output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Table générée : {output_tsv}")
    print(f"Nombre de lignes positives Borrelia : {len(out_rows)}")
    for col, n in counts.items():
        print(f"{col}\t{n}")

if __name__ == "__main__":
    main()

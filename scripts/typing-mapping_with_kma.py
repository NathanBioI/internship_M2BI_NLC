#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Typage MLST PubMLST avec KMA.

Ce script fait trois choses :
  1. mappe chaque paire FASTQ sur la base d'allèles PubMLST ;
  2. extrait les meilleurs allèles, les allèles parfaits 100/100 et le statut ST ;
  3. écrit les FASTA par locus utilisés ensuite pour MAFFT et l'arbre MLST.

Sorties principales :
  results/pubMLST_typing/<sample>/alleles.tsv
  results/pubMLST_typing/<sample>/st_sp.tsv
  results/pubMLST_typing/<sample>/perfect_alleles_specificity.tsv
  results/pubMLST_typing/mafft_input/<locus>.fasta
"""

import csv
import subprocess
from pathlib import Path

LOCI = ["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"]

FASTQ_DIR = Path("data/fastq")
PUBMLST_DB = "data/pubMLST_alleles/borrelia_burgdorferi_sensu_lato"
PROFILES_ST = Path("data/pubMLST_profile/borrelia_spp")
BIGSDB = Path("data/pubMLST_profile/BIGSdb_3429203_6349790279_07966.csv")
SPECIFICITY_FILE = Path("data/pubMLST_profile/mlst_species_specificity/all_alleles_specificity.tsv")
OUT_DIR = Path("results/pubMLST_typing")
MAFFT_DIR = OUT_DIR / "mafft_input"

MIN_ID = 97.0
MIN_COV = 99.0
KMA_IDENTITY = 97.0

# Seuils utilisés pour décider si un allèle parfait supporte une espèce
# lors de la scission des co-infections dans les FASTA MLST.
SPLIT_MIN_TOTAL_STRICT = 20
SPLIT_MIN_PURITY_STRICT = 0.95
SPLIT_MIN_TOTAL_RELAXED = 50
SPLIT_MIN_PURITY_RELAXED = 0.90

def r2_from_r1(r1: Path) -> Path:
    return Path(str(r1).replace("_R1.fastq.gz", "_R2.fastq.gz"))

def sample_name(r1: Path) -> str:
    """847_S54_L001_trimmed_R1.fastq.gz -> 847."""
    return r1.name.split("_", 1)[0]

def read_fasta(path: Path) -> dict[str, str]:
    """Lit le .fsa produit par KMA : {template: consensus}."""
    seqs, header, buf = {}, None, []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(buf).upper()
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if header is not None:
        seqs[header] = "".join(buf).upper()
    return seqs

def load_profiles() -> list[tuple[tuple[str, ...], str]]:
    """Charge le schéma PubMLST : combinaison d'allèles -> ST."""
    rows = []
    with PROFILES_ST.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = tuple((row.get(locus, "") or "").strip() for locus in LOCI)
            st = (row.get("ST", "") or "").strip()
            if st:
                rows.append((key, st))
    return rows

def load_st_species() -> dict[str, str]:
    """Charge l'espèce associée à chaque ST quand elle est disponible."""
    out = {}
    with BIGSDB.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            st = (row.get("ST (MLST)", "") or "").strip()
            species = (row.get("species", "") or "").strip()
            if st and species and st not in out:
                out[st] = species
    return out

def load_specificity() -> tuple[str, dict[tuple[str, str], dict[str, object]]]:
    """Charge all_alleles_specificity.tsv produit par loci_specificity.py."""
    parsed = {}
    with SPECIFICITY_FILE.open(encoding="utf-8") as handle:
        header = handle.readline().strip()
        cols = header.split("\t")
        idx = {name: cols.index(name) for name in ["locus", "allele", "major_species", "n_species", "total", "purity"]}

        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            locus = fields[idx["locus"]]
            allele = fields[idx["allele"]]
            parsed[(locus, allele)] = {
                "line": line.strip(),
                "major_species": fields[idx["major_species"]],
                "n_species": int(fields[idx["n_species"]]),
                "total": int(fields[idx["total"]]),
                "purity": float(fields[idx["purity"]]),
            }

    return header, parsed

def run_kma(r1: Path, r2: Path, outprefix: Path) -> None:
    """Lance KMA sur la base PubMLST."""
    subprocess.run([
        "kma",
        "-ipe", str(r1), str(r2),
        "-o", str(outprefix),
        "-t_db", PUBMLST_DB,
        "-1t1",
        "-ID", str(KMA_IDENTITY),
        "-and",
        "-ref_fsa",
        "-oa",
        "-ef",
    ], check=True)

def parse_kma_res(res_path: Path):
    """
    Retourne :
      - meilleur allèle par locus selon couverture, identité, score ;
      - liste de tous les allèles parfaits 100 % identité / 100 % couverture.
    """
    best = {locus: ("NA", -1.0, -1.0, -1.0) for locus in LOCI}
    perfect = {locus: [] for locus in LOCI}

    with res_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            cols = line.split()
            template = cols[0]
            score = float(cols[1])
            qid = float(cols[6])
            qcov = float(cols[7])

            if qid < MIN_ID or qcov < MIN_COV:
                continue

            for locus in LOCI:
                prefix = locus + "_"
                if template.startswith(prefix):
                    allele = template[len(prefix):]
                    if (qcov, qid, score) > best[locus][1:]:
                        best[locus] = (allele, qcov, qid, score)
                    if qid == 100.0 and qcov == 100.0:
                        perfect[locus].append(allele)
                    break

    return {locus: best[locus][0] for locus in LOCI}, perfect

def compatible_sts(alleles: dict[str, str], profiles: list[tuple[tuple[str, ...], str]]) -> list[str]:
    """Trouve les ST compatibles avec un profil incomplet 6/8 ou 7/8."""
    known = {locus: allele for locus, allele in alleles.items() if allele != "NA"}
    out, seen = [], set()

    for key, st in profiles:
        if all(key[LOCI.index(locus)] == allele for locus, allele in known.items()):
            if st not in seen:
                seen.add(st)
                out.append(st)

    return out

def assign_st_species(alleles: dict[str, str], profiles, st_to_species):
    """Assigne un ST complet, ou liste les candidats 7/8 ou 6/8."""
    n_known = sum(alleles[locus] != "NA" for locus in LOCI)

    if n_known == 8:
        key = tuple(alleles[locus] for locus in LOCI)
        st = next((s for k, s in profiles if k == key), "NOT_FOUND")
        species = st_to_species.get(st, "NA")
        status = "ASSIGNED" if st != "NOT_FOUND" else "NOT_FOUND"
        return status, st, species, "", ""

    if n_known in {6, 7}:
        candidates = compatible_sts(alleles, profiles)
        species_candidates = sorted({st_to_species.get(st, "") for st in candidates if st_to_species.get(st, "")})
        status = f"CANDIDATE_{n_known}of8"
        st = candidates[0] if len(candidates) == 1 else ("MULTIPLE" if candidates else "NONE")
        species = species_candidates[0] + " (candidat)" if len(species_candidates) == 1 else "NA"
        return status, st, species, ",".join(candidates), ",".join(species_candidates)

    return "INCOMPLETE", "INCOMPLETE", "NA", "", ""

def is_species_specific(data: dict[str, object]) -> bool:
    """Critères utilisés pour relier un allèle parfait à une espèce."""
    return (
        data["n_species"] == 1
        or (data["total"] >= SPLIT_MIN_TOTAL_STRICT and data["purity"] >= SPLIT_MIN_PURITY_STRICT)
        or (data["total"] >= SPLIT_MIN_TOTAL_RELAXED and data["purity"] >= SPLIT_MIN_PURITY_RELAXED)
    )

def species_split_map(perfect: dict[str, list[str]], specificity: dict[tuple[str, str], dict[str, object]]):
    """Si plusieurs espèces sont supportées, crée sample-1, sample-2, ..."""
    species = set()
    for locus in LOCI:
        for allele in perfect[locus]:
            data = specificity.get((locus, allele))
            if data and is_species_specific(data):
                species.add(data["major_species"])

    return {sp: i + 1 for i, sp in enumerate(sorted(species))} if len(species) > 1 else {}

def write_sample_outputs(sample_dir: Path, sample: str, alleles, perfect, specificity_header, specificity, st_row):
    """Écrit alleles.tsv, perfect_alleles_specificity.tsv et st_sp.tsv."""
    with (sample_dir / "alleles.tsv").open("w", encoding="utf-8") as handle:
        handle.write("sample\t" + "\t".join(LOCI) + "\n")
        handle.write(sample + "\t" + "\t".join(alleles[locus] for locus in LOCI) + "\n")

    with (sample_dir / "perfect_alleles_specificity.tsv").open("w", encoding="utf-8") as handle:
        handle.write(specificity_header + "\n")
        for locus in LOCI:
            for allele in perfect[locus]:
                data = specificity.get((locus, allele))
                if data:
                    handle.write(data["line"] + "\n")

    with (sample_dir / "st_sp.tsv").open("w", encoding="utf-8") as handle:
        handle.write("sample\tstatus\tST\tspecies\tST_candidates\tspecies_candidates\n")
        handle.write(f"{sample}\t" + "\t".join(st_row) + "\n")

def append_locus_fastas(sample: str, alleles, perfect, consensus, specificity, split_map) -> None:
    """
    Écrit les séquences utilisées pour l'arbre MLST.
    On conserve le meilleur allèle et tous les allèles parfaits.
    En cas de co-infection, les séquences sont réparties en sample-1, sample-2...
    """
    for locus in LOCI:
        alleles_to_write = sorted({alleles[locus], *perfect[locus]} - {"NA"})

        for allele in alleles_to_write:
            template = f"{locus}_{allele}"
            seq = consensus.get(template, "")
            if not seq:
                continue

            fasta_id = sample
            data = specificity.get((locus, allele))
            if split_map and data and data["major_species"] in split_map:
                fasta_id = f"{sample}-{split_map[data['major_species']]}"

            with (MAFFT_DIR / f"{locus}.fasta").open("a", encoding="utf-8") as handle:
                handle.write(f">{fasta_id}|{template}\n{seq}\n")

def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MAFFT_DIR.mkdir(parents=True, exist_ok=True)

    # Nettoyage des FASTA par locus pour éviter d'empiler les séquences à chaque relance.
    for locus in LOCI:
        (MAFFT_DIR / f"{locus}.fasta").write_text("", encoding="utf-8")

    profiles = load_profiles()
    st_to_species = load_st_species()
    specificity_header, specificity = load_specificity()

    n_done = 0

    for r1 in sorted(FASTQ_DIR.glob("*_R1.fastq.gz")):
        sample = sample_name(r1)
        r2 = r2_from_r1(r1)
        sample_dir = OUT_DIR / sample
        sample_dir.mkdir(parents=True, exist_ok=True)

        outprefix = sample_dir / sample
        run_kma(r1, r2, outprefix)

        alleles, perfect = parse_kma_res(Path(str(outprefix) + ".res"))
        consensus = read_fasta(Path(str(outprefix) + ".fsa"))
        st_row = assign_st_species(alleles, profiles, st_to_species)

        write_sample_outputs(sample_dir, sample, alleles, perfect, specificity_header, specificity, st_row)
        append_locus_fastas(sample, alleles, perfect, consensus, specificity, species_split_map(perfect, specificity))
        n_done += 1

    print(f"Échantillons PubMLST traités : {n_done}")

if __name__ == "__main__":
    main()

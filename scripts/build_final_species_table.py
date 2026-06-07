#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Construit le tableau final d'identification des espèces Borrelia.

Entrées :
  data/metadata_genaspe.csv
  results/pubMLST_typing/<sample>/
  results/ncbi_typing/<sample>/

Sortie :
  results/metadata_genaspe_species_final.tsv

Hiérarchie d'identification :
  1. PubMLST, via les allèles parfaits et leur spécificité d'espèce ;
  2. NCBI strict : Query_Identity=100, Query_Coverage=100, Depth>=5 ;
  3. NCBI relâché : Query_Identity>=95, Depth>300.
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

LOCI = {"clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"}

# Seuils PubMLST : allèle ambigu mais très majoritairement associé à une espèce.
MIN_TOTAL_FOR_MAJOR = 20
MIN_PURITY_FOR_MAJOR = 0.90

# Seuils NCBI.
NCBI_STRICT_MIN_DEPTH = 5.0
NCBI_RELAXED_MIN_QUERY_IDENTITY = 95.0
NCBI_RELAXED_MIN_DEPTH = 300.0

def clean(value: str) -> str:
    return (value or "").strip()

def is_unknown(value: str) -> bool:
    return clean(value).lower() in {"", "inconnu", "unknown", "na", "none"}

def sniff_delimiter(path: Path) -> str:
    sample = path.read_text(encoding="utf-8", errors="replace")[:4096]
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        return "\t"

def normalize_id(value: str) -> str:
    value = clean(value)
    return str(int(value)) if value.isdigit() else value

def sample_id_from_adn(value: str) -> str:
    """ADN EPPAT_23_001 -> 1."""
    match = re.search(r"(\d+)$", clean(value))
    return normalize_id(match.group(1)) if match else ""

def sample_id_from_dir(name: str) -> str:
    """847_S54... ou 847 -> 847."""
    return normalize_id(name.split("_", 1)[0])

def index_results_dirs(path: Path, skip=()) -> dict[str, Path]:
    """Indexe les dossiers résultats par numéro d'échantillon."""
    skip = set(skip)
    return {
        sample_id_from_dir(d.name): d
        for d in sorted(path.iterdir())
        if d.is_dir() and d.name not in skip
    }

def match_sample_dir(row: dict, index: dict[str, Path]) -> Path | None:
    """Apparie une ligne metadata avec un dossier de résultats."""
    for sid in (normalize_id(row.get("n°", "")), sample_id_from_adn(row.get("n°ADN EPIA", ""))):
        if sid in index:
            return index[sid]
    return None

def parse_kma_rows(res_path: Path):
    """Retourne les lignes non-commentaires d'un .res KMA sous forme de colonnes split()."""
    with res_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                yield line.split()

# ---------------------------------------------------------------------------
# PubMLST
# ---------------------------------------------------------------------------

def count_loci_100_100(res_path: Path) -> str:
    """Nombre de loci ayant au moins un hit parfait 100 % identité / 100 % couverture."""
    found = set()
    for cols in parse_kma_rows(res_path):
        template = cols[0]
        qid, qcov = float(cols[6]), float(cols[7])
        if qid == 100.0 and qcov == 100.0:
            for locus in LOCI:
                if template.startswith(locus + "_"):
                    found.add(locus)
                    break
    return str(len(found)) if found else ""

def read_st_species(path: Path) -> str:
    """Espèce indiquée par st_sp.tsv, quand elle existe."""
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            species = clean(row.get("species", ""))
            return species if species and species.upper() != "NA" else ""
    return ""

def index_to_letters(n: int) -> str:
    """0 -> A, 25 -> Z, 26 -> AA."""
    out = ""
    n += 1
    while n:
        n, rem = divmod(n - 1, 26)
        out = chr(65 + rem) + out
    return out

def species_from_breakdown(value: str) -> list[str]:
    """'Borrelia garinii (37); Borrelia afzelii (1)' -> espèces."""
    out, seen = [], set()
    for part in [p.strip() for p in clean(value).split(";") if p.strip()]:
        match = re.match(r"^(.*?)\s*\(\d+\)\s*$", part)
        species = (match.group(1) if match else part).strip()
        if species and species not in seen:
            seen.add(species)
            out.append(species)
    return out

def parse_specificity_file(path: Path):
    """
    Lit perfect_alleles_specificity.tsv.

    Une espèce reçoit :
      - un compteur numérique si l'allèle est spécifique ou très majoritaire ;
      - une lettre si l'allèle est ambigu entre plusieurs espèces.
    """
    tokens = defaultdict(lambda: {"specific_count": 0, "letters": set()})
    present = set()
    shared_to_letter = {}

    with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for i, row in enumerate(reader):
            locus = clean(row.get("locus", ""))
            allele = clean(row.get("allele", ""))
            major_species = clean(row.get("major_species", ""))
            total = float(clean(row.get("total", "0")))
            n_species = int(float(clean(row.get("n_species", "0"))))
            purity = float(clean(row.get("purity", "0")))
            breakdown = clean(row.get("species_breakdown", ""))

            if not major_species:
                continue

            # Allèle considéré comme indicatif de l'espèce majoritaire.
            if n_species == 1 or (total >= MIN_TOTAL_FOR_MAJOR and purity >= MIN_PURITY_FOR_MAJOR):
                tokens[major_species]["specific_count"] += 1
                present.add(major_species)
                continue

            # Allèle ambigu : même lettre associée aux espèces impliquées.
            involved = species_from_breakdown(breakdown) or [major_species]
            key = (locus, allele, tuple(involved))
            shared_to_letter.setdefault(key, index_to_letters(len(shared_to_letter)))
            letter = shared_to_letter[key]

            for species in involved:
                tokens[species]["letters"].add(letter)
                present.add(species)

    return tokens, present

def species_cell(token: dict) -> str:
    parts = []
    if token.get("specific_count", 0) > 0:
        parts.append(str(token["specific_count"]))
    parts.extend(sorted(token.get("letters", set())))
    return "; ".join(parts)

def choose_pubmlst_species(tokens: dict) -> str:
    """Choisit l'espèce PubMLST représentée, ou la co-infection si plusieurs espèces spécifiques."""
    specific = sorted(sp for sp, tok in tokens.items() if tok.get("specific_count", 0) > 0)
    letter_only = sorted(sp for sp, tok in tokens.items() if not tok.get("specific_count", 0) and tok.get("letters"))

    if len(specific) == 1:
        return specific[0]
    if len(specific) > 1:
        return "co-infection - " + ", ".join(specific)
    if letter_only:
        return ", ".join(letter_only)
    return "inconnu"

def pubmlst_only_letters(row: dict, species_cols: list[str]) -> bool:
    """True si PubMLST ne contient que des lettres ambiguës, sans compteur numérique."""
    has_letter = has_digit = False
    for col in species_cols:
        for token in [t.strip() for t in clean(row.get(col, "")).split(";") if t.strip()]:
            has_digit |= any(ch.isdigit() for ch in token)
            has_letter |= bool(re.fullmatch(r"[A-Z]+", token))
    return has_letter and not has_digit

# ---------------------------------------------------------------------------
# NCBI
# ---------------------------------------------------------------------------

def ncbi_species_from_template(template: str) -> str:
    """afzelii_GCF_...|clpA|... -> afzelii."""
    return clean(template).split("_", 1)[0]

def parse_ncbi_res(path: Path, mode: str) -> dict[str, int]:
    """Compte les espèces NCBI selon les seuils stricts ou relâchés."""
    counts = defaultdict(int)

    for cols in parse_kma_rows(path):
        qid = float(cols[6])
        qcov = float(cols[7])
        depth = float(cols[8])

        if mode == "strict":
            keep = qid == 100.0 and qcov == 100.0 and depth >= NCBI_STRICT_MIN_DEPTH
        else:
            keep = qid >= NCBI_RELAXED_MIN_QUERY_IDENTITY and depth > NCBI_RELAXED_MIN_DEPTH

        if keep:
            counts[ncbi_species_from_template(cols[0])] += 1

    return counts

def choose_ncbi_strict(counts: dict[str, int]) -> str:
    species = sorted(sp for sp, n in counts.items() if n > 0)
    if len(species) == 1:
        return species[0]
    if len(species) > 1:
        return "co-infection - " + ", ".join(species)
    return "inconnu"

def choose_ncbi_relaxed(counts: dict[str, int]) -> str:
    if not counts:
        return "inconnu"
    best = max(counts.values())
    winners = sorted(sp for sp, n in counts.items() if n == best)
    if len(winners) <= 2:
        return ", ".join(winners)
    return "inconnu"

# ---------------------------------------------------------------------------
# Choix final
# ---------------------------------------------------------------------------

def normalize_species(value: str) -> str:
    """Normalise les noms avant comparaison : Borrelia garinii -> garinii."""
    value = clean(value)
    if not value:
        return ""

    prefix = ""
    match = re.match(r"^\s*co[- ]infection\s*[-:]\s*(.+)$", value, flags=re.I)
    if match:
        prefix = "co-infection - "
        value = match.group(1)

    parts = []
    for part in [p.strip() for p in value.split(",") if p.strip()]:
        part = re.sub(r"\s*\(candidat\)\s*", "", part, flags=re.I)
        part = re.sub(r"^Borrelia\s+", "", part, flags=re.I).strip()
        part = re.sub(r"^burgdorferi\s+sensu\s+stricto$", "burgdorferi", part, flags=re.I)
        parts.append(part.lower())

    uniq = sorted(dict.fromkeys(parts))
    return prefix + ", ".join(uniq) if uniq else ""

def denormalize_species(value: str) -> str:
    """garinii -> Borrelia garinii ; burgdorferi -> Borrelia burgdorferi sensu stricto."""
    value = clean(value)
    if not value:
        return ""

    prefix = ""
    payload = value
    for marker in ("co-infection - ", "conflit - "):
        if value.lower().startswith(marker):
            prefix = marker
            payload = value[len(marker):]
            break

    names = []
    for part in [p.strip() for p in payload.split(",") if p.strip()]:
        names.append("Borrelia burgdorferi sensu stricto" if part == "burgdorferi" else f"Borrelia {part}")

    return prefix + ", ".join(names)

def choose_final(row: dict, pubmlst_species_cols: list[str]):
    p1_raw = clean(row.get("species_represented", ""))
    p2_raw = clean(row.get("species_represented_ncbi", ""))
    p3_raw = clean(row.get("species_represented_ncbi_95id_depth300", ""))

    p1 = normalize_species(p1_raw)
    p2 = normalize_species(p2_raw)
    p3 = normalize_species(p3_raw)

    if p1 and not is_unknown(p1) and not pubmlst_only_letters(row, pubmlst_species_cols):
        return denormalize_species(p1), "X", "", ""
    if p2 and not is_unknown(p2):
        return denormalize_species(p2), "", "X", ""
    if p3 and not is_unknown(p3):
        return denormalize_species(p3), "", "", "X"
    if any([p1_raw, p2_raw, p3_raw]):
        return "inconnu", "", "", ""
    return "", "", "", ""

def unique(items: list[str]) -> list[str]:
    return list(dict.fromkeys(items))

# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Construit metadata_genaspe_species_final.tsv.")
    parser.add_argument("--metadata", default="data/metadata_genaspe.csv")
    parser.add_argument("--pubmlst-dir", default="results/pubMLST_typing")
    parser.add_argument("--ncbi-dir", default="results/ncbi_typing")
    parser.add_argument("--out", default="results/metadata_genaspe_species_final.tsv")
    args = parser.parse_args()

    metadata = Path(args.metadata)
    pubmlst_dir = Path(args.pubmlst_dir)
    ncbi_dir = Path(args.ncbi_dir)
    output = Path(args.out)

    pubmlst_index = index_results_dirs(pubmlst_dir, skip={"mafft_input"})
    ncbi_index = index_results_dirs(ncbi_dir)

    with metadata.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=sniff_delimiter(metadata))
        base_fields = list(reader.fieldnames or [])
        metadata_rows = list(reader)

    rows = []
    pubmlst_species_used = set()
    ncbi_strict_used = set()
    ncbi_relaxed_used = set()

    matched_pubmlst = matched_ncbi = 0
    signal_pubmlst = signal_ncbi_strict = signal_ncbi_relaxed = 0

    # Première passe : calculs, puis collecte des colonnes espèces dynamiques.
    for base_row in metadata_rows:
        row = dict(base_row)

        pub_dir = match_sample_dir(row, pubmlst_index)
        pub_tokens = defaultdict(lambda: {"specific_count": 0, "letters": set()})
        strict_counts = defaultdict(int)
        relaxed_counts = defaultdict(int)

        if pub_dir:
            matched_pubmlst += 1
            row["nb_loci_100id_100cov"] = count_loci_100_100(pub_dir / f"{pub_dir.name}.res")
            row["species_identified"] = read_st_species(pub_dir / "st_sp.tsv")
            pub_tokens, present = parse_specificity_file(pub_dir / "perfect_alleles_specificity.tsv")
            row["species_represented"] = choose_pubmlst_species(pub_tokens)
            if present:
                signal_pubmlst += 1
                pubmlst_species_used.update(present)
        else:
            row["nb_loci_100id_100cov"] = ""
            row["species_identified"] = ""
            row["species_represented"] = ""

        ncbi_dir_sample = match_sample_dir(row, ncbi_index)
        if ncbi_dir_sample:
            matched_ncbi += 1
            res = ncbi_dir_sample / f"{ncbi_dir_sample.name}.res"
            strict_counts = parse_ncbi_res(res, "strict")
            relaxed_counts = parse_ncbi_res(res, "relaxed")

        row["species_represented_ncbi"] = choose_ncbi_strict(strict_counts) if ncbi_dir_sample else ""
        row["species_represented_ncbi_95id_depth300"] = choose_ncbi_relaxed(relaxed_counts) if ncbi_dir_sample else ""

        if strict_counts:
            signal_ncbi_strict += 1
            ncbi_strict_used.update(strict_counts)
        if relaxed_counts:
            signal_ncbi_relaxed += 1
            ncbi_relaxed_used.update(relaxed_counts)

        row["_pub_tokens"] = pub_tokens
        row["_strict_counts"] = dict(strict_counts)
        row["_relaxed_counts"] = dict(relaxed_counts)
        rows.append(row)

    pub_cols = sorted(pubmlst_species_used)
    strict_cols = [f"{sp} [NCBI]" for sp in sorted(ncbi_strict_used)]
    relaxed_cols = [f"{sp} [NCBI_95id_depth300]" for sp in sorted(ncbi_relaxed_used)]

    final_cols = [
        "nb_loci_100id_100cov", "species_identified",
        *pub_cols, "species_represented",
        *strict_cols, "species_represented_ncbi",
        *relaxed_cols, "species_represented_ncbi_95id_depth300",
        "species_final_retained",
        "used_pubMLST_identification",
        "used_NCBI-genomes_identification",
        "used_NCBI-genomes_lower-ID_identification",
    ]
    out_fields = unique(base_fields + final_cols)

    final_non_empty = final_unknown = 0

    # Deuxième passe : remplissage des colonnes dynamiques et choix final.
    for row in rows:
        pub_tokens = row.pop("_pub_tokens")
        strict_counts = row.pop("_strict_counts")
        relaxed_counts = row.pop("_relaxed_counts")

        for sp in pub_cols:
            row[sp] = species_cell(pub_tokens.get(sp, {}))

        for sp in sorted(ncbi_strict_used):
            value = strict_counts.get(sp, 0)
            row[f"{sp} [NCBI]"] = str(value) if value else ""

        for sp in sorted(ncbi_relaxed_used):
            value = relaxed_counts.get(sp, 0)
            row[f"{sp} [NCBI_95id_depth300]"] = str(value) if value else ""

        final, used_pub, used_ncbi, used_ncbi_low = choose_final(row, pub_cols)
        row["species_final_retained"] = final
        row["used_pubMLST_identification"] = used_pub
        row["used_NCBI-genomes_identification"] = used_ncbi
        row["used_NCBI-genomes_lower-ID_identification"] = used_ncbi_low

        if final:
            final_non_empty += 1
            final_unknown += int(final == "inconnu")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Métadonnées lues : {metadata}")
    print(f"Lignes metadata : {len(rows)}")
    print(f"Lignes appariées PubMLST : {matched_pubmlst}")
    print(f"Lignes avec signal espèce PubMLST : {signal_pubmlst}")
    print(f"Lignes appariées NCBI : {matched_ncbi}")
    print(f"Lignes avec signal NCBI strict : {signal_ncbi_strict}")
    print(f"Lignes avec signal NCBI relâché : {signal_ncbi_relaxed}")
    print(f"Espèces finales non vides : {final_non_empty}")
    print(f"Espèces finales inconnues : {final_unknown}")
    print(f"Fichier écrit : {output}")

if __name__ == "__main__":
    main()

#!/usr/bin/env bash
set -euo pipefail

# Environnement principal : MLST + NCBI + ospC + R.
# qc.sh appelle séparément qc_env et trimming_env avec conda run.
ENV_NAME="${ENV_NAME:-mlst_ospc_pipeline_env}"

THREADS="${THREADS:-AUTO}"
METADATA_IN="${METADATA_IN:-data/metadata_genaspe.csv}"

# Base PubMLST utilisée par typing-mapping_with_kma.py
PUBMLST_ALLELES_DIR="${PUBMLST_ALLELES_DIR:-data/pubMLST_alleles}"
PUBMLST_RAW_DIR="${PUBMLST_RAW_DIR:-${PUBMLST_ALLELES_DIR}/raw}"
PUBMLST_ALLELES_FASTA="${PUBMLST_ALLELES_FASTA:-${PUBMLST_ALLELES_DIR}/alleles.fasta}"
PUBMLST_DB_PREFIX="${PUBMLST_DB_PREFIX:-${PUBMLST_ALLELES_DIR}/borrelia_burgdorferi_sensu_lato}"

# Base NCBI utilisée par typing_kma_ncbi.py
NCBI_KMA_IDENTITY="${NCBI_KMA_IDENTITY:-95.0}"
NCBI_ROOT="${NCBI_ROOT:-data/ncbi}"
NCBI_DB_PREFIX="${NCBI_DB_PREFIX:-${NCBI_ROOT}/borrelia_ncbi_kma_index}"
NCBI_DB_PREP_SCRIPT="${NCBI_DB_PREP_SCRIPT:-db_prep.sh}"

RUN_IQTREE="${RUN_IQTREE:-1}"
RUN_OSPC="${RUN_OSPC:-1}"
OSPC_REF_FASTA="${OSPC_REF_FASTA:-data/ospC_ena_myannot.fasta}"

SCRIPTS_DIR="scripts"
PUBMLST_DIR="results/pubMLST_typing"
NCBI_DIR="results/ncbi_typing"
MAFFT_DIR="${PUBMLST_DIR}/mafft_input"
PHYLO_DIR="results/phylogeny"
FINAL_TSV="results/metadata_genaspe_species_final.tsv"
LOCI=(clpA clpX nifS pepX pyrG recG rplB uvrA)

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

for arg in "$@"; do
    case "$arg" in
        -no_opsc|-no_ospc|--no-opsc|--no-ospc|--no_opsc|--no_ospc)
            RUN_OSPC=0
            ;;
    esac
done

pubmlst_db_exists() {
    [[ -s "${PUBMLST_DB_PREFIX}.name" && \
       -s "${PUBMLST_DB_PREFIX}.comp.b" && \
       -s "${PUBMLST_DB_PREFIX}.length.b" && \
       -s "${PUBMLST_DB_PREFIX}.seq.b" ]]
}

prepare_pubmlst_db_if_needed() {
    if pubmlst_db_exists; then
        log "Base KMA PubMLST trouvé : ${PUBMLST_DB_PREFIX}"
    else
        log "Base KMA PubMLST absente : construction de ${PUBMLST_ALLELES_FASTA} puis indexation"
        mkdir -p "${PUBMLST_ALLELES_DIR}"
        : > "${PUBMLST_ALLELES_FASTA}"

        for locus in "${LOCI[@]}"; do
            cat "${PUBMLST_RAW_DIR}/${locus}.fas" >> "${PUBMLST_ALLELES_FASTA}"
            printf '\n' >> "${PUBMLST_ALLELES_FASTA}"
        done

        kma index -i "${PUBMLST_ALLELES_FASTA}" -o "${PUBMLST_DB_PREFIX}"
    fi
}

ncbi_db_exists() {
    [[ -s "${NCBI_DB_PREFIX}.name" && \
       -s "${NCBI_DB_PREFIX}.comp.b" && \
       -s "${NCBI_DB_PREFIX}.length.b" && \
       -s "${NCBI_DB_PREFIX}.seq.b" ]]
}

prepare_ncbi_db_if_needed() {
    if ncbi_db_exists; then
        log "Base KMA NCBI déjà présente : ${NCBI_DB_PREFIX}"
    else
        log "Base KMA NCBI absente : extraction des séquences et indexation"
        (cd "${NCBI_ROOT}/Borrelia" && bash extract_borrelia.sh)
        (cd "${NCBI_ROOT}/Borreliella" && bash extract_borreliella.sh)
        (cd "${NCBI_ROOT}" && bash "${NCBI_DB_PREP_SCRIPT}")
    fi
}

run_r() {
    Rscript "${SCRIPTS_DIR}/$1"
}

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

mkdir -p "${PHYLO_DIR}"

log "1/11 Préparation de la base KMA PubMLST si nécessaire"
prepare_pubmlst_db_if_needed

log "2/11 QC pré-trim, trimming et QC post-trim"
bash "${SCRIPTS_DIR}/qc.sh"

log "3/11 Spécificité des allèles PubMLST"
python3 "${SCRIPTS_DIR}/loci_specificity.py" --quiet

log "4/11 Typage PubMLST par KMA + extraction des FASTA par locus"
python3 "${SCRIPTS_DIR}/typing-mapping_with_kma.py"

log "5/11 Préparation de la base KMA NCBI si nécessaire"
prepare_ncbi_db_if_needed

log "6/11 Typage complémentaire NCBI par KMA"
python3 "${SCRIPTS_DIR}/typing_kma_ncbi.py" \
    --db "${NCBI_DB_PREFIX}" \
    --identity "${NCBI_KMA_IDENTITY}"

log "7/11 Construction du tableau final des espèces"
python3 "${SCRIPTS_DIR}/build_final_species_table.py" \
    --metadata "${METADATA_IN}" \
    --pubmlst-dir "${PUBMLST_DIR}" \
    --ncbi-dir "${NCBI_DIR}" \
    --out "${FINAL_TSV}"

python3 "${SCRIPTS_DIR}/make_out_table.py"

log "8/11 Alignements MAFFT MLST locus par locus"
for locus in "${LOCI[@]}"; do
    mafft --maxiterate 1000 --localpair \
        "${MAFFT_DIR}/${locus}.fasta" \
        > "${MAFFT_DIR}/${locus}.aln.fasta"
done

log "9/11 Concaténation des 8 loci MLST"
python3 "${SCRIPTS_DIR}/concatenate_loci.py" \
    "${MAFFT_DIR}/clpA.aln.fasta" \
    "${MAFFT_DIR}/clpX.aln.fasta" \
    "${MAFFT_DIR}/nifS.aln.fasta" \
    "${MAFFT_DIR}/pepX.aln.fasta" \
    "${MAFFT_DIR}/pyrG.aln.fasta" \
    "${MAFFT_DIR}/recG.aln.fasta" \
    "${MAFFT_DIR}/rplB.aln.fasta" \
    "${MAFFT_DIR}/uvrA.aln.fasta" \
    -o "${PHYLO_DIR}/MLST_concat.fasta" \
    -p "${PHYLO_DIR}/MLST_partitions.txt" \
    --missing N

if [[ "${RUN_IQTREE}" == "1" ]]; then
    log "10/11 IQ-TREE sur l'alignement MLST concaténé"
    iqtree3 \
        -st DNA \
        -s "${PHYLO_DIR}/MLST_concat.fasta" \
        -p "${PHYLO_DIR}/MLST_partitions.txt" \
        -m TESTMERGE \
        -bb 1000 \
        -nt "${THREADS}" \
        -pre "${PHYLO_DIR}/MLST_iqtree" \
        -redo
else
    log "10/11 IQ-TREE MLST ignoré"
fi

if [[ "${RUN_OSPC}" == "1" ]]; then
    log "11/11 ospC : construction des candidats"
    bash "${SCRIPTS_DIR}/ospc_build_candidates.sh" "${OSPC_REF_FASTA}"

    log "ospC : sélection, alignement et arbre"
    bash "${SCRIPTS_DIR}/run_ospc_selection_and_tree.sh" "${OSPC_REF_FASTA}"
else
    log "11/11 ospC ignoré"
fi

log "Analyses R : statistiques globales"
run_r "analyse_stats_genASPE.R"

log "Analyses R : post-hoc"
run_r "analyse_stats_genASPE_posthoc.R"

log "Figures R : répartition 2023/2024"
run_r "repartition_borrelia_year.R"

log "Figures R : répartition altitude"
run_r "repartition_borrelia_altitude.R"

log "Figures R : arbre MLST altitude/année"
run_r "tree_MLST_heatmap_alt_year.R"

log "Figures R : arbre MLST présence/absence des loci"
run_r "tree_MLST_loci.R"

if [[ "${RUN_OSPC}" == "1" ]]; then
    log "Analyses R ospC : iNEXT"
    run_r "ospc_iNEXT.R"

    log "Figures R ospC : arbre ospC"
    run_r "tree_ospC.R"

    log "Figures R ospC : MLST / ospC face à face"
    run_r "face_to_face_MLST_ospC.R"
else
    log "Analyses et figures R ospC ignorées (-no_ospc)"
fi

log "Terminé. Fichier final : ${FINAL_TSV}"

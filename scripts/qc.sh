#!/usr/bin/env bash
set -eo pipefail

# QC et trimming.
# FastQC/MultiQC et Cutadapt sont volontairement lancés dans deux environnements
# séparés, car ces outils peuvent être incompatibles dans un même env conda.

QC_ENV="${QC_ENV:-qc_env}"
TRIM_ENV="${TRIM_ENV:-trimming_env}"

INPUT_DIR="${INPUT_DIR:-data/fastq_raw}"
TRIMMED_DIR="${TRIMMED_DIR:-data/fastq}"
QC_ROOT="${QC_ROOT:-results/qc}"
QC_THREADS="${QC_THREADS:-4}"

PRE_TRIM_QC="${PRE_TRIM_QC:-${QC_ROOT}/1_pre_trim}"
POST_TRIM_QC="${POST_TRIM_QC:-${QC_ROOT}/2_post_trim}"

mkdir -p "${PRE_TRIM_QC}" "${TRIMMED_DIR}" "${POST_TRIM_QC}"

echo "-> QC pré-trim avec ${QC_ENV}"
conda run -n "${QC_ENV}" fastqc "${INPUT_DIR}"/*.fastq.gz \
    -t "${QC_THREADS}" \
    -o "${PRE_TRIM_QC}"

conda run -n "${QC_ENV}" multiqc "${PRE_TRIM_QC}" \
    -o "${PRE_TRIM_QC}" \
    -n "multiqc_pre_trim"

echo "-> Trimming avec ${TRIM_ENV}"
for R1 in "${INPUT_DIR}"/*_R1.fastq.gz; do
    R2="${R1/_R1/_R2}"
    SAMPLE="$(basename "${R1}" _R1.fastq.gz)"

    echo "trim pour ${SAMPLE}"
    conda run -n "${TRIM_ENV}" cutadapt \
        -q 20,20 \
        --minimum-length 50 \
        --cores "${QC_THREADS}" \
        -o "${TRIMMED_DIR}/${SAMPLE}_trimmed_R1.fastq.gz" \
        -p "${TRIMMED_DIR}/${SAMPLE}_trimmed_R2.fastq.gz" \
        "${R1}" "${R2}" \
        > "${TRIMMED_DIR}/${SAMPLE}_report.txt"
done

echo "-> QC post-trim avec ${QC_ENV}"
conda run -n "${QC_ENV}" fastqc "${TRIMMED_DIR}"/*.fastq.gz \
    -t "${QC_THREADS}" \
    -o "${POST_TRIM_QC}"

conda run -n "${QC_ENV}" multiqc "${POST_TRIM_QC}" \
    -o "${POST_TRIM_QC}" \
    -n "multiqc_post_trim"

echo "-> QC/trimming terminé"

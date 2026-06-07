#!/usr/bin/env bash
set -euo pipefail

# Configuration
CONDA_ENV="mlst_ospc_pipeline_env"
PRIMER_FILE="primer_mya"
DATA_DIR="data"
MISMATCHES=4

# conda
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

GENE_TABLE="$TMPDIR/gene_primers.tsv"
GENOME_LIST="$TMPDIR/genomes.list"
GLOBAL_SUMMARY="summary_missing_by_gene.tsv"

# Préparation des primers
awk '
BEGIN { OFS="\t" }
NF >= 3 {
    print $1, $2, $3
}
' "$PRIMER_FILE" > "$GENE_TABLE"

find "$DATA_DIR" -type f -name "*.fna" | sort > "$GENOME_LIST"

mapfile -t GENES < <(cut -f1 "$GENE_TABLE")

for gene in "${GENES[@]}"; do
    : > "${gene}.fa"
    : > "${gene}_missing.txt"
done

: > "$GLOBAL_SUMMARY"
echo -e "gene\tgenomes_total\tgenomes_with_sequence\tgenomes_missing" >> "$GLOBAL_SUMMARY"

total_genomes=$(wc -l < "$GENOME_LIST")

echo "Nombre de génomes détectés : $total_genomes"
echo "Gènes détectés : ${#GENES[@]}"
echo

# Extraction
while IFS= read -r genome; do
    genome_dir=$(basename "$(dirname "$genome")")
    genome_file=$(basename "$genome")

    echo "Traitement : $genome_dir / $genome_file"

    while IFS=$'\t' read -r gene fwd rev; do
        tmp_amp="$TMPDIR/${gene}_amplicon.fa"
        : > "$tmp_amp"

        seqkit amplicon \
            -F "$fwd" \
            -R "$rev" \
            -m "$MISMATCHES" \
            -r 51:-51 \
            "$genome" > "$tmp_amp" 2>/dev/null || true

        if [[ -s "$tmp_amp" ]]; then
            awk -v gid="$genome_dir" -v gene="$gene" '
            BEGIN { n=0 }
            /^>/ {
                n++
                print ">" gid "|" gene "|" n
                next
            }
            { print }
            ' "$tmp_amp" >> "${gene}.fa"
        else
            echo "$genome_dir" >> "${gene}_missing.txt"
        fi
    done < "$GENE_TABLE"
done < "$GENOME_LIST"

echo
echo "Récapitulatif :"
for gene in "${GENES[@]}"; do
    missing_count=0
    with_count=0

    if [[ -s "${gene}_missing.txt" ]]; then
        missing_count=$(wc -l < "${gene}_missing.txt")
    fi

    with_count=$(( total_genomes - missing_count ))

    echo -e "${gene}\t${total_genomes}\t${with_count}\t${missing_count}" >> "$GLOBAL_SUMMARY"
    echo "  - $gene : $with_count / $total_genomes génomes avec séquence ; $missing_count manquants"
done

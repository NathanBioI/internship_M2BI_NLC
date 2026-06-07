#!/usr/bin/env bash
set -eo pipefail
log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}
REF_FASTA_RAW="${1:-data/ospC_ena_myannot.fasta}"
FASTQ_DIR="${FASTQ_DIR:-data/fastq}"
OUT_ROOT="${OUT_ROOT:-results/ospc}"
THREADS="${THREADS:-4}"
KMA_ID="${KMA_ID:-0.70}"
REF_DIR="${REF_DIR:-data/ospc_pipeline_refs}"
REF_MIN_LEN="${REF_MIN_LEN:-450}"
REF_MAX_LEN="${REF_MAX_LEN:-1100}"
KMA_REF_FASTA="${KMA_REF_FASTA:-$REF_DIR/ospC_refs_for_kma.fasta}"
HYBPIPER_REF_FASTA="${HYBPIPER_REF_FASTA:-$REF_DIR/ospC_refs_for_hybpiper.fasta}"
REF_MAP_TSV="${REF_MAP_TSV:-$REF_DIR/ospC_reference_header_map.tsv}"
BASE_KMA_DIR="${BASE_KMA_DIR:-$OUT_ROOT/ospc_kma}"
BASE_KMA_DB_PREFIX="${BASE_KMA_DB_PREFIX:-$BASE_KMA_DIR/db/ospC_refs_kma_db}"
SEED_BUILD_DIR="${SEED_BUILD_DIR:-$OUT_ROOT/ospc_seed_recruitment}"
SEEDS_FASTA="${SEEDS_FASTA:-$SEED_BUILD_DIR/seeds/ospC_conserved_seeds.fasta}"
SEED_SPADES_DIR="${SEED_SPADES_DIR:-$OUT_ROOT/ospc_seed_spades}"
SEED_MIN_REF_LEN="${SEED_MIN_REF_LEN:-500}"
SEED_WINDOW="${SEED_WINDOW:-41}"
SEED_STEP="${SEED_STEP:-5}"
SEED_MIN_IDENT="${SEED_MIN_IDENT:-0.85}"
SEED_MAX_GAP_FRAC="${SEED_MAX_GAP_FRAC:-0.20}"
SEED_MIN_SEEDS="${SEED_MIN_SEEDS:-1}"
BBDUK_SEED_K="${BBDUK_SEED_K:-19}"
BBDUK_SEED_HDIST="${BBDUK_SEED_HDIST:-1}"
GUIDED_SPADES_DIR="${GUIDED_SPADES_DIR:-$OUT_ROOT/ospc_kma_guided_spades}"
GUIDED_CONS_MIN_LEN="${GUIDED_CONS_MIN_LEN:-300}"
GUIDED_CONS_MAX_N_FRAC="${GUIDED_CONS_MAX_N_FRAC:-0.02}"
BBDUK_GUIDED_K="${BBDUK_GUIDED_K:-19}"
BBDUK_GUIDED_HDIST="${BBDUK_GUIDED_HDIST:-1}"
SPADES_THREADS="${SPADES_THREADS:-$THREADS}"
SPADES_MEM="${SPADES_MEM:-12}"
HYBPIPER_DIR="${HYBPIPER_DIR:-$OUT_ROOT/hybpiper_ospc}"
FINAL_KMA_DIR="${FINAL_KMA_DIR:-$OUT_ROOT/ospc_kma_final}"
FINAL_KMA_DB_PREFIX_SUFFIX="${FINAL_KMA_DB_PREFIX_SUFFIX:-kma_db}"
REF_QI_MIN="${REF_QI_MIN:-60}"
REF_QC_MIN="${REF_QC_MIN:-95}"
REF_DEPTH_MIN="${REF_DEPTH_MIN:-1}"
INCLUDE_BEST_REF_IF_NONE="${INCLUDE_BEST_REF_IF_NONE:-1}"
CONTIG_MIN_LEN="${CONTIG_MIN_LEN:-1}"
mkdir -p "$OUT_ROOT"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPER_PY="${HELPER_PY:-$SCRIPT_DIR/ospc_pipeline_helpers.py}"
find_r1_files() {
  find "$FASTQ_DIR" -maxdepth 1 -type f -name "*_R1.fastq.gz" | sort
}
r2_from_r1() {
  local r1="$1"
  printf '%s\n' "${r1/_R1.fastq.gz/_R2.fastq.gz}"
}
sample_from_r1() {
  local r1="$1"
  basename "$r1" _R1.fastq.gz
}
run_recruit_and_spades() {
  local sample="$1"
  local r1="$2"
  local r2="$3"
  local bait_fasta="$4"
  local out_dir="$5"
  local k="$6"
  local hdist="$7"
  local tmp="$out_dir/tmp"
  local spades_out="$out_dir/spades"
  mkdir -p "$tmp"
  local repaired="$tmp/${sample}.repaired_interleaved.fastq.gz"
  local recruited="$tmp/${sample}.recruited_interleaved.fastq.gz"
  local recruited_r1="$tmp/${sample}.recruited_R1.fastq.gz"
  local recruited_r2="$tmp/${sample}.recruited_R2.fastq.gz"
  local singletons="$tmp/${sample}.singletons.fastq.gz"
  if repair.sh \
      in1="$r1" \
      in2="$r2" \
      out="$repaired" \
      outs="$singletons" \
      repair \
      > "$tmp/repair.stdout.log" \
      2> "$tmp/repair.stderr.log"; then
    :
  else
    log "Sample $sample : repair.sh failed"
    echo "repair_failed" > "$out_dir/status.txt"
    return 0
  fi
  if bbduk.sh \
      in="$repaired" \
      outm="$recruited" \
      ref="$bait_fasta" \
      k="$k" \
      hdist="$hdist" \
      threads=1 \
      interleaved=t \
      stats="$tmp/bbduk.stats.txt" \
      > "$tmp/bbduk.stdout.log" \
      2> "$tmp/bbduk.stderr.log"; then
    :
  else
    log "Sample $sample : BBDuk failed"
    echo "bbduk_failed" > "$out_dir/status.txt"
    return 0
  fi
  if reformat.sh \
      in="$recruited" \
      out1="$recruited_r1" \
      out2="$recruited_r2" \
      overwrite=true \
      > "$tmp/reformat.stdout.log" \
      2> "$tmp/reformat.stderr.log"; then
    :
  else
    log "Sample $sample : reformat.sh failed"
    echo "reformat_failed" > "$out_dir/status.txt"
    return 0
  fi
  local pairs=0
  if [ -s "$recruited_r1" ]; then
    pairs=$(zcat -f "$recruited_r1" 2>/dev/null | awk 'END{print int(NR/4)}')
  fi
  if [ "$pairs" -gt 0 ]; then
    rm -rf "$spades_out"
    if spades.py \
        -1 "$recruited_r1" \
        -2 "$recruited_r2" \
        -o "$spades_out" \
        --careful \
        -k 21,33,55 \
        -t "$SPADES_THREADS" \
        -m "$SPADES_MEM" \
        > "$tmp/spades.stdout.log" \
        2> "$tmp/spades.stderr.log"; then
      log "Sample $sample : SPAdes ok avec $pairs paires recrutées"
      echo "ok_recruited_pairs=$pairs" > "$out_dir/status.txt"
    else
      log "Sample $sample : SPAdes failed avec $pairs paires recrutées"
      echo "failed_recruited_pairs=$pairs" > "$out_dir/status.txt"
    fi
  else
    log "Sample $sample : zero recruited read pairs"
    echo "no_seed_recruited_reads=0" > "$out_dir/status.txt"
  fi
  return 0
}
log "[0/5] Préparation des références ospC"
mkdir -p "$REF_DIR"
    python3 "$HELPER_PY" prepare_refs \
      --input "$REF_FASTA_RAW" \
      --outdir "$REF_DIR" \
      --kma-fasta "$KMA_REF_FASTA" \
      --hybpiper-fasta "$HYBPIPER_REF_FASTA" \
      --map-tsv "$REF_MAP_TSV" \
      --min-len "$REF_MIN_LEN" \
      --max-len "$REF_MAX_LEN"
log "[1/5] KMA initial sur les références ospC"
mkdir -p "$BASE_KMA_DIR/db"
kma index -i "$KMA_REF_FASTA" -o "$BASE_KMA_DB_PREFIX"
  while IFS= read -r r1; do
    sample="$(sample_from_r1 "$r1")"
    r2="$(r2_from_r1 "$r1")"
    sample_out="$BASE_KMA_DIR/$sample"
    mkdir -p "$sample_out"
    log "Sample $sample : KMA initial"
    kma \
      -ipe "$r1" "$r2" \
      -o "$sample_out/$sample" \
      -t_db "$BASE_KMA_DB_PREFIX" \
      -1t1 \
      -ID "$KMA_ID" \
      -and \
      -ref_fsa \
      -oa \
      -ef \
      -t "$THREADS"
done < <(find_r1_files)
log "[2/5] Méthode seed + SPAdes"
  mkdir -p "$SEED_BUILD_DIR/work" "$SEED_BUILD_DIR/seeds" "$SEED_SPADES_DIR"
  FILTERED_REFS="$SEED_BUILD_DIR/work/ospC_refs.min${SEED_MIN_REF_LEN}.fasta"
  ALN_FASTA="$SEED_BUILD_DIR/work/ospC_refs.min${SEED_MIN_REF_LEN}.mafft.fasta"
  SEED_REPORT="$SEED_BUILD_DIR/seeds/ospC_conserved_seeds.tsv"
log "Construction des graines conservées ospC"
    python3 "$HELPER_PY" filter_refs \
      --input "$KMA_REF_FASTA" \
      --output "$FILTERED_REFS" \
      --min-len "$SEED_MIN_REF_LEN"
    mafft --auto --thread "$THREADS" "$FILTERED_REFS" > "$ALN_FASTA" 2> "$SEED_BUILD_DIR/work/mafft.log"
    python3 "$HELPER_PY" build_seeds \
      --alignment "$ALN_FASTA" \
      --output "$SEEDS_FASTA" \
      --report "$SEED_REPORT" \
      --window "$SEED_WINDOW" \
      --step "$SEED_STEP" \
      --min-ident "$SEED_MIN_IDENT" \
      --max-gap-frac "$SEED_MAX_GAP_FRAC" \
      --min-seeds "$SEED_MIN_SEEDS"
while IFS= read -r r1; do
    sample="$(sample_from_r1 "$r1")"
    r2="$(r2_from_r1 "$r1")"
    sample_out="$SEED_SPADES_DIR/$sample"
    mkdir -p "$sample_out"
    log "Sample $sample : seed recruitment + SPAdes"
    run_recruit_and_spades \
      "$sample" "$r1" "$r2" "$SEEDS_FASTA" "$sample_out" \
      "$BBDUK_SEED_K" "$BBDUK_SEED_HDIST"
done < <(find_r1_files)
log "[3/5] Méthode KMA-guided + SPAdes"
mkdir -p "$GUIDED_SPADES_DIR"
  find "$BASE_KMA_DIR" -mindepth 2 -maxdepth 2 -name "*.res" | sort | while read -r res; do
    sample="$(basename "$res" .res)"
    r1="$FASTQ_DIR/${sample}_R1.fastq.gz"
    r2="$FASTQ_DIR/${sample}_R2.fastq.gz"
    sample_out="$GUIDED_SPADES_DIR/$sample"
    tmp="$sample_out/tmp"
    mkdir -p "$tmp"
    bait="$tmp/${sample}.bait.fasta"
    fsa="$(dirname "$res")/${sample}.fsa"
    log "Sample $sample : choix du bait guidé par KMA"
      python3 "$HELPER_PY" choose_bait \
        --res "$res" \
        --fsa "$fsa" \
        --refs "$KMA_REF_FASTA" \
        --output "$bait" \
        --cons-min-len "$GUIDED_CONS_MIN_LEN" \
        --cons-max-n-frac "$GUIDED_CONS_MAX_N_FRAC" \
        > "$tmp/bait_choice.tsv"
    log "Sample $sample : guided recruitment + SPAdes"
    run_recruit_and_spades \
      "$sample" "$r1" "$r2" "$bait" "$sample_out" \
      "$BBDUK_GUIDED_K" "$BBDUK_GUIDED_HDIST"
  done
log "[4/5] HybPiper"
  mkdir -p "$HYBPIPER_DIR"
  ABS_HYB_REF="$(realpath "$HYBPIPER_REF_FASTA")"
  while IFS= read -r r1; do
    sample="$(sample_from_r1 "$r1")"
    r2="$(r2_from_r1 "$r1")"
    r1_abs="$(realpath "$r1")"
    r2_abs="$(realpath "$r2")"
    log "Sample $sample : HybPiper assemble"
    (
      cd "$HYBPIPER_DIR"
      hybpiper assemble \
        -r "$r1_abs" "$r2_abs" \
        -t_dna "$ABS_HYB_REF" \
        --prefix "$sample" \
        --bwa \
        --cpu "$THREADS"
    )
  done < <(find_r1_files)
  mapfile -t hyb_samples < <(find "$HYBPIPER_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
  if [ "${#hyb_samples[@]}" -gt 0 ]; then
    sample_names="$(printf '%s,' "${hyb_samples[@]}" | sed 's/,$//')"
    log "HybPiper retrieve_sequences"
    (
      cd "$HYBPIPER_DIR"
      hybpiper retrieve_sequences dna \
        -t_dna "$ABS_HYB_REF" \
        --sample_names "$sample_names" \
        --fasta_dir "final_fasta_for_tree"
    ) || log "retrieve_sequences a échoué ou n'a rien récupéré ; les contigs par échantillon restent utilisables"
fi
log "[5/5] KMA final personnalisé par échantillon"
mkdir -p "$FINAL_KMA_DIR"
  SUMMARY="$FINAL_KMA_DIR/summary.tsv"
  printf "sample\tn_refs_selected\tn_seed_contigs\tn_guided_contigs\tn_hybpiper_contigs\tdb_sequences\tkma_status\n" > "$SUMMARY"
  find "$BASE_KMA_DIR" -mindepth 2 -maxdepth 2 -name "*.res" | sort | while read -r base_res; do
    sample="$(basename "$base_res" .res)"
    r1="$FASTQ_DIR/${sample}_R1.fastq.gz"
    r2="$FASTQ_DIR/${sample}_R2.fastq.gz"
    sample_out="$FINAL_KMA_DIR/$sample"
    mkdir -p "$sample_out"
    db_fasta="$sample_out/${sample}.personalized_db.fasta"
    db_prefix="$sample_out/${sample}_${FINAL_KMA_DB_PREFIX_SUFFIX}"
    seed_contigs="$SEED_SPADES_DIR/$sample/spades/contigs.fasta"
    guided_contigs="$GUIDED_SPADES_DIR/$sample/spades/contigs.fasta"
    hybpiper_contigs="$HYBPIPER_DIR/$sample/ospC/ospC_contigs.fasta"
    log "Sample $sample : construction DB finale"
    db_stats="$(python3 "$HELPER_PY" make_final_db \
      --sample "$sample" \
      --base-res "$base_res" \
      --refs "$KMA_REF_FASTA" \
      --seed-contigs "$seed_contigs" \
      --guided-contigs "$guided_contigs" \
      --hybpiper-contigs "$hybpiper_contigs" \
      --output "$db_fasta" \
      --ref-qi-min "$REF_QI_MIN" \
      --ref-qc-min "$REF_QC_MIN" \
      --ref-depth-min "$REF_DEPTH_MIN" \
      --include-best-ref-if-none "$INCLUDE_BEST_REF_IF_NONE" \
      --contig-min-len "$CONTIG_MIN_LEN")"
    nseq="$(printf '%s\n' "$db_stats" | awk -F'\t' '{print $6}')"

    # Si aucune référence/contig n'a été récupéré pour cet échantillon,
    # la base personnalisée est vide. On ignore seulement le KMA final
    # de cet échantillon et le pipeline continue.
    if [ ! -s "$db_fasta" ] || [ "${nseq:-0}" = "0" ]; then
      log "Sample $sample : DB finale vide, KMA final ignoré"
      printf "%s	db_empty
" "$db_stats" >> "$SUMMARY"
      continue
    fi

    kma index -i "$db_fasta" -o "$db_prefix"
      log "Sample $sample : KMA final"
      kma \
        -ipe "$r1" "$r2" \
        -o "$sample_out/$sample" \
        -t_db "$db_prefix" \
        -1t1 \
        -ID "$KMA_ID" \
        -and \
        -ref_fsa \
        -oa \
        -ef \
        -t "$THREADS"
      printf "%s\tok\n" "$db_stats" >> "$SUMMARY"
done
log "Pipeline candidats terminé."
log "Références KMA : $KMA_REF_FASTA"
log "Références HybPiper : $HYBPIPER_REF_FASTA"
log "KMA initial : $BASE_KMA_DIR"
log "Seed SPAdes : $SEED_SPADES_DIR"
log "Guided SPAdes : $GUIDED_SPADES_DIR"
log "HybPiper : $HYBPIPER_DIR"
log "KMA final personnalisé : $FINAL_KMA_DIR"
log "Pour la sélection/phylogénie, utiliser PERSONALIZED_DIR=$FINAL_KMA_DIR"

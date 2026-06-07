#!/usr/bin/env bash
set -eo pipefail
log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}
REF_FASTA="${1:-data/ospC_ena_myannot.fasta}"
PERSONALIZED_DIR="${PERSONALIZED_DIR:-results/ospc/ospc_kma_final}"
OUT_DIR="${OUT_DIR:-results/ospc/ospc_selection}"
THREADS="${THREADS:-4}"
CONSENSUS_FROM="${CONSENSUS_FROM:-references}"
CONSENSUS_MIN_LEN="${CONSENSUS_MIN_LEN:-350}"
CONSENSUS_MAX_N_FRAC="${CONSENSUS_MAX_N_FRAC:-0.01}"
KMA_CONSENSUS_POLICY="${KMA_CONSENSUS_POLICY:-rescue_only}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SELECTOR="${PY_SELECTOR:-$SCRIPT_DIR/select_ospc_for_phylogeny.py}"
mkdir -p "$OUT_DIR"
log "Sélection des séquences ospC"
python3 "$PY_SELECTOR" \
  --personalized-dir "$PERSONALIZED_DIR" \
  --out-dir "$OUT_DIR" \
  --ref-fasta "$REF_FASTA" \
  --threads "$THREADS" \
  --consensus-from "$CONSENSUS_FROM" \
  --consensus-min-len "$CONSENSUS_MIN_LEN" \
  --consensus-max-n-frac "$CONSENSUS_MAX_N_FRAC" \
  --kma-consensus-policy "$KMA_CONSENSUS_POLICY"
SELECTED_CONSENSUS="$OUT_DIR/selected_ospC_consensus_oriented.fasta"
ALN="$OUT_DIR/selected_ospC_consensus_oriented.mafft.fasta"
log "Alignement MAFFT"
mafft --auto --thread "$THREADS" "$SELECTED_CONSENSUS" > "$ALN" 2> "$OUT_DIR/mafft.log"
log "Arbre IQ-TREE"
iqtree3 \
  -redo \
  -s "$ALN" \
  -m MFP \
  -bb 1000 \
  -alrt 1000 \
  -nt AUTO \
  -pre "$OUT_DIR/selected_ospC_iqtree" \
  > "$OUT_DIR/iqtree.stdout.log" \
  2> "$OUT_DIR/iqtree.stderr.log"
log "Terminé. Sorties principales :"
log "  $OUT_DIR/ospc_sample_decisions.tsv"
log "  $OUT_DIR/ospc_candidate_summary.tsv"
log "  $OUT_DIR/selected_ospC_consensus_oriented.fasta"
log "  $OUT_DIR/selected_ospC_iqtree.treefile"

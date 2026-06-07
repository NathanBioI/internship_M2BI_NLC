#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

def log(msg: str) -> None:
    print(f'[select_ospc] {msg}', file=sys.stderr)

def die(msg: str, code: int=1) -> None:
    print(f'[select_ospc] ERROR: {msg}', file=sys.stderr)
    raise SystemExit(code)

def sanitize_id(text: str) -> str:
    text = text.strip()
    text = re.sub('\\s+', '_', text)
    text = re.sub('[^A-Za-z0-9_.|:+\\-]+', '_', text)
    return text[:240]

def read_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    name: Optional[str] = None
    seq_parts: List[str] = []
    with path.open('r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield (name, ''.join(seq_parts).upper())
                name = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if name is not None:
            yield (name, ''.join(seq_parts).upper())

def write_fasta(records: Iterable[Tuple[str, str]], path: Path) -> None:
    with path.open('w', encoding='utf-8') as out:
        for name, seq in records:
            out.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + '\n')
KMA_DEFAULT_COLUMNS = ['Template', 'Score', 'Expected', 'Template_length', 'Template_Identity', 'Template_Coverage', 'Query_Identity', 'Query_Coverage', 'Depth', 'q_value', 'p_value']

def to_float(x: str) -> float:
    try:
        if x is None or x == '' or x.upper() == 'NA':
            return float('nan')
        return float(str(x).replace(',', '.'))
    except Exception:
        return float('nan')

def parse_kma_res(path: Path) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    header: Optional[List[str]] = None
    with path.open('r', encoding='utf-8', errors='replace') as handle:
        for raw in handle:
            line = raw.rstrip('\n')
            if not line:
                continue
            if line.startswith('#'):
                clean = line.lstrip('#').strip()
                if clean.lower().startswith('template'):
                    header = clean.split('\t')
                continue
            parts = line.split('\t')
            cols = header if header and len(header) <= len(parts) else KMA_DEFAULT_COLUMNS
            data = {cols[i]: parts[i] for i in range(min(len(cols), len(parts)))}
            template = str(data.get('Template', data.get('#Template', parts[0]))).strip()
            if not template:
                continue
            data['template'] = template
            data['score'] = to_float(str(data.get('Score', 'nan')))
            data['template_length'] = to_float(str(data.get('Template_length', 'nan')))
            data['template_identity'] = to_float(str(data.get('Template_Identity', 'nan')))
            data['template_coverage'] = to_float(str(data.get('Template_Coverage', 'nan')))
            data['query_identity'] = to_float(str(data.get('Query_Identity', 'nan')))
            data['query_coverage'] = to_float(str(data.get('Query_Coverage', 'nan')))
            data['depth'] = to_float(str(data.get('Depth', 'nan')))
            rows[template] = data
    return rows
CODON_TABLE = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
RC_TABLE = str.maketrans('ACGTRYSWKMBDHVNacgtryswkmbdhvn', 'TGCAYRSWMKVHDBNtgcayrswmkvhdbn')

def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()

@dataclass
class OrfInfo:
    aa_len: int = 0
    nt_start: int = 0
    nt_end: int = 0
    strand: str = '+'
    frame: int = 0
    oriented_seq: str = ''
    orf_nt: str = ''
    aa_seq: str = ''

def longest_stop_free_orf(seq: str) -> OrfInfo:
    seq = seq.upper().replace('-', '')
    best = OrfInfo(oriented_seq=seq)
    for strand, oriented in [('+', seq), ('-', revcomp(seq))]:
        for frame in (0, 1, 2):
            current_start = frame
            current_aas: List[str] = []
            for pos in range(frame, len(oriented) - 2, 3):
                codon = oriented[pos:pos + 3]
                aa = CODON_TABLE.get(codon, 'X')
                if aa == '*':
                    if len(current_aas) > best.aa_len:
                        best = OrfInfo(aa_len=len(current_aas), nt_start=current_start, nt_end=pos, strand=strand, frame=frame, oriented_seq=oriented, orf_nt=oriented[current_start:pos], aa_seq=''.join(current_aas))
                    current_start = pos + 3
                    current_aas = []
                else:
                    current_aas.append(aa)
            end = len(oriented) - (len(oriented) - frame) % 3
            if len(current_aas) > best.aa_len:
                best = OrfInfo(aa_len=len(current_aas), nt_start=current_start, nt_end=end, strand=strand, frame=frame, oriented_seq=oriented, orf_nt=oriented[current_start:end], aa_seq=''.join(current_aas))
    return best

@dataclass
class BlastHit:
    sseqid: str = ''
    pident: float = float('nan')
    length: float = float('nan')
    qlen: float = float('nan')
    slen: float = float('nan')
    qcovs: float = float('nan')
    evalue: str = ''
    bitscore: float = float('nan')
    stitle: str = ''
    species_hint: str = ''

def parse_species_hint(text: str) -> str:
    low = text.lower().replace('_', ' ')
    species_patterns = [('Borrelia afzelii', ['borrelia afzelii', 'borreliella afzelii', ' afzelii']), ('Borrelia garinii', ['borrelia garinii', 'borreliella garinii', ' garinii']), ('Borrelia burgdorferi sensu stricto', ['burgdorferi sensu stricto', 'borrelia burgdorferi', 'borreliella burgdorferi', ' burgdorferi']), ('Borrelia valaisiana', ['borrelia valaisiana', 'borreliella valaisiana', ' valaisiana']), ('Borrelia lusitaniae', ['borrelia lusitaniae', 'borreliella lusitaniae', ' lusitaniae']), ('Borrelia bavariensis', ['borrelia bavariensis', 'borreliella bavariensis', ' bavariensis']), ('Borrelia spielmanii', ['borrelia spielmanii', 'borreliella spielmanii', ' spielmanii']), ('Borrelia bissettiae', ['borrelia bissettiae', 'borrelia bissettii', ' borreliella bissettiae', ' bissettiae', ' bissettii']), ('Borrelia turdi', ['borrelia turdi', 'borreliella turdi', ' turdi']), ('Borrelia miyamotoi', ['borrelia miyamotoi', ' miyamotoi'])]
    padded = ' ' + low + ' '
    for canonical, pats in species_patterns:
        if any((pat in padded for pat in pats)):
            return canonical
    return ''

def species_compatible(mlst_species: str, candidate_species: str) -> str:
    if not mlst_species or not candidate_species:
        return 'unknown'
    mlst_hint = parse_species_hint(mlst_species)
    cand_hint = parse_species_hint(candidate_species)
    if not mlst_hint or not cand_hint:
        return 'unknown'
    return 'yes' if mlst_hint == cand_hint else 'no'

def run_blast(query_fasta: Path, ref_fasta: Path, out_dir: Path, threads: int) -> Path:
    blast_dir = out_dir / 'blast'
    blast_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = blast_dir / 'ospc_refs_db'
    blast_tsv = blast_dir / 'candidates_vs_ospc_refs.tsv'
    log('Building BLAST database')
    subprocess.run(['makeblastdb', '-in', str(ref_fasta), '-dbtype', 'nucl', '-out', str(db_prefix)], check=True, stdout=(blast_dir / 'makeblastdb.stdout.log').open('w'), stderr=(blast_dir / 'makeblastdb.stderr.log').open('w'))
    log('Running BLASTN candidate validation')
    outfmt = '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle'
    subprocess.run(['blastn', '-query', str(query_fasta), '-db', str(db_prefix), '-out', str(blast_tsv), '-outfmt', outfmt, '-max_target_seqs', '10', '-num_threads', str(max(1, threads))], check=True, stdout=(blast_dir / 'blastn.stdout.log').open('w'), stderr=(blast_dir / 'blastn.stderr.log').open('w'))
    return blast_tsv

def parse_blast_best(path: Optional[Path]) -> Dict[str, BlastHit]:
    if path is None or not path.exists():
        return {}
    best: Dict[str, BlastHit] = {}
    with path.open('r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 10:
                continue
            qseqid, sseqid = (parts[0], parts[1])
            hit = BlastHit(sseqid=sseqid, pident=to_float(parts[2]), length=to_float(parts[3]), qlen=to_float(parts[4]), slen=to_float(parts[5]), qcovs=to_float(parts[6]), evalue=parts[7], bitscore=to_float(parts[8]), stitle=parts[9], species_hint=parse_species_hint(parts[1] + ' ' + parts[9]))
            previous = best.get(qseqid)
            if previous is None or (hit.bitscore, hit.qcovs, hit.pident) > (previous.bitscore, previous.qcovs, previous.pident):
                best[qseqid] = hit
    return best

@dataclass
class Candidate:
    sample: str
    candidate_id: str
    template: str
    full_header: str
    source: str
    seq: str
    length: int
    kma: Dict[str, object]
    mlst_species: str = ''
    mlst_status: str = ''
    orf: OrfInfo = field(default_factory=OrfInfo)
    blast: BlastHit = field(default_factory=BlastHit)
    fail_reasons: List[str] = field(default_factory=list)
    eligible: bool = False
    score: float = 0.0
    selected: bool = False
    sample_decision: str = ''
    notes: str = ''

    @property
    def depth(self) -> float:
        return float(self.kma.get('depth', float('nan')))

    @property
    def qcov(self) -> float:
        return float(self.kma.get('query_coverage', float('nan')))

    @property
    def tcov(self) -> float:
        return float(self.kma.get('template_coverage', float('nan')))

    @property
    def qid(self) -> float:
        return float(self.kma.get('query_identity', float('nan')))

    @property
    def tid(self) -> float:
        return float(self.kma.get('template_identity', float('nan')))

def source_from_template(template: str) -> str:
    low = template.lower()
    if '_assembled_contig_' in low or 'hybpiper' in low:
        return 'hybpiper_contig'
    if low.startswith('guided_contig_') or ('|node_' in low and low.startswith('guided')):
        return 'guided_contig'
    if low.startswith('seed_contig_') or ('|node_' in low and low.startswith('seed')):
        return 'seed_contig'
    if 'contig' in low or 'node_' in low:
        return 'other_contig'
    if 'consensus' in low or low.endswith('.fsa'):
        return 'kma_consensus'
    return 'reference'

def source_priority(source: str) -> int:
    return {'hybpiper_contig': 9, 'guided_contig': 8, 'seed_contig': 7, 'kma_consensus': 7, 'other_contig': 5, 'reference': 0}.get(source, 0)

def read_db_fasta_map(path: Path) -> Dict[str, Tuple[str, str]]:
    seqs: Dict[str, Tuple[str, str]] = {}
    for header, seq in read_fasta(path):
        full = header.strip()
        first = full.split()[0]
        seqs[full] = (full, seq)
        seqs[first] = (full, seq)
    return seqs

def find_db_fasta_for_res(res_path: Path, sample: str) -> Optional[Path]:
    candidates = [res_path.parent / f'{sample}.personalized_db.fasta', res_path.parent / f'{sample}.personalized_db.fa', res_path.parent / f'{sample}_perso_db.fasta', res_path.parent / f'{sample}_perso_db.fa', res_path.parent / f'{sample}.perso_db.fasta', res_path.parent / f'{sample}.perso_db.fa']
    candidates.extend(sorted(res_path.parent.glob('*.personalized_db.fasta')))
    candidates.extend(sorted(res_path.parent.glob('*.personalized_db.fa')))
    candidates.extend(sorted(res_path.parent.glob('*perso_db.fasta')))
    candidates.extend(sorted(res_path.parent.glob('*perso_db.fa')))
    for p in candidates:
        if p.exists():
            return p
    return None

def find_kma_fsa_for_res(res_path: Path, sample: str) -> Optional[Path]:
    candidates = [res_path.parent / f'{sample}.fsa', res_path.parent / f'{sample}.ref.fsa', res_path.parent / f'{sample}.consensus.fasta', res_path.parent / f'{sample}.consensus.fa']
    candidates.extend(sorted(res_path.parent.glob('*.fsa')))
    candidates.extend(sorted(res_path.parent.glob('*.consensus.fasta')))
    candidates.extend(sorted(res_path.parent.glob('*.consensus.fa')))
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            return p
    return None

def clean_kma_consensus(seq: str, min_len: int, max_n_frac: float) -> Tuple[Optional[str], str]:
    raw = seq.upper().replace('-', '')
    raw = re.sub('\\s+', '', raw)
    if not raw:
        return (None, 'empty_consensus')
    n_frac_raw = raw.count('N') / len(raw)
    if n_frac_raw > max_n_frac and raw.strip('N').count('N') > 0:
        return (None, f'consensus_N_frac>{max_n_frac}')
    trimmed = raw.strip('N')
    if len(trimmed) < min_len:
        return (None, f'consensus_trimmed_len<{min_len}')
    if 'N' in trimmed:
        return (None, 'consensus_internal_N')
    bad = re.sub('[ACGT]', '', trimmed)
    if bad:
        return (None, 'consensus_ambiguous_bases')
    return (trimmed, 'ok')

def get_clean_selected_kma_consensus(selected: 'Candidate', personalized_dir: Path, min_len: int, max_n_frac: float) -> Tuple[Optional[str], str]:
    if selected.source == 'kma_consensus':
        seq = selected.seq.replace('-', '').upper()
        if selected.orf.strand == '-':
            seq = revcomp(seq)
        return (seq, 'selected_candidate_is_kma_consensus')
    sample = selected.sample
    candidate_paths: List[Path] = [personalized_dir / sample / f'{sample}.fsa', personalized_dir / sample / f'{sample}.ref.fsa']
    candidate_paths.extend(sorted((personalized_dir / sample).glob('*.fsa')) if (personalized_dir / sample).exists() else [])
    candidate_paths.extend(sorted(personalized_dir.rglob(f'{sample}.fsa')))
    seen_paths = []
    for path in candidate_paths:
        if path in seen_paths:
            continue
        seen_paths.append(path)
        if not path.exists() or path.stat().st_size == 0:
            continue
        consensus_map = read_db_fasta_map(path)
        keys = [selected.template, selected.template.split()[0]]
        for key in keys:
            if key not in consensus_map:
                continue
            full_header, cons_seq = consensus_map[key]
            clean_seq, reason = clean_kma_consensus(cons_seq, min_len=min_len, max_n_frac=max_n_frac)
            if clean_seq is None:
                return (None, f'KMA_consensus_rejected:{reason}')
            if selected.orf.strand == '-':
                clean_seq = revcomp(clean_seq)
            return (clean_seq, f'KMA_final_consensus_from={path.name};record={sanitize_id(full_header)}')
    return (None, 'no_matching_KMA_final_consensus_in_fsa')

def discover_res_files(personalized_dir: Path) -> List[Path]:
    res_files = []
    for p in personalized_dir.rglob('*.res'):
        sample = p.stem
        if find_db_fasta_for_res(p, sample) is not None:
            res_files.append(p)
    return sorted(set(res_files))

def add_reason(reasons: List[str], condition: bool, reason: str) -> None:
    if condition:
        reasons.append(reason)

def evaluate_candidate(c: Candidate, args: argparse.Namespace, blast_was_run: bool) -> None:
    c.fail_reasons = []
    seq_no_gaps = c.seq.replace('-', '').upper()
    c.length = len(seq_no_gaps)
    c.orf = longest_stop_free_orf(seq_no_gaps)
    n_count = seq_no_gaps.count('N')
    n_frac = n_count / c.length if c.length else 1.0
    if c.source == 'reference':
        c.fail_reasons.append('reference_not_allowed')
    add_reason(c.fail_reasons, c.length < args.min_len, f'len<{args.min_len}')
    add_reason(c.fail_reasons, c.length > args.max_len, f'len>{args.max_len}')
    add_reason(c.fail_reasons, not math.isnan(c.qid) and c.qid < args.min_qid, f'query_identity<{args.min_qid}')
    add_reason(c.fail_reasons, math.isnan(c.qid), 'missing_query_identity')
    if args.min_tid > 0:
        add_reason(c.fail_reasons, not math.isnan(c.tid) and c.tid < args.min_tid, f'template_identity<{args.min_tid}')
        add_reason(c.fail_reasons, math.isnan(c.tid), 'missing_template_identity')
    add_reason(c.fail_reasons, not math.isnan(c.qcov) and c.qcov < args.min_qcov, f'query_coverage<{args.min_qcov}')
    add_reason(c.fail_reasons, math.isnan(c.qcov), 'missing_query_coverage')
    add_reason(c.fail_reasons, not math.isnan(c.tcov) and c.tcov < args.min_tcov, f'template_coverage<{args.min_tcov}')
    add_reason(c.fail_reasons, math.isnan(c.tcov), 'missing_template_coverage')
    add_reason(c.fail_reasons, not math.isnan(c.depth) and c.depth < args.min_depth, f'depth<{args.min_depth}')
    add_reason(c.fail_reasons, math.isnan(c.depth), 'missing_depth')
    add_reason(c.fail_reasons, n_frac > args.max_n_frac, f'N_frac>{args.max_n_frac}')
    add_reason(c.fail_reasons, c.orf.aa_len < args.min_orf_aa, f'ORFaa<{args.min_orf_aa}')
    add_reason(c.fail_reasons, c.orf.aa_len > args.max_orf_aa, f'ORFaa>{args.max_orf_aa}')
    mlst_match = species_compatible(c.mlst_species, c.blast.species_hint)
    if blast_was_run:
        if not c.blast.sseqid:
            c.fail_reasons.append('no_BLAST_hit')
        else:
            if c.blast.qcovs < args.min_blast_qcov:
                c.fail_reasons.append(f'BLAST_qcov<{args.min_blast_qcov}')
            if c.blast.pident < args.min_blast_pident:
                c.fail_reasons.append(f'BLAST_pident<{args.min_blast_pident}')
            if c.blast.length < args.min_blast_len:
                c.fail_reasons.append(f'BLAST_aln_len<{args.min_blast_len}')
    score = 0.0
    score += source_priority(c.source)
    if args.ideal_min_len <= c.length <= args.ideal_max_len:
        score += 8
    elif args.min_len <= c.length <= args.max_len:
        score += 5
    elif c.length >= args.fragment_len:
        score += 1
    else:
        score -= 10
    if not math.isnan(c.qid):
        score += 5 if c.qid >= 98 else 3 if c.qid >= args.min_qid else -8
    if args.min_tid > 0 and (not math.isnan(c.tid)):
        score += 3 if c.tid >= args.min_tid else -5
    if not math.isnan(c.qcov):
        score += 5 if c.qcov >= 98 else 3 if c.qcov >= args.min_qcov else -5
    if not math.isnan(c.tcov):
        score += 5 if c.tcov >= 98 else 3 if c.tcov >= args.min_tcov else -5
    if not math.isnan(c.depth):
        score += 5 if c.depth >= 100 else 3 if c.depth >= args.min_depth else -5
    if args.min_orf_aa <= c.orf.aa_len <= args.max_orf_aa:
        score += 10
    else:
        score -= 20
    if n_frac <= args.max_n_frac:
        score += 2
    else:
        score -= 10
    if blast_was_run:
        if c.blast.sseqid and c.blast.qcovs >= args.min_blast_qcov and (c.blast.length >= args.min_blast_len):
            score += 8
        else:
            score -= 10
        if mlst_match == 'yes':
            score += 3
        elif mlst_match == 'no':
            score -= 3
    c.score = score
    c.eligible = len(c.fail_reasons) == 0

def candidate_sort_key(c: Candidate) -> Tuple[float, float, int, int]:
    depth = c.depth if not math.isnan(c.depth) else -1.0
    return (c.score, depth, c.length, source_priority(c.source))

def depth_value(c: Candidate) -> float:
    return c.depth if not math.isnan(c.depth) else 0.0

def decide_for_sample(sample: str, candidates: List[Candidate], args: argparse.Namespace) -> Tuple[str, Optional[Candidate], str]:
    if not candidates:
        return ('ospC_no_candidates', None, 'no candidate found in personalized DB/KMA')
    eligible_all = sorted([c for c in candidates if c.eligible], key=candidate_sort_key, reverse=True)
    fragments = [c for c in candidates if c.length >= args.fragment_len and (not c.eligible)]
    if not eligible_all:
        if fragments:
            return ('ospC_partial_or_failed_filters', None, f'{len(fragments)} fragment-like candidates, none passed strict filters')
        return ('ospC_no_reliable_sequence', None, 'no candidate passed strict filters')
    eligible = eligible_all
    ignored_consensus = 0
    if args.kma_consensus_policy == 'rescue_only':
        non_consensus = [c for c in eligible_all if c.source != 'kma_consensus']
        if non_consensus:
            eligible = sorted(non_consensus, key=candidate_sort_key, reverse=True)
            ignored_consensus = len(eligible_all) - len(eligible)
    selected = eligible[0]
    n_secondary = len(eligible) - 1
    note_parts = []
    if ignored_consensus:
        note_parts.append(f'kma_consensus_policy=rescue_only ignored {ignored_consensus} eligible KMA consensus candidate(s)')
    if n_secondary:
        note_parts.append(f'selected best eligible candidate among {len(eligible)} eligible decision candidate(s); ranking=score,depth,length,source_priority; selected_score={selected.score:.3f}; selected_depth={depth_value(selected):.4g}')
    else:
        note_parts.append(f'single eligible decision candidate; selected_score={selected.score:.3f}; selected_depth={depth_value(selected):.4g}')
    return ('ospC_selected_best_eligible', selected, '; '.join(note_parts))

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description='Select one reliable ospC sequence per sample from personalized KMA outputs.')
    p.add_argument('--personalized-dir', default='results/ospc/ospc_kma_final', type=Path)
    p.add_argument('--out-dir', default='results/ospc/ospc_selection', type=Path)
    p.add_argument('--ref-fasta', default=Path('data/ospC_ena_myannot.fasta'), type=Path, help='Original ospC reference FASTA for BLAST validation.')
    p.add_argument('--threads', default=4, type=int)
    p.add_argument('--consensus-from', choices=['references', 'all'], default='references', help='Use KMA consensus rows from reference templates only, or from all templates. Default: references.')
    p.add_argument('--consensus-min-len', default=350, type=int, help='Minimum length after trimming terminal Ns for KMA consensus candidates.')
    p.add_argument('--consensus-max-n-frac', default=0.01, type=float, help='Maximum raw N fraction for KMA consensus candidates if internal Ns are present.')
    p.add_argument('--kma-consensus-policy', choices=['all', 'rescue_only'], default='rescue_only', help='How KMA .fsa consensus candidates participate in final conflict resolution. rescue_only keeps them selectable only when no eligible assembled contig exists; all lets them compete with contigs. Default: rescue_only.')
    p.add_argument('--min-len', default=350, type=int)
    p.add_argument('--max-len', default=1100, type=int)
    p.add_argument('--ideal-min-len', default=560, type=int)
    p.add_argument('--ideal-max-len', default=780, type=int)
    p.add_argument('--fragment-len', default=350, type=int)
    p.add_argument('--min-qid', default=95.0, type=float, help='Minimum KMA Query_Identity for candidate selection.')
    p.add_argument('--min-tid', default=0.0, type=float, help='Optional minimum KMA Template_Identity; 0 disables this filter.')
    p.add_argument('--min-qcov', default=95.0, type=float)
    p.add_argument('--min-tcov', default=95.0, type=float)
    p.add_argument('--min-depth', default=20.0, type=float)
    p.add_argument('--max-n-frac', default=0.01, type=float)
    p.add_argument('--min-orf-aa', default=50, type=int)
    p.add_argument('--max-orf-aa', default=235, type=int)
    p.add_argument('--min-blast-qcov', default=70.0, type=float)
    p.add_argument('--min-blast-pident', default=70.0, type=float)
    p.add_argument('--min-blast-len', default=250.0, type=float)
    p.set_defaults(select_dominant=True, include_kma_consensus=True)
    return p

def main() -> None:
    args = build_argparser().parse_args()
    if not args.personalized_dir.exists():
        die(f'Personalized KMA directory not found: {args.personalized_dir}')
    args.out_dir.mkdir(parents=True, exist_ok=True)
    res_files = discover_res_files(args.personalized_dir)
    log(f'Found {len(res_files)} personalized KMA result files')
    all_candidates: List[Candidate] = []
    candidate_id_to_candidate: Dict[str, Candidate] = {}
    for res_path in res_files:
        sample = res_path.stem
        db_fasta = find_db_fasta_for_res(res_path, sample)
        if db_fasta is None:
            log(f'Skipping {sample}: no personalized_db FASTA found near {res_path}')
            continue
        mlst_species = ''
        mlst_status = ''
        kma_rows = parse_kma_res(res_path)
        seq_map = read_db_fasta_map(db_fasta)
        for template, kma in kma_rows.items():
            if template not in seq_map:
                first = template.split()[0]
                if first not in seq_map:
                    continue
                full_header, seq = seq_map[first]
            else:
                full_header, seq = seq_map[template]
            source = source_from_template(template)
            base_id = sanitize_id(template.split('|')[0].split()[0])
            cand_id = sanitize_id(f'{sample}__{base_id}')
            original_id = cand_id
            idx = 2
            while cand_id in candidate_id_to_candidate:
                cand_id = sanitize_id(f'{original_id}_{idx}')
                idx += 1
            c = Candidate(sample=sample, candidate_id=cand_id, template=template, full_header=full_header, source=source, seq=seq.replace('-', '').upper(), length=len(seq.replace('-', '')), kma=kma, mlst_species=mlst_species, mlst_status=mlst_status)
            all_candidates.append(c)
            candidate_id_to_candidate[c.candidate_id] = c
        if args.include_kma_consensus:
            fsa_path = find_kma_fsa_for_res(res_path, sample)
            if fsa_path is None:
                log(f'Sample {sample}: no KMA .fsa consensus file found')
            else:
                consensus_map = read_db_fasta_map(fsa_path)
                n_added = 0
                n_skipped = 0
                for template, kma in kma_rows.items():
                    original_source = source_from_template(template)
                    if args.consensus_from == 'references' and original_source != 'reference':
                        continue
                    if template in consensus_map:
                        full_header, cons_seq = consensus_map[template]
                    else:
                        first = template.split()[0]
                        if first not in consensus_map:
                            n_skipped += 1
                            continue
                        full_header, cons_seq = consensus_map[first]
                    clean_seq, reason = clean_kma_consensus(cons_seq, min_len=args.consensus_min_len, max_n_frac=args.consensus_max_n_frac)
                    if clean_seq is None:
                        n_skipped += 1
                        continue
                    base_id = sanitize_id('kma_consensus_' + template.split()[0])
                    cand_id = sanitize_id(f'{sample}__{base_id}')
                    original_id = cand_id
                    idx = 2
                    while cand_id in candidate_id_to_candidate:
                        cand_id = sanitize_id(f'{original_id}_{idx}')
                        idx += 1
                    c = Candidate(sample=sample, candidate_id=cand_id, template=template, full_header=f'KMA_consensus_from={full_header}', source='kma_consensus', seq=clean_seq, length=len(clean_seq), kma=kma, mlst_species=mlst_species, mlst_status=mlst_status)
                    all_candidates.append(c)
                    candidate_id_to_candidate[c.candidate_id] = c
                    n_added += 1
                log(f'Sample {sample}: added {n_added} clean KMA consensus candidate(s); skipped {n_skipped}')
    if not all_candidates:
        log('No candidates could be matched between KMA .res files and DB FASTA files.')
    preblast_fasta = args.out_dir / 'all_candidates.preblast.fasta'
    write_fasta(((c.candidate_id + ' ' + sanitize_id(c.template), c.seq) for c in all_candidates), preblast_fasta)
    log(f'Wrote candidate FASTA: {preblast_fasta}')
    blast_tsv: Optional[Path] = run_blast(preblast_fasta, args.ref_fasta, args.out_dir, args.threads)
    blast_was_run = True
    best_blast = parse_blast_best(blast_tsv)
    if best_blast:
        for cid, hit in best_blast.items():
            if cid in candidate_id_to_candidate:
                candidate_id_to_candidate[cid].blast = hit
    for c in all_candidates:
        evaluate_candidate(c, args, blast_was_run)
    by_sample: Dict[str, List[Candidate]] = {}
    for c in all_candidates:
        by_sample.setdefault(c.sample, []).append(c)
    sample_rows: List[Dict[str, object]] = []
    selected_contig_records: List[Tuple[str, str]] = []
    selected_orf_records: List[Tuple[str, str]] = []
    selected_consensus_records: List[Tuple[str, str]] = []
    oriented_candidate_records: List[Tuple[str, str]] = []
    for sample in sorted(by_sample):
        candidates = by_sample[sample]
        decision, selected, notes = decide_for_sample(sample, candidates, args)
        for c in candidates:
            c.sample_decision = decision
            c.notes = notes
        selected_consensus_status = ''
        selected_consensus_len = ''
        if selected is not None:
            selected.selected = True
            header = sanitize_id(f"{sample}|{decision}|{selected.candidate_id}|source={selected.source}|len={selected.length}|depth={selected.depth:.2f}|orfaa={selected.orf.aa_len}|mlst={selected.mlst_species or 'NA'}|blast={selected.blast.species_hint or selected.blast.sseqid or 'NA'}")
            selected_contig_records.append((header, selected.orf.oriented_seq or selected.seq))
            selected_orf_records.append((header, selected.orf.orf_nt or selected.orf.oriented_seq or selected.seq))
            consensus_seq, consensus_status = get_clean_selected_kma_consensus(selected, args.personalized_dir, args.consensus_min_len, args.consensus_max_n_frac)
            selected_consensus_status = consensus_status
            if consensus_seq is not None:
                selected_consensus_len = str(len(consensus_seq))
                consensus_header = sanitize_id(header + '|tree_sequence=KMA_final_consensus')
                selected_consensus_records.append((consensus_header, consensus_seq))
            else:
                log(f'Sample {sample}: selected candidate passed QC but no clean final KMA consensus was found ({consensus_status}); excluded from consensus tree FASTA')
        for c in candidates:
            oriented_candidate_records.append((c.candidate_id, c.orf.oriented_seq or c.seq))
        sample_rows.append({'sample': sample, 'mlst_species': candidates[0].mlst_species if candidates else '', 'mlst_status': candidates[0].mlst_status if candidates else '', 'n_candidates': len(candidates), 'n_eligible': sum((1 for c in candidates if c.eligible)), 'decision': decision, 'selected_candidate_id': selected.candidate_id if selected else '', 'selected_template': selected.template if selected else '', 'selected_source': selected.source if selected else '', 'selected_len': selected.length if selected else '', 'selected_depth': f'{selected.depth:.4g}' if selected else '', 'selected_orf_aa': selected.orf.aa_len if selected else '', 'selected_consensus_len': selected_consensus_len, 'selected_consensus_status': selected_consensus_status, 'selected_blast_hit': selected.blast.sseqid if selected else '', 'selected_blast_species': selected.blast.species_hint if selected else '', 'notes': notes})
    candidate_summary = args.out_dir / 'ospc_candidate_summary.tsv'
    candidate_fields = ['sample', 'mlst_species', 'mlst_status', 'candidate_id', 'template', 'source', 'length', 'score', 'query_identity', 'query_coverage', 'template_identity', 'template_coverage', 'depth', 'orf_aa_len', 'orf_strand', 'orf_frame', 'blast_hit', 'blast_pident', 'blast_qcov', 'blast_aln_len', 'blast_species_hint', 'mlst_blast_species_match', 'eligible', 'selected', 'sample_decision', 'fail_reasons', 'notes']
    with candidate_summary.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=candidate_fields, delimiter='\t')
        writer.writeheader()
        for c in sorted(all_candidates, key=lambda x: (x.sample, not x.selected, -x.score)):
            writer.writerow({'sample': c.sample, 'mlst_species': c.mlst_species, 'mlst_status': c.mlst_status, 'candidate_id': c.candidate_id, 'template': c.template, 'source': c.source, 'length': c.length, 'score': f'{c.score:.3f}', 'query_identity': f'{c.qid:.4g}' if not math.isnan(c.qid) else '', 'query_coverage': f'{c.qcov:.4g}' if not math.isnan(c.qcov) else '', 'template_identity': f'{c.tid:.4g}' if not math.isnan(c.tid) else '', 'template_coverage': f'{c.tcov:.4g}' if not math.isnan(c.tcov) else '', 'depth': f'{c.depth:.4g}' if not math.isnan(c.depth) else '', 'orf_aa_len': c.orf.aa_len, 'orf_strand': c.orf.strand, 'orf_frame': c.orf.frame, 'blast_hit': c.blast.sseqid, 'blast_pident': f'{c.blast.pident:.4g}' if not math.isnan(c.blast.pident) else '', 'blast_qcov': f'{c.blast.qcovs:.4g}' if not math.isnan(c.blast.qcovs) else '', 'blast_aln_len': f'{c.blast.length:.4g}' if not math.isnan(c.blast.length) else '', 'blast_species_hint': c.blast.species_hint, 'mlst_blast_species_match': species_compatible(c.mlst_species, c.blast.species_hint), 'eligible': 'yes' if c.eligible else 'no', 'selected': 'yes' if c.selected else 'no', 'sample_decision': c.sample_decision, 'fail_reasons': ';'.join(c.fail_reasons), 'notes': c.notes})
    sample_summary = args.out_dir / 'ospc_sample_decisions.tsv'
    with sample_summary.open('w', encoding='utf-8', newline='') as handle:
        fields = list(sample_rows[0].keys()) if sample_rows else []
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in sample_rows:
            writer.writerow(row)
    write_fasta(selected_contig_records, args.out_dir / 'selected_ospC_contigs_oriented.fasta')
    write_fasta(selected_orf_records, args.out_dir / 'selected_ospC_orf.fasta')
    write_fasta(selected_consensus_records, args.out_dir / 'selected_ospC_consensus_oriented.fasta')
    write_fasta(oriented_candidate_records, args.out_dir / 'all_candidates.oriented.fasta')
    n_selected = len(selected_consensus_records)
    n_samples = len(sample_rows)
    log(f'Wrote: {candidate_summary}')
    log(f'Wrote: {sample_summary}')
    log(f"Wrote: {args.out_dir / 'selected_ospC_consensus_oriented.fasta'}")
    log(f'Selected {n_selected}/{n_samples} samples for ospC consensus phylogeny')
    if n_selected == 0:
        log('No selected sequence. Try inspecting ospc_candidate_summary.tsv or relaxing filters carefully.')
if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import csv
import re
import sys
from pathlib import Path
from collections import Counter

def read_fasta(path):
    name = None
    seq = []
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield (name, ''.join(seq).upper())
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if name is not None:
            yield (name, ''.join(seq).upper())

def write_fasta(records, path):
    with open(path, 'w', encoding='utf-8') as out:
        for name, seq in records:
            out.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + '\n')

def safe_header(text, max_len=180):
    x = text.strip()
    x = x.replace('|', '|')
    x = re.sub('\\s+', '_', x)
    x = re.sub('[^A-Za-z0-9_.|:+\\-]+', '_', x)
    x = re.sub('_+', '_', x).strip('_')
    return x[:max_len] if x else 'unknown'

def prepare_refs(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    kma_records = []
    hyb_records = []
    kept = 0
    skipped = 0
    with open(args.map_tsv, 'w', newline='', encoding='utf-8') as mapf:
        writer = csv.writer(mapf, delimiter='\t')
        writer.writerow(['kma_id', 'hybpiper_id', 'original_header', 'length'])
        for header, seq in read_fasta(args.input):
            seq = seq.upper().replace('U', 'T').replace('-', '')
            seq = re.sub('[^ACGTN]', 'N', seq)
            ungapped_len = len(seq)
            if ungapped_len < args.min_len or ungapped_len > args.max_len:
                skipped += 1
                continue
            kept += 1
            sid = f'ref_{kept:06d}'
            kma_id = f'{sid}|{safe_header(header)}'
            hyb_id = f'{sid}-ospC'
            kma_records.append((kma_id, seq))
            hyb_records.append((hyb_id, seq))
            writer.writerow([kma_id, hyb_id, header, ungapped_len])
    write_fasta(kma_records, args.kma_fasta)
    write_fasta(hyb_records, args.hybpiper_fasta)
    print(f'kept\t{kept}')
    print(f'skipped\t{skipped}')

def read_res_rows(path):
    rows = []
    header = None
    default = ['Template', 'Score', 'Expected', 'Template_length', 'Template_Identity', 'Template_Coverage', 'Query_Identity', 'Query_Coverage', 'Depth', 'q_value', 'p_value']
    with open(path, encoding='utf-8', errors='replace') as f:
        for raw in f:
            line = raw.rstrip('\n')
            if not line:
                continue
            if line.startswith('#'):
                clean = line.lstrip('#').strip()
                if clean.lower().startswith('template'):
                    header = clean.split('\t')
                continue
            parts = line.split('\t')
            cols = header if header and len(header) <= len(parts) else default
            d = {cols[i]: parts[i] for i in range(min(len(cols), len(parts)))}
            if 'Template' not in d:
                d['Template'] = parts[0]
            rows.append(d)
    return rows

def fnum(x):
    try:
        return float(str(x).replace(',', '.'))
    except Exception:
        return float('nan')

def fasta_dict(path):
    d = {}
    for h, s in read_fasta(path):
        d[h] = (h, s)
        d[h.split()[0]] = (h, s)
    return d

def extract_best_ref(args):
    rows = read_res_rows(args.res)
    if not rows:
        return 1
    rows.sort(key=lambda r: fnum(r.get('Score', 'nan')), reverse=True)
    wanted = rows[0]['Template']
    db = fasta_dict(args.refs)
    rec = db.get(wanted) or db.get(wanted.split()[0])
    if rec is None:
        return 2
    write_fasta([rec], args.output)
    print(f'{wanted}\t{len(rec[1])}')
    return 0

def clean_consensus(seq, min_len, max_n_frac):
    raw = seq.upper().replace('-', '')
    raw = re.sub('\\s+', '', raw)
    if not raw:
        return (None, 'empty')
    trimmed = raw.strip('N')
    if len(trimmed) < min_len:
        return (None, 'too_short')
    if 'N' in trimmed:
        return (None, 'internal_N')
    n_frac = raw.count('N') / len(raw)
    if n_frac > max_n_frac and raw.strip('N').count('N') > 0:
        return (None, 'too_many_N')
    bad = re.sub('[ACGT]', '', trimmed)
    if bad:
        return (None, 'ambiguous_base')
    return (trimmed, 'ok')

def choose_bait(args):
    rows = read_res_rows(args.res)
    if not rows:
        return 1
    rows.sort(key=lambda r: fnum(r.get('Score', 'nan')), reverse=True)
    best = rows[0]
    template = best['Template']
    if args.fsa and Path(args.fsa).exists():
        cons = fasta_dict(args.fsa)
        rec = cons.get(template) or cons.get(template.split()[0])
        if rec is not None:
            clean, reason = clean_consensus(rec[1], args.cons_min_len, args.cons_max_n_frac)
            if clean is not None:
                write_fasta([(f'bait_consensus|{rec[0]}', clean)], args.output)
                print(f'consensus\t{template}\t{len(clean)}')
                return 0
    db = fasta_dict(args.refs)
    rec = db.get(template) or db.get(template.split()[0])
    if rec is None:
        return 2
    write_fasta([(f'bait_reference|{rec[0]}', rec[1])], args.output)
    print(f'reference\t{template}\t{len(rec[1])}')
    return 0

def filter_refs(args):
    records = []
    for h, s in read_fasta(args.input):
        seq = s.upper().replace('U', 'T').replace('-', '')
        if len(seq) >= args.min_len:
            records.append((h, seq))
    write_fasta(records, args.output)
    print(len(records))

def build_seeds(args):
    records = list(read_fasta(args.alignment))
    if not records:
        raise SystemExit('No aligned records')
    seqs = [s.upper() for _, s in records]
    aln_len = len(seqs[0])
    if any((len(s) != aln_len for s in seqs)):
        raise SystemExit('Alignment length mismatch')
    seeds = []
    seen = set()
    for start in range(0, aln_len - args.window + 1, args.step):
        end = start + args.window
        windows = [s[start:end] for s in seqs]
        total = len(windows) * args.window
        gap_frac = sum((w.count('-') for w in windows)) / total
        if gap_frac > args.max_gap_frac:
            continue
        informative = 0
        ident_sum = 0.0
        consensus = []
        for i in range(args.window):
            col = [w[i] for w in windows if w[i] not in '-N']
            if not col:
                consensus.append('N')
                continue
            informative += 1
            base, count = Counter(col).most_common(1)[0]
            consensus.append(base)
            ident_sum += count / len(col)
        if informative < int(args.window * 0.7):
            continue
        mean_ident = ident_sum / informative
        if mean_ident < args.min_ident:
            continue
        seed = ''.join(consensus).replace('N', '')
        if len(seed) < max(21, int(args.window * 0.7)):
            continue
        if seed in seen:
            continue
        seen.add(seed)
        seeds.append((f'seed_{len(seeds) + 1}', seed, start + 1, end, mean_ident, gap_frac))
    if len(seeds) < args.min_seeds:
        raise SystemExit(f'Only {len(seeds)} seed(s) found')
    with open(args.output, 'w', encoding='utf-8') as fa, open(args.report, 'w', encoding='utf-8') as rep:
        rep.write('seed_id\tseed_len\taln_start\taln_end\tmean_identity\tgap_fraction\tsequence\n')
        for sid, seq, a, b, ident, gapf in seeds:
            fa.write(f'>{sid}\n{seq}\n')
            rep.write(f'{sid}\t{len(seq)}\t{a}\t{b}\t{ident:.4f}\t{gapf:.4f}\t{seq}\n')
    print(len(seeds))

def selected_refs_from_base_kma(rows, qi_min, qc_min, depth_min, include_best_if_none):
    selected = []
    for r in rows:
        qi = fnum(r.get('Query_Identity', 'nan'))
        qc = fnum(r.get('Query_Coverage', 'nan'))
        depth = fnum(r.get('Depth', 'nan'))
        if qi >= qi_min and qc >= qc_min and (depth >= depth_min):
            selected.append(r['Template'])
    if not selected and include_best_if_none and rows:
        rows2 = sorted(rows, key=lambda r: fnum(r.get('Score', 'nan')), reverse=True)
        selected.append(rows2[0]['Template'])
    out = []
    seen = set()
    for x in selected:
        if x not in seen:
            out.append(x)
            seen.add(x)
    return out

def add_contigs(records, path, sample, source, min_len):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return 0
    count = 0
    for h, s in read_fasta(p):
        seq = s.upper().replace('-', '')
        if len(seq) < min_len:
            continue
        count += 1
        clean_h = safe_header(h, 120)
        name = f'{source}_contig_{count}|{sample}|{clean_h}'
        records.append((name, seq))
    return count

def make_final_db(args):
    rows = read_res_rows(args.base_res)
    ref_db = fasta_dict(args.refs)
    selected = selected_refs_from_base_kma(rows, args.ref_qi_min, args.ref_qc_min, args.ref_depth_min, args.include_best_ref_if_none)
    records = []
    n_refs = 0
    for wanted in selected:
        rec = ref_db.get(wanted) or ref_db.get(wanted.split()[0])
        if rec is not None:
            records.append(rec)
            n_refs += 1
    n_seed = add_contigs(records, args.seed_contigs, args.sample, 'seed', args.contig_min_len)
    n_guided = add_contigs(records, args.guided_contigs, args.sample, 'guided', args.contig_min_len)
    n_hyb = add_contigs(records, args.hybpiper_contigs, args.sample, 'hybpiper', args.contig_min_len)
    write_fasta(records, args.output)
    print(f'{args.sample}\t{n_refs}\t{n_seed}\t{n_guided}\t{n_hyb}\t{len(records)}')

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prepare_refs')
    p.add_argument('--input', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--kma-fasta', required=True)
    p.add_argument('--hybpiper-fasta', required=True)
    p.add_argument('--map-tsv', required=True)
    p.add_argument('--min-len', type=int, default=450)
    p.add_argument('--max-len', type=int, default=1100)
    p.set_defaults(func=prepare_refs)
    p = sub.add_parser('extract_best_ref')
    p.add_argument('--res', required=True)
    p.add_argument('--refs', required=True)
    p.add_argument('--output', required=True)
    p.set_defaults(func=extract_best_ref)
    p = sub.add_parser('choose_bait')
    p.add_argument('--res', required=True)
    p.add_argument('--fsa', default='')
    p.add_argument('--refs', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--cons-min-len', type=int, default=300)
    p.add_argument('--cons-max-n-frac', type=float, default=0.02)
    p.set_defaults(func=choose_bait)
    p = sub.add_parser('filter_refs')
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--min-len', type=int, default=500)
    p.set_defaults(func=filter_refs)
    p = sub.add_parser('build_seeds')
    p.add_argument('--alignment', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--report', required=True)
    p.add_argument('--window', type=int, default=41)
    p.add_argument('--step', type=int, default=5)
    p.add_argument('--min-ident', type=float, default=0.85)
    p.add_argument('--max-gap-frac', type=float, default=0.2)
    p.add_argument('--min-seeds', type=int, default=1)
    p.set_defaults(func=build_seeds)
    p = sub.add_parser('make_final_db')
    p.add_argument('--sample', required=True)
    p.add_argument('--base-res', required=True)
    p.add_argument('--refs', required=True)
    p.add_argument('--seed-contigs', required=True)
    p.add_argument('--guided-contigs', required=True)
    p.add_argument('--hybpiper-contigs', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--ref-qi-min', type=float, default=60)
    p.add_argument('--ref-qc-min', type=float, default=95)
    p.add_argument('--ref-depth-min', type=float, default=1)
    p.add_argument('--include-best-ref-if-none', type=int, default=1)
    p.add_argument('--contig-min-len', type=int, default=1)
    p.set_defaults(func=make_final_db)
    args = parser.parse_args()
    ret = args.func(args)
    if isinstance(ret, int):
        sys.exit(ret)
if __name__ == '__main__':
    main()

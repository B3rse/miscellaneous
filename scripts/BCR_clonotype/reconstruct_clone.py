#!/usr/bin/env python3
"""Reconstruct the observed V-region upstream of a BCR CDR3 from Archer molbar reads.

One vote per molbar family. Inputs are trusted; if something is wrong, it crashes.
"""

import argparse
import csv
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def rc(seq):
    return seq.translate(COMPLEMENT)[::-1]


def fastp_poly_g_trim(seq, qual, min_len=10):
    """Return (trimmed_seq, trimmed_qual, bases_removed) using fastp's rule."""
    n = len(seq)
    mismatches = 0
    examined = 0
    first_g = n - 1
    while examined < n:
        pos = n - examined - 1
        if seq[pos] != "G":
            mismatches += 1
        else:
            first_g = pos
        if mismatches > 5 or (mismatches > (examined + 1) // 8 and examined >= min_len - 1):
            break
        examined += 1
    if examined < min_len:
        return seq, qual, 0
    trimmed_qual = qual if qual == "*" else qual[:first_g]
    return seq[:first_g], trimmed_qual, n - first_g


def trim_poly_g_bam(flag, seq, qual, min_len=10):
    """Poly-G trim in sequencing-cycle orientation, return BAM-oriented values."""
    if seq == "*":
        return seq, qual, 0
    reverse = bool(flag & 0x10)
    machine_seq = rc(seq) if reverse else seq
    machine_qual = (qual[::-1] if qual != "*" else "*") if reverse else qual
    tseq, tqual, removed = fastp_poly_g_trim(machine_seq, machine_qual, min_len)
    if removed == 0:
        return seq, qual, 0
    if reverse:
        return rc(tseq), (tqual[::-1] if tqual != "*" else "*"), removed
    return tseq, tqual, removed


def anchor_on_cdr3(seq, qual, cdr3):
    """Return (oriented_seq, oriented_qual, cdr3_index) or None."""
    reverse_qual = qual[::-1] if qual != "*" else "*"
    for candidate, candidate_qual in ((seq, qual), (rc(seq), reverse_qual)):
        i = candidate.find(cdr3)
        if i >= 0 and candidate.find(cdr3, i + 1) < 0:
            return candidate, candidate_qual, i
    return None


def find_overlap(ref, query, query_qual, min_overlap, max_mismatch):
    """Return the unique best placement, or None."""
    best_score = None
    best = []
    reverse_qual = query_qual[::-1] if query_qual != "*" else "*"
    for q, qq in ((query, query_qual), (rc(query), reverse_qual)):
        for offset in range(-(len(q) - min_overlap), len(ref) - min_overlap + 1):
            rs, qs = max(0, offset), max(0, -offset)
            overlap = min(len(ref) - rs, len(q) - qs)
            if overlap < min_overlap:
                continue
            matches = mismatches = 0
            for k in range(overlap):
                l, r = ref[rs + k], q[qs + k]
                if l == "N" or r == "N":
                    continue
                if l == r:
                    matches += 1
                else:
                    mismatches += 1
            informative = matches + mismatches
            if informative < min_overlap or mismatches / informative > max_mismatch:
                continue
            score = 2 * matches - 3 * mismatches
            placement = (offset, q, qq)
            if best_score is None or score > best_score:
                best_score, best = score, [placement]
            elif score == best_score:
                best.append(placement)
    if len(best) != 1:
        return None
    return best[0]


def reconstruct_molbar(records, cdr3, min_overlap, max_mismatch):
    """Return {position: base} in CDR3-relative coordinates for one molbar."""
    r1 = next((r for r in records if r[0] & 0x40), None)
    r2 = next((r for r in records if r[0] & 0x80), None)
    positioned = []
    anchored = []
    unanchored = []
    for read in (r1, r2):
        if read is None:
            continue
        flag, seq, qual = read
        if seq == "*":
            continue
        result = anchor_on_cdr3(seq, qual, cdr3)
        if result:
            oseq, oqual, idx = result
            positioned.append((oseq, oqual, -idx))
            anchored.append((oseq, -idx))
        else:
            unanchored.append((seq, qual))
    if not anchored:
        return {}
    ref_seq, ref_start = max(anchored, key=lambda a: len(a[0]))
    for seq, qual in unanchored:
        placement = find_overlap(ref_seq, seq, qual, min_overlap, max_mismatch)
        if placement:
            offset, q, qq = placement
            positioned.append((q, qq, ref_start + offset))

    obs = defaultdict(list)
    for seq, qual, start in positioned:
        for i, base in enumerate(seq):
            q = ord(qual[i]) - 33 if qual != "*" and i < len(qual) else 0
            obs[start + i].append((base, max(0, q)))

    merged = {}
    for pos, calls in obs.items():
        non_n = [c for c in calls if c[0] != "N"]
        if not non_n:
            merged[pos] = "N"
            continue
        unique = {b for b, _ in non_n}
        if len(unique) == 1:
            merged[pos] = non_n[0][0]
            continue
        top = sorted(non_n, key=lambda x: -x[1])
        if len(top) < 2 or top[0][1] - top[1][1] >= 3:
            merged[pos] = top[0][0]
        else:
            merged[pos] = "N"
    return merged


def build_coverage(molecules):
    """Return upstream consensus rows; coverage counts callable A/C/G/T bases."""
    counts = defaultdict(Counter)
    for calls in molecules:
        for pos, base in calls.items():
            if pos < 0:
                counts[pos][base if base in "ACGT" else "N"] += 1
    rows = []
    for pos in sorted(counts):
        c = counts[pos]
        cov = sum(c[b] for b in "ACGT")
        if cov:
            top = max(c[b] for b in "ACGT")
            winners = [b for b in "ACGT" if c[b] == top]
            cons = winners[0] if len(winners) == 1 else "N"
            fraction = top / cov
        else:
            cons, fraction = "N", 0.0
        rows.append((pos, cons, cov, c, fraction))
    return rows


def observed_sequence(coverage_rows, min_cov):
    """Longest contiguous well-supported walk from -1 upstream."""
    by_pos = {row[0]: row for row in coverage_rows}
    bases = []
    pos = -1
    while pos in by_pos:
        _, cons, cov, _, _ = by_pos[pos]
        if cov < min_cov or cons not in "ACGT":
            break
        bases.append(cons)
        pos -= 1
    return "".join(reversed(bases))


def load_qnames(mapping_path, cdr3):
    qnames = []
    for line in Path(mapping_path).read_text().splitlines():
        if not line.strip():
            continue
        qname, seq = line.split("\t")[:2]
        if seq.strip().upper() == cdr3:
            qnames.append(qname)
    return list(dict.fromkeys(qnames))


def load_bam_records(bam_path, qname_file, trim_poly_g, poly_g_min):
    """Return {qname: [(flag, seq, qual), ...]}."""
    output = subprocess.check_output(
        ["samtools", "view", "-F", "2304", "-N", str(qname_file), str(bam_path)],
        text=True,
    )
    records = defaultdict(list)
    for line in output.splitlines():
        fields = line.split("\t")
        qname = fields[0]
        flag = int(fields[1])
        seq = fields[9].upper()
        qual = fields[10]
        if trim_poly_g:
            seq, qual, _ = trim_poly_g_bam(flag, seq, qual, poly_g_min)
        records[qname].append((flag, seq, qual))
    return records


def read_mixcr_row(path, cdr3):
    """Return the TSV row whose CDR3 matches, or None."""
    if path is None or not Path(path).exists():
        return None
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for col in ("nSeqCDR3", "cdr3nt", "cdr3", "nucleotideCDR3", "nSeqImputedCDR3"):
                if col in row and row[col].strip().upper() == cdr3:
                    return row
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--cdr3", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mixcr-results", type=Path)
    parser.add_argument("--mixcr-contigs", type=Path)
    parser.add_argument("--min-coverage", type=int, default=10)
    parser.add_argument("--min-overlap", type=int, default=20)
    parser.add_argument("--max-mismatch", type=float, default=0.10)
    parser.add_argument("--trim-poly-g", action="store_true")
    parser.add_argument("--poly-g-min-length", type=int, default=10)
    args = parser.parse_args()

    cdr3 = args.cdr3.upper()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    qnames = load_qnames(args.mapping, cdr3)
    if not qnames:
        raise SystemExit(f"No mapping entries found for CDR3 {cdr3}")
    print(f"Selected {len(qnames)} molbar families for CDR3", file=sys.stderr)

    qname_file = args.output_dir / "clone.qnames.txt"
    qname_file.write_text("\n".join(qnames) + "\n")

    records = load_bam_records(
        args.bam, qname_file, args.trim_poly_g, args.poly_g_min_length
    )

    molecules = []
    for qname in qnames:
        recs = records.get(qname, [])
        if not recs:
            continue
        molecules.append(reconstruct_molbar(recs, cdr3, args.min_overlap, args.max_mismatch))

    anchored = sum(1 for m in molecules if m)
    if not anchored:
        raise SystemExit("No selected molbar families could be anchored on the CDR3")
    print(f"Anchored {anchored} / {len(qnames)} families on the CDR3", file=sys.stderr)

    coverage = build_coverage(molecules)

    with (args.output_dir / "clone.coverage.tsv").open("w") as handle:
        handle.write("position\tconsensus\tcoverage\tA\tC\tG\tT\tN\tconsensus_fraction\n")
        for pos, cons, cov, c, frac in coverage:
            handle.write(
                f"{pos}\t{cons}\t{cov}\t{c['A']}\t{c['C']}\t{c['G']}\t{c['T']}\t{c['N']}\t{frac:.4f}\n"
            )

    v_seq = observed_sequence(coverage, args.min_coverage)
    with (args.output_dir / "clone.observed_V.fasta").open("w") as handle:
        handle.write(f">clone|CDR3={cdr3}|min_coverage={args.min_coverage}|length={len(v_seq)}\n")
        for i in range(0, len(v_seq), 80):
            handle.write(v_seq[i : i + 80] + "\n")

    print("V length at coverage thresholds:", file=sys.stderr)
    for threshold in (1, 5, 10, 20, 50):
        print(f"  >={threshold:>3}: {len(observed_sequence(coverage, threshold))} nt", file=sys.stderr)

    mixcr_row = read_mixcr_row(args.mixcr_results, cdr3)
    if mixcr_row:
        for col in ("cloneId", "cloneCount", "cloneFraction",
                    "allVHitsWithScore", "allJHitsWithScore", "allCHitsWithScore"):
            if col in mixcr_row and mixcr_row[col].strip():
                print(f"  mixcr {col}: {mixcr_row[col]}", file=sys.stderr)
    contig_row = read_mixcr_row(args.mixcr_contigs, cdr3)
    if contig_row and "nSeqImputedFR3" in contig_row:
        fr3 = contig_row["nSeqImputedFR3"]
        obs = sum(c.isupper() for c in fr3 if c.isalpha())
        imp = sum(c.islower() for c in fr3 if c.isalpha())
        print(f"  mixcr FR3: {obs} observed, {imp} imputed", file=sys.stderr)

    print(f"Wrote outputs to {args.output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()

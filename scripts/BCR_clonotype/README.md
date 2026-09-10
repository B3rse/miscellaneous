# Archer BCR clonotype reconstruction

`reconstruct_clone.py` recovers the experimentally observed V-region
sequence upstream of a specified BCR CDR3 from Archer deduplicated
molbar reads. Each Archer molecular family contributes one vote per
position at the clone-level consensus; family `Depth` is not used.

Personal analysis script — inputs are trusted, errors crash the run.

## Requirements

- Python 3.9+
- `samtools` on `$PATH`

## Usage

```bash
python reconstruct_clone.py \
  --bam sample.molbar.trimmed.deduped.merged.bam \
  --mapping sample.molbar.trimmed.deduped.Immune_repertoire.mixcr_read_ID_to_clonotype_mapping.txt \
  --cdr3 TGTGCAACGGGGGATTACCATTCCCCCTTTGACTCCTGG \
  --output-dir out/ \
  --min-coverage 10 \
  --trim-poly-g
```

Optional MiXCR comparison inputs:

```
--mixcr-results  sample.molbar.trimmed.deduped.Immune_repertoire.mixcr_results.txt
--mixcr-contigs  sample.molbar.trimmed.deduped.Immune_repertoire.mixcr_contig_results.txt
```

## Outputs

- `clone.qnames.txt` — MiXCR-selected molbar QNAMEs
- `clone.coverage.tsv` — per-position base counts and consensus
- `clone.observed_V.fasta` — contiguous well-supported V sequence at
  `--min-coverage`

Coverage-threshold profile, MiXCR summary, and family counts are
printed to stderr.

Positions in `clone.coverage.tsv` are relative to the first CDR3 base:
`-1` is immediately upstream of CDR3.
The `coverage` and `consensus_fraction` columns use callable A/C/G/T
observations; unresolved `N` calls are reported separately and do not
satisfy `--min-coverage`.

## Poly-G trimming

Off by default; enable with `--trim-poly-g`. Trimming runs in
sequencing-cycle orientation so it works correctly on reverse-strand
BAM records. The rule follows fastp's `trimPolyG`: at most one non-G
per 8 examined tail bases, up to 5 total, once the tail reaches
`--poly-g-min-length` (default 10).

## Scope limits

- CDR3 anchoring is exact-string matching.
- Mate placement is ungapped overlap.
- Tied best mate-overlap placements are treated as ambiguous and dropped.
- Mates that neither anchor on the CDR3 nor overlap the anchored mate
  are dropped.

#!/usr/bin/env python3
"""Label each pilot-set entry with what it can actually measure.

WHY THIS EXISTS
    The 18 pilot entries are not 18 comparable test cases. ICEberg records an
    element as (accession, start, end), and that accession is sometimes a whole
    chromosome and sometimes a standalone deposit of the element by itself. Those
    two support completely different measurements, and averaging them produces a
    recall figure that means nothing.

    Concretely, AB450045 is deposited as "SXT conjugative element GENE CLUSTER":
    19 kb of an element that is ~100 kb in full. It carries the integrase and the
    AMR cassette but no relaxase, no coupling protein and no T4SS. A caller that
    works by finding conjugation machinery cannot find machinery that is not in
    the file. Scoring that as a miss would be measuring the deposit, not the code.

CLASSES
    chromosomal  the accession is a genome/chromosome and the element sits inside
                 it with flanking sequence on both sides. Both DETECTION and
                 BOUNDARY OFFSET are measurable. This is the real test.
    standalone   the accession is the element itself (element length is most of
                 the record). Detection is measurable; boundary offset is not,
                 because "the whole record" is trivially the right answer and
                 there is no flank for an att site to be found in.

INPUT : pilot_set.tsv, genomes/{accession}.fna
OUTPUT: pilot_set_classified.tsv — same rows plus record_bp, flank_bp, class
        Consumed by score.py, which reports the two classes separately.
"""
import os
import sys

BENCH = os.path.dirname(os.path.abspath(__file__))

# A record counts as "standalone" when the element covers essentially all of it.
# 0.90 rather than 1.0 because a few deposits carry a kilobase or two of flank
# (GU725392 starts the element at position 77; KX077897 leaves ~1.3 kb each side).
STANDALONE_COVERAGE = 0.90


def record_length(path):
    """Total bases in a FASTA. Counted line by line: these files hold a single
    multi-megabase sequence, and slurping it whole is how the earlier
    'decontamination stripped 99% of the assembly' measurement error happened."""
    total = 0
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def main():
    # Which pilot set to classify, and where to write the result. Defaults to the
    # ICE set so existing invocations are unchanged; the IME pilot passes its own
    # pair of paths rather than overwriting the ICE benchmark's inputs.
    in_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BENCH, "pilot_set.tsv")
    out_path = (sys.argv[2] if len(sys.argv) > 2
                else os.path.join(BENCH, "pilot_set_classified.tsv"))

    rows = []
    with open(in_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))

    out_rows = []
    for row in rows:
        genome = os.path.join(BENCH, "genomes", row["accession"] + ".fna")
        if not os.path.isfile(genome):
            print(f"MISSING {genome}", file=sys.stderr)
            continue

        record_bp = record_length(genome)
        start, end = int(row["start"]), int(row["end"])
        element_bp = end - start + 1
        # How much sequence lies outside the element — the room an att search has
        # to work in, and the room a boundary error has to be visible in.
        flank_bp = record_bp - element_bp
        coverage = element_bp / record_bp if record_bp else 0.0

        row["record_bp"] = str(record_bp)
        row["flank_bp"] = str(flank_bp)
        row["element_covers_record"] = f"{coverage:.3f}"
        row["class"] = "standalone" if coverage >= STANDALONE_COVERAGE else "chromosomal"
        out_rows.append(row)

    if not out_rows:
        sys.exit(f"no usable rows in {in_path} - are the genomes fetched?")

    # Take the header from an OUTPUT row, not from rows[0]. rows[0] is an INPUT
    # row, and if its genome was missing it was skipped above and never gained
    # record_bp / flank_bp / element_covers_record / class - so the file would be
    # written without those columns and score.py would die on KeyError 'class'.
    # A silent header shift like this is exactly the failure mode that is hard to
    # notice, because the file still looks well-formed.
    columns = list(out_rows[0].keys())
    with open(out_path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in out_rows:
            handle.write("\t".join(row[c] for c in columns) + "\n")

    n_chrom = sum(1 for r in out_rows if r["class"] == "chromosomal")
    n_alone = sum(1 for r in out_rows if r["class"] == "standalone")
    print(f"chromosomal (detection + boundary): {n_chrom}")
    print(f"standalone  (detection only)      : {n_alone}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()

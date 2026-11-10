"""Tests for adapter_screen.py.

What the script does, in one sentence: it writes the files blastn needs to search the finished
genome for adapters (circular contigs extended across their start), and turns blastn's hits
into adapter loci with a counted and a reason column, plus one row of the MultiQC Adapters table.

The promises these tests pin down:
  * prepare copies every adapter record into one query file, and extends each circular contig,
    and only those, by its own first 200 bp;
  * an end hit on a linear contig counts from 20 bp, an internal hit from 28 bp, both at 90 %
    identity or more; shorter or less similar hits are dropped;
  * a circular contig has no ends: a hit near its end is judged as internal, a hit across its
    start is one locus whose end lies before its start, and a hit wholly in the appended copy
    is not counted twice;
  * a hit confined to the 50 bp phage Mu end is not counted inside a contig (a Mu prophage holds
    it) but counts at the end of a linear contig; a locus that also covers the barcode counts as
    an adapter;
  * the adapters of a locus are listed best hit first;
  * a contig end is counted once, however many loci it holds;
  * the MultiQC row carries the read-level numbers of a trimmed ONT sample, and NA otherwise.

No BLAST is run: each test writes a small genome, the query file and a BLAST table by hand, in
the column order rule blast_adapters asks blastn for. Files go to pytest's tmp_path.

Run: pytest workflow/scripts/20_qc/test_adapter_screen.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import adapter_screen as screen

ADAPTER = "CCTGTACTTCGTTCAGTTACGTATTGCT"
RB11_BLOCK = "GCTTGGGTGTTTAACC" + "GTTTCATCTATCGGAGGGAATGGA" + screen.MU_RIGHT_END


def write(path, text):
    with open(path, "w") as out:
        out.write(text)
    return str(path)


def genome_file(tmp_path, lengths):
    # Contigs of the given lengths, made of a repeated 4 bp unit (the sequence is never read
    # by the report beyond its length, except by prepare).
    text = "".join(f">{name}\n{'ACGT' * (length // 4)}\n" for name, length in lengths.items())
    return write(tmp_path / "contigs.fasta", text)


def replicons_file(tmp_path, topology):
    # The Bakta replicon table: contig, new name, type, topology, name.
    text = "".join(f"{name}\t-\tcontig\t{kind}\t-\n" for name, kind in topology.items())
    return write(tmp_path / "replicons.tsv", text)


def query_file(tmp_path):
    return write(tmp_path / "query.fasta",
                 f">ONT_adapter_top adapter\n{ADAPTER}\n"
                 f">ONT_rapid_barcode_RB11 rapid barcode RB11\n{RB11_BLOCK}\n"
                 f">ONT_rapid_barcode_3prime_flank Mu end\n{screen.MU_RIGHT_END}\n")


def blast_file(tmp_path, rows):
    # rows: (query, contig, identity, length, qstart, qend, sstart, send, bitscore, qlen, slen)
    lines = []
    for query, contig, identity, length, qstart, qend, sstart, send, bitscore, qlen, slen in rows:
        lines.append("\t".join(str(value) for value in (
            query, contig, identity, length, qstart, qend, sstart, send, "1e-5", bitscore, qlen, slen)))
    return write(tmp_path / "hits.tsv", "\n".join(lines) + "\n")


def run_report(tmp_path, lengths, rows, topology=None, choice="", check=""):
    contigs = genome_file(tmp_path, lengths)
    replicons = replicons_file(tmp_path, topology) if topology else ""
    out = str(tmp_path / "loci.tsv")
    mqc = str(tmp_path / "mqc.tsv")
    argv = ["report", "--blast", blast_file(tmp_path, rows), "--contigs", contigs,
            "--query", query_file(tmp_path), "--sample", "S1", "--out", out, "--out-mqc", mqc]
    argv += ["--replicons", replicons] if replicons else ["--replicons"]
    argv += ["--read-choice", choice] if choice else ["--read-choice"]
    argv += ["--read-check", check] if check else ["--read-check"]
    screen.main(argv)
    with open(out) as handle:
        lines = handle.read().strip().split("\n")
    header = lines[0].split("\t")
    loci = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    with open(mqc) as handle:
        mqc_lines = [line for line in handle.read().strip().split("\n") if not line.startswith("#")]
    mqc_row = dict(zip(mqc_lines[0].split("\t"), mqc_lines[1].split("\t")))
    return loci, mqc_row


# ── prepare ──

def test_prepare_joins_the_lists_and_extends_only_circular_contigs(tmp_path):
    list_one = write(tmp_path / "one.fasta", f">ONT_adapter_top adapter\n{ADAPTER}\n")
    list_two = write(tmp_path / "two.fasta", ">UniVec_X pacbio\nACGTACGTACGTACGTACGTAC\n")
    trimmed = write(tmp_path / "trimmed.fasta", f">S1_trim_1 rapid barcode RB11\n{RB11_BLOCK}\n")
    contigs = write(tmp_path / "contigs.fasta", ">chromosome\n" + "A" * 300 + "C" * 300 + "\n>linear_one\n" + "G" * 500 + "\n")
    replicons = replicons_file(tmp_path, {"chromosome": "circular", "linear_one": "linear"})
    query, subject = str(tmp_path / "query.fasta"), str(tmp_path / "subject.fasta")
    screen.main(["prepare", "--contigs", contigs, "--replicons", replicons,
                 "--lists", list_one, list_two, "--trimmed", trimmed,
                 "--out-query", query, "--out-subject", subject])
    queries = screen.read_fasta(query)
    assert list(queries) == ["ONT_adapter_top", "UniVec_X", "S1_trim_1"]
    sequences = {name: sequence for name, (header, sequence) in screen.read_fasta(subject).items()}
    assert sequences["chromosome"] == "A" * 300 + "C" * 300 + "A" * 200
    assert sequences["linear_one"] == "G" * 500


def test_prepare_without_replicons_extends_nothing(tmp_path):
    list_one = write(tmp_path / "one.fasta", f">ONT_adapter_top adapter\n{ADAPTER}\n")
    contigs = write(tmp_path / "contigs.fasta", ">NODE_1\n" + "A" * 400 + "\n")
    query, subject = str(tmp_path / "query.fasta"), str(tmp_path / "subject.fasta")
    screen.main(["prepare", "--contigs", contigs, "--replicons", "--lists", list_one,
                 "--trimmed", "--out-query", query, "--out-subject", subject])
    sequences = {name: sequence for name, (header, sequence) in screen.read_fasta(subject).items()}
    assert sequences["NODE_1"] == "A" * 400


# ── report: the length and position rules ──

def test_an_end_hit_on_a_linear_contig_counts_from_20_bp(tmp_path):
    loci, mqc = run_report(tmp_path, {"contig_1": 10000},
                           [("ONT_adapter_top", "contig_1", 100.0, 22, 1, 22, 9960, 9981, 44.1, 28, 10000)])
    assert len(loci) == 1
    assert loci[0]["position"] == "end"
    assert loci[0]["distance_to_end"] == "19"
    assert loci[0]["counted"] == "yes"
    assert loci[0]["reason"] == "adapter_at_contig_end"
    assert mqc["Contig ends with adapter"] == "1"


def test_internal_hits_need_28_bp_and_90_percent(tmp_path):
    loci, _ = run_report(tmp_path, {"contig_1": 10000}, [
        ("ONT_adapter_top", "contig_1", 100.0, 28, 1, 28, 5001, 5028, 56.0, 28, 10000),   # kept
        ("ONT_adapter_top", "contig_1", 100.0, 25, 1, 25, 7001, 7025, 50.0, 28, 10000),   # too short inside
        ("ONT_adapter_top", "contig_1", 85.0, 28, 1, 28, 3001, 3028, 30.0, 28, 10000),    # too different
    ])
    assert [(locus["start"], locus["reason"]) for locus in loci] == [("5001", "adapter_inside_contig")]


def test_a_hit_near_the_end_of_a_circular_contig_is_internal(tmp_path):
    # 22 bp near the end: an end hit on a linear contig, but too short for the internal rule
    # on a circular one; 30 bp near the end: an internal locus, with no distance to an end.
    loci, mqc = run_report(tmp_path, {"chromosome": 10000}, [
        ("ONT_adapter_top", "chromosome", 100.0, 22, 1, 22, 9960, 9981, 44.1, 28, 10200),
        ("ONT_rapid_barcode_RB11", "chromosome", 100.0, 30, 1, 30, 9900, 9929, 60.0, 90, 10200),
    ], topology={"chromosome": "circular"})
    assert len(loci) == 1
    assert loci[0]["position"] == "internal"
    assert loci[0]["distance_to_end"] == "NA"
    assert loci[0]["topology"] == "circular"
    assert mqc["Contig ends with adapter"] == "0"
    assert mqc["Adapter matches inside contigs"] == "1"


def test_a_hit_across_the_start_of_a_circular_contig_is_one_locus(tmp_path):
    # The search copy is 10,000 + 200 bp; a 60 bp hit from 9,971 to 10,030 crosses the start.
    # BLAST also reports the second half of that adapter where it really sits, at 1-30: it is
    # the same adapter and must not become a second locus. A hit wholly in the appended copy
    # (10,101-10,128) repeats one at 101-128 and is dropped; the one at 101-128 is kept.
    loci, mqc = run_report(tmp_path, {"chromosome": 10000}, [
        ("ONT_rapid_barcode_RB11", "chromosome", 100.0, 60, 1, 60, 9971, 10030, 120.0, 90, 10200),
        ("ONT_rapid_barcode_RB11", "chromosome", 100.0, 30, 31, 60, 1, 30, 60.0, 90, 10200),
        ("ONT_adapter_top", "chromosome", 100.0, 28, 1, 28, 101, 128, 56.0, 28, 10200),
        ("ONT_adapter_top", "chromosome", 100.0, 28, 1, 28, 10101, 10128, 56.0, 28, 10200),
    ], topology={"chromosome": "circular"})
    starts_and_ends = sorted((int(locus["start"]), int(locus["end"])) for locus in loci)
    assert starts_and_ends == [(101, 128), (9971, 30)]
    assert mqc["Adapter matches inside contigs"] == "2"


# ── report: the Mu end ──

def test_a_mu_end_inside_a_contig_is_listed_but_not_counted(tmp_path, capsys):
    # A 50 bp match to the Mu end alone, deep inside the chromosome: a Mu prophage end.
    loci, mqc = run_report(tmp_path, {"chromosome": 50000}, [
        ("ONT_rapid_barcode_RB11", "chromosome", 100.0, 50, 41, 90, 20001, 20050, 99.0, 90, 50000),
        ("ONT_rapid_barcode_3prime_flank", "chromosome", 100.0, 50, 1, 50, 20001, 20050, 99.0, 50, 50000),
    ])
    assert len(loci) == 1
    assert loci[0]["counted"] == "no"
    assert loci[0]["reason"] == "mu_end_inside_contig"
    assert mqc["Adapter matches inside contigs"] == "0"
    assert "Mu end inside contigs, as a Mu prophage holds (not counted)" in mqc["Verdict"]
    assert "none of the listed adapters found" in capsys.readouterr().out


def test_a_mu_end_at_a_linear_contig_end_counts(tmp_path):
    loci, mqc = run_report(tmp_path, {"contig_2": 16144}, [
        ("ONT_rapid_barcode_3prime_flank", "contig_2", 100.0, 50, 1, 50, 16033, 16082, 99.0, 50, 16144),
    ])
    assert loci[0]["reason"] == "mu_end_at_contig_end"
    assert loci[0]["counted"] == "yes"
    assert mqc["Contig ends with adapter"] == "1"


def test_a_locus_that_covers_the_barcode_is_an_adapter_not_a_mu_end(tmp_path):
    # As at the end of TUM24772's contig_2: the whole 88 bp of the RB11 block, and the Mu end
    # alone from the flank record, on the same stretch.
    loci, _ = run_report(tmp_path, {"contig_2": 16144}, [
        ("ONT_rapid_barcode_3prime_flank", "contig_2", 100.0, 50, 1, 50, 16033, 16082, 99.0, 50, 16144),
        ("ONT_rapid_barcode_RB11", "contig_2", 100.0, 88, 3, 90, 16033, 16120, 170.0, 90, 16144),
    ])
    assert len(loci) == 1
    assert loci[0]["reason"] == "adapter_at_contig_end"
    assert loci[0]["longest_hit"] == "88"


# ── report: names, counting, the MultiQC row ──

def test_adapters_are_listed_best_hit_first(tmp_path):
    loci, _ = run_report(tmp_path, {"contig_2": 16144}, [
        ("ONT_rapid_barcode_3prime_flank", "contig_2", 100.0, 50, 1, 50, 16033, 16082, 99.0, 50, 16144),
        ("ONT_rapid_barcode_RB11", "contig_2", 100.0, 88, 3, 90, 16033, 16120, 170.0, 90, 16144),
    ])
    assert loci[0]["adapters"] == "ONT_rapid_barcode_RB11,ONT_rapid_barcode_3prime_flank"


def test_a_contig_end_is_counted_once(tmp_path):
    # Two separate loci within 100 bp of the same end, and one near the other end.
    loci, mqc = run_report(tmp_path, {"contig_1": 10000}, [
        ("ONT_adapter_top", "contig_1", 100.0, 28, 1, 28, 9930, 9957, 56.0, 28, 10000),
        ("ONT_adapter_top", "contig_1", 100.0, 28, 1, 28, 9970, 9997, 56.0, 28, 10000),
        ("ONT_adapter_top", "contig_1", 100.0, 28, 1, 28, 5, 32, 56.0, 28, 10000),
    ])
    assert len(loci) == 3
    assert mqc["Contig ends with adapter"] == "2"


def test_the_mqc_row_carries_the_read_numbers_of_a_trimmed_sample(tmp_path):
    choice = write(tmp_path / "choice.tsv",
                   "percent_with_adapter\t0.00\npercent_with_barcode_without_adapter\t19.89\n"
                   "published_barcodes\tRB11 19.9 % of all reads\n"
                   "second_barcode\tnone\n"
                   "summary\trapid barcode RB11 without the adapter core in front (90 bp)\n")
    check = write(tmp_path / "check.tsv", "percent_with_adapter\t0.00\npercent_with_end_of_trimmed_sequence\t0.12\n")
    _, mqc = run_report(tmp_path, {"contig_1": 10000}, [], choice=choice, check=check)
    assert mqc["Reads with barcode, no adapter (%)"] == "19.89"
    assert mqc["ONT barcodes in the reads"] == "RB11 19.9 % of all reads"
    assert mqc["Trimmed from the reads"] == "rapid barcode RB11 without the adapter core in front (90 bp)"
    assert mqc["Barcode end left after trimming (%)"] == "0.12"
    assert mqc["Verdict"] == "none of the listed adapters found"


def test_the_mqc_row_says_na_for_a_sample_bacflux_did_not_trim(tmp_path):
    _, mqc = run_report(tmp_path, {"NODE_1": 5000}, [])
    assert mqc["Sample"] == "S1"
    assert mqc["Reads with adapter (%)"] == "NA"
    assert mqc["Trimmed from the reads"] == "reads not trimmed by BacFlux"


def test_an_empty_blast_table_gives_a_header_only_locus_table(tmp_path):
    loci, mqc = run_report(tmp_path, {"contig_1": 1000}, [])
    assert loci == []
    assert mqc["Contig ends with adapter"] == "0"

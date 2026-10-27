"""Unit tests for workflow/scripts/80_mobilome/isescan_to_table.py.

isescan_to_table turns one sample's raw ISEScan output into the three files the
rest of the module depends on: the IS table colocalise.py reads, the per-sample
QC summary that table has to be read alongside, and the audit file naming every
record dropped or flagged. These tests cover all three, plus the second CLI mode
(rule contig_lengths) that measures the assembly FASTA, and the two kinds of call
dropped for a biological reason: an IS call lying over the essential gene dnaA,
and a partial IS21 lying over the chromosome's own helicase loader, dnaC or dnaI.

No tools, no databases, no real ISEScan run: every test builds a small
ISEScan-shaped .tsv by hand and pushes it through the parser, so the whole
normalisation and QC logic is verifiable in a plain Python environment.

The fixtures follow the real ISEScan 1.7.3 output contract (24 named columns,
1-based inclusive isBegin/isEnd, 'c'/'p' completeness letter, irLen = 0 when no
terminal inverted repeat was found), so a test passing here means the parser
agrees with the tool, not just with itself.

Two facts run through most of what follows. An IS count is a floor, not a count —
identical copies collapse in a short-read assembly and never appear — which is
why the distance to the contig end travels with every row. And an empty table
has two completely different meanings: ISEScan looked and found nothing, or
ISEScan wrote nothing at all. The audit file has to tell them apart.

Run:
    python -m pytest workflow/scripts/80_mobilome/test_isescan_to_table.py -q
"""

import csv
import os
import sys

# Import the script under test regardless of where pytest was started from.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import isescan_to_table as it  # noqa: E402


# ── Fixture builders ─────────────────────────────────────────────────────────

# One realistic ISEScan record, used as the base for every test row. Values are
# internally consistent: a 1201 bp IS3-family element at 5000..6200 on contig_1,
# with 25 bp terminal inverted repeats at each end and a transposase ORF inside.
DEFAULT_FIELDS = {
    "seqID": "contig_1",
    "family": "IS3",
    "cluster": "IS3_1",
    "isBegin": "5000",
    "isEnd": "6200",
    "isLen": "1201",
    "ncopy4is": "2",
    "start1": "5000",
    "end1": "5024",
    "start2": "6176",
    "end2": "6200",
    "score": "25",
    "irId": "23",
    "irLen": "25",
    "nGaps": "1",
    "orfBegin": "5100",
    "orfEnd": "6100",
    "strand": "+",
    "orfLen": "1001",
    "E-value": "1e-30",
    "E-value4copy": "1e-30",
    "type": "c",
    "ov": "2",
    "tir": "AAGGCCTTAA:TTAAGGCCTT",
}


def isescan_row(**overrides):
    """Build one tab-separated ISEScan data line, in the tool's real column order."""
    fields = dict(DEFAULT_FIELDS)
    fields.update({key: str(value) for key, value in overrides.items()})
    return "\t".join(fields[column] for column in it.ISESCAN_COLUMNS)


def isescan_header(columns=None):
    """The header line ISEScan writes: its 24 column names, tab-separated."""
    return "\t".join(columns if columns is not None else it.ISESCAN_COLUMNS)


def write_isescan_tsv(path, rows, columns=None):
    """Write a complete ISEScan results file (header + data lines)."""
    lines = [isescan_header(columns)] + list(rows)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return str(path)


def write_lengths(path, lengths):
    """Write the two-column contig-length TSV the main mode reads."""
    it.write_contig_lengths(str(path), lengths)
    return str(path)


def read_tsv(path):
    """Read a written TSV back as (header list, list of row dicts)."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader if row]
    header = rows[0]
    return header, [dict(zip(header, row)) for row in rows[1:]]


def run_main(tmp_path, isescan_out, lengths, boundary_bp=100, extra=None, sample="S1"):
    """Run the script's main() the way the Snakemake rule will, and read back all
    three outputs. Returns (return_code, table_rows, summary_row, audit_rows).

    tmp_path is where the inputs and outputs of THIS invocation live; tests that
    run main() twice on the same fixture pass a sub-directory each time, so make
    sure it exists first.
    """
    os.makedirs(str(tmp_path), exist_ok=True)
    lengths_file = write_lengths(tmp_path / "contig_lengths.tsv", lengths)
    out_table = str(tmp_path / "is_elements.tsv")
    out_summary = str(tmp_path / "is_summary.tsv")
    out_audit = str(tmp_path / "is_discarded.tsv")

    argv = [
        "--sample", sample,
        "--isescan-out", str(isescan_out),
        "--contig-lengths", lengths_file,
        "--boundary-bp", str(boundary_bp),
        "--out-table", out_table,
        "--out-summary", out_summary,
        "--out-audit", out_audit,
    ]
    if extra:
        argv.extend(extra)

    return_code = it.main(argv)
    if return_code != 0:
        return return_code, None, None, None

    _, table_rows = read_tsv(out_table)
    _, summary_rows = read_tsv(out_summary)
    _, audit_rows = read_tsv(out_audit)
    return return_code, table_rows, summary_rows[0], audit_rows


# ── The contig-length helper ─────────────────────────────────────────────────

def test_fasta_contig_lengths_wrapped_sequence_and_described_header(tmp_path):
    # Headers carry a description; the contig ID is only the first token, because
    # that is what ISEScan puts in seqID. Sequence is wrapped over several lines.
    fasta = tmp_path / "contigs_final.fasta"
    fasta.write_text(
        ">contig_1 length=20 cov=33.1\n"
        "ACGTACGTAC\n"
        "GTACGTACGT\n"
        ">contig_2\n"
        "ACGTA\n"
    )
    lengths = it.fasta_contig_lengths(str(fasta))
    assert lengths == {"contig_1": 20, "contig_2": 5}


def test_contig_lengths_roundtrip_with_and_without_header(tmp_path):
    # write_contig_lengths writes a header; read_contig_lengths must cope with it
    # and with a hand-made headerless file (a user doing this with awk is likely).
    with_header = tmp_path / "with_header.tsv"
    it.write_contig_lengths(str(with_header), {"contig_2": 500, "contig_1": 1000})
    assert it.read_contig_lengths(str(with_header)) == {"contig_1": 1000, "contig_2": 500}
    # Sorted output keeps the file identical between runs.
    assert with_header.read_text().splitlines()[1].startswith("contig_1")

    headerless = tmp_path / "headerless.tsv"
    headerless.write_text("contig_1\t1000\ncontig_2\t500\n")
    assert it.read_contig_lengths(str(headerless)) == {"contig_1": 1000, "contig_2": 500}


def test_helper_cli_mode_writes_contig_lengths(tmp_path):
    """The second CLI mode, run by rule contig_lengths before ISEScan is parsed.

    It measures the same assembly FASTA ISEScan was given, so every
    distance-to-contig-end in the IS table is computed against the contigs the
    calls actually came off.
    """
    fasta = tmp_path / "contigs_final.fasta"
    fasta.write_text(">contig_1\nACGTACGTAC\n>contig_2\nACGT\n")
    out = tmp_path / "contig_lengths.tsv"
    return_code = it.main(["--genome-fasta", str(fasta), "--out-contig-lengths", str(out)])
    assert return_code == 0
    header, rows = read_tsv(out)
    assert header == ["contig", "length"]
    assert rows == [{"contig": "contig_1", "length": "10"},
                    {"contig": "contig_2", "length": "4"}]


# ── Normal parse ─────────────────────────────────────────────────────────────

def test_normal_multi_is_parse(tmp_path):
    # Three IS on two contigs, deliberately written out of positional order to
    # prove the output is sorted by contig then start.
    results = write_isescan_tsv(tmp_path / "contigs_final.fasta.tsv", [
        isescan_row(seqID="contig_2", family="IS200/IS605", isBegin="3000", isEnd="3800"),
        isescan_row(seqID="contig_1", family="IS3", isBegin="5000", isEnd="6200"),
        isescan_row(seqID="contig_1", family="IS6", cluster="IS6_2",
                    isBegin="12000", isEnd="12820", strand="-"),
    ])
    code, table, summary, audit = run_main(
        tmp_path, results, {"contig_1": 50000, "contig_2": 20000}
    )
    assert code == 0
    assert audit == []
    assert [row["contig"] for row in table] == ["contig_1", "contig_1", "contig_2"]
    assert [row["start"] for row in table] == ["5000", "12000", "3000"]

    first = table[0]
    # Coordinates are passed through unchanged (1-based inclusive), and the length
    # counts both end bases.
    assert first["start"] == "5000" and first["end"] == "6200"
    assert first["length_bp"] == "1201"
    assert first["strand"] == "+"
    assert first["is_family"] == "IS3"
    assert first["cluster"] == "IS3_1"
    assert first["mge_type"] == "insertion_sequence"
    # The module's ID convention: contig|mge_type-start:end.
    assert first["mge_id"] == "contig_1|insertion_sequence-5000:6200"
    assert first["sample"] == "S1"
    assert first["contig_length"] == "50000"
    # Pass-through detail columns.
    assert first["ir_len_bp"] == "25"
    assert first["ir_identical_bp"] == "23"
    assert first["ir_gaps"] == "1"
    assert first["n_copies_is"] == "2"
    assert first["evalue"] == "1e-30"
    assert first["tpase_orf_start"] == "5100"
    assert first["tpase_orf_end"] == "6100"
    assert table[1]["strand"] == "-"
    # The written header is exactly the declared schema, in order.
    assert list(first.keys()) == it.OUTPUT_COLUMNS


def test_empty_strand_becomes_na(tmp_path):
    # ISEScan leaves strand empty for an IS with no predicted transposase.
    results = write_isescan_tsv(tmp_path / "r.tsv", [isescan_row(strand="")])
    code, table, _, _ = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert table[0]["strand"] == "NA"


# ── Complete vs partial ──────────────────────────────────────────────────────

def test_complete_and_partial_flagging(tmp_path):
    results = write_isescan_tsv(tmp_path / "r.tsv", [
        isescan_row(isBegin="5000", isEnd="6200", type="c"),
        isescan_row(isBegin="20000", isEnd="20350", type="p", irLen="0", tir=""),
    ])
    code, table, summary, _ = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert table[0]["is_complete"] == "TRUE"
    assert table[0]["is_type"] == "c"
    assert table[1]["is_complete"] == "FALSE"
    assert table[1]["is_type"] == "p"
    assert summary["n_is_complete"] == "1"
    assert summary["n_is_partial"] == "1"
    assert summary["n_is_completeness_unknown"] == "0"


def test_unknown_type_letter_is_kept_and_flagged(tmp_path):
    # A completeness letter we do not recognise must not cost us a real IS: the
    # element is kept with is_complete=NA and an audit line explains it.
    results = write_isescan_tsv(tmp_path / "r.tsv", [isescan_row(type="x")])
    code, table, summary, audit = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert len(table) == 1
    assert table[0]["is_complete"] == "NA"
    assert table[0]["is_type"] == "x"
    assert summary["n_is_completeness_unknown"] == "1"
    assert len(audit) == 1
    assert audit[0]["action"] == "kept_flagged"
    assert audit[0]["reason"] == "unknown_isescan_type_value"
    assert summary["n_records_flagged"] == "1"
    assert summary["n_records_dropped"] == "0"


def test_ir_present_reflects_ir_length(tmp_path):
    results = write_isescan_tsv(tmp_path / "r.tsv", [
        isescan_row(isBegin="5000", isEnd="6200", irLen="25"),
        isescan_row(isBegin="20000", isEnd="21200", irLen="0", tir=""),
    ])
    code, table, summary, _ = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert table[0]["ir_present"] == "TRUE"
    assert table[1]["ir_present"] == "FALSE"
    assert summary["n_is_with_ir"] == "1"


# ── The short-read honesty signals ───────────────────────────────────────────
#
# An IS at a contig end sits where the assembly fell apart, so whatever was next
# to it is unknown. These columns are the whole reason the IS table can be read
# honestly: colocalise.py caps a gene's confidence when its context element
# touches an edge, and the per-sample fraction says how much of the inventory is
# in that state. The distance is always reported next to the flag, because only
# the flag depends on --boundary-bp.

def test_is_flush_against_contig_start_is_at_boundary(tmp_path):
    # An IS starting at base 1: zero bases of contig before it. This is exactly
    # where a short-read assembly breaks, so it must be flagged.
    results = write_isescan_tsv(tmp_path / "r.tsv", [isescan_row(isBegin="1", isEnd="900")])
    code, table, _, _ = run_main(tmp_path, results, {"contig_1": 10000})
    assert code == 0
    assert table[0]["dist_to_contig_start"] == "0"
    assert table[0]["dist_to_contig_end"] == "9100"
    assert table[0]["dist_to_nearest_contig_end"] == "0"
    assert table[0]["at_contig_boundary"] == "TRUE"


def test_is_flush_against_contig_end_is_at_boundary(tmp_path):
    # The mirror case: the IS runs to the last base of the contig.
    results = write_isescan_tsv(tmp_path / "r.tsv",
                                [isescan_row(isBegin="9101", isEnd="10000")])
    code, table, _, _ = run_main(tmp_path, results, {"contig_1": 10000})
    assert code == 0
    assert table[0]["dist_to_contig_start"] == "9100"
    assert table[0]["dist_to_contig_end"] == "0"
    assert table[0]["dist_to_nearest_contig_end"] == "0"
    assert table[0]["at_contig_boundary"] == "TRUE"


def test_is_in_the_middle_is_not_at_boundary(tmp_path):
    results = write_isescan_tsv(tmp_path / "r.tsv",
                                [isescan_row(isBegin="5000", isEnd="6200")])
    code, table, _, _ = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert table[0]["dist_to_nearest_contig_end"] == "4999"
    assert table[0]["at_contig_boundary"] == "FALSE"


def test_boundary_threshold_is_configurable(tmp_path):
    # 150 bp of contig before the IS: outside the default 100 bp window, inside a
    # 200 bp one. The distance itself never changes — only the flag does, which is
    # why the distance is reported alongside it.
    results = write_isescan_tsv(tmp_path / "r.tsv",
                                [isescan_row(isBegin="151", isEnd="1200")])
    lengths = {"contig_1": 50000}

    code, table, _, _ = run_main(tmp_path / "default", results, lengths, boundary_bp=100)
    assert code == 0
    assert table[0]["dist_to_nearest_contig_end"] == "150"
    assert table[0]["at_contig_boundary"] == "FALSE"

    code, table, summary, _ = run_main(tmp_path / "wide", results, lengths, boundary_bp=200)
    assert code == 0
    assert table[0]["dist_to_nearest_contig_end"] == "150"
    assert table[0]["at_contig_boundary"] == "TRUE"
    assert summary["boundary_bp"] == "200"


def test_qc_summary_counts_and_fraction(tmp_path):
    # Four IS on two contigs: two of them sit at a contig end, one is partial.
    # The fraction is the number the report must quote next to the total, because
    # the total is a floor.
    results = write_isescan_tsv(tmp_path / "r.tsv", [
        isescan_row(seqID="contig_1", isBegin="1", isEnd="900", type="c"),        # boundary
        isescan_row(seqID="contig_1", isBegin="5000", isEnd="6200", type="c"),    # interior
        isescan_row(seqID="contig_2", isBegin="19200", isEnd="20000", type="p"),  # boundary
        isescan_row(seqID="contig_2", isBegin="8000", isEnd="9200", type="c"),    # interior
    ])
    code, table, summary, audit = run_main(
        tmp_path, results, {"contig_1": 50000, "contig_2": 20000}
    )
    assert code == 0
    assert summary["sample"] == "S1"
    assert summary["n_is_total"] == "4"
    assert summary["n_is_complete"] == "3"
    assert summary["n_is_partial"] == "1"
    assert summary["n_contigs_with_is"] == "2"
    assert summary["n_is_at_contig_boundary"] == "2"
    assert summary["fraction_at_contig_boundary"] == "0.5000"
    assert summary["boundary_bp"] == "100"
    assert summary["isescan_results_file"] == results
    assert summary["n_records_dropped"] == "0"
    assert audit == []


# ── A genome with no IS, and other graceful-degradation paths ────────────────

def test_no_is_at_all_gives_empty_but_wellformed_output(tmp_path):
    # ISEScan writes NO files when it finds no IS. That is a normal result for a
    # genome, so every output must still exist, with its header, and exit 0.
    empty_isescan_dir = tmp_path / "isescan"
    empty_isescan_dir.mkdir()

    code, table, summary, audit = run_main(tmp_path, empty_isescan_dir, {"contig_1": 50000})
    assert code == 0
    assert table == []
    # The IS table is empty, so the audit file must say WHY it is empty — an
    # empty audit file next to an empty table would leave the reader unable to
    # tell a genome with no IS from a run that produced nothing.
    assert len(audit) == 1
    assert audit[0]["reason"] == "isescan_wrote_no_results_file"
    assert summary["n_is_total"] == "0"
    assert summary["n_is_complete"] == "0"
    assert summary["n_is_partial"] == "0"
    assert summary["n_is_at_contig_boundary"] == "0"
    # A fraction of nothing is reported as 0, formatted like any other value so the
    # column stays numeric when several samples' summaries are concatenated.
    assert summary["fraction_at_contig_boundary"] == "0.0000"
    assert summary["isescan_results_file"] == "NONE"

    # Headers are present even with zero rows.
    header, _ = read_tsv(tmp_path / "is_elements.tsv")
    assert header == it.OUTPUT_COLUMNS
    header, _ = read_tsv(tmp_path / "is_discarded.tsv")
    assert header == it.AUDIT_COLUMNS


def test_missing_input_path_is_graceful(tmp_path):
    # The --isescan-out path does not exist at all (ISEScan wrote nothing).
    code, table, summary, _ = run_main(
        tmp_path, tmp_path / "does_not_exist", {"contig_1": 50000}
    )
    assert code == 0
    assert table == []
    assert summary["n_is_total"] == "0"
    assert summary["isescan_results_file"] == "NONE"


def test_header_only_results_file_is_graceful(tmp_path):
    results = write_isescan_tsv(tmp_path / "r.tsv", [])
    code, table, summary, _ = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert table == []
    assert summary["n_is_total"] == "0"
    # The file WAS found and parsed, so provenance records it — unlike the
    # no-file case, where it reads NONE.
    assert summary["isescan_results_file"] == results


def test_completely_empty_results_file_is_graceful(tmp_path):
    # A zero-byte file: no header either, so there is not even a schema to check.
    # Still a normal outcome rather than a parse error.
    empty = tmp_path / "r.tsv"
    empty.write_text("")
    code, table, summary, _ = run_main(tmp_path, empty, {"contig_1": 50000})
    assert code == 0
    assert table == []
    assert summary["n_is_total"] == "0"


# ── "The tool produced nothing" vs "the genome has nothing" ──────────────────
# An empty IS table has two completely different meanings for anyone judging a
# mobility call, so the audit file must always say which one applies
# (CONTRIBUTING.md: every such decision gets an audit row with a reason).

def test_missing_results_file_and_empty_results_file_get_different_reasons(tmp_path):
    # This is the whole point: the two situations must not look alike.
    # (a) ISEScan's output directory exists but holds no results file — what a
    #     genuinely IS-free genome looks like, and also what a crashed run looks
    #     like, so the audit says both.
    no_file_dir = tmp_path / "no_file" / "isescan"
    no_file_dir.mkdir(parents=True)
    code, table, summary, audit_no_file = run_main(
        tmp_path / "no_file", no_file_dir, {"contig_1": 50000}
    )
    assert code == 0
    assert table == []
    assert len(audit_no_file) == 1
    assert audit_no_file[0]["action"] == "input_missing"
    assert audit_no_file[0]["reason"] == "isescan_wrote_no_results_file"
    # It is a statement about the sample, not about one element.
    assert audit_no_file[0]["sample"] == "S1"
    assert audit_no_file[0]["contig"] == "NA"
    assert audit_no_file[0]["start"] == "NA"
    assert audit_no_file[0]["end"] == "NA"
    # The reader is warned not to quote zero without checking the ISEScan log.
    assert "log" in audit_no_file[0]["detail"]

    # (b) ISEScan DID write a results file, it simply lists no IS. Only here is
    #     the empty table a statement about the assembly.
    (tmp_path / "empty").mkdir()
    header_only = write_isescan_tsv(tmp_path / "empty" / "r.tsv", [])
    code, table, summary, audit_empty_file = run_main(
        tmp_path / "empty", header_only, {"contig_1": 50000}
    )
    assert code == 0
    assert table == []
    assert len(audit_empty_file) == 1
    assert audit_empty_file[0]["action"] == "input_empty"
    assert audit_empty_file[0]["reason"] == "isescan_results_file_has_no_rows"

    assert audit_no_file[0]["reason"] != audit_empty_file[0]["reason"]


def test_absent_isescan_path_gets_its_own_audit_reason(tmp_path):
    # A path that is not there at all is a third case, and a different kind of
    # problem: the isescan rule always creates its output directory, so this is
    # broken wiring rather than anything about the genome. The script still exits
    # 0 (nothing downstream should die over it), but it says so in the audit.
    code, table, summary, audit = run_main(
        tmp_path, tmp_path / "does_not_exist", {"contig_1": 50000}
    )
    assert code == 0
    assert table == []
    assert len(audit) == 1
    assert audit[0]["action"] == "input_missing"
    assert audit[0]["reason"] == "isescan_output_path_missing"
    assert summary["isescan_results_file"] == "NONE"


def test_sample_level_audit_row_is_not_counted_as_a_dropped_record(tmp_path):
    # Nothing was filtered out — there was simply nothing to filter — so the
    # dropped/flagged counters in the QC summary must stay at zero, otherwise the
    # summary would suggest IS calls had been thrown away.
    empty_isescan_dir = tmp_path / "isescan"
    empty_isescan_dir.mkdir()
    code, _, summary, audit = run_main(tmp_path, empty_isescan_dir, {"contig_1": 50000})
    assert code == 0
    assert len(audit) == 1
    assert summary["n_records_dropped"] == "0"
    assert summary["n_records_flagged"] == "0"
    # The audit file keeps its declared schema, sample-level row included.
    header, _ = read_tsv(tmp_path / "is_discarded.tsv")
    assert header == it.AUDIT_COLUMNS


def test_no_sample_level_row_when_isescan_reported_is(tmp_path):
    # The explanation fires only when there is nothing to report. With real IS in
    # the file the audit must carry per-record rows only — here, one dropped row —
    # so the sample-level note never dilutes the record-level audit trail.
    results = write_isescan_tsv(tmp_path / "r.tsv", [
        isescan_row(isBegin="5000", isEnd="6200"),
        isescan_row(isBegin="6200", isEnd="5000"),   # reversed → dropped
    ])
    code, table, summary, audit = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert len(table) == 1
    assert [row["reason"] for row in audit] == ["invalid_coordinate_range"]
    assert summary["n_records_dropped"] == "1"


def test_all_records_dropped_is_not_reported_as_no_results(tmp_path):
    # ISEScan DID report IS; we dropped them ourselves. That is already explained
    # by the per-record audit rows, so no sample-level "nothing to report" row is
    # added on top of them.
    results = write_isescan_tsv(tmp_path / "r.tsv", [
        isescan_row(seqID="contig_UNKNOWN", isBegin="100", isEnd="900"),
    ])
    code, table, _, audit = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert table == []
    assert [row["reason"] for row in audit] == ["contig_not_in_contig_lengths"]


# ── Audit: everything dropped or flagged is written down with a reason ───────

def test_malformed_coordinates_go_to_audit_without_crashing(tmp_path):
    # One good IS and one with a non-numeric isBegin. The good one must survive.
    results = write_isescan_tsv(tmp_path / "r.tsv", [
        isescan_row(isBegin="5000", isEnd="6200"),
        isescan_row(isBegin="not_a_number", isEnd="6200"),
    ])
    code, table, summary, audit = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert len(table) == 1
    assert table[0]["start"] == "5000"
    assert len(audit) == 1
    assert audit[0]["action"] == "dropped"
    assert audit[0]["reason"] == "unparseable_coordinates"
    assert audit[0]["sample"] == "S1"
    assert "not_a_number" in audit[0]["detail"]
    assert summary["n_records_dropped"] == "1"
    assert summary["n_is_total"] == "1"


def test_reversed_coordinates_go_to_audit(tmp_path):
    # isEnd before isBegin: the parser drops it rather than swapping the two
    # numbers, because a repaired row would look like a perfectly good IS at
    # coordinates nothing vouches for. Dropped with a reason, as every filtering
    # decision in this module is.
    results = write_isescan_tsv(tmp_path / "r.tsv",
                                [isescan_row(isBegin="6200", isEnd="5000")])
    code, table, _, audit = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert table == []
    assert audit[0]["reason"] == "invalid_coordinate_range"


def test_ragged_row_is_dropped_with_a_reason(tmp_path):
    # A truncated line means the column-to-value mapping cannot be trusted for
    # that line, so its coordinates might be silently wrong → drop and record.
    good = isescan_row(isBegin="5000", isEnd="6200")
    ragged = "\t".join(["contig_1", "IS3", "IS3_1", "7000"])
    results = write_isescan_tsv(tmp_path / "r.tsv", [good, ragged])
    code, table, summary, audit = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert len(table) == 1
    assert len(audit) == 1
    assert audit[0]["reason"] == "unexpected_field_count"
    assert summary["n_records_dropped"] == "1"


def test_contig_missing_from_lengths_file_is_dropped_with_a_reason(tmp_path):
    # Without the contig length there is no distance-to-contig-end, i.e. no
    # honesty signal, so the IS cannot be reported as if it had one.
    results = write_isescan_tsv(tmp_path / "r.tsv", [
        isescan_row(seqID="contig_1", isBegin="5000", isEnd="6200"),
        isescan_row(seqID="contig_UNKNOWN", isBegin="100", isEnd="900"),
    ])
    code, table, _, audit = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert [row["contig"] for row in table] == ["contig_1"]
    assert audit[0]["contig"] == "contig_UNKNOWN"
    assert audit[0]["reason"] == "contig_not_in_contig_lengths"


def test_coordinates_beyond_contig_length_are_dropped(tmp_path):
    # Means the contig-lengths file describes a different assembly from the one
    # ISEScan ran on — the boundary numbers would be nonsense.
    results = write_isescan_tsv(tmp_path / "r.tsv",
                                [isescan_row(isBegin="9000", isEnd="12000")])
    code, table, _, audit = run_main(tmp_path, results, {"contig_1": 10000})
    assert code == 0
    assert table == []
    assert audit[0]["reason"] == "coordinates_beyond_contig_length"


def test_min_length_filter_is_off_by_default_and_audits_when_used(tmp_path):
    # A 351 bp partial IS: kept by default (partials at contig breaks are exactly
    # what we want to see), dropped and audited only if the user asks for it.
    results = write_isescan_tsv(tmp_path / "r.tsv", [
        isescan_row(isBegin="5000", isEnd="6200", type="c"),
        isescan_row(isBegin="20000", isEnd="20350", type="p"),
    ])
    lengths = {"contig_1": 50000}

    code, table, _, audit = run_main(tmp_path / "default", results, lengths)
    assert code == 0
    assert len(table) == 2
    assert audit == []

    code, table, summary, audit = run_main(
        tmp_path / "filtered", results, lengths, extra=["--min-length-bp", "400"]
    )
    assert code == 0
    assert len(table) == 1
    assert audit[0]["reason"] == "below_min_length_bp"
    assert summary["min_length_bp"] == "400"


# ── A call lying over dnaA is not an insertion ───────────────────────────────
# dnaA is essential, so no living isolate has an IS inside it. ISEScan still
# reports a partial IS21 over the gene (its IS21 profile matches DnaA), and the
# script drops that call when Bakta's table says where dnaA is.

# The comment lines and the header Bakta writes at the top of {sample}.tsv.
BAKTA_PREAMBLE = [
    "# Annotated with Bakta",
    "# Software: v1.12.1",
    "# Database: v6.0, full",
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs",
]


def bakta_row(contig="contig_1", feature_type="cds", start=1, stop=1440,
              locus_tag="S1_00001", gene="dnaA",
              product="chromosomal replication initiator protein DnaA", strand="+"):
    """One data line of Bakta's annotation table; by default the dnaA CDS."""
    return "\t".join([contig, feature_type, str(start), str(stop), strand,
                      locus_tag, gene, product, "COG:COG0593"])


def write_bakta_tsv(path, rows):
    """Write a Bakta-shaped annotation table: comment lines, header, rows."""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(BAKTA_PREAMBLE + list(rows)) + "\n")
    return str(path)


def dnaa_call(**overrides):
    """The call ISEScan makes over dnaA: a partial IS21 lying over most of the
    gene (30..1600 against dnaA at 1..1440, so 1,411 of its 1,571 bp inside)."""
    fields = {"family": "IS21", "cluster": "IS21_1", "isBegin": "30",
              "isEnd": "1600", "type": "p"}
    fields.update(overrides)
    return isescan_row(**fields)


def test_call_over_dnaa_is_dropped_with_a_reason(tmp_path):
    results = write_isescan_tsv(tmp_path / "r.tsv", [dnaa_call(), isescan_row()])
    bakta = write_bakta_tsv(tmp_path / "bakta.tsv", [bakta_row()])
    code, table, summary, audit = run_main(
        tmp_path, results, {"contig_1": 50000}, extra=["--bakta-tsv", bakta]
    )
    assert code == 0
    assert [row["start"] for row in table] == ["5000"]
    assert audit[0]["action"] == "dropped"
    assert audit[0]["reason"] == "overlaps_dnaA"
    assert "1411 bp shared with dnaA" in audit[0]["detail"]
    assert "90% of the call, 98% of the gene" in audit[0]["detail"]
    assert "S1_00001" in audit[0]["detail"]
    # The call sat 29 bp from the contig start. Dropped, it no longer counts as
    # an IS, nor as an IS at a contig end.
    assert summary["n_is_total"] == "1"
    assert summary["n_is_at_contig_boundary"] == "0"
    assert summary["n_records_dropped"] == "1"


def test_call_running_past_dnaa_into_the_next_gene_is_dropped(tmp_path):
    # When ISEScan's gene finder joins dnaA and dnaN into one reading frame, the
    # call runs past the end of dnaA: 72..2949 against dnaA at 1..1404. Under
    # half of the call is inside dnaA, but almost all of dnaA is inside the call.
    results = write_isescan_tsv(tmp_path / "r.tsv",
                                [dnaa_call(isBegin="72", isEnd="2949")])
    bakta = write_bakta_tsv(tmp_path / "bakta.tsv", [
        bakta_row(start=1, stop=1404),
        bakta_row(start=1409, stop=2509, locus_tag="S1_00002", gene="dnaN",
                  product="DNA polymerase III subunit beta"),
    ])
    code, table, _, audit = run_main(
        tmp_path, results, {"contig_1": 50000}, extra=["--bakta-tsv", bakta]
    )
    assert code == 0
    assert table == []
    assert audit[0]["reason"] == "overlaps_dnaA"
    assert "46% of the call, 95% of the gene" in audit[0]["detail"]


def test_without_a_bakta_table_the_dnaa_call_is_kept(tmp_path):
    # The check needs Bakta's table. Without it nothing changes, so the script
    # still runs on ISEScan's output alone.
    results = write_isescan_tsv(tmp_path / "r.tsv", [dnaa_call()])
    code, table, _, audit = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert len(table) == 1
    assert audit == []


def test_call_next_to_dnaa_is_kept(tmp_path):
    # A real IS can sit beside dnaA. This one overlaps the gene's last 141 bp,
    # under a tenth of its 1,601 bp, so it stays in the table.
    results = write_isescan_tsv(tmp_path / "r.tsv",
                                [isescan_row(isBegin="1300", isEnd="2900")])
    bakta = write_bakta_tsv(tmp_path / "bakta.tsv", [bakta_row()])
    code, table, _, audit = run_main(
        tmp_path, results, {"contig_1": 50000}, extra=["--bakta-tsv", bakta]
    )
    assert code == 0
    assert len(table) == 1
    assert audit == []


def test_call_over_another_gene_is_kept(tmp_path):
    results = write_isescan_tsv(tmp_path / "r.tsv", [dnaa_call()])
    bakta = write_bakta_tsv(tmp_path / "bakta.tsv", [
        bakta_row(gene="dnaN", product="DNA polymerase III subunit beta"),
    ])
    code, table, _, audit = run_main(
        tmp_path, results, {"contig_1": 50000}, extra=["--bakta-tsv", bakta]
    )
    assert code == 0
    assert len(table) == 1
    assert audit == []


def test_dnaa_on_another_contig_does_not_drop_the_call(tmp_path):
    results = write_isescan_tsv(tmp_path / "r.tsv", [dnaa_call()])
    bakta = write_bakta_tsv(tmp_path / "bakta.tsv", [bakta_row(contig="contig_2")])
    code, table, _, _ = run_main(
        tmp_path, results, {"contig_1": 50000, "contig_2": 50000},
        extra=["--bakta-tsv", bakta],
    )
    assert code == 0
    assert len(table) == 1


def test_dnaa_is_found_by_its_product_when_the_symbol_is_missing(tmp_path):
    path = write_bakta_tsv(tmp_path / "bakta.tsv", [bakta_row(gene="")])
    genes = it.read_dnaa_genes(path)
    assert [(gene["contig"], gene["start"], gene["end"]) for gene in genes] == [
        ("contig_1", 1, 1440)
    ]


def test_dnaa_running_across_the_contig_end_is_skipped(tmp_path):
    # Start larger than Stop: the CDS runs across the end of a circular contig.
    # One interval cannot describe it, and skipping it can only keep a call.
    path = write_bakta_tsv(tmp_path / "bakta.tsv", [bakta_row(start=49500, stop=900)])
    assert it.read_dnaa_genes(path) == []


def test_missing_bakta_table_fails_loudly(tmp_path, capsys):
    results = write_isescan_tsv(tmp_path / "r.tsv", [isescan_row()])
    code, _, _, _ = run_main(
        tmp_path, results, {"contig_1": 50000},
        extra=["--bakta-tsv", str(tmp_path / "nope.tsv")],
    )
    assert code == 1
    assert "nope.tsv" in capsys.readouterr().err


def test_bakta_table_without_its_header_fails_loudly(tmp_path, capsys):
    # Finding no dnaA in a file we cannot read would quietly switch the check
    # off, so a table without Bakta's header stops the run instead.
    path = tmp_path / "bakta.tsv"
    path.write_text(bakta_row() + "\n", encoding="utf-8")
    results = write_isescan_tsv(tmp_path / "r.tsv", [isescan_row()])
    code, _, _, _ = run_main(
        tmp_path, results, {"contig_1": 50000}, extra=["--bakta-tsv", str(path)]
    )
    assert code == 1
    assert "Sequence Id" in capsys.readouterr().err


# ── A partial IS21 lying over the chromosome's helicase loader ───────────────
# Two of ISEScan's IS21 models match the replicative helicase loader, DnaC in
# E. coli and its relatives and DnaI in the Bacilli, which is essential too. The call is
# dropped only when the loader has its chromosomal partner beside it (dnaT
# "primosomal protein DnaT" for dnaC, nrdR for dnaI) and no IS21 transposase
# within 5 kb. Everything else over a loader is kept: phage, ICE and plasmid
# loaders, a real IS21 element's IstB that Bakta calls dnaC, and the loader in a
# taxon or a draft where the partner is not on the same contig.
#
# The layout below copies the E. coli arrangement: dnaC at 10001..10738
# (738 bp) with dnaT touching it at 10739..11278, and ISEScan's partial IS21 at
# 9800..11100, which holds all of dnaC: 738 bp shared, 57% of the 1,301 bp call
# and 100% of the gene.

def dnac_row(**overrides):
    """The chromosomal dnaC CDS of the layout above."""
    fields = {"start": 10001, "stop": 10738, "locus_tag": "S1_00950",
              "gene": "dnaC", "product": "DNA replication protein DnaC"}
    fields.update(overrides)
    return bakta_row(**fields)


def dnat_row(**overrides):
    """dnaT, touching dnaC: the gap between the two CDS is 0 bp."""
    fields = {"start": 10739, "stop": 11278, "locus_tag": "S1_00951",
              "gene": "dnaT", "product": "primosomal protein DnaT"}
    fields.update(overrides)
    return bakta_row(**fields)


def loader_call(**overrides):
    """ISEScan's partial IS21 call over the dnaC of the layout above."""
    fields = {"family": "IS21", "cluster": "IS21_259", "isBegin": "9800",
              "isEnd": "11100", "type": "p"}
    fields.update(overrides)
    return isescan_row(**fields)


def run_loader_case(tmp_path, bakta_rows, calls=None, lengths=None):
    """Run main() on the given Bakta rows and IS calls (by default the partial
    IS21 over dnaC plus the ordinary IS3 at 5000..6200 from DEFAULT_FIELDS, which
    must always survive), and return (table, summary, audit)."""
    if calls is None:
        calls = [loader_call(), isescan_row(isBegin="20000", isEnd="21200")]
    if lengths is None:
        lengths = {"contig_1": 50000}
    results = write_isescan_tsv(tmp_path / "r.tsv", calls)
    bakta = write_bakta_tsv(tmp_path / "bakta.tsv", bakta_rows)
    code, table, summary, audit = run_main(
        tmp_path, results, lengths, extra=["--bakta-tsv", bakta]
    )
    assert code == 0
    return table, summary, audit


def test_partial_is21_over_dnac_with_dnat_beside_it_is_dropped(tmp_path):
    table, summary, audit = run_loader_case(tmp_path, [dnac_row(), dnat_row()])
    # Only the ordinary IS survives.
    assert [row["start"] for row in table] == ["20000"]
    assert len(audit) == 1
    assert audit[0]["action"] == "dropped"
    assert audit[0]["reason"] == "overlaps_helicase_loader"
    assert (audit[0]["contig"], audit[0]["start"], audit[0]["end"]) == (
        "contig_1", "9800", "11100")
    assert summary["n_is_total"] == "1"
    assert summary["n_records_dropped"] == "1"


def test_helicase_loader_audit_detail_names_everything(tmp_path):
    # The detail must let a reader check the decision without opening Bakta's
    # table: the loader, its locus tag and coordinates, the partner and its
    # distance, the overlap, and why the call cannot be an insertion. It is
    # longer than the 200 characters other details are cut to, and kept whole.
    _, _, audit = run_loader_case(tmp_path, [dnac_row(), dnat_row()])
    detail = audit[0]["detail"]
    assert "738 bp shared with dnaC (S1_00950, 10001..10738)" in detail
    assert "57% of the call, 100% of the gene" in detail
    assert "partner dnaT 0 bp away" in detail
    assert detail.endswith(
        "The chromosome's own helicase loader is essential; this is a profile "
        "match, not an insertion")
    assert len(detail) > it.AUDIT_DETAIL_MAX_CHARS


def test_partial_is21_over_dnai_with_nrdr_upstream_is_dropped(tmp_path):
    # The Bacillus arrangement: nrdR, then dnaB, then dnaI, with 1,500 bp of
    # sequence between nrdR and dnaI. The call holds all 903 bp of dnaI.
    bakta_rows = [
        bakta_row(start=20001, stop=20456, locus_tag="S1_01700", gene="nrdR",
                  product="transcriptional regulator NrdR"),
        bakta_row(start=20500, stop=21900, locus_tag="S1_01701", gene="dnaB",
                  product="Replicative helicase loader DnaB"),
        bakta_row(start=21957, stop=22859, locus_tag="S1_01702", gene="dnaI",
                  product="primosomal protein DnaI"),
    ]
    calls = [isescan_row(family="IS21", type="p", isBegin="21530", isEnd="23271")]
    table, _, audit = run_loader_case(tmp_path, bakta_rows, calls=calls)
    assert table == []
    assert audit[0]["reason"] == "overlaps_helicase_loader"
    assert "903 bp shared with dnaI (S1_01702, 21957..22859)" in audit[0]["detail"]
    assert "partner nrdR 1500 bp away" in audit[0]["detail"]


def test_minus_strand_loader_and_partner_are_dropped_the_same_way(tmp_path):
    # Bakta writes Start < Stop on both strands, so the strand changes nothing.
    table, _, audit = run_loader_case(
        tmp_path, [dnac_row(strand="-"), dnat_row(strand="-")]
    )
    assert [row["start"] for row in table] == ["20000"]
    assert audit[0]["reason"] == "overlaps_helicase_loader"


def test_dnac_without_a_partner_is_kept(tmp_path):
    # A phage loader: Bakta calls it dnaC, but no dnaT or nrdR sits near it.
    table, _, audit = run_loader_case(tmp_path, [dnac_row(product="Phage related protein")])
    assert [row["start"] for row in table] == ["9800", "20000"]
    assert audit == []


def test_dnat_like_protein_of_a_mobile_element_does_not_count(tmp_path):
    # Mobile elements carry DnaT-like proteins that Bakta also names dnaT, but
    # with this product. They must not vouch for the loader beside them.
    table, _, audit = run_loader_case(tmp_path, [
        dnac_row(),
        dnat_row(product="DnaT DNA-binding domain-containing protein"),
    ])
    assert len(table) == 2
    assert audit == []


def test_partner_on_another_contig_does_not_count(tmp_path):
    # On a draft assembly the partner can fall on a different contig. Nothing
    # then shows that this loader is the chromosome's, so the call is kept.
    table, _, audit = run_loader_case(
        tmp_path,
        [dnac_row(), dnat_row(contig="contig_2", start=1, stop=540)],
        lengths={"contig_1": 50000, "contig_2": 3000},
    )
    assert len(table) == 2
    assert audit == []


# The two 5,000 bp windows carry the evidence behind the rule: nrdR sits up to
# about 3.2 kb from a chromosomal dnaI, and unrelated transposases sit 2.7-3.8 kb
# from several true chromosomal loaders. These tests pin both windows and their
# inclusive edge, so a later edit to either constant fails here. The dnaC of the
# layout ends at 10738, so a CDS starting at 10738 + 1 + gap lies `gap` bp away.

def partner_dnat_at_gap(gap_bp):
    """dnaT placed `gap_bp` bp after the layout's dnaC (0 would touch it)."""
    start = 10738 + 1 + gap_bp
    return dnat_row(start=start, stop=start + 539)


def test_partner_exactly_5000_bp_away_still_counts(tmp_path):
    table, _, audit = run_loader_case(tmp_path, [dnac_row(), partner_dnat_at_gap(5000)])
    assert [row["start"] for row in table] == ["20000"]
    assert audit[0]["reason"] == "overlaps_helicase_loader"


def test_partner_5001_bp_away_does_not_count(tmp_path):
    table, _, audit = run_loader_case(tmp_path, [dnac_row(), partner_dnat_at_gap(5001)])
    assert [row["start"] for row in table] == ["9800", "20000"]
    assert audit == []


def test_nrdr_3200_bp_upstream_of_dnai_is_its_partner(tmp_path):
    # The widest gap seen between nrdR and a chromosomal dnaI was about 3.2 kb.
    dnai = dnac_row(gene="dnaI", product="primosomal protein DnaI")
    nrdr = bakta_row(start=6300, stop=10001 - 1 - 3200, locus_tag="S1_00940",
                     gene="nrdR", product="transcriptional regulator NrdR")
    table, _, audit = run_loader_case(tmp_path, [nrdr, dnai])
    assert [row["start"] for row in table] == ["20000"]
    assert audit[0]["reason"] == "overlaps_helicase_loader"


def ista_at_gap(gap_bp):
    """An IS21 transposase placed `gap_bp` bp after the layout's dnaC."""
    start = 10738 + 1 + gap_bp
    return bakta_row(start=start, stop=start + 1199, locus_tag="S1_00960",
                     gene="istA", product="IS21 family transposase")


def test_transposase_exactly_5000_bp_away_keeps_the_call(tmp_path):
    table, _, audit = run_loader_case(tmp_path, [dnac_row(), dnat_row(), ista_at_gap(5000)])
    assert [row["start"] for row in table] == ["9800", "20000"]
    assert audit == []


def test_transposase_5001_bp_away_does_not_protect_the_call(tmp_path):
    table, _, audit = run_loader_case(tmp_path, [dnac_row(), dnat_row(), ista_at_gap(5001)])
    assert [row["start"] for row in table] == ["20000"]
    assert audit[0]["reason"] == "overlaps_helicase_loader"


def test_ista_near_the_loader_keeps_the_call(tmp_path):
    # An IS21 transposase 1,261 bp away: the "loader" may be a real IS21
    # element's IstB that Bakta happens to call dnaC.
    table, _, audit = run_loader_case(tmp_path, [
        dnac_row(), dnat_row(),
        bakta_row(start=12000, stop=13500, locus_tag="S1_00953", gene="istA",
                  product="IS21 family transposase"),
    ])
    assert len(table) == 2
    assert audit == []


def test_ista_named_only_in_its_product_keeps_the_call(tmp_path):
    # No istA symbol, and a wording that says neither "IstA" nor "transposase
    # IS21": the product's "IS21" is enough. It sits 1,500 bp away on the other
    # side of dnaC, with another CDS between them.
    table, _, audit = run_loader_case(tmp_path, [
        bakta_row(start=7000, stop=8500, locus_tag="S1_00948", gene="",
                  product="Transposase for insertion sequence element IS21-like "
                          "C-terminal domain-containing protein"),
        bakta_row(start=8600, stop=9900, locus_tag="S1_00949", gene="yjjB",
                  product="threonine/serine exporter"),
        dnac_row(), dnat_row(),
    ])
    assert len(table) == 2
    assert audit == []


def test_unrelated_transposase_near_the_loader_does_not_keep_the_call(tmp_path):
    # A bare "transposase" is not IS21 evidence: transposases of other families
    # do sit a few kb from true chromosomal loaders.
    table, _, audit = run_loader_case(tmp_path, [
        bakta_row(start=6000, stop=7000, locus_tag="S1_00947", gene="",
                  product="IS3 family transposase"),
        dnac_row(), dnat_row(),
    ])
    assert [row["start"] for row in table] == ["20000"]
    assert audit[0]["reason"] == "overlaps_helicase_loader"


def wrapped_ista_row(start, stop):
    """An istA that runs across the end of the 50 kb circular contig: Bakta
    writes it with Start larger than Stop, so it covers Start..50000 and then
    1..Stop."""
    return bakta_row(start=start, stop=stop, locus_tag="S1_04999", gene="istA",
                     product="IS21 family transposase")


def test_is21_transposase_across_the_contig_end_far_from_the_loader(tmp_path):
    # The istA covers 49500..50000 and 1..900. Its nearer piece, 1..900, ends
    # 9,100 bp before dnaC (10001 - 900 - 1), outside the 5 kb window, so it is
    # no evidence of an IS21 element here and the call is dropped as usual.
    table, _, audit = run_loader_case(tmp_path, [
        dnac_row(), dnat_row(), wrapped_ista_row(49500, 900),
    ])
    assert [row["start"] for row in table] == ["20000"]
    assert audit[0]["reason"] == "overlaps_helicase_loader"


def test_is21_transposase_across_the_contig_end_near_the_loader(tmp_path):
    # The same istA (49500..50000 and 1..900), but the loader layout moved to
    # the start of the contig: dnaC 2001..2738, dnaT touching it, the call at
    # 1800..3100. The piece 1..900 now ends 1,100 bp before dnaC, inside the
    # window, so the "loader" may be an IS21 IstB and the call is kept.
    calls = [loader_call(isBegin="1800", isEnd="3100"),
             isescan_row(isBegin="20000", isEnd="21200")]
    table, _, audit = run_loader_case(tmp_path, [
        dnac_row(start=2001, stop=2738), dnat_row(start=2739, stop=3278),
        wrapped_ista_row(49500, 900),
    ], calls=calls)
    assert len(table) == 2
    assert audit == []


def test_is21_transposase_across_the_contig_end_near_the_loader_from_behind(tmp_path):
    # Here the loader sits near the END of the contig (dnaC 46001..46738, dnaT
    # touching it, the call at 45800..47100) and the istA covers 48500..50000
    # and 1..400. The piece 48500..50000 starts 1,761 bp after dnaC
    # (48500 - 46738 - 1), inside the window, so the call is kept.
    calls = [loader_call(isBegin="45800", isEnd="47100"),
             isescan_row(isBegin="20000", isEnd="21200")]
    table, _, audit = run_loader_case(tmp_path, [
        dnac_row(start=46001, stop=46738), dnat_row(start=46739, stop=47278),
        wrapped_ista_row(48500, 400),
    ], calls=calls)
    assert len(table) == 2
    assert audit == []


def test_gap_to_cds_across_contig_end_measures_the_nearer_piece():
    # The wrapped CDS covers 48500..end and 1..400 (Start 48500, Stop 400).
    wrapped = {"start": 48500, "end": 400}
    # dnaC at 10001..10738: 9,600 bp after the piece 1..400 (10001 - 400 - 1)
    # and 37,761 bp before the piece from 48500 (48500 - 10738 - 1), so the
    # nearer piece gives 9,600.
    assert it._gap_to_cds_across_contig_end({"start": 10001, "end": 10738},
                                            wrapped) == 9600
    # Touching the piece 1..400 (anchor starts at 401): nothing between them.
    assert it._gap_to_cds_across_contig_end({"start": 401, "end": 900},
                                            wrapped) == 0
    # Overlapping either piece also gives 0.
    assert it._gap_to_cds_across_contig_end({"start": 300, "end": 900},
                                            wrapped) == 0
    assert it._gap_to_cds_across_contig_end({"start": 48000, "end": 48600},
                                            wrapped) == 0


def test_complete_is21_over_a_partnered_loader_is_kept(tmp_path):
    # A complete IS21 call is a whole element, not a profile match on a gene.
    table, _, audit = run_loader_case(
        tmp_path, [dnac_row(), dnat_row()],
        calls=[loader_call(type="c")],
    )
    assert len(table) == 1
    assert audit == []


def test_partial_call_of_another_family_over_a_partnered_loader_is_kept(tmp_path):
    # The loader match is an IS21 one; an IS3 call over dnaC is something else.
    table, _, audit = run_loader_case(
        tmp_path, [dnac_row(), dnat_row()],
        calls=[loader_call(family="IS3", cluster="IS3_1")],
    )
    assert len(table) == 1
    assert audit == []


def test_call_grazing_the_loader_is_kept(tmp_path):
    # 10638..12000 shares only dnaC's last 101 bp: under half of the 738 bp
    # gene, the shorter of the two, so it does not lie over the loader.
    table, _, audit = run_loader_case(
        tmp_path, [dnac_row(), dnat_row()],
        calls=[loader_call(isBegin="10638", isEnd="12000")],
    )
    assert len(table) == 1
    assert audit == []


def test_loader_is_recognised_by_its_symbol_only(tmp_path):
    # The same product without the dnaC symbol is not treated as a loader.
    table, _, audit = run_loader_case(tmp_path, [dnac_row(gene=""), dnat_row()])
    assert len(table) == 2
    assert audit == []


def test_bakta_table_without_loaders_drops_only_the_dnaa_call(tmp_path):
    # A genome with dnaA but no dnaC or dnaI: the dnaA check works as before
    # and nothing else is dropped, the partial IS21 elsewhere included.
    calls = [dnaa_call(), loader_call(), isescan_row(isBegin="20000", isEnd="21200")]
    table, _, audit = run_loader_case(
        tmp_path, [bakta_row(), dnat_row()], calls=calls
    )
    assert [row["start"] for row in table] == ["9800", "20000"]
    assert [row["reason"] for row in audit] == ["overlaps_dnaA"]


def test_dnaa_and_loader_calls_are_dropped_with_their_own_reasons(tmp_path):
    calls = [dnaa_call(), loader_call(), isescan_row(isBegin="20000", isEnd="21200")]
    table, summary, audit = run_loader_case(
        tmp_path, [bakta_row(), dnac_row(), dnat_row()], calls=calls
    )
    assert [row["start"] for row in table] == ["20000"]
    assert [row["reason"] for row in audit] == ["overlaps_dnaA", "overlaps_helicase_loader"]
    assert summary["n_records_dropped"] == "2"


def test_what_counts_as_an_is21_transposase():
    # Symbol istA, or a product with IS21, the whole word IstA, or the domain
    # names Bakta uses for IstA's transposase family. Ignoring case.
    counts = it._is_is21_transposase
    assert counts("istA", "hypothetical protein")
    assert counts("", "IS21 family transposase")
    assert counts("", "Truncated IS232 IstA -like protein")
    assert counts("", "Mu transposase, C-terminal domain-containing protein")
    assert counts("", "Integrase catalytic domain-containing protein")
    assert counts("", "integrase core domain-containing protein")
    # Not IS21 evidence: another family, a bare transposase, and a word that
    # merely contains the letters of IstA.
    assert not counts("", "IS3 family transposase")
    assert not counts("", "transposase")
    assert not counts("", "tellurite resistance protein")


def test_judge_helicase_loaders_says_why_each_loader_was_or_was_not_used(tmp_path):
    path = write_bakta_tsv(tmp_path / "bakta.tsv", [
        dnac_row(), dnat_row(),
        bakta_row(start=30001, stop=30738, locus_tag="S1_02000", gene="dnaC",
                  product="Phage related protein"),
    ])
    judged = it.judge_helicase_loaders(it.read_bakta_genes(path))
    assert [loader["is_chromosomal"] for loader in judged] == [True, False]
    assert judged[0]["partner_gene"] == "dnaT"
    assert judged[0]["partner_gap_bp"] == 0
    assert "no dnaT or nrdR partner" in judged[1]["decision"]


# ── Loud failure when the format is not recognisable ─────────────────────────

def test_missing_required_column_fails_loudly(tmp_path):
    # A header without isEnd: guessing would give wrong coordinates that look
    # perfectly plausible, so this must stop with a message naming file + column.
    columns = [c for c in it.ISESCAN_COLUMNS if c != "isEnd"]
    truncated_row = "\t".join(
        DEFAULT_FIELDS[column] for column in columns
    )
    results = write_isescan_tsv(tmp_path / "r.tsv", [truncated_row], columns=columns)

    try:
        it.read_isescan_table(results)
    except ValueError as error:
        message = str(error)
        assert "isEnd" in message
        assert results in message
    else:
        raise AssertionError("expected a ValueError for the missing isEnd column")


def test_missing_required_column_makes_main_exit_nonzero(tmp_path, capsys):
    columns = [c for c in it.ISESCAN_COLUMNS if c != "seqID"]
    row = "\t".join(DEFAULT_FIELDS[column] for column in columns)
    results = write_isescan_tsv(tmp_path / "r.tsv", [row], columns=columns)

    code, _, _, _ = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 1
    error_text = capsys.readouterr().err
    assert "seqID" in error_text
    assert results in error_text


def test_missing_contig_lengths_file_fails_loudly(tmp_path, capsys):
    results = write_isescan_tsv(tmp_path / "r.tsv", [isescan_row()])
    code = it.main([
        "--sample", "S1",
        "--isescan-out", results,
        "--contig-lengths", str(tmp_path / "nope.tsv"),
        "--out-table", str(tmp_path / "t.tsv"),
        "--out-summary", str(tmp_path / "s.tsv"),
        "--out-audit", str(tmp_path / "a.tsv"),
    ])
    assert code == 1
    assert "nope.tsv" in capsys.readouterr().err


def test_missing_required_arguments_fail_loudly(tmp_path, capsys):
    code = it.main(["--sample", "S1"])
    assert code == 1
    assert "--isescan-out" in capsys.readouterr().err


# ── Finding the results file inside ISEScan's output directory ───────────────

def test_directory_resolution_ignores_isescan_intermediates(tmp_path):
    # ISEScan mirrors the input path: results land in {output}/{sample}/ while its
    # own proteome/ and hmm/ trees sit alongside. Only the real results count.
    output_dir = tmp_path / "isescan"
    (output_dir / "S1").mkdir(parents=True)
    (output_dir / "proteome" / "S1").mkdir(parents=True)
    (output_dir / "hmm" / "S1").mkdir(parents=True)

    real = output_dir / "S1" / "contigs_final.fasta.tsv"
    write_isescan_tsv(real, [isescan_row()])
    (output_dir / "proteome" / "S1" / "decoy.tsv").write_text("junk\n")
    (output_dir / "hmm" / "S1" / "decoy.tsv").write_text("junk\n")

    assert it.find_isescan_results(str(output_dir)) == str(real)
    # A direct path to the file works too, so the rule can pass either.
    assert it.find_isescan_results(str(real)) == str(real)
    assert it.find_isescan_results(str(tmp_path / "absent")) is None


def test_end_to_end_through_the_output_directory(tmp_path):
    output_dir = tmp_path / "isescan"
    (output_dir / "S1").mkdir(parents=True)
    write_isescan_tsv(output_dir / "S1" / "contigs_final.fasta.tsv", [
        isescan_row(seqID="contig_1", isBegin="1", isEnd="900", type="p"),
        isescan_row(seqID="contig_1", isBegin="5000", isEnd="6200", type="c"),
    ])
    code, table, summary, audit = run_main(tmp_path, output_dir, {"contig_1": 50000})
    assert code == 0
    assert len(table) == 2
    assert summary["n_is_total"] == "2"
    assert summary["n_is_at_contig_boundary"] == "1"
    assert summary["fraction_at_contig_boundary"] == "0.5000"
    assert audit == []

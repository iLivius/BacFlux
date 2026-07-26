"""Unit tests for workflow/scripts/mobilome/isescan_to_table.py.

No tools, no databases, no real ISEScan run: every test builds a small
ISEScan-shaped .tsv by hand and pushes it through the parser, so the whole
normalisation and QC logic is verifiable in a plain Python environment.

The fixtures follow the real ISEScan 1.7.3 output contract (24 named columns,
1-based inclusive isBegin/isEnd, 'c'/'p' completeness letter, irLen = 0 when no
terminal inverted repeat was found), so a test passing here means the parser
agrees with the tool, not just with itself.

Run:
    python -m pytest workflow/scripts/mobilome/test_isescan_to_table.py -q
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
    # The file WAS found and parsed, so provenance records it (unlike the no-file case).
    assert summary["isescan_results_file"] == results


def test_completely_empty_results_file_is_graceful(tmp_path):
    empty = tmp_path / "r.tsv"
    empty.write_text("")
    code, table, summary, _ = run_main(tmp_path, empty, {"contig_1": 50000})
    assert code == 0
    assert table == []
    assert summary["n_is_total"] == "0"


# ── "The tool produced nothing" vs "the genome has nothing" ──────────────────
# An empty IS table has two completely different meanings for anyone judging a
# mobility call, so the audit file must always say which one applies (CLAUDE.md:
# every such decision gets an audit row with a reason).

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
        isescan_row(isBegin="6200", isEnd="5000"),   # reversed -> dropped
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
    results = write_isescan_tsv(tmp_path / "r.tsv",
                                [isescan_row(isBegin="6200", isEnd="5000")])
    code, table, _, audit = run_main(tmp_path, results, {"contig_1": 50000})
    assert code == 0
    assert table == []
    assert audit[0]["reason"] == "invalid_coordinate_range"


def test_ragged_row_is_dropped_with_a_reason(tmp_path):
    # A truncated line means the column-to-value mapping cannot be trusted for
    # that line, so its coordinates might be silently wrong -> drop and record.
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

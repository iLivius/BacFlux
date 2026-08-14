"""Tests for isosdb_copy_number.py — turning read depth into an IS copy count.

Everywhere else the module says "on a fragmented assembly the located IS count is
a FLOOR, not a count". isosdb_copy_number.py is what turns that warning into a
number, by counting reads (immune to assembly collapse) instead of contigs.

So the failure that matters most is an INFLATED estimate: reporting that the
assembly lost copies it never had would manufacture a caveat rather than measure
one. Most of the tests below push in that direction — a partly covered database
entry, depth well under one copy, a family ISOSDB has never seen.

The fixtures are BBMap covstats rows written by hand, so no ISOSDB download and
no read mapping are needed to run any of this.

Run: pytest workflow/scripts/80_mobilome/test_isosdb_copy_number.py -q
"""

import isosdb_copy_number as ic


# ── Building BBMap covstats rows for ISOSDB and for the assembly ─────────────
# Two covstats tables go in: the sample's reads against ISOSDB (numerator) and
# the same reads against the sample's own assembly (denominator). Avg_fold is the
# depth column; Covered_percent says how much of a database entry the reads
# actually touched.

def contig(depth, length):
    return {"Avg_fold": str(depth), "Length": str(length)}


def entry(name, depth, covered=100.0, length=1200):
    return {
        "ID": name, "Avg_fold": str(depth), "Length": str(length),
        "Covered_percent": str(covered),
    }


def run(isosdb_rows, assembly_rows, families=None, located=None, **kwargs):
    return ic.estimate_copies(
        "S1", isosdb_rows, assembly_rows, families or {}, located or {}, **kwargs)


def reasons(audit_rows):
    return {row["reason"] for row in audit_rows}


def by_family(summary_rows):
    return {row["is_family"]: row for row in summary_rows}


# ── The denominator: what one copy looks like in this sample ─────────────────

def test_the_baseline_is_the_depth_of_a_typical_base():
    """Length-weighted median, not a mean over contigs.

    A bacterial assembly is one long chromosome plus a few short contigs, and the
    short ones have the least stable depth — a collapsed repeat sits at several
    times genome depth. An unweighted average lets a 2 kb outlier count as much as
    a 5 Mb chromosome and would drag the baseline up, which deflates every copy
    number derived from it.
    """
    rows = [contig(50, 5_000_000), contig(400, 2_000), contig(3, 1_500)]
    assert ic.genome_baseline_depth(rows) == 50.0


def test_no_assembly_depth_means_no_estimate_rather_than_a_guess():
    """Without a single-copy baseline a raw IS depth is meaningless — it scales
    with how deeply the sample happened to be sequenced."""
    summary, audit = run([entry("ISOSDB1", 200)], [])
    assert summary == []
    assert "no_genome_baseline_depth" in reasons(audit)


# ── Turning a depth ratio into a copy count ──────────────────────────────────
# copy number ~= depth over the IS / depth over the genome as a whole. Every test
# here is really about one of the two ways that ratio can lie: a database entry
# the reads only partly cover (a shared catalytic domain, not this element), and
# a depth so far below the baseline that it is background rather than a copy.

def test_a_five_fold_depth_ratio_reads_as_five_copies():
    summary, _audit = run(
        [entry("ISOSDB1", 250)],
        [contig(50, 5_000_000)],
        families={"ISOSDB1": "IS3"},
        located={"IS3": 2},
    )
    row = by_family(summary)["IS3"]
    assert row["estimated_copies"] == "5.0"
    assert row["located_copies"] == "2"
    # The delta IS the deliverable: three copies the assembler collapsed.
    assert row["collapse_delta"] == "3.0"


def test_family_totals_sum_across_database_entries():
    """ISOSDB is dereplicated at 95%, but one real element can still attract reads
    across several near-identical entries. The family sum is the robust number;
    the per-entry split depends on where an ambiguous read happened to land."""
    summary, _audit = run(
        [entry("ISOSDB1", 100), entry("ISOSDB2", 150)],
        [contig(50, 5_000_000)],
        families={"ISOSDB1": "IS3", "ISOSDB2": "IS3"},
    )
    row = by_family(summary)["IS3"]
    assert row["estimated_copies"] == "5.0"          # 2x + 3x
    assert row["n_db_entries_detected"] == "2"


def test_a_partially_covered_entry_is_excluded_with_a_reason():
    """Depth over a reference the reads barely touched is not that element's
    depth — it is usually a conserved domain shared with another family. Letting
    it through would inflate the estimate, which is the failure direction that
    manufactures a caveat rather than measuring one."""
    summary, audit = run(
        [entry("ISOSDB1", 500, covered=12.0)],
        [contig(50, 5_000_000)],
        families={"ISOSDB1": "IS3"},
    )
    assert by_family(summary) == {} or by_family(summary)["IS3"]["estimated_copies"] == "0.0"
    assert "database_entry_not_fully_covered" in reasons(audit)


def test_depth_far_below_the_baseline_is_not_a_copy():
    summary, audit = run(
        [entry("ISOSDB1", 2)],                       # 0.04x
        [contig(50, 5_000_000)],
        families={"ISOSDB1": "IS3"},
    )
    assert not by_family(summary)
    assert "depth_below_single_copy" in reasons(audit)


def test_a_noisy_single_copy_still_counts():
    """The threshold sits below 1.0 on purpose: sampling noise and mapping loss
    routinely push a genuine single copy to 0.6–0.8x, and discarding those would
    make the estimate read lower than the located count."""
    summary, _audit = run(
        [entry("ISOSDB1", 35)],                      # 0.7x
        [contig(50, 5_000_000)],
        families={"ISOSDB1": "IS3"},
    )
    assert by_family(summary)["IS3"]["estimated_copies"] == "0.7"


# ── Reporting the delta, and saying so when there is no delta ────────────────

def test_a_family_absent_from_isosdb_reports_NA_not_a_negative_delta():
    """The failure the first real run exposed, and the most important test here.

    ISEScan finds IS elements by profile HMM, which recognises a FAMILY. This leg
    needs nucleotide identity to a specific catalogued element. For an organism
    ISOSDB does not cover — the normal case for an environmental isolate — no
    entry attracts reads, the estimate is zero and the delta goes NEGATIVE.

    Read naively that says "the assembly collapsed nothing, we over-called", when
    the truth is "this database has nothing to say about this organism". Measured
    on hybrid sample 006 (Aquipseudomonas): 17 IS located, 1 of 22,713 ISOSDB
    entries covered end to end. Absence of a nucleotide match is not evidence of
    absence of copies, so it must report NA.
    """
    summary, audit = run(
        [], [contig(50, 5_000_000)], families={}, located={"IS110": 3})
    row = by_family(summary)["IS110"]
    assert row["located_copies"] == "3"
    assert row["estimated_copies"] == "NA"
    assert row["collapse_delta"] == "NA"
    assert row["db_informative"] == "FALSE"
    assert "family_absent_from_isosdb" in reasons(audit)


def test_a_sample_ISOSDB_cannot_speak_to_at_all_says_so_once():
    """When no family is covered the whole leg produced nothing, and the located
    count remains an UNQUANTIFIED floor. That has to be stated, or a reader sees
    a table of NAs and no explanation."""
    _summary, audit = run(
        [], [contig(50, 5_000_000)], families={}, located={"IS3": 9, "IS110": 2})
    assert "isosdb_does_not_cover_this_organism" in reasons(audit)


def test_entries_without_a_family_are_reported_not_dropped():
    """About 570 of ISOSDB's 22,713 entries carry no family annotation."""
    summary, _audit = run(
        [entry("ISOSDB99999", 100)], [contig(50, 5_000_000)], families={})
    assert "unassigned" in by_family(summary)


def test_the_summary_audit_explains_what_the_delta_means():
    _summary, audit = run(
        [entry("ISOSDB1", 250)], [contig(50, 5_000_000)],
        families={"ISOSDB1": "IS3"}, located={"IS3": 2})
    line = [a for a in audit if a["reason"] == "copy_number_estimate_complete"]
    assert line
    assert "collapse" in line[0]["detail"]
    # It must not claim to know WHERE the extra copies are.
    assert "cannot say WHERE" in line[0]["detail"]


# ── Degrading when a file is missing or empty ────────────────────────────────
# The whole leg only runs in illumina and hybrid mode (it needs reads) and only
# when an ISOSDB source is configured, so every input here can legitimately be
# absent. Nothing downstream consumes the summary — it is read by a person — so
# an absent input costs a caveat, never a run.

def test_missing_files_are_not_fatal():
    assert ic.read_covstats("/nonexistent.tsv") == []
    assert ic.read_family_map("/nonexistent.txt") == {}
    assert ic.read_located_families("/nonexistent.tsv") == {}


def test_a_header_only_covstats_is_empty_not_an_error(tmp_path):
    path = tmp_path / "cov.tsv"
    path.write_text("#ID\tAvg_fold\tLength\tCovered_percent\n")
    assert ic.read_covstats(str(path)) == []

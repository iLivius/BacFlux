"""Unit tests for traits_summary.py, the strains x traits grid of the trait layer.

traits_summary.py copies each sample's trait decisions into one table, and the same
grid as 2 / 1 / 0 into a file MultiQC draws as a heatmap. It decides nothing, so
what matters is that every call lands in the right row and the right column, and
that profiles which do not line up stop the run instead of shifting calls under
the wrong trait.

The profiles are invented and written into pytest's tmp_path with the header that
trait_profile.py writes (see PROFILE_COLUMNS there); sample names are SAMPLE1 and
so on.

Run: pytest workflow/scripts/90_report/test_traits_summary.py -q
"""

import pytest

import traits_summary as ts


# The header of {sample}_trait_profile.tsv, as trait_profile.py writes it.
PROFILE_HEADER = ["sample", "trait_id", "trait", "application", "decision", "genes",
                  "missing", "supporting", "borderline", "structure_only", "caveat"]


def write_profile(tmp_path, sample, calls, name=None):
    """Write a trait profile: calls = [(trait_id, decision), ...]; the other
    columns get plausible filler."""
    lines = ["\t".join(PROFILE_HEADER)]
    for trait_id, decision in calls:
        row = {column: "NA" for column in PROFILE_HEADER}
        row.update({"sample": sample, "trait_id": trait_id, "decision": decision,
                    "trait": f"Invented trait {trait_id}", "application": "context",
                    "caveat": "Genes present, activity not shown."})
        lines.append("\t".join(row[column] for column in PROFILE_HEADER))
    path = tmp_path / (name or f"{sample}_trait_profile.tsv")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def run(tmp_path, profiles):
    """Run main() and return (summary lines, MultiQC file lines)."""
    summary = tmp_path / "09.report" / "traits_summary.tsv"
    mqc = tmp_path / "09.report" / "traits_summary_mqc.tsv"
    exit_code = ts.main(["--profiles", *profiles, "--out-summary", str(summary),
                         "--out-mqc", str(mqc)])
    assert exit_code == 0
    return (summary.read_text(encoding="utf-8").splitlines(),
            mqc.read_text(encoding="utf-8").splitlines())


# Two invented profiles over the same three traits; CALLS_1 holds each of the
# three decisions once.
CALLS_1 = [("trait_a", "called"), ("trait_b", "partial"), ("trait_c", "not called")]
CALLS_2 = [("trait_a", "not called"), ("trait_b", "called"), ("trait_c", "called")]


def test_grid_has_one_row_per_sample_and_one_column_per_trait(tmp_path):
    """Each sample's calls land in its own row, under the right trait columns, in the
    trait order of the profiles."""
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1),
                write_profile(tmp_path, "SAMPLE2", CALLS_2)]
    summary, _ = run(tmp_path, profiles)
    assert summary == [
        "sample\ttrait_a\ttrait_b\ttrait_c",
        "SAMPLE1\tcalled\tpartial\tnot called",
        "SAMPLE2\tnot called\tcalled\tcalled",
    ]


def test_samples_keep_the_order_they_were_given_in(tmp_path):
    """Rule traits_summary passes the profiles in the order of SAMPLES (the sample
    names BacFlux read from the input file names, sorted by name), and the grid
    keeps the order it is given, so it reads like the rest of the output."""
    profiles = [write_profile(tmp_path, "SAMPLE2", CALLS_2),
                write_profile(tmp_path, "SAMPLE1", CALLS_1)]
    summary, _ = run(tmp_path, profiles)
    assert [line.split("\t")[0] for line in summary[1:]] == ["SAMPLE2", "SAMPLE1"]


def test_multiqc_file_holds_the_same_grid_as_numbers(tmp_path):
    """MultiQC colours a heatmap only from numbers: 2 called, 1 partial, 0 not
    called, under the '#' lines MultiQC reads as its instructions."""
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1),
                write_profile(tmp_path, "SAMPLE2", CALLS_2)]
    _, mqc = run(tmp_path, profiles)
    header_lines = [line for line in mqc if line.startswith("#")]
    table_lines = [line for line in mqc if not line.startswith("#")]
    assert header_lines == ts.MULTIQC_HEADER
    assert mqc[:len(header_lines)] == header_lines        # the '#' lines come first
    assert table_lines == [
        "sample\ttrait_a\ttrait_b\ttrait_c",
        "SAMPLE1\t2\t1\t0",
        "SAMPLE2\t0\t2\t2",
    ]


def test_multiqc_header_is_valid_yaml_for_a_heatmap():
    """MultiQC strips the # and reads the rest as YAML. A stray quote would make it
    drop the trait section, logging an error but still writing the report (rule
    multiqc does not pass --strict), so the header is checked here. The last #
    line keeps a ':' so that MultiQC can never take it for the column names (it
    does so for a last # line with no ':' and as many tab-separated fields as the
    first row). The colour bar is off: three categories, not a scale."""
    yaml = pytest.importorskip("yaml")
    config = yaml.safe_load("\n".join(line[1:] for line in ts.MULTIQC_HEADER))
    assert config["id"] == "bacflux_traits"
    assert config["plot_type"] == "heatmap"
    assert config["pconfig"]["min"] == 0
    assert config["pconfig"]["max"] == 2
    assert config["pconfig"]["xcats_samples"] is False
    assert config["pconfig"]["legend"] is False
    assert "2 = called (dark blue)" in config["description"]
    assert ":" in ts.MULTIQC_HEADER[-1]


def test_sample_name_comes_from_the_sample_column_not_the_file_name(tmp_path):
    """The row name is the profile's sample column, not its file name, so a renamed or
    moved file cannot mislabel a strain."""
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1, name="renamed.tsv")]
    summary, _ = run(tmp_path, profiles)
    assert summary[1].startswith("SAMPLE1\t")


def test_profiles_with_different_traits_stop_the_run(tmp_path):
    """Profiles from two different trait tables would put calls under the wrong
    column; the run stops and names the file."""
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1),
                write_profile(tmp_path, "SAMPLE2", list(reversed(CALLS_2)))]
    with pytest.raises(ValueError, match="its traits differ"):
        run(tmp_path, profiles)


def test_the_same_sample_twice_stops_the_run(tmp_path):
    """Two profiles of one sample would give the strain two rows that may disagree;
    the run stops and names both files."""
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1),
                write_profile(tmp_path, "SAMPLE1", CALLS_2, name="copy.tsv")]
    with pytest.raises(ValueError, match="sample SAMPLE1 is also in"):
        run(tmp_path, profiles)


def test_unknown_decision_stops_the_run(tmp_path):
    """A decision other than called, partial or not called has no number for the
    heatmap and means the file did not come from trait_profile.py: the run stops."""
    profiles = [write_profile(tmp_path, "SAMPLE1", [("trait_a", "maybe")])]
    with pytest.raises(ValueError, match="unknown decision 'maybe'"):
        run(tmp_path, profiles)


def test_file_that_is_not_a_profile_stops_the_run(tmp_path):
    """An empty file, or a table without the sample, trait_id and decision columns,
    stops the run instead of giving an empty or shifted grid."""
    empty = tmp_path / "empty_trait_profile.tsv"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="empty file"):
        run(tmp_path, [str(empty)])

    other = tmp_path / "other.tsv"
    other.write_text("sample\tvalue\nSAMPLE1\t3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing column"):
        run(tmp_path, [str(other)])


def test_profile_naming_two_samples_stops_the_run(tmp_path):
    """A profile holds one sample. A file naming two (two profiles pasted together,
    say) stops the run, since its calls cannot go in one row."""
    path = write_profile(tmp_path, "SAMPLE1", CALLS_1)
    text = open(path, encoding="utf-8").read().replace("SAMPLE1\ttrait_c", "SAMPLE9\ttrait_c")
    open(path, "w", encoding="utf-8").write(text)
    with pytest.raises(ValueError, match="names 2 samples"):
        run(tmp_path, [path])

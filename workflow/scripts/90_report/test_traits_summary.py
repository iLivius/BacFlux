"""Tests for traits_summary.py: the trait grid and its MultiQC plots.

traits_summary.py copies each sample's trait decisions into one table and into
the Trait profile section of the MultiQC report. It decides nothing, so what
matters is that every call lands in the right row and the right column, that the
plots label and colour the three decisions plainly, and that the catalogue it
quotes is the one the calls were made with: profiles that do not line up, or a
catalogue or profile changed since the calls were made, stop the run.

The profiles, catalogue and provenance files are invented and written into
pytest's tmp_path, with the headers trait_profile.py writes; sample names are
SAMPLE1 and so on.

Run: pytest workflow/scripts/90_report/test_traits_summary.py -q
"""

import csv
import hashlib
import json
import os
import subprocess

import pytest

import traits_summary as ts


# The header of {sample}_trait_profile.tsv, as trait_profile.py writes it.
PROFILE_HEADER = ["sample", "trait_id", "trait", "application", "decision", "genes",
                  "missing", "supporting", "borderline", "structure_only", "caveat"]
RULE_HEADER = ["trait_id", "trait", "application", "required", "caveat", "assay", "reference"]

# Two invented profiles over the same three traits; CALLS_1 holds each of the
# three decisions once.
CALLS_1 = [("trait_a", "called"), ("trait_b", "partial"), ("trait_c", "not called")]
CALLS_2 = [("trait_a", "not called"), ("trait_b", "called"), ("trait_c", "called")]


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


def write_rules(tmp_path, trait_ids, applications=None):
    """Write a small catalogue for the given traits, with a "#" line on top as
    trait_rules_v1.tsv has."""
    applications = applications or {}
    path = tmp_path / "rules.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# Invented trait rules for the tests.\n")
        writer = csv.DictWriter(handle, fieldnames=RULE_HEADER, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for trait_id in trait_ids:
            writer.writerow({
                "trait_id": trait_id, "trait": f"Invented trait {trait_id}",
                "application": applications.get(trait_id, "context"), "required": "geneA",
                "caveat": "Genes present, activity not shown.",
                "assay": "Culture assay (Doe 2026, doi:10.1234/assay).",
                "reference": "Doe et al. 2026, doi:10.1234/trait.",
            })
    return str(path)


def sha256(path):
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def write_provenance(tmp_path, sample, profile, rules):
    """Write the provenance record trait_profile.py would write with a profile:
    the fingerprints of the catalogue and of the profile."""
    record = {
        "schema_version": "1.0", "sample": sample,
        "inputs": {"rules": {"file": "rules.tsv", "sha256": sha256(rules)}},
        "outputs": {"profile": {"file": os.path.basename(profile), "sha256": sha256(profile)}},
    }
    path = tmp_path / f"{sample}_trait_provenance.json"
    path.write_text(json.dumps(record), encoding="utf-8")
    return str(path)


def run(tmp_path, profiles, applications=None, rules=None, provenance=None):
    """Run main() with a matching catalogue and provenance files, and return the
    summary lines and the plot sections."""
    if rules is None:
        trait_ids = [trait_id for trait_id, _d in ts.read_profile(profiles[0])[1]]
        rules = write_rules(tmp_path, trait_ids, applications)
    if provenance is None:
        provenance = [write_provenance(tmp_path, ts.read_profile(p)[0], p, rules)
                      for p in profiles]
    summary = tmp_path / "09.report" / "traits_summary.tsv"
    mqc = tmp_path / "09.report" / "traits_summary_mqc"
    exit_code = ts.main(["--profiles", *profiles, "--provenance", *provenance,
                         "--rules", rules, "--out-summary", str(summary),
                         "--out-mqc", str(mqc)])
    assert exit_code == 0
    sections = [json.loads(path.read_text(encoding="utf-8"))
                for path in sorted(mqc.glob("*_mqc.json"))]
    return summary.read_text(encoding="utf-8").splitlines(), sections


def heatmaps_of(sections):
    return [section for section in sections if section["plot_type"] == "heatmap"]


# ── The table ────────────────────────────────────────────────────────────────

def test_grid_has_one_row_per_sample_and_one_column_per_trait(tmp_path):
    """Each sample's calls land in its own row, under the right trait columns, in
    catalogue order."""
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1),
                write_profile(tmp_path, "SAMPLE2", CALLS_2)]
    summary, _ = run(tmp_path, profiles)
    assert summary == [
        "sample\ttrait_a\ttrait_b\ttrait_c",
        "SAMPLE1\tcalled\tpartial\tnot called",
        "SAMPLE2\tnot called\tcalled\tcalled",
    ]


def test_samples_keep_the_order_they_were_given_in(tmp_path):
    """Rule traits_summary passes the profiles in the order of SAMPLES, and the
    grid keeps the order it is given, so it reads like the rest of the output."""
    profiles = [write_profile(tmp_path, "SAMPLE2", CALLS_2),
                write_profile(tmp_path, "SAMPLE1", CALLS_1)]
    summary, _ = run(tmp_path, profiles)
    assert [line.split("\t")[0] for line in summary[1:]] == ["SAMPLE2", "SAMPLE1"]


def test_sample_name_comes_from_the_sample_column_not_the_file_name(tmp_path):
    """A renamed or moved file cannot mislabel a strain."""
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1, name="renamed.tsv")]
    summary, _ = run(tmp_path, profiles)
    assert summary[1].startswith("SAMPLE1\t")


# ── The plots ────────────────────────────────────────────────────────────────

def test_heatmap_holds_the_same_grid_with_trait_ids_on_the_axis(tmp_path):
    """Every decision keeps its row and column. The axis shows the trait IDs,
    passed as they are: one short label per trait, no name repeated beside it."""
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1),
                write_profile(tmp_path, "SAMPLE2", CALLS_2)]
    _, sections = run(tmp_path, profiles)
    heatmap = heatmaps_of(sections)[0]
    assert heatmap["data"] == [[2, 1, 0], [0, 2, 2]]
    assert heatmap["ycats"] == ["SAMPLE1", "SAMPLE2"]
    assert heatmap["xcats"] == ["trait_a", "trait_b", "trait_c"]
    assert heatmap["pconfig"]["angled_xticks"] is True
    # The names are listed in the description instead, one per ID.
    for trait_id, _decision in CALLS_1:
        assert f"<code>{trait_id}</code>: Invented trait {trait_id}" in heatmap["description"]


def test_colours_and_legend_keep_three_fixed_categories(tmp_path):
    """Three flat colour bands from 0 to 2, no colour bar, and the legend once, in
    the description of the Trait profile section, with each decision in words."""
    _, sections = run(tmp_path, [write_profile(tmp_path, "SAMPLE1", CALLS_1)])
    heatmap = heatmaps_of(sections)[0]
    config = heatmap["pconfig"]
    assert (config["min"], config["max"]) == (0, 2)
    assert config["xcats_samples"] is False
    assert config["legend"] is False
    assert config["colstops"] == ts.COLOUR_STOPS
    parent = heatmap["parent_description"]
    assert "every gene of a required set counts" in parent
    assert "a lead, not half a trait" in parent
    assert "does not show that the trait is absent" in parent
    assert "not a score" in parent
    # The grey of "not called" is visible against a white page.
    assert ts.DECISION_COLOURS["not called"] == "#d9d9d9"
    # The legend is not repeated in every panel.
    assert "a lead, not half a trait" not in heatmap["description"]


def test_application_groups_do_not_duplicate_shared_traits(tmp_path):
    """Grouping is by the exact application label, so a trait with two labels
    appears once; the panel is named after its label, without a number."""
    applications = {"trait_a": "biofertiliser", "trait_b": "biofertiliser, biocontrol"}
    _, sections = run(tmp_path, [write_profile(tmp_path, "SAMPLE1", CALLS_1)], applications)
    heatmaps = heatmaps_of(sections)
    overview = next(section for section in sections if section["plot_type"] == "bargraph")
    assert len(heatmaps) == 3
    assert sum(len(section["xcats"]) for section in heatmaps) == 3
    assert sorted(section["section_name"] for section in heatmaps) == [
        "Biofertiliser", "Biofertiliser, biocontrol", "Context"]
    assert overview["data"][0]["SAMPLE1"] == {"called": 1, "partial": 1, "not called": 1}
    assert overview["section_name"] == "Call counts"
    assert all(section["parent_id"] == "bacflux_traits" for section in sections)


def test_overview_data_is_a_list_of_one_dataset(tmp_path):
    """A plain {sample: counts} dict makes MultiQC shorten every strain's name to
    the report folder's under -d; one dataset in a list keeps the names."""
    _, sections = run(tmp_path, [write_profile(tmp_path, "SAMPLE1", CALLS_1)])
    overview = next(section for section in sections if section["plot_type"] == "bargraph")
    assert isinstance(overview["data"], list) and len(overview["data"]) == 1


def test_names_are_escaped_in_descriptions_but_labels_are_not_changed():
    """A trait name is HTML in a description, so a "<" in it is escaped there. The
    axis labels are the IDs and sample names as they are: escaping them would put
    codes such as &#x27; on the axis."""
    rules = {"a'b": {"trait": "Trait <script>bad</script> it's", "application": "context"}}
    heatmap = ts.heatmap_sections(["a'b"], [("SAMPLE1", ["called"])], rules)[0]
    assert heatmap["xcats"] == ["a'b"]
    assert "<script>" not in heatmap["description"]
    assert "&lt;script&gt;" in heatmap["description"]
    assert "it's" in heatmap["description"]


def test_rerun_removes_only_stale_owned_plot_files(tmp_path):
    folder = tmp_path / "plots"
    folder.mkdir()
    stale = folder / "bacflux_traits_application_99_mqc.json"
    stale.write_text("{}", encoding="utf-8")
    note = folder / "review-notes.txt"
    note.write_text("keep", encoding="utf-8")
    ts.write_multiqc_folder(folder, [{"id": "bacflux_traits_overview", "data": {"S": 1}}])
    assert not stale.exists()
    assert note.read_text() == "keep"


def test_references_section_links_sources_and_assays(tmp_path):
    _, sections = run(tmp_path, [write_profile(tmp_path, "SAMPLE1", CALLS_1)])
    references = next(section for section in sections if section["plot_type"] == "html")
    assert references["id"] == "bacflux_traits_references"
    assert references["section_name"] == "References and assays"
    assert "Suggested laboratory assay" in references["data"]
    assert "Caveat" in references["data"]
    assert 'href="https://doi.org/10.1234/trait"' in references["data"]
    assert 'href="https://doi.org/10.1234/assay"' in references["data"]
    assert references["data"].count("<details ") == 3
    assert all(section["doi"] == ts.METHOD_DOIS for section in sections)


def test_reference_text_is_escaped_and_doi_parentheses_are_preserved():
    linked = ts.linked_references(
        '<script>bad</script> (doi:10.1016/S0167-7012(99)00096-2).'
    )
    assert "<script>" not in linked and "&lt;script&gt;" in linked
    assert 'href="https://doi.org/10.1016/S0167-7012(99)00096-2"' in linked
    assert linked.endswith("</a>).")


# ── The catalogue and the provenance records ─────────────────────────────────

def test_a_catalogue_edited_after_the_calls_stops_the_run(tmp_path):
    """The provenance records the catalogue's fingerprint when the calls were
    made; an edit since then, even to a reference, stops the run before anything
    is written."""
    profile = write_profile(tmp_path, "SAMPLE1", CALLS_1)
    rules = write_rules(tmp_path, ["trait_a", "trait_b", "trait_c"])
    provenance = [write_provenance(tmp_path, "SAMPLE1", profile, rules)]
    text = open(rules, encoding="utf-8").read().replace("doi:10.1234/trait", "doi:10.1234/other")
    open(rules, "w", encoding="utf-8").write(text)
    with pytest.raises(ValueError, match="made with another version"):
        run(tmp_path, [profile], rules=rules, provenance=provenance)
    assert not (tmp_path / "09.report").exists()


def test_an_edited_profile_stops_the_run(tmp_path):
    profile = write_profile(tmp_path, "SAMPLE1", CALLS_1)
    rules = write_rules(tmp_path, ["trait_a", "trait_b", "trait_c"])
    provenance = [write_provenance(tmp_path, "SAMPLE1", profile, rules)]
    text = open(profile, encoding="utf-8").read().replace("\tpartial\t", "\tcalled\t")
    open(profile, "w", encoding="utf-8").write(text)
    with pytest.raises(ValueError, match="not the one rule trait_profile wrote"):
        run(tmp_path, [profile], rules=rules, provenance=provenance)


def test_every_profile_needs_its_provenance_record(tmp_path):
    one = write_profile(tmp_path, "SAMPLE1", CALLS_1)
    two = write_profile(tmp_path, "SAMPLE2", CALLS_2)
    rules = write_rules(tmp_path, ["trait_a", "trait_b", "trait_c"])
    provenance = [write_provenance(tmp_path, "SAMPLE1", one, rules)]
    with pytest.raises(ValueError, match="no provenance file for sample"):
        run(tmp_path, [one, two], rules=rules, provenance=provenance)


def test_catalogue_traits_must_be_the_profiles_traits_in_order(tmp_path):
    profile = write_profile(tmp_path, "SAMPLE1", CALLS_1)
    rules = write_rules(tmp_path, ["trait_a", "trait_c", "trait_b"])
    with pytest.raises(ValueError, match="another catalogue"):
        run(tmp_path, [profile], rules=rules)


# ── Profiles that cannot be lined up stop the run ────────────────────────────

def test_duplicate_trait_id_stops_the_run(tmp_path):
    path = write_profile(tmp_path, "SAMPLE1", [("trait_a", "called"), ("trait_a", "partial")])
    with pytest.raises(ValueError, match="duplicate trait_id"):
        ts.build_grid([path])


def test_profiles_with_different_traits_stop_the_run(tmp_path):
    """Profiles from two different trait tables would put calls under the wrong
    column; the run stops and names the file."""
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1),
                write_profile(tmp_path, "SAMPLE2", list(reversed(CALLS_2)))]
    with pytest.raises(ValueError, match="its traits differ"):
        ts.build_grid(profiles)


def test_the_same_sample_twice_stops_the_run(tmp_path):
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1),
                write_profile(tmp_path, "SAMPLE1", CALLS_2, name="copy.tsv")]
    with pytest.raises(ValueError, match="sample SAMPLE1 is also in"):
        ts.build_grid(profiles)


def test_unknown_decision_stops_the_run(tmp_path):
    profiles = [write_profile(tmp_path, "SAMPLE1", [("trait_a", "maybe")])]
    with pytest.raises(ValueError, match="unknown decision 'maybe'"):
        ts.build_grid(profiles)


def test_file_that_is_not_a_profile_stops_the_run(tmp_path):
    empty = tmp_path / "empty_trait_profile.tsv"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="empty file"):
        ts.build_grid([str(empty)])
    other = tmp_path / "other.tsv"
    other.write_text("sample\tvalue\nSAMPLE1\t3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing column"):
        ts.build_grid([str(other)])


def test_profile_naming_two_samples_stops_the_run(tmp_path):
    path = write_profile(tmp_path, "SAMPLE1", CALLS_1)
    text = open(path, encoding="utf-8").read().replace("SAMPLE1\ttrait_c", "SAMPLE9\ttrait_c")
    open(path, "w", encoding="utf-8").write(text)
    with pytest.raises(ValueError, match="names 2 samples"):
        ts.build_grid([path])


# ── Optional: render with the pinned MultiQC ─────────────────────────────────

def test_multiqc_keeps_sample_names_and_shows_the_legend(tmp_path):
    """Optional check with the real renderer, against the bug that shortened every
    strain's name to "09" under -d.

    Set BACFLUX_TEST_MULTIQC to the pinned MultiQC executable to run it; the
    ordinary tests need only Python and pytest.
    """
    executable = os.environ.get("BACFLUX_TEST_MULTIQC")
    if not executable:
        pytest.skip("set BACFLUX_TEST_MULTIQC to run the pinned renderer")
    profiles = [write_profile(tmp_path, "SAMPLE1", CALLS_1),
                write_profile(tmp_path, "SAMPLE2", CALLS_2)]
    run(tmp_path, profiles)
    result = subprocess.run(
        [executable, "09.report/traits_summary_mqc", "-d", "--strict",
         "--no-version-check", "--outdir", "rendered"],
        cwd=tmp_path, capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((tmp_path / "rendered/multiqc_data/multiqc_data.json").read_text())
    data = report["report_saved_raw_data"]["multiqc_bacflux_traits_overview_plot"]
    assert set(data) == {"SAMPLE1", "SAMPLE2"}
    assert data["SAMPLE1"] == {"called": 1, "partial": 1, "not called": 1}
    rendered = (tmp_path / "rendered/multiqc_report.html").read_text()
    assert "References and assays" in rendered
    assert 'href="https://doi.org/10.1234/trait"' in rendered
    assert "a lead, not half a trait" in rendered
    citations = (tmp_path / "rendered/multiqc_data/multiqc_citations.txt").read_text()
    assert all(doi in citations for doi in ts.METHOD_DOIS)

#!/usr/bin/env python3
"""Set every sample's trait profile side by side, and draw it for the MultiQC report.

Runs as rule traits_summary (shared/90_report.smk), once per run, after rule
trait_profile (shared/40_annotation.smk) has written every sample's
{sample}_trait_profile.tsv and {sample}_trait_provenance.json. Nothing is decided
again here: each cell is copied from that sample's profile, so the reasons stay in
the per-sample files.

What it reads
-------------
  --profiles     every sample's {sample}_trait_profile.tsv (one row per trait, in
                 the order of trait_rules_v1.tsv). The sample name is read from
                 the sample column, not from the file name.
  --provenance   every sample's {sample}_trait_provenance.json, written with the
                 profile: it holds the SHA-256 fingerprints of the catalogue and
                 of the profile as they were when the calls were made.
  --rules        trait_rules_v1.tsv, the trait catalogue. The report shows each
                 trait's name, application, references, suggested assay and
                 caveat from it, beside calls made earlier. Its fingerprint must
                 therefore match every provenance record, and every profile must
                 still be the one recorded: a catalogue edited after the calls,
                 or an edited profile, stops the rule instead of being shown.

What it writes
--------------
  --out-summary  09.report/traits_summary.tsv: a sample column, then one column
                 per trait_id in catalogue order; each cell is called, partial or
                 not called. One row per sample, in the order the profiles were
                 given (rule traits_summary passes them in the order of SAMPLES,
                 sorted by name).
  --out-mqc      09.report/traits_summary_mqc/, the plot files of the Trait
                 profile section of the MultiQC report:
                   bacflux_traits_overview_mqc.json        a bar per sample: how
                                                           many traits are called,
                                                           partial, not called;
                   bacflux_traits_application_NN_mqc.json  one heatmap per
                                                           application label of
                                                           the catalogue, samples
                                                           in rows, traits in
                                                           columns;
                   bacflux_traits_references_mqc.json      each trait's references,
                                                           suggested assay and
                                                           caveat.
                 MultiQC reads them as custom content (plot files another program
                 wrote for MultiQC, found by the _mqc in their names). Rule multiqc
                 places the overview first and the references last; the heatmaps
                 follow in alphabetical order of their application label, which is
                 how MultiQC sorts sections of one module.

Standard library only (Python 3.12 or later). Rule traits_summary runs it in
workflow/envs/python.yaml. Tests: test_traits_summary.py.
"""

import argparse
import hashlib
import html
import json
import os
import re
import sys
from pathlib import Path
from urllib.parse import quote


# ── The three decisions, their numbers and colours ───────────────────────────

# The three decisions trait_profile.py writes, and the number each becomes in a
# heatmap. MultiQC colours a heatmap only from numbers; the numbers stand for
# categories, not for a score, and the legend says so.
DECISION_SCORES = {"called": 2, "partial": 1, "not called": 0}

# Dark blue for called, light blue for partial, a mid grey for not called (a
# lighter grey nearly vanishes against the white page).
DECISION_COLOURS = {"called": "#08519c", "partial": "#9ecae1", "not called": "#d9d9d9"}

# The heatmap's colour scale, fixed from 0 to 2 in three flat bands so each
# decision always gets its own colour. MultiQC maps 0..2 onto 0..1: 0 falls in
# the first band, 1 in the middle one, 2 in the last. Without the fixed bands, a
# batch with no called trait would stretch the scale to 0..1 and draw partial
# calls in the colour of a call.
COLOUR_STOPS = [
    [0, DECISION_COLOURS["not called"]], [0.25, DECISION_COLOURS["not called"]],
    [0.25, DECISION_COLOURS["partial"]], [0.75, DECISION_COLOURS["partial"]],
    [0.75, DECISION_COLOURS["called"]], [1, DECISION_COLOURS["called"]],
]

# The legend, shown once at the top of the Trait profile section. The wording of
# each decision follows docs/analysis/traits.md.
LEGEND = (
    '<p><span style="background:#08519c;color:white;padding:1px 6px">called</span> '
    "every gene of a required set counts, and the genes that must lie together do. "
    '<span style="background:#9ecae1;padding:1px 6px">partial</span> '
    "at least half of a required set of two or more genes counts, or the set is "
    "complete but some of its genes are broken (pseudogenes) or not together where "
    "the rule asks: a lead, not half a trait. "
    '<span style="background:#d9d9d9;padding:1px 6px">not called</span> '
    "no set reached partial; this does not show that the trait is absent. "
    "The three colours are categories, not a score.</p>"
)
# The introduction of the Trait profile section. MultiQC ends a section's
# introduction with a full stop when it does not end with one, so the last
# sentence stands outside a <p> element and ends with its own full stop: after a
# closing </p>, MultiQC's full stop would show as a lone "." under the text.
PARENT_DESCRIPTION = (
    "<p>Traits predicted from each strain's genes (the BacFlux trait profile). A "
    "call says that the genes behind a trait are present; it does not show that the "
    "strain expresses the trait, and it says nothing about whether the strain is safe "
    "to use: see its taxonomy and its AMR workbook (10.workbooks).</p>"
    + LEGEND +
    "<p>The heatmaps group the traits by the catalogue's application label, with "
    "each trait's ID on the axis and its name listed above the heatmap. Traits "
    "labelled <em>context</em> describe the strain's metabolism and do not, on their "
    "own, support a product claim.</p>"
    "The genes behind each call and the reason for each decision are in "
    "04.annotation/traits/{sample}/."
)

# The prefix of every plot file this script writes; see write_multiqc_folder().
GENERATED_PREFIX = "bacflux_traits_"


# ── The methods of the trait layer, for MultiQC's citation button ────────────

# Registered with every trait section, so MultiQC's citation export lists them.
# The papers behind each trait's rule come from the catalogue instead.
METHOD_REFERENCES = (
    ("EcoFoldDB method", "Ghaly et al. (2025)", "10.1111/1462-2920.70178"),
    ("EcoFoldDB v2.0 database", "versioned database archive", "10.5281/zenodo.16908168"),
    ("Baktfold", "Bouras et al. (2026), preprint", "10.64898/2026.03.31.715528"),
    ("ProstT5", "Heinzinger et al. (2024)", "10.1093/nargab/lqae150"),
    ("Foldseek", "van Kempen et al. (2024)", "10.1038/s41587-023-01773-0"),
)
METHOD_DOIS = [doi for _method, _citation, doi in METHOD_REFERENCES]

# A DOI in the catalogue's text, written as "doi:10.xxxx/..." or as a doi.org
# link. Only these are turned into links; the rest of the text is shown as it is.
DOI_PATTERN = re.compile(r"(?:doi:\s*|https?://(?:dx\.)?doi\.org/)(10\.\d{4,9}/[^\s;<>]+)", re.I)


# ── Reading the inputs ───────────────────────────────────────────────────────

def read_profile(path):
    """Read one {sample}_trait_profile.tsv.

    Takes in: a profile written by trait_profile.py.
    Returns:  (sample, [(trait_id, decision), ...]) in the file's order.
    Stops the run when the file has no rows, lacks a needed column, has a row
    with more or fewer fields than the header, names more than one sample,
    repeats a trait, or holds a decision other than the three known ones: each
    means the file did not come from trait_profile.py as expected.
    """
    with open(path, encoding="utf-8") as handle:
        lines = [line.rstrip("\r\n") for line in handle if line.strip()]
    if not lines:
        raise ValueError(f"{path}: empty file, not a trait profile.")
    header = lines[0].split("\t")
    missing = [name for name in ("sample", "trait_id", "decision") if name not in header]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)}; not a trait profile.")

    samples = set()
    calls = []
    for number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"{path}: line {number} has {len(fields)} fields, the "
                             f"header has {len(header)}.")
        row = dict(zip(header, fields))
        if row["decision"] not in DECISION_SCORES:
            raise ValueError(f"{path}: line {number}: unknown decision '{row['decision']}'.")
        if any(trait_id == row["trait_id"] for trait_id, _decision in calls):
            raise ValueError(f"{path}: duplicate trait_id '{row['trait_id']}'.")
        samples.add(row["sample"])
        calls.append((row["trait_id"], row["decision"]))
    if not calls:
        raise ValueError(f"{path}: no trait rows.")
    if len(samples) != 1:
        raise ValueError(f"{path}: names {len(samples)} samples "
                         f"({', '.join(sorted(samples))}); a profile holds one.")
    return samples.pop(), calls


def build_grid(paths):
    """Read every profile and line them up under one list of traits.

    Takes in: the profile paths, in the order the rule passes them.
    Returns:  (trait_ids, [(sample, [decision per trait]), ...], {sample: path}).
    Stops the run when two profiles name the same sample, or list different
    traits or the same traits in a different order, because a grid of mismatched
    columns would put calls under the wrong trait.
    """
    trait_ids = None
    first_path = None
    grid = []
    path_of_sample = {}
    for path in paths:
        sample, calls = read_profile(path)
        if sample in path_of_sample:
            raise ValueError(f"{path}: sample {sample} is also in {path_of_sample[sample]}.")
        path_of_sample[sample] = path
        these_traits = [trait_id for trait_id, _decision in calls]
        if trait_ids is None:
            trait_ids, first_path = these_traits, path
        elif these_traits != trait_ids:
            raise ValueError(f"{path}: its traits differ from those of {first_path} "
                             "(a different trait table?); the grid would mix up columns.")
        grid.append((sample, [decision for _trait_id, decision in calls]))
    if trait_ids is None:
        raise ValueError("No trait profiles were supplied.")
    return trait_ids, grid, path_of_sample


def read_rules(path, trait_ids):
    """Read the trait catalogue into {trait_id: row}.

    Takes in: trait_rules_v1.tsv (its "#" lines at the top are skipped) and the
              trait IDs of the profiles.
    Stops the run when a needed column is missing, a row has the wrong number of
    fields, an ID is blank or repeated, or the catalogue's traits are not the
    profiles' traits in the same order.
    """
    with open(path, encoding="utf-8") as handle:
        lines = [line.rstrip("\r\n") for line in handle
                 if line.strip() and not line.startswith("#")]
    if not lines:
        raise ValueError(f"{path}: no trait rules.")
    header = lines[0].split("\t")
    needed = ("trait_id", "trait", "application", "required", "caveat", "assay", "reference")
    missing = [column for column in needed if column not in header]
    if missing:
        raise ValueError(f"{path}: missing rule column(s): {', '.join(missing)}.")
    rules = {}
    order = []
    for number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"{path}: rule line {number} has {len(fields)} fields, the "
                             f"header has {len(header)}.")
        row = dict(zip(header, fields))
        trait_id = row["trait_id"].strip()
        if not trait_id or trait_id in rules:
            raise ValueError(f"{path}: blank or duplicate trait_id {trait_id!r}.")
        rules[trait_id] = row
        order.append(trait_id)
    if order != list(trait_ids):
        raise ValueError(f"{path}: its traits are not those of the profiles, in the same "
                         "order; the profiles were made with another catalogue.")
    return rules


def sha256_of(path):
    """The SHA-256 fingerprint of a file's contents (see trait_profile.py's
    file_fingerprint, which wrote the ones this is compared with)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for piece in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(piece)
    return digest.hexdigest()


def check_provenance(provenance_paths, path_of_sample, rules_path):
    """Check that the catalogue and every profile are the ones the calls produced.

    Takes in: the provenance files, the profile path of every sample
              (build_grid) and the catalogue's path.
    Each provenance file names its sample; it must belong to a profile given
    here, and every profile needs one. For each, the catalogue's fingerprint must
    equal the one recorded as inputs.rules, and the profile's the one recorded as
    outputs.profile. Any difference stops the run, naming the sample, with rule
    trait_profile to be run again for it.
    """
    rules_sha256 = sha256_of(rules_path)
    checked = set()
    for path in provenance_paths:
        with open(path, encoding="utf-8") as handle:
            record = json.load(handle)
        sample = record.get("sample")
        if sample not in path_of_sample:
            raise ValueError(f"{path}: written for sample {sample!r}, which has no profile here.")
        recorded_rules = ((record.get("inputs") or {}).get("rules") or {}).get("sha256")
        if recorded_rules != rules_sha256:
            raise ValueError(f"{path}: the trait calls of {sample} were made with another "
                             f"version of {os.path.basename(rules_path)}; run rule "
                             "trait_profile again for this sample.")
        recorded_profile = ((record.get("outputs") or {}).get("profile") or {}).get("sha256")
        if recorded_profile != sha256_of(path_of_sample[sample]):
            raise ValueError(f"{path}: the profile of {sample} is not the one rule "
                             "trait_profile wrote; run it again for this sample.")
        checked.add(sample)
    unchecked = sorted(set(path_of_sample) - checked)
    if unchecked:
        raise ValueError(f"no provenance file for sample(s) {', '.join(unchecked)}.")


# ── The MultiQC sections ─────────────────────────────────────────────────────

def section(section_id, name, description):
    """The settings every trait section shares.

    All sections sit under one parent, the Trait profile entry of the report's
    navigation, whose description (with the legend) MultiQC shows once above
    them. Each registers the trait layer's method DOIs for MultiQC's citation
    button.
    """
    return {
        "id": section_id, "section_name": name, "description": description,
        "parent_id": "bacflux_traits", "parent_anchor": "bacflux_traits",
        "parent_name": "Trait profile", "parent_description": PARENT_DESCRIPTION,
        "doi": METHOD_DOIS,
    }


def linked_references(text):
    """Show catalogue text as HTML, with each DOI in it turned into a link.

    Everything else is escaped (html.escape), so a "<" in a citation is shown
    as itself. A DOI often ends a sentence in brackets, as in
    "(doi:10.1016/S0167-7012(99)00096-2).": the closing bracket of the sentence
    stays outside the link, while brackets that belong to the DOI stay inside.
    """
    pieces = []
    position = 0
    for match in DOI_PATTERN.finditer(text):
        doi = match.group(1).rstrip(".,;")
        while doi.endswith(")") and doi.count(")") > doi.count("("):
            doi = doi[:-1]
        end = match.start(1) + len(doi)
        pieces.append(html.escape(text[position:match.start()]))
        url = "https://doi.org/" + quote(doi, safe="/().:_-")
        pieces.append(f'<a href="{url}" target="_blank" rel="noopener noreferrer">'
                      + html.escape(text[match.start():end]) + "</a>")
        position = end
    pieces.append(html.escape(text[position:]))
    return "".join(pieces)


def overview_section(grid):
    """The bar chart: for each sample, how many traits are called, partial and
    not called. It helps to find strains in a batch; it is not a ranking."""
    overview = section(
        "bacflux_traits_overview", "Call counts",
        "<p>For each strain, how many of the catalogue's traits are called, partial "
        "and not called, context traits included. The counts help to find strains "
        "in a batch; they do not rank strains by benefit or safety.</p>")
    overview.update({
        "plot_type": "bargraph",
        "pconfig": {
            "id": "bacflux_traits_overview_plot", "title": "BacFlux: trait calls per strain",
            "ylab": "Traits", "cpswitch": False, "stacking": "normal",
            "sort_samples": False,
        },
        "categories": {
            decision: {"name": decision, "color": DECISION_COLOURS[decision]}
            for decision in DECISION_SCORES
        },
        # The data are given as a list holding one dataset, not as a plain
        # {sample: counts} dict. MultiQC cleans the keys of a plain dict as if
        # they were file names, and with the -d option rule multiqc uses, that
        # cleaning shortened every strain's name to "09" (the report folder).
        "data": [{
            sample: {decision: calls.count(decision) for decision in DECISION_SCORES}
            for sample, calls in grid
        }],
    })
    return overview


def heatmap_sections(trait_ids, grid, rules):
    """One heatmap per application label of the catalogue.

    Traits are grouped by their exact application label, so a trait labelled
    "biofertiliser, biocontrol" appears once, in that panel, not in two. Within
    a panel the traits keep catalogue order. The axis shows each trait's ID,
    which is short and unique; the panel's description lists every ID with the
    trait's name. Labels are passed to the plot as they are, not escaped: the
    plotting library would show an escape code (such as &#x27; for an
    apostrophe) literally, and IDs and sample names carry no HTML anyway.
    """
    groups = {}
    for index, trait_id in enumerate(trait_ids):
        groups.setdefault(rules[trait_id]["application"], []).append(index)

    sections = []
    for number, (application, indices) in enumerate(groups.items(), start=1):
        section_id = f"bacflux_traits_application_{number:02d}"
        names = "".join(
            f"<li><code>{html.escape(trait_ids[i])}</code>: "
            f"{html.escape(rules[trait_ids[i]]['trait'], quote=False)}</li>"
            for i in indices)
        heatmap = section(
            section_id, application[:1].upper() + application[1:],
            f"<p>The {len(indices)} traits with application label "
            f"<em>{html.escape(application, quote=False)}</em>:</p><ul>{names}</ul>")
        heatmap.update({
            "plot_type": "heatmap",
            "xcats": [trait_ids[i] for i in indices],
            "ycats": [sample for sample, _calls in grid],
            "data": [[DECISION_SCORES[calls[i]] for i in indices] for _sample, calls in grid],
            "pconfig": {
                "id": section_id + "_plot",
                "title": "BacFlux: " + application,
                "xlab": "Trait", "ylab": "Strain",
                "zlab": "Call (2 called, 1 partial, 0 not called)",
                "min": 0, "max": 2, "tt_decimals": 0,
                "square": False, "xcats_samples": False, "ycats_samples": True,
                # No colour bar: it would suggest a continuous scale for what are
                # three categories; the legend above names each colour.
                "legend": False, "colstops": COLOUR_STOPS,
                "display_values": False, "angled_xticks": True,
                "cluster_rows": False, "cluster_cols": False,
                # Room for the samples, plus the angled trait IDs below them.
                "height": max(420, len(grid) * 22 + 280),
            },
        })
        sections.append(heatmap)
    return sections


def references_section(trait_ids, rules):
    """The References and assays section: one expandable entry per trait.

    Each entry gives the trait's application, its required genes, the papers
    behind its rule, a suggested laboratory assay and its caveat, all from the
    catalogue, with DOIs as links. The trait layer's methods are listed first.
    """
    methods = "".join(
        f"<li><strong>{html.escape(method)}</strong>: {html.escape(citation)}, "
        + linked_references("doi:" + doi) + ".</li>"
        for method, citation, doi in METHOD_REFERENCES)
    content = [
        "<h4>Methods of the trait profile</h4><ul>", methods, "</ul>",
        "<p>Open a trait for the papers behind its rule, a suggested laboratory "
        "assay and its caveat. The papers support the genomic prediction; the "
        "assay is how to test the trait in the strain.</p>",
    ]
    for trait_id in trait_ids:
        rule = rules[trait_id]
        title = html.escape(f"{rule['trait']} ({trait_id})", quote=False)
        content.extend([
            f'<details class="mb-2"><summary><strong>{title}</strong></summary>',
            '<dl class="mt-2">',
            "<dt>Application</dt><dd>" + html.escape(rule["application"]) + "</dd>",
            "<dt>Required genes</dt><dd>" + html.escape(rule["required"]) + "</dd>",
            "<dt>References</dt><dd>" + linked_references(rule["reference"] or "none given") + "</dd>",
            "<dt>Suggested laboratory assay</dt><dd>"
            + linked_references(rule["assay"] or "none given") + "</dd>",
            "<dt>Caveat</dt><dd>" + html.escape(rule["caveat"]) + "</dd>",
            "</dl></details>",
        ])
    references = section("bacflux_traits_references", "References and assays",
                         "<p>From the trait catalogue (trait_rules_v1.tsv).</p>")
    references.update({"plot_type": "html", "data": "".join(content)})
    return references


def multiqc_sections(trait_ids, grid, rules):
    """All sections of the Trait profile part of the report, in writing order."""
    return ([overview_section(grid)]
            + heatmap_sections(trait_ids, grid, rules)
            + [references_section(trait_ids, rules)])


# ── Writing the outputs ──────────────────────────────────────────────────────

def write_summary(path, trait_ids, grid):
    """Write traits_summary.tsv: one row per sample, one column per trait."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(["sample"] + list(trait_ids)) + "\n")
        for sample, calls in grid:
            handle.write("\t".join([sample] + calls) + "\n")


def write_multiqc_folder(path, sections):
    """Write one plot file per section into the folder.

    Under Snakemake the folder is a directory() output, emptied before the rule
    runs. Run by hand into a folder used before, a plot file this script wrote
    for a section that no longer exists (a catalogue with fewer application
    labels) is removed, so it cannot linger in the report; any other file in the
    folder is left alone.
    """
    folder = Path(path)
    if folder.exists() and not folder.is_dir():
        raise ValueError(f"{path}: --out-mqc must be a folder.")
    folder.mkdir(parents=True, exist_ok=True)
    written = set()
    for plot in sections:
        filename = plot["id"] + "_mqc.json"
        (folder / filename).write_text(json.dumps(plot, ensure_ascii=False, indent=2) + "\n",
                                       encoding="utf-8")
        written.add(filename)
    for previous in folder.glob(GENERATED_PREFIX + "*_mqc.json"):
        if previous.name not in written:
            previous.unlink()


def main(argv=None):
    """Check the inputs, then write the summary table and the plot folder."""
    parser = argparse.ArgumentParser(
        description="Set every sample's trait profile side by side, for a table and MultiQC.")
    parser.add_argument("--profiles", nargs="+", required=True,
                        help="Every sample's {sample}_trait_profile.tsv.")
    parser.add_argument("--provenance", nargs="+", required=True,
                        help="Every sample's {sample}_trait_provenance.json.")
    parser.add_argument("--rules", required=True,
                        help="trait_rules_v1.tsv, the catalogue the calls were made with.")
    parser.add_argument("--out-summary", required=True,
                        help="The table of called / partial / not called.")
    parser.add_argument("--out-mqc", required=True,
                        help="The folder for the MultiQC plot files.")
    args = parser.parse_args(argv)

    # 1. Read and check everything before writing anything.
    trait_ids, grid, path_of_sample = build_grid(args.profiles)
    rules = read_rules(args.rules, trait_ids)
    check_provenance(args.provenance, path_of_sample, args.rules)
    # 2. The table, then the plots of the same calls.
    write_summary(args.out_summary, trait_ids, grid)
    write_multiqc_folder(args.out_mqc, multiqc_sections(trait_ids, grid, rules))
    print(f"Trait summary: {len(grid)} sample(s) x {len(trait_ids)} trait(s) written to "
          f"{args.out_summary} and {args.out_mqc}.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as error:
        sys.exit(f"[BacFlux] traits_summary: {error}")

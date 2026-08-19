#!/usr/bin/env python3
"""Score the ICE caller against ICEberg's curated coordinates (spec §8 Phase 7).

WHAT IS BEING MEASURED
    Two numbers per curated element, exactly as the spec asks:
      RECALL          did we call ANY element overlapping the curated interval?
      BOUNDARY OFFSET how far off were our start and end, in bp?

    Plus what we called it (ice / ime / island / conjugative_region), how the
    boundary was derived (tRNA-anchored att, de novo att, or not at all), and the
    confidence tier — because a correct call at low confidence and a correct call
    at high confidence are different results.

WHY THE TWO CLASSES ARE REPORTED SEPARATELY
    See classify_pilot.py. Seven of the eighteen accessions ARE the element
    rather than containing it, so for those the boundary answer is trivially "the
    whole record" and there is no flanking sequence for an att site to live in.
    Averaging boundary error across both classes would report a flattering number
    that means nothing. Detection is scored on all entries; boundary offset only
    on the eleven with real chromosomal context.

WHAT THIS DOES NOT MEASURE
    Agreement with ICEberg is not the same as truth. ICEberg's boundaries are
    themselves predictions in many cases — several entries are simply "the whole
    deposited record". A few hundred bp of disagreement is not automatically our
    error, and a perfect match on a standalone deposit is not automatically a
    success. Read the per-element table, not just the summary.

INPUT : pilot_set_classified.tsv, work/{sample}/ice_elements.tsv
OUTPUT: results_per_element.tsv + a summary printed to stdout
"""
import os
import sys

BENCH = os.path.dirname(os.path.abspath(__file__))

# A call counts as finding the element if it overlaps at all. Deliberately
# permissive: the question "did the machinery-based caller notice this element"
# is separate from "did it get the edges right", and the boundary offset columns
# answer the second one. A reciprocal-overlap criterion would conflate the two.
MIN_OVERLAP_BP = 1


def read_tsv(path):
    """Read a TSV into a list of dicts. Returns [] when the file is absent, which
    is a real result here — it means the caller produced no elements."""
    if not os.path.isfile(path):
        return []
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def _median(values):
    """Median of a sorted list, averaging the middle pair at even length.

    Written out rather than imported from statistics so the even-length
    behaviour is visible at the call site: the earlier sorted[len//2] form
    silently reported the UPPER of the two middle values, which overstated every
    even-n median in this benchmark.
    """
    ordered = sorted(values)
    n = len(ordered)
    if not n:
        return 0.0
    middle = n // 2
    if n % 2:
        return float(ordered[middle])
    return (float(ordered[middle - 1]) + float(ordered[middle])) / 2.0


def overlap_bp(a_start, a_end, b_start, b_end):
    """Overlap of two closed intervals, 0 when they do not touch."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def main():
    # Which classified pilot to score, and where to write the per-element
    # table. Defaults to the ICE set so existing invocations are unchanged.
    truth_path = (sys.argv[1] if len(sys.argv) > 1
                  else os.path.join(BENCH, "pilot_set_classified.tsv"))
    results_name = (sys.argv[2] if len(sys.argv) > 2 else "results_per_element.tsv")
    truth = read_tsv(truth_path)
    if not truth:
        sys.exit(f"{truth_path} not found — run classify_pilot.py first")

    results = []
    for entry in truth:
        accession = entry["accession"]
        sample = accession.replace(".", "_")
        true_start, true_end = int(entry["start"]), int(entry["end"])

        calls = read_tsv(os.path.join(BENCH, "work", sample, "ice_elements.tsv"))

        # Best call = the one overlapping the curated interval most. Ties are not
        # worth breaking carefully; if two calls overlap similarly the element was
        # split, and n_calls_overlapping below is what tells you that.
        best, best_overlap = None, 0
        n_overlapping = 0
        for call in calls:
            try:
                call_start, call_end = int(call["start"]), int(call["end"])
            except (ValueError, KeyError):
                continue
            shared = overlap_bp(true_start, true_end, call_start, call_end)
            if shared >= MIN_OVERLAP_BP:
                n_overlapping += 1
                if shared > best_overlap:
                    best, best_overlap = call, shared

        row = {
            "element": entry["element"],
            "organism": entry["organism"],
            "accession": accession,
            "class": entry["class"],
            "true_start": true_start,
            "true_end": true_end,
            "true_length_bp": true_end - true_start + 1,
            "n_calls_total": len(calls),
            "n_calls_overlapping": n_overlapping,
            "found": "TRUE" if best else "FALSE",
        }

        if best:
            call_start, call_end = int(best["start"]), int(best["end"])
            row.update({
                "called_start": call_start,
                "called_end": call_end,
                "called_length_bp": call_end - call_start + 1,
                # Signed so the direction of the error is visible: a negative
                # start offset means we started EARLIER than the curated call.
                "start_offset_bp": call_start - true_start,
                "end_offset_bp": call_end - true_end,
                "overlap_bp": best_overlap,
                # What fraction of the curated element we covered. The honest
                # summary number: 1.0 means we spanned it all.
                "recovered_fraction": round(best_overlap / row["true_length_bp"], 3),
                "mge_class": best.get("mge_class", "NA"),
                "boundary_method": best.get("boundary_method", "NA"),
                "att_length_bp": best.get("att_length_bp", "NA"),
                "att_trna": best.get("att_trna", "NA"),
                "confidence": best.get("confidence", "NA"),
                "machinery_intact": best.get("machinery_intact", "NA"),
                "anchor_classes": best.get("anchor_classes", "NA"),
            })
        else:
            for column in ("called_start", "called_end", "called_length_bp",
                           "start_offset_bp", "end_offset_bp", "overlap_bp",
                           "recovered_fraction", "mge_class", "boundary_method",
                           "att_length_bp", "att_trna", "confidence",
                           "machinery_intact", "anchor_classes"):
                row[column] = "NA"
        results.append(row)

    columns = list(results[0].keys())
    out_path = os.path.join(BENCH, results_name)
    with open(out_path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in results:
            handle.write("\t".join(str(row[c]) for c in columns) + "\n")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("=" * 78)
    print("PHASE 7 BENCHMARK — ICE detection vs ICEberg curated coordinates")
    print("=" * 78)

    for class_name, note in (("chromosomal", "detection + boundary measurable"),
                             ("standalone", "detection only; element IS the record")):
        subset = [r for r in results if r["class"] == class_name]
        if not subset:
            continue
        found = [r for r in subset if r["found"] == "TRUE"]
        print(f"\n{class_name.upper()}  (n={len(subset)}) — {note}")
        print(f"  detected: {len(found)}/{len(subset)}")

        if found:
            fractions = [float(r["recovered_fraction"]) for r in found]
            fractions.sort()
            # True median: average the middle PAIR when the count is even.
            # sorted[len//2] alone returns the upper of the two, which reported
            # the standalone-ICE median as 0.94 when the six values were
            # 0.289 0.346 0.401 0.936 0.982 0.994 and the median is 0.669 - a
            # number that then travelled into a commit message and a figure.
            median_recovered = _median(fractions)
            print(f"  median fraction of curated element recovered: {median_recovered:.2f}")

            classes = {}
            for r in found:
                classes[r["mge_class"]] = classes.get(r["mge_class"], 0) + 1
            print("  called as: " + ", ".join(f"{k}={v}" for k, v in sorted(classes.items())))

            methods = {}
            for r in found:
                methods[r["boundary_method"]] = methods.get(r["boundary_method"], 0) + 1
            print("  boundary from: " + ", ".join(f"{k}={v}" for k, v in sorted(methods.items())))

        # Boundary offsets only mean something where there is flanking context.
        if class_name == "chromosomal" and found:
            starts = sorted(abs(int(r["start_offset_bp"])) for r in found)
            ends = sorted(abs(int(r["end_offset_bp"])) for r in found)
            print(f"  |start offset| median {_median(starts):,.0f} bp  "
                  f"(min {starts[0]:,}, max {starts[-1]:,})")
            print(f"  |end offset|   median {_median(ends):,.0f} bp  "
                  f"(min {ends[0]:,}, max {ends[-1]:,})")

    missed = [r for r in results if r["found"] == "FALSE"]
    if missed:
        print(f"\nMISSED ({len(missed)}):")
        for r in missed:
            print(f"  {r['element']:<24} {r['organism']:<28} "
                  f"{r['class']:<12} calls_in_genome={r['n_calls_total']}")

    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()

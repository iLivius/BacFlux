# Benchmark definitions

**Nothing here is run by the workflow.** No Snakemake rule reads these files, and a
BacFlux run never touches this directory. It is a record, kept so that the numbers on
[Validation](../docs/mobilome/validation.md) can be traced back to what was actually
measured, and so the benchmark can be repeated later.

They are committed because they are the part that cannot be recreated. The genomes are
public accessions and can be downloaded again; the results can be recomputed. The
*design* — which genomes, chosen why, scored how — exists only here.

## `mobilome/` — the ICE/IME caller benchmark

| File | What it is |
|---|---|
| `ground_truth.tsv` | 1,677 curated elements extracted from an ICEberg 3.0 download. Every score is measured against these coordinates |
| `pilot_set.tsv` | The 18 ICE positives: a hand-made stratified draw from `ground_truth.tsv`, six of them the positive controls named in the design spec and the rest chosen to widen the set to 14 genera |
| `ime_pilot_set.tsv` | The 12 IME positives, deliberately stacked around the module's size floor — six sit **below** it and are expected misses |
| `negative_set.tsv` | The 12 genomes with no curated ICEberg entry |
| `pilot_set_classified.tsv`, `ime_pilot_classified.tsv` | The above, with the `class` column `classify_pilot.py` derives |
| `classify_pilot.py` | Splits each entry into `chromosomal` or `standalone` |
| `score.py` | The scorer: recall and boundary offset against the curated coordinates |

The `why` column in the three set files is the reasoning behind each choice, and it is
the reason these are kept rather than regenerated. Nothing derives it; it is judgement
about what would actually stress the caller, and it is recorded nowhere else.

`classify_pilot.py` exists because the 18 ICE entries are not 18 comparable tests.
ICEberg records an element as (accession, start, end), and that accession is sometimes a
whole chromosome and sometimes a standalone deposit of the element by itself. Those
support different measurements, and averaging them gives a recall figure that means
nothing — so detection is scored on all entries and boundary offset only on the ones
with real chromosomal context.

## `ebi_comparison/` — the head-to-head against the EBI pipeline

The conversion and scoring layer that made BacFlux's calls and the EBI Mobilome
Annotation Pipeline's calls comparable, so both could be scored by the same code. The
pipeline run itself is not kept — it is large and re-runnable — but this layer is where
the work went, and it is what
[the write-up](../docs/methods_ebi_comparison.md) rests on.

## Repeating a benchmark

The genomes are not here. Fetch them by accession from the set files, then run the
mobilome chain and score:

```bash
python benchmark/mobilome/classify_pilot.py      # if the classified files need rebuilding
python benchmark/mobilome/score.py benchmark/mobilome/pilot_set_classified.tsv results.tsv
```

Both scripts take paths as arguments and hold no BacFlux-specific state.

!!! note
    `ground_truth.tsv` came from ICEberg 3.0, which is unversioned and carries no reuse
    statement. A future download could differ from this copy with no way to tell, which
    is the other reason the file is kept rather than re-fetched on demand.

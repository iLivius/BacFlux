# Helper scripts — how to find the rule a script belongs to

Every script that a rule runs lives in a directory named after **the rule file that
runs it**. The directory name and the `.smk` filename are the same, so the trace-back
is by sight and needs no grepping:

| script directory | run by |
|---|---|
| `10_decontam/` | `workflow/rules/shared/10_decontam.smk` |
| `12_medaka_check/` | `workflow/rules/shared/12_medaka_check.smk` |
| `15_replicons/` | `workflow/rules/shared/15_replicons.smk` |
| `50_amr/` | `workflow/rules/shared/50_amr.smk` |
| `60_plasmid/` | `workflow/rules/shared/60_plasmid.smk` |
| `80_mobilome/` | `workflow/rules/shared/80_mobilome.smk` |

**A script sitting loose in this directory, not in a numbered one, is not part of a
workflow run.** It is a maintenance tool you run by hand. Right now that is
`clean_workdir.sh`.

## Two scripts inside `80_mobilome/` are also hand-run

`generate_gtdb_organism_table.py` and `check_gtdb_organism_table.py` are not invoked by
any rule. They build and validate the two committed tables
(`gtdb_organism_equivalences.tsv`, `gtdb_organism_genus_rules.tsv`) that
`gtdb_amrfinder_organism.py` reads at runtime, and they are re-run when a new GTDB
release comes out — not once per sample.

They stay here rather than moving in with `clean_workdir.sh` because they `import
gtdb_amrfinder_organism`, which *is* a workflow script. Splitting them off would break
that import to gain a tidier listing, which is a bad trade.

## Tests sit next to the script they test

`50_amr/card_mapping_report.py` is tested by `50_amr/test_card_mapping_report.py`, and
so on for every script. That is not only tidiness: because the two are in the same
directory, both pytest and a direct `python test_whatever.py` already have the module
on the import path, so no test needs any `sys.path` juggling to find its subject.

Run the whole suite from the repository root:

```bash
pytest workflow/scripts -q
```

## Paths are declared once

No rule spells out a script path. Each one is a constant in
`workflow/rules/shared/00_common.smk` (`SELECT_TAXONOMY_SCRIPT`, `CARD_REPORT_SCRIPT`,
`COLOCALISE_SCRIPT`, …), so moving a script means editing one line, not hunting through
the rule files. If you move something, run

```bash
python miscellaneous/check_comment_references.py
```

which fails on any path, rule name, config key or constant that a comment mentions but
that no longer exists.

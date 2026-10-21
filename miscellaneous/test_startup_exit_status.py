#!/usr/bin/env python3
"""Guard test: a start-up check must never stop the run with exit status 0.

BacFlux checks its input before any job runs: the sample files, their names,
the Illumina mates, the config values. When a check fails, the run must stop
with exit status 1, so that a wrapper script, a CI job or a cluster submission
script sees the failure. Python's sys.exit(0), and a bare sys.exit(), stop with
status 0, which those tools read as success.

The checks live in the Snakefile and the rule files, which run while Snakemake
reads the workflow. This test reads those files as plain text and fails if
either form appears outside a comment. The right form is
sys.exit("[BacFlux] message"), which prints the message and exits with 1.

Needs no Snakemake run. Run with:
    pytest miscellaneous/test_startup_exit_status.py
"""

import glob
import os
import re
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# sys.exit(0) or sys.exit() with nothing inside, allowing spaces.
STATUS_ZERO_EXIT = re.compile(r"sys\.exit\(\s*0?\s*\)")


def workflow_files():
    """The Snakefile and every rule file: the code that runs at start-up."""
    files = [os.path.join(REPO, "workflow", "Snakefile")]
    files += sorted(glob.glob(os.path.join(REPO, "workflow", "rules", "**", "*.smk"),
                              recursive=True))
    return files


def status_zero_exits(path):
    """Return "file:line: text" for every status-0 exit outside a comment."""
    found = []
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            code = line.split("#", 1)[0]
            if STATUS_ZERO_EXIT.search(code):
                relative = os.path.relpath(path, REPO)
                found.append(f"{relative}:{number}: {line.strip()}")
    return found


class TestStartupExitStatus(unittest.TestCase):

    def test_workflow_files_are_found(self):
        # Guard on the guard: if the layout moves, this test must not pass by
        # reading nothing.
        files = workflow_files()
        self.assertTrue(os.path.exists(files[0]), files[0])
        self.assertGreater(len(files), 10)

    def test_no_start_up_check_exits_with_status_zero(self):
        found = []
        for path in workflow_files():
            found += status_zero_exits(path)
        self.assertEqual(found, [], "start-up checks that exit with status 0:\n" + "\n".join(found))

    def test_the_pattern_catches_both_forms(self):
        self.assertTrue(STATUS_ZERO_EXIT.search("sys.exit(0)"))
        self.assertTrue(STATUS_ZERO_EXIT.search("sys.exit( )"))
        self.assertFalse(STATUS_ZERO_EXIT.search('sys.exit("[BacFlux] message")'))


if __name__ == "__main__":
    unittest.main()

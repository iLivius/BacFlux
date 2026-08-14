"""Unit tests for medaka_model_check.py — the pure name-handling logic, plus
main()'s branching with the two Medaka helpers monkeypatched.

Only list_available_models() and resolve_from_reads() shell out to the tool, and
every test that reaches one of them swaps in a stub, so this file runs with no
Medaka install, no basecalled reads and no model download.

Why the branching is worth pinning down: rule check_medaka_model runs the script
right after read filtering and the assembler waits on its output, so what main()
decides here is whether a long-read run stops in seconds or spends an hour in
Flye first and only then discovers the model is unusable.

Run: pytest workflow/scripts/12_medaka_check/test_medaka_model_check.py
"""

import medaka_model_check as mc


# ── The model menu a real Medaka install prints ─────────────────────────────
# A small but representative slice of a real `medaka tools list_models`, covering
# the naming shapes the tokenizer must handle: old gNNN names with/without a
# device token, new vX.Y.Z names with pore/speed, and variant/snp models.
AVAILABLE = [
    "r103_fast_g507", "r103_hac_g507", "r103_sup_g507",
    "r103_hac_variant_g507",
    "r941_min_fast_g507", "r941_min_hac_g507", "r941_min_sup_g507",
    "r941_min_hac_variant_g507", "r941_min_hac_snp_g507",
    "r941_prom_fast_g507", "r941_prom_hac_g507", "r941_prom_sup_g507",
    "r1041_e82_400bps_fast_g615", "r1041_e82_400bps_hac_v5.2.0",
    "r1041_e82_400bps_sup_v5.2.0", "r1041_e82_400bps_sup_variant_v5.0.0",
]


# ── Splitting a model name into flowcell, device, accuracy and version ──────
# Medaka model names are underscore-joined but not fixed-width: some carry a
# min/prom device token and some do not, and the version tag is either gNNN or
# vX.Y.Z. tokenize() therefore classifies token by token, and these tests cover
# each shape the real menu above contains. Everything the tokenizer does not
# recognise (pore e82, speed 400bps, network suffixes) stays in the name but is
# not matched on.

def test_tokenize_old_style_with_device():
    t = mc.tokenize("r941_min_hac_g507")
    assert t["flowcell"] == "r941"
    assert t["device"] == "min"
    assert t["accuracy"] == "hac"
    assert t["version"] == "g507"
    assert t["is_variant"] is False


def test_tokenize_no_device_token():
    # r103 names carry no min/prom token; device must come back None, not crash.
    t = mc.tokenize("r103_fast_g507")
    assert t["flowcell"] == "r103"
    assert t["device"] is None
    assert t["accuracy"] == "fast"
    assert t["version"] == "g507"


def test_tokenize_new_style_version_and_pore():
    t = mc.tokenize("r1041_e82_400bps_sup_v5.2.0")
    assert t["flowcell"] == "r1041"
    assert t["accuracy"] == "sup"
    # First version-like token (v5.2.0) wins; pore/speed are ignored for matching.
    assert t["version"] == "v5.2.0"


def test_tokenize_flags_variant_and_snp():
    assert mc.tokenize("r941_min_hac_variant_g507")["is_variant"] is True
    assert mc.tokenize("r941_min_hac_snp_g507")["is_variant"] is True


# ── Narrowing the menu when the configured model is wrong ───────────────────
# The whole value of failing early is lost if the error just says "not available"
# against a list of 80 names. suggest() keeps only the models sharing whatever
# axes the bad name did get right, so the fix is usually the single line above or
# below. variant/snp models are stripped first: they are for variant CALLING, and
# offering one as a polishing model would be a wrong answer, not a near miss.

def test_consensus_models_drops_variant_and_snp():
    kept = mc.consensus_models(AVAILABLE)
    assert "r941_min_hac_g507" in kept
    assert "r941_min_hac_variant_g507" not in kept
    assert "r941_min_hac_snp_g507" not in kept
    assert "r1041_e82_400bps_sup_variant_v5.0.0" not in kept


def test_suggest_narrows_on_a_version_typo():
    # A plausible typo: right flowcell/device/accuracy, wrong version tag.
    text = mc.suggest("r941_min_hac_g999", AVAILABLE)
    assert "r941_min_hac_g507" in text          # the obvious correct pick is shown
    assert "r941_prom_hac_g507" not in text     # different device — excluded
    assert "r103_fast_g507" not in text         # different flowcell — excluded
    assert "variant" not in text                # variant models never suggested


def test_suggest_falls_back_to_full_table_when_unrecognisable():
    text = mc.suggest("totally-bogus-name", AVAILABLE)
    # No recognisable axis — show the whole consensus menu rather than nothing.
    assert "r941_min_hac_g507" in text
    assert "r1041_e82_400bps_sup_v5.2.0" in text


def test_suggest_narrows_on_flowcell_only():
    # Only the flowcell is recognisable — narrow by flowcell alone.
    text = mc.suggest("r1041_somethingwrong", AVAILABLE)
    assert "r1041_e82_400bps_hac_v5.2.0" in text
    assert "r941_min_hac_g507" not in text


# ── main(): an explicit model out of the config ─────────────────────────────
# The user set parameters.<mode>.medaka_model to a name. main() writes the model
# to --out on success and rule long_read_consensus reads it back, so the model is
# resolved once for the run rather than a second time inside the polishing rule.
# An invalid name must exit 1 BEFORE the assembler starts, and must write no
# --out file at all: that file IS the model name Medaka will be handed, so there
# is no such thing as a half-valid one.

def test_main_explicit_valid(tmp_path, monkeypatch):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g507", "--reads", "x.fastq",
                  "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == "r941_min_hac_g507"


def test_main_explicit_local_file_accepted_as_is(tmp_path, monkeypatch):
    # A path to a real file on disk is accepted without list validation
    # (Medaka -m accepts a file). Preserves v1's escape hatch.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    model_file = tmp_path / "custom_model.tar.gz"
    model_file.write_text("x")
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", str(model_file), "--reads", "x.fastq", "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == str(model_file)


def test_main_explicit_invalid_fails_with_suggestion(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g999", "--reads", "x.fastq",
                  "--out", str(out)])
    assert rc == 1
    assert not out.exists()
    assert "r941_min_hac_g507" in capsys.readouterr().err


# ── main(): auto mode, where Medaka answers with a PATH ─────────────────────
# Auto mode means the user left medaka_model empty and Medaka infers the model
# from the basecaller tag in the filtlong FASTQ headers.
#
# The realistic thing `medaka tools resolve_model` prints: a PATH to the model
# file, NOT a bare name — so it is NEVER a member of the (bare-name) AVAILABLE
# list. Auto mode must accept it anyway (it is a valid `-m` argument). This is
# the regression the adversarial review caught: an earlier version required the
# resolved value to be in AVAILABLE, which failed every real auto run.
RESOLVED_PATH = "/opt/medaka/data/r941_min_hac_g507_model_pt.tar.gz"


def test_main_auto_success_accepts_a_path_not_in_the_name_list(tmp_path, monkeypatch):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: RESOLVED_PATH)
    assert RESOLVED_PATH not in AVAILABLE          # guard: it is a path, not a name
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "", "--reads", "x.fastq", "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == RESOLVED_PATH


def test_main_auto_failure_shows_menu(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: None)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "", "--reads", "x.fastq", "--out", str(out)])
    assert rc == 1
    assert not out.exists()
    err = capsys.readouterr().err
    assert "could not infer" in err
    assert "r941_min_hac_g507" in err


# ── main(): the opt-in fallback from a bad explicit model to auto ───────────
# parameters.<mode>.medaka_model_fallback_auto, default off. It only ever loosens
# the failure: an invalid explicit model may be replaced by an auto-inferred one,
# loudly, and if auto cannot help either the run still stops. The warning is part
# of the contract — silently polishing with a model the user did not choose is
# the thing this whole script exists to prevent.

def test_main_middle_option_falls_back_to_auto(tmp_path, monkeypatch, capsys):
    # Invalid explicit model, fallback on, auto succeeds — use auto and warn.
    # Auto returns a PATH (as real medaka does), which must be accepted.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: RESOLVED_PATH)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g999", "--reads", "x.fastq",
                  "--fallback-to-auto", "true", "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == RESOLVED_PATH
    assert "WARNING" in capsys.readouterr().err


def test_main_middle_option_but_auto_also_fails(tmp_path, monkeypatch):
    # Invalid explicit, fallback on, but auto can't help either — still fail.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: None)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g999", "--reads", "x.fastq",
                  "--fallback-to-auto", "true", "--out", str(out)])
    assert rc == 1
    assert not out.exists()

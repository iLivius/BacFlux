"""Unit tests for medaka_model_check.py — the pure string logic only (no Medaka
call), plus main()'s branching with the two Medaka helpers monkeypatched.

Run: pytest workflow/scripts/test_medaka_model_check.py
"""

import medaka_model_check as mc


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
    assert "r941_prom_hac_g507" not in text     # different device -> excluded
    assert "r103_fast_g507" not in text         # different flowcell -> excluded
    assert "variant" not in text                # variant models never suggested


def test_suggest_falls_back_to_full_table_when_unrecognisable():
    text = mc.suggest("totally-bogus-name", AVAILABLE)
    # No recognisable axis -> show the whole consensus menu.
    assert "r941_min_hac_g507" in text
    assert "r1041_e82_400bps_sup_v5.2.0" in text


def test_suggest_narrows_on_flowcell_only():
    # Only the flowcell is recognisable -> narrow by flowcell alone.
    text = mc.suggest("r1041_somethingwrong", AVAILABLE)
    assert "r1041_e82_400bps_hac_v5.2.0" in text
    assert "r941_min_hac_g507" not in text


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


def test_main_auto_success(tmp_path, monkeypatch):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: "r941_min_hac_g507")
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "", "--reads", "x.fastq", "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == "r941_min_hac_g507"


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


def test_main_middle_option_falls_back_to_auto(tmp_path, monkeypatch, capsys):
    # Invalid explicit model + fallback on + auto succeeds -> use auto, warn.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: "r941_min_hac_g507")
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g999", "--reads", "x.fastq",
                  "--fallback-to-auto", "true", "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == "r941_min_hac_g507"
    assert "WARNING" in capsys.readouterr().err


def test_main_middle_option_but_auto_also_fails(tmp_path, monkeypatch):
    # Invalid explicit + fallback on, but auto can't help either -> still fail.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: None)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g999", "--reads", "x.fastq",
                  "--fallback-to-auto", "true", "--out", str(out)])
    assert rc == 1
    assert not out.exists()

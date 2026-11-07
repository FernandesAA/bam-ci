from pathlib import Path

from bamci.config import CASE_PRESETS, Config


def test_presets_present():
    assert {"EIT", "FNL", "GFS", "ERA5"} <= set(CASE_PRESETS)
    assert CASE_PRESETS["FNL"]["kmaxq"] == 21


def test_output_name_matches_production():
    cfg = Config(case="EIT", date="2026100100", mend_out=126, kmax_out=42, data_out=Path("/tmp"))
    assert cfg.output_path().name == "GANLEIT2026100100S.unf.TQ0126L042"


def test_sigma_suffix_legacy():
    cfg = Config(case="EIT", date="2026100100", sigma_suffix=True, sigma_level=False,
                 data_out=Path("/tmp"))
    assert cfg.output_path().name.endswith("_hibrid")


def test_derived_constants():
    cfg = Config()
    assert abs(cfg.ctvr - (461.50 / 287.05 - 1.0)) < 1e-12
    assert abs(cfg.rocp - 287.05 / 1004.6) < 1e-12


def test_check_rejects_bad_date():
    assert Config(date="2026").check()

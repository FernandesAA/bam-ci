"""Resolucao configuravel: truncacao (mend) e camadas (kmax) por argumento.

Cobre os formatos que o usuario pediu explicitamente — TQ0213L64, TQ0666L42 —
e os que ja existiam (TQ0126L042), alem das grades intermediarias.
"""

from __future__ import annotations

import numpy as np
import pytest

from bamci.config import Config
from bamci.ic_format import ICHeader, ic_file_size, read_ic, write_ic
from bamci.levels import fft_friendly_size, gaussian_latitudes, grid_shape
from bamci.resolution import (
    describe,
    format_resolution,
    parse_resolution,
    resolve,
)
from bamci.spectral import analysis, legendre_columns, legendre_norm, synthesis
from bamci.sst import default_mask, model_grid

# resolucoes de referencia: (label, grid, mend, kmax, nlon, nlat)
CASES = [
    ("TQ0126L042", "Q", 126, 42, 384, 192),
    ("TQ0213L064", "Q", 213, 64, 640, 320),
    ("TQ0254L042", "Q", 254, 42, 768, 384),
    ("TQ0666L042", "Q", 666, 42, 2000, 1000),
    ("TL0095L042", "L", 95, 42, 288, 144),
]


# -- leitura dos argumentos -------------------------------------------------
@pytest.mark.parametrize("text,expected", [
    ("TQ0213L64", ("Q", 213, 64)),
    ("TQ0666L42", ("Q", 666, 42)),
    ("TQ0126L042", ("Q", 126, 42)),
    ("T213L42", ("", 213, 42)),
    ("213L42", ("", 213, 42)),
    ("q213l64", ("Q", 213, 64)),
    ("TQ0213", ("Q", 213, None)),
    ("666", ("", 666, None)),
])
def test_parse_resolution(text, expected):
    assert parse_resolution(text) == expected


@pytest.mark.parametrize("bad", ["", "TX213", "TQ", "TQ0213L", "TQ0213L0", "T-213", "TQ0213L64x"])
def test_parse_resolution_rejects(bad):
    with pytest.raises(ValueError):
        parse_resolution(bad)


def test_format_roundtrip():
    for label, grid, mend, kmax, _, _ in CASES:
        assert format_resolution(grid, mend, kmax) == label
        assert parse_resolution(label) == (grid, mend, kmax)


def test_resolve_explicit_wins_over_text():
    # -r TQ0213 --kmax 91: o kmax explicito vence o que veio no texto
    assert resolve(resolution="TQ0213", kmax=91) == ("Q", 213, 91)
    assert resolve(resolution="TQ0126L042") == ("Q", 126, 42)
    assert resolve(mend=254) == ("Q", 254, None)
    assert resolve(resolution="TL95L42")[0] == "L"


def test_resolve_rejects_bad_grid():
    with pytest.raises(ValueError):
        resolve(grid="X")


# -- grade ----------------------------------------------------------------
@pytest.mark.parametrize("label,grid,mend,kmax,nlon,nlat", CASES)
def test_grid_shape(label, grid, mend, kmax, nlon, nlat):
    assert grid_shape(mend) == (nlon, nlat)
    assert fft_friendly_size(3 * mend + 1) == nlon


def test_grid_matches_production_t126():
    assert grid_shape(126) == (384, 192)


# -- tamanho do arquivo ----------------------------------------------------
def test_t126_size_unchanged():
    assert ic_file_size(126, 42) == 11_055_828


@pytest.mark.parametrize("label,grid,mend,kmax,nlon,nlat", CASES)
def test_size_formula(label, grid, mend, kmax, nlon, nlat):
    coeff = (mend + 1) * (mend + 2)
    expected = 4 * (11 + 2 * kmax + 2) + (2 + 4 * kmax) * 4 * (coeff + 2)
    assert ic_file_size(mend, kmax) == expected


def test_describe_reports_resolution():
    info = describe("Q", 213, 64)
    assert info["label"] == "TQ0213L064"
    assert (info["nlon"], info["nlat"]) == (640, 320)
    assert info["bytes"] == ic_file_size(213, 64)


# -- Config ---------------------------------------------------------------
def test_config_default_is_t126l42():
    cfg = Config()
    assert cfg.resolution_label == "TQ0126L042"
    assert cfg.output_path().name == "GANLERA52026100100S.unf.TQ0126L042"


@pytest.mark.parametrize("label,grid,mend,kmax,nlon,nlat", CASES)
def test_config_accepts_resolution(label, grid, mend, kmax, nlon, nlat):
    del nlon, nlat
    cfg = Config.from_resolution(label, case="GFS", date="2026110500")
    assert (cfg.mend_out, cfg.kmax_out, cfg.grid) == (mend, kmax, grid)
    assert cfg.output_path().name == "GANLGFS2026110500S.unf.%s" % label
    assert not cfg.check()


def test_config_set_resolution_partial():
    cfg = Config().set_resolution("TQ0213")
    assert (cfg.mend_out, cfg.kmax_out) == (213, 42)      # kmax preservado
    cfg.set_resolution(kmax=91)
    assert (cfg.mend_out, cfg.kmax_out) == (213, 91)


# -- transformada espectral ----------------------------------------------
@pytest.mark.parametrize("mend", [21, 63, 126])
def test_columns_match_legendre_norm(mend):
    nlon, nlat = grid_shape(mend)
    mu = np.sin(np.radians(gaussian_latitudes(nlat)))
    ref = legendre_norm(mend, mu)
    for m, col in legendre_columns(mend, mu):
        np.testing.assert_allclose(col, ref[m, m:, :], rtol=0, atol=1e-12)


@pytest.mark.parametrize("mend", [126, 213, 666])
def test_normalisation_identity_high_truncation(mend):
    """sum_meia 2 w pln^2 = 1 — vale em qualquer truncacao, sem overflow."""
    from bamci.spectral import gauss_weights

    nlon, nlat = grid_shape(mend)
    mu = np.sin(np.radians(gaussian_latitudes(nlat)))[:nlat // 2]
    w = gauss_weights(nlat)[:nlat // 2]
    for m, col in legendre_columns(mend, mu):
        assert np.isfinite(col).all(), "pln estourou em m=%d" % m
        norm = (2.0 * w * col ** 2).sum(axis=1)
        np.testing.assert_allclose(norm, 1.0, atol=1e-9)


@pytest.mark.parametrize("mend", [42, 126, 213])
def test_analysis_synthesis_roundtrip(mend):
    """Roundtrip exato num campo band-limitado (a grade reproduz a truncacao)."""
    nlon, nlat = grid_shape(mend)
    coef = np.zeros((mend + 1, mend + 1), dtype=complex)
    coef[2, 7] = 1.5 - 0.25j
    coef[0, 0] = 3.0
    field = synthesis(coef, mend, nlat, nlon)
    back = synthesis(analysis(field, mend), mend, nlat, nlon)
    np.testing.assert_allclose(back, field, atol=1e-8)


# -- SST ------------------------------------------------------------------
@pytest.mark.parametrize("mend,nlat", [(126, 192), (213, 320), (254, 384), (666, 1000)])
def test_sst_grid_and_mask_name(mend, nlat):
    assert model_grid(mend)[0] == nlat
    assert default_mask(mend).endswith("ModelLandSeaMask.G%05d" % nlat)


def test_sst_mask_keeps_historic_name():
    assert default_mask(126).endswith("ModelLandSeaMask.G00192")


# -- arquivo de IC em resolucao diferente ---------------------------------
@pytest.mark.parametrize("mend,kmax", [(42, 28), (126, 64)])
def test_write_read_ic_other_resolutions(tmp_path, mend, kmax):
    ncoeff = (mend + 1) * (mend + 2)
    rng = np.random.default_rng(0)
    delta = np.full(kmax, 1.0 / kmax)
    header = ICHeader(forecast_day=0, time_of_day=0.0, date_inp=(1, 10, 2026, 0),
                      date_cur=(1, 10, 2026, 0),
                      sigma_inter=np.concatenate([[0.0], np.cumsum(delta)]),
                      sigma_layers=0.5 * (np.concatenate([[0.0], np.cumsum(delta)])[1:]
                                          + np.concatenate([[0.0], np.cumsum(delta)])[:-1]))
    path = tmp_path / ("ic_T%dL%d.unf" % (mend, kmax))
    write_ic(
        path, header,
        rng.standard_normal(ncoeff), rng.standard_normal(ncoeff),
        rng.standard_normal((ncoeff, kmax)), rng.standard_normal((ncoeff, kmax)),
        rng.standard_normal((ncoeff, kmax)), rng.standard_normal((ncoeff, kmax)),
    )
    assert path.stat().st_size == ic_file_size(mend, kmax)
    data = read_ic(path)
    assert data["kmax"] == kmax
    assert data["ncoeff"] == ncoeff
    assert data["tvir"].shape == (ncoeff, kmax)

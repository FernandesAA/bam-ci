"""Regressao do bug relatado: `bamci run --case ERA5 --grib ...` falhava.

O pipeline lia qualquer GRIB com a tabela do NCEP (shortName "r") e um arquivo
ERA5 (shortName "q") morria com "niveis ausentes para rh: [todos os niveis]".
O leitor agora vem do --case (era5 para ERA5/ERA-Interim, nl para FNL/GFS).
"""

from __future__ import annotations

import numpy as np
import pytest

import bamci.pipeline as pipeline
from bamci.config import Config
from bamci.gribio import RawField
from bamci.ic_format import ic_file_size, read_ic

NLAT, NLON = 32, 64
LEVELS = (1000.0, 850.0, 500.0)


def _era5_fields(levels=LEVELS):
    """Conjunto minimo de mensagens com os shortNames do ERA5."""
    lat = np.linspace(90.0, -90.0, NLAT)
    lon = np.linspace(0.0, 360.0, NLON, endpoint=False)
    phi = np.radians(lat)[:, None]
    lam = np.radians(lon)[None, :]

    def mk(name, level_type, level, values, units=""):
        return RawField(short_name=name, level_type=level_type, level=level,
                        values=np.asarray(values, dtype=float), latitudes=lat,
                        longitudes=lon, units=units)

    fields = []
    for k, p in enumerate(levels):
        base = 288.0 - 40.0 * k + 6.0 * np.cos(phi) * np.cos(2 * lam)
        fields.append(mk("t", "isobaricInhPa", p, base, "K"))
        fields.append(mk("q", "isobaricInhPa", p, np.full((NLAT, NLON), 0.006), "kg kg**-1"))
        fields.append(mk("u", "isobaricInhPa", p, 10.0 * np.cos(phi) + 0 * lam, "m/s"))
        fields.append(mk("v", "isobaricInhPa", p, np.full((NLAT, NLON), -1.5), "m/s"))
    sp = 1013.0 - 100.0 * np.sin(phi) ** 2
    fields.append(mk("sp", "surface", 0.0, np.broadcast_to(sp, (NLAT, NLON)) * 100.0, "Pa"))
    fields.append(mk("msl", "meanSea", 0.0, np.full((NLAT, NLON), 101500.0), "Pa"))
    fields.append(mk("z", "surface", 0.0, np.full((NLAT, NLON), 500.0 * 9.80665), "m**2 s**-2"))
    fields.append(mk("sst", "surface", 0.0, np.full((NLAT, NLON), 290.0), "K"))
    fields.append(mk("ci", "surface", 0.0, np.zeros((NLAT, NLON))))
    return fields


@pytest.fixture()
def delta_sigma(tmp_path):
    path = tmp_path / "DeltaSigma.L003"
    np.savetxt(path, np.full((1, 3), 1.0 / 3.0), fmt="%.6f")
    return path


def _run(tmp_path, delta_sigma, case, monkeypatch, source=None):
    monkeypatch.setattr(pipeline, "read_grib", lambda _p: _era5_fields())
    cfg = Config(case=case, date="2026102300")
    cfg.set_resolution("TQ0042L003")
    cfg.data_out = tmp_path
    cfg.smooth_fields = False
    return pipeline.run(cfg, grib_path="inventado.grb",
                        delta_sigma_path=str(delta_sigma), source=source)


def test_era5_case_reads_q_not_rh(tmp_path, delta_sigma, monkeypatch):
    """O caso ERA5 tem de ser lido como ERA5, nao como NCEP."""
    res = _run(tmp_path, delta_sigma, "ERA5", monkeypatch)
    assert res.path.stat().st_size == ic_file_size(42, 3)
    data = read_ic(res.path)
    assert data["kmax"] == 3
    assert np.isfinite(np.asarray(data["sphu"])).all()


def test_gfs_case_still_uses_ncep_reader(tmp_path, delta_sigma, monkeypatch):
    """Caso NCEP com arquivo NCEP continua usando o leitor 'nl'... e falha."""
    with pytest.raises(KeyError, match="umidade relativa"):
        _run(tmp_path, delta_sigma, "GFS", monkeypatch)


def test_source_override_beats_case(tmp_path, delta_sigma, monkeypatch):
    """--source era5 forca o leitor mesmo com --case GFS."""
    res = _run(tmp_path, delta_sigma, "GFS", monkeypatch, source="era5")
    assert res.path.exists()

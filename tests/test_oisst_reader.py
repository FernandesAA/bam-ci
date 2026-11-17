"""Leitura do OISST com xarray: janela de tempo sem exigir indice monotonico.

O `ds.sel(time=slice(...))` do xarray falha em arquivos cuja ordem de tempo nao
e crescente ("Cannot get left slice bound for non-monotonic index"). O reader
usa recorte posicional, entao a ordem do arquivo nao importa.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("xarray")

from bamci.sst import read_oisst  # noqa: E402


def _write_netcdf(path, times, sst_values, lat, lon, ice=True):
    """Grava um NetCDF minimo (motor scipy: nao exige netCDF4)."""
    import xarray as xr

    data = {
        "sst": (("time", "lat", "lon"), np.asarray(sst_values, dtype="float32")),
    }
    if ice:
        data["ice"] = (("time", "lat", "lon"), np.zeros_like(np.asarray(sst_values, dtype="float32")))
    ds = xr.Dataset(
        data,
        coords={"time": np.asarray(times, dtype="datetime64[ns]"),
                "lat": np.asarray(lat, dtype="float64"),
                "lon": np.asarray(lon, dtype="float64")},
    )
    ds.to_netcdf(path, engine="scipy")
    ds.close()
    return path


LAT = np.array([10.0, 0.0, -10.0])
LON = np.array([0.0, 90.0, 180.0, 270.0])


def _values(n_times, base=300.0):
    out = []
    for k in range(n_times):
        out.append(np.full((LAT.size, LON.size), base + k, dtype="float32"))
    return np.stack(out)


def test_reads_window_with_ordered_time(tmp_path):
    times = np.array(["2026-11-10", "2026-11-11", "2026-11-12"], dtype="datetime64[D]")
    p = _write_netcdf(tmp_path / "sst.nc", times, _values(3), LAT, LON)
    sst, ice, lat, lon = read_oisst(p, "20261112", days=2)
    assert sst.shape == (LAT.size, LON.size)
    # dias 11 e 12 -> media de 300+1 e 300+2
    np.testing.assert_allclose(sst, 301.5, atol=1e-4)
    assert ice.shape == sst.shape


def test_reads_window_with_non_monotonic_time(tmp_path):
    """Regressao: time fora de ordem quebrava ds.sel(time=slice(...))."""
    times = np.array(["2026-11-12", "2026-11-05", "2026-11-11"], dtype="datetime64[D]")
    values = _values(3)          # indice 0 -> 300 (12/11), 1 -> 301 (05/11), 2 -> 302 (11/11)
    p = _write_netcdf(tmp_path / "sst.nc", times, values, LAT, LON)
    sst, _, _, _ = read_oisst(p, "20261112", days=2)   # janela: 11/11 e 12/11
    # seleciona os indices 0 e 2 (300 e 302), na ordem do arquivo
    np.testing.assert_allclose(sst, (300.0 + 302.0) / 2.0, atol=1e-4)


def test_single_time_weekly_mean_uses_nearest(tmp_path):
    """Arquivo de media semanal (1 tempo): usa o tempo mais proximo, com aviso."""
    times = np.array(["2026-11-05"], dtype="datetime64[D]")
    p = _write_netcdf(tmp_path / "sst.nc", times, _values(1), LAT, LON)
    with pytest.warns(RuntimeWarning, match="mais proximo"):
        sst, _, _, _ = read_oisst(p, "20261112", days=7)
    np.testing.assert_allclose(sst, 300.0, atol=1e-4)


def test_celsius_is_converted_to_kelvin(tmp_path):
    times = np.array(["2026-11-12"], dtype="datetime64[D]")
    p = _write_netcdf(tmp_path / "sst.nc", times, _values(1, base=20.0), LAT, LON)
    sst, _, _, _ = read_oisst(p, "20261112", days=1)
    np.testing.assert_allclose(sst, 293.15, atol=1e-3)


def test_missing_ice_becomes_zeros(tmp_path):
    times = np.array(["2026-11-12"], dtype="datetime64[D]")
    p = _write_netcdf(tmp_path / "sst.nc", times, _values(1), LAT, LON, ice=False)
    sst, ice, _, _ = read_oisst(p, "20261112", days=1)
    assert np.all(ice == 0.0)
    assert sst.shape == ice.shape

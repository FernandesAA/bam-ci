import numpy as np
import pytest

from bamci.sst import (UNDEF, normalize_date, open_dataset, resolve_sst_output,
                      to_model_grid, write_sstweekly)


def test_to_model_grid_shape_and_order():
    lat = np.linspace(90, -90, 181)
    lon = np.arange(0, 360, 1.0)
    field = np.full((181, 360), 300.0)
    out = to_model_grid(field, lat, lon, nlat=48, nlon=96)
    assert out.shape == (48, 96)
    assert np.allclose(out[np.isfinite(out)], 300.0, atol=1e-6)


def test_land_becomes_undef():
    lat = np.linspace(90, -90, 19)
    lon = np.arange(0, 360, 10.0)
    field = np.full((19, 36), 300.0)
    field[:, :18] = np.nan          # metade "terra"
    out = to_model_grid(field, lat, lon, nlat=24, nlon=48)
    assert np.any(out <= UNDEF + 1)
    assert np.nanmax(out) == pytest.approx(300.0, abs=1e-6)


def test_write_two_records(tmp_path):
    ice = np.zeros((192, 384))
    sst = np.full((192, 384), 300.0)
    p = write_sstweekly(tmp_path / "SSTWeekly20261105.G00192", ice, sst)
    assert p.stat().st_size == 2 * 192 * 384 * 4
    raw = np.fromfile(p, dtype=">f4").reshape(2, 192, 384)
    assert raw[1, 0, 0] == pytest.approx(300.0)

def test_ice_from_sst_rule():
    from bamci.sst import ice_from_sst
    sst = np.array([270.0, 271.35, 272.6, 275.0])
    ice = ice_from_sst(sst)
    assert ice[0] == 1.0
    assert abs(ice[1] - 1.0) < 1e-12
    assert 0.0 < ice[2] < 1.0
    assert ice[3] == 0.0

# -- datas, saida e motor NetCDF -----------------------------------------
def test_normalize_date_accepts_common_forms():
    for text in ("20261112", "2026111200", "2026-11-12", "2026-11-12T00",
                 "202611120000", " 20261112 "):
        assert normalize_date(text) == "20261112"


@pytest.mark.parametrize("bad", ["", "2026", "abc", "20261312", "2026111"])
def test_normalize_date_rejects(bad):
    with pytest.raises(ValueError):
        normalize_date(bad)


def test_resolve_sst_output_directory(tmp_path):
    """--out como diretorio monta SSTWeekly<data>.G00<lat> dentro dele."""
    p = resolve_sst_output(str(tmp_path) + "/", "2026111200", 192)
    assert p.name == "SSTWeekly20261112.G00192"
    assert p.parent == tmp_path

    p2 = resolve_sst_output(tmp_path, "2026-11-12", 320)
    assert p2.name == "SSTWeekly20261112.G00320"


def test_resolve_sst_output_file(tmp_path):
    """--out como caminho de arquivo e respeitado (nao renomeia)."""
    target = tmp_path / "meu_sst.G00192"
    assert resolve_sst_output(target, "20261112", 192) == target


def test_resolve_sst_output_dir_requires_date(tmp_path):
    with pytest.raises(ValueError, match="--date"):
        resolve_sst_output(str(tmp_path) + "/", None, 192)


def test_open_dataset_without_xarray_is_actionable(monkeypatch):
    """Sem xarray a mensagem tem de dizer o que instalar."""
    import sys

    monkeypatch.setitem(sys.modules, "xarray", None)
    with pytest.raises(RuntimeError, match="netcdf"):
        open_dataset("qualquer.nc")


def test_open_dataset_without_engine_is_actionable(monkeypatch, tmp_path):
    """O erro real do xarray (backends sem dependencia) vira instrucao clara."""
    import sys
    import types

    fake = types.ModuleType("xarray")

    def _open(*_a, **_k):
        raise ValueError(
            "found the following matches with the input file in xarray's IO "
            "backends: ['netcdf4', 'h5netcdf']. But their dependencies may not "
            "be installed"
        )

    fake.open_dataset = _open
    monkeypatch.setitem(sys.modules, "xarray", fake)
    with pytest.raises(RuntimeError) as exc:
        open_dataset(tmp_path / "x.nc")
    msg = str(exc.value)
    assert "pip install netCDF4" in msg
    assert "h5netcdf" in msg

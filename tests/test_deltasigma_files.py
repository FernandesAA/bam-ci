"""Consistencia dos arquivos DeltaSigma.LNNN empacotados em `DeltaSigma/`.

Cada `DeltaSigma.LNNN` traz **NNN espessuras sigma** (uma por camada, formato
`(5F9.6)`), todas positivas e somando 1.

Tolerancia da soma: `1e-4` -- DeltaSigma.L064 da origem soma 0.999976
(arredondamento de 2.4e-5 no arquivo operacional). O leitor
`bamci.levels.read_delta_sigma` normaliza a soma, entao o pipeline nao e
afetado por esse residuo; o ultimo teste cobre exatamente isso.
"""

from pathlib import Path

import numpy as np
import pytest

from bamci.levels import read_delta_sigma

DELTASIGMA_DIR = Path(__file__).resolve().parents[1] / "DeltaSigma"
FILES = sorted(DELTASIGMA_DIR.glob("DeltaSigma.L*"))
EXPECTED = [9, 18, 28, 42, 64, 96, 160]
SUM_TOL = 1e-4
IDS = [p.name for p in FILES]


def _layers(name: str) -> int:
    return int(name.rsplit(".L", 1)[1])


def _values(path: Path) -> np.ndarray:
    return np.asarray([float(t) for t in path.read_text().split()], dtype=np.float64)


def test_deltasigma_dir_present():
    assert DELTASIGMA_DIR.is_dir(), DELTASIGMA_DIR
    assert IDS == [f"DeltaSigma.L{n:03d}" for n in EXPECTED]


@pytest.mark.parametrize("path", FILES, ids=IDS)
def test_n_layers_and_positive(path):
    values = _values(path)
    assert values.size == _layers(path.name), f"{path.name}: {values.size} valores"
    assert np.all(values > 0.0), f"{path.name}: valores nao positivos"


@pytest.mark.parametrize("path", FILES, ids=IDS)
def test_soma_um(path):
    total = _values(path).sum()
    assert abs(total - 1.0) <= SUM_TOL, f"{path.name}: soma {total:.6f}"


@pytest.mark.parametrize("path", FILES, ids=IDS)
def test_reader_normaliza_soma(path):
    """read_delta_sigma normaliza: a soma lida e 1 mesmo para o L064."""
    delta = read_delta_sigma(path)
    assert delta.size == _layers(path.name)
    assert abs(delta.sum() - 1.0) <= 1e-12

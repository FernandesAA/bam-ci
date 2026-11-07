"""Smoke tests da API publica de interpolacao vertical.

A implementacao fiel ao Fortran esta em tests/test_vertical_sigma.py; aqui so se
garante que a API exportada continua importavel e coerente.
"""

import numpy as np

from bamci import interp_logp, new_ps, new_sigma, vert_sigma_inter


def test_public_api_roundtrip():
    p_in = np.array([1000.0, 850.0, 700.0, 500.0])
    field = np.tile(np.log(p_in)[:, None, None], (1, 1, 1))
    p_out = np.array([950.0, 800.0, 600.0])[:, None, None]
    got = interp_logp(field, p_in, p_out)
    np.testing.assert_allclose(got[:, 0, 0], np.log([950.0, 800.0, 600.0]), atol=1e-10)


def test_new_sigma_and_new_ps_smoke():
    inter, layer = new_sigma(np.array([0.5, 0.5]))
    assert inter[0] == 1.0 and inter[-1] == 0.0
    assert layer.shape == (2,)
    ps = new_ps(np.full((1, 1), 1013.25), np.array([1000.0, 500.0]),
                np.full((2, 1, 1), 285.0), np.zeros((1, 1)))
    assert ps[0, 0] > 900.0

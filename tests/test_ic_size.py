"""Tamanho do arquivo de IC tem de fechar com o arquivo de producao.

GANLEIT2026100100S.unf.TQ126L042 medido na jaci: 11.055.828 bytes.
"""

from bamci.ic_format import ic_file_size


def test_size_matches_production_file():
    # arquivo medido na jaci: 11.055.828 bytes -> registro 1 com 95 palavras (sigma)
    assert ic_file_size(126, 42, "sigma") == 11_055_828


def test_sigma_and_hybrid_formulas():
    for mend, kmax in ((62, 28), (126, 42)):
        coeff = 4 * ((mend + 1) * (mend + 2) + 2)
        assert ic_file_size(mend, kmax, "sigma") == (11 + 2 * kmax + 2) * 4 + (2 + 4 * kmax) * coeff
        assert ic_file_size(mend, kmax, "hybrid") == (12 + 2 * kmax + 2) * 4 + (2 + 4 * kmax) * coeff

# Contribuindo

1. Crie um ambiente: `python -m venv .venv && .venv/bin/pip install -e ".[dev,grib,netcdf]"`.
2. Rode os testes: `.venv/bin/python -m pytest -q` (a suite roda sem rede).
3. Mantenha uma funcao pura por etapa numerica e um teste para cada uma.
4. Se alterar o formato dos arquivos de saida, atualize `tests/test_ic_format.py`
   (tamanho e registros) e o README (tabela de validacao).

Estilo: PT-BR nos comentarios e docstrings, nomes de codigo em ingles.

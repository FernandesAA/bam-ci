# Changelog

## 1.1.0

Resolucao configuravel. O pacote deixou de assumir TQ126L042.

- **Qualquer truncacao e qualquer numero de camadas**, por argumento:
  `-r TQ0213L64`, `-r TQ0666L42`, ou `--mend 254 --kmax 42`, ou `--grid L`.
  Novo modulo `bamci.resolution` (`parse_resolution`, `format_resolution`,
  `describe`) e `Config.set_resolution`/`Config.from_resolution`.
- **Novo subcomando `bamci res`**: mostra grade, numero de coeficientes e
  tamanho do arquivo antes de rodar (`bamci res -r TQ0213L64`).
- **Novo subcomando `bamci run`**: pipeline completo pela CLI, com
  `-r/--mend/--kmax/--grid`, `--grib`, `--delta-sigma`, `--out-dir`; confere o
  tamanho do arquivo contra a resolucao pedida.
- `--resolution` tambem em `check`, `compare`, `sst` e `pbs` (que repassa ao run).
- **SST acompanha a resolucao**: grade e nome da mascara derivados da truncacao
  (`default_mask`: T126 -> `G00192`, T213 -> `G00320`, T254 -> `G00384`,
  T666 -> `G01000`).
- **Legendre sem overflow**: a recorrencia crua estourava em T213+ (NaN). Agora
  `legendre_columns` gera uma coluna por ordem zonal com reescalonamento, o que
  tambem reduz o pico de memoria de O(mend^2 x nlat) para O(mend x nlat) por
  coluna (T666: ~3,6 GB -> ~3 MB por coluna).
- **Derivada de Legendre estavel**: `dpln_columns` usa a recorrencia
  normalizada em vez de dividir pela escala `exp(logn)` (que underflowava).
  Conferido contra `scipy.special.lpmv` (erro relativo 4e-10).
- `pipeline.run` aborta com mensagem clara quando o `DeltaSigma` nao casa com
  o `--kmax` pedido.
- 54 testes novos (`tests/test_resolution.py`), suite em 143 testes.

## 1.0.0

Primeira versao publica.

- **Condicao inicial espectral** do BAM (TQ126L042, coordenada sigma):
  leitura GRIB (ERA5/FNL/GFS) -> grade gaussiana -> interpolacao vertical ->
  transformada espectral -> arquivo UNFORMATTED big-endian.
- **SST semanal** no formato do modelo (`SSTWeeklyYYYYMMDD.G00192`), a partir de
  OISST v2.1 ou RTG/NCEP, com a mascara de terra do modelo.
- Validacao contra arquivos de producao do CPTEC/INPE (ver README).
- CLI `bamci` com `check`, `compare`, `download`, `sst`, `pbs`.

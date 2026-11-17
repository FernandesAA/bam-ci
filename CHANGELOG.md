# Changelog

## 1.1.3

`bamci sst` lia o eixo de tempo do OISST por rotulo e falhava em arquivos reais.

Sintoma: `Cannot get left slice bound for non-monotonic index with a missing
label np.datetime64('20261112-01-01')`.

Duas causas:

1. **`ds.sel(time=slice(...))` exige indice monotonico** — arquivos de media
   semanal baixados fora de ordem quebravam. Agora o recorte e **posicional**
   (`isel` sobre mascara booleana), que funciona em qualquer ordem.
2. **`np.datetime64("20261112", "D")` nao le como data**: o numpy interpreta os
   8 digitos como **ano** e devolve `20261112-01-01` (o rotulo estranho da
   mensagem de erro). A conversao passa por `datetime.strptime`, com o helper
   `_date64`.

Tambem: `time` ausente vira erro explicito ("arquivo sem dimensao de tempo") em
vez de `KeyError`; a busca do tempo mais proximo usa o mesmo caminho numerico.

- Testes: 161 -> 166 (`tests/test_oisst_reader.py` grava NetCDF de verdade via
  motor `scipy` e cobre: janela ordenada, janela fora de ordem, arquivo com um
  unico tempo, SST em Celsius e ausencia de gelo).

## 1.1.2

`bamci sst` usavel na pratica: mensagens de erro acionaveis e `--out` como diretorio.

- **Motor NetCDF ausente**: o erro do xarray ("found the following matches with
  the input file in xarray's IO backends ... their dependencies may not be
  installed") agora vira instrucao clara:
  `pip install netCDF4` (ou `h5netcdf`, ou `pip install 'bam-ci[netcdf]'`).
- **`--out` aceita diretorio**: terminando em `/` ou sendo diretorio existente,
  o arquivo e gravado como `SSTWeekly<AAAAMMDD>.G00<nlat>` dentro dele, com
  `mkdir -p` do caminho. Caminho de arquivo continua sendo respeitado.
- **`--date` normalizado**: `AAAAMMDD`, `AAAAMMDDHH` (hora descartada),
  `AAAA-MM-DD`, `AAAA-MM-DDTHH`; data inexistente (`20261312`) e recusada.
- **Arquivo de media semanal**: se a janela de 7 dias nao tiver dado (arquivo
  com um unico tempo), usa o tempo mais proximo com `RuntimeWarning` em vez de
  abortar.
- `--date` passa a ser obrigatorio apenas para `--source oisst`.
- `pyproject`: URLs do projeto apontam para o repositorio real (saíram os
  placeholders `<organizacao>`).
- Testes: 150 -> 161.

## 1.1.1

Correcao: **`bamci run` lia todo GRIB com a tabela do NCEP.**

Sintoma: com um arquivo do ERA5, `bamci run --case ERA5 --grib era5.grb ...`
abortava com `niveis ausentes para rh: [todos os niveis]`.

Causa: `pipeline.run` chamava `build_analysis(fields, "nl", levels)` com a
fonte fixa em `"nl"`. O ERA5 entrega umidade especifica (shortName `q`), nao
umidade relativa (shortName `r`), entao o campo nao era encontrado.

- `bamci.gribio.SOURCE_FOR_CASE` / `source_for_case(case)`: ERA5, ERA-40,
  ERA-Interim -> `"era5"`; FNL e GFS -> `"nl"`.
- `pipeline.run(..., source=None)` escolhe o leitor pelo `Config.case`.
- Novo `--source {era5,nl}` em `bamci run` para sobrepor o caso.
- Niveis ausentes agora sao **interpolados em log(p)** com `RuntimeWarning`
  (antes abortavam). Arquivos NCEP reais trazem menos niveis de umidade que de
  temperatura; campo completamente ausente continua sendo erro, com mensagem
  que indica o `--source` correto.
- `Config.check()` valida o `--case` e lista os validos (ex.: `ER5` -> erro
  claro em vez de falha tardia).
- Testes: 144 -> 150 (leitor por caso, regressao ERA5, interpolacao de niveis
  e `--source` sobrepondo o caso).

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

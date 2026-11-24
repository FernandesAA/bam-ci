# Changelog

## 1.2.2

- **N parametrizado**: `N_PERTURBATIONS` impar (default 7) no shell, `--n-modes`
  em `random`/`build`/`apply` no nucleo e na CLI; membros = `2N+1`.
- **Perfil de Q duplo**: `artigo_l28` (default, Tabela 1 do artigo) e
  `operacional_l42` (verbatim de `umipert_eof.nml`, zero acima da camada 17,
  exige 42 niveis).
- **`--remove-mean {yes,no}`** no build: default `yes`; o operacional usa `no`.
- **Regioes por indice de grade**: `INDEX_REGIONS` para casar com o operacional
  384x192 (`HN`, `TR`, `HS`, `NAS`, `SAS`); `NAS`/`SAS` diferem em latitude das
  caixas `NSA`/`SSA` do artigo.
- **Nomenclatura operacional de IC**: `--ic-naming operacional` grava
  `GANLCTR<data>S...`, `GANL01P<data>S...`, `GANL01N<data>S...` (mantendo o
  alias `GANLNMC...` que o modelo le).
- **`PBS_STYLE {mpp,nodes,select}`**: `mpp` para `oensMB09`
  (`mppwidth=48,mppnppn=4,mppdepth=6`), `nodes` para `bam_SeasonalClima`
  (`nodes=12:ppn=4`) e `select` para PBS Pro (`select=12:ncpus=6:mpiprocs=4`).
  Defaults do ensaio: fila `pesq`, 48 ranks, 6 threads, 2 h; previsoes longas
  `workq`/23 h via `QUEUE_FCST`/`WALLTIME_FCST`.
- **Default cientifico inalterado**: `q_profile=artigo_l28`, `remove_mean=yes`
  e preset `etsa` continuam sendo o comportamento padrao.
- **Testes**: 352 passed, 1 skipped.
## 1.2.1

- **sdist**: passa a incluir `scripts/mb09.sh` (executavel), `containers/` e `DeltaSigma/`.

## 1.2.0

Ensaio **MB09** de perturbacao de condicoes iniciais (Mendonca & Bonatti, 2009,
*Monthly Weather Review* 137, 1438-1459,
[DOI 10.1175/2008MWR2581.1](https://doi.org/10.1175/2008MWR2581.1)) completo:
nucleo, driver do modelo BAM, metricas e shell de ponta a ponta.

- **`bamci.mb09_tables`**: 6 regioes do artigo (NH, SH, TR, ETR, NSA, SSA),
  presets `oper`/`ext1`/`trop`/`ext2`/`etsa` (default `etsa`), amplitudes
  random/alvo, Tabela 1 de umidade especifica e `region_mask` com bordas
  documentadas.
- **`bamci.mb09`**: `State`, `random_deltas` (passo 1), `build_eofs`
  (diferencas, covariancia, EOF por campo, vento complexo), `apply_modes`
  (+/- com reescala e conversao para o IC), `write_ensemble` e `run_mb09`.
  O build usa explicitamente a janela **6, 9, ..., 36 h** (M=11), com
  `--window`/`--times`/`--analysis-time`; tempos fora da janela sao ignorados
  com aviso. `save_eofs`/`load_eofs` guardam o resultado em `.npz`.
- **CLI MB09**: `python -m bamci.mb09 random|build|apply` (manifestos
  `random.json` e `members.json`; 15 ICs `GANLMB09_000..014`).
- **`bamci.bamrun`**: driver do BAM 1.2 subsazonal: renderiza o `MODELIN`
  preservando o template, gera o PBS com `aprun`/OpenMP, linka
  `GANLNMC...`/`SSTMonthlyDirec...`, prepara um `datain` por membro, submete e
  espera jobs via `qsub`/`qstat`, coleta saidas e roda o `PostGrib`
  (`.unf`+`.ctl`), alem de media/espalhamento. CLI:
  `python -m bamci.bamrun {namelist,pbs,prepare,link,stage,submit,wait,collect,postin,postgrib,ensemble,shift}`.
- **`bamci.score`**: verificacao de ensemble: PAC, RMSE, spread, BSS/REL/RES
  e PECA, lendo `.ctl`/`.unf` do `PostGrib`; CLI
  `python -m bamci.score --ens ... --truth ... --var ...`.
- **`scripts/mb09.sh`**: pipeline 0-9 (IC de controle, SST, 7 ICs aleatorias,
  8 previsoes de 36 h, EOFs, 15 previsoes, pos-processamento, metricas),
  `DRY=1` por default, `--steps` para trechos e lixo em `WORK/trash`
  (sem `rm -rf`).
- **docs**: nova especificacao [`docs/mb09.md`](docs/mb09.md) com o algoritmo,
  as regioes e as decisoes **D1-D11**; secao de ensaio em [`docs/uso.md`](docs/uso.md).
- **Bug corrigido em `bamci.grads.GradsFile`**: variavel 2-D declarada como
  `0 99` no `.ctl` ocupava `nz` niveis e desalinhava todas as variaveis
  seguintes. Agora ocupa 1 nivel; teste sintetico mistura 3-D (42 niveis) com
  `pslc 0 99`/`prec 0 99` intercalados e confere offsets, shapes e leitura.
- **Testes**: **298 passed, 1 skipped** (183 -> 298), incluindo suite nova
  para tabelas, nucleo MB09, driver BAM, score, shell e o caso 2-D do GrADS.

## 1.1.4

O repositorio passa a acompanhar os arquivos `DeltaSigma.Lxxx` (espessuras
sigma) usados pelo pre-processamento, em `DeltaSigma/`.

- **7 arquivos** byte-identicos a
  `/p/projetos/ioper/fix/sazonal/bam_SeasonalClima-1.2_oper/pre/datain` (jaci),
  com os nomes originais: `L009`, `L018`, `L028`, `L042`, `L064`, `L096`,
  `L160` (83 a 1472 bytes). Nao e mais preciso buscar na jaci:
  `--delta-sigma DeltaSigma/DeltaSigma.L042`.
- **`tests/test_deltasigma_files.py`** (22 testes): para cada arquivo, o numero
  de valores e igual ao `NNN` do nome, todos positivos, soma 1 (tolerancia
  `1e-4`) e soma normalizada entregue por `read_delta_sigma` (`1.0`, atol
  `1e-12`).
- **`DeltaSigma.L064`** da origem soma **0.999976** (residuo de 2.4e-5 do
  arquivo operacional, identico ao `L064_NEW-42to64` da jaci); `read_delta_sigma`
  normaliza, entao a CI gerada nao muda. Os outros seis fecham em 1e-6.
- docs: `uso.md` (exemplo com `DeltaSigma/DeltaSigma.L064` para TQ0213L64),
  `formatos.md` (formato `(5F9.6)`, soma 1, procedencia e a nota do L064) e
  `index.md`.

- Testes: 161 -> 183 (mais 1 skip em ambos).

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

# bam-ci

Pre-processamento do **BAM 1.2 subsazonal** (modelo global do CPTEC/INPE) em Python:
transforma analises GRIB em **condicao inicial espectral** e gera a **SST semanal**
no formato que o modelo le.

Substitui a cadeia GrADS + Fortran (`GetICnAGCM`/`runPre`) por um pacote unico,
testavel e sem dependencia de executaveis externos.

```
GRIB (ERA5 / FNL / GFS)
   -> grade gaussiana da truncacao pedida   (T126 -> 384 x 192, T213 -> 640 x 320, ...)
   -> interpolacao vertical sigma L<kmax>   (NewPs + NewSigma)
   -> derivadas espectrais (divergencia / vorticidade)
   -> harmonicos esfericos (FFT + Legendre)
   -> GANL<CASE><AAAAMMDDHH>S.unf.T<Q|L><mend>L<kmax>
```

A resolucao e argumento: qualquer truncacao (126, 213, 254, 666, ...) e qualquer
numero de camadas (L42, L64, ...). TQ126L042 continua sendo o default e reproduz
o arquivo de producao de 11.055.828 B.

O leitor GRIB e escolhido pelo `--case` (ERA5/ERA-Interim -> shortNames ECMWF,
FNL/GFS -> NCEP) e pode ser sobreposto com `--source era5|nl`.

## Instalacao

```bash
python -m venv .venv && . .venv/bin/activate
pip install "bam-ci[grib,netcdf,download]"     # ou: pip install -e ".[dev,grib,netcdf]"
```

Dependencia obrigatoria: `numpy`. Extras: `eccodes` (leitura GRIB),
`xarray`+`netCDF4` (OISST em NetCDF), `cdsapi` (download ERA5 do CDS).

## Uso

```bash
# 1. baixar a analise (ERA5 e o default; grade nativa 0.25, sem regrid)
bamci download --source era5 --date 2026102300 --out DataAnl

# 2. gerar a condicao inicial (TQ126L042 e o default)
bamci run --case ERA5 --date 2026102300 -r TQ0126L042 \
          --grib DataAnl/era5.grib --delta-sigma /caminho/DeltaSigma.L042 \
          --out-dir model/datain

# 2b. qualquer truncacao / numero de camadas
bamci res -r TQ0213L64                    # grade e tamanho antes de rodar
bamci run --case GFS --date 2026110500 --mend 213 --kmax 64 \
          --grib gfs.grib --delta-sigma DeltaSigma.L064 --out-dir model/datain

# 3. SST semanal a partir do RTG do NCEP (ou --source oisst)
bamci sst --source rtg --input "/dados/rtgssthr_grb_0.083.grib2.20261105" \
          --out model/datain/SSTWeekly20261105.G00192

# 4. inspecionar / comparar arquivos
bamci check   GANLERA52026102300S.unf.TQ0126L042 -r TQ0126L042
bamci compare nossa.unf producao.unf -r TQ0213L64

# 5. submeter na jaci (PBS, fila aux)
bamci pbs --case ERA5 --date 2026102300 --walltime 02:00:00 --out icn.pbs --submit
```

## Ensaio MB09 (ensemble)

Implementa o metodo de perturbacao por EOF de Mendonca & Bonatti (2009):
7 condicoes iniciais aleatorias -> 36 h -> EOFs (janela 6-36 h) -> 7 pares
+/- + controle = **15 membros**. O preset default e `etsa`.

```bash
DRY=1 DATE=2026110500 LEAD_DAYS=30 PRESET=etsa scripts/mb09.sh   # imprime
DRY=0 DATE=2026110500 LEAD_DAYS=30 PRESET=etsa scripts/mb09.sh   # executa na jaci
```

O pipeline cobre IC de controle, SST, 7 ICs aleatorias, 8 previsoes de 36 h,
PostGrib, EOFs, 15 previsoes e metricas (`bamci.score`). O nucleo e o driver
tambem podem ser chamados avulso (`python -m bamci.mb09 random|build|apply`,
`python -m bamci.bamrun ...`). Especificacao, regioes e decisoes D1-D13 em
[docs/mb09.md](docs/mb09.md); guia de execucao em [docs/uso.md](docs/uso.md).

Resolucao alvo do ensaio no novo operacional: **TQ0213L064**
(`TRUNC=213 VERT=64`); o perfil de q default e o do artigo
(`artigo_l28`: 32 camadas perturbadas em L64 na chamada padrao, 40 com
`source_delta_sigma=DeltaSigma.L028`) e
`operacional_l42`/`remove_mean=False` sao opt-in de conformidade
historica.

## Formatos

**Condicao inicial** — Fortran UNFORMATTED sequencial, big endian, 3+4K registros:

| registro | conteudo |
|---|---|
| 1 | data + `SigI(K+1)` + `SigL(K)` (sigma) |
| 2 | topografia (coeficientes, m) |
| 3 | ln(ps/10) |
| 4..K+3 | temperatura virtual (K) |
| K+4..3K+3 | divergencia e vorticidade (s-1), alternadas por nivel |
| 3K+4..4K+3 | umidade especifica (kg/kg) |

Coeficientes na ordem do arquivo: **m externo, n interno** (`p = m(M+1) - m(m-1)/2 + (n-m)`),
normalizacao 4pi (`pln(0,0)=1/sqrt(2)`, `Cte=1/sqrt(2)`).

**SST semanal** — UNFORMATTED DIRECT, 2 registros de 192x384 float32 big endian (589.824 B):
mascara oceano/terra (+1 / -1) e SST em K (-999 sobre terra).

## Validacao (contra a producao do CPTEC/INPE)

CI de 2026-11-05 00Z gerada com GFS, comparada com o arquivo operacional:

| grandeza | resultado |
|---|---|
| tamanho / registros | identicos (11.055.828 B, 3+4x42) |
| cabecalho sigma | identico |
| correlacao da topografia | **0,994** |
| pressao de superficie | rms **9,4 hPa** (~1%) |
| temperatura virtual (k=1 / k=21) | rms **1,06 K / 0,76 K** |
| nivel do topo (k=42) | bam-ci reproduz o GFS (254,2 K); o arquivo antigo da 236,2 K |

SST de 2026-11-05 gerada com RTG/NCEP:

| grandeza | resultado |
|---|---|
| mascara oceano/terra | **identica** (24.469 pontos de terra) |
| SST | rms **0,58 K**, correlacao **0,9990**, vies +0,13 K |

## Estrutura

```
src/bamci/
  config.py       parametros (presets E40/E4W/EIT/EIH/FNL/GFS/ERA5)
  levels.py       gaussianas, HybridLevels/DeltaSigma/PressureLevels
  gribio.py       leitura GRIB e montagem da analise (ERA5/NCEP)
  download.py     ERA5 (CDS), FNL e GFS (NOMADS)
  horizontal.py   interpolacao bilinear para a grade do modelo
  vertical.py     NewPs, NewSigma e interpolacao em -ln(p)
  derivatives.py  divergencia/vorticidade por derivadas espectrais
  spectral.py     FFT + Legendre, filtro de Hoskins, ordem do arquivo
  ic_format.py    layout do arquivo de IC, writer/reader, check
  sst.py          SST semanal (OISST v2.1 / RTG) + mascara do modelo
  mb09_tables.py  regioes, presets, amplitudes e Tabela 1 do MB09
  mb09.py         nucleo EOF (random/build/apply) e CLI
  bamrun.py       driver do BAM (MODELIN, PBS, PostGrib) e CLI
  score.py        metricas de ensemble (PAC, RMSE, BSS, PECA) e CLI
  pipeline.py     orquestracao completa
  compare.py      comparacao de duas CIs
  pbs.py          geracao de job PBS
  cli.py          interface de linha de comando
scripts/
  mb09.sh         pipeline 0-9 do ensaio MB09 (DRY=1 por default)
```

## Limitacoes conhecidas

- O filtro de Hoskins (`SmthPerCut=0.12`, `Iter=10`) e aplicado como na producao; no
  fonte legado a rotina existe mas nunca e chamada.
- A grade dos campos espectrais usa **longitudes decrescentes**; a da SST, **crescentes**.
- `bamci download --source fnl` depende de um endpoint NOMADS atualmente indisponivel (403).

## Proveniencia

Implementacao independente das etapas descritas em `doc/Reanalysis2CPTECAGCM.pdf` e no
codigo `GetICnAGCM` do CPTEC/INPE (J. P. Bonatti e colaboradores). O codigo legado
**nao** e redistribuido aqui; os numeros de validacao foram medidos contra arquivos
operacionais do projeto BAM 1.2 subsazonal.

## Licenca

MIT (ver `LICENSE`). Cite com `CITATION.cff`.

## Documentacao

Guia completo em [`docs/`](docs/index.md): instalacao, uso, formatos de arquivo,
validacao contra a producao, referencia da API e problemas conhecidos.
Para publicar como site: `pip install mkdocs && mkdocs gh-deploy`.
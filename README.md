# bam-ci

Pre-processamento do **BAM 1.2 subsazonal** (modelo global do CPTEC/INPE) em Python:
transforma analises GRIB em **condicao inicial espectral** e gera a **SST semanal**
no formato que o modelo le.

Substitui a cadeia GrADS + Fortran (`GetICnAGCM`/`runPre`) por um pacote unico,
testavel e sem dependencia de executaveis externos.

```
GRIB (ERA5 / FNL / GFS)
   -> grade gaussiana TQ126 (384 x 192)
   -> interpolacao vertical sigma L42 (NewPs + NewSigma)
   -> derivadas espectrais (divergencia / vorticidade)
   -> harmonicos esfericos (FFT + Legendre)
   -> GANL<CASE><AAAAMMDDHH>S.unf.TQ0126L042   (11.055.828 B)
```

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

# 2. gerar a condicao inicial
bamci run --case ERA5 --date 2026102300 \
          --levels /caminho/DeltaSigma.L042 \
          --out model/datain

# 3. SST semanal a partir do RTG do NCEP (ou --source oisst)
bamci sst --source rtg --input "/dados/rtgssthr_grb_0.083.grib2.20261105" \
          --out model/datain/SSTWeekly20261105.G00192

# 4. inspecionar / comparar arquivos
bamci check   GANLERA52026102300S.unf.TQ0126L042 --mend 126
bamci compare nossa.unf producao.unf --mend 126

# 5. submeter na jaci (PBS, fila aux)
bamci pbs --case ERA5 --date 2026102300 --walltime 02:00:00 --out icn.pbs --submit
```

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
  pipeline.py     orquestracao completa
  compare.py      comparacao de duas CIs
  pbs.py          geracao de job PBS
  cli.py          interface de linha de comando
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

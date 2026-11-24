# Uso

## 1. Baixar a analise

```bash
bamci download --source era5 --date 2026102300 --out DataAnl   # default: ERA5, grade nativa 0.25, sem regrid
bamci download --source gfs  --date 2026102300 --out DataAnl
```

## 2. Gerar a condicao inicial

```bash
bamci run --case ERA5 --date 2026102300 -r TQ0126L042 \
          --grib DataAnl/era5.grib \
          --delta-sigma /caminho/DeltaSigma.L042 \
          --out-dir model/datain
```

Saida: `GANLERA52026102300S.unf.TQ0126L042` (11.055.828 B).

O leitor GRIB e escolhido pelo `--case`:

| caso | leitor | umidade no GRIB |
|---|---|---|
| ERA5, E40, E4W, EIT, EIH | `era5` (shortNames ECMWF) | `q` (umidade especifica) |
| FNL, GFS | `nl` (shortNames NCEP) | `r` (umidade relativa, convertida) |

Para forcar o leitor (arquivo de origem diferente do caso), use `--source`:

```bash
bamci run --case ERA5 --date 2026102300 --source nl --grib fnl.grb ...
```

Se um nivel pedido nao existir no arquivo, ele e **interpolado em log(p)** e o
programa emite um `RuntimeWarning` listando os niveis; campo ausente por
completo e erro, com a sugestao do `--source` correto.

## 2b. Outras resolucoes

A truncacao espectral e o numero de camadas sao argumentos. Aceita-se um
rotulo unico (`-r`) ou os campos separados (`--mend`, `--kmax`, `--grid`):

```bash
bamci run --case GFS --date 2026110500 -r TQ0213L64 \
          --grib gfs.grib --delta-sigma DeltaSigma/DeltaSigma.L064 --out-dir model/datain

bamci run --case ERA5 --date 2026110500 --mend 666 --kmax 42 \
          --grib era5.grib --delta-sigma DeltaSigma.L042 --out-dir model/datain

bamci run ... -r TQ0254 --kmax 91        # -r traz a truncacao, --kmax sobrescreve
```

Formatos de rotulo aceitos: `TQ0213L64`, `TQ213L64`, `T213L64`, `Q213L64`,
`213L64`, `TQ0213` (sem camadas), `TL95L42` (grade linear).

Antes de rodar, confira grade e tamanho:

```bash
bamci res -r TQ0213L64
bamci res --mend 666 --kmax 42
```

Tabela de referencia (grade gaussiana e tamanho do arquivo):

| resolucao | grade | ncoeff | bytes |
|---|---|---|---|
| TQ0126L042 | 384 x 192 | 16.256 | 11.055.828 |
| TQ0213L064 | 640 x 320 | 46.010 | 47.484.948 |
| TQ0254L042 | 768 x 384 | 65.280 | 44.392.148 |
| TQ0666L042 | 2000 x 1000 | 445.556 | 302.979.828 |

O `DeltaSigma.Lxxx` tem de ter exatamente `kmax` camadas; o pipeline aborta
com mensagem explicita se nao bater.

### Arquivos DeltaSigma que acompanham o pacote

A partir da 1.1.4 o repositorio acompanha os arquivos operacionais em
`DeltaSigma/`, com o nome original de cada um:

| arquivo | camadas | resolucao tipica |
|---|---|---|
| `DeltaSigma/DeltaSigma.L009` | 9 | - |
| `DeltaSigma/DeltaSigma.L018` | 18 | - |
| `DeltaSigma/DeltaSigma.L028` | 28 | TQ062L028 (legado) |
| `DeltaSigma/DeltaSigma.L042` | 42 | TQ0126L042 (default) |
| `DeltaSigma/DeltaSigma.L064` | 64 | TQ0213L64 |
| `DeltaSigma/DeltaSigma.L096` | 96 | - |
| `DeltaSigma/DeltaSigma.L160` | 160 | `--kmax 160` |

Entao nao e mais preciso buscar o arquivo na jaci: basta apontar para a copia
do pacote (ou para a copia local do clone):

```bash
bamci run --case GFS --date 2026110500 -r TQ0213L64 \
          --grib gfs.grib --delta-sigma DeltaSigma/DeltaSigma.L064 --out-dir model/datain

bamci run --case ERA5 --date 2026102300 -r TQ0126L042 \
          --grib DataAnl/era5.grib --delta-sigma DeltaSigma/DeltaSigma.L042 \
          --out-dir model/datain
```

Se o `bamci` estiver instalado em outro diretorio, use o caminho absoluto da
copia no clone. Formato e procedencia estao em
[Formatos de arquivo](formatos.md#deltasigma-lnnn-espessuras-sigma).

## 3. SST semanal

```bash
# RTG do NCEP (dado pontual do dia)
bamci sst --source rtg  --input "/dados/rtgssthr_grb_0.083.grib2.20261105" \
          --out model/datain/SSTWeekly20261105.G00192

# OISST v2.1 (NetCDF diario; media de 7 dias por padrao)
bamci sst --source oisst --input /dados/sst.day.mean.2026.nc --date 20261105 \
          --out model/datain/SSTWeekly20261105.G00192
```

`--out` aceita **arquivo ou diretorio**. Se for diretorio (existente ou com
`/` no final), o nome e montado a partir da data e da grade:

```bash
bamci sst --source oisst --input sst.week.mean.nc --date 2026111200 \
          --out /dados/ioper/tempo/BAM/bam-ci/GANL/
# grava /dados/.../GANL/SSTWeekly20261112.G00192
```

`--date` aceita `AAAAMMDD`, `AAAAMMDDHH` (a hora e descartada, o SST e diario),
`AAAA-MM-DD` e `AAAA-MM-DDTHH`. Se o NetCDF ja for uma **media semanal** (um
unico tempo) e a data pedida estiver fora da janela, o programa usa o tempo
mais proximo do arquivo e emite um `RuntimeWarning`.

Ler NetCDF exige um motor instalado (a mensagem de erro diz qual):

```bash
pip install netCDF4          # ou: pip install "bam-ci[netcdf]"
```

## 4. Inspecionar e comparar

```bash
bamci check   arquivo.unf --mend 126          # registros, sigma, campos
bamci check   arquivo.unf -r TQ0213L64        # idem, resolucao explicita
bamci compare nossa.unf producao.unf -r TQ0213L64   # rms/max na grade e nos coeficientes
```

## 5. Submeter na jaci (PBS, fila aux)

```bash
bamci pbs --case ERA5 --date 2026102300 --ncpus 1 --walltime 02:00:00 \
          --out icn.pbs --submit
```

O script gerado faz `cd "$PBS_O_WORKDIR"` e usa por padrao o interpretador
`.../bam-ci/.venv/bin/python` (os nos de calculo nao tem `python` no PATH).

## 6. Ensaio MB09 (perturbacao de condicoes iniciais)

O ensaio MB09 implementa o metodo de perturbacao por EOF de
**Mendonca & Bonatti (2009)** (*Monthly Weather Review* 137, 1438-1459,
[DOI 10.1175/2008MWR2581.1](https://doi.org/10.1175/2008MWR2581.1)). Para cada
data de analise:

1. gera `N` condicoes iniciais aleatorias em T, U, V, P e Q (`N` **impar**, default 7);
2. integra as `N` + o controle por 36 h, com saida a cada 3 h;
3. descarta as 6 h iniciais e calcula EOFs de 6, 9, ..., 36 h (M=11);
4. aplica `N` pares +/- sobre a analise de controle = **2N+1 membros** (15 com o default).

O preset default e `etsa` (regioes NH, SH, TR, NSA, SSA; campos P, T, Q, U, V).
A regiao e o perfil de umidade seguem o artigo; as escolhas nao especificadas
estao documentadas como **D1-D13** em [MB09](mb09.md).

**Resolucao alvo do novo operacional: TQ0213L064** (substitui T126L28 do artigo e
T126L42 do `oensMB09_bam_xc50`). Exemplo:

```bash
export TRUNC=213 VERT=64                                  # TQ0213L064
export DELTA_SIGMA=DeltaSigma/DeltaSigma.L064
export PRESET=etsa                                       # default
```

O default cientifico do perfil de q continua sendo o do artigo
(`artigo_l28`: em L64 perturba as camadas 1-32 na chamada padrao e 1-40 com
`source_delta_sigma=DeltaSigma.L028`); `operacional_l42` e `remove_mean=False`
sao opt-in de conformidade com o operacional atual.

### Execucao de ponta a ponta

O pipeline esta em [scripts/mb09.sh](../scripts/mb09.sh) e e controlado por
variaveis de ambiente:

```bash
export SUITE_ROOT=/p/projetos/ioper/fix/sazonal/bam_SeasonalClima-1.2_oper
export WORK=$PWD/mb09_work
export DATE=2026110500          # data da analise (YYYYMMDDHH)
export LEAD_DAYS=30             # previsoes finais
export PRESET=etsa              # oper | ext1 | trop | ext2 | etsa
export N_PERTURBATIONS=7        # impar; membros = 2N+1
export IC_NAMING=simple         # simple | operacional

# entradas do pre-processamento
export GRIB=DataAnl/era5.grib
export SST_INPUT=/dados/sst.day.mean.2026.nc
export DELTA_SIGMA=DeltaSigma/DeltaSigma.L042

# 1) modo DRY: so imprime os comandos
DRY=1 scripts/mb09.sh
DRY=1 scripts/mb09.sh --steps 0,1,2      # apenas um trecho

# 2) execucao real (na jaci)
DRY=0 scripts/mb09.sh
```

`N_PERTURBATIONS` tem de ser impar; o shell aborta com mensagem clara se for
par. `--ic-naming operacional` pode ser passado na linha de comando (ou
`IC_NAMING=operacional` no ambiente). `SUITE_ROOT` deve conter `model/datain`,
`model/dataout`, `pos/exec` e o `MODELIN`; os caminhos podem ser sobrepostos
com `MODELIN_TEMPLATE` e `POSTIN_TEMPLATE`. O script usa o `python` de
`PYTHON` (default `python3`) e aceita `--date`, `--work`, `--skip-ic`,
`--steps` e `--ic-naming`.

### Nomenclatura das condicoes iniciais

| `--ic-naming` | controle | membro `k` (k=1..N) + | membro `k` - |
|---|---|---|---|
| `simple` (default) | `GANLNMC<data>S.unf.<res>` | idem | idem |
| `operacional` | `GANLCTR<data>S.unf.<res>` | `GANL{kk}P<data>S.unf.<res>` | `GANL{kk}N<data>S.unf.<res>` |

Em `operacional`, `kk` tem dois digitos (`01`..`N`). O driver mantem tambem um
alias `GANLNMC<data>S.unf.<res>` apontando para o IC do membro, porque e esse o
nome que o modelo le.

### Passos do `scripts/mb09.sh`

| passo | o que faz |
|---|---|
| 0 | checagens de ambiente (binario do modelo, `aprun`, `qsub`, `python`) e criacao de `WORK` |
| 1 | condicao inicial de controle com `bamci run` (ou `--skip-ic` se ja existe) |
| 2 | SST semanal com `bamci sst` + link para `SSTMonthlyDirec<data>.G00<lat>` |
| 3 | `N` ICs aleatorias com `bamci.mb09 random --n-modes N` |
| 4 | submete controle + `N` previsoes de 36 h (`DHFCT=3`) e espera |
| 5 | pos-processa as `N+1` saidas com `PostGrib` (`.unf` + `.ctl`) |
| 6 | constroi os EOFs (`build`) e escreve os `2N+1` ICs (`apply`) |
| 7 | submete as `2N+1` previsoes ate `LEAD_DAYS` (`DHFCT=24`) e espera |
| 8 | pos-processa as `2N+1` e calcula media/espalhamento + metricas (`bamci.score`) |

O passo 9 apenas imprime o resumo e os caminhos. Cada membro usa um
`datain.<membro>` e um `dirfNameOutput` proprio, entao os jobs de 36 h e as
`2N+1` previsoes podem rodar simultaneamente.

### Recursos PBS

Os defaults do ensaio de perturbacao seguem o operacional `oensMB09_bam_xc50`
e podem ser sobrepostos por variavel de ambiente:

| parametro | default | variavel |
|---|---|---|
| fila | `pesq` | `QUEUE` |
| account | `CPTEC` | `ACCOUNT` |
| walltime | `02:00:00` | `WALLTIME` |
| estilo do recurso | `mpp` | `PBS_STYLE` |
| `mppwidth` | `48` | `MPPWIDTH` |
| `mppnppn` | `4` | `MPPNPPN` |
| `mppdepth` / OMP | `6` | `MPPDEPTH` / `OMP_THREADS` |
| ranks MPI (`aprun -n -N -d`) | `48 / 4 / 6` | `RANKS` / `RANKS_PER_NODE` / `OMP_THREADS` |
| previsoes longas | `workq` / `23:00:00` | `QUEUE_FCST` / `WALLTIME_FCST` |

O `PBS_STYLE` escolhe a sintaxe do recurso e **precisa casar com o ambiente e
com o `aprun`**:

| PBS_STYLE | linha gerada | onde usar |
|---|---|---|
| `mpp` | `#PBS -l mppwidth=48,mppnppn=4,mppdepth=6` | suite `oensMB09_bam_xc50` (Moab/Torque) |
| `nodes` | `#PBS -l nodes=12:ppn=4` | suite `bam_SeasonalClima` (Torque classico) |
| `select` | `#PBS -l select=12:ncpus=6:mpiprocs=4` | PBS Pro da jaci (`select`/`ncpus`) |

Em `nodes`/`select`, o numero de chunks vem de `NODES` (default calculado de
`RANKS/RANKS_PER_NODE`). O `aprun` default continua
`-n 48 -N 4 -d 6`; se mudar `RANKS`, `RANKS_PER_NODE` ou `OMP_THREADS`, ajuste
os dois lados para nao pedir recurso e ranks incompativeis. As rodadas curtas
de 36 h usam a fila `pesq` e `02:00:00`; as `2N+1` previsoes longas usam
`QUEUE_FCST`/`WALLTIME_FCST` (default `workq` e `23:00:00`).

### Comandos avulsos

O CLI do MB09 e um modulo (nao ha subcomando `bamci mb09`):

```bash
# N ICs aleatorias (N impar)
python -m bamci.mb09 random \
  --ic model/datain/GANLERA52026110500S.unf.TQ0126L042 \
  --out-dir work/mb09/random --n-modes 7 --seed 0 --preset etsa \
  --resolution TQ0126L042

# EOFs a partir dos .ctl do PostGrib (controle + N membros)
python -m bamci.mb09 build --work work/mb09 \
  --control work/post/ctrl/dataout/TQ0126L042/2026110500/GPOSNMC20261105002026110612C.unf.TQ0126L042.ctl \
  --members work/post/rand01/.../GPOSNMC...ctl ... \
  --preset etsa --ps-units hpa --analysis-time 2026110500 --n-modes 7 \
  --out work/mb09/eof.npz

# 2N+1 ICs finais (controle + N pares +/-)
python -m bamci.mb09 apply --eof work/mb09/eof.npz \
  --ic model/datain/GANLERA52026110500S.unf.TQ0126L042 \
  --out-dir work/mb09/ic --preset etsa --resolution TQ0126L042 \
  --prefix GANLMB09 --n-modes 7

# metricas do ensemble (PAC, RMSE, spread, BSS/REL/RES, PECA)
python -m bamci.score --ens work/post/ens*/.../GPOSNMC...ctl \
  --truth analise.ctl --var tvir --region nh --leads 1,5,10 \
  --out work/score.json
```

O driver do modelo tambem tem CLI propria (`python -m bamci.bamrun`):
`namelist`, `pbs`, `prepare`, `link`, `stage`, `submit`, `wait`, `collect`,
`postin`, `postgrib`, `ensemble` e `shift`.

### Pre-requisitos

- **jaci**: PBS Pro com fila `pesq` (curtas) e `workq` (longas), conta `CPTEC`,
  `aprun`, `qsub`/`qstat` e o ambiente
  `load_modulei_intel_2022_oneapi_mpi.bash` (intel/2022, MPI, NETCDF/PNETCDF/PIO);
- `SUITE_ROOT` da suite operacional do BAM 1.2 subsazonal e permissao de
  escrita em `SUITE_ROOT/model/datain` e `SUITE_ROOT/model/dataout`;
- `MODELIN` e `POSTIN-GRIB` da suite (ou templates com os placeholders
  `{DATAOUT}`, `{POSOUT}`, `{LABELI}`, `{LABELF}`, `{RESOLUTION}`, `{MEMBER}`).

O script **nao usa `rm -rf`**: arquivos descartados vao para `WORK/trash`.

### Aviso de frequencia de saida

As `N+1` rodadas curtas usam **`DHFCT=3`** para produzir a janela de 6-36 h do
MB09 (11 tempos). As `2N+1` previsoes finais usam **`DHFCT=24`**, como o
operacional. O build seleciona explicitamente 6, 9, ..., 36 h; tempos extras
no `.ctl` sao ignorados com aviso.

## Uso programatico

```python
from pathlib import Path
from bamci.config import Config
from bamci.gribio import read_grib, build_analysis
from bamci.download import ERA5_LEVELS
from bamci.pipeline import run

campos = read_grib("era5_pl_2026102300.grib") + read_grib("era5_sl_2026102300.grib")
analise = build_analysis(campos, "era5", ERA5_LEVELS)

cfg = Config(case="ERA5", date="2026102300", data_out=Path("model/datain"))
resultado = run(cfg, analysis=analise, delta_sigma_path="DeltaSigma.L042")
print(resultado.path)
```
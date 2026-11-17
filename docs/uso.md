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

## 2b. Outras resolucoes

A truncacao espectral e o numero de camadas sao argumentos. Aceita-se um
rotulo unico (`-r`) ou os campos separados (`--mend`, `--kmax`, `--grid`):

```bash
bamci run --case GFS --date 2026110500 -r TQ0213L64 \
          --grib gfs.grib --delta-sigma DeltaSigma.L064 --out-dir model/datain

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

## 3. SST semanal

```bash
# RTG do NCEP (dado pontual do dia)
bamci sst --source rtg  --input "/dados/rtgssthr_grb_0.083.grib2.20261105" \
          --out model/datain/SSTWeekly20261105.G00192

# OISST v2.1 (NetCDF diario; media de 7 dias por padrao)
bamci sst --source oisst --input /dados/sst.day.mean.2026.nc --date 20261105 \
          --out model/datain/SSTWeekly20261105.G00192
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

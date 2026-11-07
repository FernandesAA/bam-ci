# Uso

## 1. Baixar a analise

```bash
bamci download --source era5 --date 2026102300 --out DataAnl   # default: ERA5, grade nativa 0.25, sem regrid
bamci download --source gfs  --date 2026102300 --out DataAnl
```

## 2. Gerar a condicao inicial

```bash
bamci run --case ERA5 --date 2026102300 \
          --levels /caminho/DeltaSigma.L042 \
          --out model/datain
```

Saida: `GANLERA52026102300S.unf.TQ0126L042` (11.055.828 B).

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
bamci check   arquivo.unf --mend 126     # registros, sigma, campos, seno de m=0
bamci compare nossa.unf producao.unf --mend 126   # rms/max por campo, na grade e nos coeficientes
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

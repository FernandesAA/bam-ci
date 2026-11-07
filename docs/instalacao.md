# Instalacao

Requisitos: Python >= 3.10 e `numpy`.

```bash
python -m venv .venv
. .venv/bin/activate
pip install "bam-ci[grib,netcdf,download]"
```

## Extras

| extra | para que serve |
|---|---|
| `grib` | leitura de GRIB via `eccodes` (necessario para GFS/FNL/RTG) |
| `netcdf` | leitura de OISST em NetCDF (`xarray` + `netCDF4`) |
| `download` | download de ERA5 pelo CDS (`cdsapi`) |
| `dev` | `pytest` e `scipy` (testes) |

## Do codigo-fonte

```bash
git clone https://github.com/FernandesAA/bam-ci.git
cd bam-ci
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev,grib,netcdf]"
python -m pytest -q
```

## Credencial do CDS (ERA5)

O download de ERA5 usa a API do Copernicus CDS e precisa de uma chave:

```bash
cat > ~/.cdsapirc <<'EOF'
url: https://cds.climate.copernicus.eu/api
key: <sua-chave>
EOF
```

Alternativa: aponte `CDSAPI_RC` para um arquivo existente. Se `~/.cdsapirc` nao
existir, o pacote procura `/oper/scripts/ioper/.cdsapirc` (ambiente da jaci).

## Imagem Apptainer (autocontida)

A imagem traz o interpretador, as dependencias (numpy, eccodes, xarray, netCDF4,
cdsapi) e o `bam-ci` instalado — nao precisa de Python no host.

```bash
apptainer build --fakeroot bam-ci.sif containers/bam-ci.def     # ~5 min, ~700 MB
# ou: ./containers/build.sh

apptainer run  bam-ci.sif --help
apptainer exec --bind "$PWD:/dados" bam-ci.sif bamci check /dados/arquivo.unf --mend 126
```

Para o download de ERA5, monte a credencial:

```bash
apptainer exec --bind "$HOME/.cdsapirc:/root/.cdsapirc" bam-ci.sif \
    bamci download --source era5 --date 2026102300 --out /dados/DataAnl
```

Definicao: [`containers/bam-ci.def`](../containers/bam-ci.def) — base `python:3.12-slim`,
`libeccodes0` do apt e `pip install ".[grib,netcdf,download]"`.

### Pre-requisitos para construir

`apptainer build` precisa de **root ou fakeroot**. Em hosts sem esses privilegios:

```bash
sudo apptainer build bam-ci.sif containers/bam-ci.def     # com root
# ou, num host com Docker:
docker build -f containers/Dockerfile -t bam-ci:1.0.0 .
apptainer build bam-ci.sif docker-daemon://bam-ci:1.0.0
```

Testado nos hosts do projeto: a jaci tem `apptainer 1.5.1`, mas **sem `fakeroot`**
e sem *user namespaces*, entao a construcao precisa ser feita fora dela (root/Docker)
e o `.sif` copiado.

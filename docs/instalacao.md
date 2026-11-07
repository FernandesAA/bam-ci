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

# Problemas conhecidos

## A CI nao bate com a producao

1. **Orientacao da grade** — a CI usa longitudes **decrescentes**; a SST, crescentes.
   Se a topografia nao correlacionar > 0,9, verifique `Config.lon_reverse`.
2. **Filtro de Hoskins** — sem ele o espectro da topografia fica ate 12x alto em
   n=126. Use `Config.smooth_fields = True` (padrao).
3. **Nivel do topo (k=42)** — o arquivo antigo da um valor ~18 K abaixo do GFS.
   O bam-ci reproduz a analise; a diferenca e do legado.

## `eccodes` nao instalado

```
RuntimeError: eccodes nao instalado: use 'module load eccodes' e pip install eccodes
```

Instale o extra `grib` e, em ambientes com modulo, carregue a biblioteca antes:
`module load eccodes`.

## CDS 401 / 404

- `401 Unauthorized`: chave do CDS invalida — gere outra no portal (Perfil -> API key)
  e use a URL `https://cds.climate.copernicus.eu/api` sem prefixo de UID.
- `404`: endpoint antigo (`/api/v2`) foi descontinuado.

## NOMADS 403 (FNL/GFS/RTG)

O NOMADS bloqueia alguns caminhos de `prod` para determinadas redes. O GFS por
`filter_gfs_0p25.pl` funciona; o FNL e o RTG podem exigir espelho local
(no ambiente da jaci: `/p/projetos/ioper/data/external/ncep/AAAA/MM/DD/00/`).

## PBS: `python: command not found`

Os nos de calculo da jaci nao tem `python` no PATH. O `bamci pbs` usa por padrao
`.../bam-ci/.venv/bin/python`; ajuste com `--python`.

## PBS: job roda em `$HOME`

O script gerado faz `cd "$PBS_O_WORKDIR"`. Se você mover o script, ajuste `--workdir`.

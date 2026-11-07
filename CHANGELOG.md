# Changelog

## 1.0.0

Primeira versao publica.

- **Condicao inicial espectral** do BAM (TQ126L042, coordenada sigma):
  leitura GRIB (ERA5/FNL/GFS) -> grade gaussiana -> interpolacao vertical ->
  transformada espectral -> arquivo UNFORMATTED big-endian.
- **SST semanal** no formato do modelo (`SSTWeeklyYYYYMMDD.G00192`), a partir de
  OISST v2.1 ou RTG/NCEP, com a mascara de terra do modelo.
- Validacao contra arquivos de producao do CPTEC/INPE (ver README).
- CLI `bamci` com `check`, `compare`, `download`, `sst`, `pbs`.

# API

## `bamci.config`

- `Config` — parametros da execucao (caso, data, truncacao, niveis, caminhos).
  Presets em `CASE_PRESETS` (E40, E4W, EIT, EIH, FNL, GFS, ERA5).
  `Config.output_path()` devolve o nome do arquivo de IC.

## `bamci.gribio`

- `read_grib(path)` -> lista de `RawField` (shortName, nivel, valores, lat, lon).
- `build_analysis(campos, source, niveis)` -> `Analysis` com T, q, u, v, ps, msl,
  topografia, SST e gelo. `source`: `"era5"` ou `"nl"` (NCEP).
- `specific_humidity_from_rh(rh, T, p_hpa)` — Tetens, como no .gs do FNL.

## `bamci.horizontal`

- `interpolate_bilinear(campo, lat_src, lon_src, lat_dst, lon_dst)`.
- `regular_grid(nlat, nlon)`.

## `bamci.vertical`

- `new_sigma(delta_sigma)` -> interfaces e camadas (media de Exner).
- `new_ps(ps_ml, p_niveis, tv, topo)` -> pressao de superficie (Utils.f90).
- `interp_logp(campo, p_in, p_out)` — Lagrange cubico em -ln(p) com clamp.
- `vert_sigma_inter(campos, p_in, p_out)` — u, v, t, q + correcao de subsuperficie.

## `bamci.derivatives`

- `spectral_div_vort(u, v, nlat, nlon, mend)` — divergencia e vorticidade por
  derivadas espectrais (validado com rotacao solida).

## `bamci.spectral`

- `legendre_norm(mend, mu)`, `gauss_weights(jmax)`.
- `analysis(campo, mend)` / `synthesis(coef, mend, nlat, nlon)`.
- `hoskins_weights(mend, per_cut, iterations)`, `smooth_coefficients(...)`.
- `file_order(coef)` / `from_file_order(vetor, mend)`.

## `bamci.ic_format`

- `ICHeader`, `write_ic(...)`, `read_ic(path)`, `check_ic(path, mend)`,
  `pair_index(mend)`, `ic_file_size(mend, kmax, style)`.

## `bamci.pipeline`

- `run(cfg, analysis=None, grib_path=None, delta_sigma_path=None)` -> `PipelineResult`.

## `bamci.sst`

- `read_oisst(path, end_date, days)`, `read_rtg(paths)`.
- `to_model_grid(campo, lat, lon, reverse_lon=False)`.
- `load_land_mask(path)` — le `ModelLandSeaMask.G00192` (texto, 1=terra).
- `build_weekly(...)` / `build_weekly_rtg(...)`, `write_sstweekly(path, mascara, sst)`.

## `bamci.compare`, `bamci.pbs`

- `compare_ic(a, b, mend)`, `format_report(res)`.
- `PBSSpec`, `render`, `write_script`, `submit`.

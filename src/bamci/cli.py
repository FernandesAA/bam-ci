"""Interface de linha de comando do bamci."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from .ctl import CtlError, read_ctl
from .grads import GradsFile
from .pbs import DEFAULT_PYTHON


def _cmd_ctl(args) -> int:
    print(read_ctl(args.ctl))
    return 0


def _cmd_dump(args) -> int:
    grads = GradsFile(args.ctl)
    field = grads.read(args.var, t=args.t, z=args.z) if args.z else grads.read_volume(args.var, t=args.t)
    finite = np.isfinite(field)
    print("var   : %s" % args.var)
    print("t     : %d" % args.t)
    print("shape : %s" % (field.shape,))
    print("finitos: %d/%d" % (int(finite.sum()), field.size))
    if finite.any():
        print("min/max/mean: %.6g / %.6g / %.6g" % (np.nanmin(field), np.nanmax(field), np.nanmean(field)))
    if args.npz:
        grid = grads.ctl.grid()
        payload = {"field": field, "lat": grid["lat"], "lon": grid["lon"], "level": grid["level"]}
        np.savez(args.npz, **payload)
        print("escrito: %s" % args.npz)
    return 0


def _cmd_download(args) -> int:
    from .download import download_era5, download_ncep, parse_date

    date = parse_date(args.date)
    if args.source == "era5":
        paths = download_era5(args.out, date, flavor=args.cds_flavor)
        for name, path in paths.items():
            print("%s: %s" % (name, path))
    else:
        print(download_ncep(args.out, date, source=args.source, fnl_mode=args.fnl_mode))
    return 0


def _mend_from_args(args, default=None):
    """Resolve --mend / -r para o resto dos subcomandos."""
    from .resolution import resolve

    _, mend, _ = resolve(
        resolution=getattr(args, "resolution", None),
        grid=getattr(args, "grid", None),
        mend=getattr(args, "mend", None),
    )
    return mend if mend is not None else default


def _cmd_check(args) -> int:
    from .ic_format import check_ic, read_ic

    data = read_ic(args.file)
    print("kmax=%d ncoeff=%d data=%s" % (data["kmax"], data["ncoeff"], data["date_inp"]))
    print("sigma: a0=%.6f a-1=%.6f  camadas=%d" % (data["a_hybr"][0], data["a_hybr"][-1], data["b_hybr"].size))
    for name in ("topo", "lnps", "tvir", "divg", "vort", "sphu"):
        arr = np.asarray(data[name])
        print("%-5s rms=%.4e max=%.4e" % (name, float(np.sqrt(np.mean(arr ** 2))), float(np.max(np.abs(arr)))))
    problems = check_ic(args.file, mend=_mend_from_args(args))
    print("problemas:", problems if problems else "nenhum")
    return 0 if not problems else 1


def _cmd_compare(args) -> int:
    from .compare import compare_ic, format_report

    res = compare_ic(args.a, args.b, mend=_mend_from_args(args, default=126),
                     synth=not args.no_grid)
    print(format_report(res))
    return 0


def _cmd_sst(args) -> int:
    import glob as _glob
    from .resolution import resolve
    from .sst import (build_weekly, build_weekly_rtg, default_mask, model_grid,
                      normalize_date, resolve_sst_output)

    _, mend, _ = resolve(resolution=args.resolution, grid=args.grid, mend=args.mend)
    mend = 126 if mend is None else mend
    mask = args.topo or default_mask(mend)
    nlat, _nlon = model_grid(mend)

    date = args.date
    if date is not None:
        date = normalize_date(date)
    out_path = resolve_sst_output(args.out, date, nlat)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.source == "rtg":
        paths = sorted(_glob.glob(args.input))
        if not paths:
            print("nenhum arquivo casou com %s" % args.input, file=sys.stderr)
            return 1
        out = build_weekly_rtg(paths, out_path, topo_path=mask, days=args.days, mend=mend)
    else:
        if date is None:
            print("bamci: erro: --date e obrigatorio para --source oisst (YYYYMMDD)", file=sys.stderr)
            return 1
        out = build_weekly(args.input, date, out_path, days=args.days, mend=mend)
    print("grade  : %dx%d (T%d)" % (nlat, _nlon, mend))
    print("escrito: %s (%d bytes)" % (out, out.stat().st_size))
    return 0


def _cmd_res(args) -> int:
    from .resolution import describe, resolve

    grid, mend, kmax = resolve(
        resolution=args.resolution, grid=args.grid, mend=args.mend, kmax=args.kmax
    )
    if mend is None:
        print("informe -r TQ0213L64 ou --mend 213", file=sys.stderr)
        return 1
    info = describe(grid, mend, kmax)
    for key in ("label", "grid", "mend", "kmax", "nlon", "nlat", "ncoeff", "bytes"):
        if key in info:
            print("%-8s: %s" % (key, info[key]))
    if "kmax" in info:
        print("arquivo  : %s" % (
            "GANL<CASE><AAAAMMDDHH>S.unf." + info["label"]))
    return 0


def _cmd_run(args) -> int:
    from .config import Config
    from .ic_format import ic_file_size
    from .pipeline import run

    cfg = Config(case=args.case, date=args.date)
    cfg.set_resolution(args.resolution, grid=args.grid, mend=args.mend, kmax=args.kmax)
    if args.out_dir:
        cfg.data_out = Path(args.out_dir)
    if args.output_name:
        cfg.output_template = args.output_name
    if args.no_smooth:
        cfg.smooth_fields = False
    if args.keep_lon:
        cfg.lon_reverse = False

    problems = cfg.check()
    if problems:
        print("bamci: configuracao invalida: " + "; ".join(problems), file=sys.stderr)
        return 1

    p_levels = None
    if args.p_levels:
        p_levels = [float(v) for v in args.p_levels.replace(",", " ").split()]

    result = run(cfg, grib_path=args.grib, delta_sigma_path=args.delta_sigma,
                 p_levels=p_levels, source=args.source)
    expected = ic_file_size(cfg.mend_out, cfg.kmax_out)
    size = result.path.stat().st_size
    print("resolucao: %s" % cfg.resolution_label)
    print("leitor   : %s" % (args.source or cfg.case))
    print("escrito  : %s" % result.path)
    print("tamanho  : %d bytes (esperado %d) %s" % (size, expected, "OK" if size == expected else "DIVERGE"))
    return 0 if size == expected else 1


def _cmd_pbs(args) -> int:
    from .pbs import PBSSpec, render, submit, write_script

    command = [args.python, "-m", "bamci", "run", "--case", args.case, "--date", args.date]
    if args.resolution:
        command += ["-r", args.resolution]
    if args.mend:
        command += ["--mend", str(args.mend)]
    if args.kmax:
        command += ["--kmax", str(args.kmax)]
    if args.grid:
        command += ["--grid", args.grid]
    if args.grib:
        command += ["--grib", args.grib]
    if args.delta_sigma:
        command += ["--delta-sigma", args.delta_sigma]

    spec = PBSSpec(
        job_name=args.name,
        command=command,
        queue=args.queue,
        ncpus=args.ncpus,
        walltime=args.walltime,
        workdir=args.workdir,
        log_name=args.log,
    )
    if args.out:
        write_script(spec, args.out)
        print("escrito: %s" % args.out)
    else:
        print(render(spec))
    if args.submit:
        print(submit(args.out or "job.pbs", dry_run=args.dry_run))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bamci", description="Pre-processamento BAM 1.2 subsazonal")
    sub = parser.add_subparsers(dest="command", required=True)

    p_ctl = sub.add_parser("ctl", help="mostra o resumo de um descritor .ctl")
    p_ctl.add_argument("ctl")
    p_ctl.set_defaults(func=_cmd_ctl)

    p_dump = sub.add_parser("dump", help="le um campo do binario GrADS")
    p_dump.add_argument("ctl")
    p_dump.add_argument("var")
    p_dump.add_argument("-t", type=int, default=1, help="tempo (1-based)")
    p_dump.add_argument("-z", type=int, default=None, help="nivel (1-based); omitido = todos")
    p_dump.add_argument("--npz", default=None, help="salva o campo em .npz")
    p_dump.set_defaults(func=_cmd_dump)

    p_res = sub.add_parser("res", help="mostra grade, coeficientes e tamanho de arquivo de uma resolucao")
    p_res.add_argument("-r", "--resolution", default=None, help="ex.: TQ0213L64, T213L42")
    p_res.add_argument("--grid", default=None, choices=("Q", "L"), help="grade (padrao Q)")
    p_res.add_argument("--mend", type=int, default=None, help="truncacao espectral")
    p_res.add_argument("--kmax", type=int, default=None, help="numero de camadas")
    p_res.set_defaults(func=_cmd_res)

    p_run = sub.add_parser("run", help="gera a condicao inicial (pipeline completo)")
    p_run.add_argument("--case", default="ERA5")
    p_run.add_argument("--date", required=True, help="YYYYMMDDHH")
    p_run.add_argument("-r", "--resolution", default=None, help="ex.: TQ0213L64, T213L42")
    p_run.add_argument("--grid", default=None, choices=("Q", "L"))
    p_run.add_argument("--mend", type=int, default=None, help="truncacao espectral")
    p_run.add_argument("--kmax", type=int, default=None, help="numero de camadas")
    p_run.add_argument("--grib", required=True, help="arquivo GRIB da analise")
    p_run.add_argument("--source", default=None, choices=("era5", "nl"),
                       help="leitor GRIB (padrao: derivado do --case)")
    p_run.add_argument("--delta-sigma", required=True, help="DeltaSigma.Lxxx")
    p_run.add_argument("--p-levels", default=None, help="niveis de pressao (hPa), ex.: 1000,975,...")
    p_run.add_argument("--out-dir", default=None)
    p_run.add_argument("--output-name", default=None, help="template do nome de saida")
    p_run.add_argument("--no-smooth", action="store_true", help="desliga o filtro de Hoskins")
    p_run.add_argument("--keep-lon", action="store_true", help="nao inverte as longitudes")
    p_run.set_defaults(func=_cmd_run)

    p_chk = sub.add_parser("check", help="valida um arquivo de IC")
    p_chk.add_argument("file")
    p_chk.add_argument("-r", "--resolution", default=None, help="ex.: TQ0213L64 (sobrescreve --mend)")
    p_chk.add_argument("--grid", default=None, choices=("Q", "L"))
    p_chk.add_argument("--mend", type=int, default=None)
    p_chk.set_defaults(func=_cmd_check)

    p_cmp = sub.add_parser("compare", help="compara dois arquivos de IC")
    p_cmp.add_argument("a")
    p_cmp.add_argument("b")
    p_cmp.add_argument("-r", "--resolution", default=None, help="ex.: TQ0213L64 (sobrescreve --mend)")
    p_cmp.add_argument("--grid", default=None, choices=("Q", "L"))
    p_cmp.add_argument("--mend", type=int, default=126)
    p_cmp.add_argument("--no-grid", action="store_true")
    p_cmp.set_defaults(func=_cmd_compare)

    p_dl = sub.add_parser("download", help="baixa ERA5 (default), FNL ou GFS")
    p_dl.add_argument("--source", default="era5", choices=("era5", "fnl", "gfs"))
    p_dl.add_argument("--date", required=True, help="YYYYMMDDHH")
    p_dl.add_argument("--out", default="DataAnl")
    p_dl.add_argument("--cds-flavor", default="beta", choices=("beta", "legacy"))
    p_dl.add_argument("--fnl-mode", default="filter", choices=("filter", "direct"))
    p_dl.set_defaults(func=_cmd_download)

    p_sst = sub.add_parser("sst", help="gera o SSTWeekly.G00192 (OISST ou RTG)")
    p_sst.add_argument("--source", default="oisst", choices=("oisst", "rtg"))
    p_sst.add_argument("--input", required=True, help="arquivo NetCDF (oisst) ou padrao glob (rtg)")
    p_sst.add_argument("--date", default=None, help="YYYYMMDD (oisst)")
    p_sst.add_argument("--days", type=int, default=1,
                       help="dias a combinar (RTG e pontual: padrao 1)")
    p_sst.add_argument("-r", "--resolution", default=None, help="ex.: TQ0213L64 (define a grade)")
    p_sst.add_argument("--grid", default=None, choices=("Q", "L"))
    p_sst.add_argument("--mend", type=int, default=None, help="truncacao (padrao 126)")
    p_sst.add_argument("--topo", default=None,
                       help="mascara de terra (padrao: ModelLandSeaMask.G00<lats> da truncacao)")
    p_sst.add_argument("--out", required=True)
    p_sst.set_defaults(func=_cmd_sst)

    p_pbs = sub.add_parser("pbs", help="gera (e opcionalmente submete) um job PBS")
    p_pbs.add_argument("--case", default="ERA5")
    p_pbs.add_argument("--date", required=True)
    p_pbs.add_argument("-r", "--resolution", default=None, help="ex.: TQ0213L64")
    p_pbs.add_argument("--grid", default=None, choices=("Q", "L"))
    p_pbs.add_argument("--mend", type=int, default=None)
    p_pbs.add_argument("--kmax", type=int, default=None)
    p_pbs.add_argument("--grib", default=None, help="repassa --grib ao run")
    p_pbs.add_argument("--delta-sigma", default=None, help="repassa --delta-sigma ao run")
    p_pbs.add_argument("--name", default="bamci")
    p_pbs.add_argument("--queue", default="aux")
    p_pbs.add_argument("--ncpus", type=int, default=1)
    p_pbs.add_argument("--walltime", default="01:00:00")
    p_pbs.add_argument("--workdir", default=".")
    p_pbs.add_argument("--python", default=DEFAULT_PYTHON)
    p_pbs.add_argument("--log", default=None)
    p_pbs.add_argument("--out", default=None)
    p_pbs.add_argument("--submit", action="store_true")
    p_pbs.add_argument("--dry-run", action="store_true")
    p_pbs.set_defaults(func=_cmd_pbs)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (CtlError, FileNotFoundError, ValueError, IndexError, KeyError) as exc:
        print("bamci: erro: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

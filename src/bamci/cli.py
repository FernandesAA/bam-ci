"""Interface de linha de comando do bamci."""

from __future__ import annotations

import argparse
import sys

import numpy as np

from .ctl import CtlError, read_ctl
from .grads import GradsFile
from .pbs import DEFAULT_PYTHON
from .sst import DEFAULT_MASK


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


def _cmd_check(args) -> int:
    from .ic_format import check_ic, read_ic

    data = read_ic(args.file)
    print("kmax=%d ncoeff=%d data=%s" % (data["kmax"], data["ncoeff"], data["date_inp"]))
    print("sigma: a0=%.6f a-1=%.6f  camadas=%d" % (data["a_hybr"][0], data["a_hybr"][-1], data["b_hybr"].size))
    for name in ("topo", "lnps", "tvir", "divg", "vort", "sphu"):
        arr = np.asarray(data[name])
        print("%-5s rms=%.4e max=%.4e" % (name, float(np.sqrt(np.mean(arr ** 2))), float(np.max(np.abs(arr)))))
    problems = check_ic(args.file, mend=args.mend)
    print("problemas:", problems if problems else "nenhum")
    return 0 if not problems else 1


def _cmd_compare(args) -> int:
    from .compare import compare_ic, format_report

    res = compare_ic(args.a, args.b, mend=args.mend, synth=not args.no_grid)
    print(format_report(res))
    return 0


def _cmd_sst(args) -> int:
    import glob as _glob
    from .sst import build_weekly, build_weekly_rtg

    if args.source == "rtg":
        paths = sorted(_glob.glob(args.input))
        if not paths:
            print("nenhum arquivo casou com %s" % args.input, file=sys.stderr)
            return 1
        out = build_weekly_rtg(paths, args.out, topo_path=args.topo, days=args.days)
    else:
        out = build_weekly(args.input, args.date, args.out, days=args.days)
    print("escrito: %s (%d bytes)" % (out, out.stat().st_size))
    return 0


def _cmd_pbs(args) -> int:
    from .pbs import PBSSpec, render, submit, write_script

    spec = PBSSpec(
        job_name=args.name,
        command=[args.python, "-m", "bamci", "run", "--case", args.case, "--date", args.date],
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

    p_chk = sub.add_parser("check", help="valida um arquivo de IC")
    p_chk.add_argument("file")
    p_chk.add_argument("--mend", type=int, default=None)
    p_chk.set_defaults(func=_cmd_check)

    p_cmp = sub.add_parser("compare", help="compara dois arquivos de IC")
    p_cmp.add_argument("a")
    p_cmp.add_argument("b")
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
    p_sst.add_argument("--topo", default=DEFAULT_MASK,
                       help="mascara de terra do modelo (padrao: ModelLandSeaMask.G00192)")
    p_sst.add_argument("--out", required=True)
    p_sst.set_defaults(func=_cmd_sst)

    p_pbs = sub.add_parser("pbs", help="gera (e opcionalmente submete) um job PBS")
    p_pbs.add_argument("--case", default="ERA5")
    p_pbs.add_argument("--date", required=True)
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

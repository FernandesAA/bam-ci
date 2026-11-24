"""Testes do driver do modelo BAM e do ensaio MB09 (sem JACI, sem subprocess real)."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from bamci import bamrun


TEMPLATE = """&MODEL_RES
 trunc = 126, vert = 42, dt = 450,          ! resolucao do ensaio
 IDATEI = 00,01,01,2020,
 IDATEW = 00,01,01,2020,
 IDATEF = 00,01,01,2020,
 NMSST = 'sstmtd', DHFCT = 24,
 DHRES = 24,
 PREFX = 'NMC', PREFY = 'NMC',              ! o modelo le GANLNMC...
 TABLE = 'c',
 path_in = '/old/datain',
 dirfNameOutput = '/old/dataout',
 initlz = .TRUE., ifsst = .TRUE., fint = .TRUE.,
 UNIFIED = .TRUE., ICCON = 'ARA',
 ISCON = 'TIED', ILCON = 'HUMO', IGWD = 'YES',
 CRDCLD = 6, schemes = 3, OCFLUX = 'UKME',
 SLABOCEAN = 'SLAB', ICEMODEL = 'SSIB', co2val = 370.0
/
"""


def _template(tmp_path: Path) -> Path:
    path = tmp_path / "MODELIN.tmpl"
    path.write_text(TEMPLATE, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# namelist
# ---------------------------------------------------------------------------
def test_render_modelin_altera_so_as_chaves_pedidas(tmp_path):
    out = bamrun.render_modelin(
        _template(tmp_path),
        tmp_path / "MODELIN",
        date="2026110500",
        dhfct=3,
        trunc=126,
        vert=42,
        dt=450,
        path_in="/novo/datain",
        dirfNameOutput="/novo/dataout",
    )
    text = out.read_text(encoding="utf-8")

    # valores novos
    assert "IDATEI = 00,05,11,2026," in text
    assert "IDATEF = 00,05,11,2026," in text
    assert "DHFCT = 3," in text
    assert "path_in = '/novo/datain'," in text
    assert "dirfNameOutput = '/novo/dataout'," in text

    # nada de sobra do template
    assert "/old/datain" not in text
    assert "/old/dataout" not in text
    assert "DHFCT = 24" not in text

    # comentarios, ordem e demais chaves intactos
    assert "! resolucao do ensaio" in text
    assert "! o modelo le GANLNMC..." in text
    assert "co2val = 370.0" in text
    assert "OCFLUX = 'UKME'" in text
    assert text.index("trunc") < text.index("IDATEI") < text.index("path_in")
    assert text.count("&MODEL_RES") == 1
    assert text.rstrip().endswith("/")


def test_render_modelin_text_aceita_date_e_idatef_explicito():
    text = bamrun.render_modelin_text(
        TEMPLATE, date="2026110500", idatef="2026110612", trunc=126
    )
    assert "IDATEI = 00,05,11,2026," in text
    assert "IDATEF = 12,06,11,2026," in text


def test_render_modelin_chave_ausente_falha(tmp_path):
    with pytest.raises(KeyError):
        bamrun.render_modelin(_template(tmp_path), tmp_path / "X", nao_existe=1)


def test_data_helpers():
    assert bamrun.parse_date("2026-11-05T00") == "2026110500"
    assert bamrun.format_idate("2026110500") == "00,05,11,2026"
    assert bamrun.shift_date("2026110500", 36) == "2026110612"
    assert bamrun.run_id("2026110500", "2026110612") == "2026110500-2026110612"


# ---------------------------------------------------------------------------
# PBS
# ---------------------------------------------------------------------------
def test_render_pbs_contem_aprun_nodes_walltime():
    text = bamrun.render_pbs(
        "bam-ctrl",
        "/work",
        "/work/exec_NMC.ctrl",
        "MODELIN",
        queue="workq",
        account="CPTEC",
        walltime="23:00:00",
        nodes=3,
        ppn=40,
        omp_threads=10,
        env_source="/p/projetos/load_modules.bash",
    )
    assert "#PBS -N bam-ctrl" in text
    assert "#PBS -q workq" in text
    assert "#PBS -A CPTEC" in text
    assert "#PBS -l nodes=3:ppn=40" in text
    assert "#PBS -l walltime=23:00:00" in text
    assert "cd /work/exec_NMC.ctrl" in text
    assert "export OMP_NUM_THREADS=10" in text
    assert "export OMP_STACKSIZE=256M" in text
    assert "ulimit -s unlimited" in text
    assert "source /p/projetos/load_modules.bash" in text
    assert (
        "aprun -m1000h -n 12 -N 4 -d 10 ./ParModel_MPI < MODELIN > print.out"
        in text
    )


def test_render_pbs_aprun_args_customizado_e_validacao():
    text = bamrun.render_pbs(
        "x", "/w", "/w/e", "MODELIN", aprun_args=["-n", "4", "-d", "8"]
    )
    assert "aprun -n 4 -d 8 ./ParModel_MPI" in text

    with pytest.raises(ValueError):
        bamrun.render_pbs("x", "/w", "/w/e", nodes=0)
    with pytest.raises(ValueError):
        bamrun.render_pbs("x", "/w", "/w/e", walltime="48:00:00")


def test_write_pbs(tmp_path):
    path = bamrun.write_pbs(tmp_path / "job.pbs", "#!/bin/bash\n")
    assert path.read_text() == "#!/bin/bash\n"
    assert path.stat().st_mode & 0o111


# ---------------------------------------------------------------------------
# diretorio de execucao e entradas
# ---------------------------------------------------------------------------
def test_prepare_run_dir_copia_executavel(tmp_path):
    exe = tmp_path / "ParModel_MPI"
    exe.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    exe.chmod(0o755)

    run = bamrun.prepare_run_dir(tmp_path / "work", "m1", exe)
    assert run.run_dir.name == "exec_NMC.m1"
    assert run.executable.name == "ParModel_MPI"
    assert run.executable.is_file()
    assert run.executable.stat().st_mode & 0o111
    assert run.modelin.name == "MODELIN"


def test_prepare_run_dir_dry_run_nao_cria(tmp_path):
    run = bamrun.prepare_run_dir(
        tmp_path / "work", "m1", tmp_path / "ParModel_MPI", dry_run=True
    )
    assert not run.run_dir.exists()


def test_link_inputs_cria_nomes_ganlnmc_e_sstmonthly(tmp_path):
    suite = tmp_path / "suite"
    datain = suite / "model" / "datain"
    datain.mkdir(parents=True)
    ic = tmp_path / "GANLERA52026110500S.unf.TQ0126L042"
    sst = tmp_path / "SSTWeekly20261105.G00192"
    ic.write_bytes(b"conteudo-ic")
    sst.write_bytes(b"conteudo-sst")

    linked = bamrun.link_inputs(suite, ic, sst, "2026110500", mode="copy")
    assert linked.ic.name == "GANLNMC2026110500S.unf.TQ0126L042"
    assert linked.sst.name == "SSTMonthlyDirec20261105.G00192"
    assert linked.ic.read_bytes() == b"conteudo-ic"
    assert linked.sst.read_bytes() == b"conteudo-sst"
    assert linked.ic.parent == datain


def test_link_inputs_dry_run_nao_cria(tmp_path):
    ic = tmp_path / "GANLERA52026110500S.unf.TQ0126L042"
    sst = tmp_path / "SSTWeekly20261105.G00192"
    ic.write_bytes(b"ic")
    sst.write_bytes(b"sst")
    linked = bamrun.link_inputs(
        tmp_path / "suite", ic, sst, "2026110500", dry_run=True
    )
    assert not linked.ic.exists()
    assert not linked.sst.exists()


def test_stage_datain_espelha_estaticos(tmp_path):
    suite = tmp_path / "suite"
    datain = suite / "model" / "datain"
    datain.mkdir(parents=True)
    (datain / "ModelLandSeaMask.G00192").write_text("mascara")
    (datain / "SSTWeekly20261105.G00192").write_bytes(b"velho")
    ic = tmp_path / "GANLERA52026110500S.unf.TQ0126L042"
    sst = tmp_path / "SSTWeekly20261105.G00192"
    ic.write_bytes(b"ic")
    sst.write_bytes(b"sst")

    staged = bamrun.stage_datain(
        tmp_path / "staging", "m1", suite, ic, sst,
        labeli="2026110500", mode="copy",
    )
    assert staged.name == "datain.m1"
    assert (staged / "ModelLandSeaMask.G00192").read_text() == "mascara"
    assert (staged / "GANLNMC2026110500S.unf.TQ0126L042").read_bytes() == b"ic"
    assert (staged / "SSTMonthlyDirec20261105.G00192").read_bytes() == b"sst"
    # o SSTWeekly/IC antigos do datain da suite nao sao copiados
    assert not (staged / "SSTWeekly20261105.G00192").exists()


# ---------------------------------------------------------------------------
# PBS: qstat
# ---------------------------------------------------------------------------
def test_parse_qstat_state_formatos():
    one_line = "1234.server user workq job 1 1 1 -- 23:00 R 00:01\n"
    assert bamrun.parse_qstat_state(one_line, "1234.server") == "R"
    full = "Job Id: 1234.server\n    job_state = C\n    exit_status = 0\n"
    assert bamrun.parse_qstat_state(full) == "C"
    assert bamrun.parse_qstat_state("qstat: Unknown Job Id 1234.server") is None
    assert bamrun.parse_qstat_state("") is None


def test_wait_job_fake_qstat_r_depois_sumido(tmp_path):
    counter = tmp_path / "counter"
    script = tmp_path / "qstat"
    script.write_text(
        "#!/bin/bash\n"
        "n=$(cat %s 2>/dev/null || echo 0)\n"
        "n=$((n + 1))\n"
        "echo $n > %s\n"
        "if [ \"$n\" -eq 1 ]; then\n"
        "  echo '1234.server user workq job 1 1 1 -- 23:00 R 00:01'\n"
        "  exit 0\n"
        "fi\n"
        "echo 'qstat: Unknown Job Id 1234.server'\n"
        "exit 153\n" % (counter, counter),
        encoding="utf-8",
    )
    script.chmod(0o755)

    result = bamrun.wait_job(
        "1234.server", poll=0.01, timeout=5, qstat=str(script)
    )
    assert result.disappeared is True
    assert result.state == "C"
    assert result.polls >= 2
    assert int(counter.read_text()) >= 2


# ---------------------------------------------------------------------------
# saidas e pos-processamento
# ---------------------------------------------------------------------------
def test_collect_outputs(tmp_path):
    run_dir = tmp_path / "dataout" / "TQ0126L042" / "2026110500-2026110612"
    run_dir.mkdir(parents=True)
    (run_dir / "GFCTNMC20261105002026110612F.unf.TQ0126L042.outmdt").write_text("x")
    (run_dir / "GFCTNMC20261105002026110612F.unf.TQ0126L042.outatt").write_text("y")
    (run_dir / "leia-me.txt").write_text("z")

    found = bamrun.collect_outputs(
        tmp_path / "dataout", "2026110500", "2026110612"
    )
    assert [p.name for p in found] == [
        "GFCTNMC20261105002026110612F.unf.TQ0126L042.outatt",
        "GFCTNMC20261105002026110612F.unf.TQ0126L042.outmdt",
    ]
    with pytest.raises(FileNotFoundError):
        bamrun.collect_outputs(tmp_path / "dataout", "2026110500", "2026110700")


def test_run_postgrib_dry_run_caminhos(tmp_path):
    result = bamrun.run_postgrib(
        tmp_path / "pos",
        "2026110500",
        "2026110612",
        tmp_path / "POSTIN",
        dry_run=True,
    )
    assert result.unf.name == "GPOSNMC20261105002026110612C.unf.TQ0126L042"
    assert result.ctl.name == result.unf.name + ".ctl"
    assert result.unf.parent.name == "2026110500"


def test_render_postin(tmp_path):
    template = tmp_path / "POSTIN"
    template.write_text("arq {labeli} {labelf} @RESOLUTION@ @NAOEXISTE@")
    out = bamrun.render_postin(
        template, tmp_path / "POSTIN.1",
        labeli="2026110500", labelf="2026110612", RESOLUTION="TQ0126L042",
    )
    assert out.read_text() == "arq 2026110500 2026110612 TQ0126L042 @NAOEXISTE@"


def _write_fake_postgrib_ctl(root: Path, member: int, nx=4, ny=2, nt=2) -> Path:
    """Escreve um .ctl/.unf minimo no layout do PostGrib (ps 2-D + 4 campos 3-D)."""
    import numpy as np

    names = ["pslc", "tvir", "umes", "uvel", "vvel"]
    folder = root / ("m%d" % member)
    folder.mkdir(parents=True, exist_ok=True)
    values = np.arange(nt * len(names) * ny * nx, dtype=">f4") + member
    (folder / "out.unf").write_bytes(values.tobytes())
    (folder / "out.ctl").write_text(
        "dset ^out.unf\n"
        "undef -9.99e8\n"
        "xdef %d linear 0 90\n"
        "ydef %d linear -45 90\n"
        "zdef 1 linear 1 1\n"
        "tdef %d linear 00Z05nov2026 3hr\n"
        "vars 5\n"
        "pslc 0 99 pressao\n"
        "tvir 1 99 temperatura\n"
        "umes 1 99 umidade\n"
        "uvel 1 99 u\n"
        "vvel 1 99 v\n"
        "endvars\n" % (nx, ny, nt),
        encoding="utf-8",
    )
    return folder / "out.ctl"


def test_ensemble_mean_spread_le_ctl_do_postgrib(tmp_path):
    import numpy as np

    ctls = [_write_fake_postgrib_ctl(tmp_path, m) for m in range(3)]
    written = bamrun.ensemble_mean_spread(ctls, tmp_path / "mean", prefix="ens")
    assert len(written) == 5

    ps = np.load(tmp_path / "mean" / "ens_pslc.npz")
    assert ps["mean"].shape == (2, 1, 2, 4)    # ps 2-D lido como 1 nivel
    assert np.all(ps["spread"] > 0)

    tvir = np.load(tmp_path / "mean" / "ens_tvir.npz")
    assert tvir["mean"].shape == (2, 1, 2, 4)  # 3-D: (nt, nz, ny, nx)


# ---------------------------------------------------------------------------
# bash de ponta a ponta
# ---------------------------------------------------------------------------
def test_scripts_mb09_sh_sintaxe():
    script = Path(__file__).resolve().parents[1] / "scripts" / "mb09.sh"
    assert script.is_file()
    subprocess.run(["bash", "-n", str(script)], check=True)

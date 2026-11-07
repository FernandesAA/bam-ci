import subprocess

import pytest

from bamci.pbs import PBSSpec, render, submit, write_script


def _spec(**kw):
    base = dict(job_name="icn-EIT-2026100100", command=["python", "-m", "bamci", "run",
                                                        "--case", "EIT", "--date", "2026100100"],
                workdir=".", walltime="02:00:00")
    base.update(kw)
    return PBSSpec(**base)


def test_render_has_directives():
    txt = render(_spec(queue="aux", ncpus=4))
    assert "#PBS -q aux" in txt
    assert "#PBS -l select=1:ncpus=4" in txt
    assert "#PBS -l walltime=02:00:00" in txt
    assert "#PBS -j oe" in txt
    assert "python -m bamci run --case EIT --date 2026100100" in txt


def test_default_python_is_project_venv():
    from bamci.pbs import DEFAULT_PYTHON
    assert DEFAULT_PYTHON.endswith("/bam-ci/.venv/bin/python")


def test_default_workdir_uses_pbs_o_workdir():
    txt = render(_spec())
    assert 'cd "$PBS_O_WORKDIR"' in txt


def test_render_quotes_paths_with_space():
    txt = render(_spec(workdir="/p/projetos/ioper/my dir"))
    assert "'/p/projetos/ioper/my dir'" in txt


def test_walltime_limit():
    with pytest.raises(ValueError):
        render(_spec(walltime="48:00:00"))


def test_same_day_24h_is_allowed():
    txt = render(_spec(walltime="24:00:00"))
    assert "walltime=24:00:00" in txt


def test_ncpus_must_be_positive():
    with pytest.raises(ValueError):
        render(_spec(ncpus=0))


def test_write_script_is_executable(tmp_path):
    p = write_script(_spec(), tmp_path / "job.pbs")
    assert p.stat().st_mode & 0o111
    assert p.read_text().startswith("#!/bin/bash")


def test_submit_dry_run_does_not_call_qsub(monkeypatch):
    def boom(*a, **kw):
        raise AssertionError("qsub nao deveria ser chamado")
    monkeypatch.setattr(subprocess, "run", boom)
    assert submit("job.pbs", dry_run=True) == "qsub job.pbs"

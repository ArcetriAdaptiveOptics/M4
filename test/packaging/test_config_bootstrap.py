"""`import m4` ensures the `SYSTEM: simulated.devices:` flags exist in the opticalib config."""
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
from ruamel.yaml import YAML

PRISTINE = os.environ.get("M4_TEST_PRISTINE_CONF")
pytestmark = pytest.mark.skipif(
    not PRISTINE or not Path(PRISTINE).exists(),
    reason="M4_TEST_PRISTINE_CONF not set",
)

FLAGS = [
    "dm", "interferometer", "accelerometers", "angleRotator", "m4Exapode",
    "dp", "parSlider", "par", "rmSlider", "rm", "tempSensors",
]


def _make(tmp_path, mutate):
    y = YAML()
    with open(PRISTINE) as f:
        cfg = y.load(f)
    cfg["SYSTEM"]["data_path"] = str(tmp_path / "data")
    mutate(cfg)
    path = tmp_path / "configuration.yaml"
    with open(path, "w") as f:
        y.dump(cfg, f)
    return path


def _import_m4(path):
    env = os.environ.copy()
    env["AOCONF"] = str(path)
    return subprocess.run(
        [sys.executable, "-c", "import m4"],
        capture_output=True, text=True, timeout=120, env=env,
    )


def _load(path):
    with open(path) as f:
        return YAML().load(f)


def test_block_missing_is_created(tmp_path):
    path = _make(tmp_path, lambda c: c["SYSTEM"].pop("simulated.devices"))
    res = _import_m4(path)
    assert res.returncode == 0, res.stderr
    sim = _load(path)["SYSTEM"]["simulated.devices"]
    assert all(sim[k] is True for k in FLAGS)
    assert set(sim) == set(FLAGS)


def test_existing_false_is_preserved_and_file_untouched(tmp_path):
    def mut(c):
        c["SYSTEM"]["simulated.devices"]["dm"] = False

    path = _make(tmp_path, mut)
    before = path.read_bytes()
    res = _import_m4(path)
    assert res.returncode == 0, res.stderr
    assert path.read_bytes() == before
    assert _load(path)["SYSTEM"]["simulated.devices"]["dm"] is False


def test_partial_flags_only_missing_added(tmp_path):
    def mut(c):
        sim = c["SYSTEM"]["simulated.devices"]
        for k in list(sim):
            if k not in ("dm", "par"):
                del sim[k]
        sim["dm"] = False
        sim["par"] = False

    path = _make(tmp_path, mut)
    res = _import_m4(path)
    assert res.returncode == 0, res.stderr
    sim = _load(path)["SYSTEM"]["simulated.devices"]
    assert sim["dm"] is False and sim["par"] is False
    assert all(sim[k] is True for k in FLAGS if k not in ("dm", "par"))
    assert set(sim) == set(FLAGS)


def test_empty_block(tmp_path):
    def mut(c):
        c["SYSTEM"]["simulated.devices"] = {}

    path = _make(tmp_path, mut)
    res = _import_m4(path)
    assert res.returncode == 0, res.stderr
    sim = _load(path)["SYSTEM"]["simulated.devices"]
    assert set(sim) == set(FLAGS) and all(v is True for v in sim.values())


def test_null_block(tmp_path):
    def mut(c):
        c["SYSTEM"]["simulated.devices"] = None

    path = _make(tmp_path, mut)
    res = _import_m4(path)
    assert res.returncode == 0, res.stderr
    sim = _load(path)["SYSTEM"]["simulated.devices"]
    assert set(sim) == set(FLAGS)


def test_second_import_is_idempotent(tmp_path):
    path = _make(tmp_path, lambda c: c["SYSTEM"].pop("simulated.devices"))
    assert _import_m4(path).returncode == 0
    first = path.read_bytes()
    assert _import_m4(path).returncode == 0
    assert path.read_bytes() == first


@pytest.mark.xfail(
    reason="opticalib.core.root does `_config['SYSTEM']` at import (KeyError) "
    "before m4 gets to run; m4/__init__.py would raise the same KeyError",
    strict=False,
)
def test_missing_system_section(tmp_path):
    """Documents current behaviour when SYSTEM is absent altogether."""
    path = _make(tmp_path, lambda c: c.pop("SYSTEM"))
    res = _import_m4(path)
    assert res.returncode == 0, res.stderr[-500:]

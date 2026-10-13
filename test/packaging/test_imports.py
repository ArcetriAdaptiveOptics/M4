"""Every live module under ``m4`` must import cleanly (each in a fresh process)."""
import os
import pkgutil
import subprocess
import sys

import pytest

import m4

EXCLUDED = {
    # Script that executes at import time and needs data files.
    "m4.howto_ESO": "executes a measurement script at import, needs data files",
}


def _modules():
    names = []
    for info in pkgutil.walk_packages(m4.__path__, "m4."):
        names.append(info.name)
    return sorted(names)


MODULES = _modules()


def test_modules_discovered():
    assert len(MODULES) > 10
    assert "m4.configuration.ott" in MODULES


def test_no_deprecated_module_is_discovered():
    assert not [n for n in MODULES if "deprecated" in n.split(".")]


@pytest.mark.parametrize("name", MODULES)
def test_module_imports(name):
    if name in EXCLUDED:
        pytest.skip(EXCLUDED[name])
    res = subprocess.run(
        [sys.executable, "-c", f"import {name}"],
        capture_output=True,
        text=True,
        timeout=120,
        env=os.environ.copy(),
    )
    assert res.returncode == 0, f"import {name} failed:\n{res.stderr}"

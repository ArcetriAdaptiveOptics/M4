"""Static checks of pyproject.toml against the code base."""
import ast
import importlib
import re
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PYPROJECT = ROOT / "pyproject.toml"


@pytest.fixture(scope="module")
def cfg():
    with open(PYPROJECT, "rb") as f:
        return tomllib.load(f)


def test_legacy_packaging_files_are_gone():
    for name in ("setup.py", "setup.cfg", "requirements.txt"):
        assert not (ROOT / name).exists(), name


def test_version_is_dynamic_and_well_formed(cfg):
    proj = cfg["project"]
    assert "version" not in proj
    assert "version" in proj["dynamic"]
    assert cfg["tool"]["setuptools"]["dynamic"]["version"] == {
        "attr": "m4.__version__.__version__"
    }
    from m4.__version__ import __version__

    assert re.fullmatch(r"\d+\.\d+\.\d+", __version__)


@pytest.mark.parametrize("section", ["gui-scripts"])
def test_entry_points_resolve(cfg, section):
    eps = cfg["project"][section]
    assert eps
    for name, target in eps.items():
        mod, func = target.split(":")
        fn = getattr(importlib.import_module(mod), func)
        assert callable(fn), name


def test_expected_entry_points(cfg):
    assert cfg["project"]["gui-scripts"]["m4gui"] == "m4.gui.m4gui:main"


def test_package_discovery(cfg):
    assert cfg["tool"]["setuptools"]["packages"]["find"]["include"] == ["m4*"]


def test_certificates_are_not_package_data(cfg):
    # OPC-UA certificates are read from the source tree (see opc_ua_controller),
    # so they must not be shipped: m4 is meant to be installed with `pip install -e .`
    package_data = cfg["tool"]["setuptools"].get("package-data", {})
    for pkg, patterns in package_data.items():
        for pat in patterns:
            assert not pat.endswith((".pem", ".der")), (pkg, pat)


def test_certificates_referenced_by_opc_ua_controller_exist():
    src = (ROOT / "m4" / "devices" / "opc_ua_controller.py").read_text()
    for name in re.findall(r"[\w\-]+\.(?:der|pem)", src):
        assert (ROOT / "m4" / "devices" / name).exists(), name


def test_pytest_testpaths(cfg):
    assert cfg["tool"]["pytest"]["ini_options"]["testpaths"] == ["test"]


# ---- dependency coverage ---------------------------------------------------
# import name -> distribution name (normalised later)
DIST = {
    "zmq": "pyzmq",
    "skimage": "scikit-image",
    "ruamel": "ruamel.yaml",
    "yaml": "pyyaml",
    "serial": "pyserial",
    "playsound": "playsound",
    "PyQt5": "PyQt5",
    "IPython": "ipython",
    "opcua": "opcua",
    "cryptography": "cryptography",
    "guietta": "guietta",
    "mock": "mock",
    "pytest": "pytest",
}
# Internal lab library (Microgate), not distributable on PyPI.
ALLOWLIST = {"Microgate"}


def _norm(s):
    return re.sub(r"[-_.]+", "-", s).lower()


def _req_names(reqs):
    out = set()
    for r in reqs:
        m = re.match(r"\s*([A-Za-z0-9_.\-]+)", r)
        out.add(_norm(m.group(1)))
    return out


def _third_party_imports():
    """Return {top-level name: [(file, lineno)]} for all non-optional imports,
    and the same dict for optional ones (inside try/except ImportError or functions)."""
    hard, soft = {}, {}

    def catches_import_error(node):
        for h in node.handlers:
            t = h.type
            names = []
            if t is None:
                return True
            elts = t.elts if isinstance(t, ast.Tuple) else [t]
            for e in elts:
                names.append(getattr(e, "id", getattr(e, "attr", "")))
            if {"ImportError", "ModuleNotFoundError", "Exception"} & set(names):
                return True
        return False

    def visit(node, optional, path):
        for child in ast.iter_child_nodes(node):
            opt = optional
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                opt = True
            if isinstance(child, ast.Try) and catches_import_error(child):
                for b in child.body:
                    visit_one(b, True, path)
                for part in (child.handlers, child.orelse, child.finalbody):
                    for b in part:
                        visit_one(b, optional, path)
                continue
            visit_one(child, opt, path)

    def visit_one(node, optional, path):
        mods = []
        if isinstance(node, ast.Import):
            mods = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            mods = [node.module]
        for m in mods:
            top = m.split(".")[0]
            if top in sys.stdlib_module_names or top == "m4":
                continue
            (soft if optional else hard).setdefault(top, []).append(
                (str(path.relative_to(ROOT)), node.lineno)
            )
        visit(node, optional, path)

    for p in (ROOT / "m4").rglob("*.py"):
        tree = ast.parse(p.read_text(encoding="utf-8"), filename=str(p))
        visit(tree, False, p)
    return hard, soft


def test_third_party_imports_are_declared(cfg):
    proj = cfg["project"]
    declared = _req_names(proj["dependencies"])
    extras = set()
    for reqs in proj["optional-dependencies"].values():
        extras |= _req_names(reqs)
    # `opticalib[gui]` style names are normalised by the regex above.
    hard, soft = _third_party_imports()
    problems = []
    for top, where in {**soft, **hard}.items():
        if top in ALLOWLIST:
            continue
        dist = _norm(DIST.get(top, top))
        ok = dist in declared or dist in extras
        if not ok:
            problems.append(f"{top} (dist {dist}) used at {where[0]}")
    assert not problems, "\n".join(problems)


def test_hard_dependencies_are_in_main_dependencies(cfg):
    declared = _req_names(cfg["project"]["dependencies"])
    for must in ("opticalib", "numpy", "scipy", "pyzmq", "ruamel-yaml"):
        assert must in declared

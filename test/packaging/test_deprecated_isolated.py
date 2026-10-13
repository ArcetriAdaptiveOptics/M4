"""Nothing under m4/ may import a module that lives only in deprecated/."""
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEPRECATED = ROOT / "deprecated"
LIVE = ROOT / "m4"


def _dotted(root: Path, path: Path) -> str:
    parts = list(path.relative_to(root).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _deprecated_only_names():
    names = set()
    for p in (DEPRECATED / "m4").rglob("*.py"):
        name = _dotted(DEPRECATED, p)
        if name == "m4":
            continue
        live = ROOT / (name.replace(".", "/") + ".py")
        live_pkg = ROOT / name.replace(".", "/") / "__init__.py"
        if not live.exists() and not live_pkg.exists():
            names.add(name)
    return names


def _imports(path: Path):
    """Yield (lineno, dotted name) for every import, relative ones resolved."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    pkg = _dotted(ROOT, path).split(".")
    if path.name != "__init__.py":
        pkg = pkg[:-1]
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                yield node.lineno, a.name
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = pkg[: len(pkg) - (node.level - 1)] if node.level > 1 else pkg
                mod = ".".join(base + ([node.module] if node.module else []))
            else:
                mod = node.module or ""
            yield node.lineno, mod
            for a in node.names:
                yield node.lineno, f"{mod}.{a.name}"


def test_deprecated_folder_exists_and_has_modules():
    assert len(_deprecated_only_names()) > 5


def test_deprecated_is_not_a_package():
    assert not (DEPRECATED / "__init__.py").exists()
    assert not (DEPRECATED / "m4" / "__init__.py").exists()


def test_live_code_does_not_import_deprecated_modules():
    dead = _deprecated_only_names()
    offenders = []
    for p in LIVE.rglob("*.py"):
        for lineno, name in _imports(p):
            if name in dead:
                offenders.append(f"{p.relative_to(ROOT)}:{lineno} imports {name}")
    assert not offenders, "\n".join(offenders)


def test_live_code_does_not_import_a_deprecated_top_level_package():
    offenders = []
    for p in LIVE.rglob("*.py"):
        for lineno, name in _imports(p):
            if name == "deprecated" or name.startswith("deprecated."):
                offenders.append(f"{p.relative_to(ROOT)}:{lineno}")
    assert not offenders, "\n".join(offenders)


def test_scanner_resolves_relative_imports(tmp_path):
    f = tmp_path / "m4" / "sub" / "mod.py"
    f.parent.mkdir(parents=True)
    f.write_text("from . import x\nfrom ..y import z\n# import m4.dead\ns='import m4.dead'\n")
    global ROOT
    old, ROOT = ROOT, tmp_path
    try:
        got = {n for _, n in _imports(f)}
    finally:
        ROOT = old
    assert {"m4.sub", "m4.sub.x", "m4.y", "m4.y.z"} <= got
    assert "m4.dead" not in got

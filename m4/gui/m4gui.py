"""
Launcher of the M4 GUI.

Usage::

    m4gui -f <experiment folder or configuration.yaml>

``m4`` (and opticalib) read the configuration when they are imported, so
``AOCONF`` is set here before importing either of them.
"""

import argparse
import os


def _resolve_experiment(path: str) -> str:
    """The configuration file of an experiment folder (or the file itself)."""
    path = os.path.abspath(os.path.expanduser(path))
    if os.path.isfile(path):
        return path
    for candidate in (
        os.path.join(path, "SysConfig", "configuration.yaml"),
        os.path.join(path, "configuration.yaml"),
    ):
        if os.path.isfile(candidate):
            return candidate
    raise SystemExit(f"No configuration file found in {path}")


def main() -> None:
    """Parse the command line, select the configuration and open the GUI."""
    parser = argparse.ArgumentParser(description="M4 CalpyGUI")
    parser.add_argument(
        "-f", "--config", help="experiment folder or configuration.yaml (default: $AOCONF)"
    )
    args = parser.parse_args()
    config_path = None
    if args.config:
        config_path = _resolve_experiment(args.config)
        os.environ["AOCONF"] = config_path

    from m4.gui.calpygui import launch_gui

    launch_gui(config_path)


if __name__ == "__main__":
    main()

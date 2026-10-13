"""
The M4 main window: opticalib's GUI with the M4 procedures and environment.
"""

import os
from typing import Optional

# opticalib.gui first: it selects the Qt binding before qtpy is imported.
from opticalib.gui import PROCEDURES, CalpyGUI, Procedure
from opticalib.gui.app import launch_gui as _launch_gui

from .procedures import OttScriptsWindow

_HERE = os.path.dirname(os.path.abspath(__file__))


class M4GUI(CalpyGUI):
    """
    CalpyGUI for the M4 Optical Test Tower.

    Only class attributes are overridden: everything else (console, devices,
    data browser, viewer, ...) is opticalib's.
    """

    APP_NAME = "M4 CalpyGUI"
    # Own layout and recent experiments; theme preferences stay shared.
    SETTINGS_APP = "M4CalpyGUI"
    DESKTOP_FILE_NAME = None
    # Run in the kernel after opticalib's initCalpy.py, and after restarts.
    INIT_SCRIPTS = [os.path.join(_HERE, "initM4GUI.py")]
    # M4's cards first, then opticalib's (drop them to show only M4's).
    PROCEDURES = [
        Procedure(
            "OTT user scripts",
            "telescope",
            "Set up the OTT, move the reference mirror, align (m4.userscripts).",
            OttScriptsWindow,
        ),
    ] + PROCEDURES


def launch_gui(config_path: Optional[str] = None) -> None:
    """
    Open the M4 GUI (blocks until it is closed).

    ``m4`` reads its configuration when it is imported: set ``AOCONF`` to
    *config_path* before importing this module (``m4.gui.m4gui``
    does it).

    Parameters
    ----------
    config_path : str, optional
        The ``configuration.yaml`` file (default: ``$AOCONF``).
    """
    _launch_gui(config_path, window_class=M4GUI)

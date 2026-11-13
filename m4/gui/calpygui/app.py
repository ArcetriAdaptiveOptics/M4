"""
The M4 main window: opticalib's GUI with the M4 procedures and environment.
"""

import os
from typing import Optional

# opticalib.gui first: it selects the Qt binding before qtpy is imported.
from opticalib.gui import CalpyGUI
from opticalib.gui.app import launch_gui as _launch_gui

from .procedures import M4ACTIONS, M4CALIBRATIONS, MEASPANEL

_HERE = os.path.dirname(os.path.abspath(__file__))


class M4GUI(CalpyGUI):
    """
    CalpyGUI for the M4 Optical Test Tower.

    Only class attributes and :meth:`add_panels` are overridden: everything
    else (console, devices, data browser, viewer, ...) is opticalib's.
    """

    APP_NAME = "M4 CalpyGUI"
    # Own layout and recent experiments; theme preferences stay shared.
    SETTINGS_APP = "M4CalpyGUI"
    DESKTOP_FILE_NAME = None
    # Run in the kernel after opticalib's initCalpy.py, and after restarts.
    INIT_SCRIPTS = [os.path.join(_HERE, "initM4GUI.py")]
    # The "Procedures" panel: M4's actions only (add opticalib's PROCEDURES
    # to the list to show them too).
    PROCEDURES = M4ACTIONS + M4CALIBRATIONS

    def add_panels(self) -> None:
        """A second panel of cards, in a tab next to "Procedures"."""
        self.measurements_panel = self.add_procedure_panel(
            "meas",
            "Measurements",
            MEASPANEL,
            caption="Calibrate the OTT alignments and the DM.",
        )


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

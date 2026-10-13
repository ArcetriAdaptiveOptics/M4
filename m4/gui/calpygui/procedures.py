"""
Procedure windows of the M4 GUI
===============================

Each window lists *steps*; each step is a form (its parameters) and a code
template.  Running a step runs the generated code in the GUI console, so
what the GUI does is exactly what you would type in the shell.

:class:`OttScriptsWindow` implements a few functions of
:mod:`m4.userscripts` (``OTTScripts``), as an example to extend.
"""

from typing import List

from opticalib.gui.procedures.base import ProcedureWindow, Step
from opticalib.gui.procedures.params import (
    BoolParam,
    ExprParam,
    IntParam,
    VarParam,
)

#: Kernel variable of the OTTScripts object, as created by the setup step.
SCRIPTS = "ottscripts"


def _scripts_param() -> ExprParam:
    """The parameter naming the ``OTTScripts`` object (shared by the steps)."""
    return ExprParam(
        "scripts", "User scripts object", SCRIPTS, help="Created by 'Load the user scripts'."
    )


class OttScriptsWindow(ProcedureWindow):
    """The ``OTTScripts`` functions of ``m4.userscripts``."""

    TITLE = "OTT user scripts"
    DESCRIPTION = (
        "Set up the Optical Test Tower and run the user-level scripts of "
        "m4.userscripts: reference mirror, beam position and alignment. "
        "Every step runs, visibly, in the console."
    )
    #: Shown on top of the window while false in the kernel.
    CHECKS = [
        (
            f"'{SCRIPTS}' in globals()",
            f"No '{SCRIPTS}' in the console yet: run 'Create the OTT' and "
            "'Load the user scripts' first.",
        ),
    ]

    def steps(self) -> List[Step]:
        """The steps, grouped in sections."""
        return [
            # -- Setup ---------------------------------------------------
            Step(
                "create_ott",
                "Create the OTT",
                "Create the OTT devices, real or simulated as set in "
                "SYSTEM → simulated.devices of the configuration file "
                "(m4.create_ott).",
                [
                    VarParam("ott", "OTT variable", "ott"),
                    VarParam("dm", "DM variable", "dm"),
                    VarParam("interf", "Interferometer variable", "interf"),
                    BoolParam("no_dm", "Without the DM", False),
                ],
                lambda v: (
                    f"{v['ott']}, {v['dm']}, {v['interf']} = "
                    f"m4.create_ott(no_dm={v['no_dm']})"
                ),
                section="Setup",
            ),
            Step(
                "load_scripts",
                "Load the user scripts",
                "Create the OTTScripts object of m4.userscripts, which reads "
                "userconfiguration.yaml from the SysConfig folder.",
                [
                    ExprParam("ott", "OTT", "ott"),
                    # A plain variable name: M4's interferometer classes need
                    # not look like opticalib's to the device selector.
                    ExprParam("interf", "Interferometer", "interf"),
                    VarParam("out", "Store it in", SCRIPTS),
                ],
                lambda v: (
                    "from m4 import userscripts\n"
                    f"{v['out']} = userscripts.OTTScripts({v['ott']}, {v['interf']})"
                ),
                # Re-run the checks: the warning disappears once it exists.
                after=lambda window, results, values: window.run_checks(),
                section="Setup",
            ),
            # -- Reference mirror ------------------------------------------
            Step(
                "deploy_rm",
                "Deploy the reference mirror",
                "Move the reference mirror into the beam, at the alignment "
                "position of userconfiguration.yaml (OTTMECH → "
                "rmslider4alignment).",
                [_scripts_param()],
                lambda v: f"{v['scripts']}.deployReferenceMirror()",
                confirm="This moves the reference mirror slider.",
                section="Reference mirror",
            ),
            Step(
                "retract_rm",
                "Retract the reference mirror",
                "Move the reference mirror out of the parabola footprint "
                "(OTTMECH → rmsliderout), to measure M4.",
                [_scripts_param()],
                lambda v: f"{v['scripts']}.retractReferenceMirror()",
                confirm="This moves the reference mirror slider.",
                section="Reference mirror",
            ),
            Step(
                "positions",
                "Where are they?",
                "Read the positions of the reference mirror and of the "
                "parabola truss, in M4-centred coordinates.",
                [_scripts_param()],
                lambda v: (
                    f"rm_position = {v['scripts']}.whereisRef()\n"
                    f"truss_position = {v['scripts']}.whereisBeam()"
                ),
                # Read back after the step and shown in the window's log.
                outputs={"rm_position": "rm_position", "truss_position": "truss_position"},
                section="Reference mirror",
            ),
            # -- Alignment --------------------------------------------------
            Step(
                "align_ott",
                "Align the OTT",
                "Measure the misalignment of the OTT (tip, tilt, focus and "
                "coma through the parabola and the reference mirror) and, if "
                "asked, apply the correction (OTTScripts.alignOTT).",
                [
                    _scripts_param(),
                    IntParam("nframes", "Frames to average", 4, minimum=1, maximum=1000),
                    BoolParam("move", "Apply the correction", False),
                    BoolParam("remove_par", "Remove the parabola", True),
                ],
                lambda v: (
                    f"{v['scripts']}.alignOTT({v['nframes']}, "
                    f"move={int(v['move'] == 'True')}, removePar={v['remove_par']})"
                ),
                confirm="With 'Apply the correction' this moves the parabola and "
                "the reference mirror.",
                section="Alignment",
            ),
        ]

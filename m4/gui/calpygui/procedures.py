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
    ChoiceParam
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
                    VarParam("interf", "Interferometer variable", "interf"),
                    VarParam('myott', 'OTT scripts class variable', SCRIPTS),
                ],
                lambda v: (
                    f"{v['ott']}, {v['interf']} = "
                    f"m4.create_ott(no_dm=True)\n"
                    f"{v['myott']} = "
                    f"m4.userscripts.OTTScripts({v['ott']}, {v['interf']})"
                ),
                section="Setup",
                after=lambda window, results, values: window.run_checks(),
            ),
            # Step(
            #     "load_scripts",
            #     "Load the user scripts",
            #     "Create the OTTScripts object of m4.userscripts, which reads "
            #     "userconfiguration.yaml from the SysConfig folder.",
            #     [
            #         ExprParam("ott", "OTT", "ott"),
            #         # A plain variable name: M4's interferometer classes need
            #         # not look like opticalib's to the device selector.
            #         ExprParam("interf", "Interferometer", "interf"),
            #         VarParam("out", "Store it in", SCRIPTS),
            #     ],
            #     lambda v: (
            #         "from m4 import userscripts\n"
            #         f"{v['out']} = userscripts.OTTScripts({v['ott']}, {v['interf']})"
            #     ),
            #     # Re-run the checks: the warning disappears once it exists.
            #     after=lambda window, results, values: window.run_checks(),
            #     section="Setup",
            # ),
            # --- 4D Operations ------------------------------------------------
            Step(
                "view_ott",
                "Configure the OTT view.",
                "Set the current OTT's interferometer configuration "
                "to look either at the OTT cavity (parabola) or the DM.",
                [
                    _scripts_param(),
                    ChoiceParam(
                        "view",
                        "OTT view",
                        choices=[("cavity", "cavity"), ("dm", "DM")], 
                        default="cavity"
                    ),
                ],
                lambda v: f"{v['scripts']}.config4D4Alignment(cavity_or_dm='{v['view']}')",
                section="4D Operations",
            ),
            Step(
                "view_dm_segment",
                "Configure segment's view.",
                "Set the current OTT's interferometer configuration "
                "to look at the selected DM segment.",
                [
                    _scripts_param(),
                    IntParam(
                        "view",
                        "DM segment",
                        default=0,
                        minimum=0,
                        maximum=5
                    ),
                ],
                lambda v: f"{v['scripts']}.config4D4Segment(segment={v['view']})",
                section="4D Operations",
            ),
            ## --- OTT Mechanical Operations -----------------------------------
            Step(
                "insert_rm",
                "Insert Reference Mirror",
                "Move the reference mirror into the beam, at the alignment "
                "position of userconfiguration.yaml (OTTMECH → "
                "rmslider4alignment).",
                [_scripts_param()],
                lambda v: f"{v['scripts']}.insertReferenceMirror()",
                confirm="This moves the reference mirror slider.",
                section="Ott Mechanical Operations"
            ),
            Step(
                "rm_out",
                "Retract Reference Mirror",
                "Move the reference mirror out of the beam, at the alignment "
                "position of userconfiguration.yaml (OTTMECH → "
                "rmsliderout).",
                [_scripts_param()],
                lambda v: f"{v['scripts']}.retractReferenceMirror()",
                confirm="This moves the reference mirror slider.",
                section="Ott Mechanical Operations"
            )
        ]
# Insert RM 	configureOTT4Alignment 	rmslider4alignment 	  	 
# Retract RM 	configureOTTrefMirrorOut 	rmsliderout 	  	 
# Deploy optBeam 	configureOTT4Segment 	parslider4segment: 	  	 
# Rotate optBeam
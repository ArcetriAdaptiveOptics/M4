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
    FloatParam,
    ChoiceParam,
    TnParam,
    VarParam
)
from opticalib.gui import Procedure


_IFFUNC = "__import__('opticalib.core.config', fromlist=['_']).get_iff_config('IFFUNC')"

#: Kernel variables of the OTTScripts and M4Scripts objects.
SCRIPTS = ["ottscripts", "m4scripts", "spl", "ottmeas"]


def _scripts_param(default: str) -> ExprParam:
    """
    The parameter naming the user-scripts object of a step.

    Parameters
    ----------
    default : str
        Its kernel variable: ``SCRIPTS[0]`` (``OTTScripts``) or ``SCRIPTS[1]``
        (``M4Scripts``).  It must be a single name: a text parameter's default
        is the text shown in its field.
    """
    return ExprParam("scripts", "User scripts object", default)


def _modes(values) -> str:
    """The modes typed in, or those of the configuration when left empty."""
    return f"{_IFFUNC}['modes_list']" if values["modes"] == "None" else values["modes"]


class OttActions(ProcedureWindow):
    """The ``OTTScripts`` functions of ``m4.userscripts``."""

    TITLE = "OTT Actions"
    DESCRIPTION = (
        "Set of actions for moving and preparing the OTT for specific tasks."
    )
    #: Shown on top of the window while false in the kernel.
    CHECKS = [
        (
            f"'{SCRIPTS[0]}' in globals()",
            f"No '{SCRIPTS[0]}' in the console yet: connect the 'OttScripts' first.",
        ),
    ]

    def steps(self) -> List[Step]:
        """The steps, grouped in sections."""
        return [
            Step(
                "deploy_rm",
                "Deploy Reference Mirror",
                "Move the reference mirror into the beam, at the alignment "
                "position of userconfiguration.yaml (OTTMECH → "
                "rmslider4alignment).",
                [_scripts_param(SCRIPTS[0])],
                lambda v: f"{v['scripts']}.deployReferenceMirror()",
                confirm="This moves the reference mirror slider.",
                section="Actions"
            ),
            Step(
                "retract_rm",
                "Retract Reference Mirror",
                "Move the reference mirror out of the beam, at the alignment "
                "position of userconfiguration.yaml (OTTMECH → "
                "rmsliderout).",
                [_scripts_param(SCRIPTS[0])],
                lambda v: f"{v['scripts']}.retractReferenceMirror()",
                confirm="This moves the reference mirror slider.",
                section="Actions"
            ),
            Step(
                "view_ott",
                "Configure the OTT view.",
                "Set the current OTT's interferometer configuration "
                "to look either at the OTT cavity (parabola) or the DM.",
                [
                    _scripts_param(SCRIPTS[0]),
                    ChoiceParam(
                        "view",
                        "OTT view",
                        choices=[("cavity", "cavity"), ("dm", "DM")], 
                        default="cavity"
                    ),
                ],
                lambda v: f"{v['scripts']}.config4D4Alignment(cavity_or_dm={v['view']})",
                section="Actions",
            ),
            Step(
                "view_dm_segment",
                "Configure segment's view.",
                "Set the current OTT's interferometer configuration "
                "to look at the selected DM segment.",
                [
                    _scripts_param(SCRIPTS[0]),
                    IntParam(
                        "view",
                        "DM segment",
                        default=0,
                        minimum=0,
                        maximum=5
                    ),
                ],
                lambda v: f"{v['scripts']}.config4D4Segment(segment={v['view']})",
                section="Actions",
            ),
            Step(
                "rotate_to_Seg",
                "Rotate to Segment.",
                "Move the OTT rotator to look at the selected DM segment.",
                [
                    _scripts_param(SCRIPTS[0]),
                    IntParam(
                        "view",
                        "DM segment",
                        default=0,
                        minimum=0,
                        maximum=5
                    ),
                ],
                lambda v: f"{v['scripts']}.rotateOTT4Segment(segment={v['view']})",
                section="Actions",
            ),
        ]

class OttCalibration(ProcedureWindow):
    """The ``OTTScripts`` functions of ``m4.userscripts``."""

    TITLE = "OTT Calibration"
    DESCRIPTION = (
        "Set of actions for calibrating the optical components of the OTT, i.e. "
        "Parabola, the Reference Mirror and the M4 Hexapode."
    )
    #: Shown on top of the window while false in the kernel.
    CHECKS = [
        (
            f"'{SCRIPTS[0]}' in globals()",
            f"No '{SCRIPTS[0]}' in the console yet: connect the 'OttScripts' first.",
        ),
    ]

    def steps(self) -> List[Step]:
        """The steps, grouped in sections."""
        return [
            Step(
                "ott_calib",
                "Calibrate OTT Alignment.",
                "Runs the alignment calibration procedure for the OTT only "
                "(Parabola+RM)",
                [
                    _scripts_param(SCRIPTS[0]),
                    IntParam(
                        'nframes',
                        'Number of frames',
                        default=25,
                        minimum=1,
                        maximum=100
                    ),
                    BoolParam(
                        'save',
                        'Save Calibration',
                        default=True
                    ),
                    VarParam('out', "Store the TN in", "tn_ottcal")
                 ],
                lambda v: f"{v['out']} = {v['scripts']}.calibrateOTTAlignment(n_frames={v['nframes']}, save={v['save']})",
                section="Acquire",
            ),
            Step(
                "m4_calib",
                "Calibrate M4 Alignment.",
                "Runs the alignment calibration procedure for the M4 Hexapode.",
                [
                    _scripts_param(SCRIPTS[0]),
                    IntParam(
                        'nframes',
                        'Number of frames',
                        default=25,
                        minimum=1,
                        maximum=100
                    ),
                    BoolParam(
                        'save',
                        'Save Calibration',
                        default=True
                    ),
                    VarParam('out', "Store the TN in", "tn__m4cal")
                 ],
                lambda v: f"{v['out']} = {v['scripts']}.calibrateM4Alignment(n_frames={v['nframes']}, save={v['save']})",
                section="Acquire",
            ),
            Step(
                "focus",
                "Correct Focus.",
                "Corrects the alignment focus mode of the OTT.",
                [
                    _scripts_param(SCRIPTS[0]),
                    IntParam(
                        "nframes",
                        "Number of frames",
                        default=25,
                        minimum=1,
                        maximum=100                        
                    ),
                    BoolParam(
                        'mov',
                        'Move',
                        default=False
                    ),
                    BoolParam(
                        'rpar',
                        'Remove Parabola',
                        default=True
                    )
                ],
                lambda v: f"{v['scripts']}.alignFocus(nframes={v['nframes']}, "
                          f"move={v['mov']}, removePar={v['rpar']})",
                section="Correct",
            ),
            Step(
                "coma",
                "Correct Coma.",
                "Corrects the alignment coma mode of the OTT.",
                [
                    _scripts_param(SCRIPTS[0]),
                    IntParam(
                        "nframes",
                        "Number of frames",
                        default=25,
                        minimum=1,
                        maximum=100                        
                    ),
                    BoolParam(
                        'mov',
                        'Move',
                        default=False
                    ),
                    BoolParam(
                        'rpar',
                        'Remove Parabola',
                        default=True
                    )
                ],
                lambda v: f"{v['scripts']}.alignComa(nframes={v['nframes']}, "
                          f"move={v['mov']}, removePar={v['rpar']})",
                section="Correct",
            ),
            Step(
                "tt",
                "Correct OTT's Tip Tilt.",
                "Corrects the alignment Tip Tilt modes of the OTT.",
                [
                    _scripts_param(SCRIPTS[0]),
                    IntParam(
                        "nframes",
                        "Number of frames",
                        default=25,
                        minimum=1,
                        maximum=100                        
                    ),
                    BoolParam(
                        'mov',
                        'Move',
                        default=False
                    ),
                    BoolParam(
                        'rpar',
                        'Remove Parabola',
                        default=True
                    )
                ],
                lambda v: f"{v['scripts']}.alignOTT(nframes={v['nframes']}, "
                          f"move={v['mov']}, removePar={v['rpar']})",
                section="Correct",
            ),
            Step(
                "m4tt",
                "Correct M4's Tip Tilt.",
                "Corrects the alignment Tip Tilt modes of the M4's hexapode.",
                [
                    _scripts_param(SCRIPTS[0]),
                    IntParam(
                        "nframes",
                        "Number of frames",
                        default=25,
                        minimum=1,
                        maximum=100                        
                    ),
                    BoolParam(
                        'mov',
                        'Move',
                        default=False
                    ),
                    BoolParam(
                        'rpar',
                        'Remove Parabola',
                        default=True
                    )
                ],
                lambda v: f"{v['scripts']}.alignM4TT(nframes={v['nframes']}, "
                          f"move={v['mov']}, removePar={v['rpar']})",
                section="Correct",
            ),
        ]

class DMActions(ProcedureWindow):
    """The ``OTTScripts`` functions of ``m4.userscripts``."""

    TITLE = "DM Actions"
    DESCRIPTION = (
        "Set up the Optical Test Tower and run the user-level scripts of "
        "m4.userscripts: reference mirror, beam position and alignment. "
        "Every step runs, visibly, in the console."
    )
    #: Shown on top of the window while false in the kernel.
    CHECKS = [
        (
            f"'{SCRIPTS[1]}' in globals()",
            f"No '{SCRIPTS[1]}' in the console yet: connect the 'OttScripts' first.",
        ),
    ]

    def steps(self) -> List[Step]:
        """The steps, grouped in sections."""
        return [
            Step(
                "relax",
                "Relax DM.",
                "Sets the DM shape to the its bias shape.",
                [_scripts_param(SCRIPTS[1])],
                lambda v: f"{v['scripts']}.relax()",
                section="Actions",
            ),
            Step(
                "loadflat",
                "Load Flat Command.",
                "Loads and apply a flattening command given a Tracking Number",
                [
                    _scripts_param(SCRIPTS[1]),
                    TnParam(
                        'tn',
                        'Tracking Number',
                        folder_attr="FLAT_ROOT_FOLDER",
                        optional=True,
                    )
                ],
                lambda v: f"{v['scripts']}.loadFlatCommand(flattn={v['tn']})",
                section="Actions",
            )
        ]
        
class DMCalibration(ProcedureWindow):
    """The ``OTTScripts`` functions of ``m4.userscripts``."""

    TITLE = "DM Calibration"
    DESCRIPTION = (
        "Set up the Optical Test Tower and run the user-level scripts of "
        "m4.userscripts: reference mirror, beam position and alignment. "
        "Every step runs, visibly, in the console."
    )
    #: Shown on top of the window while false in the kernel.
    CHECKS = [
        (
            f"'{SCRIPTS[1]}' in globals()",
            f"No '{SCRIPTS[1]}' in the console yet: connect the 'OttScripts' first.",
        ),
    ]

    def steps(self) -> List[Step]:
        """The steps, grouped in sections."""
        return [
            Step(
                "iffacq_modal",
                "Acquire Modal IFF.",
                "Runs and acquire DM Influence Functions with the Mirror's modal base.",
                [
                    _scripts_param(SCRIPTS[1]),
                    ExprParam(
                        'modes',
                        'Modes',
                        "",
                        optional=True,  # empty: the modes of the configuration
                        config_default=f"{_IFFUNC}['modes_list']"
                    ),
                    IntParam(
                        "seg",
                        "Segment",
                        default=0,
                        minimum=0,
                        maximum=5
                    ),
                    IntParam(
                        "npushpull",
                        "Push-pull steps",
                        default=3,
                        minimum=1,
                        help="Length of the +1/-1 template (3 = [1, -1, 1]).",
                    ),
                    FloatParam(
                        "amp",
                        "Amplitude",
                        None,
                        optional=True,
                        config_default=f"{_IFFUNC}['amplitude']"
                    ),
                    IntParam(
                        "nrep",
                        "Number of repetitions",
                        default=1,
                        minimum=1,
                        optional=True
                    ),
                    BoolParam(
                        "shuffle",
                        "Shuffle",
                        default=False,
                        optional=True,
                    ),
                    BoolParam(
                        'view',
                        "View",
                        default=True,
                        optional=True
                    ),
                    VarParam('out', "Store the TN in", "tn_iff")
                ],
                lambda v: f"{v['out']} = {v['scripts']}.acquireModalIFF(modes={_modes(v)}, "
                    f"segment={v['seg']}, npushpull={v['npushpull']}, "
                    f"n_repetitions={v['nrep']}, amp={v['amp']}, "
                    f"shuffle={v['shuffle']}, view={v['view']})",
                section="Acquire",
            ),
            Step(
                "iffacq_zonal",
                "Acquire Zonal IFF.",
                "Runs and acquire DM Influence Functions with the Zonal modal base.",
                [
                    _scripts_param(SCRIPTS[1]),
                    ExprParam(
                        'modes',
                        'Modes',
                        "",
                        optional=True,  # empty: the modes of the configuration
                        config_default=f"{_IFFUNC}['modes_list']"
                    ),
                    IntParam(
                        "seg",
                        "Segment",
                        default=0,
                        minimum=0,
                        maximum=5
                    ),
                    IntParam(
                        "npushpull",
                        "Push-pull steps",
                        default=3,
                        minimum=1,
                        help="Length of the +1/-1 template (3 = [1, -1, 1]).",
                    ),
                    FloatParam(
                        "amp",
                        "Amplitude",
                        None,
                        optional=True,
                        config_default=f"{_IFFUNC}['amplitude']"
                    ),
                    IntParam(
                        "nrep",
                        "Number of repetitions",
                        default=1,
                        minimum=1,
                        optional=True
                    ),
                    BoolParam(
                        "shuffle",
                        "Shuffle",
                        default=False,
                        optional=True,
                    ),
                    BoolParam(
                        'view',
                        "View",
                        default=True,
                        optional=True
                    ),
                    VarParam('out', "Store the TN in", "tn_iff")
                ],
                lambda v: f"{v['out']} = {v['scripts']}.acquireZonalIFF(modes={_modes(v)}, "
                    f"segment={v['seg']}, npushpull={v['npushpull']}, "
                    f"n_repetitions={v['nrep']}, amp={v['amp']}, "
                    f"shuffle={v['shuffle']}, view={v['view']})",
                section="Acquire",
            ),
        ]

#: Cards of the "Procedures" panel.
M4ACTIONS = [
    Procedure("OTT Actions", "tower-fire", "Move and prepare the OTT.", OttActions),
    Procedure("M4 Actions", "mirror-variant", "Relax the DM, load a flat.", DMActions),
]

#: Cards of the "Calibrations" panel.
M4CALIBRATIONS = [
    Procedure(
        "OTT Calibrations",
        "tower-beach",
        "Calibrate and correct the OTT and M4 alignment.",
        OttCalibration,
    ),
    Procedure(
        "M4 Calibrations",
        "hexagon-slice-1",
        "Acquire the DM influence functions.",
        DMCalibration,
    ),
]


class SPL(ProcedureWindow):
    
    TITLE = "SPL"
    DESCRIPTION = (
        "Procedure for acquiring and analyzing the 6-chain SPL data."
    )
    #: Shown on top of the window while false in the kernel.
    CHECKS = [
        (
            f"'{SCRIPTS[2]}' in globals()",
            f"No '{SCRIPTS[2]}' in the console yet: connect the 'SPL' device first.",
        ),
    ]

    def steps(self) -> list[Step]:
        return [
            Step(
                'acquire',
                'Acquire data',
                "Acquires SPL data across the selected wavelenghts.",
                [
                    IntParam(
                        'step',
                        "Wavelenght steps [nm]",
                        default=10,
                        minimum=5,
                        maximum=50
                    ),
                    VarParam('out', "Store the TN in", "tn_spl")
                ],
                lambda v: f"{v['out']} = {SCRIPTS[2]}.acquire(step={v['step']})",
                section='Acquire'
            ),
            Step(
                'analyze',
                'Analyze data',
                "Analyzes the acquired SPL data.",
                [
                    TnParam(
                        'tn',
                        "TN to analyze",
                        "SPL_DATA_ROOT_FOLDER",
                        optional=False
                    ),
                    VarParam('out', "Store the result in", "tn_spl")
                ],
                lambda v: f"{v['out']} = {SCRIPTS[2]}.analyze(tn={v['tn']})",
                section='Processing'
            )
        ]

class TimeSeries(ProcedureWindow):
    
    TITLE = "Time Series"
    
    DESCRIPTION = "Procedure for acquiring and analyze time series data."
    
    CHECKS = [
        (
            f"'{SCRIPTS[3]}' in globals()", 
            f"No '{SCRIPTS[3]}' in the console yet.",
        )
    ]
    
    def steps(self) -> list[Step]:
        return [
            Step(
                'acquire',
                'Acquire timeseries',
                "Acquire a sequence of frames separated by a given time interval.",
                [
                    IntParam(
                        'nframes',
                        "Number of frames",
                        default=10,
                        minimum=2,
                    ),
                    IntParam(
                        'interval',
                        "Time interval [s]",
                        default=1
                    ),
                    VarParam('out', "Store the result in", "ts")
                ],
                lambda v: f"{v['out']} = {SCRIPTS[3]}.acquire(nframes={v['nframes']}, delay={v['interval']})", # TODO
                section='Acquire'
            ),
            Step(
                'noise',
                'Acquire Noise data',
                "Acquire timeseries data for noise analysis.",
                [
                    IntParam(
                        'nframes',
                        "Number of frames",
                        default=10,
                        minimum=2,
                    ),
                    VarParam('out', "Store the result in", "noise_tn")
                ],
                lambda v: f"{v['out']} = {SCRIPTS[3]}.noise(nframes={v['nframes']})", # TODO
                section='Acquire'
            ),
            # Step(
                
            # ),
            # Step(
                
            # ),
            # Step(
                
            # ),
        ]

class OptCavity(ProcedureWindow):
    
    TITLE = "OTT Cavity"
    
    DESCRIPTION = (
        "Procedure for acquiring and calibrate the OTT's optical cavity (Parabola)."
    )
    
    CHECKS = [
        (
            f"'{SCRIPTS[3]}' in globals()", 
            f"No '{SCRIPTS[3]}' in the console yet. Load the `{SCRIPTS[3]}` in devices first.",
        )
    ]
    
    def steps(self) -> list[Step]:
        return [
            Step(
                'acquire_marker',
                'Acquire Markers',
                "Acquires data for the parabola markers.",
                [
                    VarParam('out', "Store the result in", "markers")
                ],
                lambda v: f"{v['out']} = {SCRIPTS[3]}.acquire_markers()",
                section='Markers'
            ),
            Step(
                'calibrate',
                'Calibrate Parabola',
                "Calibrate the parabola markers.",
                [
                    VarParam('out', "Store the result in", "calib")
                ],
                lambda v: f"{v['out']} = {SCRIPTS[3]}.calibrate()",
                section='Markers'
            ),
            Step(
                'load',
                'Load Calibration',
                "Load the parabola calibration.",
                [
                    TnParam('tn', "Calibration tn", "CAVITY_ROOT_FODLER") # TODO
                ],
                lambda v: f"{SCRIPTS[3]}.load_cavity(tn={v['tn']})",
                section='Markers'
            ),
        ]


MEASPANEL = [
    Procedure(
        "SPL",
        "laser-pointer",
        "Acquisition and analysis procedure for the 6-chained SPL devices.",
        SPL
    ),
    Procedure(
        "Time Series",
        "timer-cog-outline",
        "Procedure for time series data acquisition, analysis and noise estimation.",
        TimeSeries
    ),
    Procedure(
        "OTT optical cavity",
        "circle-double",
        "Procedure for calibrating the OTT optical cavity though marker acquisition.",
        OptCavity
    )
]

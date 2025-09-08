"""
Author(s):
----------
- Chiara Selmi: (probably) written in 2020
- Pietro Ferraiuolo: modified in 2024
"""

from os.path import join as _join
import json as _j
import time as _t
import numpy as _np
import configparser as _cp
from opticalib.ground import osutils as _osu

statusfilename = "OTTStatus.ini"


def read_positions(tn: str) -> dict[str, float]:
    """
    Function which reads the positions of the OTT devices from the provided
    status file.

    Parameters
    ----------
    tn : str
        The track number of the test.

    Returns
    -------
    positions : dict
        A dictionary containing the positions of the OTT devices.
    """
    fold = _osu.findTracknum(tn, complete_path=True)
    f2read = _join(fold, tn, statusfilename)
    config = _cp.ConfigParser()
    print(f" Reading {f2read}...")
    config.read(f2read)
    # OTT
    pp = config["OTT"]
    positions = {
        key: (
            _np.array(_j.loads(pp[key]))
            if key in ["PAR", "RM", "M4", "DP"]
            else float(pp[key])
        )
        for key in ["PAR", "RM", "M4", "DP", "PAR_SLIDER", "RM_SLIDER", "ROT_ANGLE"]
    }
    positions = {
        "parabola": positions["PAR"],
        "referenceMirror": positions["RM"],
        "m4Exapode": positions["M4"],
        "dp": positions["DP"],
        "parabolaSlider": positions["PAR_SLIDER"],
        "referenceMirrorSlider": positions["RM_SLIDER"],
        "angleRotator": positions["ROT_ANGLE"],
    }
    return positions


def save_positions(basepath, ott):
    """
    Function which saves the current positions of the OTT devices in a status
    file.

    Parameters
    ----------
    basepath : str
        The path where the status file will be saved.
    ott : object
        The OTT object.

    Returns
    -------
    str
        A message indicating the status file has been saved.
    """
    fname = _join(basepath, statusfilename)
    f = open(fname, "w")
    f.write("[OTT]\n")
    par = ott.parabola.getPosition()
    rm = ott.referenceMirror.getPosition()
    m4 = ott.m4Exapode.getPosition()
    ps = ott.parabolaSlider.getPosition()
    rs = ott.referenceMirrorSlider.getPosition()
    ang = ott.angleRotator.getPosition()
    f.write("PAR_SLIDER = " + str(ps) + "\n")
    f.write("RM_SLIDER  = " + str(rs) + "\n")
    f.write("ROT_ANGLE  = " + str(ang) + "\n")
    f.write(
        "PAR        = "
        + _np.array2string(
            par, separator=",", formatter={"float_kind": lambda x: "%.2f" % x}
        )
        + "\n"
    )
    f.write(
        "RM         = "
        + _np.array2string(
            rm, separator=",", formatter={"float_kind": lambda x: "%.2f" % x}
        )
        + "\n"
    )
    try:
        dp = ott.dp.getPosition()
        f.write(
            "DP         = "
            + _np.array2string(
                dp, separator=",", formatter={"float_kind": lambda x: "%.2f" % x}
            )
            + "\n"
        )
    except AttributeError:
        dp = "DP not available"
        f.write("DP         = " + dp + "\n")
    f.write(
        "M4         = "
        + _np.array2string(
            m4, separator=",", formatter={"float_kind": lambda x: "%.2f" % x}
        )
        + "\n"
    )
    f.close()
    return f"Status saved in {fname}"


def go_to_geometry(tn, ott):
    """
    Reads the saved positions from the status file and moves the OTT sliders (
    Parabola slider and Reference mirror slider) to the corresponding positions.

    Parameters
    ----------
    tn : str
        The track number of the test.
    ott : object
        The OTT object.

    Returns
    -------
    str
        A message indicating the geometry positions have been set.
    """
    import time
    positions = read_positions(tn)
    ps = positions["parabolaSlider"]
    rs = positions["referenceMirrorSlider"]
    ott.parabolaSlider.setPosition(ps)
    print(f"Moving the Parabola Slider to {ps}")
    time.sleep(2.5)  # ?
    ott.referenceMirrorSlider.setPosition(rs)
    print(f"Moving the Reference Mirror Slider to {rs}")
    time.sleep(2.5)  # ?
    return "Geometry positions set."


def go_to_alignment(tn, ott):
    """
    Reads the positions from the status file and moves the OTT alignment devices
    (Parabola, Reference mirror, M4 and DP) to the corresponding positions.

    Parameters
    ----------
    tn : str
        The track number of the test.
    ott : object
        The OTT object.

    Returns
    -------
    str
        A message indicating the alignment positions have been set.
    """
    positions = read_positions(tn)
    par = positions["parabola"]
    rm = positions["referenceMirror"]
    m4 = positions["m4Exapode"]
    dp = positions["dp"]
    ott.parabola.setPosition(par)
    print(f"Moving the Parabola to {par}")
    ott.referenceMirror.setPosition(rm)
    print(f"Moving the Reference Mirror to {rm}")
    # ott.dp.setPosition(dp)
    # print(f"Moving the DP to {dp}")
    # ott.m4Exapode.setPosition(m4)
    # print(f"Moving the M4 to {m4}")
    return "Alignment positions set."

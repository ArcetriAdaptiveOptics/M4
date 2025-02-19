"""
IFF Module
==========
Author(s):
----------
- Pietro Ferraiuolo
- Runa Briguglio

Description:
------------
This module contains the necessary high/user-leve functions to acquire the IFF data, 
given a deformable mirror and an interferometer.
"""

import os
from m4.ground import read_data as rd
from m4.configuration import config_folder_names as fn
from m4.dmutils import iff_acquisition_preparation as ifa
from m4.configuration.read_iffconfig import copyConfingFile


def iffDataAcquisition(
    dm, interf, modesList=None, amplitude=None, template=None, shuffle=False
):
    """
    This is the user-lever function for the acquisition of the IFF data, given a
    deformable mirror and an interferometer.

    Except for the devices, all the arguments are optional, as, by default, the
    values are taken from the `iffConfig.ini` configuration file.

    Parameters
    ----------
    dm: object
        The inizialized deformable mirror object
    interf: object
        The initialized interferometer object to take measurements
    modesList: int | list, array like , optional
        list of modes index to be measured, relative to the command matrix to be used
    amplitude: float , optional
        command amplitude
    template: string , oprional
        template file for the command matrix
    modalBase: string , optional
        identifier of the modal base to be used
    shuffle: bool , optional
        if True, shuffle the modes before acquisition

    Returns
    -------
    tn: string
        The tracking number of the dataset acquired, saved in the OPDImages folder
    """
    ifc = ifa.IFFCapturePreparation(dm)
    tch = ifc.createTimedCmdHistory(modesList, amplitude, template, shuffle)
    info = ifc.getInfoToSave()
    dm.uploadCmdHistory(tch)
    tn = dm.runCmdHistory(interf)
    iffpath = os.path.join(fn.IFFUNCTIONS_ROOT_FOLDER, tn)
    if not os.path.exists(iffpath):
        os.mkdir(iffpath)
    try:
        for key, value in info.items():
            if key == "shuffle":
                with open(os.path.join(iffpath, f"{key}.dat"), "w") as f:
                    f.write(str(value))
            else:
                rd.saveFits_data(
                    os.path.join(iffpath, f"{key}.fits"), value, overwrite=True
                )
    except KeyError as e:
        print(f"KeyError: {key}, {e}")
    copyConfingFile(tn)
    return tn


# def iffCapture(tn):
#     """
#     This function manages the interfacing equence for collecting the IFF data
#     Parameters
#     ----------------
#     tn: string
#         the tracking number in the xxx folder where the cmd history is saved
#     Returns
#     -------
#     """

#     cmdHist = getCmdHist(tn)
#     dm.uploadCmdHist(cmdHist)
#     dm.runCmdHist()
#     print('Now launching the acquisition sequence')
#     start4DAcq(tn)
#     print('Acquisition completed. Dataset tracknum:')
#     print(tn)

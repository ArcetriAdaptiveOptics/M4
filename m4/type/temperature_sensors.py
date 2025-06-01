"""
Authors
  - C. Selmi: written in 2020


    HOW TO USE IT::

        from m4.type import temperature_sensor as ts
        folder = ts.PT_calibration(n_meas)
        ts.analyzer_PT_meas(tt)
        
"""

import os
import time
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration import config_folder_names as config
from m4.configuration.ott_parameters import OpcUaParameters


######## Sensori PT #######
def PT_calibration(n_meas):
    """
    Parameters
    ----------
        n_meas: int
            number of measurement to store

    Returns
    -------
        dove: string
            data file path of measurement
    """
    from m4.ground import tracking_number_folder
    from opcua import Client

    server = OpcUaParameters.server
    client = Client(url=server)
    client.connect()

    folder = config.PT_ROOT_FOLDER
    dove, tt = tracking_number_folder.createFolderToStoreMeasurements(folder)

    for i in range(n_meas):
        time.sleep(2)
        temp = client.get_node("ns=7;s=MAIN.i_Temperature_Sensor")
        temp_list = temp.get_value()
        temp_vector = np.array(temp_list.get_value())

        fits_file_name = os.path.join(dove, "temperature_%04.fits" % i)
        pyfits.writeto(fits_file_name, temp_vector)

        print("Misura %04d" % i)
    return dove


def analyzer_PT_meas(tt):
    """
    Parameters
    ----------
    tt: string
        tracking number folder
    """
    # tt = '20200911_142702'
    from m4.ground import smooth_function

    folder = config.PT_ROOT_FOLDER
    name = os.path.join(folder, tt)
    list = os.listdir(name)
    list.sort()

    matrix = np.zeros((len(list), OpcUaParameters.num_PT_sensor))
    matrix_s = np.zeros((len(list), OpcUaParameters.num_PT_sensor))

    i = 0
    for t in list:
        hduList = pyfits.open(os.path.join(name, t))
        temp = hduList[0].data
        matrix[i, :] = temp
        i = i + 1

    matrixDiff = matrix - matrix[0, :]
    i = 0
    for i in range(OpcUaParameters.num_PT_sensor):
        ss = smooth_function.smooth(matrixDiff[:, i], 9)
        matrix_s[:, i] = ss
        i = i + 1

    t = np.arange(0, 2 * len(list), 2)
    plt.plot(t, matrix_s / 100)
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [C]")
    plt.title("PT Calibration")

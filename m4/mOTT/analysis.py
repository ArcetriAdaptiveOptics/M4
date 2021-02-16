'''
Authors
  - C. Selmi: written in 2020
'''

import os as _os
import numpy as _np
from astropy.io import fits as _pyfits
from matplotlib import pyplot as _plt
from m4.configuration.config import fold_name as _fold_name


def testCalib():
    '''
    Returns
    -------
    intMat1: numpy array
        cube of tts1's interaction matrix
    intMat2: numpy array
        cube of tts2's interaction matrix
    '''
    tts1 = _np.array(['20201214_091212', '20201214_092528', '20201214_093842',
                    '20201214_095152', '20201214_100508', '20201214_101821',
                    '20201214_103128', '20201214_104441', '20201214_105754',
                    '20201214_111110', '20201214_112435', '20201214_113749'])
    tts2 = _np.array(['20201214_115451', '20201214_120323', '20201214_121200',
                    '20201214_122040', '20201214_122922', '20201214_123807',
                    '20201214_124640', '20201214_125504', '20201214_130327',
                    '20201214_131134', '20201214_131950', '20201214_132822'])
    intMat1 = None
    intMat2 = None
    for i in range(tts1.size):
        mat1 = pippo(tts1[i])
        mat2 = pippo(tts2[i])
        if intMat1 is None:
            intMat1 = mat1
        else:
            intMat1 = _np.stack((intMat1, mat1))
        if intMat2 is None:
            intMat2 = mat2
        else:
            intMat2 = _np.stack((intMat2, mat2))
    return intMat1, intMat2

def _readRepData(tt):
    '''
    Function to read repeatability file fits in tt folder
    '''
    file_name = _os.path.join(_fold_name.REPEATABILITY_ROOT_FOLDER, tt)
    hduList = _pyfits.open(_os.path.join(file_name, 'par.fits'))
    par = hduList[0].data
    hduList = _pyfits.open(_os.path.join(file_name, 'rm.fits'))
    rm = hduList[0].data
    #hduList = pyfits.open(os.path.join(file_name, 'images.fits'))
    #cube = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
    return par, rm

#     def analyzeOptRep(tt):
#         par, rm, cube = ._readRepData(tt)
#         z_list=[]
#         for i in range(cube.shape[2]):
#             masked_ima = cube[:,:,i]
#             coef, mat = zernike.zernikeFit(masked_ima, np.arange(2, 7))
#             z_list.append(coef)
#         return np.array(z_list)

def actsRepeatability(tt):
    '''
    Parameters
    ----------
    tt: string
        tracking number to be analyzed

    Returns
    -------
    pos01_std: numpy array
        vector of actuator's standard deviation from position 1
    pos02_std: numpy array
        vector of actuator's standard deviation from position 2
    pos01_mean: numpy array
        vector of actuator's mean from position 1
    pos02_mean: numpy array
        vector of actuator's mean from position 2
    pos0: numpy array
        vector of actuator's start position
    '''
    par, rm = _readRepData(tt)

    pos01_list_std = []
    pos02_list_std = []
    pos01_list_mean = []
    pos02_list_mean = []
    pos0_list = []
    for i in range(par.shape[2]):
        pos01 = par[:, 0, i] - par[:, 1, i]
        pos02 = par[:, 0, i] - par[:, 2, i]
        pos01_list_std.append(pos01.std())
        pos02_list_std.append(pos02.std())
        pos01_list_mean.append(pos01.mean())
        pos02_list_mean.append(pos02.mean())

        pos0 = par[:,0,i]
        pos0_list.append(pos0.std())

    pos01_std = _np.array(pos01_list_std)
    pos02_std = _np.array(pos02_list_std)
    pos01_mean = _np.array(pos01_list_mean)
    pos02_mean = _np.array(pos02_list_mean)
    pos0 = _np.array(pos0_list)
    return pos01_std, pos02_std, pos01_mean, pos02_mean, pos0

def scanAstigComa(tn):
    '''
    Parameters
    ----------
    tt: string
        tracking number to be analyzed

    Returns
    -------
    zer: numpy array
        vector of zernike
    par_pos: numpy array
        matrix containing parabola position
    rm_pos: numpy array
        matrix containing reference flat position
    '''
    dove = _os.path.join(_fold_name.CALIBRATION_ROOT_FOLDER, tn)
    name = _os.path.join(dove, 'zernike.fits')
    hduList = _pyfits.open(name)
    zer = hduList[0].data
    name = _os.path.join(dove, 'PAR_positions.fits')
    hduList = _pyfits.open(name)
    par_pos = hduList[0].data
    name = _os.path.join(dove, 'RM_positions.fits')
    hduList = _pyfits.open(name)
    rm_pos = hduList[0].data
    _plt.plot(par_pos[0:20, 3], zer[0:20, 4],'o')
    _plt.plot(par_pos[0:20, 3], zer[0:20, 5],'o')
    _plt.xlabel('Par tilt [as]')
    _plt.ylabel('Astigm. Coeff [m]')
    _plt.title(tn)
    _plt.plot(par_pos[0:20, 3], zer[0:20, 6], 'o')
    _plt.plot(par_pos[20:40, 3], zer[20:40, 7],' o')
    _plt.xlabel('Par tilt [as]')
    _plt.ylabel('Coma. Coeff [m]')
    _plt.legend(['X', 'Y'])
    _plt.title(tn)
    return zer, par_pos, rm_pos

def opticalMonitoring():
    pass

def parPistonTest():
    pass

def parTiltTest():
    pass


def alignPlot(tt):
    '''
    Parameters
    ----------
    coeff_matrix: numpy array
        zernike coefficients matrix for all the images
    tt: string
        tracking number of measurements

    Returns
    -------
    figure plot
    '''
    file_name = _os.path.join(_fold_name.REPEATABILITY_ROOT_FOLDER,
                            'Alignment', tt, 'zernike.fits')
    hdu = _pyfits.open(file_name)
    coeff_matrix = hdu[0].data

    x_old = _np.arange(coeff_matrix.shape[0])
    if coeff_matrix.shape[0] == 5:
        x = ['Start_image', 'Perturbed_image', '5_param_alignment',
             'TipTilt_alignment', 'Check_image']
    elif coeff_matrix.shape[0] == 6:
        x = ['Start_image', 'Perturbed_image', 'TipTilt pre-alignment',
             '5_param_alignment', 'TipTilt_alignment', 'Check_image']

    _plt.figure(figsize=(16, 10))
    _plt.subplot(4, 1, 1)
    _plt.plot(x_old, coeff_matrix[:, 1:3], '-o')
    _plt.grid()
    _plt.xticks(x_old, x, rotation=0)
    _plt.ylabel('TipTilt rms [nm]')
    _plt.subplot(4, 1, 2)
    _plt.plot(x_old, coeff_matrix[:, 3], '-ok')
    _plt.grid()
    _plt.xticks(x_old, x, rotation=0)
    _plt.ylabel('Focus rms [nm]')
    _plt.subplot(4, 1, 3)
    _plt.plot(x_old, coeff_matrix[:, 6:8], '-o')
    _plt.grid()
    _plt.xticks(x_old, x,rotation=0)
    _plt.ylabel('Coma rms [nm]')
    _plt.subplot(4, 1, 4)
    _plt.plot(x_old, coeff_matrix[:, 4:6]-coeff_matrix[0,4:6], '-o')
    _plt.grid()
    _plt.xticks(x_old,x, rotation=0)
    _plt.ylabel('Ast rms a.u. [nm]')

    _plt.suptitle(tt + ' Alignment', fontweight='bold', fontsize=20)
    return


###ALTRO###
def pippo(tt):
    '''Function to read interaction matrix in tt folder'''
    from m4.utils.optical_alignment import OpticalAlignment
    al = OpticalAlignment(tt)
    intMat, rec, mask = al._loadAlignmentInfo()
    return intMat

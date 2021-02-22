"""
Authors
  - C. Selmi:  written in 2019

Function for reading file fits::

    from m4.ground import read_data
    numpy_array_data = read_data.readFits_object(fits_file_path)
    or
    numpy_masked_array_data = read_data.readFits_maskedImage(fits_file_path)
    or
    amplitude, mode_vector, cmd_matrix = read_data.readTypeFromFitsName(
                                'ampName.fits', 'mvec.fits', 'cmdMatrix.fits')
"""
import os
from astropy.io import fits as pyfits
import numpy as np
import h5py
from m4.type.modalAmplitude import ModalAmplitude
from m4.type.modalBase import ModalBase
from m4.type.modesVector import ModesVector
from m4.configuration.config import fold_name


def read_phasemap(file_path, ext=0):
    ''' deve leggere sia fits che h5
    '''
    if ext == 0:
        image = readFits_maskedImage(file_path)
    else:
        hf = h5py.File(file_path, 'r')
        hf.keys()
        data1 = hf.get('dataset_1')
        data2 = hf.get('dataset_2')
        image = np.ma.masked_array(np.array(data1),
                                   np.array(data2.astype(bool)))
    return image

### Generiche
def readFits_data(fits_file_path):
    '''
    Parameters
    ----------
    fits_file_path: string
                    fits file path of numpy array to read

    Returns
    -------
            object: numpy array
    '''
    hduList = pyfits.open(fits_file_path)
    obj = hduList[0].data
    return obj

def readFits_maskedImage(fits_file_path):
    '''
    Parameters
    ----------
    fits_file_path: string
                    fits file path of masked array to read

    Returns
    -------
    object: numpy masked array
    '''
    hduList = pyfits.open(fits_file_path)
    ima = hduList[0].data
    immagine = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
    return immagine
###


def readTypeFromFitsName(amplitude_fits_file_name,
                         mode_vector_fits_file_name,
                         cmd_matrix_fits_file_name):
    '''
    Parameters
    ----------
        amplitude_fits_file_name: string
                                 vector with mode amplitude fits file name
        mode_vector_fits_file_name: string
                                    mode or actuator index vector
                                    to be applied fits file name
        cmd_matrix_fits_file_name: string
                                matrix of mode commands fits file name

    Returns
    -------
        amplitude: numpy array
                vector with mode amplitude
        modesVector: numpy array
                    mode or actuator index vector to be applied
        cmd_matrix: numpy array [nActs x nModes]
                    matrix of mode commands
                    diagonal matrix in case of zonal commands
    '''
    ma = ModalAmplitude.loadFromFits(amplitude_fits_file_name)
    amplitude = ma.getModalAmplitude()

    mv = ModesVector.loadFromFits(mode_vector_fits_file_name)
    mode_vector = mv.getModesVector()

    mb = ModalBase.loadFromFits(cmd_matrix_fits_file_name)
    cmd_matrix = mb.getModalBase()
    return amplitude, mode_vector, cmd_matrix


def _readImageFromRunaIFFs(fits_file_path):
    file_name = os.path.join(fold_name.IFFUNCTIONS_ROOT_FOLDER,
                             fits_file_path)
    hduList = pyfits.open(file_name)
    cube = hduList[0].data
    immagine = np.ma.masked_array(cube[0, :, :], mask=np.invert(cube[1, :, :].astype(bool)))
    return immagine

class InterferometerConverter():
    """ Class to use to convert H5 files to masked array"""

    @staticmethod
    def from4D(h5filename):
        """
        Parameters
        ----------
            h5filename: string
                 path of h5 file to convert

        Returns
        -------
                ima: numpy masked array
                     masked array image
        """
        file = h5py.File(h5filename, 'r')
        genraw = file['measurement0']['genraw']['data']
        data = np.array(genraw)
        mask = np.zeros(data.shape, dtype=np.bool)
        mask[np.where(data == data.max())] = True
        ima = np.ma.masked_array(data * 632.8e-9, mask=mask)
        return ima

    @staticmethod
    def fromNew4D(i4dfilename):
        """
        Parameters
        ----------
            h5filename: string
                 path of h5 file to convert

        Returns
        -------
                ima: numpy masked array
                     masked array image
        """
        file = h5py.File(i4dfilename, 'r')
        data = file.get('/Measurement/SurfaceInWaves/Data')
        meas = data[()]
        mask = np.invert(np.isfinite(meas))
        image = np.ma.masked_array(meas * 632.8e-9, mask=mask)
        return image

    @staticmethod
    def fromI4DToSimplerData(i4dname, folder, h5name):
        file = h5py.File(i4dname, 'r')
        data = file.get('/Measurement/SurfaceInWaves/Data')

        file_name = os.path.join(folder, h5name)
        hf = h5py.File(file_name, 'w')
        hf.create_dataset('Data', data=data)
        return file_name

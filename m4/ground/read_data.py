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

def read_phasemap(file_path):
    ''' per leggere i tre formati di dati interferometrici
    '''
    ext = file_path.split('.')[-1]
    if ext == 'fits':
        image = readFits_maskedImage(file_path)
    elif ext=='4D':
        image = InterferometerConverter.fromPhaseCam6110(file_path)
    elif ext=='4Ds':
        image = readFits_maskedImage(file_path)
    elif ext=='h5':
        image = InterferometerConverter.fromPhaseCam4020(file_path)
    return image

def save_phasemap(filename, masked_image, overwrite:bool=False):
    '''
    Function to save data in a standard fashion for OTT usage
    To be implemented: cube saving
    Parameters
    ----------
    filename: string
            full path
    masked_image: array
            masked array of image to be saved

    Returns
    -------
            object: numpy array
    '''
    pyfits.writeto(filename, masked_image.data, overwrite=overwrite)
    pyfits.append(filename, masked_image.mask.astype(np.uint8))

def saveFits_data(filename, data, header=None, overwrite:bool=False):
    """
    Complete function for saving simple fits data
    """
    pyfits.writeto(filename, data, header, overwrite=overwrite)

### Generiche
def readFits_data(fits_file_path):
    """
    Parameters
    ----------
    fits_file_path: str
        Complete filepath of the fits file to load.

    Returns
    -------
    object : ndarray
        The read data.
    """
    hduList = pyfits.open(fits_file_path)
    obj = hduList[0].data
    hduList.close()
    return obj

def readFits_maskedImage(fits_file_path):
    '''
    Parameters
    ----------
    fits_file_path: string
                    fits file path of masked array to read

    Returns
    -------
    masked_array: numpy masked array
    '''
    hduList = pyfits.open(fits_file_path)
    masked_array = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
    hduList.close()
    return masked_array

def readFitsSlimImage(fits_file_path):
    '''
    Parameters
    ----------
    fits_file_path: string
                    fits file path of masked array to read

    Returns
    -------
    immagine: numpy masked array
    '''
    hduList = pyfits.open(fits_file_path)
    data = hduList[0].data
    mask = np.invert(np.isfinite(data))
    immagine = np.ma.masked_array(data, mask=mask)
    hduList.close()
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


def _readImageFromRunaIFFs(file_name):
    #file_name = os.path.join(fold_name.IFFUNCTIONS_ROOT_FOLDER,
    #                        fits_file_path)
    hduList = pyfits.open(file_name)
    cube = hduList[0].data
    immagine = np.ma.masked_array(cube[0, :, :], mask=np.invert(cube[1, :, :].astype(bool)))
    return immagine

class InterferometerConverter():
    """ Class to use to convert H5 files to masked array"""

    @staticmethod
    def fromPhaseCam4020(h5filename):
        """
        Function for PhaseCam4020
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
        mask = np.zeros(data.shape, dtype=bool)
        mask[np.where(data == data.max())] = True
        ima = np.ma.masked_array(data * 632.8e-9, mask=mask)
        return ima

    @staticmethod
    def fromPhaseCam6110(i4dfilename):
        """
        Function for PhaseCam6110
        Parameters
        ----------
            h5filename: string
                 path of h5 file to convert

        Returns
        -------
                ima: numpy masked array
                     masked array image
        """
        with h5py.File(i4dfilename, 'r') as ff:
            data = ff.get('/Measurement/SurfaceInWaves/Data')
            meas = data[()]
            mask = np.invert(np.isfinite(meas))

        image = np.ma.masked_array(meas * 632.8e-9, mask=mask)
        
        return image

    @staticmethod
    def fromFakeInterf(filename):
        """
        Function for fake interferometer
        Parameters
        ----------
            file: string
                 path name for data

        Returns
        -------
                ima: numpy masked array
                     masked array image
        """
        masked_ima = readFits_maskedImage(filename)
        return masked_ima

    @staticmethod
    def fromI4DToSimplerData(i4dname, folder, h5name):
        ''' Function for converting files from 4d 6110 files to H5 files
        Parameters
        ----------
        i4dname: string
            file name path of 4d data
        folder: string
            folder path for new data
        h5name: string
            name for h5 data

        Returns
        -------
        file_name: string
            finale path name
        '''
        file = h5py.File(i4dname, 'r')
        data = file.get('/Measurement/SurfaceInWaves/Data')

        file_name = os.path.join(folder, h5name)
        hf = h5py.File(file_name, 'w')
        hf.create_dataset('Data', data=data)
        return file_name

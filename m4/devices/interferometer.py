'''
Authors
  - C. Selmi: written in 2020
'''
import logging
import os
import time
import shutil
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from m4.configuration import config_folder_names as fold_name
from m4.configuration.ott_parameters import Interferometer
from m4.ground.read_data import InterferometerConverter
from m4.devices.base_interferometer import BaseInterferometer

class I4d4020(BaseInterferometer):
    ''' Class for i4d interferometer

    HOW TO USE IT::

        from m4.devices.interferometer import *
        i4d4020 = I4d4020()
        or
        i4d6110 = I4d6110()
    '''

    def __init__(self):
        """The constructor """
        from oaautils import i4d
        self._ic = InterferometerConverter()
        self._interf = i4d.I4D()
        self._logger = logging.getLogger('4D')

    def acquire_phasemap(self, nframes=1, show=0):
        """
        Parameters
        ----------
            nframes: int
                number of frames
            show: int
                0 to not show the image

        Returns
        -------
            masked_ima: numpy masked array
                    interferometer image
        """
        if nframes == 1:
            masked_ima = self._getMeasurementOnTheFly(self._interf)
        else:
            cube_images = None
            for i in range(nframes):
                ima = self._getMeasurementOnTheFly(self._interf)
                if cube_images is None:
                    cube_images = ima
                else:
                    cube_images = np.ma.dstack((cube_images, ima))
            masked_ima = np.ma.mean(cube_images, axis=2)
        if show != 0:
            plt.clf()
            plt.imshow(masked_ima, origin='lower')
            plt.colorbar()
        return masked_ima

    def save_phasemap(self, location, file_name, masked_image):
        """
        Parameters
        ----------
        location: string
            measurement file path
        file_name: string
            measuremnet fits file name
        masked_image: numpy masked array
            data to save
        """
        fits_file_name = os.path.join(location, file_name)
        pyfits.writeto(fits_file_name, masked_image.data)
        pyfits.append(fits_file_name, masked_image.mask.astype(int))

    def _getMeasurementOnTheFly(self, interf):
        '''
        Parameters
        ----------
            interf: object
                interferometer

        Returns
        -------
            masked_image: numpy masked image
                interferogram
        '''
        filename = '/tmp/prova4d'

        nMeasure = 1
        interf.connect()
        interf.capture(1, name='DM_temp')
        interf.produce('DM_temp')
        interf.disconnect()
        time.sleep(1.0)
        fName = os.path.join(fold_name.PHASECAM_ROOT_FOLDER, 'DM_temp')
        #fName = '/home/m4/4d/Zcopy/DM_temp'

        for i in range(nMeasure):
            shutil.move(fName + '/hdf5/img_%04d.h5' %i,
                        filename + "_m%02d" %i + ".h5")

        shutil.rmtree(fName + '/hdf5')
        shutil.rmtree(fName + '/raw')

        return self._ic.fromPhaseCam4020('/tmp/prova4d_m00.h5')


class I4d6110(BaseInterferometer):
    ''' Class for i4d 6110 interferometer
    '''

    def __init__(self):
        """The constructor """
        from m4.devices.i4d import I4D
        self._i4d = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)
        self._ic = InterferometerConverter()
        self._logger = logging.getLogger('4D')

    def acquire_phasemap(self, nframes=1, delay=0):
        """
        Parameters
        ----------
            nframes: int
                number of frames
            delay: int [s]
                delay between images

        Returns
        -------
            masked_ima: numpy masked array
                    interferometer image
        """
        if nframes == 1:
            width, height, pixel_size_in_microns, data_array = self._i4d.takeSingleMeasurement()
            masked_ima = self._fromDataArrayToMaskedArray(width, height, data_array*632.8e-9)
        else:
            image_list = []
            for i in range(nframes):
                width, height, pixel_size_in_microns, data_array = self._i4d.takeSingleMeasurement()
                masked_ima = self._fromDataArrayToMaskedArray(width, height, data_array*632.8e-9)
                image_list.append(masked_ima)
                time.sleep(delay)
            images = np.ma.dstack(image_list)
            masked_ima = np.ma.mean(images, 2)

#         if show != 0:
#             plt.clf()
#             plt.imshow(masked_ima, origin='lower')
#             plt.colorbar()
        return masked_ima

    def _fromDataArrayToMaskedArray(self, width, height, data_array):
        data = np.reshape(data_array, (width, height))
        idx, idy = np.where(np.isnan(data))
        mask = np.zeros((data.shape[0], data.shape[1]))
        mask[idx, idy] = 1
        masked_ima = np.ma.masked_array(data, mask=mask.astype(bool))
        return masked_ima

    def save_phasemap(self, location, file_name, masked_image):
        """
        Parameters
        ----------
        location: string
            measurement file path
        file_name: string
            measuremnet fits file name
        masked_image: numpy masked array
            data to save
        """
        fits_file_name = os.path.join(location, file_name)
        pyfits.writeto(fits_file_name, masked_image.data)
        pyfits.append(fits_file_name, masked_image.mask.astype(int))

    def burstFramesToSpecificDirectory(self, directory, numberOfFrames):
        '''
        Parameters
        ----------
        directory: string
            directory where to save files
        numberOfFrames: int
            number of frames to acquire
        '''
        self._i4d.burstFramesToSpecificDirectory(directory, numberOfFrames)
        return

    def burstFramesToDisk(self, numberOfFrames):
        '''
        Parameters
        ----------
        numberOfFrames: int
            number of frames to acquire

        '''
        self._i4d.burstFramesToDisk(numberOfFrames)
        return

    def convertRawFrames(self,measurementsDirectory, rawFramesDirectory):
        '''
        Parameters
        ----------
        tn: int
            number of frames to acquire

        '''
        
        self._i4d.convertRawFramesInDirectoryToMeasurementsInDestinationDirectory(measurementsDirectory, rawFramesDirectory)
        return

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
from m4.ground import timestamp
from m4.configuration import config_folder_names as fold_name
from m4.configuration.ott_parameters import Interferometer
from m4.ground.read_data import InterferometerConverter
from m4.devices.base_interferometer import BaseInterferometer
import playsound
from m4.configuration.ott_parameters import Sound
from m4.ground.read_4DConfSettingFile import ConfSettingReader #modRB 20231027 to implement frame2OTTframe here

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
        self._ts = timestamp.Timestamp()

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

    def acquire_detector(self, nframes=1, delay=0):
        """
        Parameters
        ----------
            nframes: int
                number of frames
            delay: int [s]
                delay between images

        Returns
        -------
            data2d: numpy masked array
                    detector interferometer image 
        """
        self.acquire_phasemap()
        
        if nframes == 1:
            data, height, pixel_size_in_microns, width=self._i4d.getFringeAmplitudeData()
            data2d = np.reshape(data, (width, height))
        else:
            image_list = []
            for i in range(nframes):
                data, height, pixel_size_in_microns, width=self._i4d.getFringeAmplitudeData()
                data2d_t = np.reshape(data, (width, height))
                image_list.append(data2d_t)
                time.sleep(delay)
            images = np.ma.dstack(image_list)
            data2d = np.ma.mean(images, 2)

        return data2d

    def _fromDataArrayToMaskedArray(self, width, height, data_array):
       # data = np.reshape(data_array, (width, height))
        data = np.reshape(data_array, (height,width)) #mod20231002, rectangular frames were bad. now fixed

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
        pyfits.append(fits_file_name, masked_image.mask.astype(np.uint8))  #was int, then uint16

    # def burstAndConvertFrom4DPCTom4OTTpc(self, n_frames):
    #     '''
    #     Attention: check if 4d is mounted
    #
    #     Parameters
    #     ----------
    #     n_frames: int
    #          number of frames for burst
    #     '''
    #     fold_capture = 'D:/M4/Capture' #directory where to save files
    #     tn = self._ts.now()
    #     print(tn)
    #     self._i4d.burstFramesToSpecificDirectory(os.path.join(fold_capture, tn+'/'), n_frames)
    #
    #     # convert the frames 
    #     fold_produced ='D:/M4/Produced'
    #     self._i4d.convertRawFramesInDirectoryToMeasurementsInDestinationDirectory(os.path.join(fold_produced, tn), os.path.join(fold_capture, tn))
    #
    #     shutil.move(os.path.join('/home/m4/4d/M4/Produced', tn), fold_name.OPD_IMAGES_ROOT_FOLDER)

    def capture(self, numberOfFrames, folder_name=None):
        '''
        Parameters
        ----------
        numberOfFrames: int
            number of frames to acquire

        Other parameters
        ---------------
        folder_name: string
            if None a tacking number is generate
        
        Returns
        -------
        folder_name: string
            name of folder measurements
        '''
        if folder_name is None:
            folder_name = self._ts.now()
        print(folder_name)
        
        self._i4d.burstFramesToSpecificDirectory(os.path.join(Interferometer.CAPTURE_FOLDER_NAME_4D_PC,
                                                              folder_name), numberOfFrames)
        if Sound.PLAY is True:
            playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH, 'Capture-completed.mp3'))
        return folder_name
    
    def produce(self, folder_name):
        '''
        Parameters
        ----------
        folder_name: string
            name of folder measurements to convert
        '''
        self._i4d.convertRawFramesInDirectoryToMeasurementsInDestinationDirectory(
            os.path.join(Interferometer.PRODUCE_FOLDER_NAME_4D_PC, folder_name),
            os.path.join(Interferometer.CAPTURE_FOLDER_NAME_4D_PC, folder_name))
        
        shutil.move(os.path.join(Interferometer.PRODUCE_FOLDER_NAME_M4OTT_PC, folder_name),
                    fold_name.OPD_IMAGES_ROOT_FOLDER)
        self._rename4D(folder_name)
        if Sound.PLAY is True:
            playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH,'produce-completed.mp3'))
            
    def _rename4D(self, folder):
        """
        Renames the produced 'x.4D' files into '000x.4D'
        """
        fold = os.path.join(fold_name.OPD_IMAGES_ROOT_FOLDER, folder)
        files = os.listdir(fold)
        for file in files:
            if file.endswith('.4D'):
                num_str = file.split('.')[0]
                if num_str.isdigit():
                    num = int(num_str)
                    new_name = f"{num:05d}.4D"
                    old_file = os.path.join(fold, file)
                    new_file = os.path.join(fold, new_name)
                    os.rename(old_file, new_file)

    def getCameraSettings(self):
        '''
        Return 
        ----------
        output: list
        the output is a 4 elements list with width_pixel, height_pixel, offset_x, offset_y, as read from the local copy of the 4D camera settings file 
        '''

        file_path = Interferometer.SETTINGS_CONF_FILE_M4OTT_PC
        setting_reader = ConfSettingReader(file_path)
        width_pixel = setting_reader.getImageWidhtInPixels()
        height_pixel = setting_reader.getImageHeightInPixels()
        offset_x = setting_reader.getOffsetX()
        offset_y = setting_reader.getOffsetY()
        return [width_pixel, height_pixel, offset_x, offset_y]

    def getFrameRate(self):
        '''
        Return 
        ----------
        frame_rate: float
        frame rate of the interferometer
        '''

        file_path = Interferometer.SETTINGS_CONF_FILE_M4OTT_PC
        setting_reader = ConfSettingReader(file_path)
        frame_rate = setting_reader.getFrameRate()
        return frame_rate

    def intoFullFrame(self, img):
        '''
        The function fits the passed frame (expected cropped) into the full interferometer frame (2048x2048), after reading the cropping parameters.

        Parameters
        ----------
        img: masked_array

        Return 
        ----------
        output: masked_array
        the output is the interferometer full frame
        '''

        off = (self.getCameraSettings())[2: 4]
        off = np.flip(off)
        nfullpix = np.array([2048, 2048])
        fullimg = np.zeros(nfullpix)
        fullmask = np.ones(nfullpix)
        offx = off[0]
        offy = off[1]
        sx   = np.shape(img)[0] #croppar[2] 
        sy   = np.shape(img)[1] #croppar[3]
        fullimg[offx: offx+sx, offy: offy + sy]  = img.data
        fullmask[offx: offx+sx, offy: offy + sy] = img.mask
        fullimg = np.ma.masked_array(fullimg, fullmask)
        return fullimg


'''
@author: cs
'''

import os
from scipy import ndimage
import numpy as np
import pyfits
from m4.ground.configuration import Configuration
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.flattening import Flattenig
from m4.analyzer_iffunctions import AnalyzerIFF
from m4.ground import tracking_number_folder
from m4.ground.interferometer_converter import InterferometerConverter

class ZernikeCommand():

    def __init__(self, segment_roi, tt_list_for_an):
        self._roi = segment_roi
        self._ttListForAn = tt_list_for_an
        self._ic = InterferometerConverter()
        self._pupilXYRadius = Configuration.M4_PUPIL_XYRADIUS
        self._zg = ZernikeGenerator(2*self._pupilXYRadius[2])
        self._segmentImaDiameter = 512
        self._measure = None
        self._totalModeCommand = None
        self._anList = []


    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "ZernikeCommandTest")

    def zernikeCommandTest(self, zernike_modes_vector_amplitude):
        '''
        args:
            zernike_modes_vector_amplitude = vector of amplitude of zernike
                                            modes to be test, starting from z=2
        '''
        #ci metto tutto il procedimento
        store_in_folder = self._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        self._dove, self._tt = save._createFolderToStoreMeasurements()

        for tt in self._ttListForAn:
            an = self._createAnalyzer(tt)
            self._anList.append(an)

        self._nModes = zernike_modes_vector_amplitude.shape[0]
        for j in range(self._nModes):
            amp = zernike_modes_vector_amplitude[j]
            zernike_coeff_array = np.zeros(zernike_modes_vector_amplitude.shape[0])
            zernike_coeff_array[j] = amp
            number_of_zernike_mode = j + 2

            total_mode_command = self.singleZernikeCommandTest(zernike_coeff_array,
                                                                   self._anList,
                                                                   number_of_zernike_mode)
            #single_zernike_cube = self._singleZernikeModeCubeMeasurementCreator(number_of_zernike_mode)
            single_zernike_cube = self._TESTFUCTIONFORZERNIKECUBE(number_of_zernike_mode)
            cube_name = 'cube_mode%04d.fits' %number_of_zernike_mode
            self._saveSingleZernikeCube(single_zernike_cube, cube_name)

            total_mode_image = self.imageReconstructor(single_zernike_cube)
            self.differentialPistonMeasurement(total_mode_command)
        return self._tt, total_mode_image


    def imageReconstructor(self, singleZernikeCube):
        #vanno rimesse insieme (rotare e attaccare)
        total_mode_image = 0
        return total_mode_image

    def _readCubes(self):
        cubeList = []
        for i in range(self._nModes):
            j= i+2
            cube_name = 'cube_mode%04d.fits' %j
            fits_file_name = os.path.join(self._dove, cube_name)
            hduList = pyfits.open(fits_file_name)
            cube = np.ma.masked_array(hduList[0].data,
                                      hduList[1].data.astype(bool))
            cubeList.append(cube)
        return cubeList

    def differentialPistonMeasurement(self, total_mode_command):
        #applico questo comando allo specchio (pos e neg) e misuro il pistone
        pass


    def _singleZernikeModeCubeMeasurementCreator(self, number_of_zernike_mode):
        fits_file_path = os.path.join(self._storageFolder(), self._tt)
        singleZernikeCube = None
        for i in range(Configuration.N_SEG):
            pos_file_path = os.path.join(fits_file_path,
                                         'mode%04d_measure_segment%02d_pos.h5' %(number_of_zernike_mode, i))
            positive_image = self._ic.from4D(pos_file_path)
            neg_file_path = os.path.join(fits_file_path,
                                         'mode%04d_measure_segment%02d_neg.h5' %(number_of_zernike_mode, i))
            negative_image = self._ic.from4D(neg_file_path)

            segment_image = positive_image - negative_image / 2 # * amp
            if singleZernikeCube is None:
                singleZernikeCube = segment_image
            else:
                singleZernikeCube = np.ma.dstack((singleZernikeCube, segment_image))

        return singleZernikeCube

    def _TESTFUCTIONFORZERNIKECUBE(self, number_of_zernike_mode):
        fits_file_path = os.path.join(self._storageFolder(), self._tt)
        singleZernikeCube = None
        for i in range(2):
            pos_file_path = os.path.join(fits_file_path,
                                         'mode%04d_measure_segment%02d_pos.fits' %(number_of_zernike_mode, i))
            positive_image = self._readImage(pos_file_path)
            neg_file_path = os.path.join(fits_file_path,
                                         'mode%04d_measure_segment%02d_neg.fits' %(number_of_zernike_mode, i))
            negative_image = self._readImage(neg_file_path)

            seg_ima = positive_image.data - negative_image.data / 2 # * amp
            segment_image = np.ma.masked_array(seg_ima, mask= positive_image.mask)
            if singleZernikeCube is None:
                singleZernikeCube = segment_image
            else:
                singleZernikeCube = np.ma.dstack((singleZernikeCube, segment_image))
        return singleZernikeCube

    def _readImage(self, file_path):
        hduList = pyfits.open(file_path)
        immagine = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
        return immagine

###
    def singleZernikeCommandTest(self, zernike_coeff_array, an_list, number_of_zernike_mode):
        ''' scelto un solo modo di zernike lo applico a tutti i segmenti e
        restituisco il tt delle misure
        '''
        commandsList = []
        self._totalModeCommand = None

        #for i in range(Configuration.N_SEG):
        for i in range(2):
            command_for_segment = self.zernikeCommandForSegment(zernike_coeff_array, i, an_list[i])
            commandsList.append(command_for_segment)
            measure_name = 'mode%04d_measure_segment%02d' %(number_of_zernike_mode, i)
            self.applySegmentCommand(command_for_segment, measure_name, self._dove)

        totalModeCommand = self._totalCommandCreator(commandsList)
        name = 'total_command_mode%04d.fits' %number_of_zernike_mode
        fits_file_name = os.path.join(self._dove, name)
        pyfits.writeto(fits_file_name, totalModeCommand)
        return totalModeCommand

    def _saveSingleZernikeCube(self, single_zernike_cube, cube_name):
        fits_file_path = os.path.join(self._storageFolder(), self._tt)
        fits_file_name = os.path.join(fits_file_path, cube_name)
        pyfits.writeto(fits_file_name, single_zernike_cube.data)
        pyfits.append(fits_file_name, single_zernike_cube.mask.astype(int))

###

    def _totalCommandCreator(self, commandsList):
        for i in range(len(commandsList)):
            if self._totalModeCommand is None:
                self._totalModeCommand = commandsList[i]
            else:
                self._totalModeCommand = np.concatenate((self._totalModeCommand,
                                                     commandsList[i]),
                                                     axis=None)
        return self._totalModeCommand

    def zernikeCommandForSegment(self, zernike_coeff_array, segment_number, an):
        ''' Calcola il comando del modo scelto da dare al segmento
        '''
        surface_map = self._createZernikeModesOnM4(zernike_coeff_array)
        theta, r = self._moveSegmentView(segment_number)
        cropped_image, rotate_image, final_image = self._cropImage(theta, r, surface_map)
        fl = Flattenig(an)
        command_for_segment = - fl.flatCommand(final_image)
        return command_for_segment

    def applySegmentCommand(self, command_for_segment, measure_name, dove):
        ''' Non potendo applicare il comando allo specchio e misurare il wf
        salvo elementi dal cubo di dati in /IFFunctions/20170216_123645
        '''
        #applica al segmento il comando positivo e negativo
        # dm.setShape(command_for_segment)
        self._measure = self._TESTIMAGETOSAVE(17)
        name = measure_name + '_pos.fits'
        self._saveMeasurement(dove, name)
        # dm.setShape(-command_for_segment)
        self._measure = self._TESTIMAGETOSAVE(36)
        name = measure_name + '_neg.fits'
        self._saveMeasurement(dove, name)
        return

    def _createAnalyzer(self, tt):
        ''' Per ora crea an con il cubo di 25 misure che ho usato per i test
        '''
        #file_name = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
        #fits_file_name = os.path.join(file_name, 'Cube.fits')
        an = AnalyzerIFF.loadCubeFromIFFMeasureToCreateAn(tt)
        return an

    def _createZernikeModesOnM4(self, zernike_coeff_array):
        ''' crea una superficie circolare con il diametro in pixel dello specchio
            ho scleto 1024x1024
        '''
        zernike_surface_map = 0.0
        first_zern_mode_index = 2
        last_zern_mode_index = 2 + len(zernike_coeff_array)
        index_zernike_modes = np.arange(first_zern_mode_index, last_zern_mode_index)
        zd = self._zg.getZernikeDict(index_zernike_modes)

        for i in index_zernike_modes:
            zernike_surface_map = zernike_surface_map + \
                                    zernike_coeff_array[i-2] * zd[i]

        return zernike_surface_map

    def _cropImage(self, theta, r, zernike_surface_map):
        ''' Taglia l'immagine dello zernike grande quanto m4 in un'imagine
        512x512 centrata nel centro del semento scelto
        '''
        centrex = np.int(self._pupilXYRadius[0] + r * np.cos(theta))
        centrey = np.int(self._pupilXYRadius[1] + r * np.sin(theta))
        segment_radius = self._segmentImaDiameter // 2
        cropped_image = zernike_surface_map[centrey - segment_radius:
                                            centrey + segment_radius,
                                            centrex - segment_radius:
                                            centrex + segment_radius]
        rotate_image = ndimage.rotate(cropped_image, -theta - 90, reshape=False)
        final_image = np.ma.masked_array(rotate_image.data, mask= self._roi)

        return cropped_image, rotate_image, final_image

    def _moveSegmentView(self, segment_number):
        ''' Dovrà comandare il movimento della parabola.
            Per ora restituisce l'angolo e il raggio del centrodel segmento scelto
            '''
        tetha = Configuration.REFERENCE_ANGLE * segment_number
        r = Configuration.SEGMENT_DISTANCE_FROM_CENTRE
        return tetha, r


    def _TESTIMAGETOSAVE(self, n_cube_image):
        fold = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, '20170216_123645')
        cube_path = os.path.join(fold, 'CubeMeasure.fits')
        hduList = pyfits.open(cube_path)
        cube = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        test_image = cube[:,:,n_cube_image]
        return test_image

    def _saveMeasurement(self, dove, name):
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, self._measure.data)
        pyfits.append(fits_file_name, self._measure.mask.astype(int))
        
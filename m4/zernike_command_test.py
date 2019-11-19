'''
@author: cs
'''

from scipy import ndimage
import os
import numpy as np
import pyfits
from m4.ground.configuration import Configuration
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.flattening import Flattenig
from m4.analyzer_iffunctions import AnalyzerIFF
from m4.ground import tracking_number_folder
from m4.ground import object_from_fits_file_name as obj

class ZernikeCommand():

    def __init__(self, segment_roi):
        self._roi = segment_roi
        self._pupilXYRadius = Configuration.M4_PUPIL_XYRADIUS
        self._zg = ZernikeGenerator(2*self._pupilXYRadius[2])
        self._segmentImaDiameter = 512
        self._measure = None
        self._ttList = None
        self._commandsList = None
        self._finalSegmentImage = None
        self._totalCommand = None

    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "ZernikeCommandTest")

    def zernikeCommandTest(self):
        #ci metto tutto il procedimento
        pass


    def imageReconstructor(self):
        final_segments_images_list = self._segmentsImagesFromMeasurement()
        #vanno rimesse insieme (rotare e attaccare)
        pass
        
    def differentialPistonMeasurement(self):
        totalCommand = self._totalCommandCreator()
        #applico questo comando allo specchio (pos e neg) e misuro il pistone
        pass

    def _segmentsImagesFromMeasurement(self):
        for tt in self._ttList:
            fits_file_path = os.path.join(self._storageFolder(), self._ttList)
            pos_file_path = os.path.join(fits_file_path, 'Measure_pos')
            positive_image = obj.readImageFromFitsFileName(pos_file_path)
            neg_file_path = os.path.join(fits_file_path, 'Measure_neg')
            negative_image = obj.readImageFromFitsFileName(neg_file_path)

            segment_image = positive_image - negative_image / 2 # * amp
            self._finalSegmentImage.appen(segment_image)
        return self._finalSegmentImage 

    def _totalCommandCreator(self):
        for i in range(len(self._commandsList)):
            if self._totalCommand is None:
                self._totalCommand = self._commandsList[i]
            else:
                self._totalCommand = np.concatenate((self._totalCommand,
                                                     self._commandsList[i]),
                                                     axis=None)
        return self._totalCommand



###
    def singleZernikeCommandTest(self, zernike_coeff_array, tt_list_for_an):
        ''' scelto un solo modo di zernike lo applico a tutti i segmento e
        restituisco i tt delle misure
        '''
        for i in range(Configuration.N_SEG):
            an = self._createAnalyzer(tt_list_for_an[i])
            command_for_segment = self.zernikeCommandForSegment(zernike_coeff_array, i, an)
            self._commandsList.append(command_for_segment)
            measure_name = 'Measure'
            tt = self.applySegmentCommand(command_for_segment, measure_name)
            self._ttList.append(tt)
            return self._ttList

###


    def zernikeCommandForSegment(self, zernike_coeff_array, segment_number, an):
        ''' Calcola il comando del modo scelto da dare al segmento
        '''
        surface_map = self._createZernikeModesOnM4(zernike_coeff_array)
        theta, r = self._moveSegmentView(segment_number)
        cropped_image, rotate_image, final_image = self._cropImage(theta, r, surface_map)
        fl = Flattenig(an)
        command_for_segment = - fl.flatCommand(final_image)
        return command_for_segment

    def applySegmentCommand(self, command_for_segment, measure_name):
        ''' Non potendo applicare il comando allo specchio e misure il wf
        salvo il comando stesso
        '''
        #applica al segmento il comando positivo e negativo
        store_in_folder = self._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        dove, tt = save._createFolderToStoreMeasurements()
        # dm.setShape(command_for_segment)
        self._measure = command_for_segment
        name = measure_name + '_pos.fits'
        self._saveMeasurement(dove, name)
        # dm.setShape(-command_for_segment)
        self._measure = - command_for_segment
        name = measure_name + '_neg.fits'
        self._saveMeasurement(dove, name)
        return tt

    def _createAnalyzer(self, tt):
        ''' Per ora crea an con il cubo di 25 misure che ho usato per i test
        '''
        file_name = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
        #fits_file_name = os.path.join(file_name, 'Cube.fits')
        an = AnalyzerIFF.loadTestMeasureFromFits(tt)
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

    def _saveMeasurement(self, dove, name):
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, self._measure)
        
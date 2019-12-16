'''
@author: cs
'''

import os
import logging
from scipy import ndimage
import numpy as np
from astropy.io import fits as pyfits
from m4.ground.configuration import Configuration
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.flattening import Flattenig
from m4.analyzer_iffunctions import AnalyzerIFF
from m4.ground import tracking_number_folder
from m4.ground.interferometer_converter import InterferometerConverter

class ZernikeCommand():
    ''' Class created to control Zernike modes at all m4
        (6 segments at the same time)
    '''

    def __init__(self, segment_roi, tt_list_for_an):
        """The constructor:
            segment_roi: central roi for the segment
            tt_list_for_an: list of tracking number to use for an creation
        """
        self._logger = logging.getLogger('ZER_CMD_TEST:')
        self._roi = segment_roi
        self._ttListForAn = tt_list_for_an
        self._ic = InterferometerConverter()
        self._pupilXYRadius = Configuration.M4_MECHANICAL_PUPIL_XYRADIUS
        self._zg = ZernikeGenerator(2*self._pupilXYRadius[2])
        self._diameterInPixelForSegmentImages = Configuration.DIAMETER_IN_PIXEL_FOR_SEGMENT_IMAGES
        self._bigDiameter = Configuration.BIG_IMAGE_DIAMETER
        self._measure = None
        self._totalModeCommand = None
        self._anList = []
        self._zernikeSurfaceCube = None
        self._m4ImagesCube = None


    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "ZernikeCommandTest")

    def _saveMeasurementInfo(self):
        fits_file_name = os.path.join(self._dove, 'info.fits')
        header = pyfits.Header()
        header['ANTTSEG0'] = self._ttListForAn[0]
        header['ANTTSEG1'] = self._ttListForAn[1]
        header['ANTTSEG2'] = self._ttListForAn[2]
        header['ANTTSEG3'] = self._ttListForAn[3]
        header['ANTTSEG4'] = self._ttListForAn[4]
        header['ANTTSEG5'] = self._ttListForAn[5]
        pyfits.writeto(fits_file_name, self._ampVector, header)
        pyfits.append(fits_file_name, self._roi.astype(int))

    def zernikeCommandTest(self, zernike_modes_vector_amplitude):
        '''
        args:
            zernike_modes_vector_amplitude = vector of amplitude of zernike
                                            modes to be test, starting from z=2

        returns:
                tt = tracking number for the measurement
                zernikeSurfaceCube = all zernike surface generated on M4 [512, 512, n_modes]
                m4ImagesCube = cube consisting of the 6 images of the segments [512, 512, n_modes]
        '''
        self._ampVector = zernike_modes_vector_amplitude
        store_in_folder = self._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        self._dove, self._tt = save._createFolderToStoreMeasurements()

        self._saveMeasurementInfo()
        self._logger.debug('Info file creation')

        self._logger.info('Creation of an list')
        for tt in self._ttListForAn:
            an = self._createAnalyzer(tt)
            self._anList.append(an)

        self._nModes = zernike_modes_vector_amplitude.shape[0]
        for j in range(self._nModes):
            amp = zernike_modes_vector_amplitude[j]
            if amp == 0:
                pass
            else:
                zernike_coeff_array = np.zeros(zernike_modes_vector_amplitude.shape[0])
                zernike_coeff_array[j] = amp
                number_of_zernike_mode = j + 2

                surface_map, total_mode_command = self.singleZernikeCommandTest(zernike_coeff_array,
                                                                                self._anList,
                                                                                number_of_zernike_mode)
                if self._zernikeSurfaceCube is None:
                    self._zernikeSurfaceCube = surface_map
                else:
                    self._zernikeSurfaceCube = np.dstack((self._zernikeSurfaceCube, surface_map))
                #single_zernike_cube = self._singleZernikeModeCubeMeasurementCreator(number_of_zernike_mode)
                single_zernike_cube = self._TESTFUCTIONFORZERNIKECUBE(number_of_zernike_mode, amp)
                cube_name = 'cube_mode%04d.fits' %number_of_zernike_mode
                self._saveSingleZernikeCube(single_zernike_cube, cube_name)

                total_mode_image = self.imageReconstructor(single_zernike_cube)
                piston = self._differentialPistonMeasurement(total_mode_command)
                final_total_mode_image = total_mode_image - piston
                if self._m4ImagesCube is None:
                    self._m4ImagesCube = final_total_mode_image
                else:
                    self._m4ImagesCube = np.ma.dstack((self._m4ImagesCube, final_total_mode_image))


        fits_file_name = os.path.join(self._dove, 'surfCube.fits')
        pyfits.writeto(fits_file_name, self._zernikeSurfaceCube)
        fits_file_name = os.path.join(self._dove, 'm4ImageCube.fits')
        pyfits.writeto(fits_file_name, self._m4ImagesCube.data)
        pyfits.append(fits_file_name, self._m4ImagesCube.mask.astype(int))
        return self._tt, self._zernikeSurfaceCube, self._m4ImagesCube

    def readSurfM4ImageCubes(self, tt):
        dove = os.path.join(self._storageFolder(), tt)

        fits_file_name = os.path.join(dove, 'surfCube.fits')
        hduList = pyfits.open(fits_file_name)
        self._zernikeSurfaceCube = hduList[0].data
        fits_file_name = os.path.join(dove, 'm4ImageCube.fits')
        hduList = pyfits.open(fits_file_name)
        self._m4ImagesCube = hduList[0].data
        return self._zernikeSurfaceCube, self._m4ImagesCube

    def analyzerResults(self, tt):
        surf_cube, m4_images_cube = self.readSurfM4ImageCubes(tt)
        diff_list = []
        rms_list = []
        a = m4_images_cube[0][2].shape
        for i in range(a[0]):
            diff = m4_images_cube[: , :, i] - surf_cube[:, :, i]
            rms = diff.std()
            diff_list.append(diff)
            rms_list.append(rms)
        rms = np.array(rms_list)
        return diff_list, rms


    def imageReconstructor(self, singleZernikeCube):
        """
        args:
            singleZernikeCube = cube of one zernike mode
                                (shape [pixels, pixels, number of segments measured])

        returns:
                final_total_mode_image = reconstructed image joining the six segments
                imaList = image list of the six segments
        """
        #vanno rimesse insieme (rotare e attaccare))
        final_total_mode_image = np.zeros((self._bigDiameter, self._bigDiameter))
        final_total_mask = np.ones((self._bigDiameter, self._bigDiameter))
        imaList = []
        maskList = []
        for i in range(Configuration.N_SEG):
            total_mode_image = np.zeros((self._bigDiameter, self._bigDiameter))
            total_image_masks = np.ones((self._bigDiameter, self._bigDiameter))
            seg_img = singleZernikeCube[:, :, i]
            seg_mask = singleZernikeCube[:, :, i].mask.astype(int)
            theta = Configuration.REFERENCE_ANGLE_RAD * i
            theta_degrees = Configuration.REFERENCE_ANGLE_DEGREES * i
            r = Configuration.SEGMENT_DISTANCE_FROM_CENTRE
            seg_img_rot = ndimage.rotate(seg_img, - theta_degrees + 90, reshape=False)
            seg_mask_rot = ndimage.rotate(seg_mask, - theta_degrees + 90, reshape=False, cval=1)

            centerx = np.int(self._bigDiameter/2 + r * np.cos(theta))
            centery = np.int(self._bigDiameter/2 + r * np.sin(theta))
            segment_radius = self._diameterInPixelForSegmentImages // 2
            total_mode_image[centery - segment_radius:
                             centery + segment_radius,
                             centerx - segment_radius:
                             centerx + segment_radius] = seg_img_rot

            total_image_masks[centery - segment_radius:
                             centery + segment_radius,
                             centerx - segment_radius:
                             centerx + segment_radius] = seg_mask_rot

            imaList.append(total_mode_image)
            final_total_mode_image = final_total_mode_image + total_mode_image
            final_total_mask = final_total_mask + total_image_masks
            maskList.append(total_image_masks)

        mask = np.zeros(final_total_mask.shape, dtype=np.bool)
        mask[np.where(final_total_mask == final_total_mask.max())] = 1

        true_final_image = np.ma.masked_array(final_total_mode_image, mask=mask)

        return true_final_image

    def readCubes(self, tt=None):
        """Reads all the cubes present in the last generated tracking number
        and puts them in a list.

        returns:
                cubeList = list of cubes [n_mode][pixels, pixels, segment]
        """
        if tt is None:
            dove = self._dove
        else:
            dove = os.path.join(self._storageFolder(), tt)

        cubeList = []
        list = os.listdir(dove)
        list.sort()
        cube_name = list[0:2]
        for i in range(len(cube_name)):
            #j= i+2
            #cube_name = 'cube_mode%04d.fits' %j
            fits_file_name = os.path.join(dove, cube_name[i])
            hduList = pyfits.open(fits_file_name)
            cube = np.ma.masked_array(hduList[0].data,
                                      hduList[1].data.astype(bool))
            cubeList.append(cube)
        return cubeList

    def _differentialPistonMeasurement(self, total_mode_command):
        """
        args:
            total_mode_command = np.array(N_ACTS_TOT)

        returns:
                piston = piston value measured with SPL
        """
        #applico questo comando allo specchio (pos e neg) e misuro il pistone
        piston = 0
        return piston


    def _singleZernikeModeCubeMeasurementCreator(self, number_of_zernike_mode, amp):
        """
        args:
            number_of_zernike_mode = int number of zernike mode
            amp = int number of the mode amplitude applied

        returns:
                singleZernikeCube = cube of one zernike mode
                                    (shape [pixels, pixels, number of segments measured])
        """
        self._logger.info('Single cube creation: for mode %s', number_of_zernike_mode)
        fits_file_path = os.path.join(self._storageFolder(), self._tt)
        singleZernikeCube = None
        for i in range(Configuration.N_SEG):
            pos_file_path = os.path.join(fits_file_path,
                                         'mode%04d_measure_segment%02d_pos.h5' %(number_of_zernike_mode, i))
            positive_image = self._ic.from4D(pos_file_path)
            neg_file_path = os.path.join(fits_file_path,
                                         'mode%04d_measure_segment%02d_neg.h5' %(number_of_zernike_mode, i))
            negative_image = self._ic.from4D(neg_file_path)

            seg_ima = positive_image.data - negative_image.data / (2 * amp)
            segment_image = np.ma.masked_array(seg_ima, mask=positive_image.mask)
            if singleZernikeCube is None:
                singleZernikeCube = segment_image
            else:
                singleZernikeCube = np.ma.dstack((singleZernikeCube, segment_image))

        return singleZernikeCube

    def _TESTFUCTIONFORZERNIKECUBE(self, number_of_zernike_mode, amp):
        fits_file_path = os.path.join(self._storageFolder(), self._tt)
        singleZernikeCube = None
        #for i in range(2):
        for i in range(Configuration.N_SEG):
            pos_file_path = os.path.join(fits_file_path,
                                         'mode%04d_measure_segment%02d_pos.fits' %(number_of_zernike_mode, i))
            positive_image = self._readImage(pos_file_path)
            neg_file_path = os.path.join(fits_file_path,
                                         'mode%04d_measure_segment%02d_neg.fits' %(number_of_zernike_mode, i))
            negative_image = self._readImage(neg_file_path)

            seg_ima = positive_image.data - negative_image.data / (2 * amp)
            segment_image = np.ma.masked_array(seg_ima, mask=positive_image.mask)
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
        ''' Chosen only one mode of Zernike I apply it to all segments and
        I return the total mode command from the measurement.

        args:
            zernike_coeff_array = array containing the amplitude of the mode
                            to create located in the its corresponding position
                            exemple: np.array([z2, z3, z4....])
            an_list = list of the 6 object an
            number_of_zernike_mode = (int) zernike mode to measure

        returns:
                surface_map = zernike mode surface map
                totalModeCommand = array of 5253 elements
        '''
        self._logger.info('Total command calculation for zernike mode %s', number_of_zernike_mode)
        commandsList = []
        self._totalModeCommand = None
        surface_map = np.zeros((self._bigDiameter, self._bigDiameter))
        surface_map_on_m4 = self._createZernikeModesOnM4(zernike_coeff_array)
        center_big_image = int(self._bigDiameter/2)
        raggio_zer = self._pupilXYRadius[2]
        surface_map[center_big_image - raggio_zer:center_big_image + raggio_zer,
                    center_big_image - raggio_zer:center_big_image + raggio_zer] = surface_map_on_m4

        for i in range(Configuration.N_SEG):
        #for i in range(2):
            self._logger.debug('Single segment command calculation: segment number %s', i)
            command_for_segment = self.zernikeCommandForSegment(surface_map, i, an_list[i])
            commandsList.append(command_for_segment)
            measure_name = 'mode%04d_measure_segment%02d' %(number_of_zernike_mode, i)
            self._applySegmentCommand(command_for_segment, an_list[i], measure_name, self._dove)

        totalModeCommand = self._totalCommandCreator(commandsList)
        name = 'total_command_mode%04d.fits' %number_of_zernike_mode
        fits_file_name = os.path.join(self._dove, name)
        pyfits.writeto(fits_file_name, totalModeCommand)
        return surface_map, totalModeCommand

    def _saveSingleZernikeCube(self, single_zernike_cube, cube_name):
        """ Function to save the cube of a single zernike mode. It is formed
        with the 6 mode measurement on the  segments.

        args:
            single_zernike_cube = cube of 6 masked arrays elements
            cube_name = fits file name for the cube
        """
        fits_file_path = os.path.join(self._storageFolder(), self._tt)
        fits_file_name = os.path.join(fits_file_path, cube_name)
        pyfits.writeto(fits_file_name, single_zernike_cube.data)
        pyfits.append(fits_file_name, single_zernike_cube.mask.astype(int))

###

    def _totalCommandCreator(self, commandsList):
        """ Use the command list to create the total command for the mirror.
        args:
            commandsList = list of 6 elements

        returns:
            totalModeCommand = array of 5253 elements
        """
        for i in range(len(commandsList)):
            if self._totalModeCommand is None:
                self._totalModeCommand = commandsList[i]
            else:
                self._totalModeCommand = np.concatenate((self._totalModeCommand,
                                                         commandsList[i]),
                                                         axis=None)
        return self._totalModeCommand

    def zernikeCommandForSegment(self, surface_map, segment_number, an):
        ''' Calculate the command of the chosen mode to give to the segment.
        args:
            surface_map = zernike mode surface map on big image
            segment_number = number of the chosen segment
            an = object analyzer

        returns:
                command_for_segment = vector of 892 element containing the
                                segment's command
        '''
        theta, theta_degrees, r = self._moveSegmentView(segment_number)
        cropped_image, rotate_image, final_image = \
                        self._cropImage(theta, theta_degrees, r, surface_map)
        fl = Flattenig(an)
        command_for_segment = - fl.flatCommand(final_image)
        return command_for_segment

    def _applySegmentCommand(self, command_for_segment, an, measure_name, dove):
        ''' Non potendo applicare il comando allo specchio e misurare il wf
        salvo il wf sintetico da me ricostruito

        args:
            command_for_segment = vector of 892 element containing the
                                segment's command
            an = analyzer object of the segment
            measure_name = fits file name of the measure
            dove = fits file path to save
        '''
        #applica al segmento il comando positivo e negativo
        # dm.setShape(command_for_segment)
        # passo in argomento le ampiezze invece del comando così
        #nella funzione test creo il wf sitetico
        self._measure = self._TESTIMAGESAVEFROMCMD(command_for_segment, an)
        name = measure_name + '_pos.fits'
        self._saveMeasurement(dove, name)
        # dm.setShape(-command_for_segment)
        self._measure = self._TESTIMAGESAVEFROMCMD(-command_for_segment, an)
        name = measure_name + '_neg.fits'
        self._saveMeasurement(dove, name)
        return

    def _createAnalyzer(self, tt):
        ''' Usavo il cubo di 25 misure che ho usato per i test
            Poi usando 892 oggetti modi in tt = 20170216_123645
            Ora gli 811 modi in tt = 20170630_105105
        '''
        #file_name = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
        #fits_file_name = os.path.join(file_name, 'Cube.fits')
        an = AnalyzerIFF.loadCubeFromIFFMeasureToCreateAn(tt)
        return an

    def _createZernikeModesOnM4(self, zernike_coeff_array):
        ''' Creates the Zernike mode on a circular surface with the diameter,
            in pixels, of the mirror
            (I have chosen 1024x1024)

            args:
            zernike_coeff_array = array containing the amplitude of the mode
                to create located in the its corresponding position
                exemple: np.array([z2, z3, z4....])

            returns:
                    zernike_surface_map_on_m4 = surface map for the zernike mode
        '''
        zernike_surface_map_on_m4 = 0.0
        first_zern_mode_index = 2
        last_zern_mode_index = 2 + len(zernike_coeff_array)
        index_zernike_modes = np.arange(first_zern_mode_index, last_zern_mode_index)
        zd = self._zg.getZernikeDict(index_zernike_modes)

        for i in index_zernike_modes:
            zernike_surface_map_on_m4 = zernike_surface_map_on_m4 + \
                                    zernike_coeff_array[i-2] * zd[i]

        return zernike_surface_map_on_m4

    def _cropImage(self, theta, theta_degrees, r, zernike_surface_map):
        ''' Cut the image of Zernike, as big as m4 pixels dimension, into an image
        512x512 centered in the center of the chosen segment

        args:
            theta = theta angle in radians of the segment
            theta_degrees = theta angle in degrees of the segment
            r = radial distance of segment's center
            zernike_surface_map = 1236x1236 surface map for the zernike mode

        returns:
            cropped_image = image cut around the center of the segment
            rotate_image = image rotated to obtain the standard view of
                            the segment, corresponding to the mask used
                            to create the object
            final_image = image with its mask
        '''
        centerx = np.int(self._bigDiameter/2 + r * np.cos(theta))
        centery = np.int(self._bigDiameter/2 + r * np.sin(theta))
        segment_radius = self._diameterInPixelForSegmentImages // 2
        cropped_image = zernike_surface_map[centery - segment_radius:
                                            centery + segment_radius,
                                            centerx - segment_radius:
                                            centerx + segment_radius]
        rotate_image = ndimage.rotate(cropped_image, theta_degrees - 90, reshape=False)
        final_image = np.ma.masked_array(rotate_image.data, mask=self._roi)

        return cropped_image, rotate_image, final_image

    def _moveSegmentView(self, segment_number):
        ''' Dovrà comandare il movimento della parabola.
            Per ora restituisce l'angolo e il raggio del centro del segmento scelto
            '''
        theta = Configuration.REFERENCE_ANGLE_RAD * segment_number
        theta_degrees = Configuration.REFERENCE_ANGLE_DEGREES * segment_number
        r = Configuration.SEGMENT_DISTANCE_FROM_CENTRE
        return theta, theta_degrees, r


    def _TESTIMAGETOSAVEFROMCUBE(self, n_cube_image):
        """ Test functions: use cube measurements saved in fold
            to create test image to save

            args:
                n_cube_image = number of cube image to use

            returns:
                    test_image = image from the cube
        """
        fold = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, '20170216_123645')
        cube_path = os.path.join(fold, 'CubeMeasure.fits')
        hduList = pyfits.open(cube_path)
        cube = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        test_image = cube[:, :, n_cube_image]
        return test_image

    def _TESTIMAGESAVEFROMCMD(self, segment_command, an):
        """ Test functions: create the image using the synthetic wf
        """
        fl = Flattenig(an)
        wf = fl.syntheticWfCreator(self._roi, segment_command)
        return wf

    def _saveMeasurement(self, dove, name):
        """ Save measurement data.
            args:
                dove = path where to save the measurement file
                name = measurement fits file name
        """
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, self._measure.data)
        pyfits.append(fits_file_name, self._measure.mask.astype(int))
        
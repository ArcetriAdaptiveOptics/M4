'''
Authors
  - G. Pariani, R.Briguglio: written in 2016
  - C. Selmi: ported to Python in 2020
  - A. Puglisi and C. Selmi: python code debugging in 2021
'''

import os
import numpy as np
import glob
import logging
from astropy.io import fits as pyfits
from photutils import centroids
import matplotlib.pyplot as plt
from m4.ground import smooth_function as sf
from m4.configuration.ott_parameters import OttParameters
from m4.configuration import config_folder_names as fold_name


class SplAnalyzer():
    '''
    Class used to analyze images, whit SPL camera, at different wavelengths.

    HOW TO USE IT::
        from m4.analyzers.SPL_data_analyzer import SplAnalyzer
        an = SplAnalyzer()
        piston, piston_smooth = an.analyzer(tt)
    '''

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('SPL_AN:')
        self.tn_fringes = OttParameters.TN_FRINGES
        self._Qm = None
        self._QmSmooth = None
        self._matrix = None
        self._matrixSmooth = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        #return Path(__file__).parent/'data'
        #return fold_name.SPL_ROOT_FOLDER
        lift_path = '/home/labot/LIFT/SPL/data'
        return lift_path

    def analyzer(self, tt):
        '''
        Analyze measurement data and compare it with synthetic data.

        Parameters
        ----------
        tt: string
            tracking number of the measurement data
        Returns
        -------
        piston: int
                piston value
        piston_smooth: int
                piston value after smoothing data
        '''
        dove = os.path.join(self._storageFolder(), tt)

        self._logger.info('Analysis of tt = %s', tt)

        lambda_path = os.path.join(dove, 'lambda_vector.fits')
        hduList = pyfits.open(lambda_path)
        lambda_vector = hduList[0].data

        cube, cube_normalized = self._readCameraFrames(tt)
        img = np.sum(cube_normalized, 2)
        pick = self._newThr(img)
        matrix = np.zeros((pick[3]-pick[2] + 1, lambda_vector.shape[0])) # 150 + 1 pixel
        matrix_smooth = np.zeros((pick[3]-pick[2] + 1, lambda_vector.shape[0]))
        crop_frame_cube = None
        for i in range(lambda_vector.shape[0]):
            frame = cube[:, :, i]
            crop_frame = frame[pick[0]:pick[1], pick[2]:pick[3] + 1]
            #import code
            #code.interact(local=dict(globals(), **locals()))
            if np.max(crop_frame) > 4000:
                print("**************** WARNING: saturation detected!")
            if crop_frame_cube is None:
                crop_frame_cube = crop_frame
            else:
                crop_frame_cube = np.dstack((crop_frame_cube, crop_frame))

            y = np.sum(crop_frame, 0)
            area = np.sum(y[:])
            y_norm = y / area
            if i == 0:
                mm = 1.2 * np.max(y_norm)
            matrix[:, i] = y_norm

            w = sf.smooth(y_norm, 4)
            w = w[:pick[3]-pick[2] + 1]
            matrix_smooth[:, i] = w

        matrix[np.where(matrix == np.nan)] = 0
        self._matrix = matrix
        self._matrixSmooth = matrix_smooth
        #self._saveMatrix(matrix)
        piston, piston_smooth = self._templateComparison(matrix,
                                                         matrix_smooth,
                                                         lambda_vector)

        self._savePistonResult(tt, piston, piston_smooth)

        return piston, piston_smooth
    ### PLOT MATRIX ###
    # x = lambda_vector
    # y = np.arange(151)* spl._pix2um
    # imshow(spl._matrix, extent = [x[0],x[19],y[0], y[150]], origin= 'lower');
    #colorbar(); plt.xlabel('lambda [nm]'); plt.ylabel('position [um]');
    #plt.title('TN = %s' %tt)

    def _saveMatrix(self, matrix, tt):
        destination_file_path = self._storageFolder()
        fits_file_name = os.path.join(destination_file_path,
                                      'fringe_result.fits', tt)
        pyfits.writeto(fits_file_name, matrix)


    def _readCameraFrames(self, tt):
        ''' Read images in a specific tracking number and return the cube
        args:
            tt = tracking number
        returns:
            cube = [pixels, pixels, n_frames=lambda]
            cube_normalized = creating with normalized images,
                                [pixels, pixels, n_frames=lambda]
        '''
        dove = os.path.join(self._storageFolder(), tt)

        path_list = []
        for f in  glob.iglob(os.path.join(dove, 'image_*.fits')):
            path_list.append(f)

        path_list.sort()
        cube = []
        cube_normalized = []
        for i in range(len(path_list)):
            #print('Reading ',path_list[i])
            hduList = pyfits.open(path_list[i])
            image = hduList[0].data
            image_normalized = image / np.sum(image)
            cube.append(image)
            cube_normalized.append(image_normalized)
        cube = np.dstack(cube)
        cube_normalized = np.dstack(cube_normalized)
        return cube, cube_normalized

    def _newThr(self, img):
        ''' Calculate the peak position of the image '''
        #histogram thr to be discussed
#        counts, bin_edges = np.histogram(img)
#        edges = (bin_edges[2:] + bin_edges[1:len(bin_edges)-1]) / 2
#        thr = 5 * edges[np.where(counts == max(counts))]
#        idx = np.where(img < thr)
#        img[idx] = 0
        # cy, cx = scin.measurements.center_of_mass(img)
        cx, cy = centroids.centroid_2dg(img)
        baricenterCoord = [np.int(round(cy)), np.int(round(cx))]
        pick = [baricenterCoord[0]-25, baricenterCoord[0]+25,
                baricenterCoord[1]-75, baricenterCoord[1]+75]
        return pick

    def _templateComparison(self, matrix, matrix_smooth, lambda_vector):
        '''
        Compare the matrix obtained from the measurements with
        the one recreated with the synthetic data in tn_fringes.
        Parameters
        ----------
        matrix: [pixels, lambda]
                measured matrix
        matrix_smooth: [pixels, lambda]
                        measured matrix smooth
        lambda_vector: numpy array
                        vector of wavelengths (between 400/700 nm)
        Returns
        -------
        piston: int
                piston value
        '''
        #tn_fringes = '20181026'
        self._logger.debug('Template Comparison with data in tn_fringes = %s',
                           self.tn_fringes)
        dove = os.path.join(self._storageFolder(), 'Fringes',
                            self.tn_fringes)
        dove_delta = os.path.join(dove, 'Differential_piston.fits')
        hduList = pyfits.open(dove_delta)
        delta = hduList[0].data
        dove_lambda_synth = os.path.join(dove, 'Lambda.fits')
        hduList = pyfits.open(dove_lambda_synth)
        lambda_synth_from_data = hduList[0].data * 1e9
        lambda_synth = self._myLambaSynth(lambda_synth_from_data)

        idx = np.isin(lambda_synth, lambda_vector)

        Qm = matrix - np.mean(matrix)
        self._Qm = Qm
        Qm_smooth = matrix_smooth - np.mean(matrix_smooth)
        self._QmSmooth = Qm_smooth
        #creare la matrice di Giorgio della giusta dimensione
        F = []
        for i in range(1, delta.shape[0]):
            file_name = os.path.join(dove, 'Fringe_%05d.fits' %i)
            #print('Reading ',file_name)
            hduList = pyfits.open(file_name)
            fringe = hduList[0].data
            fringe_selected = fringe[:, idx]
            F.append(fringe_selected)
        F = np.dstack(F)
        Qt = F - np.mean(F)
        self._Qt = Qt

        R = np.zeros(delta.shape[0]-1)
        R_smooth = np.zeros(delta.shape[0]-1)
        for i in range(delta.shape[0]-1):
            # print('Processing',i)
            R[i] = np.sum(Qm[:, :]*Qt[:, :, i]) / (np.sum(Qm[:, :]**2)**5
                                                   * np.sum(Qt[:, :, i]**2)**5)
            R_smooth[i] = np.sum(Qm_smooth[:, :]*Qt[:, :, i]) / \
                        (np.sum(Qm_smooth[:, :]**2)**5 * np.sum(Qt[:, :, i]**2)**5)

        plt.plot(R)
        plt.show()
        idp = np.where(R == max(R))
        idp_smooth = np.where(R_smooth == max(R_smooth))
        piston = delta[idp]
        piston_smooth = delta[idp_smooth]
        return piston, piston_smooth

    def _myLambaSynth(self, lambda_synth_from_data):
        ''' Transforms its values into integers
        '''
        my_lambda = np.zeros(lambda_synth_from_data.shape[0])
        for j in range(lambda_synth_from_data.shape[0]):
            bb = np.int(round(lambda_synth_from_data[j]))
            my_lambda[j] = bb
        return my_lambda

    def _savePistonResult(self, tt, piston, piston_smooth):
        dove = os.path.join(self._storageFolder(), tt,
                            'piston_result_prova.txt')
        file = open(dove, 'w+')
        file.write('%4e, %4e' %(piston[0], piston_smooth[0]))
        file.close()

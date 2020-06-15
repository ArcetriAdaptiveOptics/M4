'''
Autors
  - G. Pariani, R.Briguglio: written in 2016
  - C. Selmi: ported to Python in 2020
'''

import os
import logging
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.utils import image_extender
from m4.utils.img_redux import TipTiltDetrend
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.ground.configuration import Configuration
from m4.ground.interferometer_converter import InterferometerConverter

class Caliball():

    def __init__(self):
        """The constructor """
        self._maskthreshold = 5
        self._rmsthreshold = 3
        self._logger = logging.getLogger('CALIBALL:')
        self._ic = InterferometerConverter()
        self._zOnM4 = ZernikeOnM4()
        self._ttd = TipTiltDetrend()


    @staticmethod
    def _storageFolder():
        """ Creates the path for measurement data"""
        return os.path.join(Configuration.OPD_DATA_FOLDER,
                            "Caliball")

    def doStat(self):
        cube = self._readCube()
        #cube_ttr = self._readCube(1)
        rs_img = self._readRsImg()
        mask_point = np.zeros(cube.shape[2])
        for i in range(mask_point.shape[0]):
            mask_point[i] = np.sum(np.invert(cube[:,:,i].mask))
        mask_ord = np.sort(mask_point)
        aa = np.where(mask_point == min(mask_point))[0][0]
        bad_dataset = cube[:, :, aa]
        bad_mask = cube[:, :, aa].mask

        plt.plot(np.arange(mask_point.shape[0]), mask_ord, 'o'); plt.xscale('log')
        plt.ylabel('# valid points'); plt.xlabel('# frames')
        plt.title('Cumulative plot of mask valid points')

        rs = image_extender.imageExtender(rs_img)
        coef, mat = self._zOnM4.zernikeFit(rs, np.array([2, 3]))
        rs_ttr = self._ttd.ttRemoverFromCoeff(coef, rs)
        r0 = rs_ttr.std()

        cube_ttr = self.readCubeTTR()
        rs_vect = np.zeros(cube_ttr.shape[0])
        for j in range(cube_ttr.shape[0]):
            rs_vect[j] = cube_ttr[j,:,:].std()

        plt.figure()
        plt.plot(np.arange(cube_ttr.shape[0]), rs_vect, label='Data'); plt.yscale('log')
        plt.plot(np.zeros(cube_ttr.shape[0]) + r0, label='Average')
        plt.ylabel('m RMS'); plt.xlabel('# frames')
        plt.title('Images WFE'); plt.legend()
        return bad_dataset, bad_mask
    #plot(x, mask_ord, 'o'); pyplot.xscale('log'); plt.ylabel('# valid points'); plt.xlabel('# frames'); plt.title('Cumulative plot of mask valid points')

    def pixel_std(self):
        cube_ttr = self._readCube(1)
        std_image = cube_ttr.std(axis=2)
        mean_image = cube_ttr.mean(axis=2)

        plt.imshow(std_image, origin='lower'); plt.colorbar()
        plt.title('Pixel StDev')
        plt.figure()
        plt.imshow(mean_image, origin='lower'); plt.colorbar()
        plt.title('RS average value')
        return std_image, mean_image

    def doselaverage(self):
        cube_ttr = self._readCube(1)
        mask_point = np.zeros(cube_ttr.shape[2])
        for i in range(mask_point.shape[0]):
            mask_point[i] = np.sum(cube_ttr[:,:,i].mask)
        idx = np.where(mask_point <= (min(mask_point) + self._maskthreshold))

        rs_std = np.zeros(idx[0].shape[0])
        for i in range(idx[0].shape[0]):
            rs_std[i] = np.std(cube_ttr[:,:,idx[0][i]])

        rthresh = np.mean(rs_std)+self._rmsthreshold*np.std(rs_std)
        idr = np.where(rs_std < rthresh)
        r_img = np.mean(cube_ttr[:,:,idr], axis=3)
        std_img = np.std(cube_ttr[:,:,idr], axis=3)
        r_image = r_img[:,:,0]
        std_image = std_img[:,:,0]

        plt.imshow(r_image, origin='lower'); plt.colorbar()
        plt.title('RS measurement error')
        plt.figure()
        plt.imshow(std_image, origin='lower'); plt.colorbar()
        plt.title('Pixel StDev, thresh = %d' %self._rmsthreshold)
        return r_image, std_image
    # non deve starci idr ma gli idx[idr]

###
    def createImgCubeFile(self):
        path = Caliball._storageFolder()
        total_cube = None
        for j in range(1,4):
            fold = os.path.join(path, 'test%d' %j)
            cube = self._createMeasurementCube(fold)
            if total_cube is None:
                total_cube = cube
            else:
                total_cube = np.ma.dstack((total_cube, cube))
        self._saveCube(total_cube, 'Total_Cube.fits')
        return total_cube

    def createCubeTTrFromCube(self):
        # ci mette un'eternità a fare l'estenzione dell'immagine
        cube = self._readCube()
        cube_ttr = None
        for i in range(cube.shape[2]):
            image = image_extender.imageExtender(cube[:,:,i])
            coef, mat = self._zOnM4.zernikeFit(image, np.array([2, 3]))
            image_ttr = self._ttd.ttRemoverFromCoeff(coef, image)
            if cube_ttr is None:
                cube_ttr = image_ttr
            else:
                cube_ttr = np.ma.dstack((cube_ttr, image_ttr))
        self._saveCube(cube_ttr, 'Total_Cube_ttr.fits')
        return cube_ttr

    def createCubeTTrFromRuna(self):
        file_name = '/Users/rm/imgcubefit.fits'
        hduList = pyfits.open(file_name)
        cube_runa = hduList[0].data #(3000, 500, 496)
        cube = self._readCube()
        cube_ttr = None
        for i in range(cube_runa.shape[0]):
            image = np.ma.masked_array(cube_runa[i, :, :], mask=cube.mask[:,:,i])
            if cube_ttr is None:
                cube_ttr = image
            else:
                cube_ttr = np.ma.dstack((cube_ttr, image))
        return cube_ttr

    def createRsImgFile(self):
        cube = self._readCube()
        mask_point = np.zeros(cube.shape[2])
        for i in range(mask_point.shape[0]):
            mask_point[i] = np.sum(cube[:,:,i].mask)
        idx = np.where(mask_point <= (min(mask_point) + self._maskthreshold))
        tot = idx[0].shape[0]
        image = np.sum(cube[:,:,idx], 3) / tot
        mask = np.prod(cube[:,:,idx].mask, 3)
        rs_img = np.ma.masked_array(image[:,:,0], mask=mask)

        fits_file_name = os.path.join(Caliball._storageFolder(), 'rs_img.fits')
        pyfits.writeto(fits_file_name, rs_img.data)
        pyfits.append(fits_file_name, rs_img.mask.astype(int))
        return rs_img

    def _createMeasurementCube(self, fold):
        cube = None
        list = os.listdir(fold)
        for i in range(len(list)-1):
            name = 'img_%04d.h5' %i
            file_name = os.path.join(fold, name)
            ima = self._ic.from4D(file_name)
            if cube is None:
                cube = ima
            else:
                cube = np.ma.dstack((cube, ima))
        return cube

    def _saveCube(self, total_cube, name):
        fits_file_name = os.path.join(Caliball._storageFolder(), name)
        pyfits.writeto(fits_file_name, total_cube.data)
        pyfits.append(fits_file_name, total_cube.mask.astype(int))

    def _readCube(self, ttr=None):
        if ttr is None:
            file_name = os.path.join(Caliball._storageFolder(), 'Total_Cube.fits')
        else:
            #file_name = os.path.join(Caliball._storageFolder(), 'Total_Cube_ttr_runa.fits')
            file_name = os.path.join(Caliball._storageFolder(), 'Total_Cube_ttr.fits')
        hduList = pyfits.open(file_name)
        cube = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        return cube

    def readCubeTTR(self):
        file_name = '/Users/rm/imgcubefit.fits'
        #file_name = '/mnt/cargo/data/M4/OTT-Review/CaliBallTest/imgcubefit.fits'
        hduList = pyfits.open(file_name)
        cube_ttr = hduList[0].data #(3000, 500, 496)
        return cube_ttr

    def _readRsImg(self):
        file_name = os.path.join(Caliball._storageFolder(), 'rs_img.fits')
        #file_name = '/mnt/cargo/data/M4/OTT-Review/CaliBallTest/RS-img.fits'
        hduList = pyfits.open(file_name)
        rs_img = np.ma.masked_array(hduList[0].data,
                                    hduList[1].data.astype(bool))
        return rs_img

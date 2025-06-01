"""
Authors
  - G. Pariani, R.Briguglio: written in 2016
  - C. Selmi: ported to Python in 2020
"""

import os
import logging
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.utils.image_reducer import TipTiltDetrend
from m4.ground import zernike
from m4.configuration import config_folder_names as fold_name
from m4.ground.read_data import InterferometerConverter
from m4.misc import tip_tilt_interf_fit


class Caliball:
    """
    Class for data analysis.

    HOW TO USE IT::

        from m4.caliball_average import Caliball
        folder_name = 'tt_rotBall_num'
        cal = Caliball(folder_name)
        cal.createDataForAnalysis()
        r_image, std_image = cal.dataSelectionAndAnalysis()

    """

    def __init__(self, folder_name):
        """The constructor"""
        self._folderName = folder_name
        self._maskthreshold = 1000  # 5
        self._rmsthreshold = 3
        self._logger = logging.getLogger("CALIBALL:")
        self._ic = InterferometerConverter()
        self._ttd = TipTiltDetrend()

        self._cube = None
        self._cube_ttr = None

    @staticmethod
    def _storageFolder():
        """Creates the path for measurement data"""
        return fold_name.CALIBALL_ROOT_FOLDER

    def createDataForAnalysis(self):
        """Create file fits for data analysis

        Returns
        ----------
        self._folderName: string
                        folder name for data
        """
        self._cube = self._createMeasurementCube()
        self._cube_ttr = self._createCubeTTrFromCube()
        rs_ima = self._createRsImgFile()
        return self._folderName

    def validMaskPoint(self):
        """Create the cumulative plot of mask valid points
        (valid points within each frames)

        Returns
        ----------
        bad_dataset: numpy array
                    worst measurement data
        bad_mask: numpy array
                    worst mask
        """
        cube = self._readCube()
        # cube_ttr = self._readCube(1)
        mask_point = np.zeros(cube.shape[2])
        for i in range(mask_point.shape[0]):
            mask_point[i] = np.sum(np.invert(cube[:, :, i].mask))
        mask_ord = np.sort(mask_point)
        aa = np.where(mask_point == min(mask_point))[0][0]
        bad_dataset = cube[:, :, aa]
        bad_mask = cube[:, :, aa].mask

        plt.plot(np.arange(mask_point.shape[0]), mask_ord, "o")
        plt.xscale("log")
        plt.ylabel("# valid points")
        plt.xlabel("# frames")
        plt.title("Cumulative plot of mask valid points")
        return bad_dataset, bad_mask

    def validFrames(self):
        """Create the plot of individual RMS of each frames
        (without tip/tilt) together with RMS of the average frame
        """
        rs_img = self._readRsImg()
        coef, mat = zernike.zernikeFit(rs_img, np.array([2, 3]))
        surf = zernike.zernikeSurface(rs_img, coef, mat)
        rs_ttr = rs_img - surf
        r0 = rs_ttr.std()

        cube_ttr = self._readCube(1)
        rs_vect = np.zeros(cube_ttr.shape[2])
        for j in range(cube_ttr.shape[2]):
            rs_vect[j] = cube_ttr[:, :, j].std()

        plt.figure()
        plt.plot(np.arange(cube_ttr.shape[2]), rs_vect, label="Data")
        plt.yscale("log")
        plt.plot(np.zeros(cube_ttr.shape[2]) + r0, label="Average")
        plt.ylabel("m RMS")
        plt.xlabel("# frames")
        plt.title("Images WFE")
        plt.legend()

    # plot(x, mask_ord, 'o'); pyplot.xscale('log'); plt.ylabel('# valid points');
    # plt.xlabel('# frames'); plt.title('Cumulative plot of mask valid points')

    def pixel_std(self):
        """
        Create the plot of single pixel stdev and mean along the sequence

        Returns
        -------
        std_image: masked array
                single pixel stdev along the sequence
        mean_image: masked array
                mean along the sequence
        """
        cube_ttr = self._readCube(1)
        std_image = cube_ttr.std(axis=2)
        mean_image = cube_ttr.mean(axis=2)

        plt.imshow(std_image, origin="lower")
        plt.colorbar()
        plt.title("Pixel StDev")
        plt.figure()
        plt.imshow(mean_image, origin="lower")
        plt.colorbar()
        plt.title("RS average value")
        return std_image, mean_image

    def dataSelectionAndAnalysis(self):
        """
        Create the plot of pixel stdev and mean after
        rejecting the frames whit bad mask

        Returns
        -------
        std_image: masked array
                stdev along the sequence
        mean_image: masked array
                mean along the sequence
        """
        cube_ttr = self._readCube(1)
        mask_point = np.zeros(cube_ttr.shape[2])
        for i in range(mask_point.shape[0]):
            mask_point[i] = np.sum(np.invert(cube_ttr[:, :, i].mask))
        idx = np.where(mask_point >= (max(mask_point) - self._maskthreshold))

        rs_std = np.zeros(idx[0].shape[0])
        for i in range(idx[0].shape[0]):
            rs_std[i] = np.std(cube_ttr[:, :, idx[0][i]])

        rthresh = np.mean(rs_std) + self._rmsthreshold * np.std(rs_std)
        idr = np.where(rs_std < rthresh)

        idxr = idr
        for i in range(idr[0].shape[0]):
            idxr[0][i] = idx[0][idr[0][i]]
        r_img = np.mean(cube_ttr[:, :, idxr], axis=2)
        std_img = np.std(cube_ttr[:, :, idxr], axis=2)
        mask = np.prod(cube_ttr[:, :, idxr].mask, 2)
        r_image = np.ma.masked_array(r_img[:, :, 0], mask=mask[:, :, 0])
        std_image = np.ma.masked_array(std_img[:, :, 0], mask=mask[:, :, 0])

        plt.imshow(r_image, origin="lower")
        plt.colorbar()
        plt.title("RS measurement error")
        plt.figure()
        plt.imshow(std_image, origin="lower")
        plt.colorbar()
        plt.title("Pixel StDev, thresh = %d" % self._rmsthreshold)
        return r_image, std_image

    # non deve starci idr ma gli idx[idr]

    ###

    def _createCubeTTrFromCube(self, fitEx=None):
        # ci mette un eternit a fare l estenzione dell immagine
        # cube = self._readCube()
        cube_ttr = None
        if fitEx is None:
            for i in range(self._cube.shape[2]):
                image = self._cube[:, :, i]
                coef, mat = zernike.zernikeFit(image, np.array([2, 3]))
                surf = zernike.zernikeSurface(image, coef, mat)
                image_ttr = image - surf
                if cube_ttr is None:
                    cube_ttr = image_ttr
                else:
                    cube_ttr = np.ma.dstack((cube_ttr, image_ttr))
        else:
            for i in range(self._cube.shape[2]):
                coef, interf_coef = tip_tilt_interf_fit.fit(self._cube[:, :, i])
                image_ttr = self._ttd.ttRemoverFromCoeff(coef, self._cube[:, :, i])
                if cube_ttr is None:
                    cube_ttr = image_ttr
                else:
                    cube_ttr = np.ma.dstack((cube_ttr, image_ttr))

        self._saveCube(cube_ttr, "Total_Cube_ttr.fits")
        return cube_ttr

    def _createRsImgFile(self):
        cube = self._readCube()
        mask_point = np.zeros(cube.shape[2])
        for i in range(mask_point.shape[0]):
            mask_point[i] = np.sum(cube[:, :, i].mask)
        idx = np.where(mask_point <= (min(mask_point) + self._maskthreshold))
        tot = idx[0].shape[0]
        image = np.sum(cube[:, :, idx], 3) / tot
        mask = np.prod(cube[:, :, idx].mask, 3)
        rs_img = np.ma.masked_array(image[:, :, 0], mask=mask)

        fits_file_name = os.path.join(
            Caliball._storageFolder(), self._folderName, "rs_img.fits"
        )
        pyfits.writeto(fits_file_name, rs_img.data)
        pyfits.append(fits_file_name, rs_img.mask.astype(int))
        return rs_img

    def _createMeasurementCube(self):
        cube = None
        fold = os.path.join(Caliball._storageFolder(), self._folderName, "hdf5")
        list = os.listdir(fold)
        for i in range(len(list) - 1):
            name = "img_%04d.h5" % i
            file_name = os.path.join(fold, name)
            ima = self._ic.fromPhaseCam4020(file_name)
            if cube is None:
                cube = ima
            else:
                cube = np.ma.dstack((cube, ima))
        self._saveCube(cube, "Total_Cube.fits")
        return cube

    def _saveCube(self, total_cube, name):
        fits_file_name = os.path.join(Caliball._storageFolder(), self._folderName, name)
        pyfits.writeto(fits_file_name, total_cube.data)
        pyfits.append(fits_file_name, total_cube.mask.astype(int))

    def _readCube(self, ttr=None):
        if ttr is None:
            file_name = os.path.join(
                Caliball._storageFolder(), self._folderName, "Total_Cube.fits"
            )
        else:
            # file_name = os.path.join(Caliball._storageFolder(), 'Total_Cube_ttr_runa.fits')
            file_name = os.path.join(
                Caliball._storageFolder(), self._folderName, "Total_Cube_ttr.fits"
            )
        hduList = pyfits.open(file_name)
        cube = np.ma.masked_array(hduList[0].data, hduList[1].data.astype(bool))
        return cube

    def _readRsImg(self):
        file_name = os.path.join(
            Caliball._storageFolder(), self._folderName, "rs_img.fits"
        )
        # file_name = '/mnt/cargo/data/M4/OTT-Review/CaliBallTest/RS-img.fits'
        hduList = pyfits.open(file_name)
        rs_img = np.ma.masked_array(hduList[0].data, hduList[1].data.astype(bool))
        return rs_img

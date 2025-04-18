"""
Authors
  - C. Selmi: written in 2020
"""

import os
import glob
import numpy as np
import logging
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration import config_folder_names as fold_name
from m4.configuration.ott_parameters import Interferometer
from m4.ground import read_data
from m4.ground import zernike


class Analysis:
    """
    HOW TO USE IT::

    from m4.mini_OTT.analysis import Analysis
    an = Analysis(tt)
    """

    def __init__(self, tt):
        """The constructor"""
        self._logger = logging.getLogger("mOTT_ANALYIS:")
        self.tt = tt

    @staticmethod
    def _storageFolder():
        """Return the path where to save measurement data"""
        return fold_name.OUT_FOLDER()

    # Results = self._storageFolder()

    def _readRepData(self):
        """
        Function to read repeatability file fits in tt folder
        """
        file_name = os.path.join(fold_name.REPEATABILITY_ROOT_FOLDER, self.tt)
        hduList = pyfits.open(os.path.join(file_name, "par.fits"))
        par = hduList[0].data
        hduList = pyfits.open(os.path.join(file_name, "rm.fits"))
        rm = hduList[0].data
        # hduList = pyfits.open(os.path.join(file_name, 'images.fits'))
        # cube = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
        return par, rm

    def actsRepeatability(self):
        """
        Parameters
        ----------
        tt: string
            tracking number to be analyzed

        Returns
        -------
        pos01_std: numpy array
            vector of actuator's standard deviation from position 1
        pos02_std: numpy array
            vector of actuator's standard deviation from position 2
        pos01_mean: numpy array
            vector of actuator's mean from position 1
        pos02_mean: numpy array
            vector of actuator's mean from position 2
        pos0: numpy array
            vector of actuator's start position
        """
        par, rm = self._readRepData()

        pos01_list_std = []
        pos02_list_std = []
        pos01_list_mean = []
        pos02_list_mean = []
        pos0_list = []
        for i in range(par.shape[2]):
            pos01 = par[:, 0, i] - par[:, 1, i]
            pos02 = par[:, 0, i] - par[:, 2, i]
            pos01_list_std.append(pos01.std())
            pos02_list_std.append(pos02.std())
            pos01_list_mean.append(pos01.mean())
            pos02_list_mean.append(pos02.mean())

            pos0 = par[:, 0, i]
            pos0_list.append(pos0.std())

        pos01_std = np.array(pos01_list_std)
        pos02_std = np.array(pos02_list_std)
        pos01_mean = np.array(pos01_list_mean)
        pos02_mean = np.array(pos02_list_mean)
        pos0 = np.array(pos0_list)
        return pos01_std, pos02_std, pos01_mean, pos02_mean, pos0

    def scanAstigComa(self):
        """
        Parameters
        ----------
        tt: string
            tracking number to be analyzed

        Returns
        -------
        zer: numpy array
            vector of zernike
        par_pos: numpy array
            matrix containing parabola position
        rm_pos: numpy array
            matrix containing reference flat position
        """
        dove = os.path.join(fold_name.CALIBRATION_ROOT_FOLDER, self.tt)
        name = os.path.join(dove, "zernike.fits")
        hduList = pyfits.open(name)
        zer = hduList[0].data
        name = os.path.join(dove, "PAR_positions.fits")
        hduList = pyfits.open(name)
        par_pos = hduList[0].data
        name = os.path.join(dove, "RM_positions.fits")
        hduList = pyfits.open(name)
        rm_pos = hduList[0].data
        plt.plot(par_pos[0:20, 3], zer[0:20, 4], "o")
        plt.plot(par_pos[0:20, 3], zer[0:20, 5], "o")
        plt.xlabel("Par tilt [as]")
        plt.ylabel("Astigm. Coeff [m]")
        plt.title(self.tt)
        plt.plot(par_pos[0:20, 3], zer[0:20, 6], "o")
        plt.plot(par_pos[20:40, 3], zer[20:40, 7], " o")
        plt.xlabel("Par tilt [as]")
        plt.ylabel("Coma. Coeff [m]")
        plt.legend(["X", "Y"])
        plt.title(self.tt)
        return zer, par_pos, rm_pos

    def opticalMonitoring(self):
        pass

    def parPistonTest(self):
        pass

    def parTiltTest(self):
        pass

    def alignPlot(self):
        """
        Parameters
        ----------
        coeff_matrix: numpy array
            zernike coefficients matrix for all the images
        tt: string
            tracking number of measurements

        Returns
        -------
        figure plot
        """
        file_name = os.path.join(
            fold_name.REPEATABILITY_ROOT_FOLDER, "Alignment", self.tt, "zernike.fits"
        )
        hdu = pyfits.open(file_name)
        coeff_matrix = hdu[0].data

        x_old = np.arange(coeff_matrix.shape[0])
        if coeff_matrix.shape[0] == 5:
            x = [
                "Start_image",
                "Perturbed_image",
                "5_param_alignment",
                "TipTilt_alignment",
                "Check_image",
            ]
        elif coeff_matrix.shape[0] == 6:
            x = [
                "Start_image",
                "Perturbed_image",
                "TipTilt pre-alignment",
                "5_param_alignment",
                "TipTilt_alignment",
                "Check_image",
            ]

        plt.figure(figsize=(16, 10))
        plt.subplot(4, 1, 1)
        plt.plot(x_old, coeff_matrix[:, 1:3], "-o")
        plt.grid()
        plt.xticks(x_old, x, rotation=0)
        plt.ylabel("TipTilt rms [nm]")
        plt.subplot(4, 1, 2)
        plt.plot(x_old, coeff_matrix[:, 3], "-ok")
        plt.grid()
        plt.xticks(x_old, x, rotation=0)
        plt.ylabel("Focus rms [nm]")
        plt.subplot(4, 1, 3)
        plt.plot(x_old, coeff_matrix[:, 6:8], "-o")
        plt.grid()
        plt.xticks(x_old, x, rotation=0)
        plt.ylabel("Coma rms [nm]")
        plt.subplot(4, 1, 4)
        plt.plot(x_old, coeff_matrix[:, 4:6] - coeff_matrix[0, 4:6], "-o")
        plt.grid()
        plt.xticks(x_old, x, rotation=0)
        plt.ylabel("Ast rms a.u. [nm]")

        plt.suptitle(self.tt + " Alignment", fontweight="bold", fontsize=20)
        return

    def _convert(self, seconds):
        mmin, sec = divmod(seconds, 60)
        hour, mmin = divmod(mmin, 60)
        #   return "%d:%02d:%02d" % (hour, mmin, sec)
        return "%d:%02d" % (hour, mmin)

    def longTerm_analysis(self, zc1, zc2=None, ntick=6, figure=True):
        where = fold_name.OPD_SERIES_ROOT_FOLDER
        #'/mnt/m4storage/Data/M4Data/OPTData/OPD_series/'
        path = os.path.join(where, self.tt)
        D = sorted(glob.glob(os.path.join(path, self.tt[0:-6]) + "*"))
        zern = pyfits.open(os.path.join(path, "zernike.fits"))[0].data
        temp = pyfits.open(os.path.join(path, "temperature.fits"))[0].data
        t0 = str(D[0][1 + D[0].rfind("_") : D[0].find(".")])
        hs = float(t0[0:2]) * 3600
        ms = float(t0[2:4]) * 60
        s = float(t0[4::])
        t0s = hs + ms + s
        time = zern[:, 0]
        times = time + t0s
        ttime = np.copy(times)
        timeh = ttime / 3600
        taxis = []
        for ii in range(len(ttime)):
            if ttime[ii] > 24 * 3600:
                ttime[ii] -= 24 * 3600
            taxis.append(self._convert(ttime[ii]))
        PtV = np.max(zern[:, 4] * 10**9) - np.min(zern[:, 4] * 10**9)
        sstd = np.std(zern[:, 4] * 10**9)

        plt.figure()
        start = 0
        stop = len(zern[:, 0])  # - 2000
        plt.plot(
            timeh[start:stop], zern[start:stop, zc1] * 10**9, "k", label="Z" + str(zc1)
        )
        if zc2 is not None:
            plt.plot(
                timeh[start:stop],
                zern[start:stop, zc2] * 10**9,
                "r",
                label="Z" + str(zc2),
            )
        #    plt.plot(timeh[start:stop],zern[start:stop,3]*10**9,'k',label='Z2')
        #    plt.plot(timeh[start:stop],zern[start:stop,2]*10**9,'r',label='Z3')
        if figure is True:
            plt.xticks(
                timeh, taxis[:: len(taxis) / ntick], rotation=0
            )  # ,fontsize=fsize)
            plt.locator_params(nbins=ntick, axis="x", tight=None)
            plt.xlabel("Time [hrs:min]")
        else:
            plt.xlabel("Time [h]")
        plt.ylabel("rms [nm]")
        plt.title(self.tt)
        plt.legend()
        plt.grid()

        results_path = os.path.join(fold_name.OUT_FOLDER, "LongTermStability")
        dove = os.path.join(results_path, self.tt)
        if os.path.exists(dove):
            dove = dove
        else:
            os.makedirs(dove)

        name = os.path.join(dove, "%s-tiptiltstability.png" % self.tt)
        if os.path.isfile(name):
            os.remove(name)
        plt.savefig(name)

    #    return timeh, zern, taxis, PtV, sstd, temp

    def longTerm_rmsConvection(self):
        where = fold_name.OPD_SERIES_ROOT_FOLDER
        # where = '/mnt/m4storage/Data/M4Data/OPTData/OPD_series/'
        # where = '/home/labot/data/M4/Data/M4Data/OPTData/OPD_series'
        path = os.path.join(where, self.tt)
        D1 = sorted(glob.glob(os.path.join(path, "*.fits")))
        # D = D1[0:-int(len(D1)/8)]
        D = D1[0:-2]

        cube = None
        for name in D:
            print(name)
            image = read_data.readFits_maskedImage(name)
            if cube is None:
                cube = image
            else:
                cube = np.ma.dstack((cube, image))
        mean = np.ma.mean(cube, axis=2)

        rms_list = []
        for i in range(cube.shape[2]):
            # print(i)
            ima = cube[:, :, i] - mean
            coef, mat = zernike.zernikeFit(ima, np.arange(8) + 1)
            surf = zernike.zernikeSurface(ima, coef, mat)
            new_ima = ima - surf
            rms = new_ima.std()
            rms_list.append(rms)

        rms = np.array(rms_list)
        x = np.arange(rms.size)
        x_list = []
        for i in range(rms.size):
            aa = D[i].split("_")[-1]
            tt = aa.split(".")[0]
            x_list.append(tt)
        x_time = np.array(x_list)
        ntick = 11
        plt.figure(figsize=(10, 6))
        plt.plot(x, rms, "-")
        plt.xticks(x, x_time[:: len(x_time) / ntick], rotation=45)
        plt.locator_params(nbins=ntick, axis="x", tight=True)
        plt.ylabel("rms[m]")
        plt.title("%s" % self.tt)

        results_path = os.path.join(fold_name.OUT_FOLDER, "LongTermStability")
        dove = os.path.join(results_path, self.tt)
        if os.path.exists(dove):
            dove = dove
        else:
            os.makedirs(dove)
        name = os.path.join(dove, "rmsMeanDiff.fits")
        pyfits.writeto(name, rms, overwrite=True)

        name = os.path.join(dove, "%s-rms.png" % self.tt)
        if os.path.isfile(name):
            os.remove(name)
        plt.savefig(name)
        return rms

    def tipTiltVSTime(self, data_file_path, start, stop):
        """
        Attenzione: questa funzione vuole il file path
        """
        last_name = data_file_path.split("/")[-1]
        if last_name == "hdf5":
            list_tot = glob.glob(os.path.join(data_file_path, "*.h5"))
            tt = data_file_path.split("/")[-2]
            # ext = 1
        else:
            list_tot = glob.glob(os.path.join(data_file_path, "*.fits"))
            tt = data_file_path.split("/")[-1]
            # ext = 0

        list_tot.sort()
        list = list_tot[start:stop]

        if last_name == "hdf5":
            time = np.arange(0, len(list)) * 1 / Interferometer.BURST_FREQ
        else:
            time = "calcolare dal tt"

        tip_tilt = []
        for name in list:
            image = read_data.read_phasemap(name)
            coef, mat = zernike.zernikeFit(image, np.array([1, 2, 3]))
            tip_tilt.append(coef)

        tip = np.array(tip_tilt)[:, 1]
        tilt = np.array(tip_tilt)[:, 2]
        # PLOT
        plt.figure(figsize=(10, 6))
        plt.plot(time, tip * 1e9, label="TIP")
        plt.plot(time, tilt * 1e9, label="TILT")
        plt.legend()
        plt.title(tt)
        plt.ylabel("zernike coef [nm]")
        plt.xlabel("time [s]")

        # SAVE
        results_path = os.path.join(fold_name.OUT_FOLDER, "Noise")
        x = data_file_path.split("/")
        if last_name == "hdf5":
            dove = os.path.join(results_path, x[len(x) - 2])
        else:
            dove = os.path.join(results_path, x[len(x) - 1])
        if os.path.exists(dove):
            dove = dove
        else:
            os.makedirs(dove)
        name = os.path.join(dove, "tipTilt_time.png")
        if os.path.isfile(name):
            os.remove(name)
        plt.savefig(name)

        return np.array(tip_tilt), time


# ###ALTRO###
# def pippo(tt):
#     '''Function to read interaction matrix in tt folder'''
#     from m4.utils.optical_alignment import OpticalAlignment
#     al = OpticalAlignment(tt)
#     intMat, rec, mask = al._loadAlignmentInfo()
#     return intMat
#
# def testCalib():
#     '''
#     Returns
#     -------
#     intMat1: numpy array
#         cube of tts1's interaction matrix
#     intMat2: numpy array
#         cube of tts2's interaction matrix
#     '''
#     tts1 = np.array(['20201214_091212', '20201214_092528', '20201214_093842',
#                     '20201214_095152', '20201214_100508', '20201214_101821',
#                     '20201214_103128', '20201214_104441', '20201214_105754',
#                     '20201214_111110', '20201214_112435', '20201214_113749'])
#     tts2 = np.array(['20201214_115451', '20201214_120323', '20201214_121200',
#                     '20201214_122040', '20201214_122922', '20201214_123807',
#                     '20201214_124640', '20201214_125504', '20201214_130327',
#                     '20201214_131134', '20201214_131950', '20201214_132822'])
#     intMat1 = None
#     intMat2 = None
#     for i in range(tts1.size):
#         mat1 = pippo(tts1[i])
#         mat2 = pippo(tts2[i])
#         if intMat1 is None:
#             intMat1 = mat1
#         else:
#             intMat1 = np.stack((intMat1, mat1))
#         if intMat2 is None:
#             intMat2 = mat2
#         else:
#             intMat2 = np.stack((intMat2, mat2))
#     return intMat1, intMat2
#
#     def analyzeOptRep(tt):
#         par, rm, cube = ._readRepData(tt)
#         z_list=[]
#         for i in range(cube.shape[2]):
#             masked_ima = cube[:,:,i]
#             coef, mat = zernike.zernikeFit(masked_ima, np.arange(2, 7))
#             z_list.append(coef)
#         return np.array(z_list)

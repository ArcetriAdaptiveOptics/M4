"""
Authors
  - C. Selmi:  written in 2020
"""

import os
import logging
import numpy as np
from astropy.io import fits as pyfits
from m4.ground import tracking_number_folder
from m4.configuration import config_folder_names as fold_name
from m4.utils import osutils as osu
from m4.configuration.ott_parameters import Interferometer
from m4.ground.read_data import InterferometerConverter
from m4.analyzers.analyzer_iffunctions import AnalyzerIFF
from m4.ground import zernike
from m4.ground import read_data

class Noise:
    """
    Class for noise evaluation

    HOW TO USE IT::

        from m4.noise_functions import Noise
        n = Noise()
        #acquisizione dati e analisi dalla cartella hdf5
        tt = n.noise_analysis_from_hdf5_folder(tidy_or_shuffle, template)
        #analisi di pi cartelle di dati
        rms_medio, tip, tilt, n_temp = n.different_template_analyzer(tt_list)
    """

    def __init__(self):
        """The constructor"""
        self._logger = logging.getLogger("NOISE:")
        self._ic = InterferometerConverter()
        self._numberOfNoiseIma = None
        self._cubeFromAnalysis = None

    @staticmethod
    def _storageFolder():
        """Creates the path where to save measurement data"""
        return fold_name.NOISE_ROOT_FOLDER

    def _defAnalyzer(
        self,
        data_file_path,
        tidy_or_shuffle,
        template,
        n_push_pull=None,
        actsVector=None,
    ):
        """
        arg:
            tidy_or_shuffle = (int) 0 per tidy, 1 per shuffle
            template = np.array composed by 1 and -1
            actsVector = vector of actuators or modes
            n_push_pull = (int) number of push pull

        returns:
            an = analyzer object
        """
        an = AnalyzerIFF()
        # lista = glob.glob(os.path.join(data_file_path,'*.h5'))
        # if len(lista)==0:
        #    lista = glob.glob(os.path.join(data_file_path,'*.fits'))
        #    if len(lista)==0:
        #        lista = glob.glob(os.path.join(data_file_path,'*.4D'))
        lista = self._createOrdListFromFilePath(data_file_path)
        n_tot = len(lista)

        an._template = template
        if n_push_pull is None:
            an._nPushPull = 6
            # an._nPushPull = int(n_tot/an._template.size) #cosi il cubo diventa una immagine
        else:
            an._nPushPull = n_push_pull
        if actsVector is None:
            n_acts = int(n_tot / (an._template.size * an._nPushPull))
            an._actsVector = np.arange(n_acts)
            an._modeVector = np.copy(an._actsVector)
        else:
            an._actsVector = actsVector
            an._modeVector = np.copy(an._actsVector)
        an._cmdAmplitude = np.ones(an._actsVector.shape[0])

        indexingList = []
        if tidy_or_shuffle == 0:
            for i in range(an._nPushPull):
                indexingList.append(an._modeVector)
        elif tidy_or_shuffle == 1:
            for j in range(an._nPushPull):
                np.random.shuffle(an._modeVector)
                indexingList.append(an._modeVector)
        an._indexingList = np.array(indexingList)
        return an

    def noise_analysis_from_hdf5_folder(
        self,
        data_file_path,
        tidy_or_shuffle,
        template,
        n_push_pull=None,
        actsVector=None,
    ):
        """
        Parameters
        ----------
            data_file_path: string
                            measurement data folder
            tidy_or_shuffle: int
                            0 for tidy, 1 for shuffle
            template: numpy array
                    vector composed by 1 and -1
        Other Parameters
        ----------------
            actsVector: numpy array, optional
                        vector of actuators or modes
            n_push_pull: int
                        number of push pull

        Returns
        -------
                tt: string
                    tracking number of measurements made
        """
        an = self._defAnalyzer(
            data_file_path, tidy_or_shuffle, template, actsVector, n_push_pull
        )

        dove, tt = tracking_number_folder.createFolderToStoreMeasurements(
            self._storageFolder()
        )

        self._logger.info("Creating analysis in %s", tt)
        self._cubeFromAnalysis = an.createCubeFromImageFolder(data_file_path)
        fits_file_name = os.path.join(dove, "Cube.fits")
        pyfits.writeto(fits_file_name, self._cubeFromAnalysis.data)
        pyfits.append(fits_file_name, self._cubeFromAnalysis.mask.astype(int))

        self._saveInfo(
            dove, tidy_or_shuffle, an._template, an._actsVector, an._nPushPull
        )

        rms_mean, quad_mean, tilt_mean, ptv_mean = self._rmsFromCube(
            self._cubeFromAnalysis
        )
        self._saveResults(rms_mean, quad_mean, ptv_mean, dove)
        return tt

    def _rmsFromCube(self, cube_to_process):
        """
        Parameters
        ----------
            cube_to_process: [pixel, pixel, number of imagescube]
                            cube generated by the analyzer_iffunctions

        Returns
        -------
            rms_mean: numpy array
                    rms averaged on the number of modes used in the iff's acquisition
            tip: numpy array
                tip averaged on the number of modes used in the iff's acquisition
            tilt: numpy array
                tilt averaged on the number of modes used in the iff's acquisition
        """
        self._logger.debug("Calculation of rms, tip and tilt")
        rms_list = []
        coef_tilt_list = []
        coef_tip_list = []
        quad_list = []
        ptv_list = []
        for i in range(cube_to_process.shape[2]):
            image = cube_to_process[:, :, i]
            coef, mat = zernike.zernikeFit(image, np.array([1, 2, 3]))
            sur = zernike.zernikeSurface(image, coef, mat)
            image_ttr = image - sur
            ptv = np.max(image_ttr) - np.min(image_ttr)
            ptv_list.append(ptv)
            rms = image_ttr.std()
            rms_list.append(rms)
            coef_tip_list.append(coef[1])
            coef_tilt_list.append(coef[2])
            quad = np.sqrt(coef[1] ** 2 + coef[2] ** 2)
            quad_list.append(quad)
        ptv_vector = np.array(ptv_list)
        rms_vector = np.array(rms_list)
        tip = np.array(coef_tip_list).mean()
        tilt = np.array(coef_tilt_list).mean()
        quad_tt = np.array(quad_list).mean()
        rms_mean = np.mean(rms_vector)
        ptv_mean = np.mean(ptv_vector)
        tipvec  = np.array(coef_tip_list)
        tiltvec = np.array(coef_tilt_list)
        ttvec = np.array([tiltvec,tipvec])
        return rms_mean, quad_tt, tilt, ptv_mean#,ttvec

    def spectrumAllData(self, data_file_path):
        """
        Parameters
        ----------
        data_file_path: string
                measurement data folder

        Returns
        -------
        spe_tip: numpy array
            spectrum obtained by fft of tip data
        freq_tip: numpy array
            frequency obtained by fft of tip data
        spe_tilt: numpy array
            spectrum obtained by fft of tilt data
        freq_tilt: numpy array
            frequency obtained by fft of tilt data
        """
        lista = self._createOrdListFromFilePath(data_file_path)

        coef_tilt_list = []
        coef_tip_list = []
        for i in range(len(lista)):
            name = "img_%04d" % i
            print(name)
            # file_name = os.path.join(data_file_path, name)
            start_image = image = read_data.read_phasemap(lista[i])
            coef, mat = zernike.zernikeFit(start_image, np.array([2, 3]))
            coef_tip_list.append(coef[0])
            coef_tilt_list.append(coef[1])
        tip = np.array(coef_tip_list)
        tilt = np.array(coef_tilt_list)

        spe_tip, freq_tip = self._fft(tip)
        spe_tilt, freq_tilt = self._fft(tilt)
        return spe_tip, freq_tip, spe_tilt, freq_tilt

    def _fft(self, vector):
        dt = 35e-3
        n = vector.size
        T = n * dt

        spe = np.fft.fftshift(np.fft.rfft(vector, norm="ortho"))
        freq = np.fft.fftshift(np.fft.rfftfreq(vector.size, d=dt))
        # freq = np.fft.fftshift(np.fft.rfftfreq(spe.size, d=dt))
        #         res = np.fft.fft(rms_vector)
        #         qq = res.real**2+res.imag**2
        return spe, freq

    def _saveResults(self, rms_mean, quad_mean, ptv_mean, destination_file_path):
        """Save results as text file"""
        fits_file_name = os.path.join(destination_file_path, "results.txt")
        file = open(fits_file_name, "w+")
        file.write("%e %e %e" % (rms_mean, quad_mean, ptv_mean))
        file.close()

    def _saveInfo(self, dove, tidy_or_shuffle, template, actsVector, n_push_pull):
        """Save measurement data as file fits"""
        fits_file_name = os.path.join(dove, "Info.fits")
        header = pyfits.Header()
        header["NPUSHPUL"] = n_push_pull
        header["TIDYSHUF"] = tidy_or_shuffle
        pyfits.writeto(fits_file_name, template)
        pyfits.append(fits_file_name, actsVector)
        return

    def _readCube(self, tt):
        """
        args:
            tt = tracking number of measurement

        return:
            _cubeFromAnalysis = cube obtained after iff analysis
        """
        store_in_folder = Noise._storageFolder()
        file_path = os.path.join(store_in_folder, tt)
        fits_file_name = os.path.join(file_path, "Cube.fits")
        hduList = pyfits.open(fits_file_name)
        self._cubeFromAnalysis = np.ma.masked_array(
            hduList[0].data, hduList[1].data.astype(bool)
        )
        return self._cubeFromAnalysis

    def _readTempFromInfoFile(self, tt):
        """
        args:
            tt = tracking number of measurement

        return:
            n_temp = (int) length of template vector
        """
        store_in_folder = Noise._storageFolder()
        file_path = os.path.join(store_in_folder, tt)
        fits_file_name = os.path.join(file_path, "Info.fits")
        hduList = pyfits.open(fits_file_name)
        n_temp = hduList[0].data.shape[0]
        return n_temp

    def different_template_analyzer(self, tt_list):
        """
        Parameters
        ----------
            tt_list: lista
                    lista of tracking number to analyze

        Returns
        -------
            rms_medio: numpy array
                     vector of mean rms (one for each data folder)
            n_tempo: numpy array
                    vector of the length of the templates used
        """
        self._logger.info("Analysis whit different template")
        self._logger.debug("tt_list used: %s", tt_list)
        rms_list = []
        ptv_list = []
        quad_list = []
        tilt_list = []
        n_temp_list = []
        for tt in tt_list:
            cube = self._readCube(tt)
            n_temp = self._readTempFromInfoFile(tt)
            rms, quad, tilt, ptv = self._rmsFromCube(cube)
            rms_list.append(rms)
            quad_list.append(quad)
            tilt_list.append(tilt)
            ptv_list.append(ptv)
            n_temp_list.append(n_temp)

        rms_medio = np.array(rms_list)
        quad = np.array(quad_list)
        tilt = np.array(tilt_list)
        n_temp = np.array(n_temp_list)
        ptv_medio = np.array(ptv_list)
        return rms_medio, quad, n_temp, ptv_medio

    ### tt_list ###
    # measurementFolder ='/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Noise'
    # lista= os.listdir(measurementFolder); lista.sort()
    # tt_list = lista[1:len(lista)-2]
    ### PLOT ###
    # plot(n_temp, rms, '-o'); plt.xlabel('n_temp'); plt.ylabel('rms_medio')

    ### FUNZIONE DI STRUTTURA ###
    def comp_spatial_stf(self, image_ma, pixsc=1e-3, nbins=None):
        # lista = self._createOrdListFromFilePath(data_file_path)

        img = image_ma.data
        mask = image_ma.mask == 0
        sy = img.shape[0]
        sx = img.shape[1]

        if nbins is None:
            nbins = int(sy / 10)

        sy, sx = img.shape

        vec = img[0: int(sx / 2) + 1, int(sy / 2)]
        for posiz in np.arange(int(sx / 2), 0, -1):
            val = vec[posiz]
            if ~np.isnan(val):
                # print(posiz, val)
                break
        xx, yy = np.meshgrid(np.arange(sx) - posiz, np.arange(sy) - sy / 2)
        rr = np.sqrt(xx**2 + yy**2)
        data2hist = (img[mask] - val) ** 2
        r2hist = rr[mask]
        hist, bin_edges = np.histogram(r2hist, bins=nbins)

        idd = np.digitize(r2hist, bin_edges)
        idd_uniq = np.unique(idd)
        stf = hist * 0.0

        for ii in np.arange(len(stf)):
            data2mean = data2hist[idd == ii]
            if len(data2mean) == 0:
                stf[ii] = 0
            else:
                stf[ii] = np.mean(data2hist[idd == ii])

        return (bin_edges[1:] - (bin_edges[1] - bin_edges[0])) * pixsc, stf

    def analysis_with_structure_function(
        self, data_file_path, tau_vector, h5_or_fits=None, nzern=None
    ):
        """
        .. 4000 = total number of image in hdf5

        Parameters
        ----------
            data_file_path: string
                            measurement data folder
            tau_vector: numpy array
                        vector of tau to use

        Other Parameters
        ----------------
        h5_or_fits: None
            if it is none the .h5 or .4D data analysis is performed
            else the fits analysis is performed
        nzern: int
            if None first three zernike are are subtracted from the differential image
            else the number specified

        Returns
        -------
            rms_medio: numpy array
                    calculated on the difference of the images (distant tau)
            quad_med: numpy array
                     squaring sum of tip and tilt calculated on the difference
                    of the images
        """
        # if h5_or_fits is None:
        #    lista = glob.glob(os.path.join(data_file_path, '*.h5'))
        #    if len(lista)==0:
        #        lista = glob.glob(os.path.join(data_file_path,'*.4D'))
        #    lista.sort()
        # else:
        #    listtot = glob.glob(os.path.join(data_file_path, '*.fits'))
        #    listtot.sort()
        #    lista = listtot[0:-2]
        # print(lista)

        lista = self._createOrdListFromFilePath(data_file_path)
        if nzern is None:
            zv = np.arange(3) + 1
        else:
            zv = np.arange(nzern) + 1

        image_number = len(lista)
        i_max = int(
            (image_number - tau_vector[tau_vector.shape[0] - 1])
            / (tau_vector[tau_vector.shape[0] - 1] * 2)
        )
        if i_max <= 10:
            print("Warning low sampling...")
        #    raise OSError('tau = %s too large. i_max = %d' %(tau_vector[tau_vector.shape[0]-1], i_max))
        rms_medio_list = []
        quad_med_list = []
        for j in range(tau_vector.shape[0]):
            dist = tau_vector[j]
            # print(dist)
            rms_list = []
            quad_list = []
            for i in range(i_max):
                k = i * dist * 2
                if h5_or_fits is None:
                    # k = i * dist * 2
                    # name = 'img_%04d' %k
                    # file_name = os.path.join(data_file_path, name)
                    image_k = read_data.read_phasemap(lista[k])
                    # name = 'img_%04d' %(k+dist)
                    # file_name = os.path.join(data_file_path, name)

                    image_dist = read_data.read_phasemap(lista[k + dist])

                else:
                    image_k = read_data.readFits_maskedImage(lista[k])
                    image_dist = read_data.readFits_maskedImage(
                        lista[k + dist])

                image_diff = image_k - image_dist
                zernike_coeff_array, mat = zernike.zernikeFit(image_diff, zv)
                sur = zernike.zernikeSurface(
                    image_diff, zernike_coeff_array, mat)
                image_ttr = image_diff - sur
                quad = np.sqrt(
                    zernike_coeff_array[0] ** 2 + zernike_coeff_array[1] ** 2
                )

                rms = image_ttr.std()
                rms_list.append(rms)
                quad_list.append(quad)
            rms_vector = np.array(rms_list)
            aa = rms_vector.mean()
            rms_medio_list.append(aa)
            quad_med_list.append(np.array(quad_list).mean())
        rms_medio = np.array(rms_medio_list)
        quad_med = np.array(quad_med_list)

        # per calcolo statistical amplitude of convention
        n_meas = rms_vector.shape[0] * 2 * tau_vector.shape[0]

        return rms_medio, quad_med, n_meas

    # plot(tau_vector, rms, '-o'); plt.xlabel('tau'); plt.ylabel('rms_medio')

    def piston_noise(self, data_file_path):
        """Remove tip and tilt from image and average the results
        .. dovrei vedere una variazione nel tempo

        Parameters
        ----------
            data_file_path: string
                            measurement data folder

        Returns
        -------
            mean: numpy array
                vector containing images's mean
            time: numpy array
                vector of the time at which the image were taken
        """
        lista = self._createOrdListFromFilePath(data_file_path)

        image_number = len(lista)
        time = np.arange(image_number) * (1 / Interferometer.BURST_FREQ)

        mean_list = []
        for j in range(image_number):
            # name = 'img_%04d' %j
            # file_name = os.path.join(data_file_path, name)
            image = image = read_data.read_phasemap(lista[j])
            zernike_coeff_array, mat = zernike.zernikeFit(
                image, np.array([2, 3]))
            sur = zernike.zernikeSurface(image, zernike_coeff_array, mat)
            image_ttr = image - sur
            mean = image_ttr.mean()
            mean_list.append(mean)

        spe, freq = self._fft(np.array(mean_list))
        return np.array(mean_list), time, spe, freq

    def tiptilt_series(self, data_file_path):
        """Remove tip and tilt from image and average the results
        .. dovrei vedere una variazione nel tempo

        Parameters
        ----------
            data_file_path: string
                            measurement data folder

        Returns
        -------
            mean: numpy array
                vector containing images's mean
            time: numpy array
                vector of the time at which the image were taken
        """
        lista = self._createOrdListFromFilePath(data_file_path)
        image_number = len(lista)
        time = np.arange(image_number) * (1 / Interferometer.BURST_FREQ)

        tt_list = []
        for j in range(image_number):
            # name = 'img_%04d.h5' %j
            # file_name = os.path.join(data_file_path, name)
            image = image = read_data.read_phasemap(lista[j])
            coeff, mat = zernike.zernikeFit(image, np.array([1, 2, 3]))

            tt_list.append(coeff)

        tt = np.array(tt_list)
        return tt

    # def _createOrdListFromFilePath(self, data_file_path):

    #     #        lista = glob.glob(os.path.join(data_file_path,'*.h5'))
    #     #        if len(lista)==0:
    #     #            lista = glob.glob(os.path.join(data_file_path,'*.fits'))
    #     #            if len(lista)==0:
    #     #                lista = glob.glob(os.path.join(data_file_path,'*.4D'))
    #     #        lista.sort()
    #     # da controllare ordinamento nel caso di file .4D
    #     lista = th.fileList(None, fold=data_file_path, name="*.h5")
    #     if len(lista) == 0:
    #         lista = th.fileList(None, fold=data_file_path, name="*.4D")
    #         if len(lista) == 0:
    #             lista = th.fileList(None, fold=data_file_path, name="20*.fits")
    #     return lista

    def _createOrdListFromFilePath(self, data_file_path):
        """
        Returns an ordered list of file inside the input path.

        Parameters
        ----------
        data_file_path : str
            File path of the data.

        Returns
        -------
        lista : list of str
            List of data obtained through the file path.

        """
        lista = osu.getFileList(fold=data_file_path, key=".h5")
        if len(lista) == 0:
            lista = osu.getFileList(fold=data_file_path, key=".4D")
            if len(lista) == 0:
                lista = osu.getFileList(fold=data_file_path, key=".fits")
        return lista

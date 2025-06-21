"""
Authors
  - C. Selmi: written in 2020

    Class to be used for alignment of the optical tower
    and the deformable mirror

    HOW TO USE IT::

        from m4 import noise
        noise.noise_vibrations(data_file_path, numbers_array, tidy_or_shuffle)
        or
        noise.spectrumFromData(data_file_path)
        or
        noise.convection_noise(data_file_path, tau_vector)
"""

import os
import time
import glob
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration import folders as config
from m4.configuration.ott_parameters import Interferometer
from m4.analyzers.noise_data_analyzer import Noise


def _path_noise_results(data_file_path, h5_or_fits=None):
    """Function to get tt"""
    results_path = os.path.join(config.NOISE_ROOT_FOLDER, "Noise")
    x = data_file_path.split("/")
    if h5_or_fits is None:
        dove = os.path.join(results_path, x[len(x) - 1])  # da controllare il path
    else:
        dove = os.path.join(results_path, x[len(x) - 1])
    if os.path.exists(dove):
        dove = dove
    else:
        os.makedirs(dove)
    return dove


def _createTemplateList(numbers_array):
    """
    Parameters
    ----------
        numbers_array: numpy array
                    vector containing integers numbers for
                    template creation
    Returns
    -------
        template_list: list
                    list of template to use
    """
    template_list = []
    vec = np.array([1, -1])
    for i in numbers_array:
        if i % 2 == 0:
            # pari
            k = i - 2
            temp = np.tile(vec, int(i / 2))
        elif i % 2 == 1:
            # dispari
            k = i - 2
            if k == 1:
                temp_pari = vec
            else:
                temp_pari = np.tile(vec, int((i - 1) / 2))
            temp = np.append(temp_pari, 1)
        template_list.append(temp)
    return template_list


def noise_vibrations(data_file_path, numbers_array, tidy_or_shuffle=0, show=True):
    """
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
        numbers_array: numpy array
                    vector containing integers numbers for
                    template creation
        tidy_or_shuffle: int
                        0 for tidy, 1 for shuffle

    Returns
    ------
    The output of this function is the plot of results and save this results
    """
    n = Noise()
    dove = _path_noise_results(data_file_path)
    template_list = _createTemplateList(numbers_array)
    tt_list = []
    for temp in template_list:
        tt = n.noise_analysis_from_hdf5_folder(data_file_path, tidy_or_shuffle, temp)
        time.sleep(1)
        tt_list.append(tt)
    fits_file_name = os.path.join(dove, "trackingnumbers_%d.txt" % tidy_or_shuffle)
    file = open(fits_file_name, "w+")
    file.write("Tidy or shuffle = %d \n" % tidy_or_shuffle)
    for tt in tt_list:
        file.write("%s \n" % tt)
    file.close()
    rms_medio, quad_medio, n_temp, ptv_medio = n.different_template_analyzer(tt_list)
    if show:
        pyfits.writeto(
            os.path.join(dove, "rms_vector_%d.fits" % tidy_or_shuffle),
            rms_medio,
            overwrite=True,
        )
        pyfits.writeto(
            os.path.join(dove, "tiptilt_vector_%d.fits" % tidy_or_shuffle),
            quad_medio,
            overwrite=True,
        )
        pyfits.writeto(
            os.path.join(dove, "n_temp_vector_%d.fits" % tidy_or_shuffle),
            n_temp,
            overwrite=True,
        )
        pyfits.writeto(
            os.path.join(dove, "ptv_%d.fits" % tidy_or_shuffle),
            ptv_medio,
            overwrite=True,
        )
        plt.clf()
        # WFE = 2*rms_medio
        plt.plot(n_temp, rms_medio * 1e9, "-o")
        plt.xlabel("n_temp")
        plt.ylabel("rms [nm]")
        plt.title("%s" % tt)
        plt.grid()
        name = os.path.join(dove, "rms_ntemp_%d.png" % tidy_or_shuffle)
        if os.path.isfile(name):
            os.remove(name)
        plt.savefig(name)
        plt.figure()
        plt.plot(n_temp, quad_medio * 1e9, "-o")
        plt.xlabel("n_temp")
        plt.ylabel("TipTilt [nm]")
        plt.title("%s" % tt)
        plt.grid()
        name = os.path.join(dove, "tiptilt_ntemp_%d.png" % tidy_or_shuffle)
        if os.path.isfile(name):
            os.remove(name)
        plt.savefig(name)
        plt.figure()
        plt.plot(n_temp, ptv_medio * 1e9, "-o")
        plt.xlabel("n_temp")
        plt.ylabel("PtV [nm]")
        plt.title("%s" % tt)
        plt.grid()
        name = os.path.join(dove, "ptv_ntemp_%d.png" % tidy_or_shuffle)
        if os.path.isfile(name):
            os.remove(name)
        plt.savefig(name)
    tt = data_file_path.split("/")[-1]
    return [rms_medio, quad_medio, n_temp, ptv_medio]


def spectrumFromData(data_file_path):
    """
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
    """
    print("Spectrum analysis")
    n = Noise()
    dove = _path_noise_results(data_file_path)

    spe_tip, freq_tip, spe_tilt, freq_tilt = n.spectrumAllData(data_file_path)
    plt.figure()
    plt.clf()
    plt.plot(freq_tip, np.absolute(spe_tip), "o")
    plt.xlabel("Freq[HZ]")
    plt.ylabel("|FFT(sig)|")
    plt.title("tip_spectrum")
    name = os.path.join(dove, "tip_spectrum.png")
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)
    plt.figure()
    plt.plot(freq_tilt, np.absolute(spe_tilt), "o")
    plt.xlabel("Freq[HZ]")
    plt.ylabel("|FFT(sig)|")
    plt.title("tilt_spectrum")
    name = os.path.join(dove, "tilt_spectrum.png")
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)


def convection_noise(
    data_file_path,
    tau_vector,
    freq=Interferometer.BURST_FREQ,
    fits_analysis=False,
    nzern=None,
    show=True,
):
    """
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
        tau_vector: numpy array
                    vector of tau to use

    Other Parameters
    ----------------
        fits_analysis: Boolean
            if False the h5 or 4D data analysis is performed
    """
    # last_name = data_file_path.split('/')[-1]
    if fits_analysis is False:
        h5_or_fits = None
    else:
        h5_or_fits = 7
    n = Noise()
    dove = _path_noise_results(data_file_path, h5_or_fits)
    rms, quad, n_meas = n.analysis_whit_structure_function(
        data_file_path, tau_vector, h5_or_fits, nzern=nzern
    )
    rms_nm = rms * 1e9
    if h5_or_fits is None:
        x = tau_vector * (1 / freq)
        param = [5, 0.5, 32]
        try:
            pp, fit = _curvFit(param, x, rms_nm)
            decorr_time = 1 / pp[0] + pp[1]
        except:
            pp = np.array([0, 0, rms[-1] * 1e9])
            decorr_time = -1
            fit = rms_nm.copy() * 0
        if show:
            pyfits.writeto(
                os.path.join(dove, "rms_vector_conv.fits"), rms, overwrite=True
            )
            pyfits.writeto(
                os.path.join(dove, "tiptilt_vector_conv.fits"), quad, overwrite=True
            )
            pyfits.writeto(
                os.path.join(dove, "tau_vector.fits"), tau_vector, overwrite=True
            )
            pyfits.writeto(
                os.path.join(dove, "time_vector_conv.fits"), x, overwrite=True
            )
            plt.figure()
            plt.clf()
            plt.plot(x, rms * 1e9, "-o", label="meas")
            plt.xlabel("time [s]")
            plt.ylabel("rms [nm]")
            plt.plot(x, fit, "-", label="fit")
            plt.grid()
            plt.plot(
                [x[0], x[-1]],
                [pp[2], pp[2]],
                "--r",
                linewidth=3,
                label="%.2f [nm]" % pp[2],
            )
            plt.legend()
            tt = data_file_path.split("/")[-1]
            plt.title("%s" % tt)
            plt.show()
            name = os.path.join(dove, "rms_tau.png")
            if os.path.isfile(name):
                os.remove(name)
            plt.savefig(name)
        return [rms, quad, n_meas, pp[2], x, fit]
    else:
        time_diff = _time_for_plot(data_file_path)
        x = tau_vector * time_diff
        if show:
            pyfits.writeto(
                os.path.join(dove, "time_vector_conv.fits"), x, overwrite=True
            )
            plt.figure()
            plt.clf()
            plt.plot(x, rms * 1e9, "-o", label="time_diff = %d" % time_diff)
            plt.xlabel("time [s]")
            plt.ylabel("rms [nm]")
            plt.grid()
            plt.legend()
            tt = dove.split("/")[-1]
            plt.title("%s" % tt)
            plt.show()
            name = os.path.join(dove, "rms_tau.png")
            if os.path.isfile(name):
                os.remove(name)
            plt.savefig(name)
        return [rms, quad, n_meas]


# stimare tc dal grafico e usare 2*tau_c = epsilon_c / np.sqrt(n) n = 4000
#     tau_c = 30 * (1/27.58)
#     epsilon_c = 2 * tau_c * np.sqrt(n_meas)
#     fits_file_name = os.path.join(dove, 'epsilon_c.txt')
#     file = open(fits_file_name, 'w+')
#     file.write('Epsilon_c = %e' %epsilon_c)
#     file.close()


def _time_for_plot(stab_path):
    listtot = glob.glob(os.path.join(stab_path, "*.fits"))
    listtot.sort()
    aa = listtot[0].split("/")
    t0 = aa[-1].split("_")[1].split(".")[0]
    bb = listtot[1].split("/")
    t1 = bb[-1].split("_")[1].split(".")[0]
    hs = float(t0[0:2]) * 3600
    ms = float(t0[2:4]) * 60
    s = float(t0[4::])
    t0s = hs + ms + s
    hs = float(t1[0:2]) * 3600
    ms = float(t1[2:4]) * 60
    s = float(t1[4::])
    t1s = hs + ms + s
    time_diff = t1s - t0s
    return time_diff


def _funFit(x, a, b, c):
    fun = -np.exp(-a * (x - b)) + c
    return fun


def _curvFit(param, x, rms_nm):
    from scipy.optimize import curve_fit

    pp, pcov = curve_fit(_funFit, x, rms_nm, param)
    fit = _funFit(x, *pp)
    return pp, fit


def piston_noise(data_file_path):
    """
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
    """
    n = Noise()
    dove = _path_noise_results(data_file_path)

    mean, time, spe, freq = n.piston_noise(data_file_path)
    pyfits.writeto(os.path.join(dove, "piston_vector.fits"), mean)
    pyfits.writeto(os.path.join(dove, "time_vector.fits"), time)

    plt.clf()
    plt.plot(time, mean)
    plt.xlabel("time[s]")
    plt.ylabel("mean_image")
    plt.savefig(os.path.join(dove, "piston_noise.png"))
    plt.figure()
    plt.plot(freq, np.absolute(spe), "o")
    plt.xlabel("Freq[HZ]")
    plt.ylabel("|FFT(sig)|")
    plt.title("piston_power_spectrum")
    plt.savefig(os.path.join(dove, "piston_spectrum.png"))

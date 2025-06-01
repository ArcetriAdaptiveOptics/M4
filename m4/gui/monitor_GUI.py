"""
Authors
  - C. Selmi: written in 2023
"""

from guietta import Gui, G, _, ___, HB, MA, III
from m4.analyzers.noise_data_analyzer import Noise
from m4.configuration.ott_parameters import Interferometer
from scipy.optimize import curve_fit
import sys
import time
import numpy as np
from m4.ground import read_data
from m4.ground import zernike


class Runner:
    """
    Class for creating interactive GUI

    HOW TO USE IT::

        conf = 'myConfGlobalPath.yaml'
        from m4.configuration import start
        ott, interf, dm = start.create_ott(conf)
        from m4.gui import monitor_GUI
        g = monitor_GUI.Runner()
        g.run()
    """

    def __init__(self):
        """The constructor"""

    # file_path = '/Users/cselmi/Desktop/Arcetri/M4/Data/M4Data/20201215_140200_noise/hdf5'
    def _setUp(self):

        def Start_Mean_Burst(gui, *args):
            """current alignment value (Zernike: TipTilt focus, coma): media burst e fit modi di allineamento."""
            n = Noise()
            lista = n._createOrdListFromFilePath(np.str(gui.file_path))
            # cube_list = []
            for i in range(len(lista)):
                name = "img_%04d" % i
                print(name)
                image = read_data.read_phasemap(lista[i])
                if i == 0:
                    ima_sum = image.data
                    mask = image.mask
                else:
                    ima_sum += image.data
                    mask = np.ma.mask_or(mask, image.mask)
                # cube_list.append(image)
            # cube = np.ma.dstack(cube_list)
            # ima_mean = np.ma.mean(cube, axis=2)
            ima_mean = np.ma.masked_array(ima_sum / len(lista), mask=mask)
            gui.stdBurst = np.std(ima_mean)
            gui.plot3 = ima_mean
            gui.plot3.set_clim(vmin=ima_mean.min(), vmax=ima_mean.max())
            gui.plot3.set_title("Media burst")
            zernike_coeff_array, mat = zernike.zernikeFit(ima_mean, np.arange(11) + 1)
            gui.tip_value = zernike_coeff_array[0]
            gui.tilt_value = zernike_coeff_array[1]
            gui.focus_value = zernike_coeff_array[2]
            gui.coma1_value = zernike_coeff_array[6]
            gui.coma2_value = zernike_coeff_array[7]

        def Start_Differential_TipTilt(gui, *args):
            """current differential tip/tilt (at the acquisition frequency): noise.noise_vibrations with fixed template."""
            data_file_path = np.str(gui.file_path)
            data_tt = data_file_path.split("/")[-2]
            n = Noise()
            template_list = [
                np.array([1, -1, 1]),
                np.array([1, -1, 1, -1, 1]),
                np.array([1, -1, 1, -1, 1, -1, 1]),
            ]

            tt_list = []
            for temp in template_list:
                tt = n.noise_analysis_from_hdf5_folder(data_file_path, 0, temp)
                time.sleep(1)
                tt_list.append(tt)
            rms_medio, quad_medio, n_temp, ptv_medio = n.different_template_analyzer(
                tt_list
            )
            # noise.noise_vibrations(data_file_path, np.array([23,25,27]), 0)
            setTipTiltPlot(gui, data_tt, n_temp, rms_medio)

        def setTipTiltPlot(gui, tt, n_temp, rms_medio):
            ax = gui.plot.ax
            ax.clear()
            ax.set_title(tt)

            ax.plot(n_temp, rms_medio * 1e9, "-o")
            ax.set_xlabel("n_temp")
            ax.set_ylabel("rms [nm]")
            ax.figure.canvas.draw()
            # ax.grid()

        def Start_Decorrelation_Time(gui, *args):
            """current decorrelation time: noise.convection_noise. Output: grafico."""
            data_file_path = np.str(gui.file_path)
            data_tt = data_file_path.split("/")[-2]
            n = Noise()

            tau_vector = np.arange(1, 60)
            rms, quad, n_meas = n.analysis_whit_structure_function(
                data_file_path, tau_vector, h5_or_fits=None, nzern=None
            )

            rms_nm = rms * 1e9
            x = tau_vector * (1 / Interferometer.BURST_FREQ)
            param = [5, 0.5, 32]
            try:
                pp, pcov = curve_fit(_funFit, x, rms_nm, param)
                fit = _funFit(x, *pp)
                decorr_time = 1 / pp[0] + pp[1]
            except:
                pp = np.array([0, 0, rms[-1] * 1e9])
                decorr_time = -1
                fit = rms_nm.copy() * 0
            setDecorrelationTimePlot(gui, data_tt, rms, x, fit, pp)

        def _funFit(x, a, b, c):
            fun = -np.exp(-a * (x - b)) + c
            return fun

        def setDecorrelationTimePlot(gui, tt, rms, x, fit, pp):
            ax = gui.plot2.ax
            ax.clear()
            ax.set_title(tt)
            ax.plot(x, rms * 1e9, "-o", label="meas")
            ax.set_xlabel("time [s]")
            ax.set_ylabel("rms [nm]")
            ax.plot(x, fit, "-", label="fit")
            ax.plot(
                [x[0], x[-1]],
                [pp[2], pp[2]],
                "--r",
                linewidth=3,
                label="%.2f [nm]" % pp[2],
            )
            ax.legend()
            ax.figure.canvas.draw()

        def Start_Decorrelation_Noise(gui, *args):
            """current decorrelated noise: diff delle immagini dopo tau (prima - ultima). Output: image and RMS"""
            n_start_ima = 0
            tau = int(gui.tau or 3)
            n_after_tau_ima = int(tau * Interferometer.BURST_FREQ)

            n = Noise()
            lista = n._createOrdListFromFilePath(np.str(gui.file_path))
            start_ima = read_data.read_phasemap(lista[n_start_ima])
            after_tau_ima = read_data.read_phasemap(lista[n_after_tau_ima])
            diff_ima = start_ima - after_tau_ima
            gui.stdDiff = np.std(diff_ima)
            gui.plot4 = diff_ima
            gui.plot4.set_clim(vmin=diff_ima.min(), vmax=diff_ima.max())

        def Start_Convection_Footprint_PSD(gui, *args):
            """convection footprint PSD: codice di Marchino che restituisce PSD da immagine del punto 4 (decorrelation noise)"""
            n_start_ima = 0
            n_after_tau_ima = int(int(gui.tau) * Interferometer.BURST_FREQ)

            n = Noise()
            lista = n._createOrdListFromFilePath(np.str(gui.file_path))
            start_ima = read_data.read_phasemap(lista[n_start_ima])
            after_tau_ima = read_data.read_phasemap(lista[n_after_tau_ima])
            diff_ima = start_ima - after_tau_ima

        def Start_Convective_Regions(gui, *args):
            """convective regions (pixel-wise stdev): cubo delle immagini senza tip/tilt, std nella direzione giusta che restituisce un'immagine."""
            # troppo lento
            # n = Noise()
            # lista = n._createOrdListFromFilePath(np.str(gui.file_path))
            # cube_list = []
            # for i in range(len(lista)):
            #     name = 'img_%04d' %i
            #     print(name)
            #     image = read_data.read_phasemap(lista[i])
            #     cube_list.append(image)
            #
            # ima_list = []
            # for image in cube_list:
            #     print(i)
            #     coef, mat = zernike.zernikeFit(image,
            #                                    np.arange(11)+1)
            #     surf = zernike.zernikeSurface(image, coef[:2], mat[:,:2])
            #     image_ttr = image - surf
            #     ima_list.append(image_ttr)
            #
            # cube = np.ma.dstack(ima_list)
            # imm = np.std(cube, axis=2)
            # gui.plot6 = imm
            # gui.plot6.set_clim(vmin=imm.min(), vmax=imm.max())
            pass

        def Start_All_The_Analysis(gui, *args):
            """Start all the analysis"""
            # non funziona per via del tau in ingresso
            Start_Mean_Burst(gui, *args)
            Start_Differential_TipTilt(gui, *args)
            Start_Decorrelation_Time(gui, *args)
            Start_Decorrelation_Noise(gui, *args)
            Start_Convection_Footprint_PSD(gui, *args)
            Start_Convective_Regions(gui, *args)

        gui_analysis = Gui(
            [
                "Set data file path",
                "__file_path__",
                ___,
                ___,
                Start_All_The_Analysis,
                ___,
                ___,
                ___,
            ],
            [Start_Mean_Burst, ___, ___, ___, MA("plot3"), ___, ___, ___],
            ["RMS value:", "stdBurst", _, " ", III, III, III, III],
            ["Zernike:", _, "tip_value", _, "tilt_value", " ", "focus_value", " "],
            [_, _, "coma1_value", _, "coma2_value", _, _, _],
            [
                Start_Differential_TipTilt,
                ___,
                ___,
                ___,
                Start_Decorrelation_Time,
                ___,
                ___,
                ___,
            ],
            [MA("plot"), ___, ___, ___, MA("plot2"), ___, ___, ___],
            ["Set tau [s]", _, "__tau__:3", _, _, _, _, _],
            [
                Start_Decorrelation_Noise,
                ___,
                ___,
                ___,
                Start_Convection_Footprint_PSD,
                ___,
                ___,
                ___,
            ],
            [MA("plot4"), ___, ___, ___, MA("plot5"), ___, ___, ___],
            ["RMS value:", "stdDiff", _, _, _, _, _, _],
            [Start_Convective_Regions, ___, ___, ___, MA("plot6"), ___, ___, ___],
        )

        self.gui = Gui([G("OTT_Monitor_GUI")])
        self.gui.OTT_Monitor_GUI = gui_analysis
        for i in range(8):
            gui_analysis.layout().setColumnStretch(i, 10)

    def run(self):
        """Run the GUI"""
        self._setUp()
        self.gui.run()


def main():
    from m4.configuration import ott

    ott, interf, dm = ott.create_ott()
    runner = Runner()
    sys.exit(runner.run())


if __name__ == "__main__":
    main()

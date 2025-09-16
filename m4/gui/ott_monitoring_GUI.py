"""
Author(s)
---------
- Pietro Ferraiuolo : written in 2024

Description
-----------
Monitoring application for the M4's Optical Test Tower, complete with GUI.

How to Use it
-------------
To run the GUI, simply initialize the class and use the 'rungui()' method.
We need an interferometer object initialized, which we will have if the initialization
of the tower has ben executed:

>>> from m4.gui.monitor import SystemMonitoring
>>> mm = systemMonitoring(interf)
Monitoring interferometer configuration loaded.
>>> mm.rungui()

Alternatively, to run a single sequence of monitoring, the method can be launched::

    >>> mm.monitoring()
"""

import time, os, shutil, threading, numpy as np
from matplotlib.pyplot import tight_layout
from guietta import Gui, _, ___, III, M, G, P, execute_in_main_thread
from m4 import noise
from m4.configuration.root import folders as fn
from m4.configuration import userconfig as uc
from opticalib.ground.osutils import newtn as ts
from opticalib.ground import zernike as zern
from opticalib.ground import osutils as osu


class SystemMonitoring:
    """
    Class for the continuous monitoring of the OTT.
    The only public method is the actual monitoring task.

    Methods
    -------
    monitoring() :
        Sets of actions, acquisition and analysis, which monitors the OTT

    How to Use it
    -------------
    To run a single monitoring event, simply use the 'monitoring' method

    >>> from m4.gui.monitor import SystemMonitoring
    >>> # As we have initialized the interferometer device...
    >>> mm = SystemMonitoring(interf)
    Monitoring interferometer configuration loaded.
    >>> mm.monitoring() # single loop operation

    As for the continuous loop monitoring, with integrated GUI, run the 'rungui'
    method:

    >>> mm.rungui()
    """

    def __init__(self, interferometer):
        """The Constructor"""
        self._llog = os.path.join(fn.MONITORING_ROOT_FOLDER, "Monitor_CompleteLog.txt")
        self._slog = os.path.join(fn.MONITORING_ROOT_FOLDER, "Monitor_ShortLog.txt")
        self.template = np.array([3, 5, 7, 9])
        self.delay = 2
        self.n_frames = 5
        self.interf = interferometer
        self.cam_info = None
        self.freq = None  # Hz - Gets updated every measure
        self.fast_data_path = None
        self.slow_data_path = None
        self.__load_monitor_config()
        self.__update_interf_settings()

        self.is_monitoring = False
        self.monitoring_thread = None

    def rungui(self):
        """
        Method to run the Graphical User Interface of the monitoring app for the
        OTT.

        GUI Structure
        =============
        Plots :
            The left column is populated with the plots resulting from the data
            analysis.
        Results :
            A sub-gui at the top-right, labeled as 'RESULTS', is where data analysis
            numerical results are shown.
        Camera info :
            A sub-gui at the center-right, labeled as 'CAMERA INFO', contains the
            information relative to the camera settings, such as frequency, frame
            pixels and offsets.
        Control :
            A sub-gui at the bottom-right, labeled as 'CONTROL', contains the
            live logging messages of the events, a progress bar showing the state
            of the monitoring loop and three buttons:
                start : A start monitoring button. Once pressed, the script will
                    the looping event of monitoring
                stop : A stop button, which will stop the script from going to
                    the next loop event.
                    Note: if the monitoring script is still running when the
                    stop button is pressed, this will cause the GUI to freeze
                    until the threads are done. Just wait until it finishes.
                close: the close button, upload to the interferometer its default
                    configuration and exit the program by closing the GUI.
        """
        gui = Gui(
            [_, _, M("plot1"), ___, ___, G("RESULTS"), ___],
            [M("plot1"), ___, III, III, III, III, III],
            [III, III, M("plot2"), ___, ___, G("CAMERA INFO"), ___],
            [M("plot1"), ___, III, III, III, III, III],
            [III, III, M("plot3"), ___, ___, G("CONTROL"), ___],
            [_, _, III, III, III, III, III],
        )
        results_gui = Gui(
            ["Slow Mean RMS", "res1", "nm"],
            ["Slow STD", "res2", "nm"],
            ["Fast TT rms", "res3", "nm"],
            ["Slow TT rms", "res4", "nm"],
        )
        camera_gui = Gui(
            ["Frequency", "cam1", "Hz"],
            ["Frame Width", "cam2", "px"],
            ["Frame Height", "cam3", "px"],
            ["X-Offset", "cam4", "px"],
            ["Y-Offset", "cam5", "px"],
        )
        control_gui = Gui(
            ["message", ___, ___],
            [P("progress"), ___, ___],
            [["start"], ["stop"], ["close"]],
        )
        gui.RESULTS = results_gui
        gui.CAMERAINFO = camera_gui
        gui.CONTROL = control_gui
        control_gui.main_gui = gui
        control_gui.results = results_gui
        control_gui.camera = camera_gui

        @execute_in_main_thread(gui)
        def show_results():
            results_gui.res1 = f"{self.slow_results[0]*1e9:.2f}"
            results_gui.res2 = f"{self.slow_results[1]*1e9:.2f}"
            results_gui.res3 = f"{self.fast_tt[:,1].std()*1e9:.2f}"
            results_gui.res4 = f"{self.slow_tt[:,1].std()*1e9:.2f}"
            camera_gui.cam1 = f"{self.freq:.1f}"
            camera_gui.cam2 = f"{self.cam_info[0]:d}"
            camera_gui.cam3 = f"{self.cam_info[1]:d}"
            camera_gui.cam4 = f"{self.cam_info[2]:d}"
            camera_gui.cam5 = f"{self.cam_info[3]:d}"
            plot1(gui)
            plot2(gui)
            plot3(gui)
            plot4(gui)
            plot5(gui)

        @execute_in_main_thread(gui)
        def _show_progress(percentage, text):
            control_gui.widgets["progress"].setValue(percentage)
            control_gui.widgets["message"].setText(text)

        def start(gui, *args):
            if not self.is_monitoring:
                self.is_monitoring = True
                self.monitoring_thread = threading.Thread(
                    target=monitoring_loop, args=(gui,)
                )
                self.monitoring_thread.start()
                gui.widgets["start"].setEnabled(False)

        def stop(gui, *args):
            self.is_monitoring = False
            if self.monitoring_thread:
                self.monitoring_thread.join()
                gui.widgets["start"].setEnabled(True)

        def close(sub_gui, *args):
            self.is_monitoring = False
            if self.monitoring_thread:
                self.monitoring_thread.join()
            sub_gui.main_gui.close()
            self.__load_default_config()

        def monitoring_loop():
            while self.is_monitoring:
                self.monitoring(progress=_show_progress)
                show_results()
                time.sleep(3.5)

        def plot1(gui, *args):
            ax = gui.plot1.ax
            ax.clear()
            ax.set_title("RMS")
            ax.set_xlabel("Frames per Template")
            ax.set_ylabel("Root Mean Square [nm]")
            ax.plot(
                self.fast_results["vn"]["ntemp"],
                self.fast_results["vn"]["rms"],
                "-o",
                c="black",
            )
            ax.grid()
            ax.set_in_layout(tight_layout)
            ax.figure.canvas.draw()

        def plot2(gui, *args):
            ax = gui.plot2.ax
            ax.clear()
            ax.set_title("Tip-Tilt Quadratic Residual")
            ax.set_xlabel("Frames per Template")
            ax.set_ylabel("Tip - Tilt [nm]")
            ax.plot(
                self.fast_results["vn"]["ntemp"],
                self.fast_results["vn"]["tt"],
                "-o",
                c="black",
            )
            ax.grid()
            ax.set_in_layout(tight_layout)
            ax.figure.canvas.draw()

        def plot3(gui, *args):
            ax = gui.plot3.ax
            ax.clear()
            ax.set_title("Decorrelation Noise Fit")
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("Root Mean Square [nm]")
            ax.plot(
                self.fast_results["cn"]["x"],
                self.fast_results["cn"]["rms"] * 1e9,
                "-o",
                c="black",
                label="Measurements",
            )
            x = self.fast_results["cn"]["x"]
            pp = self.fast_results["cn"]["pp"]
            ax.plot(
                [x[0], x[-1]],
                [pp, pp],
                "--",
                linewidth=3,
                color="red",
                label=f"{pp:.2f} [nm]",
            )
            ax.plot(x, self.fast_results["cn"]["fit"])
            ax.grid()
            ax.legend(loc="best", fontsize="small")
            ax.set_in_layout(tight_layout)
            ax.figure.canvas.draw()

        def plot4(gui, *args):
            ax = gui.plot4.ax
            ax.clear()
            ax.set_title("Fast Tip-Tilt")
            ax.set_xlabel("")
            ax.set_ylabel("Tip-Tilt [nm]")
            ax.scatter(self.fast_tt[:, 0], self.fast_tt[:, 1], s=7.5)
            ax.set_in_layout(tight_layout)
            ax.set_xlim(-1e-6, 1e-6)
            ax.figure.canvas.draw()

        def plot5(gui, *args):
            ax = gui.plot5.ax
            ax.clear()
            ax.set_title("Slow Tip-Tilt")
            ax.set_xlabel("")
            ax.set_ylabel("Tip-Tilt [nm]")
            ax.scatter(self.slow_tt[:, 0], self.slow_tt[:, 1], s=7.5)
            ax.set_in_layout(tight_layout)
            # ax.set_xlim(-1e-6, 1e-6)
            ax.figure.canvas.draw()

        gui.events(
            [_, _, plot1, _, _, _, _],
            [plot4, _, _, _, _, _, _],
            [_, _, plot2, _, _, _, _],
            [plot5, _, _, _, _, _, _],
            [_, _, plot3, _, _, _, _],
            [_, _, _, _, _, _, _],
        )
        control_gui.events(
            [_, _, _],
            [_, _, _],
            [start, stop, close],
        )
        gui.title("OTT Monitoring")
        gui.window().resize(1250, 950)
        gui.run()

    def monitoring(self, progress=None):
        """
        Monitoring task for the OTT. It performs two types of acquisitions, a
        'fast acquisition' which does a 'capture' of images at the current
        frequency for 3 seconds, and a 'slow acquisition', which acquires 5
        frames delayed by 2 seconds. The corresponding 'slow' and 'fast' analysis
        are performed on these sets of data.
        This will produce plots and print out results in two monitoring logs, as
        well as ultimately remove all the acquired data, as they are not usefull
        anymore.
        """
        if progress:
            progress(0, "Restarting Acquisition")
            time.sleep(0.25)
        self.__update_interf_settings()
        if progress:
            progress(10, "Interferometer information updated")
            time.sleep(0.5)
            progress(20, "Fast data acquisition started")
        self._fast_acquisition()
        if progress:
            progress(
                30, "Fast data acquisition Completed\nInitiating Slow data acquisition"
            )
        self._slow_acquisition()
        if progress:
            progress(50, "Slow data acquisition completed")
            progress(55, "Slow data analysis in progress...")
        self._slow_analysis()
        if progress:
            progress(60, "Slow data analysis completed")
            progress(65, "Fast data analysis in progress...")
        self._fast_analysis()
        if progress:
            progress(80, "Fast data analysis completed")
        self.__write_log_message()
        if progress:
            progress(85, "Log report written")
            time.sleep(0.25)
        self.__clear_data_folder(self.slow_data_path)
        self.__clear_data_folder(self.fast_data_path)
        if progress:
            progress(95, "Data cleaning completed")
            time.sleep(0.5)
            progress(100, "Monitoring loop completed")

    def _fast_analysis(self):
        """
        Fast analysis performed on the fast-acquired data. It basically performs
        the differential algorythm on the set of data, with different templates,
        and returns the mean RMS for each template, as well as the quadratic
        resitual mean tip-tilt.

        Returns
        -------
        results : dict
            A dictionary of dictionaries for each analysis results. It contains
            res1 : dict
                Results of the 'noise_vibration' function, which are:
                    rms_medio : Root mean square
                    quad_medio : tip-tilt quadratic residual
                    n_temp : number of frames per template
                    ptv_medio : peak-to-valley mean value
            res2 : dict
                Results of the 'convection_noise' function, which are:
                    rms :
                    quad :
                    n_meas :
                    pp :
                    decor_time :
        """
        tos = 0
        numbers_array = self.template
        tau_vector = np.arange(1, 21, 2)
        par1 = noise.noise_vibrations(
            self.fast_data_path, numbers_array, tidy_or_shuffle=tos, show=False
        )
        par2 = noise.convection_noise(
            self.fast_data_path, tau_vector, freq=self.freq, show=False
        )
        keys1 = ["rms", "tt", "ntemp", "ptv"]
        keys2 = ["rms", "tt", "nmeas", "pp", "x", "fit"]
        keys3 = ["vn", "cn"]
        res1 = dict(zip(keys1, par1))
        res2 = dict(zip(keys2, par2))
        tt = self._abs_tilt_analysis(self.fast_data_path)
        self.fast_tt = tt
        self.fast_results = dict(zip(keys3, (res1, res2)))

    def _slow_analysis(self):
        """
        Analysis performed on the slow-acquired data. It takes the images in couples
        and subtracts them removing zernike modes, and computes the results
        standard deviation, which, for each couple, is stored in a list.

        Returns
        -------
        svec : list of floats
            List of standard deviations of the zernike-removed subtracted images.
        """
        gap = 2
        fl = osu.getFileList(fold=self.slow_data_path)
        nfile = len(fl)
        npoints = int(nfile / gap)
        tt = self._abs_tilt_analysis(self.slow_data_path)
        slist = []
        for i in range(0, npoints):
            q0 = osu.read_phasemap(fl[i * gap])
            q1 = osu.read_phasemap(fl[i * gap + 1])
            diff = zern.removeZernike(q1 - q0)
            slist.append(diff.std())
        svec = np.array(slist)
        mean = np.mean(svec)
        std = np.std(svec)
        self.slow_tt = tt
        self.slow_results = [mean, std]

    def _abs_tilt_analysis(self, path):
        """


        Parameters
        ----------
        path : TYPE
            DESCRIPTION.

        Returns
        -------
        tt : TYPE
            DESCRIPTION.

        """
        fl = osu.getFileList(fold=path)
        nf = len(fl)
        tt = np.zeros([nf, 2])
        for i in range(0, nf):
            q = osu.read_phasemap(fl[i])
            coeff, _ = zern.zernikeFit(q, [1, 2, 3])
            tt[i, :] = coeff[1:]
        return tt

    def _fast_acquisition(self):
        """
        Fast acquisition procedure. Capture of freq*3 images, with subsequent
        produce. The tn path is stored in the 'fast_data_path' variable.
        """
        n_frames = int(self.freq * 3)
        tn = self.interf.capture(n_frames)
        self.interf.produce(tn)
        self.fast_data_path = os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn)
        print("Fast acquisition completed.")

    def _slow_acquisition(self):
        """
        Slow acquisition procedure. Acquires a phasemap every 2 seconds, for a
        total of 5 images. The tn path is stored in the 'slow_data_path' variable.
        """
        n_frames = self.n_frames
        tn = ts.now()
        self.slow_data_path = os.path.join(fn.OPD_SERIES_ROOT_FOLDER, tn)
        os.mkdir(self.slow_data_path)
        i = 0
        while i < (n_frames + 1):
            img = self.interf.acquire_phasemap()
            self.interf.save_phasemap(self.slow_data_path, ts.now() + ".fits", img)
            time.sleep(self.delay)
            i += 1
        print("Slow acquisition completed.")

    def __update_interf_settings(self):
        """
        Reads and loads the current interferometer settings.
        """
        self.freq = self.interf.getFrameRate()
        self.cam_info = self.interf.getCameraSettings()

    def __clear_data_folder(self, path):
        """
        Clear the input datapath with all its content.

        Parameters
        ----------
        path : str
            Complete datapath to be removed, with its content.
        """
        if os.path.exists(path):
            shutil.rmtree(path)

    def __load_monitor_config(self):
        """
        Updates the interferometer settings with the monitoring acquisition
        configuration file, found in 'm4.configuration.userconfig'.
        """
        self.interf.loadConfiguration(uc.phasecam_monitorconfig)
        print("Monitoring interferometer configuration loaded.")

    def __load_default_config(self):
        """
        Updates the interferometer settings with a default configuration file,
        found in 'm4.configuration.userconfig'.
        """
        self.interf.loadConfiguration(uc.phasecam_baseconfig)
        print("Default interferometer configuration loaded.")

    def __write_log_message(self):
        """
        Function which writes both the complete and the short monitoring log files.
        """
        tn = ts.now()
        long_msg = f"""#______________________________________________________________________________
{tn}
[Diagnosis]
Mean : {self.slow_results[0]*1e9:.2f}nm ; Std = {self.slow_results[1]*1e9:.2f}nm ; res3 = () ; res4 = ()
[Camera Settings]
FrameRate = {self.freq:.1f} ; Width = {self.cam_info[0]}px ; Height={self.cam_info[1]}px
x-offset = {self.cam_info[2]}px ; y-offset={self.cam_info[3]}px
"""
        short_msg = f"""{tn}    {self.slow_results[0]*1e9:.2f}nm    {self.slow_results[1]*1e9:.2f}nm    res3    res4    res5"""
        with open(self._llog, "a", encoding="utf-8") as log:
            log.write(long_msg)
        with open(self._slog, "a", encoding="utf-8") as log:
            log.write(short_msg)


def main():
    """
    Application Launcher
    """
    from m4.configuration.ott import create_ott

    _, interf, _ = create_ott()
    monitoring = SystemMonitoring(interf)
    monitoring.rungui()


if __name__ == "__main__":
    main()

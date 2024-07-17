"""
Author(s)
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import time, schedule, threading, os, shutil, numpy as np
#import tkinter as tk # GUI
import matplotlib.pyplot as plt
from m4.noise import _curvFit, _createTemplateList
from m4.analyzers.noise_data_analyzer import Noise
from m4.ground.timestamp import Timestamp
from m4.configuration import update_folder_paths as ufp
ts = Timestamp()
fn = ufp.folders

class SystemMonitoring():
    """
    """
    def __init__(self, interferometer, timing=1):
        """The Constructor"""
        self._log       = os.path.join(fn.MONITORING_ROOT_FOLDER, 'MonitorLog.txt')
        self._timing    = timing
        self.delay      = 3
        self.n          = Noise()
        self.interf     = interferometer
        self.cam_info   = None
        self.freq       = 20.0 #Hz - Gets updated every measure
        self.fast_data_path = None
        self.slow_data_path = None

    def start_monitoring(self):
        self.__update_interf_settings()
        # Plan task for every 20 seconds
        schedule.every(20).seconds.do()
        # Execute
        stop_run_continuous = self._run_continuously()
        stop_run_continuous.set()

    def _convectionNoiseMonitoring(self):
        tau_vector = np.arange(1, 400, 1/55) #??? how to decide
        rms, quad, n_meas = self.n.analysis_whit_structure_function(
                                               self.fast_data_path, tau_vector)
        rms_nm = rms * 1e9
        x = tau_vector * (1 / self.freq)
        param = [5, 0.5, 32]
        try:
            pp, fit = _curvFit(param, x, rms_nm)
            decorr_time = 1 / pp[0] + pp[1]
        except:
            pp = np.array([0, 0, rms[-1] * 1e9])
            decorr_time = -1
            fit = rms_nm.copy() * 0
        plt.figure()
        plt.plot(x, rms_nm, "-o", label="Measured data")
        plt.plot(x, fit, "-", label="Curve fit")
        plt.plot([x[0], x[-1]], [pp[2], pp[2]], "--r", linewidth=3,
                                                     label=f"{pp[2]:.2f} [nm]")
        plt.grid()
        plt.xlabel("Time [s]")
        plt.ylabel("RMS [nm]")
        plt.legend()
        plt.show()
        out_parameters = [pp[2], rms_nm, decorr_time]
        return out_parameters

    def _vibrationNoiseMonitoring(self):
        tidy_or_shuffle = 0
        numbers_array = np.array([3,11,25,37,51]) #template
        template_list = _createTemplateList(numbers_array) #??? should work
        tt_list = []
        for temp in template_list:
            tt = self.n.noise_analysis_from_hdf5_folder(self.slow_data_path,
                                                        temp, tidy_or_shuffle)
            time.sleep(1)
            tt_list.append(tt)
        rms_medio,quad_medio,n_temp,ptv_medio = self.n.different_template_analyzer(tt_list)
        # 1
        plt.figure()
        plt.plot(n_temp, rms_medio*1e9, "-o")
        plt.xlabel("n_temp")
        plt.ylabel("RMS [nm]")
        plt.grid()
        # 2
        plt.figure()
        plt.plot(n_temp, quad_medio * 1e9, "-o")
        plt.xlabel("n_temp")
        plt.ylabel("TipTilt [nm]")
        plt.grid()
        # 3
        plt.figure()
        plt.plot(n_temp, ptv_medio*1e9, "-o")
        plt.xlabel("n_temp")
        plt.ylabel("Piston [nm]")
        plt.grid()
        out_parameters = [rms_medio, quad_medio, ptv_medio]
        return out_parameters

    def _fast_acquisition(self):
        tn = self.interf.capture(self.freq*2)
        self.interf.produce(tn)
        self.fast_data_path = os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn)
        print("Fast acquisition completed.")

    def _slow_acquisition(self):
        n_frames = 4
        tn = ts.now()
        self.slow_data_path = os.path.join(fn.OPD_SERIES_ROOT_FOLDER, tn)
        for i in range(n_frames):
            img = self.interf.acquire_phasemap()
            self.inter.save_phasemap(self.slow_data_path, ts.now()+'.fits', img)
            time.sleep(self.delay)
        print("Slow acquisition completed.")

    def __update_interf_settings(self):
        self.freq = self.interf.getFrameRate()
        self.cam_info = self.interf.getCameraSettings()

    def __clear_data_folder(self, path):
        if os.path.exists(path):
            shutil.rmtree(path)

    def __write_log_message(self):
        log_msg = \
f"""#______________________________________________________________________________
{ts.now()}
[Diagnosis]
res1 = () ; res2 = () ; res3 = () ; res4 = ()
[Camera Settings]
FrameRate = {self.freq:.1f} ; Width = {self.cam_info[0]}px ; Height={self.cam_info[1]}px
x-offset = {self.cam_info[2]}px ; y-offset={self.cam_info[3]}px
"""
        with open(self._log, 'a', encoding='utf-8') as log:
            log.write(log_msg)

#______________________________________________________________________________
def _run_continuously(self):
    cease_continuous_run = threading.Event()

    class ScheduleThread(threading.Thread):
        @classmethod
        def run(cls):
            while not cease_continuous_run.is_set():
                schedule.run_pending()
                time.sleep(self._timing)

    continuous_thread = ScheduleThread()
    continuous_thread.start()
    return cease_continuous_run

"""
Author(s)
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import time, os, shutil, numpy as np #schedule, threading
from matplotlib import pyplot as plt
from guietta import Gui, _, ___, III, M
from m4 import noise
from m4.analyzers.timehistory import frame
from m4.ground.timestamp import Timestamp
from m4.ground import zernike as zern, read_data as rd
from m4.configuration import update_folder_paths as ufp
from m4.configuration import userconfig as uc
ts = Timestamp()
fn = ufp.folders

class SystemMonitoring():
    """
    """
    def __init__(self, interferometer, timing=1):
        """The Constructor"""
        self._llog      = os.path.join(fn.MONITORING_ROOT_FOLDER, 'Monitor_CompleteLog.txt')
        self._slog      = os.path.join(fn.MONITORING_ROOT_FOLDER, 'Monitor_ShortLog.txt')
        self._timing    = timing
        self.template   = np.array([3,5,7,9])
        self.delay      = 2
        self.n_frames   = 5
        self.interf     = interferometer
        self.cam_info   = None
        self.freq       = None # Hz - Gets updated every measure
        self.fast_data_path = None
        self.slow_data_path = None
        self.__load_monitor_config()
        self.__update_interf_settings()

    def monitoring(self):
        start = time.time()
        self.__update_interf_settings()
        self._fast_acquisition()
        self._slow_acquisition()
        svec = self._slow_analysis()
        spars = []
        spars.append(np.mean(svec))
        spars.append(np.std(svec))
        fpar1, fpar2 = self._fast_analysis()
        plt.figure()
        plt.plot(fpar1[2], fpar1[0], '-x', c='black')
        plt.title('RMS')
        plt.grid()
        plt.show()
        plt.figure()
        plt.plot(fpar1[2], fpar1[1], '-x', c='black')
        plt.title('TipTilt')
        plt.grid()
        plt.show()
        self.__write_log_message(spars)
        self.__clear_data_folder(self.slow_data_path)
        self.__clear_data_folder(self.fast_data_path)
        self.__load_default_config() # da fare alla fine del monitoring, non della funzione
        end = time.time()
        print(f'Execution time: {end-start:.2f} seconds')

    def start_monitoring(self):
        self.__update_interf_settings()
        # Plan task for every 20 seconds
        schedule.every(20).seconds.do()
        # Execute
        stop_run_continuous = self._run_continuously()
        stop_run_continuous.set()

    def _fast_analysis(self):
        tos = 0
        numbers_array = self.template
        tau_vector = np.arange(1,21,2)
        par1 = noise.noise_vibrations(self.fast_data_path, numbers_array,
                                                tidy_or_shuffle=tos, show=False)
        par2 = noise.convection_noise(self.fast_data_path, tau_vector, 
                                                freq=self.freq, show=True)
        return par1, par2

    def _slow_analysis(self):
        gap = 2
        fl = sorted([os.path.join(self.slow_data_path, file) for file in os.listdir(self.slow_data_path)])
        nfile = len(fl)
        npoints = int(nfile / gap)# - 2
        slist   = []
        for i in range(0, npoints):
            q0 = rd.read_phasemap(fl[i*gap])
            q1 = rd.read_phasemap(fl[i*gap+1])
            diff = zern.removeZernike(q1-q0)
            slist.append(diff.std())
        svec = np.array(slist)
        return svec

    def _fast_acquisition(self):
        n_frames = int(self.freq*3)
        tn = self.interf.capture(n_frames)
        self.interf.produce(tn)
        self.fast_data_path = os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn)
        print("Fast acquisition completed.")

    def _slow_acquisition(self):
        n_frames = self.n_frames
        tn = ts.now()
        self.slow_data_path = os.path.join(fn.OPD_SERIES_ROOT_FOLDER, tn)
        os.mkdir(self.slow_data_path)
        for i in range(n_frames):
            img = self.interf.acquire_phasemap()
            self.interf.save_phasemap(self.slow_data_path, ts.now()+'.fits', img)
            time.sleep(self.delay)
        print("Slow acquisition completed.")

    def __update_interf_settings(self):
        self.freq = self.interf.getFrameRate()
        self.cam_info = self.interf.getCameraSettings()

    def __clear_data_folder(self, path):
        if os.path.exists(path):
            shutil.rmtree(path)

    def __load_monitor_config(self):
        self.interf.loadConfiguration(uc.phasecam_monitorconfig)

    def __load_default_config(self):
        self.interf.loadConfiguration(uc.phasecam_baseconfig)

    def __write_log_message(self, results):
        tn = ts.now()
        long_msg = \
f"""#______________________________________________________________________________
{tn}
[Diagnosis]
Mean : {results[0]*1e9:.2f}nm ; Std = {results[1]*1e9:.2f}nm ; res3 = () ; res4 = ()
[Camera Settings]
FrameRate = {self.freq:.1f} ; Width = {self.cam_info[0]}px ; Height={self.cam_info[1]}px
x-offset = {self.cam_info[2]}px ; y-offset={self.cam_info[3]}px
"""
        short_msg = \
f"""{tn}    {results[0]*1e9:.2f}nm    {results[1]*1e9:.2f}nm    res3    res4    res5"""
        with open(self._llog, 'a', encoding='utf-8') as log:
            log.write(long_msg)
        with open(self._slog, 'a', encoding='utf-8') as log:
            log.write(short_msg)

#______________________________________________________________________________
gui = Gui(
    [ M('plot1') ,   ___  ,   ___   ,   ___  ,  M('plot2') ,   ___  ,   ___  ,   ___  ],
    [    III     ,   III  ,   III   ,   III  ,     III     ,   III  ,   III  ,   III  ],
    [    III     ,   III  ,   III   ,   III  ,     III     ,   III  ,   III  ,   III  ],
    [    III     ,   III  ,   III   ,   III  ,     III     ,   III  ,   III  ,   III  ],
    [ 'RESULTS'  , 'res1' , 'res2'  , 'res3' ,   'res4'    , 'res5' , 'res6' , 'res7' ],
    [     _      ,    _   ,['start'],    _   ,      _      ,['stop'],    _   ,    _   ],
    )

def start(gui, *args):
    gui.res1 = "freq = 20.0"
    plot1(gui)
    plot2(gui)
    return

def stop(gui, *args):
    gui.res1 = 'res1'
    return

def plot1(gui, *args):
    t = np.linspace(0, 2*np.pi, 1000)
    gui.plot1 = np.sin(t)

def plot2(gui, *args):
    t = np.linspace(0, 2*np.pi, 1000)
    gui.plot2 = np.cos(t)

gui.events(
    [ plot1  ,   _   ,   _   ,   _   , plot2 ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   , start ,   _   ,   _   , stop  ,   _   ,   _   ],
    )

# def _run_continuously(self):
#     cease_continuous_run = threading.Event()

#     class ScheduleThread(threading.Thread):
#         @classmethod
#         def run(cls):
#             while not cease_continuous_run.is_set():
#                 schedule.run_pending()
#                 time.sleep(self._timing)

#     continuous_thread = ScheduleThread()
#     continuous_thread.start()
#     return cease_continuous_run

"""
Author(s)
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import time, schedule, threading, os, shutil, numpy as np
from guietta import Gui, _, ___, III, M
from m4 import noise
from m4.ground.timestamp import Timestamp
from m4.configuration import update_folder_paths as ufp
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
        self.template   = np.array([3]) #,11,25,37,51])
        self.delay      = 3
        self.n_frames   = 4
        # self.n          = Noise()
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
        tau_vector = np.arange(1, 40, 2) #??? Test on real data
        out_parameters = noise.convection_noise(self.fast_data_path, tau_vector,
                                                freq=self.freq, show=False)
        return out_parameters

    def _vibrationNoiseMonitoring(self):
        tos = 0
        numbers_array = self.template
        out_parameters = noise.noise_vibrations(self.slow_data_path, numbers_array,
                                                tidy_or_shuffle=tos, show=False)
        return out_parameters

    def _fast_acquisition(self):
        tn = self.interf.capture(self.freq*2)
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

    def __write_log_message(self):
        tn = ts.now()
        long_msg = \
f"""#______________________________________________________________________________
{tn}
[Diagnosis]
res1 = () ; res2 = () ; res3 = () ; res4 = ()
[Camera Settings]
FrameRate = {self.freq:.1f} ; Width = {self.cam_info[0]}px ; Height={self.cam_info[1]}px
x-offset = {self.cam_info[2]}px ; y-offset={self.cam_info[3]}px
"""
        short_msg = \
f"""{tn}    res1    res2    res3    res4    res5"""
        with open(self._llog, 'a', encoding='utf-8') as log:
            log.write(long_msg)
        with open(self._slog, 'a', encoding='utf-8') as log:
            log.write(short_msg)

#______________________________________________________________________________
gui = Gui(
    [     _      ,  _  ,  _  ,  _  ,  _  ,         'RESULTS'       ,  'res1'  , 'res2'  , 'res3' , 'res4' , 'res5' ],
    [ M('plot1') , ___ , ___ , ___ , ___ ,             _           ,     _    ,    _    ,   _    ,    _   ,    _   ],
    [    III     , III , III , III , III ,             _           ,  'res1'  , 'res2'  , 'res3' , 'res4' , 'res5' ],
    [    III     , III , III , III , III ,    'CAMERA SETTINGS'    ,  'freq'  , 'Width' ,'heigth', 'x-off', 'y-off'],
    [    III     , III , III , III , III ,             _           ,  'freq'  , 'Width' ,'heigth', 'x-off', 'y-off'],
    [    III     , III , III , III , III ,             _           ,     _    ,    _    ,   _    ,    _   ,    _   ],
    [ M('plot2') , ___ , ___ , ___ , ___ , 'Acquisition Parameters',     _    ,    _    ,   _    ,    _   ,    _   ],
    [    III     , III , III , III , III ,          'SLOW'         ,'N Frames','__frm__',   _    ,    _   ,    _   ],
    [    III     , III , III , III , III ,             _           ,'Delay'   ,'__dly__',   _    ,    _   ,    _   ],
    [    III     , III , III , III , III ,          'FAST'         ,'Duration','__sec__',   _    ,    _   ,    _   ],
    [    III     , III , III , III , III ,             _           ,     _    ,    _    ,   _    ,    _   ,    _   ],
    [    III     , III , III , III , III ,         ['start']       , ['stop'] ,    _    ,   _    ,    _   ,    _   ],
    )
gui.column_stretch([0,1,1,1,1,1,1,1,1,1,1])
gui.row_stretch([1,1,1,1,1,1,1,0,1,1,1])

def start(gui, *args):
    gui.freq = 20.0

    return

def stop(gui, *args):
    gui.freq = 0.0
    return

def plot1(gui, *args):
    t = np.linspace(1, 10, 1000)
    gui.plot1 = np.sin(t)

def plot2(gui, *args):
    t = np.linspace(1, 10, 1000)
    gui.plot2 = np.cos(t)

gui.events(
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [ plot1  ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [ plot2  ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
    [   _    ,   _   ,   _   ,   _   ,   _   , start ,  stop ,   _   ,   _   ,   _   ,   _   ],
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

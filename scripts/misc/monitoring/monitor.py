"""
Author(s)
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import time, os, shutil, numpy as np #schedule, threading
from guietta import Gui, _, ___, III, M, G, P
from m4 import noise
from m4.configuration import update_folder_paths as ufp
from m4.configuration import userconfig as uc
from m4.ground.timestamp import Timestamp
from m4.ground import zernike as zern, read_data as rd
from m4.utils import osutils as osu
ts = Timestamp()
fn = ufp.folders

class SystemMonitoring():
    """
    Class for the continuous monitoring of the OTT.
    The only public method is the actual monitoring task.

    Methods
    -------
    monitoring() :
        Sets of actions, acquisition and analysis, which monitors the OTT

    How to Use it
    -------------
    >>> from m4.gui.monitor import SystemMonitoring
    >>> # As we have initialized the interferometer device...
    >>> mm = SystemMonitoring(interf)
    >>> mm.monitoring() # single loop operation
    """
    def __init__(self, interferometer):
        """The Constructor"""
        self._llog      = os.path.join(fn.MONITORING_ROOT_FOLDER, 'Monitor_CompleteLog.txt')
        self._slog      = os.path.join(fn.MONITORING_ROOT_FOLDER, 'Monitor_ShortLog.txt')
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

    def rungui(self):
        gui = Gui(
            [ M('plot1')  ,     ___   , ___ ,   G('RESULTS')   , ___ ],
            [    III      ,     III   , III ,       III        , III ],
            [ M('plot2')  ,     ___   , ___ , G('CAMERA INFO') , ___ ],
            [    III      ,     III   , III ,       III        , III ],
            [ M('plot3')  ,     ___   , ___ ,   G('CONTROL')   , ___ ],
            [    III      ,     III   , III ,       III        , III ],
            )
        results_gui = Gui(
            [    'res1'   , 'numero' ,  'nm'  ],
            [    'res2'   , 'numero' ,  'nm'  ],
            [    'res3'   , 'numero' , 'unit1'],
            [    'res4'   , 'numero' , 'unit2'],
            )
        camera_gui = Gui(
            [   'Frequency'  , 'numer1' , 'Hz' ],
            [  'Frame Width' , 'numer2' , 'px' ],
            [ 'Frame Height' , 'numer3' , 'px' ],
            [   'X-Offset'   , 'numer4' , 'px' ],
            [   'Y-Offset'   , 'numer5' , 'px' ],
            )
        control_gui = Gui(
            [ P('progress') ,    ___     ,    ___   ],
            [   ['start']   ,  ['stop']  , ['close']],
            )
        gui.RESULTS     = results_gui
        gui.CAMERAINFO  = camera_gui
        gui.CONTROL     = control_gui
        def show_results(gui, *args):
            gui.res1 = f"Slow mean rms: {self.slow_results[0]*1e9:.1f}nm"
            gui.res2 = f"Slow std: {self.slow_results[1]*1e9:.1f}nm"
            gui.res3 = "OK"
            gui.res4 = "Ok"
            gui.res5 = "Ok"
            gui.res6 = "Ok"
            gui.Frequency  = f"Frequency:  {self.freq:.1f}Hz"
            gui.FrameWidth = f"Frame Width: {self.cam_info[0]:d}px"
            gui.FrameHeight= f"Frame Height: {self.cam_info[1]:d}px"
            gui.XOffset    = f"X-Offset: {self.cam_info[2]:d}px"
            gui.YOffset    = f"Y-Offset: {self.cam_info[3]:d}px"
            plot1(gui)
            plot2(gui)
            plot3(gui)
        def start(gui, *args):
            gui.execute_in_background(self.monitoring,\
                        args=(lambda x: gui.widgets['progress'].setValue(x),),\
                        callback=show_results)
        def stop(gui, *args):
            pass
        def close(gui, *args):
            gui.close()
        def plot1(gui, *args):
            ax = gui.plot1.ax
            ax.clear()
            ax.set_title('RMS')
            ax.set_xlabel('Frames per Template')
            ax.set_ylabel('Root Mean Square [nm]')
            ax.plot(self.fast_results['vn']['ntemp'], self.fast_results['vn']['rms'], '-o', c='black')
            ax.grid()
            ax.figure.canvas.draw()
        def plot2(gui, *args):
            ax = gui.plot2.ax
            ax.clear()
            ax.set_title('Tip-Tilt Quadratic Residual')
            ax.set_xlabel('Frames per Template')
            ax.set_ylabel('Tip - Tilt [nm]')
            ax.plot(self.fast_results['vn']['ntemp'], self.fast_results['vn']['tt'], '-o', c='black')
            ax.grid()
            ax.figure.canvas.draw()
        def plot3(gui, *args):
            ax = gui.plot3.ax
            ax.clear()
            ax.set_title('Decorrelation Noise Fit')
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Root Mean Square [nm]')
            ax.plot(self.fast_results['cn']['x'], self.fast_results['cn']['rms']*1e9, '-o',\
                    c='black', label='Measurements')
            x = self.fast_results['cn']['x']
            pp = self.fast_results['cn']['pp']
            ax.plot([x[0], x[-1]],[pp, pp], "--", linewidth=3, color='red', \
                    label=f"{pp:.2f} [nm]")
            ax.plot(x, self.fast_results['cn']['fit'])
            ax.grid()
            ax.legend()
            ax.figure.canvas.draw()

        gui.events(
            [ plot1  ,   _   ,   _   ,   _   ,   _   ],
            [   _    ,   _   ,   _   ,   _   ,   _   ],
            [ plot2  ,   _   ,   _   ,   _   ,   _   ],
            [   _    ,   _   ,   _   ,   _   ,   _   ],
            [ plot3  ,   _   ,   _   ,   _   ,   _   ],
            [   _    ,   _   ,   _   ,   _   ,   _   ],
            )
        control_gui.events(
            [   _  ,   _  ,   _   ],
            [ start, stop , close ],
            )
        gui.title('OTT Monitoring')
        gui.window().resize(600, 800)
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
        self.__update_interf_settings()
        if progress:
            progress(10)
        self._fast_acquisition()
        if progress:
            progress(30)
        self._slow_acquisition()
        if progress:
            progress(50)
        self._slow_analysis()
        self._fast_analysis()
        if progress:
            progress(80)
        self.__write_log_message()
        self.__clear_data_folder(self.slow_data_path)
        self.__clear_data_folder(self.fast_data_path)
        if progress:
            progress(100)
#        self.__load_default_config() # da fare alla fine del monitoring, non della funzione

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
        tau_vector = np.arange(1,21,2)
        par1 = noise.noise_vibrations(self.fast_data_path, numbers_array,
                                                tidy_or_shuffle=tos, show=False)
        par2 = noise.convection_noise(self.fast_data_path, tau_vector,
                                                freq=self.freq, show=False)
        keys1 = ['rms', 'tt', 'ntemp', 'ptv']
        keys2 = ['rms', 'tt', 'nmeas', 'pp', 'x', 'fit']
        keys3 = ['vn', 'cn']
        res1 = dict(zip(keys1, par1))
        res2 = dict(zip(keys2, par2))
        self.fast_results = dict(zip(keys3,(res1,res2)))

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
        slist   = []
        for i in range(0, npoints):
            q0 = rd.read_phasemap(fl[i*gap])
            q1 = rd.read_phasemap(fl[i*gap+1])
            diff = zern.removeZernike(q1-q0)
            slist.append(diff.std())
        svec = np.array(slist)
        mean = np.mean(svec)
        std = np.std(svec)
        self.slow_results = [mean, std]

    def _fast_acquisition(self):
        """
        Fast acquisition procedure. Capture of freq*3 images, with subsequent
        produce. The tn path is stored in the 'fast_data_path' variable.
        """
        n_frames = int(self.freq*3)
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
        for i in range(n_frames):
            img = self.interf.acquire_phasemap()
            self.interf.save_phasemap(self.slow_data_path, ts.now()+'.fits', img)
            time.sleep(self.delay)
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

    def __load_default_config(self):
        """
        Updates the interferometer settings with a default configuration file,
        found in 'm4.configuration.userconfig'.
        """
        self.interf.loadConfiguration(uc.phasecam_baseconfig)

    def __write_log_message(self):
        """
        Function which writes both the complete and the short monitoring log files.
        """
        tn = ts.now()
        long_msg = \
f"""#______________________________________________________________________________
{tn}
[Diagnosis]
Mean : {self.slow_results[0]*1e9:.2f}nm ; Std = {self.slow_results[1]*1e9:.2f}nm ; res3 = () ; res4 = ()
[Camera Settings]
FrameRate = {self.freq:.1f} ; Width = {self.cam_info[0]}px ; Height={self.cam_info[1]}px
x-offset = {self.cam_info[2]}px ; y-offset={self.cam_info[3]}px
"""
        short_msg = \
f"""{tn}    {self.slow_results[0]*1e9:.2f}nm    {self.slow_results[1]*1e9:.2f}nm    res3    res4    res5"""
        with open(self._llog, 'a', encoding='utf-8') as log:
            log.write(long_msg)
        with open(self._slog, 'a', encoding='utf-8') as log:
            log.write(short_msg)











#______________________________________________________________________________
# =============================================================================
# Debugging of the GUI and relative applications
# =============================================================================

# gui = Gui(
#     [ M('plot1')  ,     ___   ,     ___     ,  M('plot2')  ,   ___    ,    ___   ,  M('plot3') ,   ___   ,   ___  ],
#     [    III      ,     III   ,     III     ,     III      ,   III    ,    III   ,     III     ,   III   ,   III  ],
#     [    III      ,     III   ,     III     ,     III      ,   III    ,    III   ,     III     ,   III   ,   III  ],
#     [ 'RESULTS'   ,   'res1'  ,   'res2'    ,    'res3'    ,  'res4'  ,  'res5'  ,    'res6'   , 'res7'  , 'res8' ],
#     ['CAMERA INFO','Frequency','Frame Width','Frame Height','X-Offset','Y-Offset',      _      ,    _    ,    _   ],
#     [     _       ,     _     ,  ['start']  ,      _       ,    _     , ['stop'] ,      _      ,    _    ,['close']],
#     )

# def start(gui, *args):
#     init = time.time()
#     mm.monitoring()
#     gui.res1 = f"Slow mean rms: {mm.slow_results[0]*1e9:.1f}nm"
#     gui.res2 = f"Slow std: {mm.slow_results[1]*1e9:.1f}nm"
#     gui.res3 = "OK"
#     gui.res4 = "Ok"
#     gui.res5 = "Ok"
#     gui.res6 = "Ok"
#     gui.res7 = "Ok"
#     gui.Frequency  = f"Frequency:  {mm.freq:.1f}Hz"
#     gui.FrameWidth = f"Frame Width: {mm.cam_info[0]:d}px"
#     gui.FrameHeight= f"Frame Height: {mm.cam_info[1]:d}px"
#     gui.XOffset    = f"X-Offset: {mm.cam_info[2]:d}px"
#     gui.YOffset    = f"Y-Offset: {mm.cam_info[3]:d}px"
#     plot1(gui)
#     plot2(gui)
#     plot3(gui)
#     end = time.time()
#     gui.res8 = f"Elapsed Time: {end-init:.1f}s"

# def stop(gui, *args):
#     pass

# def close(gui, *args):
#     gui.close()

# def plot1(gui, *args):
#     ax = gui.plot1.ax
#     ax.clear()
#     ax.set_title('RMS')
#     ax.set_xlabel('Frames per Template')
#     ax.set_ylabel('Root Mean Square [nm]')
#     ax.plot(mm.fast_results['vn']['ntemp'], mm.fast_results['vn']['rms'], '-o', c='black')
#     ax.grid()
#     ax.figure.canvas.draw()

# def plot2(gui, *args):
#     ax = gui.plot2.ax
#     ax.clear()
#     ax.set_title('Tip-Tilt Quadratic Residual')
#     ax.set_xlabel('Frames per Template')
#     ax.set_ylabel('Tip - Tilt [nm]')
#     ax.plot(mm.fast_results['vn']['ntemp'], mm.fast_results['vn']['tt'], '-o', c='black')
#     ax.grid()
#     ax.figure.canvas.draw()

# def plot3(gui, *args):
#     ax = gui.plot3.ax
#     ax.clear()
#     ax.set_title('Decorrelation Noise Fit')
#     ax.set_xlabel('Time [s]')
#     ax.set_ylabel('Root Mean Square [nm]')
#     ax.plot(mm.fast_results['cn']['x'], mm.fast_results['cn']['rms']*1e9, '-o',\
#             c='black', label='Measurements')
#     x = mm.fast_results['cn']['x']
#     pp = mm.fast_results['cn']['pp']
#     ax.plot([x[0], x[-1]],[pp, pp], "--", linewidth=3, color='red', \
#             label=f"{pp:.2f} [nm]")
#     ax.plot(x, mm.fast_results['cn']['fit'])
#     ax.grid()
#     ax.legend()
#     ax.figure.canvas.draw()

# gui.events(
#     [ plot1  ,   _   ,   _   , plot2 ,   _   ,   _   , plot3 ,   _   ,   _   ],
#     [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
#     [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
#     [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
#     [   _    ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ,   _   ],
#     [   _    ,   _   , start ,   _   ,   _   , stop  ,   _   ,   _   , close ],
#     )

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

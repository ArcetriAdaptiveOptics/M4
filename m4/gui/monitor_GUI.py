'''
Authors
  - C. Selmi: written in 2023
'''
from guietta import Gui, G, _, ___, HB, MA
from m4.analyzers.noise_data_analyzer import Noise
from m4.configuration.ott_parameters import Interferometer
from scipy.optimize import curve_fit
import sys
import time
import numpy as np
from m4.ground import read_data
from m4.ground import zernike


class Runner():
    '''
    Class for creating interactive GUI

    HOW TO USE IT::

        conf = 'myConfGlobalPath.yaml'
        from m4.configuration import start
        ott, interf, dm = start.create_ott(conf)
        from m4.gui import monitor_GUI
        g = monitor_GUI.Runner(ott)
        g.run()
    '''

    def __init__(self):
        '''The constructor
        '''


    #file_path = '/Users/cselmi/Desktop/Arcetri/M4/Data/M4Data/20201215_140200_noise/hdf5'
    def _setUp(self):
        
        def Start_Mean_Burst(gui, *args):
            ''' current alignment value (Zernike: TipTilt focus, coma): media burst e fit modi di allineamento.
            '''
            n = Noise()
            lista = n._createOrdListFromFilePath(np.str(gui.file_path))
            cube_list = []
            for i in range(len(lista)):
                name = 'img_%04d' %i
                print(name)
                image = read_data.read_phasemap(lista[i])
                cube_list.append(image)
            cube = np.ma.dstack(cube_list)
            ima_mean = np.ma.mean(cube, axis=2)
            gui.stdBurst = np.std(ima_mean)
            gui.plot3 = ima_mean
            gui.plot3.set_clim(vmin=ima_mean.min(), vmax=ima_mean.max())
            #gui.set_title('Media burst') #comando sbagliato
            gui.zernike_coeff_array, mat = zernike.zernikeFit(ima_mean,
                                                              np.arange(11)+1)
            
        

        def Start_Differential_TipTilt(gui, *args):
            ''' current differential tip/tilt (at the acquisition frequency): noise.noise_vibrations with fixed template.
            '''
            data_file_path = np.str(gui.file_path)
            data_tt = data_file_path.split('/')[-2]
            n = Noise()
            template_list = [np.array([1,-1,1]), np.array([1,-1,1,-1,1]), np.array([1,-1,1,-1,1,-1,1])]

            tt_list = []
            for temp in template_list:
                tt = n.noise_analysis_from_hdf5_folder(data_file_path, 0,
                                                       temp)
                time.sleep(1)
                tt_list.append(tt)
            rms_medio, quad_medio, n_temp, ptv_medio = n.different_template_analyzer(tt_list)
            #noise.noise_vibrations(data_file_path, np.array([23,25,27]), 0)
            setTipTiltPlot(gui, data_tt, n_temp, rms_medio)   

        def setTipTiltPlot(gui, tt, n_temp, rms_medio):
            ax = gui.plot.ax
            ax.clear()
            ax.set_title(tt)
            
            ax.plot(n_temp, rms_medio * 1e9, '-o')
            ax.set_xlabel('n_temp')
            ax.set_ylabel('rms [nm]')
            ax.figure.canvas.draw()
            #ax.grid()
        
        def Start_Decorrelation_Time(gui, *args):
            ''' current decorrelation time: noise.convection_noise. Output: grafico.
            '''
            data_file_path = np.str(gui.file_path)
            data_tt = data_file_path.split('/')[-2]
            n = Noise()
            
            tau_vector = np.arange(1, 60)
            rms, quad, n_meas = n.analysis_whit_structure_function(data_file_path,
                                                           tau_vector,
                                                           h5_or_fits=None, nzern=None)

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
            ax.plot(x, rms * 1e9, '-o', label='meas')
            ax.set_xlabel('time [s]')
            ax.set_ylabel('rms [nm]')
            ax.plot(x, fit, '-', label='fit')
            ax.plot([x[0], x[-1]], [pp[2], pp[2]], '--r', linewidth=3,
                     label='%.2f [nm]' % pp[2])
            ax.legend()
            ax.figure.canvas.draw()

        def Start_Decorrelation_Noise(gui, *args):
            ''' current decorrelated noise: diff delle immagini dopo tau (prima - ultima). Output: image and RMS
            '''
            pass

        def Start_Convection_Footprint_PSD(gui, *args):
            ''' convection footprint PSD: codice di Marchino che restituisce PSD da immagine del punto 4 (decorrelation noise)
            '''
            pass

        def Start_Convective_Regions(gui, *args):
            ''' convective regions (pixel-wise stdev): cubo delle immagini senza tip/tilt, std nella direzione giusta che restituisce un'immagine.
            '''
            pass

        def Start_Analysis(gui, *args):
            ''' Start all the analysis
            '''
            pass

        gui_analysis = Gui(['Set data file path', '__file_path__', ___, ___, ___, ___, ___, ___],
                           [Start_Mean_Burst, ___, ___, ___, MA('plot3'), ___, ___, ___],
                           [ 'RMS value:',___, ___, 'stdBurst', _, _, _, _],
                           [ 'Zernike:',___, ___, 'zernike_coeff_array', ___, ___, ___, ___],
                           [Start_Differential_TipTilt, ___, ___, ___, Start_Decorrelation_Time, ___, ___, ___],
                           [MA('plot'), ___, ___, ___, MA('plot2'), ___, ___, ___])
        
        self.gui = Gui([G('OTT_Monitor_GUI')])
        self.gui.OTT_Monitor_GUI = gui_analysis
        
    def run(self):
        ''' Run the GUI '''
        self._setUp()
        self.gui.run()

def main():
    from m4.configuration import start
    conf = '/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
    ott, interf, dm = start.create_ott(conf)

    runner = Runner(ott)
    sys.exit(runner.run())

if __name__ == '__main__':
    main()
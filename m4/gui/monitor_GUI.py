'''
Authors
  - C. Selmi: written in 2023
'''
from guietta import Gui, G, _, ___, HB, MA
from m4 import noise
from m4.analyzers.noise_data_analyzer import Noise
import sys
import time
import numpy as np


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
        
        def Start_Alignment_Value(gui, *args):
            ''' current alignment value (Zernike: TipTilt focus, coma): media burst e fit modi di allineamento.
            '''
            pass

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
            pass

        def Start_Decorrelation_Noise(gui, *args):
            ''' current decorrelated noise: diff delle immagini dopo tau (prima - ultima). Output: image and RMS
            '''
            pass

        def Start_Convection_Footprint_PSD(gui, *args):
            ''' convection footprint PSD: codice di Marchino che restituisce PSD da immagine del punto 4
            '''
            pass

        def Start_Convective_Reions(gui, *args):
            ''' convective regions (pixel-wise stdev): cubo delle immagini senza tip/tilt, std nella direzione giusta che restituisce un'immagine.
            '''
            pass

        def Start_Analysis(gui, *args):
            ''' Start all the analysis
            '''
            pass

        gui_analysis = Gui(['Set data file path', '__file_path__', _, _, _],
                              [Start_Differential_TipTilt, _, _, _, _],
                              [MA('plot'), _, _, _, _])
        
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
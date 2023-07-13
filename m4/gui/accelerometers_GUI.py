'''
Authors
  - C. Selmi: written in 2023
'''
import numpy as np
from m4.ground.package_data import data_root_dir
from m4.analyzers.accelerometers_data_analyzer import AccelerometersDataAnalyzer
from guietta import Gui, G, _, ___, HB
import sys

class Runner():
    '''
    Class for creating interactive GUI

    HOW TO USE IT::

        conf = 'myConfGlobalPath.yaml'
        from m4.configuration import start
        ott, interf, dm = start.create_ott(conf)
        from m4.gui import accelerometers_GUI
        g = accelerometers_GUI.Runner(ott)
        g.run()
    '''

    def __init__(self, ott):
        '''The constructor
        ott: object
            tower object
        '''
        self.ott = ott

    def _setUp(self):

        def setPlot(gui):
            '''
            image: numpy array
                image of the tower
            '''
            gui.heart_empty_30.beat()

        def Start_Acquisition(gui, *args):
            ''' Function to reload tower data and image
            '''
            recording_seconds = np.int(gui.rec_sec) or 5.
            if self.ott:
                gui.tn = self.ott.accelerometers.acquireData(recording_seconds)
        
        def Start_Analysis(gui, *args):
            track_num = np.str(gui.track_num)
            an = AccelerometersDataAnalyzer(track_num)
            an.readAndShow()

        gui_acquisition = Gui(['Set recording time', '__rec_sec__', _, '[s]'],
                              [Start_Acquisition, _, _, _],
                              [_, _, 'Tracking Number generated:', 'tn'],
                              [ _, _, _, HB('heart_empty_30.png', 'heart_full_30.png')], images_dir=data_root_dir())
        
        gui_analysis = Gui(['Set tracking number', '__track_num__'],
                              [Start_Analysis, _])

        self.gui = Gui([ G('Acc_acquisitions'), G('Acc_analysis')])
        self.gui.Acc_acquisitions = gui_acquisition
        self.gui.Acc_analysis = gui_analysis

        
    def run(self):
        ''' Run the GUI '''
        self._setUp()
        self.gui.run()

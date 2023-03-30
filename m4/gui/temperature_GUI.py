'''
Authors
  - C. Selmi: written in 2023
'''
import numpy as np
from m4.ground.package_data import data_root_dir
from guietta import Gui, G, MA, _, ___, III, HB, PGI

class Runner():
    '''
    Class for creating interactive GUI

    HOW TO USE IT::

        conf = 'myConfGlobalPath.yaml'
        from m4.configuration import start
        ott, interf, dm = start.create_ott(conf)
        from m4.gui import temperature_GUI
        g = temperature_GUI.Runner(ott)
        g.run()
    '''

    def __init__(self, ott):
        '''The constructor
        ott: object
            tower object
        '''
        self.ott = ott

    def _setUp(self):

        def setPlot(gui, image):
            '''
            image: numpy array
                image of the tower
            '''
            gui.plot = image
            gui.plot.set_clim(vmin=image.min(), vmax=image.max())
            gui.heart_empty_30.beat()
            
        
        gui_image = Gui([ PGI('heart_full_30'), ___    , ___, ___ ],
                        [   _ ,  HB('heart_empty_30.png', 'heart_full_30.png'), _, _], images_dir=data_root_dir())
        #gui_image.plot.set_title('Temperatures')
        
        self.gui = Gui([ G('OTT_temperatures')])
        self.gui.OTT_temperature = gui_image
        self.gui.heart_full_30 = np.arange(100).reshape((10,10))
        
    def run(self):
        ''' Run the GUI '''
        self._setUp()
        self.gui.run()
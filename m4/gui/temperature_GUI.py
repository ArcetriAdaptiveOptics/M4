'''
Authors
  - C. Selmi: written in 2023
'''
import numpy as np
from m4.ground.package_data import data_root_dir
from guietta import Gui, G, MA, _, ___, III, HB, PGI
import sys

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
            gui.heart_empty_30.beat()
        
        def getstatus(gui):
            ''' Function to reload tower data and image
            '''
            temp = self.ott.temperature.getTemperature()
            gui.t0 = temp[0]
            gui.t1 = temp[1]
            gui.t2 = temp[2]
            gui.t3 = temp[3]
            gui.t4 = temp[4]
            gui.t5 = temp[5]
            gui.t6 = temp[6]
            gui.t7 = temp[7]
            gui.t8 = temp[8]
            gui.t9 = temp[9]
            gui.t10 = temp[10]
            gui.t11 = temp[11]
            gui.t12 = temp[12]
            gui.t13 = temp[13]
            gui.t14 = temp[14]
            gui.t15 = temp[15]
            gui.t16 = temp[16]
            gui.t17 = temp[17]
            gui.t18 = temp[18]
            gui.t19 = temp[19]
            gui.t20 = temp[20]
            gui.t21 = temp[21]
            gui.t22 = temp[22]
            gui.t23 = temp[23]


        gui_control = Gui(['T0, T1', 't0', 't1', '[deg]'],
                        ['T2, T3', 't2', 't3', '[deg]' ],
                        ['T4, T5', 't4', 't5', '[deg]' ],
                        ['T6, T7', 't6', 't7', '[deg]' ],
                        ['T8, T9', 't8', 't9', '[deg]' ],
                        ['T10, T11', 't10', 't11', '[deg]' ],
                        ['T12, T13', 't12', 't13', '[deg]' ],
                        ['T14, T15', 't14', 't15', '[deg]' ],
                        ['T16, T17', 't16', 't17', '[deg]' ],
                        ['T18, T19', 't18', 't19', '[deg]' ],
                        ['T20, T21', 't20', 't21', '[deg]' ],
                        ['T22, T23', 't22', 't23', '[deg]' ],
                        [ _, _, _, HB('heart_empty_30.png', 'heart_full_30.png')], images_dir=data_root_dir())

        gui_image = Gui(['OttTemps.png', _], images_dir=data_root_dir())

        gui_control.timer_start(getstatus, 1)
        
        self.gui = Gui([ G('OTT_temperature'), G('OTT_view')])
        self.gui.OTT_temperature = gui_control
        self.gui.OTT_view = gui_image

        
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


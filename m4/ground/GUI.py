'''
Authors
  - C. Selmi: written in 2021
'''

import sys
import time
import numpy as np
from m4.ott_sim.ott_images import OttImages
from guietta import Gui, G, MA, _, ___, III, HB
from m4.ground.package_data import data_root_dir
#from guietta import Empty, Exceptions



class Runner():
    '''
    Class for creating interactive GUI

    HOW TO USE IT::

        conf = 'myConfGlobalPath.yaml'
        from m4.configuration import start
        ott, interf = start.create_ott(conf)
        from m4.ground import GUI
        g = GUI.Runner(ott)
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

        def getstatus(gui):
            ''' Function to reload tower data and image
            '''
            if self.ott:
                gui.parpos = self.ott.parabola.getPosition()
                gui.rmpos = self.ott.referenceMirror.getPosition()
                gui.m4pos = self.ott.m4Exapode.getPosition()
                gui.pslider = self.ott.parabolaSlider.getPosition()
                gui.anglepos = self.ott.angleRotator.getPosition()
                gui.rslider = self.ott.referenceMirrorSlider.getPosition()
            ottIma = OttImages(self.ott)
            image = ottIma.ott_view()
            setPlot(gui, image)

        def Set_Parabola(gui, *args):
            ''' Function to move the parabola'''
            vec = self.ott.parabola.getPosition()
            vec[2] = np.float(gui.parpist or vec[2])
            vec[3] = np.float(gui.partip or vec[3])
            vec[4] = np.float(gui.partilt or vec[4])
            if self.ott:
                self.ott.parabola.setPosition(vec)

        def Set_RefMirror(gui, *args):
            ''' Function to move the reference mirror '''
            vec = self.ott.referenceMirror.getPosition()
            vec[2] = np.float(gui.rmpist or vec[2])
            vec[3] = np.float(gui.rmtip or vec[3])
            vec[4] = np.float(gui.rmtilt or vec[4])
            if self.ott:
                self.ott.referenceMirror.setPosition(vec)

        def Set_M4(gui, *args):
            ''' Function to move the exapode '''
            vec = self.ott.m4Exapode.getPosition()
            vec[2] = np.float(gui.m4pist or vec[2])
            vec[3] = np.float(gui.m4tip or vec[3])
            vec[4] = np.float(gui.m4tilt or vec[4])
            if self.ott:
                self.ott.m4Exapode.setPosition(vec)

        def Set_ReferenceMirrorSlider(gui, *args):
            ''' Function to move the reference mirror slider '''
            pos = np.int(gui.rmslider) or self.ott.referenceMirrorSlider.getPosition()
            if self.ott:
                self.ott.referenceMirrorSlider.setPosition(pos)

        def Set_AngleRotator(gui, *args):
            ''' Function to rotate the tower angle'''
            pos = np.int(gui.angle) or self.ott.angleRotator.getPosition()
            if self.ott:
                self.ott.angleRotator.setPosition(pos)

        def Set_ParabolaSlider(gui, *args):
            ''' Function to move the parabola slider '''
            pos = np.int(gui.parslider) or self.ott.parabolaSlider.getPosition()
            if self.ott:
                self.ott.parabolaSlider.setPosition(pos)


        gui_image = Gui([ MA('plot')            , ___    , ___, ___ ],
                        ['Par position:'        ,'parpos', ___, '[-, -, mm, arcsec, arcsec, -]'],
                        ['Rm position:'        , 'rmpos', ___, '[-, -, mm, arcsec, arcsec, -]'],
                        ['M4 position:'        , 'm4pos', ___, 'mm'],
                        ['Par slider position:', 'pslider', ___,'mm'],
                        ['Rm slider position:', 'rslider', ___,'mm'],
                        ['Ang rot position:', 'anglepos', ___, 'deg'],
                        [   _               ,    _      ,  HB('heart_empty_30.png',
                                                              'heart_full_30.png'), _], images_dir=data_root_dir())

        control_gui = Gui(['New Par position', '0', '0', '__parpist__', '__partip__', '__partilt__', '0', '[-, -, mm, arcsec, arcsec, -]'],
                          [Set_Parabola, ___, ___, ___, ___, ___, ___, ___],
                          ['New Rm position', '0', '0', '__rmpist__', '__rmtip__', '__rmtilt__', '0', '[-, -, mm, arcsec, arcsec, -]'],
                          [Set_RefMirror, ___, ___, ___, ___, ___, ___, ___],
                          ['New M4 position', '0', '0', '__m4pist__', '__m4tip__', '__m4tilt__', '0', 'mm'],
                          [Set_M4, ___, ___, ___, ___, ___, ___, ___],
                          ['New Par Slider position', '__parslider__', _, _, _, _, _, 'mm'],
                          [Set_ParabolaSlider, ___, ___, ___, ___, ___, ___, ___],
                          ['New Rm Slider position', '__rmslider__', _, _, _, _, _, 'mm'],
                          [Set_ReferenceMirrorSlider, ___, ___, ___, ___, ___, ___, ___],
                          ['New Angle Rot position', '__angle__', _, _, _, _, _, 'deg'],
                          [Set_AngleRotator, ___, ___, ___, ___, ___, ___, ___]) #exceptions=Exceptions.OFF)

        ottIma = OttImages(self.ott)
        image = ottIma.ott_view()

        gui_image.plot = image
        gui_image.plot.set_title('OTT geometry')
        gui_image.plot.colorbar()
        #gui_image.plot.text(30, 50, 'ciao')

        gui_image.timer_start(getstatus, 1)

        self.gui = Gui(
             [ G('OTT') , G('Control') ]
             )
        self.gui_control = control_gui
        self.gui_image = gui_image
        self.gui.OTT = gui_image
        self.gui.Control = control_gui


    def runImage(self):
        self._setUp()
        self.gui_image.run()

    def run(self):
        ''' Run the GUI '''
        self._setUp()
        self.gui.run()


def main():
    from m4.configuration import start
    conf = '/home/labot/Desktop/labotConfig.yaml' #modificare all'occorrenza
    ott, interf = start.create_ott(conf)

    runner = Runner(ott)
    sys.exit(runner.run())

if __name__ == '__main__':
    main()

'''
Authors
  - C. Selmi: written in 2021
'''

import sys
import numpy as np

from m4.ott_sim.ott_images import OttImages

from guietta import Gui, G, MA, _, ___, III
from guietta import Empty, Exceptions



class Runner():

    def __init__(self, ott, interf):
        self.ott = ott
        self.interf = interf

    def _setUp(self):

        def setPlot(gui, image):
            gui.plot = image

        def getstatus(gui):
            if self.ott:
                gui.parpos = self.ott.parabola.getPosition()
                gui.rmpos = self.ott.referenceMirror.getPosition()
                gui.m4pos = self.ott.m4.getPosition()
                gui.pslider = self.ott.parabolaSlider.getPosition()
                gui.anglepos = self.ott.angleRotator.getPosition()
                gui.rslider = self.ott.referenceMirrorSlider.getPosition()
            image = self.interf.acquire_phasemap()
            setPlot(gui, image)

        def movepar(gui, *args):
            pist = np.float(gui.parpist or 0) #mettere pos precedente
            tip = np.float(gui.partip or 0)
            tilt = np.float(gui.partilt or 0)
            vec = np.array([0, 0, pist, tip, tilt, 0])
            if self.ott:
                self.ott.parabola.setPosition(vec)

        def moverm(gui, *args):
            pist = np.array(gui.rmpist or 0)
            tip = np.array(gui.rmtip or 0)
            tilt = np.array(gui.rmtilt or 0)
            vec = np.array([0, 0, pist, tip, tilt, 0])
            if self.ott:
                self.ott.referenceMirror.setPosition(vec)

        def movem4(gui, *args):
            pist = np.array(gui.m4pist or 0)
            tip = np.array(gui.m4tip or 0)
            tilt = np.array(gui.m4tilt or 0)
            vec = np.array([0, 0, pist, tip, tilt, 0])
            if self.ott:
                self.ott.m4.setPosition(vec)

        def moverslider(gui, *args):
            pos = np.int(gui.rmslider)
            if self.ott:
                self.ott.referenceMirrorSlider.setPosition(pos)

        def moveangle(gui, *args):
            pos = np.int(gui.angle)
            if self.ott:
                self.ott.angleRotator.setPosition(pos)

        def movepslider(gui, *args):
            pos = np.int(gui.parslider)
            if self.ott:
                self.ott.parabolaSlider.setPosition(pos)


        image = self.interf.acquire_phasemap()
        #getstatus()
        gui_image = Gui([ MA('plot'), ___, ___, ___, ___, ___, ___, ___ ],
                        [ III, III, III, III, III, III, III, III ],
                        [ III, III, III, III, III, III, III, III ],
                        [ III, III, III, III, III, III, III, III ],
                        [ III, III, III, III, III, III, III, III ],
                        [ III, III, III, III, III, III, III, III ],
                        ['Par position:','parpos', ___, ___, ___, ___, ___, 'mm'],
                        ['Rm position:', 'rmpos', ___, ___, ___, ___, ___, 'mm'],
                        ['M4 position:', 'm4pos', ___, ___, ___, ___, ___, 'mm'],
                        ['Par slider position:', 'pslider', ___, ___, ___, ___, ___,'mm'],
                        ['Rm slider position:', 'rslider', ___, ___, ___, ___, ___,'mm'],
                        ['Ang rot position:', 'anglepos', ___, ___, ___, ___, ___, 'deg'])
        control_gui = Gui(['New Par position', '0', '0', '__parpist__', '__partip__', '__partilt__', '0', 'mm'],
                          [['Set_Parabola'], ___, ___, ___, ___, ___, ___, ___],
                          ['New Rm position', '0', '0', '__rmpist__', '__rmtip__', '__rmtilt__', '0', 'mm'],
                          [['Set_RefMirror'], ___, ___, ___, ___, ___, ___, ___],
                          ['New M4 position', '0', '0', '__m4pist__', '__m4tip__', '__m4tilt__', '0', 'mm'],
                          [['Set_M4'], ___, ___, ___, ___, ___, ___, ___],
                          ['New Par Slider position', '__parslider__', _, _, _, _, _, 'mm'],
                          [['Set_ParabolaSlider'], ___, ___, ___, ___, ___, ___, ___],
                          ['New Rm Slider position', '__rmslider__', _, _, _, _, _, 'mm'],
                          [['Set_ReferenceMirrorSlider'], ___, ___, ___, ___, ___, ___, ___],
                          ['New Angle Rot position', '__angle__', _, _, _, _, _, 'deg'],
                          [['Set_AngleRotator'], ___, ___, ___, ___, ___, ___, ___]) #exceptions=Exceptions.OFF)

        setPlot(gui_image, image)


        control_gui.Set_Parabola = movepar
        control_gui.Set_RefMirror = moverm
        control_gui.Set_M4 = movem4
        control_gui.Set_ParabolaSlider = movepslider
        control_gui.Set_ReferenceMirrorSlider = moverslider
        control_gui.Set_AngleRotator = moveangle

        gui_image.timer_start(getstatus, 1)

        self.gui = Gui(
             [ G('OTT') , G('Control') ]
             )

        self.gui.Control = control_gui
        self.gui.OTT = gui_image

#         while 0:
#             try:
#                 name, event = gui_image.get(timeout=1)
# 
#             except Empty:
#                 oi = OttImages(self.ott)
#                 smap1, smask = oi.ott_smap()
#                 setPlot(gui_image, image)
#                 continue
# 
#             if name is None:
#                 break
#             if name in ['viridis', 'hot']:
#                 cmap = name

    def run(self):
        self._setUp()
        self.gui.run()

#non funziona
if __name__ == '__main__':
    from m4.configuration import start
    conf = '/Users/rm/eclipse-workspace/M4/m4/configuration/myConfig.yaml'
    ott, interf = start.create_ott(conf)                                    

    runner = Runner(ott, interf)
    sys.exit(runner.run())


'''
Authors
  - C. Selmi: written in 2020
'''

import numpy as np
import matplotlib.pyplot as plt
from m4.ott_sim.ott_images import OttImages

from guietta import Gui, M, _, ___, III, VS
from guietta import Ax
from guietta import Empty, Exceptions
from m4.configuration.create_ott import OTT


def main(ott, interf):
#     a = M('plot')
#     print(isinstance(a,QWidget))
    image = np.arange(512**2).reshape((512, 512))
    ott.referenceMirror.setPosition(np.array([0.e+00, 0.e+00, 0.e+00, 1.e-07, -4.e-07, 0.e+00]))
    ott.parabola.setPosition(np.array([0, 0, 9.9999997e-06, 1.0836526e-07, 2.2898718e-09, 0]))
    ott.parabolaSlider.setPosition(0.75)
    ott.angleRotator.setPosition(90.)
    ott.referenceMirrorSlider.setPosition(0.6)
    oi = OttImages(ott)
    smap1, smask = oi.ott_smap()

    cmap = 'viridis'
    gui = Gui([ ['hot'], ___, ['viridis'], ___ ],
               [M('plot'), M('plot2'), M('plot3'), M('plot4')])
    setPlot(gui, image)

    while True:
        try:
            name, event = gui.get(timeout=1)

        except Empty:
            setPlot(gui, image**2, cmap)
            setPlot2(gui, oi, smap1, smask, cmap)
            setPlot3(gui, oi, cmap)
            setPlot4(gui, oi, smap1, cmap)
            continue

        if name is None:
            break
        if name in ['viridis', 'hot']:
            cmap = name



def setPlot(gui, image, cmap='viridis'):
    with Ax(gui.plot) as ax:
        ax.clear()
        ax.set_title('Titolo')
        im = ax.imshow(image, origin='lower', cmap=cmap)
        #ax.figure.colorbar(im, ax=ax)

def setPlot2(gui, oi, smap1, smask, cmap='viridis'):
    image = oi.pwrap(smap1, smask)
    with Ax(gui.plot2) as ax:
        ax.clear()
        ax.set_title('Titolo')
        im = ax.imshow(image, origin='lower', cmap=cmap)
        #ax.figure.colorbar(im, ax=ax)

def setPlot3(gui, oi, cmap='viridis'):
    with Ax(gui.plot3) as ax:
        ax.clear()
        ax.set_title('Titolo')
        im = ax.imshow(oi.ott_view(), origin='lower', cmap=cmap)
        #ax.figure.colorbar(im, ax=ax)

def setPlot4(gui, oi, smap1, cmap='hot'):
    with Ax(gui.plot4) as ax:
        ax.clear()
        ax.set_title('Titolo')
        im = ax.imshow(smap1, origin='lower', cmap=cmap)
        #ax.figure.colorbar(im, ax=ax

class Runnes():

    def __init__(self, ott, interf):
        self.ott = ott
        self.interf = interf

    def _setUp(self):

        def getstatus(gui):
            if self.ott:
                gui.parpos = self.ott.parabola.getPosition()
                gui.rmpos = self.ott.referenceMirror.getPosition()
                gui.m4pos = self.ott.m4.getPosition()
                gui.pslider = self.ott.parabolaSlider.getPosition()
                gui.anglepos = self.ott.angleRotator.getPosition()
                gui.rslider = self.ott.referenceMirrorSlider.getPosition()

        def movepar(gui):
            vec = np.array(gui.parpos)
            if self.ott:
                self.ott.parabola.setPosition(vec)

        def moverm(gui):
            vec = np.array(gui.rmpos)
            if self.ott:
                self.ott.referenceMirror.setPosition(vec)

        def movem4(gui):
            vec = np.array(gui.m4pos)
            if self.ott:
                self.ott.m4.setPosition(vec)

        def moverslider(gui):
            pos = np.int(gui.rslider)
            if self.ott:
                self.ott.referenceMirrorSlider.setPosition(pos)

        def moveangle(gui):
            pos = np.int(gui.angle)
            if self.ott:
                self.ott.angleRotator.setPosition(pos)

        def movepslider(gui):
            pos = np.int(gui.pslider)
            if self.ott:
                self.ott.parabolaSlider.setPosition(pos)


        image = self.interf.acquire_phasemap()
        cmap = 'viridis'
        #getstatus()
        gui = Gui([ ['hot'], ___, ['viridis'], ___ , _, 'Par position:', 'parpos', 'mm?'],
                   [M('plot2'), ___, ___, ___,['Set_Parabola'], ___, '__parpos__', 'mm?'],
                   [ III, III, III, III, _, 'Rm position:', 'rmpos', 'mm?'],
                   [ III, III, III, III, ['Set_RefMirror'], ___, '__rmpos__', 'mm?'],
                   [ III, III, III, III, _, 'M4 position:', 'm4pos', 'mm?'],
                   [ III, III, III, III, ['Set_M4'], ___, '__m4pos__', 'mm?'],
                   [ III, III, III, III, _, 'Par slider position:', 'pslider', 'mm?'],
                   [ III, III, III, III, ['Set_ParabolaSlider'], ___, '__pslider__', 'mm'],
                   [ III, III, III, III, _, 'Rm slider position:', 'rslider', 'mm?'],
                   [ III, III, III, III, ['Set_ReferenceMirrorSlider'], ___, '__rslider__', 'mm'],
                   [ III, III, III, III, _, 'Ang rot position:', 'anglepos', 'mm?'],
                   [ III, III, III, III, ['Set_AngleRotator'], ___, '__angle__', 'deg']) #exceptions=Exceptions.OFF)
        oi = OttImages(self.ott)
        smap1, smask = oi.ott_smap()
        setPlot2(gui, oi, smap1, smask, cmap)

        gui.Set_Parabola = movepar
        gui.Set_RefMirror = moverm
        gui.Set_M4 = movem4
        gui.Set_ParabolaSlider = movepslider
        gui.Set_ReferenceMirrorSlider = moverslider
        gui.Set_AngleRotator = moveangle

        while True:
            try:
                name, event = gui.get(timeout=1)

            except Empty:
                oi = OttImages(self.ott)
                smap1, smask = oi.ott_smap()
                setPlot2(gui, oi, smap1, smask, cmap)
                continue

            if name is None:
                break
            if name in ['viridis', 'hot']:
                cmap = name
        self.gui = gui

    def run(self):
        self._setUp()
        self.gui.run()

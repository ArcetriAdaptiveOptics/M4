'''
Authors
  - C. Selmi: written in 2020
'''

import numpy as np
import matplotlib.pyplot as plt
from m4.ott_sim.ott_images import OttImages

from guietta import Gui, M, _, ___, III, VS
from guietta import Ax
from guietta import Empty


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


# class GuiCreator():
#
#     def __init__(self):
#         """The constructor """
#         self._gui = Gui(
#               [  M('plot') , ___ ,  ___ , ___ ],
#               [     III    , III ,  III ,     III      ],
#               [     III    , III ,  III ,     III      ],
#               [     III    , III ,  III ,     III],
#              )
#
#     def create_one_image(self, image):
#         with Ax(self._gui.plot) as ax:
#             #ax.clear()
#             ax.set_title('Titolo')
#             im = ax.imshow(image, origin='lower', cmap='Reds')
#             ax.figure.colorbar(im, ax=ax)
#
#     def plot_image(self, image):
#         self._gui.events(
#             [  _            ,  _ , _ ,   self.create_one_image ],
#             [  _            ,  _ , _ ,   _          ],
#             [  _            ,  _ , _ ,   _          ], )
#         self.create_one_image(image)
#         self._gui.run()

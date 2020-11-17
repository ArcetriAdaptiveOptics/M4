'''
Autors
  - C. Selmi: written in 2019
'''

import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive
from guietta import Gui, M, _, ___, III, VS
from guietta import Ax

class GuiCreator():

    def __init__(self):
        """The constructor """
        self._gui = Gui(
              [  M('plot') , ___ ,  ___ , ___ ],
              [     III    , III ,  III ,     III      ],
              [     III    , III ,  III ,     III      ],
              [     III    , III ,  III ,     III],
             )


    def create_one_image(self, image):
        with Ax(self._gui.plot) as ax:
            ax.clear()
            ax.set_title('Titolo')
            im = ax.imshow(image, origin='lower', cmap='Reds')
            ax.figure.colorbar(im, ax=ax)

    def plot_image(self, image):
        self._gui.events(
            [  _            ,  _ , _ ,   self.create_one_image ],
            [  _            ,  _ , _ ,   _          ],
            [  _            ,  _ , _ ,   _          ], )
        self.create_one_image(image)
        self._gui.run()

    def test(self, color='Blues'):
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        ax.set_title('Titolo')
        im = ax.imshow(self._image, origin='lower', cmap=color)
        ax.figure.colorbar(im, ax=ax)

    def piùColor(self, image, color='Blues'):
        self._image = image
        interactive(self.test, {'manual': True}, color=['Blues', 'Reds', 'Greens'])

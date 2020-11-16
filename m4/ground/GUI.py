'''
Autors
  - C. Selmi: written in 2019
'''

import numpy as npax 
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

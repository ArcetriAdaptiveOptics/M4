'''
Autors
  - C. Selmi: written in 2019
'''

import numpy as np
from guietta import Gui, M, _, ___, III, VS
from guietta import Ax

class GuiCreator():

    def __init__(self):
        """The constructor """
        self._gui = Gui(
              [  M('plot') , ___ ,  ___ , VS('slider') ],
              [     III    , III ,  III ,     III      ],
              [     III    , III ,  III ,     III      ],
              [     III    , III ,  III ,  '^^^ Move the slider'  ],
             )


    def plot_image(self, image):
        with Ax(self._gui.plot) as ax:
            ax.set_title('Titolo')
            ax.imshow(image, origin='lower')
            #ax.colorbar()

    def replot(self, image):
        self.plot_image(image)
        self._gui.events(
            [  _            ,  _ , _ ,   self.plot_image ],
            [  _            ,  _ , _ ,   _          ],
            [  _            ,  _ , _ ,   _          ], )
        self._gui.run()

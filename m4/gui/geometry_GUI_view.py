"""
Authors
  - C. Selmi: written in 2021
"""
import os
import sys
import time
import numpy as np
from m4.ott_sim.ott_images import OttImages
from guietta import Gui, G, MA, _, ___, III, HB
from m4.ground.package_data import data_root_dir
from m4.devices.opt_beam import Parabola, ReferenceMirror, AngleRotator
conf = os.environ['PYOTTCONF']

# from guietta import Empty, Exceptions


class Runner:
    """
    Class for creating interactive GUI

    HOW TO USE IT::

        conf = 'myConfGlobalPath.yaml'
        from m4.configuration import start
        ott, interf, dm = start.create_ott(conf)
        from m4.gui import geometry_GUI
        g = geometry_GUI.Runner(ott)
        g.run()
    """

    def __init__(self, ott):
        """The constructor
        ott: object
            tower object
        """
        self.ott    = ott
        self.truss  = Parabola(ott, conf)
        self.rm     = ReferenceMirror(ott, conf)
        self.angrot = AngleRotator(ott, conf)

    def _setUp(self):

        def setPlot(gui, image):
            """
            image: numpy array
                image of the tower
            """
            gui.plot = image
            gui.plot.set_clim(vmin=image.min(), vmax=image.max())
            gui.heart_empty_30.beat()

        def getstatus(gui):
            """Function to reload tower data and image"""
            if self.truss and self.rm and self.angrot:
                gui.parpos = self.ott.parabola.getPosition()
                gui.rmpos = self.ott.referenceMirror.getPosition()
                gui.m4pos = self.ott.m4Exapode.getPosition()
                #gui.pslider = self.ott.parabolaSlider.getPosition()+1
                gui.pslider = self.truss.trussGetPosition()
                #gui.anglepos = self.ott.angleRotator.getPosition()
                gui.anglepos = self.angrot.getPosition()
                #gui.rslider = self.ott.referenceMirrorSlider.getPosition()+1
                gui.rslider = self.rm.rmGetPosition()
                #gui.psliderM4 = self.ott.parabolaSlider.getPosition()+self.paroffset*1000
                #gui.rsliderM4 = self.ott.referenceMirrorSlider+self.rmoffset*1000
            ottIma = OttImages(self.ott)
            image = ottIma.ott_view()
            setPlot(gui, image)

        def Set_Parabola(gui, *args):
            """Function to move the parabola"""
            vec = self.ott.parabola.getPosition()
            vec[2] = float(gui.parpist or vec[2])
            vec[3] = float(gui.partip or vec[3])
            vec[4] = float(gui.partilt or vec[4])
            if self.ott:
                self.ott.parabola.setPosition(vec)

        def Set_RefMirror(gui, *args):
            """Function to move the reference mirror"""
            vec = self.ott.referenceMirror.getPosition()
            vec[2] = float(gui.rmpist or vec[2])
            vec[3] = float(gui.rmtip or vec[3])
            vec[4] = float(gui.rmtilt or vec[4])
            if self.ott:
                self.ott.referenceMirror.setPosition(vec)

        def Set_M4(gui, *args):
            """Function to move the exapode"""
            vec = self.ott.m4Exapode.getPosition()
            vec[2] = float(gui.m4pist or vec[2])
            vec[3] = float(gui.m4tip or vec[3])
            vec[4] = float(gui.m4tilt or vec[4])
            if self.ott:
                self.ott.m4Exapode.setPosition(vec)

        def Set_ReferenceMirrorSlider(gui, *args):
            """Function to move the reference mirror slider"""
            pos = self.rm.rmGetPosition()
            if gui.rslider:
                pos = float(gui.rslider)
            if self.rm:
                self.rm.moveRmTo(pos)

        def Set_AngleRotator(gui, *args):
            """Function to rotate the tower angle"""
            pos = self.angrot.getPosition()
            if gui.anglepos:
                pos = int(gui.anglepos)
            if self.angrot:
                self.angrot.setPosition(pos)

        def Set_ParabolaSlider(gui, *args):
            """Function to move the parabola slider"""
            pos = self.truss.trussGetPosition()
            if gui.pslider:
                pos = float(gui.pslider)
            if self.truss:
                self.truss.moveTrussTo(pos)

        gui_image = Gui(
            [MA("plot"), ___, ___, ___],
            ["Par position:", "parpos", ___, "[-, -, mm, arcsec, arcsec, -]"],
            ["Rm position:", "rmpos", ___, "[-, -, mm, arcsec, arcsec, -]"],
            ["M4 exapode position:", "m4pos", ___, "mm"],
            ["Par slider position:", "pslider", ___, "m"],
            ["Rm slider position:", "rslider", ___, "m"],
            ["Ang rot position:", "anglepos", ___, "deg"],
            [_, _, HB("heart_empty_30.png", "heart_full_30.png"), _],
            images_dir=data_root_dir(),
        )

        control_gui = Gui(
            [
                "New Par position",
                "0",
                "0",
                "__parpist__",
                "__partip__",
                "__partilt__",
                "0",
                "[-, -, mm, arcsec, arcsec, -]",
            ],
            [Set_Parabola, ___, ___, ___, ___, ___, ___, ___],
            [
                "New Rm position",
                "0",
                "0",
                "__rmpist__",
                "__rmtip__",
                "__rmtilt__",
                "0",
                "[-, -, mm, arcsec, arcsec, -]",
            ],
            [Set_RefMirror, ___, ___, ___, ___, ___, ___, ___],
            [
                "New M4 exapode position",
                "0",
                "0",
                "__m4pist__",
                "__m4tip__",
                "__m4tilt__",
                "0",
                "mm",
            ],
            [Set_M4, ___, ___, ___, ___, ___, ___, ___],
            ["New Par Slider position", "__pslider__", _, _, _, _, _, "m"],
            [Set_ParabolaSlider, ___, ___, ___, ___, ___, ___, ___],
            ["New Rm Slider position", "__rslider__", _, _, _, _, _, "m"],
            [Set_ReferenceMirrorSlider, ___, ___, ___, ___, ___, ___, ___],
            ["New Angle Rot position", "__anglepos__", _, _, _, _, _, "deg"],
            [Set_AngleRotator, ___, ___, ___, ___, ___, ___, ___],
        )  # exceptions=Exceptions.OFF)

        ottIma = OttImages(self.ott)
        image = ottIma.ott_m4view()

        gui_image.plot = image
        gui_image.plot.set_title("OTT geometry")
        gui_image.plot.colorbar()
        # gui_image.plot.text(30, 50, 'ciao')

        gui_image.timer_start(getstatus, 1)

        self.gui = Gui([G("OTT"), G("Control")])
        self.gui_control = control_gui
        self.gui_image = gui_image
        self.gui.OTT = gui_image
        self.gui.Control = control_gui

    def runImage(self):
        self._setUp()
        self.gui_image.run()

    def run(self):
        """Run the GUI"""
        self._setUp()
        self.gui.run()


def main():
    from m4.configuration import start
    import os
    conf = os.environ["PYOTTCONF"]
    ott, interf, dm = start.create_ott(conf)

    runner = Runner(ott)
    sys.exit(runner.run())


if __name__ == "__main__":
    main()

"""
Authors
  - C. Selmi: written in 2023
"""

import numpy as np
from m4.ground.package_data import data_root_dir
from m4.analyzers.accelerometers_data_analyzer import AccelerometersDataAnalyzer
from guietta import Gui, G, _, ___, HB, MA
import sys
from m4.configuration.ott_parameters import OpcUaParameters


class Runner:
    """
    Class for creating interactive GUI

    HOW TO USE IT::

        conf = 'myConfGlobalPath.yaml'
        from m4.configuration import start
        ott, interf, dm = start.create_ott(conf)
        from m4.gui import accelerometers_GUI
        g = accelerometers_GUI.Runner(ott)
        g.run()
    """

    def __init__(self, ott):
        """The constructor
        ott: object
            tower object
        """
        self.ott = ott

    def _setUp(self):

        def setPlot(gui, an):
            """
            image: numpy array
                image of the tower
            """
            spe1 = an._spe[:, 1:]
            freq1 = an._freq[1:]

            ax = gui.plot.ax
            ax.clear()
            ax.set_title(an.tt)
            label_list = []
            for i in OpcUaParameters.accelerometers_plc_id:
                ss = (
                    "Ch-"
                    + str(i)
                    + " "
                    + OpcUaParameters.accelerometrs_directions[i - 1]
                )
                label_list.append(ss)
            for i in range(spe1.shape[0]):
                prova = ax.plot(freq1, np.abs(spe1[i, :]), "-")
                ax.legend(prova, label_list[i], loc=0)
                # gui.proxy('label%d'%i).set('%s' % label_list[i])
            ax.set_xlim([0, 100])
            # gui.plot.set_clim(vmin=0, vmax=100)
            ax.set_xlabel("Freq[Hz]")
            ax.set_ylabel("Amplitude Spectrum |m/s2|")
            ax.figure.canvas.draw()
            # gui.heart_empty_30.beat()

        def Start_Acquisition(gui, *args):
            """Function to reload tower data and image"""
            recording_seconds = int(gui.rec_sec) or 5.0
            if self.ott:
                gui.tn = self.ott.accelerometers.acquireData(recording_seconds)

        def Start_Analysis(gui, *args):
            track_num = np.str(gui.track_num)
            an = AccelerometersDataAnalyzer(track_num)
            an._spe, an._freq = an._acc.power_spectrum()
            an.datah5 = an._acc.datah5
            setPlot(gui, an)

        gui_acquisition = Gui(
            ["Set recording time", "__rec_sec__", _, "[s]"],
            [Start_Acquisition, _, _, _],
            [_, _, "Tracking Number generated:", "tn"],
            [_, _, _, HB("heart_empty_30.png", "heart_full_30.png")],
            images_dir=data_root_dir(),
        )

        gui_analysis = Gui(
            ["Set tracking number", "__track_num__"],
            [Start_Analysis, _],
            [MA("plot"), _],
        )

        self.gui = Gui([G("Acc_acquisitions"), G("Acc_analysis")])
        self.gui.Acc_acquisitions = gui_acquisition
        self.gui.Acc_analysis = gui_analysis

    def run(self):
        """Run the GUI"""
        self._setUp()
        self.gui.run()


def main():
    from m4.configuration import start

    conf = "/mnt/m4storage/Data/SYSCONFData/m4Config.yaml"
    ott, interf, dm = start.create_ott(conf)

    runner = Runner(ott)
    sys.exit(runner.run())


if __name__ == "__main__":
    main()

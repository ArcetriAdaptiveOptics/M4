"""
Authors
  - C. Selmi: written in 2019
"""

import logging
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from m4.configuration import config_folder_names as fn


imgFold = fn.OPD_IMAGES_ROOT_FOLDER
ifFold = fn.IFFUNCTIONS_ROOT_FOLDER
intMatFold = fn.INTMAT_ROOT_FOLDER
confFold = fn.CONFIGURATION_ROOT_FOLDER


class ComputeReconstructor:
    """
    This class analyzes the measurements made through the IFF class
    and calculates the reconstructor to be used in the control loop.

    HOW TO USE IT::

        tn = "YYYYMMDD_HHMMSS"
        cr = ComputeReconstructor.loadInteractionMatrix(tn)
        rec = cr.run()

    OR

        cr = ComputeReconstructor(interation_matrix_cube) 
        rec = cr.run(Interactive=True)

        where the interaction_matrix_cube is a masked_array dstack of
        shape [pixels, pixels, n_images]
    """

    def __init__(self, interation_matrix_cube=None, tn=None):
        """The constructor"""
        self._logger = logging.getLogger("COMPUTE_REC:")
        self._intMatCube = interation_matrix_cube
        self._intMat = None
        self._analysisMask = None
        self._intMat_U = None
        self._intMat_S = None
        self._intMat_Vt = None
        self._threshold = None
        self._filtered_sv = None
        self._tn = tn
        if tn is not None:
            self._logger.info(
                "Loaded interaction matrix from TN  %s", self._tn)

        if self._intMatCube is not None:
            analysis_mask = self._intMatCube[:, :, 0].mask
            # compute logical "or" mask of all the images in the cube
            for i in range(1, self._intMatCube.shape[2]):
                analysis_mask = np.logical_or(
                    analysis_mask, self._intMatCube[:, :, i].mask)

            self.setAnalysisMask(analysis_mask)
        else:
            self.analysisMask = None

    def run(self, Interactive=False, sv_threshold=None):

        self._logger.info("Computing reconstructor")
        self._computeIntMat()

        self._logger.info("Computing singular values")
        self._intMat_U, self._intMat_S, self._intMat_Vt = \
            np.linalg.svd(self._intMat, full_matrices=False)
        if Interactive:
            self._threshold = self.make_interactive_plot(self._intMat_S)
        else:
            if sv_threshold is None:
                self._threshold = {
                    "y": np.finfo(np.float32).eps,
                    "x": np.argmin(np.abs(self._intMat_S - np.finfo(np.float32).eps)),
                }
            else:
                self._threshold = {
                    "y": sv_threshold,
                    "x": np.argmin(np.abs(self._intMat_S - sv_threshold)),
                }
        sv_threshold = self._intMat_S.copy()
        sv_threshold[self._threshold["x"]:] = 0
        self._filtered_sv = sv_threshold
        self._logger.info("Assembling reconstructor")
        return self._intMat_Vt.T @ np.diag(sv_threshold) @ self._intMat_U.T

    @staticmethod
    def loadIntMatFromFolder(tn):
        """
        Creates the object using saved information about path measurements
        """
        # load cube from file
        try:
            intmat_fullpath = os.path(
                os.path.join(intMatFold, tn), "IMCube.fits")
            hdu = pyfits.open(intmat_fullpath)
            _intMatCube = hdu[0].data
            hdu.close()
        except Exception as e:
            raise e

        return ComputeReconstructor(interation_matrix_cube=_intMatCube, tn=tn)

    def _computeIntMat(self):
        # compute the interaction matrix from the cube by selecting from each
        #  masked arrey all the pixels given by the analysis mask
        self._logger.info(
            "Computing interaction matrix from cube of size %s", self._intMatCube.shape)
        try:
            self._intMat = np.array(
                [
                    self._intMatCube[:, :, i][self._analysisMask == 0]
                    for i in range(self._intMatCube.shape[2])
                ]
            )
        except Exception as e:
            self._logger.error(
                "Error in computing interaction matrix from cube:%s", e)
            raise e

    def setAnalysisMask(self, analysis_mask):
        """Set the analysis mask chosen

        Parameters
        ----------
                analysis_mask: numpy array [pixels, pixels]
        """
        self._analysisMask = analysis_mask

    def analysisMask(self):
        """Set the analysis mask chosen

        Parameters
        ----------
                analysis_mask: numpy array [pixels, pixels]
        """
        return self._analysisMask

    @staticmethod
    def make_interactive_plot(singular_values, current_threshold=None):

        # Creare il grafico
        fig, ax = plt.subplots()
        modelist = np.arange(len(singular_values))
        ax.loglog(modelist, singular_values, "b-o")
        ax.set_xlabel("Mode number")
        ax.set_ylabel("Singular value")
        ax.grid()
        ax.autoscale(tight=False)
        ax.set_ylim([min(singular_values), max(singular_values)])
        ax.title.set_text("Singular values")
        # if current_threshold is None:
        threshold = dict()

        threshold["y"] = np.finfo(np.float32).eps  # 0.01
        threshold["x"] = 0  # 5 % len(singular_values)
        ax.axhline(threshold["y"], color="g", linestyle="-", linewidth=1)
        ax.axvline(threshold["x"], color="g", linestyle="-", linewidth=1)
        ax.plot(
            modelist[singular_values < threshold["y"]],
            singular_values[singular_values < threshold["y"]],
            "r-x",
        )

        # Funzione per aggiornare la posizione della croce rossa
        def update_crosshair(event):
            if event.inaxes:
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                ax.cla()
                ax.loglog(modelist, singular_values, "b-o")
                ax.plot(
                    modelist[singular_values < threshold["y"]],
                    singular_values[singular_values < threshold["y"]],
                    "r-x",
                )
                ax.autoscale(tight=False)
                ax.set_xlabel("Mode number")
                ax.set_ylabel("Singular value")
                ax.grid()
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                ax.axhline(threshold["y"], color="g",
                           linestyle="-", linewidth=1)
                ax.axvline(threshold["x"], color="g",
                           linestyle="-", linewidth=1)
                ax.title.set_text("Singular values")
                x_mouse, y_mouse = event.xdata, event.ydata
                ax.axhline(y_mouse, color="r", linestyle="-", linewidth=1)
                ax.axvline(
                    np.argmin(np.abs(singular_values - y_mouse)),
                    color="r",
                    linestyle="-",
                    linewidth=1,
                )
                ax.figure.canvas.draw()
                print(
                    f"Current threshold: X = {
                        threshold['x']:.2f}, Y = {threshold['y']:.5f}"
                )

        # Connettere la funzione all'evento di movimento del mouse
        ax.add_artist(ax.patch)
        fig.canvas.mpl_connect(
            "motion_notify_event", lambda event: update_crosshair(event)
        )

        # Funzione per registrare le coordinate del click

        def record_click(event):
            if event.inaxes and event.dblclick:
                x_click, y_click = event.xdata, event.ydata
                threshold["y"] = y_click
                threshold["x"] = np.argmin(np.abs(singular_values - y_click))
                ax.axhline(threshold["x"], color="g",
                           linestyle="-", linewidth=1)
                ax.axvline(threshold["y"], color="g",
                           linestyle="-", linewidth=1)
                print(
                    f"New eigenvalues threshold: X = {
                        x_click:.2f}, Y = {y_click:.2f}"
                )
                return threshold["x"], threshold["y"]

        # Connettere la funzione all'evento di click del tasto sinistro
        fig.canvas.mpl_connect("button_press_event",
                               lambda event: record_click(event))

        # Visualizzare il grafico
        plt.show()
        return threshold


# if __name__ == "__main__":
#    import random
#    random.seed(0)
#    mat = np.array([
#        [random.random() for i in range(100)] for j in range(100)])
#    sqmat = mat @ mat.T
#    rec = ComputeReconstructor.make_interactive_plot_bokeh(sqmat)

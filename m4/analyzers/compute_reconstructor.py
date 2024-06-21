'''
Authors
  - C. Selmi: written in 2019
'''

import os
import logging
import h5py
from astropy.io import fits as pyfits
import numpy as np
from m4.ground import read_data
from m4.ground.read_data import InterferometerConverter
from m4.utils.influence_functions_maker import IFFunctionsMaker
from m4.utils.roi import ROI
from m4.utils.image_reducer import TipTiltDetrend
from m4.configuration import config_folder_names as fold_name
import matplotlib.pyplot as plt


class ComputeReconstructor():
    '''
    This class analyzes the measurements made through the IFF class
    and calculates the reconstructor to be used in the control loop.

    HOW TO USE IT::

        from m4.analyzers.analyzer_iffunctions import AnalyzerIFF
        from m4.analyzers.compute_reconstructor import ComputeReconstructor
        fileName = os.path.join(".../IFFunctions", tt)
        an = AnalyzerIFF.loadInfoFromTtFolder(fileName)
        rec = ComputeReconstructor(an)
        an.setRec() = rec.getReconstructor()
    '''

    def __init__(self, analyzer_iff_obj):
        """The constructor """
        self._logger = logging.getLogger('COMPUTE_REC:')
        self._cube = analyzer_iff_obj.getCube()
        self._analysisMask = None

    @staticmethod
    def loadReconstructorFromFolder(tt):
        """ Creates the object using information about path measurements

            TO BE IMPLEMENTED
        """
        return tt

    def reconstructor(self):
        '''
        Returns
        -------
                cube: masked array [pixels, pixels, number of images]
                    cube from analysis
        '''
        return self._cube

    def analysisMask(self):
        ''' Set the analysis mask chosen

        Parameters
        ----------
                analysis_mask: numpy array [pixels, pixels]
        '''
        return self._analysisMask

    @staticmethod
    def make_interactive_plot(singular_values, current_threshold=None):

        # Creare il grafico
        fig, ax = plt.subplots()
        modelist = np.arange(len(singular_values))
        ax.loglog(modelist, singular_values, 'b-o')
        ax.set_xlabel('Mode number')
        ax.set_ylabel('Singular value')
        ax.grid()
        ax.autoscale(tight=True)
        ax.set_ylim([min(singular_values), max(singular_values)])
        ax.title.set_text('Singular values')
        # if current_threshold is None:
        threshold = dict()

        threshold['y'] = 0.01  # np.finfo(np.float32).eps
        threshold['x'] = 5 % len(singular_values)
        ax.axhline(threshold['y'], color='g',
                   linestyle='-', linewidth=1)
        ax.axvline(threshold['x'],  color='g',
                   linestyle='-', linewidth=1)
        ax.plot(modelist[singular_values < threshold['y']],
                singular_values[singular_values < threshold['y']], 'r-x')

        # Funzione per aggiornare la posizione della croce rossa
        def update_crosshair(event):
            if event.inaxes:
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                ax.cla()
                ax.loglog(modelist, singular_values, 'b-o')
                ax.plot(modelist[singular_values < threshold['y']],
                        singular_values[singular_values < threshold['y']], 'r-x')
                ax.autoscale(tight=True)
                ax.set_xlabel('Mode number')
                ax.set_ylabel('Singular value')
                ax.grid()
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                ax.axhline(threshold['y'], color='g',
                           linestyle='-', linewidth=1)
                ax.axvline(threshold['x'],  color='g',
                           linestyle='-', linewidth=1)
                ax.title.set_text('Singular values')
                x_mouse, y_mouse = event.xdata, event.ydata
                ax.axhline(y_mouse, color='r',
                           linestyle='-', linewidth=1)
                ax.axvline(x_mouse,  color='r',
                           linestyle='-', linewidth=1)
                ax.figure.canvas.draw()
                print(f"Current threshold: X = {
                      threshold['x']:.2f}, Y = {threshold['y']:.5f}")

        # Connettere la funzione all'evento di movimento del mouse
        ax.add_artist(ax.patch)
        fig.canvas.mpl_connect('motion_notify_event',
                               lambda event: update_crosshair(event))

        # Funzione per registrare le coordinate del click

        def record_click(event):
            if event.inaxes and event.dblclick:
                x_click, y_click = event.xdata, event.ydata
                threshold['x'] = x_click
                threshold['y'] = y_click
                ax.axhline(threshold['x'], color='g',
                           linestyle='-', linewidth=1)
                ax.axvline(threshold['y'],  color='g',
                           linestyle='-', linewidth=1)
                print(f"New eigenvalues threshold: X = {
                      x_click:.2f}, Y = {y_click:.2f}")
                return x_click, y_click

        # Connettere la funzione all'evento di click del tasto sinistro
        fig.canvas.mpl_connect('button_press_event',
                               lambda event: record_click(event))

        # Visualizzare il grafico
        plt.show()


if __name__ == '__main__':

    import random
    import unittest
    import mock
    # generate sample data
    random.seed(0)
    random_matrix = [
        [random.random() for i in range(100)] for j in range(1000)]
    random_matrix_dp = np.linalg.inv(
        np.array(random_matrix).T @ np.array(random_matrix))
    us, s, vsT = np.linalg.svd(random_matrix_dp)

    print(s)
    th1, th2 = ComputeReconstructor.make_interactive_plot(s)
    print(th1, th2)

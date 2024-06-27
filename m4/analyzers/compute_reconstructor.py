'''
Authors
  - C. Selmi: written in 2019
'''

# import os
import logging
# import h5py
# from astropy.io import fits as pyfits
import numpy as np
# from m4.ground import read_data
# from m4.ground.read_data import InterferometerConverter
# from m4.utils.influence_functions_maker import IFFunctionsMaker
from m4.utils.roi import ROI
# from m4.utils.image_reducer import TipTiltDetrend
# from m4.configuration import config_folder_names as fold_name
import matplotlib.pyplot as plt
from m4.analyzers.analyzer_iffunctions import AnalyzerIFF


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
-------------------        an.setRec() = rec.getReconstructor()       -------------------
    '''

    def __init__(self, AnalyzerIFF_obj: AnalyzerIFF):
        """The constructor """
        self._logger = logging.getLogger('COMPUTE_REC:')
        self._intMat = AnalyzerIFF_obj.getInteractionMatrix()
        self._analysisMask = None
        self._intMat_U = None
        self._intMat_S = None
        self._intMat_Vt = None
        self._threshold = None
        self._filtered_sv = None

    def run(self, Interactive=False, sv_threshold=None):
        self._intMat_U, self._intMat_S, self._intMat_Vt = np.linalg.svd(
            self._intMat)
        if Interactive:
            self._threshold = self.make_interactive_plot(self._intMat_S)
        else:
            if sv_threshold is None:
                self._threshold = {'y': np.finfo(np.float32).eps, 'x':
                                   np.argmin(np.abs(self._intMat_S -
                                                    np.finfo(np.float32).eps))}
            else:
                self._threshold = {'y': sv_threshold, 'x':
                                   np.argmin(np.abs(self._intMat_S -
                                                    sv_threshold))}
        sv_threshold = self._intMat_S.copy()
        sv_threshold[self._threshold['x']:] = 0
        self._filtered_sv = sv_threshold
        return self._intMat_Vt.T @ np.diag(sv_threshold) @ self._intMat_U.T

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
        ax.autoscale(tight=False)
        ax.set_ylim([min(singular_values), max(singular_values)])
        ax.title.set_text('Singular values')
        # if current_threshold is None:
        threshold = dict()

        threshold['y'] = np.finfo(np.float32).eps  # 0.01
        threshold['x'] = 0  # 5 % len(singular_values)
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
                ax.autoscale(tight=False)
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
                ax.axvline(np.argmin(np.abs(singular_values - y_mouse)),  color='r',
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
                threshold['y'] = y_click
                threshold['x'] = np.argmin(np.abs(singular_values - y_click))
                ax.axhline(threshold['x'], color='g',
                           linestyle='-', linewidth=1)
                ax.axvline(threshold['y'],  color='g',
                           linestyle='-', linewidth=1)
                print(f"New eigenvalues threshold: X = {
                      x_click:.2f}, Y = {y_click:.2f}")
                return threshold['x'], threshold['y']

        # Connettere la funzione all'evento di click del tasto sinistro
        fig.canvas.mpl_connect('button_press_event',
                               lambda event: record_click(event))

        # Visualizzare il grafico
        plt.show()

    def make_interactive_plot_bokeh(singular_values, current_threshold=None):

        # still to be debugged
        from bokeh.plotting import figure, show
        from bokeh.models import ColumnDataSource
        from bokeh.models.tools import HoverTool
        from bokeh.layouts import column
        from bokeh.events import DoubleTap
        from bokeh.io import curdoc
        from bokeh.models import Span
        from bokeh.models import CustomJS

        # Create the plot
        fig = figure(title='Singular values', x_axis_label='Mode number',
                     y_axis_label='Singular value', tools='pan,box_zoom,reset,save')
        modelist = np.arange(len(singular_values))
        fig.line(modelist, singular_values, line_width=2)
        fig.circle(modelist, singular_values, fill_color="white", size=8)
        fig.grid.grid_line_alpha = 0.3
        fig.title.text_font_size = '16pt'
        fig.xaxis.axis_label_text_font_size = "14pt"
        fig.yaxis.axis_label_text_font_size = "14pt"
        fig.xaxis.major_label_text_font_size = "12pt"
        fig.yaxis.major_label_text_font_size = "12pt"

        # Add a threshold line
        threshold = dict()
        threshold['y'] = np.finfo(np.float32).eps  # 0.01
        threshold['x'] = 0  # 5 % len(singular_values)
        fig.line([0, len(singular_values)], [threshold['y'], threshold['y']],
                 line_color='green', line_width=2, line_dash='dashed')
        fig.line([threshold['x'], threshold['x']], [0, max(singular_values)],
                 line_color='green', line_width=2, line_dash='dashed')
        fig.line(modelist[singular_values < threshold['y']],
                 singular_values[singular_values < threshold['y']], line_color='red', line_width=2)

        # Function to update the position of the red cross
        def update_crosshair(event):
            if event.x and event.y:
                x_mouse, y_mouse = event.x, event.y
                threshold['y'] = y_mouse
                threshold['x'] = np.argmin(np.abs(singular_values - y_mouse))
                fig.line([0, len(singular_values)], [threshold['y'], threshold['y']],
                         line_color='green', line_width=2, line_dash='dashed')
                fig.line([threshold['x'], threshold['x']], [0, max(singular_values)],
                         line_color='green', line_width=2, line_dash='dashed')
                fig.line(modelist[singular_values < threshold['y']],
                         singular_values[singular_values < threshold['y']], line_color='red', line_width=2)
                print(f"Current threshold: X = {
                      threshold['x']:.2f}, Y = {threshold['y']:.5f}")

        # Connect the function to the mouse movement event
        fig.on_event(DoubleTap, update_crosshair)

        # Function to record the coordinates of the click
        def record_click(event):
            if event.x and event.y:
                x_click, y_click = event.x, event.y
                threshold['y'] = y_click
                threshold['x'] = np.argmin(np.abs(singular_values - y_click))
                fig.line([0, len(singular_values)], [threshold['y'], threshold['y']],
                         line_color='green', line_width=2, line_dash='dashed')
                fig.line([threshold['x'], threshold['x']], [0, max(singular_values)],
                         line_color='green', line_width=2, line_dash='dashed')
                print(f"New eigenvalues threshold: X = {
                      x_click:.2f}, Y = {y_click:.2f}")
                return threshold['x'], threshold['y']

        # Connect the function to the left button click event
        fig.on_event(DoubleTap, record_click)

        # Show the plot
        show(fig)

        return threshold['x'], threshold['y']


# if __name__ == "__main__":
#    import random
#    random.seed(0)
#    mat = np.array([
#        [random.random() for i in range(100)] for j in range(100)])
#    sqmat = mat @ mat.T
#    rec = ComputeReconstructor.make_interactive_plot_bokeh(sqmat)

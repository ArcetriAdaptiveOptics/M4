'''
@author: cs
'''

import numpy as np
from m4.ground import fit_ellipse
from m4.utils.roi import ROI
from astropy.io import fits as pyfits

class ParabolIdent():
    '''
    '''

    def __init__(self):
        """The constructor """
        self._r = ROI()

    def _fiduciali(self, image):
        roiList = self._r.roiGenerator(image)

    #xx = np.matrix('1 2 5 7 9 6 3 8; 7 6 8 7 5 7 2 4')
    def ellipse(self, punti_fiduciali):
        z,a,b,alpha = fit_ellipse.fitellipse(punti_fiduciali, 'linear')

        Q_matrix = np.zeros((2,2))
        Q_matrix[0,0] = np.cos(alpha)
        Q_matrix[0,1] = -np.sin(alpha)
        Q_matrix[1,0] = np.sin(alpha)
        Q_matrix[1,1] = np.cos(alpha)
        C_matrix = np.zeros((2,1))
        C_matrix[0,0] = a * np.cos(alpha)
        C_matrix[1,0] = b * np.sin(alpha)
        ellipse = z + Q_matrix * C_matrix

    #ss = np.array((x1,x2,x3,x4,x5,x6))
    def ellipse2(self, punti_fiduciali):
        import cv2
        ellipse = cv2.fitEllipse(punti_fiduciali)

    def imaTest(self):
        file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/20161226_122557/mode_0569.fits'
        hduList = pyfits.open(file_name)
        ima = hduList[0].data
        immagine = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
        return immagine

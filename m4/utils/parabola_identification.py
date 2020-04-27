'''
@author: cs
'''

import numpy as np
#import cv2
from numpy.linalg import eig, inv
from m4.ground import fit_ellipse
from skimage.draw import circle as draw_circle
from m4.utils.roi import ROI
from skimage.draw import ellipse as draw_ellipse
from astropy.io import fits as pyfits

class ParabolIdent():
    '''
    '''

    def __init__(self):
        """The constructor """
        self._r = ROI()

    def parable(self, image):
        x, y = self._fiduciali(image)
        centro, axs, raggio = self._fitEllipse(x, y)
        ellipse = self._drawEllipse(centro, axs, image)
        circle = self._drawCircle(centro, raggio, image)
        return ellipse, circle, centro, axs, raggio

    def _fiduciali(self, image):
        roiList = self._r.roiGenerator(image)
        x = np.array([512, 359, 512, 665])
        y = np.array([359, 512, 665, 512])
        return x, y

    def prova(self, image, centro, raggio):
        size = np.array([image.shape[0], image.shape[1]])
        ima_x = np.arange(size[0], dtype = float)
        ima_y = np.arange(size[1], dtype = float)
        xx = np.tile(ima_x, (size[0], 1))-centro[0] / raggio
        yy = np.tile(ima_y, (size[1], 1)).T-centro[1] / raggio
        return xx, yy

#     def ellipse(self, punti_fiduciali):
#         #xx = np.matrix('1 2 5 7 9 6 3 8; 7 6 8 7 5 7 2 4')
#         z,a,b,alpha = fit_ellipse.fitellipse(punti_fiduciali, 'linear')
# 
#         Q_matrix = np.zeros((2,2))
#         Q_matrix[0,0] = np.cos(alpha)
#         Q_matrix[0,1] = -np.sin(alpha)
#         Q_matrix[1,0] = np.sin(alpha)
#         Q_matrix[1,1] = np.cos(alpha)
#         C_matrix = np.zeros((2,1))
#         C_matrix[0,0] = a * np.cos(alpha)
#         C_matrix[1,0] = b * np.sin(alpha)
#         ellipse = z + Q_matrix * C_matrix
#         return ellipse

#     def ellipse2(self, punti_fiduciali):
#         #ss = np.array((x1,x2,x3,x4,x5,x6))
#         ellipse = cv2.fitEllipse(punti_fiduciali)
#         img = np.zeros((15,15))
#         aa = cv2.ellipse(img, ellipse, (255,0, 255), 1, cv2.LINE_AA)
#         return aa

    def _fitEllipse(self, x, y):
        #fit ellipse
        x = x[:,np.newaxis]
        y = y[:,np.newaxis]
        D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
        S = np.dot(D.T,D)
        C = np.zeros([6,6])
        C[0,2] = C[2,0] = 2; C[1,1] = -1
        E, V =  eig(np.dot(inv(S), C))
        n = np.argmax(np.abs(E))
        a_vect = V[:,n]
        #center
        b,c,d,f,g,a = a_vect[1]/2, a_vect[2], a_vect[3]/2, a_vect[4]/2, a_vect[5], a_vect[0]
        num = b*b-a*c
        x0=(c*d-b*f)/num
        y0=(a*f-b*d)/num
        centro = np.array([x0,y0])
        #lunghezza degli assi        
        up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
        down1=(b*b-a*c)*((c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        down2=(b*b-a*c)*((a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        res1=np.sqrt(up/down1)
        res2=np.sqrt(np.abs(up/down2)) #attenzione
        axs = np.array([res1, res2])
        raggio = axs.mean()
        return centro, axs, raggio

    def _drawEllipse(self, centro, axs, image):
        ellipse = np.zeros((image.shape[0],image.shape[1]))
        rr, cc = draw_ellipse(centro[1], centro[0], axs[1], axs[0])
        ellipse[rr, cc]=1
        return ellipse

    def _drawCircle(self, centro, raggio, image):
        circle = np.zeros((image.shape[0],image.shape[1]))
        rr, cc = draw_circle(centro[1], centro[0], raggio)
        circle[rr, cc]=1
        return circle

    def imaTest(self):
        file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/20161226_122557/mode_0569.fits'
        hduList = pyfits.open(file_name)
        ima = hduList[0].data
        immagine = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
        return immagine

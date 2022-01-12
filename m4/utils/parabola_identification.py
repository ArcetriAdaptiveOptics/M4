'''
Authors
  - C. Selmi: written in 2020
'''

import numpy as np
#import cv2
from numpy.linalg import eig, inv
from skimage.draw import circle as draw_circle
from skimage.draw import ellipse as draw_ellipse
from skimage.measure import label
from m4.ground import read_data
import sklearn.feature_extraction as skf_e
import sklearn.cluster as skc
from m4.configuration.ott_parameters import OttParameters

class ParabolIdent():
    ''' Class to be used to determine the position of the parable

    HOW TO USE IT::

        from m4.utils.parabola_identification import ParabolIdent
        pi = ParabolIdent()
    '''

    def __init__(self):
        """The constructor """
        self._rFiducialPoint = OttParameters.RADIUS_FIDUCIAL_POINT

    def testZernikeOnPar(self, image):
        '''
        Parameters
        ----------
        image: numpy masked array
            image to test

        Returns
        -------
        coef1:
        coef2:
        '''
        from m4.ground import zernike
        x, y = self.fiduciali(image)
        centro, axs, raggio = self._fitEllipse(x, y)
        pxs = raggio / OttParameters.RADIUS_FIDUCIAL_POINT
        par_radius = pxs * OttParameters.parab_radius
#         import pdb
#         pdb.set_trace()
        circle = self._drawCircle(centro, par_radius, image)
        ima = np.ma.masked_array(image, mask=np.invert(circle.astype(bool)))
        coef1, mat = zernike.zernikeFit(ima, np.arange(10)+1)

        xx, yy = self.coord(image, centro, raggio)
        mm = (circle==1)
        coef2 = zernike._surf_fit(xx[mm], yy[mm], image.data[mm], np.arange(10)+1)
        return coef1, coef2

    def parable(self, image):
        '''
        Parameters
        ----------
            image: numpy masked array
                image to use

        Returns
        -------
            circle: numpy array
                    circumference of 1 to be plotted
            centro: numpy array
                    coordinates of the center
            axs: numpy array
                    major and minor axis coming from the fit of the ellipse
            raggio: int
                radius of the parabola circumference
        '''
        x, y = self.fiduciali(image)
        centro, axs, raggio = self._fitEllipse(x, y)
        #ellipse = self._drawEllipse(centro, axs, image)
        circle = self._drawCircle(centro, raggio, image)
        return circle, centro, axs, raggio

    def _fiducialiAMano(self):
        ''' Guardati da imshow di image
        '''
        x = np.array([512, 359, 512, 665])
        y = np.array([359, 512, 665, 512])
        return x, y

    def fiduciali(self, ima):
        ''' Calculates the coordinates of the fiducial points of the parabola and
        return it in a single vector of x and y

        Parameters
        ----------
            image: numpy masked array
                image from which to extract the fiducial point

        Returns
        -------
            x: numpy array
                vector of x coordinates of fiducial points
            y: numpy array
                vector of y coordinates of fiducial points
        '''
        graph = skf_e.image.img_to_graph(ima.data, mask=ima.mask)
        labels = skc.spectral_clustering(graph, n_clusters=4, eigen_solver='arpack')
        labels = label(ima.mask.astype(int))
        roiList = []
        for i in range(1, 13):
            maski = np.zeros(labels.shape, dtype=np.bool)
            maski[np.where(labels == i)] = 1
            roiList.append(maski)
        x_list = []
        y_list = []
        for i in np.array([4,6,7,8]):
            aa = np.where(roiList[i]==1)
            x = aa[1].mean()
            y = aa[0].mean()
            x_list.append(x)
            y_list.append(y)
        return np.array(x_list), np.array(y_list)

    def coord(self, image, centro, raggio):
        '''
        Parameters
        ----------
            image: numpy masked array
            centro: numpy array
                     coordinates of the center
            raggio: int
                 radius of the parabola circumference

        Returns
        -------
            xx: numpy array
                coordinate x della parabola
            yy: numpy array
                coordinate y della parabola
        '''
        raggio = raggio / self._rFiducialPoint
        size = np.array([image.shape[0], image.shape[1]])
        ima_x = np.arange(size[0], dtype = float)
        ima_y = np.arange(size[1], dtype = float)
        xx = (np.tile(ima_x, (size[0], 1))-centro[0]) / raggio
        yy = (np.tile(ima_y, (size[1], 1)).T-centro[1]) / raggio
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
        '''
        args:
            x = vector of the x coordinates of the points to be used in the fit
            y = vector of the y coordinates of the points to be used in the fit

        returns:
            centro = coordinates of the center
            axs = major and minor axis coming from the fit of the ellipse
            raggio = radius of the parabola circumference
        '''
        #fit ellipse
        x = x[:, np.newaxis]
        y = y[:, np.newaxis]
        D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
        S = np.dot(D.T, D)
        C = np.zeros([6, 6])
        C[0, 2] = C[2, 0] = 2
        C[1, 1] = -1
        E, V =  eig(np.dot(inv(S), C))
#         import pdb
#         pdb.set_trace()
        n = np.argmax(np.abs(E))
        a_vect = V[:, n]
        #center
        b, c, d, f, g, a = a_vect[1]/2, a_vect[2], a_vect[3]/2, a_vect[4]/2, a_vect[5], a_vect[0]
        num = b*b - a*c
        x0 = (c*d - b*f)/num
        y0 = (a*f - b*d)/num
        centro = np.array([x0, y0])
        #lunghezza degli assi        
        up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
        down1 = (b*b-a*c)*((c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        down2 = (b*b-a*c)*((a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        res1 = np.sqrt(np.abs(up/down1)) #prima non c'era il valore assoluto
        res2 = np.sqrt(np.abs(up/down2)) #attenzione
        axs = np.array([res1, res2])
        raggio = axs.mean()
        return centro, axs, raggio

    def _drawEllipse(self, centro, axs, image):
        '''
        args:
            centro = coordinates of the center
            axs = major and minor axis coming from the fit of the ellipse
            image = masked array

        returns:
            ellipse = ellipse of one to display
        '''
        ellipse = np.zeros((image.shape[0],image.shape[1]))
        rr, cc = draw_ellipse(centro[1], centro[0], axs[1], axs[0])
        ellipse[rr, cc]=1
        return ellipse

    def _drawCircle(self, centro, raggio, image):
        '''
        args:
            centro = coordinates of the center
            raggio = radius of the parabola circumference
            image = masked array

        returns:
            circle = circle of one to display
        '''
        circle = np.zeros((image.shape[0], image.shape[1]))
        rr, cc = draw_circle(centro[1], centro[0], raggio/self._rFiducialPoint)
        circle[rr, cc] = 1
        return circle

    def _imaTest(self, file_path):
        from astropy.io import fits as pyfits
        #file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/20161226_122557/mode_0569.fits'
        #file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/Seg/img_0000.fits'
        #file_name = '/Users/rm/eclipse-workspace/M4/test/utils/img_0000.fits'
        hduList = pyfits.open(file_path)
        immagine = np.ma.masked_array(hduList[0].data[0,:,:],
                                      mask=np.invert(hduList[0].data[1,:,:].astype(bool)))  
        #immagine = read_data.read_phasemap(file_path)
        return immagine

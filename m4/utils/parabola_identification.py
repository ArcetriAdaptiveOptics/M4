'''
Authors
  - C. Selmi: written in 2020
'''

import numpy as np
import pandas as pd
from skimage.measure import label, regionprops_table
#from m4.ground import zernike
from m4.ground import geo
from numpy.linalg import eig, inv
from skimage.draw import disk as draw_circle
from m4.configuration.ott_parameters import OttParameters

class ParabolaCirculaPupil():
    ''' Class to be used to determine the position of the parable

    HOW TO USE IT::

        from m4.utils.parabola_identification import ParabolaCirculaPupil
        pz = ParabolaCirculaPupil()
        circle_mask = pz.par_mask_on_ott(image)
    '''

    def __init__(self):
        """The constructor """
        self._rFiducialPoint = OttParameters.RADIUS_FIDUCIAL_POINT
        self._nmarkers = None

    def par_mask_on_ott(self, image):
        image_masked_central_fid = geo.draw_mask(image, np.int(image.shape[0]/2), np.int(image.shape[1]/2),
                                                OttParameters.INNER_MARKERS_REJECTION_RADIUS)
        imaf = self.fiduciali(image_masked_central_fid)
        centro, axs, raggio = self._fitEllipse(imaf[0], imaf[1])
        pxs = raggio / OttParameters.RADIUS_FIDUCIAL_POINT
        par_radius = pxs * OttParameters.parab_radius
        circle_mask = self._drawCircle(centro, par_radius, image)
        return circle_mask


    def fiduciali(self, image):
        image = image.mask
        image_not_blobs = label(image == 0)
        image_blobs = label(image != 0)
        properties =['area','centroid',
                     'major_axis_length', 'minor_axis_length',
                     'eccentricity','bbox']
        image_df = pd.DataFrame(regionprops_table(image_blobs, properties = properties))
        sel_image = (image_df['area'] < 1000) &  (image_df['area'] > np.median(image_df['area'])/2)
        imaf = np.array([image_df['centroid-0'][sel_image],image_df['centroid-1'][sel_image]])
        nmarkers = np.max(imaf.shape)
        self._nmarkers = nmarkers
        print(nmarkers)
        return imaf

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
        rr, cc = draw_circle((centro[1], centro[0]), raggio)
        circle[rr, cc] = 1
        return circle

    def _imaTest(self, file_path):
        from astropy.io import fits as pyfits
        #file_path = '/Users/rm/eclipse-workspace/M4/test/utils/img_0000.fits'
        hduList = pyfits.open(file_path)
        immagine = np.ma.masked_array(hduList[0].data[0,:,:],
                                      mask=np.invert(hduList[0].data[1,:,:].astype(bool)))  
        return immagine

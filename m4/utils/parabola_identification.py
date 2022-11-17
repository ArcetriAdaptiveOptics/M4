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

from m4.ground import tracking_number_folder
from m4.configuration import config_folder_names as fold_name
import time

class ParabolaActivities():
    ''' Class to be used to determine the position of the parable

    HOW TO USE IT::

        from m4.utils.parabola_identification import ParabolaActivities
        pa = ParabolaActivities()

        circle_mask = pa.par_mask_on_ott(image)
        or
        tt, cube = pa.parab_cgh_measure(interf, n_frames, delay)
        or
        pa.check_concentricity(image)
    '''

    def __init__(self):
        """The constructor """
        self._rFiducialPoint = OttParameters.RADIUS_FIDUCIAL_POINT
        self._nmarkers = None

    def par_mask_on_ott(self, image):
        ''' NOTE: this function uses the parameter INNER_MARKERS_REJECTION_RADIUS
            to hide the central markers

        Parameters
        ----------
        image: numpy masked array

        Returns
        -------
        circle_mask: numpy ndarray
            circolar mask representing the parabola (zeros for masked point)
        '''
        image_masked_central_fid = geo.draw_mask(image, np.int(image.shape[0]/2), np.int(image.shape[1]/2),
                                       OttParameters.INNER_MARKERS_REJECTION_RADIUS)
        imaf = self.fiduciali(image_masked_central_fid)
        centro, axs, raggio = self._fitEllipse(imaf[0], imaf[1])
        pxs = raggio / OttParameters.RADIUS_FIDUCIAL_POINT
        par_radius = pxs * OttParameters.parab_radius
        circle_mask = self._drawCircle(centro, par_radius, image)
        return circle_mask


    def fiduciali(self, image):
        ''' Function to obtain the position of parabola markers
        Parameters
        ----------
        image: numpy masked array
            image on which to search for markers

        Returns
        -------
        imaf: numpy array [2, n_markers]
            x and y coordinates for markers
        '''
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

    def parab_cgh_measure(self, interf, n_frames, delay=0):
        ''' Function for data acquisition and saving
        Parameters
        ----------
        interf: object
            interferometer object create whit the start up
        n_frames: int
            numbers of frame to acquire
        delay: int [s]
            delay between measurements

        Returns
        -------
        tt: string
            tracking number folder
        cube: numpy masked array
            cube of measurements save
        '''
        dove, tt = tracking_number_folder.createFolderToStoreMeasurements(fold_name.PARABOLA_CGH_FOLDER)

        cube_list = []
        for i in range(n_frames):
            masked_image = interf.acquire_phasemap()
            file_name = 'image_%04d' %i
            interf.save_phasemap(dove, file_name, masked_image)
            cube_list.append()
            time.sleep(delay)
        cube = np.dstack(cube_list)
        return tt, cube

    def check_concentricity(self, image):
        '''
        Parameters
        ----------
        image: numpy masked array

        Returns
        -------
        print the centre, axis and radius (mean of axis)
        '''
        min_inner_markers_rejection_radius = np.array(['n_pixel1',
                                                       'n_pixel2',
                                                       'n_pixel3'])
        max_inner_markers_rejection_radius = np.array(['n_pixel1',
                                                       'n_pixel2',
                                                       'n_pixel3'])
        for i in range(min_inner_markers_rejection_radius.shape[0]):
            mask1 = self._draw_mask(image, np.int(image.shape[0]/2),
                                                np.int(image.shape[1]/2),
                                                min_inner_markers_rejection_radius[i])
            mask2 = self._draw_mask(image, np.int(image.shape[0]/2),
                                                np.int(image.shape[1]/2),
                                                max_inner_markers_rejection_radius[i])
            new_image_mask = np.ma.mask_or(mask1, np.invert(mask2))
            image_masked_central_fid = np.ma.masked_array(image, new_image_mask)

            imaf = self.fiduciali(image_masked_central_fid)
            centro, axs, raggio = self._fitEllipse(imaf[0], imaf[1])
            print(centro, axs, raggio)

    ## Funtions from Runa modified##
    def _draw_mask(self, img, cx, cy, r, out=0):
        """ Function to create circular mask
        Created by Runa

        Parameters
        ----------
        img: numpy array
            image to mask
        cx: int [pixel]
            center x of the mask
        cy: int [pixel]
            center y of the mask
        r: int [pixel]
            radius of the mask

        Returns
        -------
        img1: numpy array
            start image mask whit circular new mask
        """
        ss = np.shape(img)
        x = np.arange(ss[0])
        x = np.transpose(np.tile(x, [ss[1], 1]))
        y = np.arange(ss[1])
        y = np.tile(y, [ss[0], 1])
        x = x - cx
        y = y - cy
        nr = np.size(r)
        if nr == 2:
            rr = x*x/r[0]**2+y*y/r[1]**2
            r1 = 1
        else:
            rr = x*x+y*y
            r1 = r**2
        pp = np.where(rr < r1)
        img1 = img.mask.copy()
        if out == 1:
            img1[pp] = 0
        else:
            img1[pp] = 1
        #plt.imshow(img1)
        return img1

'''
@author: cs
'''
import os
import numpy as np
from astropy.io import fits as pyfits
from photutils.centroids import fit_2dgaussian
from m4.ground.configuration import Configuration

class GeomTransf():

    def __init__(self):
        """The constructor """
        self._pixelCut = 15

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return os.path.join(Configuration.OPD_DATA_FOLDER,
                            "GeomTransf")


    def _readMeasure(self):
        tt1 = '20200511_134609'
        tt2 = '20200511_134613'
        file_path1 = os.path.join(GeomTransf._storageFolder(), tt1, 'act-marker.fits')
        file_path2 = os.path.join(GeomTransf._storageFolder(), tt2, 'act-marker.fits')
        hduList = pyfits.open(file_path1)
        ima = hduList[0].data
        immagine1 = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
        hduList = pyfits.open(file_path2)
        ima = hduList[0].data
        immagine2 = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
        return immagine1, immagine2

    def fiduciali(self, ima):
        par_ima_list = []
        for i in range(5):
            ima_cut, y_min, x_min = self._imaCutter(ima, self._pixelCut)
            par_cut = fit_2dgaussian(ima_cut)._parameters
            par_ima = np.copy(par_cut)
            par_ima[2] = np.int(par_cut[2] + x_min)
            par_ima[3] = np.int(par_cut[3] + y_min)
            coord = np.array([par_ima[2], par_ima[3]])
            ima[par_ima[3].astype(int)-self._pixelCut: par_ima[3].astype(int)+self._pixelCut,
                par_ima[2].astype(int)-self._pixelCut: par_ima[2].astype(int)+self._pixelCut] = 0
            par_ima_list.append(coord)
        return np.array(par_ima_list)
    #ha scambiato gli ultimi due punti

    def _imaCutter(self, image, px):
        y_peak, x_peak = np.where(image == np.max(image))
        y_min = y_peak[0]-px
        x_min = x_peak[0]-px
        image_cut = image[y_min:y_peak[0]+px+1, x_min:x_peak[0]+px+1]
        return image_cut, y_min, x_min

    def transf_one_point(self, ima, xy_coord):
        x_coord = xy_coord[0]
        y_coord = xy_coord[1]
        size = np.array([ima.shape[0], ima.shape[1]])
        ima_x = np.arange(size[0], dtype = float)
        ima_y = np.arange(size[1], dtype = float)
        xx = np.tile(ima_x, (size[0], 1))-x_coord
        yy = np.tile(ima_y, (size[1], 1)).T-y_coord
        return xx, yy

    def main_one_point(self):
        im0, im1 = self._readMeasure()
        par_list0 = self.fiduciali(im0)
        par_list1 = self.fiduciali(im1)
        xy_coord0 = par_list0[0]
        xy_coord1 = par_list1[0]
        xx0, yy0 = self.transf_one_point(im0, xy_coord0)
        xx1, yy1 = self.transf_one_point(im1, xy_coord1)

        xx0_inv = np.linalg.pinv(xx0)
        k = np.dot(xx1, xx0_inv)
        return k, xx0, xx1

    def prova(self):
        x0 = np.zeros((5,2))
        x0[0,0] = 115 ; x0[0,1] = 322
        x0[1,0] = 284 ; x0[1,1] = 361
        x0[2,0] = 387 ; x0[2,1] = 322
        x0[3,0] = 217 ; x0[3,1] = 224
        x0[4,0] = 262 ; x0[4,1] = 145

        x1 = np.zeros((5,2))
        x1[0,0] = 102 ; x1[0,1] = 323
        x1[1,0] = 272 ; x1[1,1] = 368
        x1[2,0] = 375 ; x1[2,1] = 333
        x1[3,0] = 208 ; x1[3,1] = 229
        x1[4,0] = 256 ; x1[4,1] = 152

        x0_inv = np.linalg.pinv(x0)
        k = np.dot(x1,x0_inv)
        return x0, x1, k
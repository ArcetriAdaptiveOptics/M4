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

    def principal_main(self, point_to_use):
        im0, im1 = self._readMeasure()
        if point_to_use == 1:
            k, xx0, xx1 = self.main_one_point(im0, im1)
            return k, xx0, xx1
        elif point_to_use == 5:
            nc, x0, x1 = self.main_five_points(im0, im1)
            return nc, x0, x1

    def fiduciali(self, ima, point_to_use):
        par_ima_list = []
        for i in range(point_to_use):
            ima_cut, y_min, x_min = self._imaCutter(ima, self._pixelCut)
            par_cut = fit_2dgaussian(ima_cut)._parameters
            par_ima = np.copy(par_cut)
            par_ima[2] = np.int(par_cut[2] + x_min)
            par_ima[3] = np.int(par_cut[3] + y_min)
            coord = np.array([par_cut[2] + x_min, par_cut[3] + y_min])
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

    def main_five_points(self, im0, im1):
        point_to_use = 5
        x1_no = self.fiduciali(im1, point_to_use)
        x1 = self._invertRighe(x1_no).T
        par_x0 = self.fiduciali(im0, point_to_use)
        x0 = par_x0.T
        mat = self._grid_five_points(par_x0)
        cmat = self._calc_cmat_five_points(x1)
        rec = np.linalg.pinv(mat)
        tran = np.dot(x0, rec)
        nc = np.dot(tran, cmat)
        return nc, x0, x1

    def _invertRighe(self, par):
        new_par = np.copy(par)
        new_par[3,:] = par[4,:]
        new_par[4,:] = par[3,:]
        return new_par

    def _grid_five_points(self, par_list):
        x = par_list.T[0]
        y = par_list.T[1]
        mat = self._calc_mat_five_points(x, y)
        return mat

    def _calc_mat_five_points(self, x, y):
        mat = np.zeros((6, x.shape[0]))
        for i in range(x.shape[0]):
            mat[0,i]  = x[i]
            mat[1,i]  = y[i]
            mat[2,i]  = x[i]**2
            mat[3,i]  = y[i]**2
            mat[4,i]  = x[i]*y[i]
            mat[5,i]  =   1.
        return mat

    def _calc_cmat_five_points(self, x1):
        cmat = np.zeros((6, x1.shape[1]))
        thex = x1[0] #non facendo il reshape non cambia nulla da x
        they = x1[1]
        cmat[0,:] = thex
        cmat[1,:] = they
        cmat[2,:] = thex**2
        cmat[3,:] = they**2
        cmat[4,:] = thex*they
        cmat[5,:] = np.ones(thex.shape[0])
        return cmat



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

###
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
###

    def _calc_grid_one_point(self, ima, xy_coord):
        x_coord = xy_coord[0]
        y_coord = xy_coord[1]
        size = np.array([ima.shape[0], ima.shape[1]])
        ima_x = np.arange(size[0], dtype = float)
        ima_y = np.arange(size[1], dtype = float)
        xx = np.tile(ima_x, (size[0], 1))-x_coord
        yy = np.tile(ima_y, (size[1], 1)).T-y_coord
        return xx, yy

    def _grid_one_point(self, ima, point_to_use):
        par_list = self.fiduciali(ima, point_to_use)
        xy_coord = par_list[0]
        xx, yy = self._calc_grid_one_point(ima, xy_coord)
        return xx, yy

    def main_one_point(self, im0, im1):
        point_to_use = 1
        xx0, yy0 = self._grid_one_point(im0, point_to_use)
        xx1, yy1 = self._grid_one_point(im1,point_to_use)

        xx0_inv = np.linalg.pinv(xx0)
        k = np.dot(xx1, xx0_inv)
        return k, xx0, xx1
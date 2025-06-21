"""
Authors
  - C. Selmi: written in 2022
"""

import os
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration import folders as fold_name
from photutils.centroids import fit_2dgaussian


class ShellActuatorFootprintRegistration:
    """ """

    def __init__(self):
        """The constructor"""
        self._pixelCut = 15

    @staticmethod
    def _storageFolder():
        """Creates the path where to save measurement data"""
        return fold_name.GEOTRANSF_ROOT_FOLDER

    def peak_coord_finder(self, image, n_point_to_find):
        """
        Parameters
        ----------
            image: numpy array
                IFFs image
            n_point_to_find: int
                number of actuator moved to obtain IFFs image

        Returns
        -------
            xy_coord: numpy array [n_point_to_find, 2]
                list of x, y coord
        """
        ima = np.copy(image)
        shell_act_coord_list = []
        for i in range(n_point_to_find):
            ima_cut, y_min, x_min = self._imaCutter(ima, self._pixelCut)
            par_cut = fit_2dgaussian(ima_cut)._parameters
            par_ima = np.copy(par_cut)
            par_ima[2] = int(par_cut[2] + x_min)
            par_ima[3] = int(par_cut[3] + y_min)
            coord = np.array([par_cut[2] + x_min, par_cut[3] + y_min])
            ima[
                par_ima[3].astype(int)
                - self._pixelCut : par_ima[3].astype(int)
                + self._pixelCut,
                par_ima[2].astype(int)
                - self._pixelCut : par_ima[2].astype(int)
                + self._pixelCut,
            ] = 0
            shell_act_coord_list.append(coord)
        return np.array(shell_act_coord_list)

    def _imaCutter(self, image, px):
        """
        args:
            image = np.masked_array of image
            px = number of pixel with which to cut the image

        returns:
            image_cut = cut image
            y_min = y min used for the cut
            x_min = x min udes for the cut
        """
        y_peak, x_peak = np.where(image == np.max(image))
        y_min = y_peak[0] - px
        x_min = x_peak[0] - px
        image_cut = image[y_min : y_peak[0] + px + 1, x_min : x_peak[0] + px + 1]
        return image_cut, y_min, x_min

    def move_peak_whit_shift_function(self, image, img_coord, ref_coord):
        """ """
        coord2 = ref_coord.copy()
        coord2[3] = ref_coord[4]
        coord2[4] = ref_coord[3]
        x = "?"
        y = "?"

        from scipy.ndimage import shift

        data = shift(image, [x, y])
        new_image = np.ma.masked_array(data, mask=image.mask)
        return new_image

    ### test ###
    def _readMeasure(self):
        """Legge due immagini di prova"""
        tt1 = "20200511_134609"
        tt2 = "20200511_134613"
        file_path1 = os.path.join(
            ShellActuatorFootprintRegistration._storageFolder(), tt1, "act-marker.fits"
        )
        file_path2 = os.path.join(
            ShellActuatorFootprintRegistration._storageFolder(), tt2, "act-marker.fits"
        )
        hduList = pyfits.open(file_path1)
        ima = hduList[0].data
        immagine1 = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
        hduList = pyfits.open(file_path2)
        ima = hduList[0].data
        immagine2 = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
        return immagine1, immagine2

"""
Authors
  - M. Xompero: written in 2022
"""

import numpy as _np
import os as _os
import scipy.linalg
from matplotlib import pyplot as plt
from astropy.io import fits
import pandas as pd
from skimage.measure import label, regionprops, regionprops_table
from m4.configuration import folders as fold_name
from opticalib.ground.osutils import newtn as _newtn


class ParabolaFootprintRegistration:
    """
    HOW TO USE IT::

            from m4.utils.parabola_footprint_registration import ParabolaFootprintRegistration
            pfr = ParabolaFootprintRegistration()
            cgh_on_ott, difference = pfr.main()
    """

    def __init__(self):
        """The constructor"""
        self._ottf_df = None
        self._sel_ott = None
        self._nmarkers = None
        self._grid_x = None
        self._grid_y = None
        self._params_cgh = None
        self._params_ott = None
        # results
        self._rms_ott = None
        self._rms_cgh = None
        self._rms_diff = None

    @staticmethod
    def _storageFolder():
        """Folder for data"""
        return fold_name.PARCGH_ROOT_FOLDER
        # return '/Volumes/My Passport/M4/Data/M4Data'

    def main(self):
        cgh_image, ott_image = self._readTestImages()
        cghf, ottf, cgh_not_blobs = self.marker_finder(cgh_image, ott_image)
        polycoeff = self.fit_trasformation_parameter(cghf, ottf)
        cgh_on_ott, mask_float, difference = (
            self.cgh_coord_tranform_and_trig_interpolatio(
                cgh_image, ott_image, cgh_not_blobs, polycoeff
            )
        )
        self.show_results(cgh_image, cgh_on_ott, ott_image, mask_float, difference)
        self._saveData(cgh_on_ott)
        return cgh_on_ott, difference

    def marker_finder(self, cgh_image, ott_image):
        cgh_blobs = label(cgh_image == 0)
        cgh_not_blobs = label(cgh_image != 0)
        ott_not_blobs = label(ott_image != 0)
        ott_blobs = label(ott_image == 0)
        properties = [
            "area",
            "centroid",
            "major_axis_length",
            "minor_axis_length",
            "eccentricity",
            "bbox",
        ]
        cghf_df = pd.DataFrame(regionprops_table(cgh_blobs, properties=properties))
        ottf_df = pd.DataFrame(regionprops_table(ott_blobs, properties=properties))
        self._ottf_df = ottf_df
        sel_cgh = (cghf_df["area"] < 1000) & (
            cghf_df["area"] > _np.median(cghf_df["area"]) / 2
        )
        sel_ott = (ottf_df["area"] < 1000) & (
            ottf_df["area"] > _np.median(ottf_df["area"]) / 2
        )
        self._sel_ott = sel_ott
        cghf = _np.array(
            [cghf_df["centroid-0"][sel_cgh], cghf_df["centroid-1"][sel_cgh]]
        )
        ottf = _np.array(
            [ottf_df["centroid-0"][sel_ott], ottf_df["centroid-1"][sel_ott]]
        )
        nmarkers = _np.max(cghf.shape)
        self._nmarkers = nmarkers
        print(nmarkers)
        # plot footprint
        fig, axs = plt.subplots(1, 2, figsize=(15, 5))
        for ii in _np.arange(nmarkers):
            axs[0].plot([cghf[0, ii], ottf[0, ii]], [cghf[1, ii], ottf[1, ii]])
        axs[0].title.set_text("footprint")

        # check the markers are correctly identified or reorder them
        try:
            from arte.utils.shape_fitter import ShapeFitter

            print("\nARTE Module was installed. Using it.")

            sf = ShapeFitter(cgh_not_blobs)
            sf.fit_circle_correlation()  # era sf.fit
            params_cgh = sf.parameters()
            sf = ShapeFitter(ott_not_blobs)
            sf.fit_circle_correlation()  # era sf.fit
            params_ott = sf.parameters()
        except ImportError:
            print("\nThere was no ARTE module installed. Used pre-computed values")
            params_cgh = _np.array([1022.99196625, 1023.9987704, 571.2935786])
            params_ott = _np.array([1023.88384835, 1022.79593859, 969.3681117])

        self._params_cgh = params_cgh
        self._params_ott = params_ott
        zoom = params_ott[2] / params_cgh[2]
        print("Parameters pupil cgh: ", params_cgh)
        print("Parameters pupil ott: ", params_ott)
        print("Zoom ott/cgh: ", zoom)

        from scipy.interpolate import interp2d

        npix = _np.max(cgh_image.shape)
        fit2d = interp2d(
            _np.arange(npix),
            _np.arange(npix),
            _np.array(cgh_image != 0, dtype="float"),
            kind="linear",
        )
        new_x = (_np.arange(npix) - params_cgh[0]) / zoom + params_cgh[0]
        new_y = (_np.arange(npix) - params_cgh[1]) / zoom + params_cgh[1]
        cgh_resampled = _np.array(fit2d(new_x, new_y) > 0.5, dtype="float")

        cgh_resampled_blobs = label(cgh_resampled == 0)
        properties = [
            "area",
            "centroid",
            "major_axis_length",
            "minor_axis_length",
            "eccentricity",
            "bbox",
        ]

        cghf_resampled_df = pd.DataFrame(
            regionprops_table(cgh_resampled_blobs, properties=properties)
        )
        cghf_resampled = _np.array(
            [cghf_resampled_df["centroid-0"], cghf_resampled_df["centroid-1"]]
        )
        sel_cgh_resampled = (cghf_df["area"] < 1000) & (
            cghf_df["area"] > _np.median(cghf_df["area"]) / 2
        )
        cghf_resampled = _np.array(
            [
                cghf_resampled_df["centroid-0"][sel_cgh_resampled],
                cghf_resampled_df["centroid-1"][sel_cgh_resampled],
            ]
        )
        # plot cgh reshaped
        axs[1].imshow(cgh_resampled)
        axs[1].plot(cghf_resampled[0, :], cghf_resampled[1, :], "ro")
        axs[1].title.set_text("cgh_resampled")

        order = _np.zeros(nmarkers)
        for ii in _np.arange(nmarkers):
            order[ii] = _np.argmin(
                _np.sqrt(
                    (cghf_resampled[0, :] - ottf[0, ii]) ** 2
                    + (cghf_resampled[1, :] - ottf[1, ii]) ** 2
                )
            )

        if _np.sum((_np.arange(nmarkers) - order) ** 2) < 1e-5:
            print("Check ok, no reordering need")
        else:
            print("reorder needed (not implemented yet)")
            assert False
        return cghf, ottf, cgh_not_blobs

    def fit_trasformation_parameter(self, cghf, ottf, forder=10):
        base_cgh = self._expandbase(cghf[0, :], cghf[1, :], forder=forder)
        base_ott = self._expandbase(ottf[0, :], ottf[1, :], forder=forder)

        base_cgh_plus = scipy.linalg.pinv(base_cgh)
        polycoeff = _np.matmul(ottf, base_cgh_plus)
        print(polycoeff)
        return polycoeff

    def _expandbase(self, cx, cy, forder=10):
        print("Fitting order %i" % forder)
        if forder == 3:
            zz = _np.stack((cx, cy, _np.ones(cx.size)), axis=0)
        if forder == 6:
            zz = _np.stack((cx**2, cy**2, cx * cy, cx, cy, _np.ones(cx.size)), axis=0)
        if forder == 5:
            zz = _np.stack((cx**2, cy**2, cx, cy, _np.ones(cx.size)), axis=0)
        if forder == 10:
            zz = _np.stack(
                (
                    cx**3,
                    cy**3,
                    cx**2 * cy,
                    cy**2 * cx,
                    cx**2,
                    cy**2,
                    cx * cy,
                    cx,
                    cy,
                    _np.ones(cx.size),
                ),
                axis=0,
            )

        return zz

    def cgh_coord_tranform_and_trig_interpolatio(
        self, cgh_image, ott_image, cgh_not_blobs, polycoeff, forder=10
    ):
        from scipy.interpolate import LinearNDInterpolator
        from scipy.spatial import Delaunay

        npix = _np.max(cgh_image.shape)
        grid_x, grid_y = _np.mgrid[0:npix:1, 0:npix:1]
        idw = cgh_not_blobs == 1
        coord_mat = self._expandbase(grid_x[idw], grid_y[idw], forder=forder)
        tf_coord_mat = _np.dot(_np.transpose(coord_mat), _np.transpose(polycoeff))
        tri = Delaunay(tf_coord_mat)
        fitND = LinearNDInterpolator(tri, cgh_image[idw], fill_value=0)
        cgh_on_ott = fitND(grid_x, grid_y)
        self._grid_x = grid_x
        self._grid_y = grid_y
        mask_not = ott_image == 0
        mask = ~mask_not
        mask_float = _np.array(mask, dtype="float")
        difference = (ott_image - 2 * cgh_on_ott) * mask_float * 632.8e-9
        self._rms_diff = _np.std(difference[mask])
        self._rms_ott = _np.std(ott_image[mask] * 632.8e-9)
        self._rms_cgh = _np.std(cgh_image[mask] * 632.8e-9)
        print(
            "Ott_image rms = %g, gch_image rms = %g, Difference rms = %g"
            % (self._rms_ott, self._rms_cgh, self._rms_diff)
        )
        return cgh_on_ott, mask_float, difference

    def cgh_tf(
        self, cgh_image, ott_image, cgh_not_blobs, polycoeff, forder=10, display=False
    ):
        from scipy.interpolate import LinearNDInterpolator
        from scipy.spatial import Delaunay

        npix = _np.max(cgh_image.shape)
        grid_y, grid_x = _np.mgrid[0:npix:1, 0:npix:1]
        idw = cgh_not_blobs == 1
        coord_mat = self._expandbase(grid_x[idw], grid_y[idw], forder=forder)
        tf_coord_mat = _np.matmul(_np.transpose(coord_mat), _np.transpose(polycoeff))
        tf_coord_mat = tf_coord_mat[:, [1, 0]]

        # tri = Delaunay(tf_coord_mat[:,[1,0]])
        tri = Delaunay(tf_coord_mat)
        fitND = LinearNDInterpolator(tri, cgh_image[idw], fill_value=0)
        cgh_on_ott = fitND(grid_x, grid_y)
        fitND = LinearNDInterpolator(tri, cgh_not_blobs[idw], fill_value=0)
        cgh_on_ott_mask = fitND(grid_x, grid_y) == 1

        if display:
            plt.figure()
            plt.imshow(cgh_on_ott_mask)
            plt.plot(grid_x[idw], grid_y[idw], ".r")
            plt.plot(tf_coord_mat[:, 0], tf_coord_mat[:, 1], ".g")

        self._grid_x = grid_x
        self._grid_y = grid_y
        mask_not = ott_image == 0
        mask = ~mask_not
        mask_float = _np.array(mask, dtype="float")
        difference = (ott_image - 2 * cgh_on_ott) * mask_float  # *632.8e-9
        self._rms_diff = _np.std(difference[mask])
        self._rms_ott = _np.std(ott_image[mask])  # *632.8e-9)
        self._rms_cgh = _np.std(cgh_image[mask])  # *632.8e-9)
        print(
            "Ott_image rms = %g, gch_image rms = %g, Difference rms = %g"
            % (self._rms_ott, self._rms_cgh, self._rms_diff)
        )
        return cgh_on_ott, mask_float, difference

    def image_transformation(self, cgh_image, ott_image, cghf, ottf, forder=10):
        from scipy.interpolate import LinearNDInterpolator
        from scipy.spatial import Delaunay

        polycoeff = self.fit_trasformation_parameter(cghf, ottf, forder=forder)
        npix = _np.max(cgh_image.shape)
        # npix = _np.max(ott_image.shape)  #only for test
        grid_x, grid_y = _np.mgrid[0:npix:1, 0:npix:1]
        idw = cgh_image.mask == False
        # idw = cgh_not_blobs == 1
        coord_mat = self._expandbase(grid_x[idw], grid_y[idw], forder=forder)
        tf_coord_mat = _np.matmul(_np.transpose(coord_mat), _np.transpose(polycoeff))
        tri = Delaunay(tf_coord_mat)
        fitND = LinearNDInterpolator(tri, cgh_image.data[idw], fill_value=0)
        cgh_on_ott = fitND(grid_x, grid_y)
        self._grid_x = grid_x
        self._grid_y = grid_y
        mask_not = ott_image == 0
        mask = ~mask_not
        mask = -1 * mask + 1
        mask_float = _np.array(mask, dtype="float")
        difference = (ott_image.data - 2 * cgh_on_ott) * (
            -1 * mask_float + 1
        )  # *632.8e-9
        # self._rms_diff = _np.std(difference[mask])
        # self._rms_ott = _np.std(ott_image[mask]*632.8e-9)
        # self._rms_cgh = _np.std(cgh_image[mask]*632.8e-9)
        # print("Ott_image rms = %g, gch_image rms = %g, Difference rms = %g" %(self._rms_ott, self._rms_cgh, self._rms_diff))
        cgh_tra = _np.ma.masked_array(cgh_on_ott, mask_float)
        return cgh_on_ott, mask_float, cgh_tra, difference

    def show_results(self, cgh_image, cgh_on_ott, ott_image, mask_float, difference):
        # show results
        fig, axs = plt.subplots(2, 2, figsize=(10, 10))
        axs[0, 0].imshow(cgh_image)
        axs[0, 0].title.set_text("cgh")
        axs[0, 1].imshow(cgh_on_ott * mask_float)
        axs[0, 1].title.set_text("cgh mapped on ott")
        axs[1, 0].imshow(ott_image)
        axs[1, 0].title.set_text("ott")
        axs[1, 1].imshow(difference)
        axs[1, 1].title.set_text("difference")

    def _readTestImages(self):
        """
        load data and manually adjust coordinate flip
        """
        f1 = "CCD_PAR_test_21markers_5mm_grid.txt"
        f2 = "CCD_OTT_21markers_5mm_grid.txt"

        # f1 = 'CCD_PAR_test_21markers_5mm_grid_shell_0deg.txt'

        cgh_path = _os.path.join(self._storageFolder(), f1)
        ott_path = _os.path.join(self._storageFolder(), f2)
        cgh_image = _np.loadtxt(cgh_path)
        ott_image = _np.loadtxt(ott_path)

        cgh_flip = _np.flip(_np.flip(cgh_image, axis=0), axis=1)
        return cgh_flip, ott_image

    def plot_difference_and_marker_view(self, difference):
        def crop_on_radius(npix):
            return (
                (self._grid_x - self._params_ott[0]) ** 2
                + (self._grid_y - self._params_ott[1]) ** 2
            ) < (self._params_ott[2] - npix) ** 2

        plt.figure(figsize=(10, 10))
        npix2crop = 5
        cropmask = crop_on_radius(npix2crop)
        cropmask_float = _np.array(cropmask, dtype="float")
        difference_crop = difference * cropmask_float
        plt.imshow(difference_crop)
        plt.title(
            "difference: edge crop %d, std [nm]: %g"
            % (npix2crop, _np.std(difference_crop[cropmask]))
        )
        plt.colorbar()

        nrow = _np.int32(_np.ceil(self._nmarkers / 3))
        fig, axs = plt.subplots(nrow, 3, figsize=(15, 5 * nrow))
        ottbbox = _np.array(
            [
                self._ottf_df["bbox-0"][self._sel_ott],
                self._ottf_df["bbox-1"][self._sel_ott],
                self._ottf_df["bbox-2"][self._sel_ott],
                self._ottf_df["bbox-3"][self._sel_ott],
            ]
        )
        for ii in _np.arange(self._nmarkers):
            caxs = axs[int(ii / 3), int(ii % 3)]
            caxs.imshow(
                difference_crop[
                    ottbbox[0, ii] - 50 : ottbbox[2, ii] + 50,
                    ottbbox[1, ii] - 50 : ottbbox[3, ii] + 50,
                ]
            )

    def _saveData(self, cgh_on_ott):
        tn = _newtn()
        dove = _os.path.join(self._storageFolder(), tn)

        # info
        fits_file_name = _os.path.join(dove, "final rms.txt")
        file = open(fits_file_name, "w+")
        file.write(
            " Ott_image rms = %g\n gch_image rms = %g\n Difference rms = %g"
            % (self._rms_ott, self._rms_cgh, self._rms_diff)
        )
        file.close()
        # data
        fits_file_name = _os.path.join(dove, "cgh_on_ott.fits")
        fits.writeto(fits_file_name, cgh_on_ott)

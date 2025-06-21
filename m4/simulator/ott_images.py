# -*- coding: utf-8 -*-
"""
Created on Thu May 21 18:06:09 2020

@author: Runa

The module creates simulated image of the OTT elements (M4, refmirror and parabola)
according to the current OTT configuration.
Input (loaded from create_ott) is sensitivity matrices for the various elements
(from Zemax) and mechanical dimensions of items
pixel scale is selected to have output images of 512,512import numpy as np
"""
import numpy as np
from os.path import join as _join
from m4.configuration import folders as conf
from m4.configuration.ott_parameters import Interferometer, OttParameters
from opticalib.ground import zernike, geo
from opticalib import load_fits as _lf
import matplotlib.pyplot as plt


class OttImages:
    """
    Class for creating OTT geometry and simulated images of the OTT elements
    (M4, refmirror, and parabola) according to the current OTT configuration.
    """

    def __init__(self, ott):
        """
        Initializes the OttImages class with the given OTT (simulated) system.

        Parameters
        ----------
        ott : object
            The simulated OTT object.
        """
        self._ott = ott
        self.m4offset = 0.0
        self.offset = 0.0
        self.smap = np.zeros((Interferometer.N_PIXEL[0], Interferometer.N_PIXEL[1]))
        self.rmap = np.zeros(
            (
                (2 * OttParameters.rflat_radius * OttParameters.pscale).astype(int),
                (2 * OttParameters.rflat_radius * OttParameters.pscale).astype(int),
            )
        )
        self.m4pupil = _lf(
            _join(conf.MIRROR_FOLDER, conf.mirror_conf, "m4_mech_pupil-bin2.fits")
        )
        self.m4ima = self.m4pupil * 0.0
        self.mask = _lf(_join(conf.MIRROR_FOLDER, conf.mirror_conf, "ott_mask.fits"))
        self.parmask = np.ma.make_mask(
            _lf(_join(conf.OPTICAL_FOLDER, conf.optical_conf, "ottmask.fits"))
        )
        self.zmat = _lf(_join(conf.OPTICAL_FOLDER, conf.optical_conf, "Zmat.fits"))

    def ott_smap(self, offset=None, quant=None, show: bool = False):
        """
        Creates a simulated image of the OTT elements with the given offset and
        quantization.

        Parameters
        ----------
        offset : float, optional
            The offset to be added to the simulated image.
        quant : int, optional
            The quantization level to be applied to the simulated image.
        show : int, optional
            If non-zero, displays the simulated image.

        Returns
        -------
        smap1 : numpy array
            The simulated image of the OTT elements.
        smask : numpy array
            The mask applied to the simulated image.
        """
        npix = Interferometer.N_PIXEL
        pmap = self.ott_parab_ima()
        m4map, m4smask = self.ott_m4_ima()
        mmask = m4smask * self.mask
        smap = (pmap + m4map) * mmask
        if offset is not None:
            smap = smap + offset
        roffset = (
            self._ott.referenceMirrorSlider.getPositionInM()
            - self._ott.parabolaSlider.getPositionInM()
        ) * OttParameters.pscale
        draw_ref = (
            abs(roffset) - OttParameters.rflat_radius * OttParameters.pscale
            < npix[1] / 2
        )
        if draw_ref == 1:
            rmap, rmask = self.ott_rflat_ima()
        rmap = (rmap + pmap) * rmask
        romask = geo.draw_mask(
            self.mask * 0,
            npix[0] / 2,
            npix[1] / 2 + roffset,
            OttParameters.rflat_radius * OttParameters.pscale
            - OttParameters.rflat_cell * OttParameters.pscale,
        )
        romask = romask + rmask  # ???
        irmask = -rmask + 1
        smap1 = smap * irmask + rmap
        smap1 = smap1 + self.offset
        smask = rmap + mmask
        smask[romask != 0] = 1  # ??? original: [smask !=0] = 1
        smask[romask == 1] = 0
        smap1 = smap1 * smask
        if quant == 1:
            print("add quantization")
        if show:
            plt.clf()
            plt.subplot(131)
            plt.imshow(self.ott_view())
            plt.subplot(132)
            plt.imshow(smap1, cmap="hot")
            plt.colorbar(shrink=0.3)
            plt.subplot(133)
            plt.imshow(self.pwrap(smap1, smask), cmap="gray")
            plt.tight_layout()
        return smap1, smask

    def pwrap(self, img, mask):
        """
        Applies phase wrapping to the given image.

        Parameters
        ----------
        img : numpy array
            The image to be phase wrapped.
        mask : numpy array
            The mask to be applied to the image.

        Returns
        -------
        img1 : numpy array
            The phase-wrapped image.
        """
        wav = Interferometer.WAVEL
        img1 = img.copy()
        optfact = 1
        img1[mask != 0] = np.sin(2 * np.pi * img1[mask != 0] * optfact / (wav))
        return img1

    def ott_parab_ima(self):
        """
        Creates a simulated image of the parabola element.

        Returns
        -------
        smap : numpy array
            The simulated image of the parabola element.
        """
        npix = Interferometer.N_PIXEL
        smap = (self.smap).copy()
        ww = np.dot(self.zmat, self.zmx_parpos2z())

        for i in range(0, 5):
            smap[self.parmask == True] = (
                smap[self.parmask == True]
                + ww[:, i] * (self._ott.parabola.getPositionInM())[i]
            )

        mask = geo.draw_mask(
            self.mask,
            npix[0] / 2,
            npix[1] / 2,
            OttParameters.fold_radius * OttParameters.pscale,
        )
        smap = smap * mask
        return smap

    def ott_rflat_ima(self, deshape=0):
        """
        Creates a simulated image of the reference flat element.

        Parameters
        ----------
        deshape : int, optional
            If non-zero, applies deformation to the reference flat element.

        Returns
        -------
        rmap : numpy array
            The simulated image of the reference flat element.
        rmask : numpy array
            The mask applied to the reference flat element.
        """
        npix = Interferometer.N_PIXEL
        pscale = OttParameters.pscale
        rmap = self.smap.copy()
        roffset = (
            self._ott.referenceMirrorSlider.getPositionInM()
            - self._ott.parabolaSlider.getPositionInM()
        ) * pscale
        rflat_radius = OttParameters.rflat_radius
        if abs(roffset) - rflat_radius * pscale < npix[1] / 2:
            rmask = np.ones([npix[0], npix[1]])
            rmask = geo.draw_mask(
                self.mask * 0,
                npix[0] / 2,
                npix[1] / 2 + roffset,
                rflat_radius * pscale - 1,
            )
            if deshape == 1:
                pass

            ww = np.dot(self.zmat, self.zmx_refflatpos2z())
            for i in range(0, 5):
                rmap[self.parmask == True] = (
                    rmap[self.parmask == True]
                    + ww[:, i] * (self._ott.referenceMirror.getPositionInM())[i]
                )
            rmap = rmap * rmask

        return rmap, rmask

    def ott_m4_ima(self):
        """
        Creates a simulated image of the M4 element.

        Returns
        -------
        m4img : numpy array
            The simulated image of the M4 element.
        m4smask : numpy array
            The mask applied to the M4 element.
        """
        theta = self._ott.angleRotator.getPosition() * np.pi / 180.0
        rmat = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
        ss = self.m4pupil.shape

        m4ima = self.m4ima + self.m4offset
        usepupil = self.m4pupil.copy()

        m4smask = self.ott_map2ima(usepupil)
        m4img = self.ott_map2ima(m4ima)
        mask = m4smask * self.mask

        m4smap = self.smap.copy()
        ww = np.dot(self.zmat, self.zmx_m4pos2z())
        for i in range(0, 5):
            m4smap[self.parmask == True] = (
                m4smap[self.parmask == True]
                + ww[:, i] * (self._ott.m4Exapode.getPositionInM())[i]
            )

        m4img = (m4img + m4smap) * m4smask
        return m4img, m4smask

    def ott_map2ima(self, w):
        """
        Maps the given array to the image coordinates.

        Parameters
        ----------
        w : numpy array
            The array to be mapped.

        Returns
        -------
        simg : numpy array
            The mapped array in image coordinates.
        """
        npix = Interferometer.N_PIXEL
        theta = self._ott.angleRotator.getPosition() * np.pi / 180.0
        rmat = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
        ss = w.shape
        theangle = -30.0 - self._ott.angleRotator.getPosition()
        simg = geo.rotate(w, theangle)
        parxy = self._ott.parabolaSlider.getPositionInM() * OttParameters.pscale
        x0 = np.fix(ss[0] / 2 - npix[0] / 2)
        x1 = np.fix(ss[0] / 2 + npix[0] / 2 - 1)
        y0 = np.fix(ss[1] / 2 + parxy - npix[1] / 2)
        y1 = np.fix(ss[1] / 2 + parxy + npix[1] / 2 - 1)
        x0 = x0.astype(int)
        x1 = x1.astype(int)
        y0 = y0.astype(int)
        y1 = y1.astype(int)
        simg = simg[x0 : x1 + 1, y0 : y1 + 1]
        return simg

    def ott_view(self, show: bool = False):
        """
        Creates a view of the OTT elements.

        Parameters
        ----------
        show : int, optional
            If non-zero, displays the OTT view.

        Returns
        -------
        ottimg : numpy array
            The view of the OTT elements.
        """
        m4 = self.m4pupil.copy()
        pixscale = OttParameters.PIXEL_SCALE
        parxy = [self._ott.parabolaSlider.getPositionInM() * pixscale, 0]
        refmxy = [self._ott.referenceMirrorSlider.getPositionInM() * pixscale, 0]
        ang = (-self._ott.angleRotator.getPosition()) * np.pi / 180
        rmat = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])
        parxy = rmat.dot(parxy)
        refmxy = rmat.dot(refmxy)
        ss = np.array(np.shape(m4))
        m4c = (ss - 1) / 2
        parcircle = geo.draw_mask(
            m4 * 0,
            parxy[0] + m4c[0],
            parxy[1] + m4c[1],
            OttParameters.parab_radius * pixscale,
        )
        refmcircle = geo.draw_mask(
            m4 * 0,
            refmxy[0] + m4c[0],
            refmxy[1] + m4c[1],
            OttParameters.rflat_radius * pixscale,
        )
        ottimg = m4 + parcircle + refmcircle
        if show:
            plt.imshow(ottimg)
        return ottimg

    def ott_m4view(self, show=None):
        """
        Creates a view of the M4 element.

        Parameters
        ----------
        show : int, optional
            If non-zero, displays the M4 view.

        Returns
        -------
        ottimg : numpy array
            The view of the M4 element.
        """
        from m4.devices.opt_beam import Parabola, ReferenceMirror

        par = Parabola(self._ott)
        rm = ReferenceMirror(self._ott)

        m4 = self.m4pupil.copy()
        pixscale = OttParameters.PIXEL_SCALE
        parxy = [par.trussGetPosition() * pixscale, 0]
        refmxy = [rm.rmGetPosition() * pixscale, 0]
        ang = (-self._ott.angleRotator.getPosition()) * np.pi / 180
        rmat = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])
        parxy = rmat.dot(parxy)
        refmxy = rmat.dot(refmxy)
        ss = np.array(np.shape(m4))
        m4c = (ss - 1) / 2

        parcircle = geo.draw_mask(
            m4 * 0,
            parxy[0] + m4c[0],
            parxy[1] + m4c[1],
            OttParameters.parab_radius * pixscale,
        )
        refmcircle = geo.draw_mask(
            m4 * 0,
            refmxy[0] + m4c[0],
            refmxy[1] + m4c[1],
            OttParameters.rflat_radius * pixscale,
        )

        ottimg = m4 + parcircle + refmcircle

        if show is not None:
            plt.imshow(ottimg)

        return ottimg

    def iff_images(self, zonal_modal):
        """
        Creates interferometric images based on the zonal or modal approach.

        Parameters
        ----------
        zonal_modal : int
            If 0, uses the zonal approach. If 1, uses the modal approach.

        Returns
        -------
        numpy array
            The interferometric image.
        """
        a = np.zeros((21, 21))
        b = geo.draw_mask(a, 10, 10, 10)

        segmask1 = np.ma.make_mask(self.segmask1)

        if zonal_modal == 0:
            s1 = segmask1.copy()
            s1 = s1.astype(float)
            b[segmask1 == True] = self.ifmat[3, :]
            return b
        elif zonal_modal == 1:
            comm = np.dot(self.ifmat.T, self.vmat[:, 3])
            s1 = segmask1.copy()
            s1 = s1.astype(float)
            s1[segmask1 == True] = comm
            return s1

    def _readMatFromTxt(self, file_name):
        """
        Reads a matrix from a text file.

        Parameters
        ----------
        file_name : str
            The path to the text file containing the matrix.

        Returns
        -------
        mat : numpy array
            The matrix read from the text file.
        """
        file = open(file_name, "r")
        triplets = file.read().split()
        x = np.array(triplets)
        mat = x.reshape(11, 6)
        return mat.astype(float)

    def zmx_parpos2z(self):
        """
        Reads the matrix for parable positions to Zernike coefficients.

        Returns
        -------
        mat : numpy array
            The matrix for parable positions to Zernike coefficients.
        """
        file_name = _join(conf.OPTICAL_FOLDER, conf.optical_conf, "PAR_pos2z.txt")
        mat = self._readMatFromTxt(file_name)
        return mat

    def zmx_refflatpos2z(self):
        """
        Reads the matrix for reference flat positions to Zernike coefficients.

        Returns
        -------
        mat : numpy array
            The matrix for reference flat positions to Zernike coefficients.
        """
        file_name = _join(conf.OPTICAL_FOLDER, conf.optical_conf, "M4_pos2z.txt")
        mat = self._readMatFromTxt(file_name)
        return mat

    def zmx_m4pos2z(self):
        """
        Reads the matrix for deformable mirror positions to Zernike coefficients.

        Returns
        -------
        mat : numpy array
            The matrix for deformable mirror positions to Zernike coefficients.
        """
        file_name = _join(conf.OPTICAL_FOLDER, conf.optical_conf, "M4_pos2z.txt")
        mat = self._readMatFromTxt(file_name)
        return mat

    def create_zmat(self, file_name):
        """
        Creates a Zernike matrix from a FITS file.

        Parameters
        ----------
        file_name : str
            The path to the FITS file.

        Returns
        -------
        zmat : numpy array
            The Zernike matrix.
        """
        load = _lf(file_name)
        final_mask = np.invert(load.astype(bool))
        prova = np.ma.masked_array(
            np.ones(
                (
                    (2 * OttParameters.parab_radius * OttParameters.pscale).astype(int),
                    (2 * OttParameters.parab_radius * OttParameters.pscale).astype(int),
                )
            ),
            mask=final_mask,
        )
        zernike_mode = np.arange(10) + 1
        mm = np.where(final_mask == False)
        _, _, _, xx, yy = geo.qpupil(final_mask)
        zmat = zernike._getZernike(xx[mm], yy[mm], zernike_mode)
        return zmat


# class _Old_OttImages:
#     """
#         Class for creating OTT geometry

#     HOW TO USE IT::
#         from m4.ott_sim.ott_images import OttImages
#         geo = OttImages(ott)
#         imshow(geo.ott_view())
#     """

#     def __init__(self, ott):
#         """The constructor"""
#         self._ott = ott
#         self.m4offset = 0.0
#         self.offset = 0.0
#         self.smap = np.zeros((Interferometer.N_PIXEL[0], Interferometer.N_PIXEL[1]))
#         self.rmap = np.zeros(
#             (
#                 (2 * OttParameters.rflat_radius * OttParameters.pscale).astype(int),
#                 (2 * OttParameters.rflat_radius * OttParameters.pscale).astype(int),
#             )
#         )
#         self.m4pupil = read_data.readFits_data(
#             os.path.join(
#                 conf.MIRROR_FOLDER, conf.mirror_conf, "m4_mech_pupil-bin2.fits"
#             )
#         )
#         self.m4ima = self.m4pupil * 0.0
#         self.mask = read_data.readFits_data(
#             os.path.join(conf.MIRROR_FOLDER, conf.mirror_conf, "ott_mask.fits")
#         )
#         self.parmask = np.ma.make_mask(
#             read_data.readFits_data(
#                 os.path.join(conf.OPTICAL_FOLDER, conf.optical_conf, "ottmask.fits")
#             )
#         )
#         #         self.segmask1 = np.ma.make_mask(obj.readFits_object(
#         #                                         os.path.join(conf.path_name.MIRROR_FOLDER,
#         #                                         conf.mirror_conf, 'py-sect4-mask.fits')))
#         self.zmat = read_data.readFits_data(
#             os.path.join(conf.OPTICAL_FOLDER, conf.optical_conf, "Zmat.fits")
#         )

#     ### Functions to create images ###
#     def ott_smap(self, offset=None, quant=None, show=0):
#         npix = Interferometer.N_PIXEL
#         pmap = self.ott_parab_ima()
#         m4map, m4smask = self.ott_m4_ima()
#         mmask = m4smask * self.mask
#         smap = (pmap + m4map) * mmask
#         if offset is not None:
#             smap = smap + offset

#         # fullmask  = maskm
#         roffset = (
#             self._ott.referenceMirrorSlider.getPositionInM()
#             - self._ott.parabolaSlider.getPositionInM()
#         ) * OttParameters.pscale
#         draw_ref = (
#             abs(roffset) - OttParameters.rflat_radius * OttParameters.pscale
#             < npix[1] / 2
#         )
#         if draw_ref == 1:
#             rmap, rmask = self.ott_rflat_ima()
#             # smap[ott.idx[idr]] = rmap[ott.idx[idr]]+pmap[ott.idx[idr]]
#             # smap[idc] = 0
#             # fullmask = maskm+maskr
#             # fullmask[np.where(fullmask)] = 1
#             # fullmask[idc] = 0

#         rmap = (rmap + pmap) * rmask
#         romask = geo.draw_mask(
#             self.mask * 0,
#             npix[0] / 2,
#             npix[1] / 2 + roffset,
#             OttParameters.rflat_radius * OttParameters.pscale
#             - OttParameters.rflat_cell * OttParameters.pscale,
#         )
#         romask = romask + rmask
#         irmask = -rmask + 1

#         smap1 = smap * irmask + rmap
#         smap1 = smap1 + self.offset  # this is to be verified
#         smask = rmap + mmask
#         smask[smask != 0] = 1
#         smask[romask == 1] = 0
#         smap1 = smap1 * smask
#         if quant == 1:
#             print("add quantization")
#         # ToDo smap = add_quantization(smap, bad=fullmask, quant=ott.interf.quantization):

#         # mask = fullmask
#         # cir = qpupil(fullmask, xx=xx, yy=yy,nocircle=1)
#         """
#         ToDo
#         if keyword_set(show) then begin
#           smap1 = smap

#           out = where(fullmask eq 0)
#           smap1[out] = max(smap1[where(fullmask)])
#           wset, 0
#           image_show, /As,  smap1;, titl='WF (?) map'
#           wset, 1
#           gg = ott_geometry(/sh)
#           wset, 3
#           display_data,  mirror.gpos[*], mirror.coord, spot_mag=2, /no_n,/Sh, back=max(mirror.gpos);image_Show, /As, /sh, mirror.m4ima
#           if keyword_set(fringes) then begin
#             wset, 2
#             loadct, 0,/silent
#             interf = pwrap(smap, lambda=ott.interf.lambda, optfact=1, bad=mask, detector_mask=detmask)*fullmask

#             out = where(fullmask eq 0)

#             interf[0:ott.interf.npix[0]/5,4.5*ott.interf.npix[1]/5:ott.interf.npix[1]-1]=-1
#             interf[4*ott.interf.npix[0]/5:ott.interf.npix[0]-1,4.5*ott.interf.npix[1]/5:ott.interf.npix[1]-1]=1

#             image_show, /As, interf, min_v=-1, max_v=1
#             loadct, 3,/silent
#           endif
#       endif
#     """
#         if show != 0:
#             plt.clf()
#             plt.subplot(131)
#             plt.imshow(self.ott_view())
#             plt.subplot(132)
#             plt.imshow(smap1, cmap="hot")
#             plt.colorbar()
#             plt.subplot(133)
#             plt.imshow(self.pwrap(smap1, smask), cmap="gray")
#         return smap1, smask

#     def pwrap(self, img, mask):
#         wav = Interferometer.WAVEL
#         img1 = img.copy()
#         optfact = 1
#         img1[mask != 0] = np.sin(2 * np.pi * img1[mask != 0] * optfact / (wav))
#         return img1

#     def ott_parab_ima(self):
#         npix = Interferometer.N_PIXEL
#         smap = (self.smap).copy()
#         ww = np.dot(self.zmat, self.zmx_parpos2z())

#         for i in range(0, 5):
#             smap[self.parmask == True] = (
#                 smap[self.parmask == True]
#                 + ww[:, i] * (self._ott.parabola.getPositionInM())[i]
#             )

#         mask = geo.draw_mask(
#             self.mask,
#             npix[0] / 2,
#             npix[1] / 2,
#             OttParameters.fold_radius * OttParameters.pscale,
#         )
#         smap = smap * mask
#         return smap

#     def ott_rflat_ima(self, deshape=0):
#         npix = Interferometer.N_PIXEL
#         pscale = OttParameters.pscale
#         rmap = self.smap.copy()
#         roffset = (
#             self._ott.referenceMirrorSlider.getPositionInM()
#             - self._ott.parabolaSlider.getPositionInM()
#         ) * pscale
#         rflat_radius = OttParameters.rflat_radius
#         if abs(roffset) - rflat_radius * pscale < npix[1] / 2:
#             rmask = np.ones([npix[0], npix[1]])
#             rmask = geo.draw_mask(
#                 self.mask * 0,
#                 npix[0] / 2,
#                 npix[1] / 2 + roffset,
#                 rflat_radius * pscale - 1,
#             )
#             if deshape == 1:
#                 # rima = read_refflat_deformation(1.)
#                 # rmap[npix[0]/2-rflat_radius*pscale-1, npix[1]/2-rflat_radius*pscale-1] = rima
#                 pass

#             # idr   = where(rmask[ott.idx])
#             ww = np.dot(self.zmat, self.zmx_refflatpos2z())
#             for i in range(0, 5):
#                 rmap[self.parmask == True] = (
#                     rmap[self.parmask == True]
#                     + ww[:, i] * (self._ott.referenceMirror.getPositionInM())[i]
#                 )
#                 # rmap[ott.idx[idr]] += ww[i,idr]* ott.refflat[i]

#             rmap = rmap * rmask
#             # cm = edge_erode(erode=2, rmask)
#             # idc = where((cm+rmask) eq 1)

#         return rmap, rmask

#     def ott_m4_ima(self):  # bad=mask, idm=idm, idtot=idtot, optpupil=optpupil

#         theta = self._ott.angleRotator.getPosition() * np.pi / 180.0
#         rmat = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
#         ss = self.m4pupil.shape

#         m4ima = self.m4ima + self.m4offset
#         usepupil = self.m4pupil.copy()
#         # if keyword_set(optpupil) then usepupil = usepupil*ott.m4optpupil   ;mod 20161010

#         m4smask = self.ott_map2ima(usepupil)
#         m4img = self.ott_map2ima(m4ima)
#         mask = m4smask * self.mask
#         # idtot = where(mask[ott.idx])
#         m4smap = self.smap.copy()
#         ww = np.dot(self.zmat, self.zmx_m4pos2z())
#         for i in range(0, 5):
#             m4smap[self.parmask == True] = (
#                 m4smap[self.parmask == True]
#                 + ww[:, i] * (self._ott.m4Exapode.getPositionInM())[i]
#             )
#             # m4smap[ott.idx[idtot]] = m4smap[ott.idx[idtot]] + ww[i,idtot]* ott.m4[i]

#         # m4img[ott.idx[idtot]] += m4smap[ott.idx[idtot]]
#         m4img = (m4img + m4smap) * m4smask
#         return m4img, m4smask

#     def ott_map2ima(self, w):  # debugged
#         npix = Interferometer.N_PIXEL
#         theta = self._ott.angleRotator.getPosition() * np.pi / 180.0
#         rmat = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
#         ss = w.shape
#         theangle = -30.0 - self._ott.angleRotator.getPosition()
#         simg = geo.rotate(w, theangle)
#         parxy = self._ott.parabolaSlider.getPositionInM() * OttParameters.pscale
#         x0 = np.fix(ss[0] / 2 - npix[0] / 2)
#         x1 = np.fix(ss[0] / 2 + npix[0] / 2 - 1)
#         y0 = np.fix(ss[1] / 2 + parxy - npix[1] / 2)
#         y1 = np.fix(ss[1] / 2 + parxy + npix[1] / 2 - 1)
#         x0 = x0.astype(int)
#         x1 = x1.astype(int)
#         y0 = y0.astype(int)
#         y1 = y1.astype(int)
#         simg = simg[x0 : x1 + 1, y0 : y1 + 1]
#         return simg

#     def ott_view(self, show=None):
#         # pixscale = 200. #pix/m
#         # parod = 1.44
#         # rmod = 0.6
#         m4 = self.m4pupil.copy()
#         pixscale = OttParameters.PIXEL_SCALE
#         parxy = [self._ott.parabolaSlider.getPositionInM() * pixscale, 0]
#         refmxy = [self._ott.referenceMirrorSlider.getPositionInM() * pixscale, 0]
#         # ang = (-30-self._ott.angleRotator.getPosition())*np.pi/180
#         ang = (-self._ott.angleRotator.getPosition()) * np.pi / 180
#         rmat = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])
#         parxy = rmat.dot(parxy)
#         refmxy = rmat.dot(refmxy)
#         ss = np.array(np.shape(m4))
#         m4c = (ss - 1) / 2
#         parcircle = geo.draw_mask(
#             m4 * 0,
#             parxy[0] + m4c[0],
#             parxy[1] + m4c[1],
#             OttParameters.parab_radius * pixscale,
#         )
#         refmcircle = geo.draw_mask(
#             m4 * 0,
#             refmxy[0] + m4c[0],
#             refmxy[1] + m4c[1],
#             OttParameters.rflat_radius * pixscale,
#         )
#         ottimg = m4 + parcircle + refmcircle
#         if show is not None:
#             plt.imshow(ottimg)

#         return ottimg

#     def ott_m4view(self, show=None):
#         from m4.devices.opt_beam import Parabola, ReferenceMirror
#         import os

#         par = Parabola(self._ott)
#         rm = ReferenceMirror(self._ott)

#         m4 = self.m4pupil.copy()
#         pixscale = OttParameters.PIXEL_SCALE
#         parxy = [par.trussGetPosition() * pixscale, 0]
#         refmxy = [rm.rmGetPosition() * pixscale, 0]
#         ang = (-self._ott.angleRotator.getPosition()) * np.pi / 180
#         rmat = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])
#         parxy = rmat.dot(parxy)
#         refmxy = rmat.dot(refmxy)
#         ss = np.array(np.shape(m4))
#         m4c = (ss - 1) / 2

#         parcircle = geo.draw_mask(
#             m4 * 0,
#             parxy[0] + m4c[0],
#             parxy[1] + m4c[1],
#             OttParameters.parab_radius * pixscale,
#         )
#         refmcircle = geo.draw_mask(
#             m4 * 0,
#             refmxy[0] + m4c[0],
#             refmxy[1] + m4c[1],
#             OttParameters.rflat_radius * pixscale,
#         )

#         ottimg = m4 + parcircle + refmcircle

#         if show is not None:
#             plt.imshow(ottimg)

#         return ottimg

#     def iff_images(self, zonal_modal):
#         a = np.zeros((21, 21))
#         b = geo.draw_mask(a, 10, 10, 10)

#         segmask1 = np.ma.make_mask(self.segmask1)

#         if zonal_modal == 0:
#             s1 = segmask1.copy()
#             s1 = s1.astype(float)
#             b[segmask1 == True] = self.ifmat[3, :]
#             return b
#         elif zonal_modal == 1:
#             comm = np.dot(self.ifmat.T, self.vmat[:, 3])
#             s1 = segmask1.copy()
#             s1 = s1.astype(float)
#             s1[segmask1 == True] = comm
#             return s1

#     ### Sensitivity matrices ###
#     def _readMatFromTxt(self, file_name):
#         """Function to read matrix of 11 Zernike x 6 displacements,
#         m RMS, per 1 m displacement - or 1 radiant rotation

#         Parameters
#         ----------
#         file_name: string
#                     matrix file path

#         Returns
#         -------
#                 mat: numpy array [11,6]
#                     matrix from txt file
#         """
#         file = open(file_name, "r")
#         triplets = file.read().split()
#         x = np.array(triplets)
#         mat = x.reshape(11, 6)
#         return mat.astype(float)

#     def zmx_parpos2z(self):
#         """
#         Returns
#         -------
#                 mat: numpy array [11,6]
#                     matrix parable positions to zernike
#         """
#         file_name = os.path.join(
#             conf.OPTICAL_FOLDER, conf.optical_conf, "PAR_pos2z.txt"
#         )
#         # file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_PAR_pos2z.txt'
#         mat = self._readMatFromTxt(file_name)
#         return mat

#     def zmx_refflatpos2z(self):
#         """
#         Returns
#         -------
#                 mat: numpy array [11,6]
#                     matrix reference flat positions to zernike
#         """
#         file_name = os.path.join(conf.OPTICAL_FOLDER, conf.optical_conf, "M4_pos2z.txt")
#         # file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_FM_pos2z.txt'
#         mat = self._readMatFromTxt(file_name)
#         return mat

#     def zmx_m4pos2z(self):
#         """
#         Returns
#         -------
#                 mat: numpy array [11,6]
#                     matrix deformable mirror positions to zernike
#         """
#         file_name = os.path.join(conf.OPTICAL_FOLDER, conf.optical_conf, "M4_pos2z.txt")
#         # file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_M4_pos2z.txt'
#         mat = self._readMatFromTxt(file_name)
#         return mat

#     # Zmat
#     def create_zmat(self, file_name):
#         """
#         Returns
#         -------
#             zmat: numpy array

#         """
#         # file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ott_mask.fits'
#         hduList = pyfits.open(file_name)
#         final_mask = np.invert(hduList[0].data.astype(bool))

#         prova = np.ma.masked_array(
#             np.ones(
#                 (
#                     (2 * OttParameters.parab_radius * OttParameters.pscale).astype(int),
#                     (2 * OttParameters.parab_radius * OttParameters.pscale).astype(int),
#                 )
#             ),
#             mask=final_mask,
#         )
#         zernike_mode = np.arange(10) + 1
#         mm = np.where(final_mask == False)
#         x, y, r, xx, yy = geo.qpupil(final_mask)
#         zmat = zernike._getZernike(xx[mm], yy[mm], zernike_mode)
#         return zmat

"""
Authors
  - L. Oggioni
  - 20221216
"""

import glob
import os
import numpy as np
from matplotlib import pyplot as plt
from m4.utils import osutils as osu
from m4.ground import zernike


class PistMeas:
    """
    HOW TO USE IT::

    from m4.utils import PistonMeasures
    from m4.configuration import start

    ## conf='---/MyConf.yaml'
    conf='G:\il mio Drive\Lavoro_brera\M4\LucaConfig.yaml'

    ott, interf, dm = start.create_ott(conf)
    pm=PistonMeasures.PistMeas(ott, interf, dm)

    a='G:/Il mio Drive/Lavoro_brera/PROGETTI/M4/Data'
    tn='20210327_151400_noise'

    """

    def __init__(self, ott, interf, dm):
        """The constructor"""
        self._ott = ott
        self._interf = interf
        self._dm = dm

    def main(self):
        a = "G:/Il mio Drive/Lavoro_brera/PROGETTI/M4/Data"
        tn = "20210118_233600_noise"
        c = self.load_datacube(a, tn, Nmax=500)
        cc = self.remZern_datacube(c, modes=np.array([1, 2, 3, 4]))
        return c, cc

    def load_datacube(
        self,
        a="G:/Il mio Drive/Lavoro_brera/PROGETTI/M4/Data",
        tn="20210327_151400_noise",
        Nmax=0,
    ):
        """
        load a datacube of images acquired with the interferometers
        Parameters
        ----------
            a: string
                path

            tn: string
                file name
        """
        fl = self.fileList(tn, a)
        if Nmax > 0 & Nmax < len(fl):
            fl = fl[0:Nmax]
        imgcube = osu.createCube(fl)
        return imgcube

    def remZern_datacube(self, imgcube, modes=np.array([1, 2, 3, 4])):
        """
        Remove the chosen zernike from the images in the datacube,
        plot the removed coefficients

        Parameters
        ----------
            imgcube: cube of images

            mode: numpy array
                zernike coefficient to be removed
        """
        cube = np.zeros(np.shape(imgcube))
        coeff = np.zeros([len(cube), len(modes)])
        for j in range(len(imgcube)):
            cube[j, :, :], coeff[j, :] = self.removeZernike(imgcube[j, :, :], modes)
        plt.figure()
        for j in range(len(modes)):
            plt.plot(coeff[:, j])
        return cube, coeff

    def showMovie(self, imgcube, diff=0):
        plt.figure()
        for j in range(len(imgcube)):
            plt.pause(0.1)
            plt.show()
            plt.imshow(imgcube[j, :, :])

    def script_singlezone(self, imgcube, N_diffalg=5, mask=None):
        """
        script for evaluate the performance of the unwrapping algorithm and the differential algorithm.
        The piston value of each image is averaged on the full aperture

        Parameters
        ----------
            imgcube: cube of images

            N_diffalg: int
                        number of point for the differential algorithm
        """
        f, ([ax1, ax2], [ax3, ax4]) = plt.subplots(2, 2)
        f.tight_layout()
        if mask is None:
            pist = np.ma.mean(imgcube[:, :, :], axis=1)
            pist = np.ma.mean(pist, axis=1)
        else:
            pist = np.zeros(imgcube.shape[0])
            for ii in np.arange(1, imgcube.shape[0]):
                c = imgcube[ii, :, :].data
                temp = c[mask]
                pist[ii] = np.ma.mean(temp)
        thr = 632e-9 / 4
        pu = self.signal_unwrap(pist, thr)
        t = np.arange(len(pist)) / 28.57
        marker_style1 = dict(color="tab:blue", linestyle=":", marker="o")
        marker_style2 = dict(color="tab:red", linestyle=":", marker="*")
        ax1.plot(t, pist - pist[0], **marker_style1)
        ax1.plot(t, pu, **marker_style2)
        ax1.set_title("Piston unwrapping")
        ax1.set_xlabel("time (s)")
        ax1.set_ylabel("Piston (m)")
        dd = pist[1:] - pist[0:-1]
        ax2.plot(t[0:-1], dd, **marker_style1)
        ax2.plot(t[0:-1], np.zeros(len(dd)) + thr)
        ax2.plot(t[0:-1], np.zeros(len(dd)) - thr)
        ax2.set_title("differential vector")
        ax2.set_xlabel("time (s)")
        ax2.set_ylabel("Delta Piston (m)")
        template = np.array([1, 3, 5, 7, 9, 11, 13, 15])
        m, sig = self.find_best_diffalg(pu, template)
        ax3.plot(template, m, "-o")
        ax3.set_title("DiffAlg, find param: average")
        ax3.set_xlabel("Length DiffAlg")
        ax3.set_ylabel("Average (m)")
        ax4.plot(template, sig, "-o")
        ax4.set_title("DiffAlg, find param: std")
        ax4.set_xlabel("Length DiffAlg")
        ax4.set_ylabel("STD (m)")
        ##
        f2, (f2ax1, f2ax2, f2ax3) = plt.subplots(3, 1)
        f2.tight_layout()
        f2ax1.plot(t, pist - pist[0], **marker_style1)
        f2ax1.set_xlabel("time (s)")
        f2ax1.set_ylabel("Piston (m)")
        f2ax1.set_title("raw data")
        f2ax2.plot(t[0:-1], dd, **marker_style1)
        f2ax2.plot(t[0:-1], np.zeros(len(dd)) + thr)
        f2ax2.plot(t[0:-1], np.zeros(len(dd)) - thr)
        f2ax2.set_title("differential vector")
        f2ax2.set_xlabel("time (s)")
        f2ax2.set_ylabel("Delta Piston (m)")
        f2ax3.plot(t, pu, **marker_style1)
        f2ax3.set_xlabel("time (s)")
        f2ax3.set_ylabel("Piston (m)")
        f2ax3.set_title("unwrapped data")
        ##
        f3 = plt.figure()
        s = self.diffalg(pu, N_diffalg)
        plt.plot(s, "-o")
        np.disp(s.std())
        plt.title("DiffAlg with length %i" % N_diffalg)
        plt.xlabel("")
        plt.ylabel("m")

    # def script_confronto_zone(self,imgcube):
    #
    #     #close('all')
    #     #tn='20210118_233600_noise'
    #
    #     mx=np.array([245,245,300,200])
    #     my=np.array([190,290,240,240])
    #
    #
    #     L=5
    #
    #     plt.imshow(imgcube[1,:,:])
    #     ax = plt.gca()
    #
    #     plt.figure(5)
    #     (ax1, ax2) = plt.subplots(2,1)
    #
    #     for j in range(len(mx)):
    #
    #         disp(j)
    #
    #         rect = matplt.plotlib.patches.Rectangle((mx[j]-L, my[j]-L), 2*L, 2*L, linewidth=1, edgecolor='r', facecolor='none')
    #         ax.add_patch(rect)
    #         plt.show()
    #
    #         pist=np.ma.mean(imgcube[:,mx[j]-L:mx[j]+L,my[j]-L:my[j]+L], axis=1)
    #         pist=np.ma.mean(pist, axis=1)
    #
    #         thr=632e-9/4
    #         pu = self.signal_unwrap(pist,thr)
    #         t=np.arange(len(pist))/28.57
    #
    #         plt.figure(2)
    #         plt.plot(t,pist-pist[0],'-o')
    #
    #         plt.figure(3)
    #         plt.plot(t,pu,'-o')
    #
    #         plt.figure(4)
    #         dd= pist[1:]-pist[0:-1]
    #         plt.plot(t[0:-1],dd,'-o')
    #         plt.plot(t[0:-1],np.zeros(len(dd))+thr)
    #         plt.plot(t[0:-1],np.zeros(len(dd))-thr)
    #
    #
    #         template=np.array([1,3,5,7,11,13,15])
    #         m,sig=self.find_best_diffalg(pu,template)
    #         ax1.plot(template,m,'-o')
    #         ax2.plot(template,sig,'-o')
    #
    #
    #         plt.figure(6)
    #         s=self.diffalg(pu,N=5)
    #         plt.plot(s,'-o')
    #
    #        plt.figure(5)
    #        dd= pist[1:]-pist[0:-1]
    #        plt.plot(t[0:-1],np.absolute(dd.data),'-o')
    #        plt.plot(t[0:-1],np.zeros(len(dd))+thr)
    def signal_unwrap(self, x, thr=632e-9 / 4, phase=632e-9 / 2):
        """
        unwrap the signal coming from an interferometric measure

        Parameters
        ----------
            x: numpy array
                vector of values to be unwrupped
            thr: float
                threshold in m
            phase: float
                    step added by the interferometer, in m
        """
        v = x - x[0]
        npx = np.size(v)
        for i in np.arange(1, npx):
            dv = v[i] - v[i - 1]
            if dv > thr:
                v[i] = v[i] - np.abs(phase * (round(dv / phase)))
            if dv < -thr:
                v[i] = v[i] + np.abs(phase * (round(dv / phase)))
        return v

    def signal_unwrap_single(self, x, thr=632e-9 / 4, phase=632e-9 / 4):
        v = x.copy()
        v = v - v[1:]
        npx = np.size(v)
        for i in np.arange(1, npx):
            # dv = v[i]-v[i-1]
            if v[i] > thr:
                v[i] = v[i] - np.abs(phase * (round(v[i] / phase)))
            if v[i] < -thr:
                v[i] = v[i] + np.abs(phase * (round(v[i] / phase)))
        return v

    def diffalg(self, x, N=3):
        """
        differential algorithm, used to filter out low frequency noise

        Parameters
        ----------
            x: numpy array
                vector of values to be filtered
            N: integer
                must be odd
        """
        dim = int(np.floor(len(x) / N))
        s = np.zeros(dim)
        for j in range(dim):
            for k in range(N - 1):
                s[j] = s[j] + ((x[j * N + k + 1] - x[j * N + k]) * (-1) ** k) / 2
            s[j] = s[j] / (N - 1)
        return s

    def find_best_diffalg(self, x, template):
        """
        Apply the differential algorithm with different N specified in the template to find the best one

        Parameters
        ----------
            x: numpy array
                vector of values to be filtered
            template: numpy array
                the elements must be odd
        """
        m = np.zeros(len(template))
        sig = np.zeros(len(template))
        for i in range(len(template)):
            s = self.diffalg(x, template[i])
            m[i] = np.ma.mean(s)
            sig[i] = np.ma.std(s)
        return m, sig

    def fileList(self, tn, fold=None):
        """
        Parameters
        ----------
        tn: string
            tracking number where to search for the images file list

        Returns
        -------
        lsdir: the list of image files
        fold1: the path where the file is found
        """
        if fold is not None:
            name = "*.h5"
            addfold = "/hdf5/"
        else:
            folds = osu.findTracknum(tn)
            addfold = "/"
            name = "20*"
            for folder in folds:
                if folder == "OPDImages":
                    addfold = "/hdf5/"
                    name = "img*"
                    fold = folder
        fold1 = fold + "/" + tn + addfold  # to be re-checked at OTT!!
        lsdir = sorted(
            glob.glob(fold1 + name),
            key=lambda x: (os.path.basename(x).split("/")[-1].split(".")[0]),
        )
        # lsdir = lsdir[0]
        return lsdir

    def removeZernike(self, ima, modes=np.array([1, 2, 3, 4])):
        """
        remove the zernike from an image

        """

        coeff, mat = zernike.zernikeFit(ima, modes)
        surf = zernike.zernikeSurface(ima, coeff, mat)
        new_ima = ima - surf
        return new_ima, coeff

    def create_circular_mask(self, h, w, center=None, radius=None):
        if center is None:  # use the middle of the image
            center = (int(w / 2), int(h / 2))
        if (
            radius is None
        ):  # use the smallest distance between the center and image walls
            radius = min(center[0], center[1], w - center[0], h - center[1])
        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)
        mask = dist_from_center <= radius
        return mask

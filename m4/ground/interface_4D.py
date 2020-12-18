# -*- coding: utf-8 -*-
"""
Authors
  - C. Selmi:  written in March 2020

Function for simulation or data acquisition from the interferometer::

    from m4.configuration import start
    ott = start.create_ott()
    from m4.utils.interface_4D import comm4d
    c4d = comm4d()
    image = c4d.acq4d(ott, nframes, 0)
    
"""
import os
import numpy as np
from matplotlib import pyplot as plt
# a='D:\Astro\ARCETRI\Python\M4-master'
from astropy.io import fits as pyfits
from m4.configuration import config as conf
from m4.ott_sim.ott_images import OttImages
from m4.ground.read_data import InterferometerConverter
from m4.configuration.config import fold_name


class comm4d():

    phcam_ip = ''
    phcam_datapath = ''
    phcam_conf = ''

    def acq4d(self, ott, nframes=1, show=0):
        """
        Parameters
        ----------
            ott: object
                the ott object
            nframes: int
                number of frames
            show: int
                0 to not show the image

        Returns
        -------
            masked_ima: numpy masked array
                    interferometer image
        """
        ottIma = OttImages(ott)

        if conf.simulated ==1:
            opd, mask = ottIma.ott_smap(show=show)
            masked_ima = np.ma.masked_array(opd.T, mask=np.invert(mask.astype(bool)).T)

        else:
            print('Frame acquisition')
            from oaautils import i4d
            interf = i4d.I4D()
            if nframes == 1:
                masked_ima = self._getMeasurementOnTheFly(interf)
            else:
                cube_images = None
                for i in range(nframes):
                    ima = self._getMeasurementOnTheFly(interf)
                    if cube_images is None:
                        cube_images = ima
                    else:
                        cube_images = np.ma.dstack((cube_images, ima))
                masked_ima = np.ma.mean(cube_images, axis=2)

            if show==0:
                pass
            else:
                plt.clf()
                plt.imshow(masked_ima)
                plt.colorbar()

        return masked_ima

    def _getMeasurementOnTheFly(self, interf):

        def _measureAndStoreH5(interf, filename):
            import time
            import shutil
            nMeasure=1
            interf.connect()
            interf.capture(1, name= 'DM_temp')
            interf.produce('DM_temp')
            interf.disconnect()
            time.sleep(1.0)
            fName = os.path.join(fold_name.PHASECAM_ROOT_FOLDER, 'DM_temp')
            #fName = '/home/m4/4d/Zcopy/DM_temp'

            for i in range(nMeasure):
                shutil.move(fName + '/hdf5/img_%04d.h5' %i,
                        filename + "_m%02d" %i + ".h5")

            shutil.rmtree(fName + '/hdf5')
            shutil.rmtree(fName + '/raw')


        _measureAndStoreH5(interf, '/tmp/prova4d')
        return InterferometerConverter.from4D('/tmp/prova4d_m00.h5')


    def save_phasemap(self, dove, name, image):
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, image.data)
        pyfits.append(fits_file_name, image.mask.astype(int))

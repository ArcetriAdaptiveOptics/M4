# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 06:20:03 2020

@author: Runa
"""
# from importlib import reload
# import os
# from matplotlib import pyplot as plt
# a='D:\Astro\ARCETRI\Python\M4-master'
# os.chdir(a)
# from m4.configuration.create_ott import *
# from m4.configuration import start
# ott=start.create_ott()
#from m4.ott_sim.ott_images import *
from m4.configuration import config as conf
from m4.configuration.ott_parameters import *
from m4.ott_sim.ott_images import OttImages
from m4.ground.interferometer_converter import InterferometerConverter


class comm4d():

    phcam_ip = ''
    phcam_datapath = ''
    phcam_conf = ''

    def acq4d(self, ott, nframes=1, show=0):
        ottIma = OttImages(ott)

        if conf.simulated ==1:
            opd, mask = ottIma.ott_smap(show=show)
            masked_ima = np.ma.masked_array(opd.T, mask=np.invert(mask.astype(bool)).T)

        else:
            print('some function to acquire the interferometer....')
            from oaautils import i4d
            interf = i4d.I4D()
            masked_ima = self._getMeasurementOnTheFly(interf)

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
            fName = '/home/labot/4d/Zcopy/DM_temp'
    
            for i in range(nMeasure):
                shutil.move(fName + '/hdf5/img_%04d.h5' %i,
                        filename + "_m%02d" %i + ".h5")
    
            shutil.rmtree(fName + '/hdf5')
            shutil.rmtree(fName + '/raw')
    
    
        _measureAndStoreH5(interf, '/tmp/prova4d')
        return InterferometerConverter.from4D('/tmp/prova4d_m00.h5')

'''
Authors
  - C. Selmi: written in 2019
'''

import os
import time
import glob
import numpy as np
from astropy.io import fits as pyfits
from m4.ground.read_data import InterferometerConverter
ic = InterferometerConverter()

path2 = '/home/labot/testdir/' #SSD SDA
path0 = '/mnt/data/M4/Data/testdir/' #HDD
path1 = '/media/labot/4799bd14-ebd0-4072-b5b1-5dae7c3bbd8b/data/testdir/' #NVME
path1 = '/mnt/testlabot/data/testdir/' #NVME

#mypath = '/Users/rm/Desktop/Arcetri/M4/4D/files4D'
path_list = [path0, path1, path2]
#path_list = [mypath, mypath, mypath]

cube_name = 'Cube.fits'

def read_data_from_differents_paths(n_frames):
    dt_list = []
    for path in path_list:
        list = glob.glob(os.path.join(path,'*.4D'))
        t0 = time.time()
        for i in range(n_frames):
            ima = ic.fromNew4D(list[i])
        t1 = time.time()
        dt_list.append(t1-t0)

    print('HDD, NVMe, SSD, %d files read' %n_frames)
    print(dt_list)

def create_test_cubes(n_frames):
    cube = None
    dt_list = []
    for path in path_list:
        t0 = time.time()
        list = glob.glob(os.path.join(path,'*.4D'))
        for i in range(n_frames):
            ima = ic.fromNew4D(list[i])
            if cube is None:
                cube = ima
            else:
                cube = np.ma.dstack((cube, ima))
        file_name = os.path.join(path, cube_name)
        pyfits.writeto(file_name, cube.data, overwrite=True)
        pyfits.append(file_name, cube.mask.astype(int), overwrite=True)
        t1 = time.time()
        dt_list.append(t1-t0)
        cube = None

    print('HDD, NVMe, SSD, %d files written in cube' %n_frames)
    print(dt_list)

def read_test_cubes():
    dt_list = []
    for path in path_list:
        t0 = time.time()
        cube_path = os.path.join(path, cube_name)
        hduList = pyfits.open(cube_path)
        cube = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        t1 = time.time()
        dt_list.append(t1-t0)

    n_frames = cube.shape[2]
    print('HDD, NVMe, SSD, %d files in cube read' %n_frames)
    print(dt_list)

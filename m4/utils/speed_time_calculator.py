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

path0 = '/home/runa/testdir/'
path1 = '/mnt/cargo/data/testdir/'
path2 = '/mnt/magnum/data/testdir/'
#mypath = '/Users/rm/Desktop/Arcetri/M4/4D/files4D'
path_list = [path0, path1, path2]
#path_list = [mypath, mypath, mypath]

cube_name = 'Cube.fits'

def read_data_from_differents_paths(n_frames):
    dt_list = []
    for path in path_list:
        t0 = time.time()
        list = glob.glob(os.path.join(path,'*.4D'))
        for i in range(n_frames):
            ima = ic.fromNew4D(list[i])
        t1 = time.time()
        dt_list.append(round(t1-t0))

    print('SSD, NVMe, HHD, %d files read' %n_frames)
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
        dt_list.append(round(t1-t0))
        cube = None

    print('SSD, NVMe, HHD, %d files written in cube' %n_frames)
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
        dt_list.append(round(t1-t0))

    n_frames = cube.shape[2]
    print('SSD, NVMe, HHD, %d files in cube read' %n_frames)
    print(dt_list)

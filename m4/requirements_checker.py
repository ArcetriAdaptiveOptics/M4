'''
Authors
  - C. Selmi:  written in 2021
'''

import os
import numpy as np
from matplotlib import pyplot as plt
from m4.configuration import config_folder_names as config
from m4.analyzers import requirement_analyzer as req_check

'''
HOW TO USE IT::

    from m4.main import requirements_checker as rc
    rc.analysis_req(data_file_path, zernike_vector_to_subtract, step=None, offset=None)
'''

def analysis_req(data_file_path, zernike_vector_to_subtract, step=None, offset=None):
    ''' Analysis of noise requirements for a tn that create and use 3 different robust images.
    Function also crate a plot for results.

    Parameters
    ----------
    data_file_path: string
        total path for data analysis
    zernike_vector_to_subtract: numpy array
        vector of zernike to be subtracted from robust image

    Other Parameters
    ----------------
    step: int
        distance between patches
    offset: if it is None data analysis is made by split n_images in two
            else re-reads the offset image saved in the tt folder and subtracts it
            to each image during cube creation
    '''
    last_name = data_file_path.split('/')[-1]
    if last_name == 'hdf5':
        tt = data_file_path.split('/')[-2]
    else:
        tt = data_file_path.split('/')[-1]

    results_path = os.path.join(config.OUT_FOLDER, 'Req')
    dove = os.path.join(results_path, tt)
    if os.path.exists(dove):
        dove = dove
    else:
        os.makedirs(dove)
    fits_file_name = os.path.join(dove, 'info.txt')
    file = open(fits_file_name, 'w+')
    if offset is None:
        file.write('Data produced without offset optic image')
    else:
        file.write('Data produced with offset optic image')
    file.close()

    print('Creating cube 50')
    image50 = req_check.robustImageFromDataSet(50, data_file_path, zernike_vector_to_subtract, offset)
    print('Creating cube 100')
    image100 = req_check.robustImageFromDataSet(100, data_file_path, zernike_vector_to_subtract, offset)
    print('Creating cube 300')
    image300 = req_check.robustImageFromDataSet(300, data_file_path, zernike_vector_to_subtract, offset)
#     print('Creating cube 600')
#     image600 = req_check.robustImageFromDataSet(600, data_file_path, offset)

    image_list = [image50, image100, image300]  # , image600]
    slop_list, diff_piston_list, roc_list, rms31, rms500 = fromImagesToReq(image_list, None, step)

    x = np.array([50, 100, 300])  # ,600])
    # GRAFICO STD IMAGES
    y = np.array([image50.std(), image100.std(), image300.std()])  # ,image600.std()])
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'rms_image [m]', dove, 'std.png')
    # GRAFICO SLOPE
    y = np.array(slop_list)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'rms_slope [arcsec]', dove, 'slope.png')
    # GRAFICO DIFF PISTON
    y = np.array(diff_piston_list)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'diff_piston [m]', dove, 'diff_piston.png')
    # GRAFICO ROC
    y = np.array(roc_list)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'roc [m]', dove, 'roc.png')
    # GRAFICO RMS 31 MM
    y = np.array(rms31)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'rms_31mm [m]', dove, 'rms_31mm.png')
    # GRAFICO RMS 500 MM
    y = np.array(rms500)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'rms_500mm [m]', dove, 'rms_500mm.png')


def fromImagesToReq(image_list, pscale=None, step=None, n_patches=None):
    ''' Function that, given a list of images, calculates the outputs
    of all the required requirements

    Parameters
    ----------
    image_list: list
        list of image to be analyzed

    Other Parameters
    ----------------
    pscale: int
            value of image's pixel scale [px/m]
    step: int
        distance between patches
    offset: if it is None data analysis is made by split n_images in two
            else re-reads the offset image saved in the tt folder and subtracts it
            to each image during cube creation
    '''
    slop_list = []
    diff_piston_list = []
    roc_list = []
    rms31 = []
    rms500 = []
    print(pscale)
    for image in image_list:
        print('Producing slope')
        slop_list.append(req_check.test242(image, pscale))
        print('Producing differential piston')
        diff_piston_list.append(req_check.diffPiston(image))
        print('Producing roc')
        roc_list.append(req_check.test283(image, pscale, step))
        print('Producing rms31')
        rms31.append(req_check.test243(image, 0.015, pscale, step, n_patches))
        print('Producing rms51')
        rms500.append(req_check.test243(image, 0.1, pscale, step, n_patches))
    return slop_list, diff_piston_list, roc_list, rms31, rms500


def plotAndSaveForReqAnalysis(x, y, xlabel, ylabel, dove, image_name):
    ''' Function for plotting results
    '''
    plt.figure(figsize=(10, 6))
    plt.plot(np.sqrt(x), y, '-o')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    tt = dove.split('/')[-1]
    plt.title('%s' % tt)
    name = os.path.join(dove, image_name)
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)



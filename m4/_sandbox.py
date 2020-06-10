'''
@author: cs
'''

import os
import shutil
from astropy.io import fits as pyfits
import numpy as np
from m4.ground.configuration import Configuration
from m4.ground import object_from_fits_file_name
from m4.ground import object_from_fits_file_name as obj


### FUNZIONI PER TEST IFF ###
def testIFF_shuffleMeasureCreator(device, cmd_matrix_tag, mode_vect_tag,
                                  amp_tag, n_push_pull, template=None):
    ''' Al posto di n_push_pull passo in argomento il n_push_pull
    '''
    from m4.type.modesVector import ModesVector
    mv = ModesVector.loadFromFits(mode_vect_tag)
    mode_vect_input = mv.getModesVector()
    from m4.influence_functions_maker import IFFunctionsMaker
    IF = IFFunctionsMaker(device)

    tt = IF.acq_IFFunctions(mode_vect_tag, n_push_pull, \
                            amp_tag, cmd_matrix_tag, 1, template)

    folder = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
    who, tt_cmdH, acts_vector, cmd_matrix, \
    amplitude, n_push_pull, indexingList, template = IF.loadInfoFromFits(folder)

    cube = IF._testIFFunctions_createCube25fromFileFitsMeasure()
    from m4.type.commandHistory import CmdHistory
    cmdH = CmdHistory(device)
    ampl_reorg = cmdH._amplitudeReorganization(mode_vect_input, indexingList,
                                               amplitude, n_push_pull)

    misure = None
    for i in range(indexingList.shape[0]):
        for j in range(indexingList.shape[1]):
            mask = np.invert(cube[:,:,indexingList[i][j]].mask)
            k = i * indexingList.shape[1] + j
            if misure is None:
                misure = np.ma.masked_array(cube[:,:,indexingList[i][j]].data * template[0] *
                                            ampl_reorg[k], mask=mask)
                for l in range(1, template.shape[0]):
                    misure = np.ma.dstack((misure,
                                           np.ma.masked_array(cube[:,:,indexingList[i][j]].data *
                                                              ampl_reorg[k] * template[l],
                                                              mask=mask)))
            else:
                for l in range(template.shape[0]):
                    misure = np.ma.dstack((misure,
                                           np.ma.masked_array(cube[:,:,indexingList[i][j]].data *
                                                              ampl_reorg[k] * template[l],
                                                              mask=mask)))

    fits_file_name = os.path.join(folder, 'misure.fits')
    header = pyfits.Header()
    header['NPUSHPUL'] = n_push_pull
    header['WHO'] = who
    header['TT_CMDH'] = tt_cmdH
    pyfits.writeto(fits_file_name, acts_vector, header)
    pyfits.append(fits_file_name, cmd_matrix, header)
    pyfits.append(fits_file_name, amplitude, header)
    pyfits.append(fits_file_name, indexingList, header)
    pyfits.append(fits_file_name, misure.data, header)
    pyfits.append(fits_file_name, misure.mask.astype(int), header)
    pyfits.append(fits_file_name, template, header)
    return tt

def testIFF_tidyMeasureCreator(device, cmd_matrix_tag, mode_vect_tag,
                               amp_tag, n_push_pull, template=None):
    from m4.influence_functions_maker import IFFunctionsMaker
    IF = IFFunctionsMaker(device)

    tt = IF.acq_IFFunctions(mode_vect_tag, n_push_pull, amp_tag, cmd_matrix_tag,
                            None, template)

    folder = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
    who, tt_cmdH, acts_vector, cmd_matrix, \
            amplitude, n_push_pull, indexingList, template = IF.loadInfoFromFits(folder)

    cube = IF._testIFFunctions_createCube25fromFileFitsMeasure()
    ampl = np.tile(amplitude, n_push_pull)

    misure = None
    for i in range(indexingList.shape[0]):
        for j in range(indexingList.shape[1]):
            mask = np.invert(cube[:,:,indexingList[i][j]].mask)
            k = i * indexingList.shape[1] + j
            if misure is None:
                misure = np.ma.masked_array(cube[:,:,indexingList[i][j]].data * template[0] *
                                            ampl[k], mask=mask)
                for l in range(1, template.shape[0]):
                    misure = np.ma.dstack((misure,
                                       np.ma.masked_array(cube[:,:,indexingList[i][j]].data *
                                                              ampl[k] * template[l],
                                                              mask=mask)))
            else:
                for l in range(template.shape[0]):
                    misure = np.ma.dstack((misure,
                                           np.ma.masked_array(cube[:,:,indexingList[i][j]].data *
                                                              ampl[k] * template[l],
                                                              mask=mask)))

    fits_file_name = os.path.join(folder, 'misure.fits')
    header = pyfits.Header()
    header['NPUSHPUL'] = n_push_pull
    header['WHO'] = who
    header['TT_CMDH'] = tt_cmdH
    pyfits.writeto(fits_file_name, acts_vector, header)
    pyfits.append(fits_file_name, cmd_matrix, header)
    pyfits.append(fits_file_name, amplitude, header)
    pyfits.append(fits_file_name, indexingList, header)
    pyfits.append(fits_file_name, misure.data, header)
    pyfits.append(fits_file_name, misure.mask.astype(int), header)
    pyfits.append(fits_file_name, template, header)
    return tt


def testIFF_an(tt, ttD=None):
    from m4.analyzer_iffunctions import AnalyzerIFF
    file_name = os.path.join \
                (Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
    an = AnalyzerIFF.loadTestMeasureFromFits(file_name)
    cube = an.createCube(ttD)
    an.saveCubeAsFits('Cube.fits')
    int_mat = an.getInteractionMatrix()
    rec = an.getReconstructor()
    prod = np.dot(rec, int_mat)
    return an, prod, cube

def testIFF_spiano(an):
    ampr = np.random.randn(25)
    wf = np.dot(an._cube, ampr)

    an.setDetectorMask(wf.mask | an.getMasterMask())
    rec = an.getReconstructor()
    wf_masked = np.ma.masked_array(wf.data,
                                   mask=np.ma.mask_or(wf.mask, an.getMasterMask()))

    amp = np.dot(rec, wf_masked.compressed())
    surf = np.dot(an._cube, amp)
    sp_wf = wf - surf
    return amp, sp_wf

### FINE FUNZIONI TEST IFF ###


### Code test for data: /mnt/shared/M4Data/OPTData/OPDImages/20191210_143625/hdf5 ###

def testOPDImages_TestFormat():
    ''' In questo modo sto portando tutti i dati che sono in OPDImages nella forma
        che servono a me per far girare la funzione di creazione cubo fatta per il test
    '''
    from m4.analyzer_iffunctions import AnalyzerIFF
    an = AnalyzerIFF()
    # metterci an.createTestCubeMeasure_forTestDataInOpdImages ma dal pc giù
    cubeMeasure = an.createTestCubeMeasure_forTestDataInOpdImages() #meglio se lo leggo che ci mette una vita
    cube = an.createCube()
    return cube

def testOPDImages_madeToOrder():
    from m4.analyzer_iffunctions import AnalyzerIFF
    an = AnalyzerIFF()
    cube = an.createCube_forTestDataInOpdImages()
    return cube

### Fine test ###

def immaginiprova():
    fits_root = "/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova"
    fits_file_name = os.path.join(fits_root, 'mode_0005.fits')
    hduList = pyfits.open(fits_file_name)
    ima = hduList[0].data
    m4 = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
    fits_file_name = os.path.join(fits_root, 'mode_0006.fits')
    hduList = pyfits.open(fits_file_name)
    ima = hduList[0].data
    segment = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
    return m4, segment

def immaginiProvaTTDetrendSeg():
    push = object_from_fits_file_name.readImageFromFitsFileName \
            ('Seg/img_0000.fits')
    pull = object_from_fits_file_name.readImageFromFitsFileName \
            ('Seg/img_0001.fits')
    mode_0 = np.ma.masked_array(pull.data - push.data, mask=push.mask)

#push= object_from_fits_file_name.readImageFromFitsFileName('Seg/img_0002.fits')
#pull= object_from_fits_file_name.readImageFromFitsFileName('Seg/img_0003.fits')
#mode1= np.ma.masked_array(pull.data - push.data, mask= np.invert(push.mask))
    return mode_0

def immaginiProvaTTDetrendAll():
    push = object_from_fits_file_name.readImageFromFitsFileName \
            ('All/img_0000.fits')
    pull = object_from_fits_file_name.readImageFromFitsFileName \
            ('All/img_0001.fits')
    mode_0 = np.ma.masked_array(pull.data - push.data, mask=push.mask)

    push = object_from_fits_file_name.readImageFromFitsFileName \
            ('All/img_0002.fits')
    pull = object_from_fits_file_name.readImageFromFitsFileName \
            ('All/img_0003.fits')
    mode_1 = np.ma.masked_array(pull.data - push.data, mask=push.mask)
    return mode_0, mode_1

def loadTestCubeFromRunaIFF(tt, n_measure):
    import m4.ground.object_from_fits_file_name as obj
    cube = None

    for i in range(n_measure):
        name = 'mode_%04d.fits' %i
        file_name = os.path.join(tt, name)
        ima = obj.readImageFromRunaIFFs(file_name)
        if cube is None:
            cube = ima
        else:
            cube = np.ma.dstack((cube, ima))
    return cube


def immaginiDaIFFRuna():
    doveSeg = '20170630_105105/mode_0197.fits'
    #doveM4 = '20161226_122557/mode_0267.fits'

    seg = object_from_fits_file_name.readImageFromRunaIFFs(doveSeg)
    #m4 = object_from_fits_file_name.readImageFromRunaIFFs(doveM4)
    return seg

# Per fare la tt list sto usando sempre tt = '20170630_105105'
def provaZernike(seg, zernike_modes_vector_amplitude, tt_list_for_an):
    from m4.utils.roi import ROI
    r = ROI()
    roi = r.roiGenerator(seg)
    from m4.zernike_command_test import ZernikeCommand
    zc = ZernikeCommand(roi[3], tt_list_for_an)
    tt, surf_cube, m4_images_cube = zc.zernikeCommandTest(zernike_modes_vector_amplitude)
    return zc, tt, surf_cube, m4_images_cube

# TEST noise #
def provaAcquisitionNoise(template):
    cmd_matrix_tag = 'matTestIF.fits'
    mode_vect_tag = 'vectTestIF.fits'
    #amp_tag = 'ZerosAmp.fits'
    amp_tag = 'OnesAmp.fits'
    n_push_pull = 3
    #template = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1])
    from m4.utils import create_device
    device = create_device.myDevice("segment")

    from m4.influence_functions_maker import IFFunctionsMaker
    IF = IFFunctionsMaker(device)

    tt = IF.acq_IFFunctions(mode_vect_tag, n_push_pull, amp_tag, cmd_matrix_tag,
                            None, template)
    source = os.path.join('/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions', tt)
    destination = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Noise'
    destination_file_path = shutil.move(source, destination)

    from m4.noise_functions import Noise
    n = Noise()
    destination_file_path = os.path.join(destination, tt)
    data_file_path = os.path.join(destination, 'hdf5')
    cube_measure = n._createAndSaveCubeFromH5Data(data_file_path, destination_file_path, device)

    return destination_file_path

### prova ROI ###
def provaAutRoi():
    file_name = 'ProvaCodice/Immagini_prova/mode_0005.fits'
    central_view = obj.readImageFromFitsFileName(file_name)
    file_name = 'ProvaCodice/Immagini_prova/mode_0006.fits'
    segment_view = obj.readImageFromFitsFileName(file_name)
    file_name = 'ProvaCodice/Immagini_prova/Allineamento/20191007_134908/img.fits'
    rm_in = obj.readImageFromFitsFileName(file_name) #for seg
    return central_view, segment_view, rm_in

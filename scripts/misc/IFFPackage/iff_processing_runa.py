"""
Authors
    - P. Ferriauolo
    - R. Briguglio

    Written in june 2024

------------------------

Module which contains functions for data processing and analysis of M4's IFF.
"""
import os
import glob
import numpy as np
from astropy.io import fits as pyfits
from m4.ground import timestamp
from m4.ground import zernike as zern
from m4.configuration import read_iffconfig
from m4.mini_OTT import timehistory as th
from m4.ground import read_data as rd
from m4.configuration import config_folder_names as fn
#from m4.devices import deformable_mirror as dm
#from scripts.misc.IFFPackage import iff_acquisition_preparation

#ifa = iff_acquisition_preparation.IFFCapturePreparation(dm.M4AU())
regEnd = 0
trigFrame = 0
imgFold     = fn.OPD_IMAGES_ROOT_FOLDER
ifFold      = fn.IFFUNCTIONS_ROOT_FOLDER
intMatFold  = fn.INTMAT_ROOT_FOLDER
confFold    = fn.CONFIGURATION_ROOT_FOLDER

ampVecFile = 'ampVector.fits'
modesVecFile = 'modesVector.fits'
templateFile = 'Template.fits'
regisActFile = 'registrationActs.fits'
shuffleFile  = 'shuffle.dat'
indexListFile= 'indexList.fits'

frameCenter = [200,200]

def process(tn, register = False):
    """
    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    infoT, infoR, infoIF = getAcqInfo(tn)
    ampVector, modesVector, template,indexList, registrationActs, shuffle = getAcqPar(tn)

    regMat, modesMat = findTriggerStartFrame(tn)
    actImgList = registrationRedux(regMat, template)
    modesMatReorg = modesReorganization(modesMat)
  
    iffRedux(modesMatReorg)
    if register is not False:
        dx = findFrameOffset(actImgList, actlist)
    else:
        dx = register
    createCube(tn, dx)

def iffRedux(fileMat, ampVect, modeList, template, shuffle=0):
    """
    Parameters
    ----------
    fileMat : string, matrix
        is already reorganized for the shuffle, if requested

    Returns
    -------
    None.

    """
    fold = os.path.join(ifFold, tn)
    nmodes = len(ampVect)
    nframes= len(template)
    if shuffle != 0:
        nframes = shuffle*3
    for i in range(0, nmodes-1):
        img = pushPullRedux(fileMat[i,:], template, shuffle)
        img = img / (2*ampVect[i])
        img_name = os.path.join(fold, 'mode_{:04d}.fits'.format(modeList[i]))
        rd.save_phasemap(img_name, norm_image)

def registrationRedux(fileMat, template):
    """ 
    Parameters
    ----------
    regFrames : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    nAct = (fileMat.shape)[0]
    imglist = []
    for i in range(0, nAct-1):
        img = pushPullRedux(fileMat[i,:], template)
        imglist.append(img)
    return imglist

def findFrameOffset(imglist, actlist):
    """
    This function computes the position difference between the current frame and a reference
    Parameters
    ----------
    imglist : list | masked arrays
        list of the actuator images to be used
    actlist: int | array
        list of actuators (index)

    Returns
    -------
    dp: float
        position difference
    """
    actcoordFile =  os.path.join(ifFold, tn,coordfile)
    actcoord = rd.readFits_data(actcoordfile)
    xy = fa.findFrameCoord(imglist, actlist, actcoord)
    dp = xy - frameCenter
    return dp


def pushPullRedux(fileVec, template, shuffle=0):
    """
    This function performs the basic operation of processing PushPull data. It has been extracted from the higher level scripts to allow re-using at a lower level where requested (e.g. processing of PushPull for registration frames).
   Parameters
    ----------
    fileVec : string | array
        it is a row in the fileMat (the organized matrix of the images filename), corresponding to all the realizations of the same mode (or act), with given template. If shuffle option has been used, the fileMat (and fileVec) shall be reorganized before running this script
    template: int | array
        template for the PushPull acquisition
    shuffle: int 
        keyword for shuffle option used or not. if 0, shuffle not used. if !=0 it is the repetitions of the templated samplin

    Returns
    image: masked array
        final processed image
    -------
    None.

    """
    image_list = []
    for l in range(0, template.shape[0]-1):
        ima = rd.read_phasemap(fileVec[l])
        image_list.append(ima)
    image = np.zeros((ima.shape[0], ima.shape[1]))
    if shuffle == 0:
        for x in range(1, len(image_list)):
            opd2add = image_list[x]*template[x] + image_list[x-1]*template[x-1]
            master_mask2add = np.ma.mask_or(image_list[x].mask, image_list[x-1].mask)
            if x==1:
                master_mask = master_mask2add
            else:
                master_mask = np.na.mask_or(master_mask, master_mask2add)
            image += opd2add
    else:
        print('Shuffle option')
        for i in range(0, shuffle-1):
            for x in range(1,2):
                opd2add = image_list[i*3+x]*template[x] + image_list[i*3+x-1]*template[x-1]
                master_mask2add = np.ma.mask_or(image_list[i*3+x].mask, image_list[i*3+x-1].mask)
                if (i == 0 and x == 1):
                    master_mask = master_mask2add
                else:
                    master_mask = np.na.mask_or(master_mask, master_mask2add)
                image += opd2add
    
    image = np.ma.masked_array(image, mask=master_mask)/ (template.shape[0]-1)
    return image

def createCube(tn, register):
    """
    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    cube : TYPE
        DESCRIPTION.

    """
    filelist = th.fileList(tn)
    cube_list = []
    for imgfits in filelist:
        image = rd.readFits_maskedImage(imgfits)
        if register is not False:
             image= np.roll(image, register)
        #qui aggiungere la registrazione?
        cube_list.append(image)
    cube = np.ma.dstack(cube_list)
    path = os.path.join(fn.INTMAT_ROOT_FOLDER, tn)
    save_name = os.path.join(path, 'IMCube.fits')
    rd.save_phasemap(save_name, cube, isCube=True)
    print("Cube saved in '{}'".format(path))
    return cube #!!! return?

def stackCubes(tnlist):
    """
    Parameters
    ----------
    tnlist : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    search_fold = fn.IFFUNCTIONS_ROOT_FOLDER
    cubelist = []
    for tn in tnlist:
        fold = os.path.join(search_fold, tn)
        cube_name = os.path.join(fold, 'IMCube.fits')
        hdul = pyfits.open(cube_name)
        cubelist.append(hdul[0].data)
        hdul.close()
        #!!! Add the matrix stacking
    stacked_cube = np.ma.dstack(cubelist)
    new_tn = timestamp.Timestamp.now()
    save_fold = os.path.join(fn.INTMAT_ROOT_FOLDER, new_tn)
    rd.save_phasemap(os.path.join(save_fold, 'IMCube.fits'), stacked_cube, isCube=True)
    print("Stacked cube saved in {}".format(save_fold))
    return stacked_cube #!!! Return?
    

def findTriggerStartFrame(tn, amplitude=None):
    """
    This function identifies the triggering frame and the registration frames in
    the data frame sequence.

    Parameters
    ----------
    tn : str / list
        Tracking number of the folder containing the images.
    amplitude: float
        Amplitude of the trigger command applied. If no value is passed, it will
        be loaded from the 'iffconfig.ini' configuration file.

    Returns
    -------
    trigFrame : int
        Index of the frame containing the trigger. The subsequent frame (triggerId + 1) 
        is the first 'useful' frame in the time history.
    regFrames : int | ArrayLike
        List containing the first and last indices of the registration frames. 
        The subsequent frame is the start of the IFF acquisition commanded matrix
        history.
    """
    # !!! Readable way
    # flist = os.listdir(os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn))
    # flist.sort()
    # filelist = [os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, (tn+'/'+image)) for image in flist]
    # Concise way
    ampVector, modesVector, template,indexList, registrationActs, shuffle = getAcqPar(tn)
    nmodes = len(modesVector)
    nPP = len(template)
    infoT, infoR, infoIF = getAcqInfo(tn)
    nTzeros = infoT[0]
    triggAmp = infoT[2]
    nRzeros = infoR[0]
    nRact = len(registrationActs)  #len(infoR[1])
    nPPR = len(infoR[3])
    nIfzeros = infoIF[0]
    fileList = getFileList(tn)
    img0 = rd.read_phasemap(fileList[0])
    go = 1
    i = 1
    while go !=0:
        if go > nTzeros:
            print('Heading Zeros exceeded, error!!')
        img1 = rd.read_phasemap(fileList[i])
        rr2check = zern.removeZernike(img1-img0,[1,2,3]).std()
        if rr2check > triggAmp/3:
            go=0
        else:
            i = i+1
            img0 = img1
    
    print('Trigger frame found at position:' + str(i))
    firstReg = i+1+nRzeros
    regFrameNames = fileList[firstReg:firstReg+nRact*nPPR]
    regMat = np.reshape(regFrameNames,(nRact,nPPR))
    #ifFuncFirst = i+1+nRact*nPPR+1+nIfzeros  #to be checked!!!
    firstIFF = firstReg+nRact*nPPR+nIfzeros
    ifFramenames = fileList[firstIFF:firstIFF+nmodes*nPP+1]
    if len(ifFramenames) != len(fileList[firstIFF:]):
        print('Warning! Too much files in folder')
    ifMat = np.reshape(ifFramenames,(nmodes,nPP))
    print('Registration file matrix:')
    print(regMat.shape)
    print('IFFunc file matrix:')
    print(ifMat.shape)

    return regMat, ifMat #imgMode0 = ifMat[0,:]

def getAcqInfo(tn):
    '''
    '''
    infoT = read_iffconfig.getConfig('TRIGGER')
    infoR = read_iffconfig.getConfig('REGISTRATION')
    infoIF = read_iffconfig.getConfig('IFFUNC')
    return infoT, infoR, infoIF

def getAcqPar(tn):
    '''
    Utility function to read all the sampling parameters, necessary for the processing
    Parameters
    ----------
    tn : string
        Tracking number of the dataset

    Returns
    -------
    ampVector: float, array
        vector of the command amplitudes
    modesVector: int, array
        list of modes sampled
    template: int, array
        acquisition template, typically a sequence of +1;-1 (es: [+1,-1,+1])
    indexList: int, array
        positions index of the modes sampling, in case of shuffle option
    registrationActs: int, array
        list of actuators sampled for the frames registration
    shuffle
    TYPE
        DESCRIPTION.

    """

    '''
    base = os.path.join(ifFold, tn)
    ampVector = rd.readFits_data(os.path.join(base, ampVecFile))
    template = rd.readFits_data(os.path.join(base,templateFile))
    modesVector = rd.readFits_data(os.path.join(base, modesVecFile))
    indexList = rd.readFits_data(os.path.join(base,indexListFile))
    registrationActs = rd.readFits_data(os.path.join(base, regisActFile))
    shf = open(os.path.join(base, shuffleFile), 'r')
    shuffle =int(shf.read())
    shf.close()
    return ampVector, modesVector, template,indexList, registrationActs, shuffle

def getFileList(tn):
    bb = os.path.join(imgFold, tn)
    name  = '/*.4Ds'
    print('Replace search key')
    lsdirs = glob.glob(bb + name)
    lsdirs.sort(key=_sortFunc4D)
    lsdirs.sort(key=_sortFunc4D)
    return lsdirs


def _sortFunc4D(elem):
    iid = os.path.basename(elem)[:-4]
    print('Replace search key')
    iis = "%5.5i" % int(iid)
    return iis


def _ampReorganization():
    """
    

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    indexList = ifa.getIndexingList()
    return indexList

def _indexReorganization():
    """
    

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    amps = ifa.getAmplitude()
    return amps

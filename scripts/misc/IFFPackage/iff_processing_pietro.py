"""
Authors
    - P. Ferriauolo
    - R. Briguglio

    Written in june 2024

------------------------

Module which contains functions for data processing and analysis of M4's IFF.
"""
import os
import numpy as np
from m4.ground import timestamp
from m4.configuration import read_iffconfig
from m4.ground import read_data as rd
from m4.configuration import config_folder_names as fn
from m4.ground import zernike as zern

imgFold     = fn.OPD_IMAGES_ROOT_FOLDER
ifFold      = fn.IFFUNCTIONS_ROOT_FOLDER
intMatFold  = fn.INTMAT_ROOT_FOLDER
confFold    = fn.CONFIGURATION_ROOT_FOLDER
frameCenter = [200,200]
ampVecFile     = 'ampVector.fits'
modesVecFile   = 'modesVector.fits'
templateFile   = 'Template.fits'
regisActFile   = 'registrationActs.fits'
shuffleFile    = 'shuffle.dat'
indexListFile  = 'indexList.fits'
coordfile      = '' #!!! Da definire?

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
    infoT, infoR, infoIF = _getAcqInfo()
    ampVector, modesVector, template, indexList, registrationActs, shuffle = _getAcqPar(tn)
    _, regMat = getRegFileMatrix(tn)
    modesMat = getIffFileMatrix(tn)
    # regMat, modesMat = findTriggerStartFrame(tn)
    actImgList = registrationRedux(regMat, template)
    modesMatReorg = _modesReorganization(modesMat) # ???
    iffRedux(tn, modesMatReorg, ampVector, modesVector, template, shuffle)
    if register is not False:
        dx = findFrameOffset(tn, actImgList, registrationActs)
    else:
        dx = register
    saveCube(tn)
    
def saveCube(tn):
    """
    Creates and save a cube from the fits files contained in the tn folder, 
    along with the command matrix and the modes vector fits.

    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    cube : TYPE
        DESCRIPTION.

    """
    filelist = _getFileList(tn)
    cube_list = []
    for imgfits in filelist:
        image = rd.readFits_maskedImage(imgfits)
        cube_list.append(image)
    cube = np.ma.dstack(cube_list)
    # Saving the cube
    cube_path = os.path.join(intMatFold, tn, 'IMCube.fits')
    rd.save_phasemap(cube_path, cube, isCube=True)
    # Copying the cmdMatrix and the ModesVector into the INTMAT Folder
    cmat = rd.readFits_data(os.path.join(ifFold, tn, 'cmdMatrix.fits'))
    mvec = rd.readFits_data(os.path.join(ifFold, tn, 'modesvector.fits'))
    rd.save_phasemap(os.path.join(intMatFold, tn, 'cmdMatrix.fits'), cmat)
    rd.save_phasemap(os.path.join(intMatFold, tn, 'modesVector.fits'), mvec)
    print("Cube saved in '{}'".format(cube_path))

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
    new_tn = timestamp.Timestamp.now()
    cubelist = []
    matrixList = []
    modesVectList = []
    for tn in tnlist:
        fold = os.path.join(ifFold, tn)
        cube_name = os.path.join(fold, 'IMCube.fits')
        matrix_name = os.path.join(fold, 'cmdMatrix.fits')
        modesVec_name = os.path.join(fold, 'modesVector.fits')
        cubelist.append(rd.readFits_data(cube_name))
        matrixList.append(rd.readFits_data(matrix_name))
        modesVectList.append(rd.readFits_data(modesVec_name))
    stacked_cube = np.ma.dstack(cubelist)
    stacked_cmat = np.dstack(matrixList)
    stacked_mvec = np.dstack(modesVectList)
    save_cube = os.path.join(fn.INTMAT_ROOT_FOLDER, new_tn, 'IMCube.fits')
    save_cmat = os.path.join(fn.INTMAT_ROOT_FOLDER, new_tn, 'cmdMatrix.fits')
    save_mvec = os.path.join(fn.INTMAT_ROOT_FOLDER, new_tn, 'modesVector.fits')
    rd.save_phasemap(save_cube, stacked_cube, isCube=True)
    rd.save_phasemap(save_cmat, stacked_cmat)
    rd.save_phasemap(save_mvec, stacked_mvec)
    print("Stacked cube and matrices saved in {}".format(new_tn))

# def iffRedux(tn):
#     """
#     Processes 4D images to produce FITS images
#     This function retrieves the 4D images acquired and pre-processed from the
#     interferometer, applies the differential algorythm, and produces a fits image
#     for each measured mode.
#     Parameters
#     ----------
#     tn : str 
#         Tracking number in the OPDImages folder where the acquired data is stored
#     Returns
#     -------
#     """
#     _,modesList,_,template,_  = read_iffconfig.getConfig('IFFUNC')
#     nPushPull = len(template)
#     #indexingList = _indexReorganization()    # to be implemented
#     #amplitude = _ampReorganization()      # to be implemented
#     indexingList = rd.readFits_data(os.path.join(fn.IFFUNCTIONS_ROOT_FOLDER, tn)+'/indexList.fits')
#     amplitude = rd.readFits_data(os.path.join(fn.IFFUNCTIONS_ROOT_FOLDER, tn)+'/ampVector.fits')
#     filelist = sorted([os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, (tn+'/'+image)) for image in os.listdir(os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn))])
#     filelist = filelist[regEnd:]
#     for i in range(len(modesList)):
#         print('Mode', i) #!!! Debugging, no need to print
#         for k in range(nPushPull):
#             p = nPushPull * i + k
#             n = indexingList[i]
#             # Crea il pacchetto di immagini del modo 'i', contenente nPushPull images
#             file_name = filelist[p]
#             image_list = []
#             for l in range(0, template.shape[0]-1):
#                 file_name = filelist[p+l]
#                 ima = rd.read_phasemap(file_name)
#                 image_list.append(ima)
#             image = np.zeros((ima.shape[0], ima.shape[1]))
#             # Algorimo differenziale
#             for x in range(1, len(image_list)):
#                 opd2add = image_list[x]*template[x] + image_list[x-1]*template[x-1]
#                 master_mask2add = np.ma.mask_or(image_list[x].mask, image_list[x-1].mask)
#                 if x==1:
#                     master_mask = master_mask2add
#                 else:
#                     master_mask = np.na.mask_or(master_mask, master_mask2add)
#                 image += opd2add
#             image = np.ma.masked_array(image, mask=master_mask)
#             norm_image = image / (2*amplitude[n] * (template.shape[0]-1))
#         fold = os.path.join(ifFold, tn)
#         img_name = os.path.join(fold, 'mode_{:04d}.fits'.format(i))
#         rd.save_phasemap(img_name, norm_image)

def iffRedux(tn, fileMat, ampVect, modeList, template, shuffle=0):
    """
    

    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.
    fileMat : TYPE
        DESCRIPTION.
    ampVect : TYPE
        DESCRIPTION.
    modeList : TYPE
        DESCRIPTION.
    template : TYPE
        DESCRIPTION.
    shuffle : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    None.

    """
    fold = os.path.join(ifFold, tn)
    nmodes = len(ampVect)
    # nframes= len(template)
    # if shuffle != 0:
    #     nframes = shuffle*3
    for i in range(0, nmodes):
        img = pushPullRedux(fileMat[i,:], template, shuffle)
        norm_img = img / (2*ampVect[i])
        img_name = os.path.join(fold, 'mode_{:04d}.fits'.format(modeList[i]))
        rd.save_phasemap(img_name, norm_img)

def pushPullRedux(fileVec, template, shuffle=0):
    """
    This function performs the basic operation of processing PushPull data. 
    It has been extracted from the higher level scripts to allow re-using at a 
    lower level where requested (e.g. processing of PushPull for registration frames).
    Parameters
    ----------
    fileVec : string | array
        it is a row in the fileMat (the organized matrix of the images filename), 
        corresponding to all the realizations of the same mode (or act), with 
        given template. If shuffle option has been used, the fileMat (and fileVec) 
        shall be reorganized before running this script
    template: int | array
        template for the PushPull acquisition
    shuffle: int 
        keyword for shuffle option used or not. if 0, shuffle not used. if !=0 
        it is the number of repetitions of the templated sampling

    Returns
    image: masked array
        final processed image
    -------
    None.

    """
    image_list = []
    template = np.array(template)
    for l in range(0, template.shape[0]):
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
    nAct = fileMat.shape[0]
    imglist = []
    for i in range(0, nAct-1):
        img = pushPullRedux(fileMat[i,:], template)
        imglist.append(img)
    return imglist

def findFrameOffset(tn, imglist, actlist):
    """
    This function computes the position difference between the current frame and
    a reference
    
    Parameters
    ----------
    tn : str
        Tracking number
    imglist : list | masked arrays
        List of the actuator images to be used
    actlist: int | array
        List of actuators (index)

    Returns
    -------
    dp: float
        Position difference
    """
    actCoordFile = os.path.join(ifFold, tn, coordfile)
    actCoord = rd.readFits_data(actCoordFile)
    xy = fa.findFrameCoord(imglist, actlist, actCoord) # ???
    dp = xy - frameCenter
    return dp

# !!! Deprecated - divided into the three functions below
# def getTriggerAndRegFrame(tn, amplitude=None):
#     """
#     This function identifies the triggering frame and the registration frames in
#     the data frame sequence.

#     Parameters
#     ----------
#     tn : str / list
#         Tracking number of the folder containing the images.
#     amplitude: float, optional
#         Amplitude of the trigger command applied. If no value is passed, it will
#         be loaded from the 'iffconfig.ini' configuration file.

#     Returns
#     -------
#     trigFrame : int
#         Index of the frame containing the trigger. The subsequent frame (triggerId + 1) 
#         is the first 'useful' frame in the time history.
#     regFrames : int | ArrayLike
#         List containing the first and last indices of the registration frames. 
#         The subsequent frame is the start of the IFF acquisition commanded matrix
#         history.
#     """
#     infoT = read_iffconfig.getConfig('TRIGGER')
#     infoR = read_iffconfig.getConfig('REGISTRATION')
#     infoIF= read_iffconfig.getConfig('IFFUNC')
#     timing = read_iffconfig.getTiming()
#     if amplitude is not None:
#         infoT['amplitude'] = amplitude
#     fileList = _getFileList(tn)
#     img0 = rd.read_phasemap(fileList[0])
#     go = 0
#     i = 1
#     while go !=0:
#         if go > infoT['zeros']:
#             raise RuntimeError('Heading Zeros exceeded, error!!')
#         img1 = rd.read_phasemap(fileList[i])
#         rr2check = zern.removeZernike(img1-img0,[1,2,3]).std()
#         if rr2check > infoT['amplitude']/3:
#             go=0
#         else:
#             i+=1
#             img0 = img1
#     trigFrame = i
#     regStart  = trigFrame + infoR['zeros']*timing +1
#     regEnd    = regStart + len(infoR['modes'])*len(infoR['template'])*timing
#     regList   = fileList[regStart:regEnd]
#     regMat    = np.reshape(regList, (len(infoR['modes']), len(infoR['template'])))
#     iffList   = fileList[regEnd+infoIF['zeros']:]
#     iffMat    = np.reshape(iffList, (len(infoIF['modes']), len(infoIF['template'])))
#     return regMat, iffMat

def getTriggerFrame(tn, amplitude=None):
    infoT, _, _ = _getAcqInfo()
    if amplitude is not None:
        infoT['amplitude'] = amplitude
    fileList = _getFileList(tn)
    img0 = rd.read_phasemap(fileList[0])
    go = 0
    i = 1
    while go !=0:
        if go > infoT['zeros']:
            raise RuntimeError('Heading Zeros exceeded, error!!')
        img1 = rd.read_phasemap(fileList[i])
        rr2check = zern.removeZernike(img1-img0,[1,2,3]).std()
        if rr2check > infoT['amplitude']/3:
            go=0
        else:
            i+=1
            img0 = img1
    trigFrame = i
    return trigFrame
    
def getRegFileMatrix(tn):
    fileList    = _getFileList(tn)
    _, infoR, _ = _getAcqInfo()
    timing      = read_iffconfig.getTiming()
    trigFrame   = getTriggerFrame(tn)
    regStart    = trigFrame + infoR['zeros']*timing +1
    regEnd      = regStart + len(infoR['modes'])*len(infoR['template'])*timing
    regList     = fileList[regStart:regEnd]
    regMat      = np.reshape(regList, (len(infoR['modes']), len(infoR['template'])))
    return regEnd, regMat
    
def getIffFileMatrix(tn):
    fileList    = _getFileList(tn)
    infoIF      = _getAcqInfo()
    regEnd, _   = getRegFileMatrix(tn)
    iffList     = fileList[regEnd+infoIF['zeros']:]
    iffMat      = np.reshape(iffList, (len(infoIF['modes']), len(infoIF['template'])))
    return iffMat

def _getFileList(tn, fold=None):
    """
    Returns the file list of a given tn datapath.

    Parameters
    ----------
    tn : str
        DESCRIPTION.
    fold : str, optional
        DESCRIPTION. If None, the folder is OPD_IMAGES_ROOT_FOLDER 

    Returns
    -------
    fl : list
        List of files inside folder.

    """
    if fold is None:
        fl = sorted([os.path.join(imgFold, (tn+'/'+image)) for image in os.listdir(os.path.join(imgFold, tn))])
    else:
        fl = sorted([os.path.join(fold, (tn+'/'+image)) for image in os.listdir(os.path.join(fold, tn))])
    return fl

def _getAcqPar(tn):
    """
    

    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    ampVector : TYPE
        DESCRIPTION.
    modesVector : TYPE
        DESCRIPTION.
    template : TYPE
        DESCRIPTION.
    indexList : TYPE
        DESCRIPTION.
    registrationActs : TYPE
        DESCRIPTION.
    shuffle : TYPE
        DESCRIPTION.

    """
    base = os.path.join(ifFold, tn)
    ampVector        = rd.readFits_data(os.path.join(base, ampVecFile))
    template         = rd.readFits_data(os.path.join(base, templateFile))
    modesVector      = rd.readFits_data(os.path.join(base, modesVecFile))
    indexList        = rd.readFits_data(os.path.join(base, indexListFile))
    registrationActs = rd.readFits_data(os.path.join(base, regisActFile))
    shf = open(os.path.join(base, shuffleFile), 'r')
    shuffle = int(shf.read())
    shf.close()
    return ampVector, modesVector, template,indexList, registrationActs, shuffle

def _getAcqInfo():
    """
    

    Returns
    -------
    infoT : TYPE
        DESCRIPTION.
    infoR : TYPE
        DESCRIPTION.
    infoIF : TYPE
        DESCRIPTION.

    """
    infoT = read_iffconfig.getConfig('TRIGGER')
    infoR = read_iffconfig.getConfig('REGISTRATION')
    infoIF = read_iffconfig.getConfig('IFFUNC')
    return infoT, infoR, infoIF

def _ampReorganization(ampVector):
    pass

def _modesReorganization(modesVector):
    pass


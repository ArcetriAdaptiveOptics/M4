"""
Authors
    - P. Ferraiuolo
    - R. Briguglio
        Written in june 2024
------------------------
Module containing functions to process the data acquired for IFF measurements.
"""
import os
import numpy as np
from m4.ground import timestamp
from m4.configuration import read_iffconfig
from m4.ground import read_data as rd
from m4.configuration import config_folder_names as fn
from m4.ground import zernike as zern
from astropy.io import fits as pyfits
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
    High level function with processes the data contained in the given tracking
    number OPDimages folder, performing the differential algorythm and saving 
    the final cube.

    Parameters
    ----------
    tn : str
        Tracking number of the data in the OPDImages folder.
    register : int, optional
        Parameter which enables the registration option. The default is False.
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
    saveCube(tn, dx)
    
def saveCube(tn, register=False):
    """
    Creates and save a cube from the fits files contained in the tn folder, 
    along with the command matrix and the modes vector fits.

    Parameters
    ----------
    tn : str
        Tracking number of the IFFunctions data folder from which create the cube.
        
    Returns
    -------
    cube : masked_array
        Data cube of the images, with shape (npx, npx, nmodes).
    """
    filelist = _getFileList(tn, fold=ifFold)
    filelist = [file for file in filelist if 'mode_' in file]
    cube_list = []
    for imgfits in filelist:
        image = rd.readFits_maskedImage(imgfits)
        if register is not False:
            image= np.roll(image, register)
        cube_list.append(image)
    cube = np.ma.dstack(cube_list)
    # Saving the cube
    new_fold = os.path.join(intMatFold, tn)
    os.mkdir(new_fold)
    cube_path = os.path.join(new_fold, 'IMCube.fits')
    rd.save_phasemap(cube_path, cube, isCube=True)
    # Copying the cmdMatrix and the ModesVector into the INTMAT Folder
    cmat = rd.readFits_data(os.path.join(ifFold, tn, 'cmdMatrix.fits'))
    mvec = rd.readFits_data(os.path.join(ifFold, tn, 'modesVector.fits'))
    pyfits.writeto(os.path.join(intMatFold, tn, 'cmdMatrix.fits'), cmat)
    pyfits.writeto(os.path.join(intMatFold, tn, 'modesVector.fits'), mvec)
    print(f"Cube saved in '{cube_path}'")
    return cube

def stackCubes(tnlist):
    """
    Stack the cubes sontained in the corresponding tracking number folder, creating
    a new cube, along with stacked command matrix and modes vector.

    Parameters
    ----------
    tnlist : list of str
        List containing the tracking number of the cubes to stack.

    Returns
    -------
    stacked_cube : masked_array
        Final cube, stacked along the 3th axis.
    """
    new_tn = timestamp.Timestamp.now()
    stacked_cube_fold = os.path.join(fn.INTMAT_ROOT_FOLDER, new_tn)
    os.mkdir(stacked_cube_fold)
    cubelist = []
    matrixList = []
    modesVectList = []
    for tn in tnlist:
        fold = os.path.join(intMatFold, tn)
        cube_name = os.path.join(fold, 'IMCube.fits')
        matrix_name = os.path.join(fold, 'cmdMatrix.fits')
        modesVec_name = os.path.join(fold, 'modesVector.fits')
        cubelist.append(rd.readFits_data(cube_name))
        matrixList.append(rd.readFits_data(matrix_name))
        modesVectList.append(rd.readFits_data(modesVec_name))
    stacked_cube = np.ma.dstack(cubelist)
    stacked_cmat = np.dstack(matrixList)
    stacked_mvec = np.dstack(modesVectList)
    save_cube = os.path.join(stacked_cube_fold, 'IMCube.fits')
    save_cmat = os.path.join(stacked_cube_fold, 'cmdMatrix.fits')
    save_mvec = os.path.join(stacked_cube_fold, 'modesVector.fits')
    rd.save_phasemap(save_cube, stacked_cube, isCube=True)
    pyfits.writeto(save_cmat, stacked_cmat)
    pyfits.writeto(save_mvec, stacked_mvec)
    print(f"Stacked cube and matrices saved in {new_tn}")

def iffRedux(tn, fileMat, ampVect, modeList, template, shuffle=0):
    """
    Reduction function that performs the push-pull analysis on each mode, saving
    out the final processed image for each mode.

    Parameters
    ----------
    tn : str
        Tracking number of the data in the OPDImages folder.
    fileMat : ndarray
        A matrix of images in string format, in which each row is a mode and the
        columns are its template realization.
    ampVect : float | ArrayLike
        Vector containing the amplitude for each commanded mode.
    modeList : int | ArrayLike
        Vector conaining the list of commanded modes.
    template : int | ArrayLike
        Template for the push-pull command actuation.
    shuffle : int, optional
        A value different from 0 activates the shuffle option, and the imput 
        value is the number of repetition for each mode's push-pull packet. The
        default is 0, which means the shuffle is OFF.
    """
    fold = os.path.join(ifFold, tn)
    nmodes = len(ampVect)
    # nframes= len(template)
    # if shuffle != 0:
    #     nframes = shuffle*3
    for i in range(0, nmodes):
        img = pushPullRedux(fileMat[i,:], template, shuffle)
        norm_img = img / (2*ampVect[i])
        img_name = os.path.join(fold, f"mode_{modeList[i]:04d}.fits")
        rd.save_phasemap(img_name, norm_img)

def pushPullRedux(fileVec, template, shuffle=0):
    r"""
    Performs the basic operation of processing PushPull data.
    
    Packs all mode's push-pull into a list and then performs the differential
    algorithm
    
    > $\sum_i \dfrac{img_i \cdot template_i - img_{i-1}\cdot template_{i-1}}{}$
    
    Parameters
    ----------
    fileVec : string | array
        It is a row in the fileMat (the organized matrix of the images filename), 
        corresponding to all the realizations of the same mode (or act), with a
        given template. If shuffle option has been used, the fileMat (and fileVec) 
        shall be reorganized before running the script.
    template: int | ArrayLike
        Template for the PushPull acquisition.
    shuffle: int, optional
        A value different from 0 activates the shuffle option, and the imput 
        value is the number of repetition for each mode's templated sampling. 
        The default value is 0, which means the shuffle option is OFF.

    Returns
    -------
    image: masked_array
        Final processed mode's image.
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
                master_mask = np.ma.mask_or(master_mask, master_mask2add)
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

def registrationRedux(fileMat, template=None):
    """
    Reduction function that performs the push-pull analysis on the registration
    data.
    
    Parameters
    ----------
    fileMat : ndarray
        A matrix of images in string format, in which each row is a mode and the
        columns are its template realization.

    Returns
    -------
    imgList : ArrayLike
        List of the processed registration images.
    """
    if template is None:
        _,infoR,_ = _getAcqInfo()
        template = infoR['template']
    nAct = fileMat.shape[0]
    imglist = []
    for i in range(0, nAct-1):
        img = pushPullRedux(fileMat[i,:], template)
        imglist.append(img)
    return imglist

def findFrameOffset(tn, imglist, actlist):
    """
    This function computes the position difference between the current frame and
    a reference one.
    
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

def getTriggerFrame(tn, amplitude=None):
    """
    Analyze the tracking number's images list and search for the trigger frame.

    Parameters
    ----------
    tn : str
        Tracking number of the data in the OPDImages folder.
    amplitude : int os float, optional
        Amplitude of the commanded trigger mode, which serves as the check value
        for finding the frame. If no value is passed it is loaded from the iffConfig.ini
        file.

    Returns
    -------
    trigFrame : int
        Index which identifies the trigger frame in the images folder file list.
        
    Raises
    ------
    RuntimeError
        Error raised if the file iteration goes beyon the expected trigger frame
        wich can be inferred through the number of trigger zeros in the iffConfig.ini
        file.
    """
    infoT, _, _ = _getAcqInfo()
    if amplitude is not None:
        infoT['amplitude'] = amplitude
    fileList = _getFileList(tn)
    img0 = rd.read_phasemap(fileList[0])
    go = i = 1
    while go !=0:
        if go > infoT['zeros']:
            raise RuntimeError('Heading Zeros exceeded, error!!')
        img1 = rd.read_phasemap(fileList[i])
        rr2check = zern.removeZernike(img1-img0,[1,2,3]).std()
        if rr2check > infoT['amplitude']/3:
            go=0
        else:
            i+=1
            go+=1
            img0 = img1
    trigFrame = i
    return trigFrame
    
def getRegFileMatrix(tn):
    """
    Search for the registration frames in the images file list, and creates the
    registration file matrix.

    Parameters
    ----------
    tn : str
        Tracking number of the data in the OPDImages folder.

    Returns
    -------
    regEnd : int
        Index which identifies the last registration frame in the images file
        list.
    regMat : ndarray
        A matrix of images in string format, containing the registration frames.
        It has shape (registration_modes, n_push_pull).
    """
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
    """
    Creates the iffMat

    Parameters
    ----------
    tn : str
        Tracking number of the data in the OPDImages folder.

    Returns
    -------
    iffMat : ndarray
        A matrix of images in string format, conatining all the images for the
        IFF acquisition, that is all the modes with each push-pull realization.
        It has shape (modes, n_push_pull)
    """
    fileList    = _getFileList(tn)
    _,_,infoIF  = _getAcqInfo()
    regEnd, _   = getRegFileMatrix(tn)
    iffList     = fileList[regEnd+infoIF['zeros']:]
    iffMat      = np.reshape(iffList, 
                             (len(infoIF['modes']), len(infoIF['template'])))
    return iffMat

def _getFileList(tn, fold=None):
    """
    Returns the file list of a given tracking number datapath.

    Parameters
    ----------
    tn : str
        Tracking number of the data in the OPDImages folder.
    fold : str, optional
        Folder in which searching for the tracking number. If None, the default
        folder is OPD_IMAGES_ROOT_FOLDER.

    Returns
    -------
    fl : list
        List of sorted files inside the folder.
    """
    if fold is None:
        fl = sorted([os.path.join(imgFold, (tn+'/'+image)) \
                     for image in os.listdir(os.path.join(imgFold, tn))])
    else:
        fl = sorted([os.path.join(fold, (tn+'/'+image)) \
                     for image in os.listdir(os.path.join(fold, tn))])
    return fl

def _getAcqPar(tn):
    """
    Reads ad returns the acquisition parameters from fits files.

    Parameters
    ----------
    tn : str
        Tracking number of the data in the OPDImages folder.

    Returns
    -------
    ampVector : float | ArrayLike
        Vector containg the amplitude of each commanded mode.
    modesVector : int | ArrayLike
        Vector containing the list of commanded modes.
    template : int | ArrayLike
        Sampling template ampplied on each mode.
    indexList : int | ArrayLike
        Indexing of the modes inside the commanded matrix.
    registrationActs : int | ArrayLike
        Vector containing the commanded actuators for the registration.
    shuffle : int
        Shuffle information. If it's nor 0, the values indicates the number of
        template sampling repetition for each mode.
    """
    base = os.path.join(ifFold, tn)
    ampVector        = rd.readFits_data(os.path.join(base, ampVecFile))
    template         = rd.readFits_data(os.path.join(base, templateFile))
    modesVector      = rd.readFits_data(os.path.join(base, modesVecFile))
    indexList        = rd.readFits_data(os.path.join(base, indexListFile))
    registrationActs = rd.readFits_data(os.path.join(base, regisActFile))
    with open(os.path.join(base, shuffleFile), 'r') as shf:
        shuffle = int(shf.read())
    return ampVector, modesVector, template,indexList, registrationActs,shuffle

def _getAcqInfo():
    """
    Returns the information read from the iffConfig.ini file.

    Returns
    -------
    infoT : dict
        Information read about the TRIGGER options.
    infoR : dict
        Information read about the REGISTRATION options.
    infoIF : dict
        Information read about the IFFUNC option.
    """
    infoT = read_iffconfig.getConfig('TRIGGER')
    infoR = read_iffconfig.getConfig('REGISTRATION')
    infoIF = read_iffconfig.getConfig('IFFUNC')
    return infoT, infoR, infoIF

# def _ampReorganization(ampVector):
#     pass

# def _modesReorganization(modesVector):
#     pass
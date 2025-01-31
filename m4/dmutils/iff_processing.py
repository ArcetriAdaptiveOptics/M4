"""
Author(s):
----------
    - Pietro Ferraiuolo
    - Runa Briguglio
    
Written in June 2024

Description
-----------
Module containing all the functions necessary to process the data acquired for 
the Influence Function measurements done on M4.

High-level Functions
--------------------
process(tn, registration=False) 
    Function that processes the data contained in the OPDImages/tn folder. by p
    erforming the differential algorithm, it procudes fits images for each comm
    anded mode into the IFFunctions/tn folder, and creates a cube from these in
    to INTMatrices/tn. If 'registration is not False', upon createing the cube,
    the registration algorithm is performed.

stackCubes(tnlist)
    Function that, given as imput a tracking number list containing cubes data,
    will stack the found cubes into a new one with a new tracking number, into 
    INTMatrices/new_tn. A 'flag.txt' file will be created to give more informat
    ion on the process.

Notes
-----
In order for the module to work properly, the tower initialization must be run
so that the folder names configuration file is populated. 
From the IPython console

>>> run '/path/to/m4/initOTT.py'
>>> from m4.dmutils import iff_processing as ifp

Example
-------
>>> tn1 = '20160516_114916'
>>> tn2 = '20160516_114917' # A copy of tn1 (simulated) data
>>> ifp.process(tn1)
Cube saved in '.../m4/data/M4Data/OPTData/INTMatrices/20160516_114916/IMcube.fits'
>>> ifp.process(tn2)
Cube saved in '.../m4/data/M4Data/OPTData/INTMatrices/20160516_114917/IMcube.fits'
>>> tnlist = [tn1, tn2]
>>> ifp.stackCubes(tnlist)
Stacekd cube and matrices saved in '.../m4/data/M4Data/OPTData/INTMatrices/'new_tn'/IMcube.fits'
"""
import os, shutil
import configparser
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration import read_iffconfig, update_folder_paths as ufp
from m4.ground import read_data as rd, zernike as zern, timestamp
from m4.utils import osutils as osu
from scripts.misc.IFFPackage import actuator_identification_lib as fa
config = configparser.ConfigParser()
fn = ufp.folders
imgFold     = fn.OPD_IMAGES_ROOT_FOLDER
ifFold      = fn.IFFUNCTIONS_ROOT_FOLDER
intMatFold  = fn.INTMAT_ROOT_FOLDER
confFold    = fn.CONFIGURATION_ROOT_FOLDER
frameCenter = [200,200]
ampVecFile     = 'ampVector.fits'
modesVecFile   = 'modesVector.fits'
templateFile   = 'template.fits'
regisActFile   = 'regActs.fits'
shuffleFile    = 'shuffle.dat'
indexListFile  = 'indexList.fits'
coordfile      = '' #TODO
cmdMatFile     = 'cmdMatrix.fits'
cubeFile       = 'IMCube.fits'
flagFile       = 'flag.txt'

def process(tn, register:bool=False, save_cube:bool=False):
    """
    High level function with processes the data contained in the given tracking
    number OPDimages folder, performing the differential algorithm and saving 
    the final cube.

    Parameters
    ----------
    tn : str
        Tracking number of the data in the OPDImages folder.
    register : bool, optional
        Parameter which enables the registration option. The default is False.
    """
    ampVector,modesVector,template,_,registrationActs,shuffle = _getAcqPar(tn)
    _,regMat = getRegFileMatrix(tn)
    modesMat = getIffFileMatrix(tn)
    actImgList = registrationRedux(regMat, template) #FIXME
    modesMatReorg = _modesReorganization(modesMat)
    iffRedux(tn, modesMatReorg, ampVector, modesVector, template, shuffle)
    if register:
        dx = findFrameOffset(tn, actImgList, registrationActs)
    else:
        dx = register
    if save_cube:
        saveCube(tn, dx)

def saveCube(tn, register=False):
    """
    Creates and save a cube from the fits files contained in the tn folder, 
    along with the command matrix and the modes vector fits.

    Parameters
    ----------
    tn : str
        Tracking number of the IFFunctions data folder from which create the cu
        be.
    register : int or tuple, optional
        If not False, and int or a tuple of int must be passed as value, and 
        the registration algorithm is performed on the images before stacking them
        into the cube. Default is False.
        
    Returns
    -------
    cube : masked_array
        Data cube of the images, with shape (npx, npx, nmodes).
    """
    old_fold = os.path.join(ifFold, tn)
    filelist = osu.getFileList(fold=old_fold, key="mode_")
    cube = osu.createCube(filelist, register=register)
    # Saving the cube
    new_fold = os.path.join(intMatFold, tn)
    os.mkdir(new_fold)
    cube_path = os.path.join(new_fold, cubeFile)
    rd.save_phasemap(cube_path, cube)
    # Copying the cmdMatrix and the ModesVector into the INTMAT Folder
    cmat = rd.readFits_data(os.path.join(ifFold, tn, 'cmdMatrix.fits'))
    mvec = rd.readFits_data(os.path.join(ifFold, tn, 'modesVector.fits'))
    pyfits.writeto(os.path.join(intMatFold, tn, 'cmdMatrix.fits'), cmat)
    pyfits.writeto(os.path.join(intMatFold, tn, 'modesVector.fits'), mvec)
    with open(os.path.join(intMatFold, tn, flagFile), 'w', encoding='utf-8') as f:
        f.write(f"Cube created from '{old_fold}' data.\n \n \n")
    print(f"Cube saved in '{cube_path}'")
    print(f"Shape: {cube.shape}")
    return cube

def stackCubes(tnlist):
    """
    Stack the cubes sontained in the corresponding tracking number folder, creating
    a new cube, along with stacked command matrix and modes vector.

    Parameters
    ----------
    tnlist : list of str
        List containing the tracking numbers of the cubes to stack.

    Returns
    -------
    stacked_cube : masked_array
        Final cube, stacked along the 3th axis.
    """
    new_tn              = timestamp.Timestamp.now()
    stacked_cube_fold   = os.path.join(fn.INTMAT_ROOT_FOLDER, new_tn)
    os.mkdir(stacked_cube_fold)
    cube_parameters     = _getCubeList(tnlist)
    flag                = _checkStackedCubes(tnlist)
    # Stacking the cube and the matrices
    stacked_cube        = np.ma.dstack(cube_parameters[0])
    stacked_cmat        = np.dstack(cube_parameters[1])
    stacked_mvec        = np.dstack(cube_parameters[2])
    # Saving everithing to a new file into a new tn
    save_cube           = os.path.join(stacked_cube_fold, cubeFile)
    save_cmat           = os.path.join(stacked_cube_fold, 'cmdMatrix.fits')
    save_mvec           = os.path.join(stacked_cube_fold, 'modesVector.fits')
    rd.save_phasemap(save_cube, stacked_cube)
    pyfits.writeto(save_cmat, stacked_cmat)
    pyfits.writeto(save_mvec, stacked_mvec)
    with open(os.path.join(stacked_cube_fold, flagFile), 'w', encoding='UTF-8') as file:
        flag.write(file)
    print(f"Stacked cube and matrices saved in {new_tn}")
    
def filterZernikeCube(tn, zern_modes:list=None, save:bool=True):
    """
    Function which filters out the desired zernike modes from a cube.

    Parameters
    ----------
    tn : str
        Tracking number of the cube to filter.
    zern_modes : list, optional
        List of zernike modes to filter out. The default is [1,2,3] 
        (piston, tip and tilt).

    Returns
    -------
    ffcube : masked array
        Filtered cube.
    """
    new_tn = os.path.join(intMatFold, timestamp.Timestamp.now())
    os.mkdir(new_tn)
    oldCube = os.path.join(intMatFold, tn, cubeFile)
    ocFlag  = os.path.join(intMatFold, tn, flagFile)
    newCube = os.path.join(new_tn, cubeFile)
    ocFlag  = os.path.join(intMatFold, tn, flagFile)
    newFlag = os.path.join(new_tn, flagFile)
    CmdMat = os.path.join(intMatFold, tn, cmdMatFile)
    ModesVec = os.path.join(intMatFold, tn, modesVecFile)
    cube = rd.readFits_maskedImage(oldCube)
    zern2filter = zern_modes if zern_modes is not None else [1,2,3]
    fcube = []
    for i in range(cube.shape[-1]):
        filtered = zern.removeZernike(cube[:,:,i], zern2filter)
        fcube.append(filtered)
    ffcube = np.ma.dstack(fcube)
    if save:
        rd.save_phasemap(newCube, ffcube)
        shutil.copyfile(CmdMat, os.path.join(new_tn, cmdMatFile))
        shutil.copyfile(ModesVec, os.path.join(new_tn, modesVecFile))
        with open(ocFlag, 'r', encoding='utf-8') as oflag:
            flag = oflag.readlines()
        flag.pop(-1)
        flag += f"Zernike modes filtered = {zern2filter}"
        with open(newFlag, 'w', encoding='utf-8') as nflag:
            nflag.writelines(flag)
        print(f"Filtered cube saved at {new_tn}")
    return ffcube, new_tn.split('/')[-1]

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
    nmodes = len(modeList)
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
    for i in range(0, template.shape[0]):
        ima = rd.read_phasemap(fileVec[i])
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
                opd2add = image_list[i*3+x]*template[x] + \
                            image_list[i*3+x-1]*template[x-1]
                master_mask2add = np.ma.mask_or(image_list[i*3+x].mask,
                                                image_list[i*3+x-1].mask)
                if (i == 0 and x == 1):
                    master_mask = master_mask2add
                else:
                    master_mask = np.na.mask_or(master_mask, master_mask2add)
                image += opd2add
    image = np.ma.masked_array(image, mask=master_mask) / (template.shape[0]-1)#!!!
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
    nActs = fileMat.shape[0]
    imglist = []
    for i in range(0, nActs-1):
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
    xy = fa.findFrameCoord(imglist, actlist, actCoord)
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
    fileList = osu.getFileList(tn)
    img0 = rd.read_phasemap(fileList[0])
    go = i = 1
    # add the condition where if there are not trigger frames the code is skipped and the 
    # the rest is handled with care
    while go !=0:
        thresh = infoT['amplitude']/3
        img1 = rd.read_phasemap(fileList[i])
        rr2check = zern.removeZernike(img1-img0,[1,2,3]).std()
        if go > infoT['zeros']:
            raise RuntimeError(f"Frame {go}. Heading Zeros exceeded: std= {rr2check:.2e} < {thresh:.2e} =Amp/3")
        if rr2check > thresh:
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
    fileList    = osu.getFileList(tn)
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
    fileList    = osu.getFileList(tn)
    _,_,infoIF  = _getAcqInfo()
    regEnd, _   = getRegFileMatrix(tn)
    iffList     = fileList[regEnd+infoIF['zeros']:]
    iffMat      = np.reshape(iffList,
                             (len(infoIF['modes']), len(infoIF['template'])))
    return iffMat

def _getCubeList(tnlist):
    cubeList = []
    matrixList = []
    modesVectList = []
    for tn in tnlist:
        fold = os.path.join(intMatFold, tn)
        cube_name = os.path.join(fold, 'IMCube.fits')
        matrix_name = os.path.join(fold, 'cmdMatrix.fits')
        modesVec_name = os.path.join(fold, 'modesVector.fits')
        cubeList.append(rd.readFits_data(cube_name))
        matrixList.append(rd.readFits_data(matrix_name))
        modesVectList.append(rd.readFits_data(modesVec_name))
    return cubeList, matrixList, modesVectList

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
    with open(os.path.join(base, shuffleFile), 'r', encoding='UTF-8') as shf:
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

def _checkStackedCubes(tnlist):
    """
    Inspect the cubes to stack, to check whether there are shared modes, or not.

    Parameters
    ----------
    tnlist : list of str
        List containing the tracking number of the cubes to stack.

    Returns
    -------
    flag : dict
        Dictionary containing the flagging information about the stacked cube, 
        to be later dump into the 'flag.txt' file.
    """
    _,_,modesVectList = _getCubeList(tnlist)
    nmodes = len(modesVectList[0])
    nvects = len(modesVectList)
    for i in range(nvects):
        for j in range(i + 1, nvects):
            common_modes = set(modesVectList[i]).intersection(modesVectList[j])
    c_nmodes = len(common_modes)
    if c_nmodes in range(1,nmodes):
        flag = __shared_modes(tnlist,modesVectList)
    elif c_nmodes==nmodes:
        flag = __averaged(tnlist,modesVectList)
    else:
        flag = __stacked(tnlist,modesVectList)
    return flag

def __stacked(tnlist,modesVectList):
    """
    Creates the dictionary to dump into the 'flag.txt' file accordingly to 
    sequentially stacked cubes with no repeated modes.

    Parameters
    ----------
    tnlist : list of str
        List containing the tracking number of the cubes to stack.
    modesVectList : list of ndarray
        A list containing the modes vectors for each cube.

    Returns
    -------
    config : dict
        Dictionary containing the flagging information about the stacked cube.
    """
    c_type = 'Sequentially stacked cubes'
    text=''
    for i, tn in enumerate(tnlist):
        text += f"""{tn}, modes {list(modesVectList[i])} \\
           """
    flag = {
        'Flag': {
            'Cube type': c_type,
            'Source cubes': text,
            }
        }
    config['Flag'] = {}
    for key,value in flag['Flag'].items():
        config['Flag'][key] = value
    return config

def __averaged(tnlist,modesVectList):
    """
    Creates the dictionary to dump into the 'flag.txt' file accordingly to 
    averaged cubes with same commanded modes.

    Parameters
    ----------
    tnlist : list of str
        List containing the tracking number of the cubes to stack.
    modesVectList : list of ndarray
        A list containing the modes vectors for each cube.

    Returns
    -------
    config : dict
        Dictionary containing the flagging information about the stacked cube.
    """
    c_type = 'Mean of cubes'
    text=''
    for i, tn in enumerate(tnlist):
        text += f"""{tn}, modes {list(modesVectList[i])} \\
           """
    flag = {
        'Flag': {
            'Cube type': c_type,
            'Source cubes': text,
            }
        }
    config['Flag'] = {}
    for key,value in flag['Flag'].items():
        config['Flag'][key] = value
    return config

def __shared_modes(tnlist,modesVectList):
    """
    Creates the dictionary to dump into the 'flag.txt' file accordingly to 
    stacked cubes with some shared mode, which should be treated carefully.

    Parameters
    ----------
    tnlist : list of str
        List containing the tracking number of the cubes to stack.
    modesVectList : list of ndarray
        A list containing the modes vectors for each cube.

    Returns
    -------
    config : dict
        Dictionary containing the flagging information about the stacked cube.
    """
    c_type = '!!!Warning: repeated modes in stacked cube'
    text=''
    for i, tn in enumerate(tnlist):
        text += f"""{tn}, modes {list(modesVectList[i])} \\
           """
    flag = {
        'Flag': {
            'Cube type': c_type,
            'Source cubes': text,
            }
        }
    config['Flag'] = {}
    for key,value in flag['Flag'].items():
        config['Flag'][key] = value
    return config

#TODO
def _ampReorganization(ampVector):
    reorganizaed_amps = ampVector
    return reorganizaed_amps

def _modesReorganization(modesVector):
    reorganizaed_modes = modesVector
    return reorganizaed_modes

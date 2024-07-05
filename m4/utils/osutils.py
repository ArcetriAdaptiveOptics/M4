"""
Author(s) 
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
Utils module, containing all the useful interaction between the software and the
machine OS.
"""
import os
from m4.configuration import config_folder_names as fn
from m4.ground.read_4DConfSettingFile import ConfSettingReader

optdata = fn.OPT_DATA_FOLDER
opdimg = fn.OPD_IMAGES_ROOT_FOLDER

def findTracknum(tn):
    """
    Search for the tracking number given in input within all the data path subf
    olders

    Parameters
    ----------
    tn : str
        Tracking number to be searched.

    Returns
    -------
    tn_path : list of str
        List containing all the folders (within the OPTData path) in which the 
        tracking number is present.
    """
    tn_path = []
    for fold in os.listdir(optdata):
        search_fold = os.path.join(optdata, fold)
        if tn in os.listdir(search_fold):
            tn_path.append(fold)
    return tn_path

def getFileList(tn, fold=None, key:str=None):
    """
    Returns the file list of a given tracking number datapath.

    Parameters
    ----------
    tn : str
        Tracking number of the data in the OPDImages folder.
    fold : str, optional
        Folder in which searching for the tracking number. If None, the default
        folder is the OPD_IMAGES_ROOT_FOLDER.
    key : str, optional
        A key which identify specific files to return

    Returns
    -------
    fl : list os str
        List of sorted files inside the folder.

    Examples
    --------
    Here are some examples regarding the use of the 'key' argument. Let's say w
    e need a list of files inside ''tn = '20160516_114916' '' in the IFFunction
    s folder. 

        >>> iffold = '.../M4/m4/data/M4Data/OPTData/IFFunctions'
        >>> tn = '20160516_114916'
        >>> getFileList(tn, fold=iffold)
        ['.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/cmdMatrix.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0000.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0001.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0002.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0003.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/modesVector.fits']
        
    Let's suppose we want only the list of 'mode_000x.fits' files:
    
        >>> getFileList(tn, fold=fold, key='mode_')
        ['.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0000.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0001.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0002.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0003.fits']
    
    Notice that, in this specific case, it was necessary to include the undersc
    ore after 'mode' to exclude the 'modesVector.fits' file from the list.
    """
    if fold is None:
        fl = sorted([os.path.join(opdimg, (tn+'/'+file)) \
                     for file in os.listdir(os.path.join(opdimg, tn))])
    else:
        fl = sorted([os.path.join(fold, (tn+'/'+file)) \
                     for file in os.listdir(os.path.join(fold, tn))])
    if key is not None:
        try:
            selected_list = []
            for file in fl:
                if key in file:
                    selected_list.append(file)
        except TypeError:
            raise TypeError("'key' argument must be a string")
        fl = selected_list
    return fl

def tnRange(tn0, tn1):
    """
    Returns the list of tracking numbers between tn0 and tn1, within the same f
    older, if they both exist in it.

    Parameters
    ----------
    tn0 : str
        Starting tracking number.
    tn1 : str
        Finish tracking number.

    Returns
    -------
    tnMat : list of str
        A list or a matrix of tracking number in between the start and finish 
        ones.

    Raises
    ------
    Exception
        An exception is raised if the two tracking numbers are not found in the
        same folder
    """
    tn0_fold = findTracknum(tn0)
    tn1_fold = findTracknum(tn1)
    if len(tn0_fold)==1 and len(tn1_fold)==1:
        if tn0_fold[0]==tn1_fold[0]:
            fold = os.path.join(optdata, tn0_fold[0])
            tn_folds = sorted(os.listdir(fold))
            id0 = tn_folds.index(tn0)
            id1 = tn_folds.index(tn1)
            tnMat = [os.path.join(fold, tn) for tn in tn_folds[id0:id1+1]]
        else:
            raise Exception("The tracking numbers are in different foldes")
    else:
        tnMat = []
        for ff in tn0_fold:
            if ff in tn1_fold:
                fold = os.path.join(optdata, ff)
                tn_folds = sorted(os.listdir(fold))
                id0 = tn_folds.index(tn0)
                id1 = tn_folds.index(tn1)
                tnMat.append([os.path.join(fold, tn) for tn in tn_folds[id0:id1+1]])
    return tnMat

def rename4D(folder):
    """
    Renames the produced 'x.4D' files into '0000x.4D'
    
    Parameters
    ----------
    folder : str
        The folder where the 4D data is stored.
    """
    fold = os.path.join(opdimg, folder)
    files = os.listdir(fold)
    for file in files:
        if file.endswith('.4D'):
            num_str = file.split('.')[0]
            if num_str.isdigit():
                num = int(num_str)
                new_name = f"{num:05d}.4D"
                old_file = os.path.join(fold, file)
                new_file = os.path.join(fold, new_name)
                os.rename(old_file, new_file)

def getConf4DSettingsPath(tn):
    """
    Returns the path of the '4DSettings.ini' camera configuration file.

    Parameters
    ----------
    tn : str
        The tracking number of the configuration to search.

    Returns
    -------
    config : str
        Complete file path of the configuration file.
    """
    path = getFileList(tn, key='4DSetting')
    return path

def getCameraSettings(tn):
    """
    Return a dictionary with all the camera settings loaded from the '4DSetting.ini'
    file

    Parameters
    ----------
    tn : str
        Tracking number where to search the camera configuration file.

    Returns
    -------
    camera_settings : dict
        Dictionary with the useful camera parameters at the moment of acquisition.
        These are:
            width (px) : The camera width in pixel scale
            height (px) : The camera height in pixel scale
            x-offset : The offset of the camera in the x (width) axis
            y-offset : The offset of the camera in the y (height) axis
            frame rate : The frame rate of acquisition of the camera
    """
    conf_path = getConf4DSettingsPath(tn)[0]
    setting_reader = ConfSettingReader(conf_path)
    camera_settings = {}
    camera_settings['width (px)'] = setting_reader.getImageWidhtInPixels()
    camera_settings['height (px)'] = setting_reader.getImageHeightInPixels()
    camera_settings['x-offset'] = setting_reader.getOffsetX()
    camera_settings['y-offset'] = setting_reader.getOffsetY()
    camera_settings['frame rate'] = setting_reader.getFrameRate()
    return camera_settings

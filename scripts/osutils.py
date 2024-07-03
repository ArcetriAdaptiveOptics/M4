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
    fl : list
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
        
    Let's suppose we want only the list of mode fits files.
    
        >>> getFileList(tn, fold=fold, key='mode_')
        ['.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0000.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0001.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0002.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0003.fits']
    
    It was necessary, in this case, to include the underscore to not include th
    e 'modesVector.fits' file in the list
    """
    if fold is None:
        fl = sorted([os.path.join(opdimg, (tn+'/'+image)) \
                     for image in os.listdir(os.path.join(opdimg, tn))])
    else:
        fl = sorted([os.path.join(fold, (tn+'/'+image)) \
                     for image in os.listdir(os.path.join(fold, tn))])
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
    

    Parameters
    ----------
    tn0 : TYPE
        DESCRIPTION.
    tn1 : TYPE
        DESCRIPTION.

    Returns
    -------
    tnMat : TYPE
        DESCRIPTION.

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
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
import numpy as np
from arte.dataelab.base_file_walker import AbstractFileNameWalker
from m4.ground import read_data as rd
from m4.ground.timestamp import Timestamp
from m4.configuration import config_folder_names as fn
from m4.ground.read_4DConfSettingFile import ConfSettingReader
from pathlib import Path

OPTDATA = fn.OPT_DATA_FOLDER
OPDIMG = fn.OPD_IMAGES_ROOT_FOLDER

class FileWalker(AbstractFileNameWalker):

    def __init__(self, data_root_dir=OPTDATA):
        self._data_root_dir = Path(data_root_dir)

    def snapshot_dir(self, tn=None):
        if isinstance(tn, str):
            tn = Timestamp(tn)
        return Path(self._data_root_dir, tn)

    # def rsi_measure_path(self, tn):
    #    return self.snapshot_dir(tn) / 'imas.fits'

    def findTracknum(self, tn):
        tn_path = []

        for root, dirs, files in Path(self._data_root_dir).walk(on_error=print):
            found_tn_dirs = list(filter(lambda l: tn in l, dirs))
            found_tn_files = list(filter(lambda l: tn in l, files))

#            if len(found_tn_files) > 0:
#                # tn_path.append(str(root)[-15:])
#                tn_path.append(root)
#            if len(found_tn_dirs) > 0:
#                tn_path.append(Path(root, found_tn_dirs[0]))
            if len(found_tn_files) > 0:
                [tn_path.append(Path(root, ff)) for ff in found_tn_files]

            if len(found_tn_dirs) > 0:
                [tn_path.append(Path(root, dd)) for dd in found_tn_dirs]

        return tn_path

    def find_tag_between_dates(self, tn_start, tn_stop):
        if isinstance(tn_start, Timestamp):
            tn_start = tn_start.asNowString()
        if isinstance(tn_stop, Timestamp):
            tn_stop = tn_stop.asNowString()

        tn_path = []
        for root, dirs, files in Path(self._data_root_dir).walk(on_error=print):

            found_tn_dirs = list(
                filter(lambda l: tn_stop >= l, filter(lambda l: tn_start <= l, dirs)))
            found_tn_files = list(
                filter(lambda l: tn_stop >= l, filter(lambda l: tn_start <= l, files)))

            if len(found_tn_files) > 0:
                [tn_path.append(Path(root, ff)) for ff in found_tn_files]

            if len(found_tn_dirs) > 0:
                [tn_path.append(Path(root, dd)) for dd in found_tn_dirs]

        return sorted(tn_path, key=lambda l: str(l))

#        tn_list = tnRange(tn_start.asNowString(), tn_stop.asNowString())
#        return sorted([tn for tn in tn_list if tn >= str(tn_start) and tn <= str(tn_stop)])



def findTracknum(tn, complete_path:bool=False):

    """
    Search for the tracking number given in input within all the data path subfolders.

    Parameters
    ----------
    tn : str
        Tracking number to be searched.
    complete_path : bool, optional
        Option for wheter to return the list of full paths to the folders which
        contain the tracking number or only their names.

    Returns
    -------
    tn_path : list of str
<<<<<<< Updated upstream

        List containing all the folders (within the OPTData path) in which the 
=======
        List containing all the folders (within the OPTData path) in which the
>>>>>>> Stashed changes
        tracking number is present, sorted in alphabetical order.

    """
    tn_path = []
    for fold in os.listdir(OPTDATA):
        search_fold = os.path.join(OPTDATA, fold)
        if tn in os.listdir(search_fold):
            if complete_path:
                tn_path.append(search_fold)
            else:
                tn_path.append(fold)
    return sorted(tn_path)


def getFileList(tn=None, fold=None, key:str=None):

    """
    Search for files in a given tracking number or complete path, sorts them and
    puts them into a list.

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

    How to Use it
    -------------
    If the complete path for the files to retrieve is available, then this function
    should be called with the 'fold' argument set with the path, while 'tn' is
    defaulted to None.

    In any other case, the tn must be given: it will search for the tracking
    number into the OPDImages folder, but if the search has to point another
    folder, then the fold argument comes into play again. By passing both the
    tn (with a tracking number) and the fold argument (with only the name of the
    folder) then the search for files will be done for the tn found in the
    specified folder. Hereafter there is an example, with the correct use of the
    key argument too.

    Examples
    --------
<<<<<<< Updated upstream
    Here are some examples regarding the use of the 'key' argument. Let's say w
    e need a list of files inside ''tn = '20160516_114916' '' in the IFFunction
=======
    Here are some examples regarding the use of the 'key' argument. Let's say we
    need a list of files inside ''tn = '20160516_114916' '' in the IFFunction
>>>>>>> Stashed changes
    s folder.

        >>> iffold = 'IFFunctions'
        >>> tn = '20160516_114916'
        >>> getFileList(tn, fold=iffold)
        ['.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/cmdMatrix.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0000.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0001.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0002.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0003.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/modesVector.fits']

    Let's suppose we want only the list of 'mode_000x.fits' files:

<<<<<<< Updated upstream
        >>> getFileList(tn, fold=fold, key='mode_')
=======
        >>> getFileList(tn, fold=iffold, key='mode_')
>>>>>>> Stashed changes
        ['.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0000.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0001.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0002.fits',
         '.../M4/m4/data/M4Data/OPTData/IFFunctions/20160516_114916/mode_0003.fits']

<<<<<<< Updated upstream
    Notice that, in this specific case, it was necessary to include the undersc
    ore after 'mode' to exclude the 'modesVector.fits' file from the list.
=======
    Notice that, in this specific case, it was necessary to include the underscore
    after 'mode' to exclude the 'modesVector.fits' file from the list.
>>>>>>> Stashed changes
    """

    if tn is None and fold is not None:
        fl = sorted([os.path.join(fold, file) \
                     for file in os.listdir(fold)])
    else:
        if fold is None:
            fl = sorted([os.path.join(OPDIMG, tn, file) \
                         for file in os.listdir(os.path.join(OPDIMG, tn))])
        else:
<<<<<<< Updated upstream
            fl = sorted([os.path.join(fold, tn, file) \
                         for file in os.listdir(os.path.join(fold, tn))])

=======
            try:
                for the_folds in findTracknum(tn, complete_path=True):
                    if fold in the_folds:
                        path = os.path.join(the_folds, tn)
                        fl = sorted([os.path.join(path, file) \
                                                 for file in os.listdir(path)])
                    else: 
                        raise Exception
            except Exception as exc:
                raise FileNotFoundError(f"Invalid Path: no data found for '.../{fold}/{tn}'") from exc
>>>>>>> Stashed changes
    if key is not None:
        try:
            selected_list = []
            for file in fl:
                if key in file:
                    selected_list.append(file)
        except TypeError as err:
            raise TypeError("'key' argument must be a string") from err
        fl = selected_list
    return fl


def tnRange(tn0, tn1):
    """
    Returns the list of tracking numbers between tn0 and tn1, within the same
    folder, if they both exist in it.

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
    if len(tn0_fold) == 1 and len(tn1_fold) == 1:
        if tn0_fold[0] == tn1_fold[0]:
            fold = os.path.join(OPTDATA, tn0_fold[0])
            tn_folds = sorted(os.listdir(fold))
            id0 = tn_folds.index(tn0)
            id1 = tn_folds.index(tn1)
            tnMat = [os.path.join(fold, tn) for tn in tn_folds[id0:id1+1]]
        else:
            raise FileNotFoundError(
                "The tracking numbers are in different foldes")
    else:
        tnMat = []
        for ff in tn0_fold:
            if ff in tn1_fold:
                fold = os.path.join(OPTDATA, ff)
                tn_folds = sorted(os.listdir(fold))
                id0 = tn_folds.index(tn0)
                id1 = tn_folds.index(tn1)
                tnMat.append([os.path.join(fold, tn)
                             for tn in tn_folds[id0:id1+1]])
    return tnMat


def rename4D(folder):
    """
    Renames the produced 'x.4D' files into '0000x.4D'

    Parameters
    ----------
    folder : str
        The folder where the 4D data is stored.
    """
    fold = os.path.join(OPDIMG, folder)
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
            pixel format : The format of the interferometer pixels
    """
    conf_path = getConf4DSettingsPath(tn)[0]
    setting_reader = ConfSettingReader(conf_path)
    camera_settings = {}
    camera_settings['width (px)'] = setting_reader.getImageWidhtInPixels()
    camera_settings['height (px)'] = setting_reader.getImageHeightInPixels()
    camera_settings['x-offset'] = setting_reader.getOffsetX()
    camera_settings['y-offset'] = setting_reader.getOffsetY()
    camera_settings['frame rate'] = setting_reader.getFrameRate()
    camera_settings['pixel format'] = setting_reader.getPixelFormat()
    return camera_settings


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



def createCube(filelist, register=False):

    """
    Creates a cube of images from an images file list

    Parameters
    ----------
    filelist : list of str
        List of file paths to the images/frames to be stacked into a cube.
    register : int or tuple, optional
        If not False, and int or a tuple of int must be passed as value, and
        the registration algorithm is performed on the images before stacking them
        into the cube. Default is False.

    Returns
    -------
    cube : ndarray
        Data cube containing the images/frames stacked.
    """
    cube_list = []
    for imgfits in filelist:
        image = rd.readFits_maskedImage(imgfits)
        if register is not False:
            image = np.roll(image, register)
        cube_list.append(image)
    cube = np.ma.dstack(cube_list)
    return cube

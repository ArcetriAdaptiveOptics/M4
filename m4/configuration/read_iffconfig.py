'''
Author(s)
    -P.Ferraiuolo
    -R.Briguglio
    
Written in 06/2024
'''
import os
import configparser
import json
import numpy as np
import shutil
import m4.configuration.config_folder_names as fn

config=configparser.ConfigParser()
cfoldname       = fn.CONFIGURATION_ROOT_FOLDER
iff_configFile  = 'iffConfig.ini'
nzeroName       = 'numberOfZeros'
modeIdName      = 'modeId'
modeAmpName     = 'modeAmp'
templateName    = 'template'
modalBaseName   = 'modalBase'

def getConfig(key, bpath=cfoldname):
    """
    Reads the configuration file for the IFF acquisition.
    The key passed is the block of information retrieved

    Parameters
    ----------
    key : str
        Key value of the block of information to read. Can be
            - 'TRIGGER'
            - 'REGISTRATION'
            - 'IFFUNC'
    bpath : str, OPTIONAL
        Base path of the file to read. Default points to the Configuration root
        folder
            
    Returns
    -------
    info : dict
        A dictionary containing all the configuration file's info:
            - nzeros
            - modeId
            - modeAmp 
            - template
            - modalBase 
    """
    fname = os.path.join(bpath, iff_configFile)
    config.read(fname)
    cc = config[key]
    nzeros      = int(cc[nzeroName])
    modeId_str  = cc[modeIdName]
    try:
        modeId = np.array(json.loads(modeId_str))
    except json.JSONDecodeError:
        modeId = np.array(eval(modeId_str))
    modeAmp     = float(cc[modeAmpName])
    modalBase   = cc[modalBaseName]
    template    = np.array(json.loads(cc[templateName]))
    info = {'zeros': nzeros,
            'modes': modeId,
            'amplitude': modeAmp,
            'template': template,
            'modalBase': modalBase
        }
    return info

def copyConfingFile(tn, old_path=cfoldname):
    """
    This function copies the configuration file to the new folder created for the
    IFF data, to keep record of the configuration used on data acquisition.

    Parameters
    ----------
    tn : str
        Tracking number of the new data.
    old_path : str, OPTIONAL
        Base path of the file to read. Default points to the Configuration root
        folder.

    Returns
    -------
    res : str
        String containing the path where the file has been copied
    """
    fname = os.path.join(old_path, iff_configFile)
    nfname= os.path.join(fn.IFFUNCTIONS_ROOT_FOLDER, tn, iff_configFile)
    res = shutil.copy2(fname, nfname)
    return f"{iff_configFile} copied to {res}"


def updateConfigFile(key, info, bpath=cfoldname):
    """
    Updates the configuration file for the IFF acquisition.
    The key passed is the block of information to update

    Parameters
    ----------
    key : str
        Key value of the block of information to update. Can be
            - 'TRIGGER'
            - 'REGISTRATION'
            - 'IFFUNC'
    info : dict
        A dictionary containing all the configuration file's info:
            - nzeros
            - modeId
            - modeAmp 
            - template
            - modalBase 
    bpath : str, OPTIONAL
        Base path of the file to read. Default points to the Configuration root
        folder
    """
    fname = os.path.join(bpath, iff_configFile)
    fnameBck = os.path.join(bpath, 'iffConfig_backup.ini')
    shutil.copyfile(fname, fnameBck)  # Create a backup of the original file
    config.read(fname)
    cc = config[key]
    cc[nzeroName]       = str(info['zeros'])
    cc[modeIdName]      = json.dumps(info['modes'].tolist())
    cc[modeAmpName]     = str(info['amplitude'])
    cc[modalBaseName]   = info['modalBase']
    cc[templateName]    = json.dumps(info['template'].tolist())
    with open(fname, 'w') as configfile:
        config.write(configfile)

def getNActs_fromConf(bpath=cfoldname):
    """
    Retrieves the number of actuators from the iffConfig.ini file. 
    DEPRECATED

    Parameters
    ----------
    bpath : str, OPTIONAL
        Base path of the file to read. Default points to the Configuration root\
        folder

    Returns
    -------
    nacts : int
        Number of DM's used actuators

    """
    fname = os.path.join(bpath, iff_configFile)
    config.read(fname)
    cc = config['DM']
    nacts = int(cc['NActs'])
    return nacts

def getTiming(bpath=cfoldname):
    """
    Retrieves the timing information from the iffConfig.ini file
    DEPRECATED??

    Parameters
    ----------
    bpath : str, OPTIONAL
        Base path of the file to read. Default points to the Configuration root\
        folder

    Returns
    -------
    timing : int
        Timing for the synchronization with the mirrors working frequency
    """
    fname = os.path.join(bpath, iff_configFile)
    config.read(fname)
    cc = config['DM']
    timing = int(cc['Timing'])
    return timing

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
    modeId      = np.array(json.loads(cc[modeIdName]))
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

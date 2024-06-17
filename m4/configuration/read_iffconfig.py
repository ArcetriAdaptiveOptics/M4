import os
import configparser
import json
import numpy as np
config=configparser.ConfigParser()
import m4.configuration.config_folder_names as fn
iff_configFile   = 'iffConfig.ini'
cfoldname = fn.CONFIGURATION_ROOT_FOLDER
nzeroName    = 'numberOfZeros'
modeIdName   = 'modeId'
modeAmpName  = 'modeAmp'
templateName = 'template'
modalBaseName= 'modalBase'

def getConfig(key, bpath=cfoldname):
    '''
    Reads the configuration file for the IFF acquisition. The key passed is the block of information retrieved

    Parameters
    ----------
    bpath : str
        Base path of the file to read
    key : str
        Key value of the block of information to read. Can be
            - 'TRIGGER'
            - 'REGISTRATION'
            - 'IFFUNC'

    Returns
    -------
    nzeros : int
        Number of zero columns that preceed the trigger mode
    modeId : ArrayLike
        Mode(s) to be applied
    modeAmp : float
        Amplitude of the applied mode(s)
    template : int | ArrayLike
        Template of the mode(s) to apply
    '''
    #if bpath is None:
    #    bpath = os.path.dirname(os.environ['PYOTTCONF'])
    #if os.path.isdir(bpath):
    #    fname = os.path.join(bpath, iff_configFile)
    #else:
    #    fname = bpath
    fname=os.path.join(cfoldname,iff_configFile)
    config.read(fname)
    cc = config[key]

    nzeros = int(cc[nzeroName])
    modeId = np.array(json.loads(cc[modeIdName]))
    modeAmp = float(cc[modeAmpName])
    modalBase = cc[modalBaseName]
    template = np.array(json.loads(cc[templateName]))

    return nzeros, modeId, modeAmp, template, modalBase

def getNActs_fromConf(bpath=cfoldname):
    """

    """
    #if bpath is None:
    #    bpath = os.path.dirname(os.environ['PYOTTCONF'])
    #if os.path.isdir(bpath):
    #    fname = os.path.join(bpath, iff_configFile)
    #else:
    #    fname = bpath
    fname=os.path.join(cfoldname,iff_configFile)
    config.read(fname)
    cc = config['DM']
    nacts = int(cc['NActs'])

    return nacts

def getTiming(bpath=cfoldname):
    """

    """
    #if bpath is None:
    #    bpath = os.path.dirname(os.environ['PYOTTCONF'])
    #if os.path.isdir(bpath):
    #    fname = os.path.join(bpath, iff_configFile)
    #else:
    #    fname = bpath
    #fname=os.path.join(cfoldname,iff_configFile)

    fname=os.path.join(cfoldname,iff_configFile)
    config.read(fname)
    cc = config['DM']
    timing = int(cc['Timing'])

    return timing

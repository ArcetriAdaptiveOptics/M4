import os
import configparser
import json
import numpy as np
config=configparser.ConfigParser()

bpath = os.environ['M4CONF']
iff_configFile   = 'iffConfig.ini'

nzeroName    = 'numberOfZeros'
modeIdName   = 'modeId'
modeAmpName  = 'modeAmp'
templateName = 'template'

def getConfig(key):
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
    fname = os.path.join(bpath, iff_configFile)
    print(fname)

    config.read(fname)
    cc = config[key]

    nzeros = int(cc[nzeroName])
    modeId = np.array(json.loads(cc[modeIdName]))
    modeAmp = float(cc[modeAmpName])
    template = np.array(json.loads(cc[templateName]))

    return nzeros, modeId, modeAmp, template

def getNActs_fromConf():
    fname = os.path.join(bpath, iff_configFile)
    config.read(fname)
    cc = config['DM']

    nacts = int(cc['NActs'])

    return nacts

def getTiming():
    fname = os.path.join(bpath, iff_configFile)
    config.read(fname)
    cc = config['DM']

    timing = int(cc['Timing'])

    return timing

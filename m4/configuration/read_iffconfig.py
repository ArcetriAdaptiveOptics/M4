import os
import configparser
import json
import numpy as np
config=configparser.ConfigParser()

iff_configFile   = 'iffConfig.ini'
triggerName      = 'TRIGGER'
registrationName = 'REGISTRATION'
iffuncName       = 'IFFUNC'

nzeroName    = 'numberOfZeros'
modeIdName   = 'modeId'
modeAmpName  = 'modeAmp'
templateName = 'template'

def getConfig(bpath,key):
    fname = os.path.join(bpath, iff_configFile)
    print(fname)
    config.read(fname)
    cc = config[key]
    nzeros = int(cc[nzeroName])
    modeId = np.array(json.loads(cc[modeIdName]))
    modeAmp = float(cc[modeAmpName])
    template = np.array(json.loads(cc[templateName]))
    return nzeros, modeId,modeAmp, template


import configparser
import numpy as np
config=configparser.ConfigParser()

hw4d_configfile = '4DSettings.ini'
cameraname = 'ACA2440'
namewidth = 'ImageWidthInPixels'
nameheigth    = 'ImageHeightInPixels'
nameoffx  = 'OffsetX'
nameoffy  = 'OffsetY'
namefreq  = 'FrameRate'

def getCameraConfig(path):
    config.read(path+'/'+hw4d_configfile)
    cc = config[cameraname]
    w = cc[namewidth]
    h = cc[nameheigth]
    offx = cc[nameoffx]
    offy = cc[nameoffy]
    freq = cc[namefreq]
    return np.array([w,h,offx, offy,freq])

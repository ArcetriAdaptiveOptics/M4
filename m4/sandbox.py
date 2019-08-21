'''
@author: cs
'''


import numpy as np
import os
from matplotlib import pyplot as plt


def main1908_createFileInfo():
    from m4.utils import createDevice 
    device= createDevice.myDevice("segment") 
    from m4.influenceFunctions import IFFunctions 
    IF= IFFunctions(device) 
    cmdMatrix=np.ones([800,700]) 
    cmdMatrix.T[30]=np.zeros(800) 
    indexing=np.array([10,11,12,13,14,15,16,17,18,19,30,31,32,33,34]) 
    amplitude= np.ones(15)
    save="/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions" 
    dove= IF.acq_IFFunctions(save,indexing,3, amplitude,cmdMatrix)
    #dove= IF.acq_IFFunctions(save,indexing,3,amplitude)
    
    return dove

def main1908_analyzer(tt):
    from m4.utils import createDevice 
    device= createDevice.myDevice("segment") 
    from m4.analyzerIFF import AnalyzerIFF
    #save="/Users/rm/Desktop/Arcetri/M4/ProvaCodice"
    save=os.path.join("/Users/rm/Desktop/Arcetri/M4/ProvaCodice/", tt)
    an= AnalyzerIFF.loadModalIFFInfoFromH5Folder(device,save)
    
    return an

def main2108():
    from m4.type.commandHistory import CmdHistory
    from m4.utils import createDevice 
    from m4.influenceFunctions import IFFunctions
    device= createDevice.myDevice("segment")
    IF= IFFunctions(device)
    cmdH= CmdHistory(device)
    
    indexing= np.array([11,12,13,14,15])
    cmdMatrix=np.ones([800,700]) 
    cmdMatrix.T[30]=np.zeros(800)
    amplitude=np.array([1,2,3,4,5])
    matrixToApply, indexingList= cmdH.createCmdHistory(
                                            indexing, 3, cmdMatrix)
    
    #vect= IF._amplitudeReorganization(indexing, indexingList, amplitude, 3) 
    
    return matrixToApply, indexingList
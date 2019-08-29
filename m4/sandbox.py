'''
@author: cs
'''


import numpy as np
import os
import copy
from matplotlib import pyplot as plt


def main1908_createFileInfo():
    from m4.utils import createDevice 
    device= createDevice.myDevice("segment") 
    from m4.influenceFunctionsMaker import IFFunctionsMaker 
    IF= IFFunctionsMaker(device) 
    cmdMatrix=np.ones([892,700]) 
    cmdMatrix.T[30]=np.zeros(892) 
    indexing=np.array([10,11,12,13,14,15,16,17,18,19,30,31,32,33,34]) 
    amplitude= np.ones(15)
    save="/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions" 
    #dove= IF.acq_IFFunctions(save,indexing,3, amplitude,cmdMatrix)
    dove= IF.acq_IFFunctions(save,indexing,3,amplitude)
    
    return dove

def main1908_analyzer(tt):
    from m4.analyzerIFFunctions import AnalyzerIFF
    save=os.path.join("/Users/rm/Desktop/Arcetri/M4/ProvaCodice/", tt)
    #an= AnalyzerIFF.loadModalIFFInfoFromH5Folder(device,save)
    an= AnalyzerIFF.loadZonalIFFInfoFromH5Folder(save)
    
    return an

def main2108_amplitudeReaorgTEST():
    from m4.type.commandHistory import CmdHistory
    from m4.utils import createDevice 
    from m4.influenceFunctionsMaker import IFFunctionsMaker
    device= createDevice.myDevice("segment")
    IF= IFFunctionsMaker(device)
    cmdH= CmdHistory(device)
    
    indexing= np.array([11,12,13,14,15])
    cmdMatrix=np.ones([892,700]) 
    cmdMatrix.T[30]=np.zeros(892)
    amplitude=np.array([1,2,3,4,5])
    indexingImput= copy.copy(indexing)
    matrixToApply, indexingList= cmdH.createShuffleCmdHistory(
                                            indexing, 3, cmdMatrix)
    
    vect= IF._amplitudeReorganization(indexingImput, indexingList, amplitude, 3) 
    
    return indexingList, vect


def main2808_commandHistoryTest():
    from m4.utils import createDevice 
    from m4.type.commandHistory import CmdHistory 
    device= createDevice.myDevice("segment") 
    cmdH= CmdHistory(device)
    modeVector=np.array([11,12,13,14,15])
    cmdMatrix=np.ones([892,700]) 
    cmdMatrix.T[30]=np.zeros(892)
    amplitude=np.array([1,2,3,4,5])
    matrixToApply= cmdH.shuffleCommandHistoryMaker(modeVector, 
                                                   amplitude, cmdMatrix, 3)
    #matrixToApply2= cmdH.tidyCommandHistoryMaker(modeVector, 
                                                    #amplitude, cmdMatrix, 3)
    
    return matrixToApply

def main2908_provaIFF():
    from m4.utils import createDevice 
    device= createDevice.myDevice("segment") 
    from m4.influenceFunctionsMaker import IFFunctionsMaker
    IF= IFFunctionsMaker(device) 
    
    #cmdMatrix=np.ones([892,800]) 
    #cmdMatrix.T[30]=np.zeros(892) 
    cmdMatrix=np.zeros((892,800))  
    for i in range(800):
        cmdMatrix[i][i]=1
    modeVect=np.array([10,11,12,13,14,15,16,17,18,19,30,31,32,33,34]) 
    amp= np.array([0,1,2,3,4,5,6,7,8,9,0.30,0.31,0.32,0.33,0.34])
    
    IF.acq_IFFunctions(modeVect, 3, amp, cmdMatrix, 1)
    
    
    
    
    
    
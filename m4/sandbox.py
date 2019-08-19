'''
@author: cs
'''


import numpy as np
from matplotlib import pyplot as plt


def main1908_createFileInfo():
    from m4.utils import createDevice 
    device= createDevice.myDevice("segment") 
    from m4.influenceFunctions import IFFunctions 
    IF= IFFunctions(device) 
    cmdMatrix=np.ones([800,700]) 
    cmdMatrix.T[30]=np.zeros(800) 
    indexing=np.array([10,11,12,13,14,15,16,17,18,19,30,31,32,33,34]) 
    save="/Users/rm/Desktop/Arcetri/M4/ProvaCodice" 
    mat, ind= IF.pushPull(save,indexing,3,0.1,cmdMatrix)
    #mat, ind= IF.pushPull(save,indexing,3,0.1)
    
    return mat, ind

def main1908_analyzer():
    from m4.utils import createDevice 
    device= createDevice.myDevice("segment") 
    from m4.analyzerIFF import AnalyzerIFF
    save="/Users/rm/Desktop/Arcetri/M4/ProvaCodice"
    an= AnalyzerIFF.fromH5Folder(device,save)
    
    return an
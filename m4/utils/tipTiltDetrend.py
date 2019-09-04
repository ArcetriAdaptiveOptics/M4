'''
@author: cs
'''

import numpy as np
from m4.utils.configuration import Configuration


class TipTiltDetrend():
    
    def __init__(self):
        self._pupillXYR= Configuration.ParabolaPupilXYRadius
    
    
    def tipTiltRemover(self, image, roi, zm4):
        imaList=[]
        for i in range(roi.shape[1]):
            a= np.ma.masked_array(image.data, mask= roi[:,i])
            imaList.append(a)
        
        zm4.setPupilCenterAndRadiusInIFCoords(self._pupillXYR[0], self._pupillXYR[1], self._pupillXYR[2])    
            
        tipTiltList=[]
        for ima in imaList:
            tt= zm4.zernikeFit(ima, np.array([2,3]))
            tipTiltList.append(tt)
        
        ttMean= np.mean(tt, 0)
        zernikeCmd, surfaceMap= zm4.zernikeToDMCommand(ttMean)
        ttImage= image - surfaceMap
            
        return ttImage
        
'''
@author: cs
'''

import numpy as np
from m4.utils.roi import ROI


class TipTiltDetrend():
    
    
    def tipTiltRemover(self, image, roi, zm4):
        imaList=[]
        roiCoord=[]
        for i in range(roi.shape[1]):
            a= np.ma.masked_array(image.data, mask= roi[:,i])
            pupillXYR= ROI.getZernikeCoord(roi[i])
            roiCoord.append(pupillXYR)
            imaList.append(a)
            
        tipTiltList=[]
        for ima, pupillXYR in imaList, roiCoord:
            zm4.setPupilCenterAndRadiusInIFCoords(pupillXYR[0], pupillXYR[1], pupillXYR[2])
            tt= zm4.zernikeFit(ima, np.array([2,3]))
            tipTiltList.append(tt)
        
        ttMean= np.mean(tt, 0)
        
        cmdList=[]
        surfaceMapList=[]
        for pupillXYR in roiCoord:
            zm4.setPupilCenterAndRadiusInIFCoords(pupillXYR[0], pupillXYR[1], pupillXYR[2])
            zernikeCmd, surfaceMap= zm4.zernikeToDMCommand(ttMean)
            cmdList.append(zernikeCmd)
            surfaceMapList.append(surfaceMap)
        
        ttImage=image
        for sur in surfaceMapList:
            ttImage= ttImage + sur
            
        return ttImage
        
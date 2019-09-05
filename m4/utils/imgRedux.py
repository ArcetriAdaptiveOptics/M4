'''
@author: cs
'''

import numpy as np
from m4.utils.configuration import Configuration
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.utils.roi import ROI


class TipTiltDetrend():
    
    def __init__(self):
        self._pupillXYR= Configuration.ParabolaPupilXYRadius
    
    
    def tipTiltRemover(self, image, roi, analyzerIFFunctions):
        imaList=[]
        for r in roi:
            a= np.ma.masked_array(image.data, mask= r)
            imaList.append(a)
        
        zm4= ZernikeGenerator(analyzerIFFunctions)
        zm4.setPupilCenterAndRadiusInIFCoords(self._pupillXYR[0], self._pupillXYR[1], self._pupillXYR[2])    
            
        tipTiltList=[]
        for ima in imaList:
            tt= zm4.zernikeFit(ima, np.array([2,3]))
            tipTiltList.append(tt)
        
        ttMean= np.mean(tt, 0)
        zernikeCmd, surfaceMap= zm4.zernikeToDMCommand(ttMean)
        ttImage= image - surfaceMap
            
        return ttImage
        
        

class PhaseSolve():
    
    def __init__(self):
        self._r=ROI()
        self._lambda= Configuration.Lambda
        self._n= None
    
    
    def n_calculator(self, splValues): 
        n=np.array(splValues.shape[0])   
        for i in range(splValues.shape[0]):
            n[i]= (2.* splValues[i]) / self._lambda
        self._n= n
        return self._n
    
    
    def m4PhaseSolver(self, m4Ima): 
        roiList= self._r._ROIonM4(m4Ima)
        
        media=[]
        imgList=[]
        for roi in roiList:
            img= np.ma.masked_array(m4Ima, mask=roi)
            m= img.mean()
            media.append(m)
            imgList.append(img)
               
        aa= np.arange(self._n.shape[0])   
        zipped= zip(aa, imgList)
        img_phaseSolveList=[]
        for i, img in zipped:
            img_phaseSolve= img - self._n[i]
            img_phaseSolveList.append(img_phaseSolve)
        
        img_phaseSolveList[0]= imgList[0]
        
        for j in range(img_phaseSolveList.shape[0]):
            m4NewImage= img_phaseSolveList[0] + img_phaseSolveList[j]
            
        return m4NewImage
        
        
    def masterRoiPhaseSolver(self, segIma):
        pass
            
            
            
            
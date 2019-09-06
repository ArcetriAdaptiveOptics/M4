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
        n=np.zeros(splValues.shape[0])   
        for i in range(splValues.shape[0]):
            n[i]= (2.* splValues[i]) / self._lambda
        self._n= n
        return self._n
    
    
    def m4PhaseSolver(self, m4Ima, splValues): 
        self.n_calculator(splValues)
        roiList= self._r._ROIonM4(m4Ima)
        m4NewImage= None
        
        media=[]
        imgList=[]
        for roi in roiList:
            img= np.zeros(m4Ima.shape)
            img[np.where(roi== True)]= np.ma.compress(roi.ravel(), m4Ima)
            imgg= np.ma.masked_array(img, mask= roi)
            m= img.mean()
            media.append(m)
            imgList.append(imgg)
               
        aa= np.arange(self._n.shape[0])   
        zipped= zip(aa, imgList)
        img_phaseSolveList=[]
        for i, imgg in zipped:
            img_phaseSolve= np.ma.masked_array(imgg.data - self._n[i], mask= np.invert(imgg.mask))
            img_phaseSolveList.append(img_phaseSolve)
        
        #img_phaseSolveList[len(img_phaseSolveList)-1].data= imgList[len(imgList)-1].data
          
          
        for j in range(1, len(img_phaseSolveList)):
            if m4NewImage is None:
                m4NewImage= np.ma.array(img_phaseSolveList[0].filled(1)* img_phaseSolveList[j].filled(1), 
                                         mask=(img_phaseSolveList[0].mask * img_phaseSolveList[j].mask))
            else:
                m4NewImage = np.ma.array(m4NewImage.filled(1) * img_phaseSolveList[j].filled(1), 
                                         mask=(m4NewImage.mask * img_phaseSolveList[j].mask))
            
        return m4NewImage, img_phaseSolveList, imgList
    
    
    
    
    
        
        
    def masterRoiPhaseSolver(self, segIma):
        pass
            
            
            
            
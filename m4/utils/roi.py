'''
@author: cs
'''

import numpy as np


class ROI():
    
    def __init__(self):
        self._DetectorROI= None
        self._AnalysisROI= None
        self._pupilXYRadius= None
        
    
    def getDetectorROI(self, image):
        self._DetectorROI= image.mask()
        return self._DetectorROI
    
    def getAnalysisROI(self):
        return self._AnalysisROI 
        
    def setAnalysisROI(self, analysisMask):
        self._AnalysisROI= analysisMask
        
        
    def getZernikeCoord(self, roi):
        return self._pupilXYRadius 
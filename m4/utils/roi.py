'''
@author: cs
'''

import numpy as np
import sklearn.feature_extraction as skf_e
import sklearn.cluster as skc
import skimage.segmentation as sks


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
    
    def _ROIonSegment(self, ima):
        graph = skf_e.image.img_to_graph(ima.data, mask=ima.mask)
        labels = skc.spectral_clustering(graph, n_clusters=4, eigen_solver='arpack')
        label_im = -np.ones(ima.mask.shape)
        label_im[ima.mask] = labels
        
        markers = ima.mask.astype('int')*0
        markers[0,0]=1
        markers[180,79]=2
        markers[296,258]=3 
        markers[170,436]=4
        markers[34,258]=5
        roi_mask = sks.random_walker(ima.mask, markers) 
        
        roiList=[]
        for i in range(2,5):
            maski=np.zeros(roi_mask.shape, dtype=np.bool)
            maski[np.where(roi_mask==i)]=1 
            roiList.append(maski)
            
        return roiList
    
'''
@author: cs
'''

import numpy as np


class TipTiltDetrend():
    
    
    def tipTiltRemover(self, image, roi):
        imaList=[]
        for i in range(roi.shape[1]):
            a= np.ma.masked_array(image.data, mask= roi[:,i])
            imaList.append(a)
            
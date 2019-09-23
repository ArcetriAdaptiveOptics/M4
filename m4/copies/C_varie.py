'''
@author: cs
'''


# TEST
        folder= os.path.join(Configuration.CALIBRATION_ROOT_FOLDER, "IFFunctions", tt)
        who, tt_cmdH, actsVector, cmdMatrix, cmdAmpl, nPushPull, indexingList= IF.loadInfo(folder)
        from m4.type.commandHistory import CmdHistory
        cmdH= CmdHistory.load(device, tt)
        commandHistory= cmdH._cmdHToApply
        
        
        
        
# ROI

        def tipTiltRemover(self, image, roi, finalIndex, analysisInd):   
        ''' 
            arg:
                image= immagine da analizzare
                roi= roi dell'immagine
                finalIndex= indice della roi finale
                analysisInd= indice delle roi da utilizzare per l'analisi
        '''
        coefList=[]
        for r in roi:
            ima= np.ma.masked_array(image.data, mask=r)
            coef= self._zernikeFit(ima, np.array([2,3]))
            coefList.append(coef)
            
        tip= (coefList[0][0] + coefList[2][0])/2.
        tilt= (coefList[0][1] + coefList[2][1])/2.
        
        surfcoef= np.array([tip, tilt]) 
        surfaceMap=self._zernikeSurface(surfcoef)
        
        cx= self._pupilXYRadius[0]
        cy= self._pupilXYRadius[1]
        r= self._pupilXYRadius[2]
        imaCut=image[cy-r:cy+r, cx-r:cx+r]
        imageTTR= np.ma.masked_array(imaCut.data - surfaceMap, mask=roi[1])
            
        return surfaceMap, imageTTR
    
            
    def centralViewTipTilRemover(self, m4image, m4roi):
        a= np.arange(len(m4roi))-1
        a[0]=5
        b= np.arange(len(m4roi))
        c= np.arange(len(m4roi))+1
        c[5]=0
        
        zipped= zip(a, b, c) 
        finalImaList=[]
        for  a,b,c in zipped:
            roiList=[]
            roiList.append(m4roi[a])
            roiList.append(m4roi[b])
            roiList.append(m4roi[c])
            
            surfaceMap, imageTTR= self.tipTiltRemover(m4image, roiList)
            finalImaList.append(imageTTR)
        
        return finalImaList
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
    
    
    
    
    
    
# ZERNIKE

        self._an = analyzerIFFunctions
        self._shapeIFs= self._an.getCube()[:,:,0].shape
        self._nPixelOnDiameter= None
        self._nActs= self._an.getCube()[:,0,0].shape
        self._pupilXYRadius= None
        self._zg= None 
        
    def _zernikeFit(self, img, zernikeMode):
        '''
        zernikeMode= vector of Zernike modes to remove
        '''
        mat= np.zeros((img.compressed().shape[0], zernikeMode.size)) 
        for i in range(0, zernikeMode.size):
            z=self._zg.getZernike(zernikeMode[i])
            aa= np.ma.masked_array(z, mask= img.mask)
            mat.T[i]= aa.compressed()
            
        self._mat= mat
        inv= np.linalg.pinv(mat)   
        a= np.dot(inv, img.compressed())  
        return a
    
    
    def _zernikeSurface(self, surfaceZernikeCoeffArray, roi, ind):
        zernikeSurfaceMap= np.dot(self._totalMatList[ind], surfaceZernikeCoeffArray)
        mask= np.invert(roi[ind])
        surf= np.ma.masked_array(np.zeros((2*self._pupilXYRadius[2],2*self._pupilXYRadius[2])), mask=mask)               
        surf[mask]= zernikeSurfaceMap 
        return surf
  
            
            
            
# PHASE SOLVE

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
            img_phaseSolve= np.ma.masked_array(imgg.data - self._n[i], mask= imgg.mask)
            img_phaseSolveList.append(img_phaseSolve)
        
        img_phaseSolveList[len(img_phaseSolveList)-1]= np.ma.masked_array(imgList[len(imgList)-1].data, 
                                                                mask= np.invert(imgList[len(imgList)-1].mask))
          
          
        for j in range(1, len(img_phaseSolveList)):
            if m4NewImage is None:
                m4NewImage= np.ma.array(img_phaseSolveList[0].filled(1)* img_phaseSolveList[j].filled(1), 
                                         mask=(img_phaseSolveList[0].mask * img_phaseSolveList[j].mask))
            else:
                m4NewImage = np.ma.array(m4NewImage.filled(1) * img_phaseSolveList[j].filled(1), 
                                         mask=(m4NewImage.mask * img_phaseSolveList[j].mask))
            
        return m4NewImage, img_phaseSolveList, imgList
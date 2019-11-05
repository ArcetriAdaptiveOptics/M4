'''
@author: cs
'''


import numpy as np
import os
import copy
import pyfits
from matplotlib import pyplot as plt
from m4.utils import roi
from m4.ground.configuration import Configuration
from m4.ground import objectFromFitsFileName



def main1908_createFileInfo():
    from m4.utils import createDevice 
    device= createDevice.myDevice("segment") 
    from m4.influenceFunctionsMaker import IFFunctionsMaker 
    IF= IFFunctionsMaker(device) 
    cmdMatrix=np.ones([892,700]) 
    cmdMatrix.T[30]=np.zeros(892) 
    indexing=np.array([10,11,12,13,14,15,16,17,18,19,30,31,32,33,34]) 
    amplitude= np.ones(15)
    save="/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions" 
    #dove= IF.acq_IFFunctions(save,indexing,3, amplitude,cmdMatrix)
    dove= IF.acq_IFFunctions(save,indexing,3,amplitude)
    
    return dove

def main1908_analyzer(tt):
    from m4.analyzerIFFunctions import AnalyzerIFF
    save=os.path.join("/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions", tt)
    an= AnalyzerIFF.loadInfoFromh5Folder(save)
    
    return an

def main2108_amplitudeReaorgTEST():
    from m4.type.commandHistory import CmdHistory
    from m4.utils import createDevice 
    from m4.influenceFunctionsMaker import IFFunctionsMaker
    device= createDevice.myDevice("segment")
    IF= IFFunctionsMaker(device)
    cmdH= CmdHistory(device)
    
    indexing= np.array([11,12,13,14,15])
    cmdMatrix=np.ones([892,700]) 
    cmdMatrix.T[30]=np.zeros(892)
    amplitude=np.array([1,2,3,4,5])
    indexingImput= copy.copy(indexing)
    matrixToApply, indexingList= cmdH.createShuffleCmdHistory(
                                            indexing, 3, cmdMatrix)
    
    vect= IF._amplitudeReorganization(indexingImput, indexingList, amplitude, 3) 
    
    return indexingList, vect


def main2808_commandHistoryTest():
    from m4.utils import createDevice 
    from m4.type.commandHistory import CmdHistory 
    device= createDevice.myDevice("segment") 
    cmdH= CmdHistory(device)
    modeVector=np.array([11,12,13,14,15])
    cmdMatrix=np.ones([892,700]) 
    cmdMatrix.T[30]=np.zeros(892)
    amplitude=np.array([1,2,3,4,5])
    matrixToApply= cmdH.shuffleCommandHistoryMaker(modeVector, 
                                                   amplitude, cmdMatrix, 3)
    #matrixToApply2= cmdH.tidyCommandHistoryMaker(modeVector, 
                                                    #amplitude, cmdMatrix, 3)
    
    return matrixToApply

def main2908_provaIFF():
    from m4.utils import createDevice 
    device= createDevice.myDevice("segment") 
    from m4.influenceFunctionsMaker import IFFunctionsMaker
    IF= IFFunctionsMaker(device) 
    
    #cmdMatrix=np.ones([892,800]) 
    #cmdMatrix.T[30]=np.zeros(892) 
    cmdMatrix=np.zeros((892,800))  
    for i in range(800):
        cmdMatrix[i][i]=1
    modeVect=np.array([10,11,12,13,14,15,16,17,18,19,30,31,32,33,34]) 
    amp= np.array([0,1,2,3,4,5,6,7,8,9,0.30,0.31,0.32,0.33,0.34])
    
    IF.acq_IFFunctions(modeVect, 3, amp, cmdMatrix, 1)
    
    
def main0409_imageM4():
    fitsRoot= "/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova"
    fitsFileName= os.path.join(fitsRoot, 'mode_0005.fits')
    hduList= pyfits.open(fitsFileName)
    ima= hduList[0].data
    m4= ima[1]
    fitsFileName= os.path.join(fitsRoot, 'mode_0006.fits')
    hduList= pyfits.open(fitsFileName)
    ima= hduList[0].data
    segment= ima[1]
    
    #mask=np.zeros(segment.shape, dtype=np.bool) 
    #mask[np.where(segment==segment.max())]=1  
    imaseg=np.ma.masked_array(segment * 632.8e-9, mask=segment)
    
    mask=np.zeros(m4.shape, dtype=np.bool) 
    mask[np.where(m4==m4.max())]=1  
    imam4=np.ma.masked_array(m4 * 632.8e-9, mask=mask)
    return m4, segment, imaseg, imam4


def main0409_ROI(ima):
    import sklearn.feature_extraction as skf_e
    import sklearn.cluster as skc
    import skimage.segmentation as sks
    graph = skf_e.image.img_to_graph(ima.data, mask=ima.mask)
    labels = skc.spectral_clustering(graph, n_clusters=7, eigen_solver='arpack')
    label_im = -np.ones(ima.mask.shape)
    label_im[ima.mask] = labels
    
    markers = ima.mask.astype('int')*0
    markers[0,0]=1
    markers[83,257]=2
    markers[164,407]=3 
    markers[339,409]=4
    markers[433,254]=5
    markers[338,103]=6
    markers[157,101]=7
    markers[232,270]=8
    roi_mask = sks.random_walker(ima.mask, markers)  
    
    return label_im, markers, roi_mask
    
  
def main0509_makeM4ImgWhitPhase(m4):
    mask=np.zeros(m4.shape, dtype=np.bool) 
    mask[np.where(m4==m4.max())]=1  
    imam4=np.ma.masked_array(m4 * 10, mask=mask)
    
    from m4.utils.roi import ROI
    r=ROI()
    roi= r._ROIonM4(imam4)

    finalIma= None
    a=np.array((10,20,30,40,50,60))
    aa= np.arange(a.shape[0])
    cc= np.arange(len(roi)-1)
    zipped= zip(aa, cc)
    
    bbList=[]
    for i,j in zipped:
        bb=np.zeros(imam4.shape)
        bb[np.where(roi[j]== True)]=a[i]
        bbList.append(bb)
    
    for bb in bbList:   
        if finalIma is None: 
            finalIma= imam4.data + bb.data
        else:
            finalIma= finalIma + bb.data
    
    finalImag= np.ma.masked_array(finalIma, mask=np.invert(imam4.mask))
    
    fitsFileName= os.path.join('/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova', 'imgProvaM4.fits')
    header= pyfits.Header()
    pyfits.writeto(fitsFileName, finalImag.data, header)
    pyfits.append(fitsFileName, finalImag.mask.astype(int), header)
    pyfits.append(fitsFileName, a, header)
    return finalImag

def main1009_makeSegImgWhitPhase(seg):
    mask=np.zeros(seg.shape, dtype=np.bool) 
    mask[np.where(seg==seg.max())]=1  
    imaseg=np.ma.masked_array(seg * 10, mask=mask)
    
    from m4.utils.roi import ROI
    r=ROI()
    roi= r._ROIonSegment(imaseg)
    finalIma= None
    a=np.array((0,20,40))
    aa= np.arange(a.shape[0])
    cc= np.arange(len(roi))
    zipped= zip(aa, cc)
    bbList=[]
    for i,j in zipped:
        bb=np.zeros(imaseg.shape)
        bb[np.where(roi[j]== True)]=a[i]
        bbList.append(bb)
    for bb in bbList:   
        if finalIma is None: 
            finalIma= imaseg.data + bb.data
        else:
            finalIma= finalIma + bb.data
    
    finalImag= np.ma.masked_array(finalIma, mask=imaseg.mask)
    fitsFileName= os.path.join('/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova', 'imgProvaSeg.fits')
    header= pyfits.Header()
    pyfits.writeto(fitsFileName, finalImag.data, header)
    pyfits.append(fitsFileName, finalImag.mask.astype(int), header)
    pyfits.append(fitsFileName, a, header)
    return finalImag
    
    
    
def main0609_TESTpsM4():  
    fitsFileName= os.path.join('/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova', 'imgProva.fits')
    header= pyfits.getheader(fitsFileName)
    hduList= pyfits.open(fitsFileName)
    finalIma= np.ma.masked_array(hduList[0].data, hduList[1].data.astype(bool))
    #m4, seg, imas, ima= main0409_imageM4()
    #finalIma= main0509_makeImgWhitPhase(m4) 
    from m4.utils.imgRedux import PhaseSolve
    ps= PhaseSolve()
    spl=np.array((20,30,40,50,60,10))
    m4NewImage, img_phaseSolveList, imgList= ps.m4PhaseSolver(finalIma, spl)
    return m4NewImage, img_phaseSolveList, imgList
    

def main1009_TESTpsSeg():
    fitsFileName= os.path.join('/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova', 'imgProvaSeg.fits')
    header= pyfits.getheader(fitsFileName)
    hduList= pyfits.open(fitsFileName)
    finalIma= np.ma.masked_array(hduList[0].data, hduList[1].data.astype(bool))
    from m4.utils.imgRedux import PhaseSolve
    ps= PhaseSolve()
    spl=np.array([5])
    img_phaseSolve= ps.masterRoiPhaseSolver(finalIma, spl)
    return img_phaseSolve


def main1009_tiptiltImage():
    fitsFileName= os.path.join('/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova', 'imgProvaSeg.fits')
    header= pyfits.getheader(fitsFileName)
    hduList= pyfits.open(fitsFileName)
    finalIma= np.ma.masked_array(hduList[0].data, hduList[1].data.astype(bool))
    from m4.utils.roi import ROI
    r=ROI()
    roi= r._ROIonSegment(finalIma)
    from m4.utils.imgRedux import TipTiltDetrend 
    tt= TipTiltDetrend()
    surfaceMap, imageTTR= tt.tipTiltRemover(finalIma, roi[3])
    return finalIma, imageTTR

def main1109_tiptiltZernike(imas):
    from m4.utils.imgRedux import TipTiltDetrend 
    tt= TipTiltDetrend()
    sur= tt._zernikeSurface(np.array([10,10, 3]))
    aa= np.ma.masked_array(sur, mask= np.invert(imas.mask))
    from m4.utils.roi import ROI
    r=ROI()
    roi= r._ROIonSegment(imas)
    
    surfaceMap, imageTTR= tt.tipTiltRemover(aa, roi, 3)
    return aa, imageTTR


### FUNZIONI PER TEST IFF ###    
def testIFF_shuffleMeasureCreator(device, cmdMatrixTag, modeVectTag, ampTag, nPushPull):
    from m4.type.modesVector import ModesVector
    mv= ModesVector.loadFromFits(modeVectTag)
    modeVectInput= mv.getModesVector()
    from m4.influenceFunctionsMaker import IFFunctionsMaker
    IF= IFFunctionsMaker(device) 
    
    tt= IF.acq_IFFunctions(modeVectTag, nPushPull, ampTag, cmdMatrixTag, 1)
    
    folder= os.path.join(Configuration.CALIBRATION_ROOT_FOLDER, "IFFunctions", tt)
    who, tt_cmdH, actsVector, cmdMatrix, amplitude, nPushPull, indexingList= IF.loadInfoFromFits(folder)
    
    cube= IF._testIFFunctions_createCube25fromFileFitsMeasure()
    from m4.type.commandHistory import CmdHistory
    cmdH= CmdHistory(device)
    amplReorg= cmdH._amplitudeReorganization(modeVectInput, indexingList, amplitude, nPushPull)
     
    misure= None 
    for i in range(indexingList.shape[0]):
        for j in range(indexingList.shape[1]):
            mask= np.invert(cube[:,:,indexingList[i][j]].mask)
            k= i * indexingList.shape[1] + j
            if misure is None:
                misure= np.ma.masked_array(cube[:,:,indexingList[i][j]].data * amplReorg[k], mask= mask)
                misure= np.ma.dstack((misure, np.ma.masked_array(cube[:,:,indexingList[i][j]].data * amplReorg[k] * -1, 
                                                                 mask= mask)))
            else:
                misure= np.ma.dstack((misure, np.ma.masked_array(cube[:,:,indexingList[i][j]].data * amplReorg[k] * 1,
                                                                  mask= mask)))
                misure= np.ma.dstack((misure, np.ma.masked_array(cube[:,:,indexingList[i][j]].data * amplReorg[k] * -1,
                                                                  mask= mask)))
    
    fitsFileName= os.path.join(folder, 'misure.fits')
    header= pyfits.Header()
    header['NPUSHPUL']= nPushPull
    header['WHO']= who
    header['TT_CMDH']= tt_cmdH
    pyfits.writeto(fitsFileName, actsVector, header)
    pyfits.append(fitsFileName, cmdMatrix, header)
    pyfits.append(fitsFileName, amplitude, header)
    pyfits.append(fitsFileName, indexingList, header)
    pyfits.append(fitsFileName, misure.data, header)
    pyfits.append(fitsFileName, misure.mask.astype(int), header)
    return tt

def testIFF_tidyMeasureCreator(device, cmdMatrixTag, modeVectTag, ampTag, nPushPull):
    from m4.influenceFunctionsMaker import IFFunctionsMaker
    IF= IFFunctionsMaker(device) 
    
    tt= IF.acq_IFFunctions(modeVectTag, nPushPull, ampTag, cmdMatrixTag)
    
    folder= os.path.join(Configuration.CALIBRATION_ROOT_FOLDER, "IFFunctions", tt)
    who, tt_cmdH, actsVector, cmdMatrix, amplitude, nPushPull, indexingList= IF.loadInfoFromFits(folder)
    
    cube= IF._testIFFunctions_createCube25fromFileFitsMeasure()
    ampl= np.tile(amplitude, nPushPull)
    
    misure= None 
    for i in range(indexingList.shape[0]):
        for j in range(indexingList.shape[1]):
            mask= np.invert(cube[:,:,indexingList[i][j]].mask)
            k= i * indexingList.shape[1] + j
            if misure is None:
                misure= np.ma.masked_array(cube[:,:,indexingList[i][j]].data * ampl[k], mask= mask)
                misure= np.ma.dstack((misure, np.ma.masked_array(cube[:,:,indexingList[i][j]].data * ampl[k] * -1, 
                                                                 mask= mask)))
            else:
                misure= np.ma.dstack((misure, np.ma.masked_array(cube[:,:,indexingList[i][j]].data * ampl[k] * 1,
                                                                  mask= mask)))
                misure= np.ma.dstack((misure, np.ma.masked_array(cube[:,:,indexingList[i][j]].data * ampl[k] * -1,
                                                                  mask= mask)))
    
    fitsFileName= os.path.join(folder, 'misure.fits')
    header= pyfits.Header()
    header['NPUSHPUL']= nPushPull
    header['WHO']= who
    header['TT_CMDH']= tt_cmdH
    pyfits.writeto(fitsFileName, actsVector, header)
    pyfits.append(fitsFileName, cmdMatrix, header)
    pyfits.append(fitsFileName, amplitude, header)
    pyfits.append(fitsFileName, indexingList, header)
    pyfits.append(fitsFileName, misure.data, header)
    pyfits.append(fitsFileName, misure.mask.astype(int), header)    
    return tt
    

def testIFF_an(tt, ttD= None):
    from m4.analyzerIFFunctions import AnalyzerIFF
    fileName= os.path.join("/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions", tt)
    an= AnalyzerIFF.loadTestMeasureFromFits(fileName)
    cube= an.createCube(ttD)
    intMat= an.getInteractionMatrix()
    rec= an.getReconstructor()
    prod= np.dot(rec, intMat)
    return an, prod, cube

def testIFF_spiano(an):
    ampr= np.random.randn(25)
    wf= np.dot(an._cube, ampr)
    
    an.setDetectorMask(wf.mask | an.getMasterMask())
    rec= an.getReconstructor()
    wf_masked = np.ma.masked_array(wf.data, 
                                       mask=np.ma.mask_or(wf.mask, an.getMasterMask()))
        
    amp= np.dot(rec, wf_masked.compressed())
    surf= np.dot(an._cube, amp)
    spWf= wf - surf
    return amp, spWf
    
### FINE FUNZIONI TEST IFF ###

def main1709_ttDetrend(image):
    from m4.utils.roi import ROI
    r=ROI()
    roi= r._ROIonSegment(image)
    from m4.utils.imgRedux import TipTiltDetrend 
    tt= TipTiltDetrend()
    surfaceMap, imageTTR= tt.tipTiltRemover(image, roi, 3)
    return surfaceMap, imageTTR

def main2609_ttDetrend(image):
    from m4.utils.roi import ROI
    r=ROI()
    roi= r._ROIonSegment(image)
    from m4.utils.imgRedux import TipTiltDetrend 
    tt= TipTiltDetrend()
    return roi, tt

def immaginiprova():
    fitsRoot= "/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova"
    fitsFileName= os.path.join(fitsRoot, 'mode_0005.fits')
    hduList= pyfits.open(fitsFileName)
    ima= hduList[0].data
    m4= np.ma.masked_array(ima[0], mask= np.invert(ima[1]))
    fitsFileName= os.path.join(fitsRoot, 'mode_0006.fits')
    hduList= pyfits.open(fitsFileName)
    ima= hduList[0].data
    segment= np.ma.masked_array(ima[0], mask= np.invert(ima[1]))
    return m4, segment

def immaginiProvaTTDetrendSeg():
    push= objectFromFitsFileName.readImageFromFitsFileName('Seg/img_0000.fits')
    pull= objectFromFitsFileName.readImageFromFitsFileName('Seg/img_0001.fits')
    mode0= np.ma.masked_array(pull.data - push.data, mask= push.mask)
    
#     push= objectFromFitsFileName.readImageFromFitsFileName('Seg/img_0002.fits')
#     pull= objectFromFitsFileName.readImageFromFitsFileName('Seg/img_0003.fits')
#     mode1= np.ma.masked_array(pull.data - push.data, mask= np.invert(push.mask))
    return mode0
    
def immaginiProvaTTDetrendAll():
    push= objectFromFitsFileName.readImageFromFitsFileName('All/img_0000.fits')
    pull= objectFromFitsFileName.readImageFromFitsFileName('All/img_0001.fits')
    mode0= np.ma.masked_array(pull.data - push.data, mask= push.mask)
    
    push= objectFromFitsFileName.readImageFromFitsFileName('All/img_0002.fits')
    pull= objectFromFitsFileName.readImageFromFitsFileName('All/img_0003.fits')
    mode1= np.ma.masked_array(pull.data - push.data, mask= push.mask)
    return mode0, mode1


def main1001_calib():
    ampVect= np.ones(5)*5.0e-06 
    from m4.utils.opticalCalibration import Opt_Calibration
    a= Opt_Calibration()
    a._commandAmpVector= ampVect 
    
    a.createCube()
    cub= a.getCube() 
    ima= cub[:,:,0]
    from m4.utils.roi import ROI 
    r= ROI()
    roi= r._ROIonAlignmentImage(ima)  
    mask= roi[2] 
    mat= a.getInteractionMatrix(mask)
    rec= a.getReconstructor(mask)
    return mat, rec






    
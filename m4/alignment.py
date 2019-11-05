'''
@author: cs
'''


from m4.utils.opticalAlignment import Opt_Alignment
from m4.utils.opticalCalibration import Opt_Calibration
from m4.ground import objectFromFitsFileName as obj
from m4.utils.roi import ROI
from m4.utils.zernikeOnM4 import ZernikeOnM4
import numpy as np

class Alignment():
    
    def __init__(self):
        self._cal= Opt_Calibration()
        self._r= ROI()
        self._zOnM4= ZernikeOnM4()
        
        
    def OTT_Calibration(self, commandAmpVector, nPushPull, maskIndex):
        '''
            arg:
                commandAmpVector= vettore contenente i valori dei movimenti dei 5 gradi di libertà
                nPushPull= numero di push pull per ogni grado di libertà
                maskIndex= int (3 per la maschera dell'RM)
        '''
        self._moveRM()
        self._tt= self._cal.measureCalibrationMatrix(0, commandAmpVector, nPushPull)
        intMat, rec= self._cal.analyzerCalibrationMeasurement(self._tt, maskIndex)
        return self._tt
    
    def OTT_Alignement(self, tt= None):
        if tt is None:
            a= Opt_Alignment(self._tt)
        else:
            a= Opt_Alignment(tt)
        cmd= a.opt_align()
        self._applyCmd()
        return cmd
       
    
    def M4_Calibration(self, commandAmpVector_ForPARRMAlignement, nPushPull_ForPARRMAlignement, 
                       maskIndex_ForPARRMAlignement, commandAmpVector_ForM4Calibration, 
                       nPushPull_ForM4Calibration, maskIndex_ForM4Alignement):
        self._moveSegmentView()
        tt= self.OTT_Calibration(commandAmpVector_ForPARRMAlignement, 
                                       nPushPull_ForPARRMAlignement, maskIndex_ForPARRMAlignement)
        cmd= self.OTT_Alignement(tt)
        self._moveRM()
        zernikeCoefComa, comaSurface= self._measureComaOnSegmentMask()
        self._tt= self._cal.measureCalibrationMatrix(3, commandAmpVector_ForM4Calibration, nPushPull_ForM4Calibration)
        intMat, rec= self._cal.analyzerCalibrationMeasurement(self._tt, maskIndex_ForM4Alignement)
        return self._tt, zernikeCoefComa, comaSurface
    
    def M4_Alignement(self, zernikeCoefComa, tt= None):
        if tt is None:
            a= Opt_Alignment(self._tt)
        else:
            a= Opt_Alignment(tt)
        cmd= a.opt_align(zernikeCoefComa)
        return cmd
        
    def _moveRM(self):
        pass
    
    def _moveSegmentView(self):
        pass
    
    def _applyCmd(self):
        pass
    
    def _measureComaOnSegmentMask(self):
        ima= obj.readImageFromFitsFileName('Allineamento/20191001_081344/img.fits')
        roi= self._r.roiGenerator(ima)
        segmentIma= np.ma.masked_array(ima.data, mask= roi[11])
        
        coef, mat= self._zOnM4.zernikeFit(segmentIma, np.arange(2,11))
        coma= coef[5]
        comaSurface= self._zOnM4.zernikeSurface(np.array([coma]), segmentIma.mask, mat, np.array([5]))
        return coma, comaSurface
        
        
    
    

    
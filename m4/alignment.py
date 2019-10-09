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
        
        
    def OTT_Alignment(self, commandAmpVector, nPushPull):
        '''
            arg:
                commandAmpVector= vettore contenente i valori dei movimenti dei 5 gradi di libertà
                nPushPull= numero di push pull per ogni grasi di libertà
        '''
        self._moveRM()
        tt= self._cal.measureCalibrationMatrix(0, commandAmpVector, nPushPull)
        intMat, rec= self._cal.analyzerCalibrationMeasurement(tt)
        a= Opt_Alignment(tt)
        cmd= a.opt_align()
        self._applyCmd()
        return cmd
    
    def M4_Alignment(self, commandAmpVector_ForPARRMAlignement, nPushPull_ForPARRMAlignement,
                        commandAmpVector_ForM4Calibration, nPushPull_ForM4Calibration):
        self._moveSegmentView()
        self.OTT_Alignment(commandAmpVector_ForPARRMAlignement, nPushPull_ForPARRMAlignement)
        self._moveRM()
        piston= self._measureComaOnSegmentMask()
        tt= self._cal.measureCalibrationMatrix(3, commandAmpVector_ForM4Calibration, nPushPull_ForM4Calibration)
        a= Opt_Alignment(tt)
        cmd= a.opt_align(piston[0])
        return cmd
        
    def _moveRM(self):
        pass
    
    def _moveSegmentView(self):
        pass
    
    def _applyCmd(self):
        pass
    
    def _measureComaOnSegmentMask(self):
        ima= obj.readImageFromFitsFileName('Allineamento/20191001_081344/img.fits')
        roi= self._r.ROIonAlignmentImage(ima)
        segmentIma= np.ma.masked_array(ima.data, mask= roi[1])
        
        coef, mat= self._zOnM4.zernikeFit(segmentIma, np.arange(2,11))
        coma= coef[5]
        comaSurface= self._zOnM4.zernikeSurface(np.array([coma]), segmentIma.mask, mat, np.array([5]))
        return coma, comaSurface
        
        
    
    

    
'''
@author: cs
'''


from m4.utils.opticalAlignment import Opt_Alignment
from m4.utils.opticalCalibration import Opt_Calibration


class Alignment():
    
    def __init__(self):
        self._cal= Opt_Calibration()
        
        
    def OTT_Alignment(self, commandAmpVector, nPushPull):
        self._moveRM()
        tt= self._cal.measureCalibrationMatrix(0, commandAmpVector, nPushPull)
        intMat, rec= self._cal.analyzerCalibrationMeasurement(tt)
        self._all= Opt_Alignment(tt)
        cmdf, cmdt= self._all.opt_align()
        self._applyCmd()
        return cmdf, cmdt
    
    def M4_Alignment(self):
        self._moveSegmentView()
        pass
        
    def _moveRM(self):
        pass
    
    def _moveSegmentView(self):
        pass
    
    def _applyCmd(self):
        pass
    

    
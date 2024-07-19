from m4.configuration import userconfig as myconf
from m4 import main
import numpy as np

class OTTScripts:
    """
    xxx
    
    Methods
    =======
    alignTT
    alignFocus
    alignComa
    
    """

    def __init__(self, ott, interf, dm):
        """The Constructor"""
        print('Tell me, Master')
        self._ott = ott
        self._interf = interf
        self._dm = dm
        
    def _checkAlignmInfo(self, move, removePar):
        if move == 0:
            print('Running in test mode, specify move=1 to actually align the OTT')
            doit = False
        else:
            doit = True
        if removePar == True:
            tnPar = myconf.remappedpar_tn
        else:
            tnPar = None
        return doit, tnPar
    
    def alignTT(self, nframes, move=0, removePar=True):
        doit, tnPar = self._checkAlignmInfo(move, removePar)
        zern2corrf = np.array([0,1])
        dofidf = np.array([3,4])
        tna = main.align_PARAndRM(self._ott, self._interf, myconf.alignmentCalibration_tn, zern2corrf, dofidf, nframes, doit, tnPar)

    def alignFocus(self, nframes, move=0, removePar=True):
        doit, tnPar = _checkAlignmInfo(move, removePar)
        zern2corrf = np.array([2])
        dofidf = np.array([0])
        tna = main.align_PARAndRM(self._ott, self._interf, myconf.alignmentCalibration_tn, zern2corrf, dofidf, nframes, doit, tnPar)

    def alignComa(self, nframes, move=0, removePar=True):
        if nframes < 5:
            print('The requested number of frames is not enough to guarantee a robust alignment. Trying...')
        doit, tnPar = _checkAlignmInfo(move, removePar)
        zern2corrf = np.array([3,4])
        dofidf = np.array([1,2,3,4])
        tna = main.align_PARAndRM(self._ott, self._interf, myconf.alignmentCalibration_tn, zern2corrf, dofidf, nframes, doit, tnPar)

    def calibrateAlignment(self, nPushPull = 2, n_frames=20):
        par_pist = myconf.alignCal_parPist
        par_tip, par_tilt = myconf.alignCal_parTip, myconf.alignCal_parTilt
        rm_tip, rm_tilt = myconf.alignCal_rmTip, myconf.alignCal_rmTilt
        command_amp_vector = np.array([par_pist, par_tip, par_tilt, rm_tip, rm_tilt])
        print('Calibrating the OTT alignment, with the following command amplitudes')
        print(command_amp_vector)
        tnc = 'xxxx'
        #tnc = main.calibrate_PARAndRM(self._ott, self._interf, command_amp_vector, nPushPull, n_frames, delay=0)

        return tnc

    def config4D4Alignment(self):
        print('Applying 4D configuration file: '+myconf.phasecam_alignmentconfig)
        self._interf.loadConfiguration(myconf.phasecam_alignmentconfig)

    def config4D4Markers(self):
        print('Applying 4D configuration file: '+myconf.phasecam_markerconfig)
        self._interf.loadConfiguration(myconf.phasecam_markerconfig)


    def shiftAndTrackFlat(self):
        '''
        moves the RM and adjusts the alignment at each steps
        '''
        pass

    def shiftAndTrackTruss(self):
        pass

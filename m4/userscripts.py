import numpy as np
from m4.configuration import userconfig as myconf
from m4 import main, noise
from m4.configuration import update_folder_paths as ufp
from m4.utils import markers as mrk
from m4.mini_OTT import measurements
fn = ufp.folders

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
#       print('Tell me, Master')
        self._ott = ott
        self._interf = interf
        self._dm = dm
        self._meas = measurements.Measurements(ott, interf)

    def alignTT(self, nframes, move=0, removePar=True):
        doit, tnPar = self._checkAlignmInfo(move, removePar)
        zern2corrf = np.array([0,1])
        dofidf = np.array([3,4])
        tna = main.align_PARAndRM(self._ott, self._interf,
                                  myconf.alignmentCalibration_tn, zern2corrf,
                                  dofidf, nframes, doit, tnPar)

    def alignFocus(self, nframes, move=0, removePar=True):
        doit, tnPar = self._checkAlignmInfo(move, removePar)
        zern2corrf = np.array([2])
        dofidf = np.array([0])
        tna = main.align_PARAndRM(self._ott, self._interf,
                                  myconf.alignmentCalibration_tn, zern2corrf,
                                  dofidf, nframes, doit, tnPar)

    def alignComa(self, nframes, move=0, removePar=True):
        if nframes < 5:
            print('The requested number of frames is not enough to guarantee a\
                  robust alignment. Trying...')
        doit, tnPar = self._checkAlignmInfo(move, removePar)
        zern2corrf = np.array([3,4])
        dofidf = np.array([1,2,3,4])
        tna = main.align_PARAndRM(self._ott, self._interf,
                                  myconf.alignmentCalibration_tn, zern2corrf,
                                  dofidf, nframes, doit, tnPar)

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

    def acquireNoise(self):
        self._interf.loadConfiguration(myconf.phasecam_noiseconfig)
        tn = self._interf.capture(myconf.noise_nframes)
        self._interf.produce(tn)
        self._interf.loadConfiguration(myconf.phasecam_baseconfig)
        dfpath = fn.OPD_IMAGES_ROOT_FOLDER+'/'+tn+'/'
        noise.convection_noise(dfpath, myconf.noise_tau_vector)
        noise.noise_vibrations(dfpath,myconf.noise.difftemplate)
        return tn

    def acquireTimeAverage(self, nframes, delay=2):
        tn = self._meas.opticalMonitoring(nframes, delay)
        return tn

    def acquireCurrentFootprint(self):
        c0 = mrk.measureMarkerPos(None, self._interf)
        return c0

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
    

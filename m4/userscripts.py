import numpy as np
from m4.configuration import userconfig as myconf
#from m4 import main, noise #main is no longer required for alignment
from m4 import noise
from m4 import main  #this is to be removed
from m4.configuration import update_folder_paths as ufp
from m4.utils import markers as mrk
from m4.mini_OTT import measurements
from m4.analyzers import timehistory as th
from m4.utils.alignment import Alignment
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
        self.alignment = Alignment(ott, interf)

    def generalAlignment(self, zz, dof,nframes, move=0, removePar=True):
        self.config4D4Alignment()
        doit, tnPar = self._checkAlignmInfo(move, removePar)
        #zern2corrf = np.array([0,1])
        #dofidf = np.array([3,4])
        if removePar == True:  #qui bisogna aggiungere il Tn dell'allineamento!!
            self.alignment.reload_calibrated_parabola(tnPar)
        self.alignment.correct_alignment(dof, zz, move, nframes)

    def alignTT(nframes, move=0, removePar=True):
        zz = np.array([0,1])
        dd = np.array([3,4])
        self.generalAlignment(zz, dd,nframes, move, removePar)

    def alignComa(nframes, move=0, removePar=True):
        zz = np.array([0,1,6,7])
        dd = np.array([1,2,3,4])
        self.generalAlignment(zz, dd,nframes, move, removePar)

    def alignFocus(nframes, move=0, removePar=True):
        zz = np.array([2])
        dd = np.array([0])
        self.generalAlignment(zz, dd,nframes, move, removePar)  
    '''
    def alignTT(self, nframes, move=0, removePar=True):
        self.config4D4Alignment()
        doit, tnPar = self._checkAlignmInfo(move, removePar)
        zern2corrf = np.array([0,1])
        dofidf = np.array([3,4])
        if removePar == True: 
            self.alignment.reloadCalibratedParabola(tnPar)
        self.alignment.correct_alignment(dofidf, zern2corrf, move, nframes)

        #tna = main.align_PARAndRM(self._ott, self._interf,                                  myconf.alignmentCalibration_tn, zern2corrf,                                  dofidf, nframes, doit, tnPar)  #obsolete

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
    '''
    def calibrateAlignment(self, nPushPull = 2, n_frames=20, removePar = True):
        doit, tnPar = self._checkAlignmInfo(1, removePar)
        par_pist = myconf.alignCal_parPist
        par_tip, par_tilt = myconf.alignCal_parTip, myconf.alignCal_parTilt
        rm_tip, rm_tilt = myconf.alignCal_rmTip, myconf.alignCal_rmTilt
        command_amp_vector = np.array([par_pist, par_tip, par_tilt, rm_tip, rm_tilt])
        print('Calibrating the OTT alignment, with the following command amplitudes and Parabola removal option:')
        print(command_amp_vector)
        print(tnPar)
        tnc = main.calibrate_PARAndRM(self._ott, self._interf, command_amp_vector, nPushPull, n_frames, delay=0, tnPar=tnPar)
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

    def analyzeTimeAverage(self, tn, zern2remove=[1,2,3]):
        img = th.averageFrames(tn)
        img = th.zernike.removeZernike(img, zern2remove)
        return img

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

    def prepareFlat(self, tniff):
        #self._flat=flattening(tniff)
        #return self._flat
        pass

    def flattening(self, tniff, nmodes, nframes,zern2remove=[1,2,3]):
        #flatclass=flattening(tniff)
        #tn = self.acquireTimeAverage(self, nframes)
        #img = self.analyzeTimeAverage(self, tn, zern2remove)
        #flatcmd = flatclass.(tniff, img, nmodes, zern2remove)
        #flatcmd is the ampl vector referred to the IFF cmdmatrix
        #cmdMat2use = flatclass._cmdMat[:,0:nmodes]....

        #command = cmdmat2use @ flatcmd
        #now command is just nmodes. ---> fill the entire command vector???
        #dm.mirrorCommand(-flatcmd)
        pass

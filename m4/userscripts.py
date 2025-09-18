"""
Author(s):
----------
- Runa Briguglio - INAF - Osservatorio Astrofisico di Arcetri : written in 2024

Description
-----------
This module provides "user-level" scripts for using the M4 Optical Test Tower (OTT) and its devices, in the frame of the M4 Optical Calibration and Verification.
The scripts consist of high level functions called/initialized with pre-defined parameters and configurations, which are specifically set for that specific, user-level action.
The configurations are defined in a userconfiguration.py file.
Such user-scripts include:
    - Scripts for setting up the OTT for specific actions such as alignment or DM measurement.
    - Scripts for configuration of the OTT interferometer with pre-defined configuration file.
    - Scripts for the calibration and correction of the optical alignment.
    - Scripts for collecting noise and monitoring measurements and to analyze them (TBC)
    - Scripts for commanding at high level the DM for common operations

How to Use it
-------------
Hardware Requirement:
    - The OTT shall be fully functional, including PLC connected, PAR and RefMirr actuators powered up and connected, PAR and RefMirr translation stages homed,...
    - The 4SightFocus Sw shall be up and running on the 4D PC and the shared folders mounted in the calibration workstation
    - The DM is not mandatory, as the scripts may be initialized also with no DM defined,. This configuration allows full control of the OTT and interferometer for alignment and monitoring.

Usage Example
-------------
>>> import m4
>>> ott,_,interf=m4.create_ott()
>>> from m4 import userscripts
>>> ottuser = OTTScripts(ott, interf)
>>> # everything is ready
>>> 

"""
import numpy as np
from m4.configuration import userconfig as myconf
#from m4 import main, noise  # main is no longer required for alignment
#from m4 import noise
#from m4 import main  # this is to be removed
#from m4.configuration import update_folder_paths as ufp

#from m4.utils import markers as mrk
#from m4.mini_OTT import measurements
#from m4.analyzers import timehistory as th
#from m4.utils.alignment import Alignment
from opticalib import alignment
from m4.devices import opt_beam

#fn = ufp.folders


class OTTScripts:
    """
    xxx

    Methods
    =======
    alignTT
    alignFocus
    alignComa

    """

    def __init__(self, ott=None, interf=None, dm=None):
        """The Constructor"""
        #       print('Tell me, Master')
        self._ott = ott
        self._interf = interf
        self._dm = dm
        #self._meas = measurements.Measurements(ott, interf)
        self.alignment = Alignment(ott, interf)
        self.collimator= opt_beam.Parabola(ott)
        self.refMirror = opt_beam.ReferenceMirror(ott)

    def configureOTT4Alignment(self):
        """
        This function moves the Reference Mirror to the defined position in the OTT to allow the optical alignment. 
        Parameters
        ----------

        Returns
        -------
        """
        print('Moving Reference Mirror to '+str(myconf.rmslider4alignment))
        self.refMirror.moveRmsTo(myconf.rmslider4alignment)

    def configureOTTrefMirrorOut(self):
        """
        This function moves the Reference Mirror to the defined position outside the PAR footprint when measuring M4.
        Parameters
        ----------

        Returns
        -------
        """
        print('Moving Reference Mirror outside the beam to'+str(myconf.rmsliderout))
        self.refMirror.moveRmsTo(myconf.rmsliderout)

    def configureOTT4Segment(self):
        """
        This function moves the Parabolic Mirror to the defined position in the OTT to view the M4 segment.
        Parameters
        ----------

        Returns
        -------
        """
        print('Moving Truss to '+str(myconf.parslider4segment))
        self.collimator.moveTrussTo(myconf.parslider4segment)

    def generalAlignment(self, zz, dof, nframes, move=0, removePar=True):
        """
        This function is a general script for managing the alignment, to be used in different configurations to align different items and different Zernikes in the OTT
        Parameters
        ----------
        zz   :  list
            List of Zernike indexes to be corrected. Zernike modes are organized in this scope with the following indexing:
            0 - Tip
            1 - Tilt
            2 - Focus
            3 - ComaX
            4 - ComaY
        dof  : list
            List of indexes of the Degrees of Freedom to be corrected. DoF are organized according with the following convention, which in turn is that adopted in the alignment configuration file:
            0 - PAR dZ == PAR Piston
            1 - PAR rX == PAR Tip
            2 - PAR rY == PAR Tilt
            3 - RM  rX == RefMirr Tip
            4 - RM  rY == RefMirr Tilt
            5 - M4  rX == M4 Tip
            6 - M4  rY == M4 Tip
        nframes   :  int
            Number of frames to be averaged to acquire the image for computing the Zernike to be corrected
        move      :  int
           Flag to indicate whether to apply the correction movement (1) or not (0)
        removePar :  bool
            Flag to indicate whether the PAR figure shall be removed or not.
        Returns
        -------
        """
        self.config4D4Alignment()
        doit, tnPar = self._checkAlignmInfo(move, removePar)
        # zern2corrf = np.array([0,1])
        # dofidf = np.array([3,4])
        if removePar == True:  # qui bisogna aggiungere il Tn dell'allineamento!!
            self.alignment.reload_calibrated_parabola(tnPar)
        self.alignment.correct_alignment(dof, zz, move, nframes)

    def alignM4TT(self, nframes, move=0, removePar=True):
        """
        This function moves the M4 Hexapod to correct the TipTilt on the M4 visible segment. TBD!! This part shall be expanded to manage the different rotation angles in the OTT,and also the CenterView.
        Parameters
        ----------
        nframes   :  int
            Number of frames to be averaged to acquire the image for computing the Zernike to be corrected
        move      :  int
           Flag to indicate whether to apply the correction movement (1) or not (0)
        removePar :  bool
            Flag to indicate whether the PAR figure shall be removed or not.
        Returns
        -------
        """    
        zz = np.array([0, 1])
        dd = np.array([5, 6])
        self.generalAlignment(zz, dd, nframes, move, removePar)

    def alignOTT(self, nframes, move=0, removePar=True):
        """
        This function moves the Parabolic Mirror to the defined position in the OTT to view the M4 segment.
        Parameters
        ----------
        nframes   :  int
            Number of frames to be averaged to acquire the image for computing the Zernik
e to be corrected
        move      :  int
           Flag to indicate whether to apply the correction movement (1) or not (0)
        removePar :  bool
            Flag to indicate whether the PAR figure shall be removed or not.           
        Returns
        -------
        """  
        zz = np.array([0, 1])
        dd = np.array([3, 4])
        self.generalAlignment(zz, dd, nframes, move, removePar)

    def alignComa(self, nframes, move=0, removePar=True):
        """
        This function moves the Parabolic Mirror to the defined position in the OTT to view the M4 segment.
        Parameters
        ----------
        nframes   :  int
            Number of frames to be averaged to acquire the image for computing the Zernik
e to be corrected
        move      :  int
           Flag to indicate whether to apply the correction movement (1) or not (0)
        removePar :  bool
            Flag to indicate whether the PAR figure shall be removed or not.           
        Returns
        -------
        """
        zz = np.array([0, 1, 6, 7])
        dd = np.array([1, 2, 3, 4])
        self.generalAlignment(zz, dd, nframes, move, removePar)

    def alignFocus(self, nframes, move=0, removePar=True):
        """
        This function moves the Parabolic Mirror to the defined position in the OTT to view the M4 segment.
        Parameters
        ----------
        nframes   :  int
            Number of frames to be averaged to acquire the image for computing the Zernik
e to be corrected
        move      :  int
           Flag to indicate whether to apply the correction movement (1) or not (0)
        removePar :  bool
            Flag to indicate whether the PAR figure shall be removed or not.          
        Returns
        -------
        """
        zz = np.array([2])
        dd = np.array([0])
        self.generalAlignment(zz, dd, nframes, move, removePar)

    """
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
    """

    def calibrateOTTAlignment(self, cmdAmp, n_frames, save=True):
        """
        xxx
        Parameters
        ----------

        Returns
        -------
        """
        doit, tnPar = self._checkAlignmInfo(1, removePar)
        par_pist = myconf.alignCal_parPist
        par_tip, par_tilt = myconf.alignCal_parTip, myconf.alignCal_parTilt
        rm_tip, rm_tilt = myconf.alignCal_rmTip, myconf.alignCal_rmTilt
        command_amp_vector = np.array([par_pist, par_tip, par_tilt, rm_tip, rm_tilt])
        print(
            "Calibrating the OTT alignment, with the following command amplitudes and Parabola removal option:"
        )
        print(command_amp_vector)
        print(tnPar)
        print('Add here the command for OTT alignment calibration == PAR + RM')
        #tncal = al.calibrate_PARAndRM( command_amp_vector,  n_frames   )
        return tncal

    def config4D4Alignment(self):
        """
        xxx
        Parameters
        ----------

        Returns
        -------
        """
        print("Applying 4D configuration file: " + myconf.phasecam_alignmentconfig)
        self._interf.loadConfiguration(myconf.phasecam_alignmentconfig)

    def config4D4Markers(self):
        """
        xxx
        Parameters
        ----------

        Returns
        -------
        """
        print("Applying 4D configuration file: " + myconf.phasecam_markerconfig)
        self._interf.loadConfiguration(myconf.phasecam_markerconfig)

    def shiftAndTrackFlat(self):
        """
        moves the RM and adjusts the alignment at each steps
        """
        pass

    def shiftAndTrackTruss(self):
        pass

    def acquireNoise(self):
        self._interf.loadConfiguration(myconf.phasecam_noiseconfig)
        tn = self._interf.capture(myconf.noise_nframes)
        self._interf.produce(tn)
        self._interf.loadConfiguration(myconf.phasecam_baseconfig)
        dfpath = fn.OPD_IMAGES_ROOT_FOLDER + "/" + tn + "/"
        noise.convection_noise(dfpath, myconf.noise_tau_vector)
        noise.noise_vibrations(dfpath, myconf.noise.difftemplate)
        return tn

    def acquireTimeAverage(self, nframes, delay=2):
        tn = self._meas.opticalMonitoring(nframes, delay)
        return tn

    def analyzeTimeAverage(self, tn, zern2remove=[1, 2, 3]):
        img = th.averageFrames(tn)
        img = th.zernike.removeZernike(img, zern2remove)
        return img

    def acquireCurrentFootprint(self):
        c0 = mrk.measureMarkerPos(None, self._interf)
        return c0

    def _checkAlignmInfo(self, move, removePar):
        if move == 0:
            print("Running in test mode, specify move=1 to actually align the OTT")
            doit = False
        else:
            doit = True
        if removePar == True:
            tnPar = myconf.remappedpar_tn
        else:
            tnPar = None
        return doit, tnPar

    def prepareFlat(self, tniff):
        # self._flat=flattening(tniff)
        # return self._flat
        pass

    def flattening(self, tniff, nmodes, nframes, zern2remove=[1, 2, 3]):
        # flatclass=flattening(tniff)
        # tn = self.acquireTimeAverage(self, nframes)
        # img = self.analyzeTimeAverage(self, tn, zern2remove)
        # flatcmd = flatclass.(tniff, img, nmodes, zern2remove)
        # flatcmd is the ampl vector referred to the IFF cmdmatrix
        # cmdMat2use = flatclass._cmdMat[:,0:nmodes]....

        # command = cmdmat2use @ flatcmd
        # now command is just nmodes. ---> fill the entire command vector???
        # dm.mirrorCommand(-flatcmd)
        pass

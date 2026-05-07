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
import opticalib
import os
from opticalib.core.read_config import load_yaml_config as lya
#from m4.configuration import userconfig as myconf
## patch to work from MicWs
#import Microgate.utils.setupLog as setupLog
#setupLog.consoleProfile()


# from m4 import main, noise  # main is no longer required for alignment
# from m4 import noise
# from m4 import main  # this is to be removed
# from m4.configuration import update_folder_paths as ufp

# from m4.utils import markers as mrk
# from m4.mini_OTT import measurements
# from m4.analyzers import timehistory as th
# from m4.utils.alignment import Alignment
# from opticalib import alignment
from m4 import alignment
from m4.devices import opt_beam
from opticalib.dmutils import (
    iff_module as ifm,
    iff_processing as ifp,
    iff_acquisition_preparation as ifa,
)
from opticalib.dmutils.flattening import Flattening
from opticalib import analyzer as imgaz
from opticalib.ground import modal_decomposer as mdl
from scripts.misc import ott_measurements as measurements

configfilename = 'userconfiguration.yaml'

# fn = ufp.folders

def read_userconfig(masterkey,label="BASE"):
    myconf = lya(os.path.join(opt.folders.CONFIGURATION_FOLDER,configfilename))
    theconf = myconf[masterkey]
    z = theconf.get(label)
    if z == None:
        return theconf
    else:
        return z

class OTTScripts:
    """
    xxx

    Methods
    =======
    configureOTT4Alignment
    configureOTTrefMirrorOut
    configureOTT4Segment
    whereisRef
    whereisBeam
    alignM4TT
    alignOTT
    alignComa
    alignFocus
    calibrateOTTAlignment
    calibrateM4Alignment
    config4D4Alignment
    config4D4Markers

    acquireNoise
    acquireTimeAverage
    analyzeTimeAverage
    acquireCurrentFootprint

    """

    def __init__(self, ott=None, interf=None):
        """The Constructor"""
        print('Tell me, Master')
        self._ott    = ott
        self._interf = interf
        self.meas    = measurements.Measurements(interf, ott)
        self.alignment  = alignment.OttAligner(ott, interf)
        self.collimator = opt_beam.Parabola(ott)
        self.refMirror  = opt_beam.ReferenceMirror(ott)
        myconf4d = read_userconfig('CONFIGURATION4D')
        myconfott= read_userconfg('OTTMECH')
        myconfottcal = read_userconf('OTTCAL')

    def configureOTT4Alignment(self):
        """
        This function moves the Reference Mirror to the defined position in the OTT to allow the optical alignment.
        Parameters
        ----------

        Returns
        -------
        """
        print("Moving Reference Mirror to " + str(myconf4d['rmslider4alignment']))
        self.refMirror.moveRmsTo(myconf4d['rmslider4alignment'])

    def configureOTTrefMirrorOut(self):
        """
        This function moves the Reference Mirror to the defined position outside the PAR footprint when measuring M4.
        Parameters
        ----------

        Returns
        -------
        """
        conf = myconfott['rmsliderout']
        print("Moving Reference Mirror outside the beam to" + str(conf))
        self.refMirror.moveRmsTo(conf)

    def configureOTT4Segment(self):
        """
        This function moves the Parabolic Mirror to the defined position in the OTT to view the M4 segment.
        Parameters
        ----------

        Returns
        -------
        """
        conf = myconfott['parslider4segment']
        print("Moving Truss to " + str(conf)
        self.collimator.moveTrussTo(conf)

    def whereisRef(self):
        """
        This function returns the position of the Ref Mirror in the OTT in M4-centered coordinates.
        Parameters
        ----------

        Returns
        -------
        pp    :  float
            RefMirror position
        """
        pp = self.refMirror.rmsGetPosition()
        print(str(pp))
        return pp

    def whereisBeam(self):
        """
        This function returns the position of the Parabolic mirror mount (Truss) in the OTT in M4-centered coordinates.

        Parameters
        ----------

        Returns
        pp    : float
            PAR mirror position
        -------
        """
        pp = self.collimator.trussGetPosition()
        print(str(pp))
        return pp

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
        conf = myconfottcal['alignmentCalibration_tn']
        cavity_or_dm = 'cavity'
        if dof is [5,6]:
            cavity_or_dm = 'dm'
        self.config4D4Alignment(cavity_or_dm)
        self.alignment.load_calibration(conf)
        doit, tnPar = self._checkAlignmInfo(move, removePar)
        if removePar == True:  # qui bisogna aggiungere il Tn dell'allineamento!!
            print("Reload fitting_surface, ToBeChecked!")
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
                This function correct the RefMirror orientation to minimioze the tilt as seen on the RefMirror footprint.
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
                This function corrects the Parabola and Ref Mirror orientations to minimize Tilt and Coma, as measured over the RefMirror footprint.
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
        zz = np.array([0, 1, 3, 4])
        dd = np.array([1, 2, 3, 4])
        self.generalAlignment(zz, dd, nframes, move, removePar)

    def alignFocus(self, nframes, move=0, removePar=True):
        """
                This function corrects the vertical (piston) position of the Parabola to minimize the focus as measured over the Ref Mirror footprint.

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

    def calibrateOTTAlignment(self,  n_frames= 25, save=True):
        """
                This function measures the alignment calibration for the OTT, as the Zernike amplitude (fotted over the RefMirror area and with the PAR as pupil, corresponding to the displacement of each relevant DoF in the OTT. In particular, the PAR dZ, rX and rY and RefMirr rX and rY are moved, while measuring tip, tilt, focus and both comaX and comaY.
                Parameters
                ----------
                cmdAmp     :   list
                    Amplitude, in natural dimensions for the PLC [mm and arcsec], for each DoF. M4 tip tilt are zeroed to skip the execution of M4 hexapod commands.
                n_frames   :   int
                    Number of frames to be averaged to acquire the image for computing the Zernik
        e to be corrected
                save      :  bool
                    Flag to indicate whether to save or not the data


                Returns
                -------
        """
        
        self.config4D4Alignment(cavity_or_dm = 'cavity')
        doit, tnPar = self._checkAlignmInfo(1, removePar)
        par_pist = myconfott['alignCal_parPist']
        par_tip, par_tilt = myconfott['alignCal_parTip'], myconfott['alignCal_parTilt']
        rm_tip, rm_tilt = myconfott['alignCal_rmTip'], myconfott['alignCal_rmTilt']
        m4_tip, m4_tilt = myconfott['alignCal_m4Tip'] * 0, myconfott['alignCal_m4Tilt'] * 0

        command_amp_vector = np.array(
            [par_pist, par_tip, par_tilt, rm_tip, rm_tilt, m4_tip, m4_tilt]
        )
        print(
            "Calibrating the OTT alignment, with the following command amplitudes and Parabola removal option:"
        )
        print(command_amp_vector)
        print(tnPar)
        tncal = al.calibrate_alignment(command_amp_vector, n_frames,save=True)
        return tncal

    def calibrateM4Alignment(self, n_frames, save=True):
        """
                This function measures the alignment calibration for the OTT, as the Zernike amplitude (fotted over the RefMirror area and with the PAR as pupil, corresponding to the displacement of each relevant DoF in the OTT. In particular, the PAR dZ, rX and rY and RefMirr rX and rY are moved, while measuring tip, tilt, focus and both comaX and comaY.
                Parameters
                ----------
                cmdAmp     :   list
                    Amplitude, in natural dimensions for the PLC [mm and arcsec], for each DoF. PAR and RefMirr Dof are zeroed to skip the execution of the associated commands.
                n_frames   :   int
                    Number of frames to be averaged to acquire the image for computing the Zernik
        e to be corrected
                save      : bool
                    Flag to indicate whether to save or not the data

                Returns
                -------
        """
        self.config4D4Alignment(cavity_or_dm = 'dm')
        doit, tnPar = self._checkAlignmInfo(1, removePar)

        par_pist = myconfott['alignCal_parPist']*0
        par_tip, par_tilt = myconfott['alignCal_parTip']*0, myconfott['alignCal_parTilt']*0
        rm_tip, rm_tilt = myconfott['alignCal_rmTip']*0, myconfott['alignCal_rmTilt']*0
        m4_tip, m4_tilt = myconfott['alignCal_m4Tip'] , myconfott['alignCal_m4Tilt'] 


)

        command_amp_vector = np.array(
            [par_pist, par_tip, par_tilt, rm_tip, rm_tilt, m4_tip, m4_tilt]
        )
        print(
            "Calibrating the OTT alignment, with the following command amplitudes and Parabola removal option:"
        )
        print(command_amp_vector)
        print(tnPar)
        tncal = al.calibrate_alignment(cmdAmp, n_frames,save=True)
        return tncal

    def currentMisalignment(self, cavity_or_dm = 'cavity', nframes = 16, removePar = True):
        '''
        '''
        self.config4D4Alignment(cavity_or_dm)
        doit, tnPar = self._checkAlignmInfo(0, removePar)
        zernres = []
        for in range(nframes):
            fullimg = self.interf.acquireFullFrame(nframes)
            zernres.append(al._zern_routine(fullimg))
        zernres = np.array(zernres)
        zernm = np.mean(zernres,0)
        zernerr(np.std(zernres,0))
        print('Residual alignment (Tip, Tilt, Focus, ComaX, ComaY) & Measurement error:')
        print(zernm)
        print(zernerr)
        return zernm, zernerr


    def config4D4Alignment(self, cavity_or_dm = 'cavity'):
        """
        The function communicates with the 4D Focus SW to load the interferometer configuration file to enable the optical alignment.
        Parameters
        ----------

        Returns
        -------
        """
        
        if cavity_or_dm == 'cavity':
            cfile = myconf4d['OTTalignmentconfig']
        if cavity_or_dm == 'dm':
            cfile = myconf4d['M4alignmentconfig']
        print("Applying 4D configuration file: " + cfile)
        self._interf.loadConfiguration(cfile)


    def config4D4Segment(self, segment = None):
        """
        The function communicates with the 4D Focus SW to load the interferometer configuration file to enable the optical alignment.
        Parameters
        ----------

        Returns
        -------
        """
        if segment is None:
            cfile = myconf4d['segmentconfig']
        print("Applying 4D configuration file: " + cfile)
        self._interf.loadConfiguration(cfile)

    def config4D4Markers(self):
        """
        The function communicates with the 4D Focus SW to load the interferometer configuration file to enable the acquisition of the markers.

        Parameters
        ----------

        Returns
        -------
        """
        conf = myconf4d['markerconfig']
        print("Applying 4D configuration file: " + conf)
        self._interf.loadConfiguration(conf)


    def _checkAlignmInfo(self, move, removePar):
        if move == 0:
            print("Running in test mode, specify move=1 to actually align the OTT")
            doit = False
        else:
            doit = True
        if removePar == True:
            tnPar = myconfottcal['remappedpar_tn']
        else:
            tnPar = None
        return doit, tnPar

    def shiftAndTrackFlat(self):
        """
        moves the RM and adjusts the alignment at each steps
        """
        pass

    def shiftAndTrackTruss(self):
        pass



class MeasurementScripts:
    """
    xxx

    Methods
    =======
    acquireNoise
    acquireTimeAverage
    analyzeTimeAverage
    acquireCurrentFootprint

    """

    def __init__(self, ott=None, interf=None, dm=None):
        """The Constructor"""
        print('Tell me, Master')
        self._ott    = ott
        self._interf = interf
        self._dm     = dm
        self.meas    = measurements.Measurements(interf, ott)
        self.alignment  = alignment.OttAligner(ott, interf)
        #self.collimator = opt_beam.Parabola(ott)
        #self.refMirror  = opt_beam.ReferenceMirror(ott)
        myconf4d = read_userconfig('CONFIGURATION4D')
        myconfott= read_userconfg('OTTMECH')
        myconfottcal = read_userconf('OTTCAL')
        myconfmeas   = read_userconf('MEASUREMENT')

    def acquireNoise(self):
        '''
        '''
        self._interf.loadConfiguration(myconf4d['phasecam_noiseconfig'])
        tn = self._interf.capture(myconfmeas['noise_nframes'])
        self._interf.produce(tn)
        self._interf.loadConfiguration(myconf4d['baseconfig'])
        dfpath = fn.OPD_IMAGES_ROOT_FOLDER + "/" + tn + "/"
        noise.convection_noise(dfpath, myconf.noise_tau_vector)
        noise.noise_vibrations(dfpath, myconf.noise.difftemplate)
        return tn

    def acquireTimeSeries(self, nframes, delay=2):
        '''
        '''
        tn = self.meas.opticalMonitoring(nframes, delay)
        return tn

    def analyzeTimeSeries(self, tn, zern2remove=[1, 2, 3],fitmode = 'global'):
        img = imgaz.averageFrames(tn)
        zernfit = mdl.ZernikeFitter(img)
        img = zernfit.removeZernike(img,zern2remove,fitmode)
        return img

    def acquireCurrentFootprint(self):
        c0 = mrk.measureMarkerPos(None, self._interf)
        return c0



class M4Scripts:
    """
    xxx

    Methods
    =======
    initReconstructor
    loadFlatCommand
    relax
    opticalFlat
    acquireModalIFF
    

    """

    def __init__(self, dm=None, interf=None):
        """The Constructor"""
        #       print('Tell me, Master')

        self.interf = interf
        self.dm = dm
        self.ifa = opticalib.dmutils.iff_module
        from opticalib.dmutils.iff_acquisition_preparation import IFFCapturePreparation as ifa
        #self.ifc = self.ifa.IFFCapturePreparation(dm)
        self.flattening = None
        myconf4d = read_userconfig('CONFIGURATION4D')
        myconfott= read_userconfg('OTTMECH')
        myconfottcal = read_userconf('OTTCAL')
        myconfmeas   = read_userconf('MEASUREMENT')
        myconfdmconf   = read_userconf('DM_CONFIG')
        myconfdmmeas   = read_userconf('DM_MEAS')



    def initReconstructor(tn):
        self.flattening = opticalib.dmutils.flattening.Flattening(tn)

    def loadFlatCommand(self,flattn=None, incremental=10):
        """
        The function loads the actuator positions corresponding to a save flattening vector and applies the command to the DM after checking the bias vectors.
        Parameters
        ----------
        tn     :   str
            the Tracking number where the data are saved. If not passed (None), the default flat (tn saved in userconfig) will be applied

        Returns
        -------
        """
        if flattn is None:
            flattn = myconfdmconf['dm_defaultFlatCmd']
        flatfold = opticalib.folders.FLAT_ROOT_FOLDER
        flatfile = os.path.join(flatfold, flattn, 'flatCommand.fits')
        fcmd = opticalib.load_fits(flatfile)
        self.dm.set_shape(fcmd, incremental=incremental)
    
    def relax(self):
        self.dm.set_shape(-self.dm._last_cmd, differential = True, incremental = 10)
    
    def opticalFlat(self,nmodes, segmentId=[0,1], tn=None):
        if tn is None:
            tn = myconfdmcoinf['dm_defaultIFF']
        #if segmentId == [0,1]:
        mid = np.stack((np.arange(nmodes),np.arange(111,111+nmodes)),axis=0).flatten()
        f = opticalib.dmutils.flattening.Flattening(tn)
        f.applyFlatCommand(self.dm, self.interf, mid,modes2discard=2,nframes=4,incremental=10)
        

    def acquireModalIFF(modes, segment, amp=None):
        ampvec = opticalib.load_fits(os.path.join(opticalib.folders.IFFUNCTIONS_ROOT_FOLDER,usr.myconfdmmeas['iff_modal_ampTN'],'ampVector.fits')) if (amp is None) else amp
        tn = self.generalIffAcquisition(modes, segment, amp, npushpull)

        pass

    def generalIffAcquisition(modes, amp, modalbase,npushpull = 3, segment = None, interferometer = self.interf):
        
        ampvec = opticalib.load_fits(os.path.join(opticalib.folders.IFFUNCTIONS_ROOT_FOLDER,usr.myconf.iff_modal_ampTN,'ampVector.fits'))
        template = np.ones(npushpull)
        template[1::2]=-1
        if segment is None:
            mlist = modes
        else:
            mlist = np.arange(modes)+self.dm.nActsPerSegment*segment
        #self.ifa._updateModalBase(modalbase)
        print(mlist)
        print('Interf 2 use')
        print(useinterf)
        tn = opticalib.dmutils.iff_module.iffDataAcquisition(self.dm, interferometer,mlist, amplitude, template, modalbase, shuffle)
        return tn
       
    def fitZernCommand(tn, nmodes, tid, roiid=None,n2discard = 2):
    #print(tn)
        f = self.initReconstructor(tn)
        #f = Flattening(tn)
        m = f._getMasterMask()
        ss = m.shape
        img = np.zeros(ss)
        mm1 = np.ones(ss)
        mm1[m == False] = 0
        img = np.ma.masked_array(img,mm1)
        zfit = md.ZernikeFitter(img)
        roiimg = roi.roiGenerator(img)
        cc, zmat=zfit.fit(img,[1,2,3])
        img.data[img.mask == 0] = zmat[tid,:]
        img.data[roiimg[roiid] == 1] = 0
        img1 = img.copy()
        img1 = np.ma.masked_array(img.data, roiimg[roiid])
        if tid != 0:
            print('Requested mode is not Piston, removing the mean and normalizing to 1 std')
            img1 = img1 - np.ma.mean(img1)
            img1 = img1/img1.std()
        print('Mean, StD')
        print(np.mean(img1));    print(img1.std())
        img.data[roiimg[roiid] == 0] = img1.data[roiimg[roiid] == 0]
        img.data[roiimg[roiid] == 1] = 0
        f.loadImage2Shape(img)
        f.computeRecMat(n2discard)  #nmodes 2 remove
        fcmd = f.computeFlatCmd(nmodes)
        return fcmd



class RequirementScripts:
    """
    xxx

    Methods
    =======
    sampleSlope
    sampleHF
    sampleCurvature
    
    """

    def __init__(self):
        """The Constructor"""
        from m4.analyzers import requirement_analyzer as raz
        self.pixscaleott = 0.00076  # mm/pix in full res 4D
        self.hf_samplerad = 0.015  # mm
        self.rebinfactor = 4
        

    def sampleSlope(self,tn, rebfact = None):
        """
        This function computes the surface slope of a frame against the slope error specification. The average of a time series measurements is computed and the slope is evaluated as the difference between the frame and its rolled one (shifted 1pix in both X and Y), normalized with the pixelscale. The result is rebinned to get rid of the pixel-by-pixel scatter
        Parameters
        ----------
        tn     :   str
            the Tracking number where the data are saved.

        Returns
        -------
        val    : float
            The std of the surface slope error (in arcsec)
        """
        if rebfact is None:
            rebfact = self.rebinfactor
        img = opticalib.analyzer.average_frames(tn)
        dimg = img - _np.roll(img,(1,1),axis=(0,1))
        dimg = az.modeRebinner(dimg,rebfact)
        dimg = dimg/ (pixscaleott*rebfact) * 206265 #conversion surf --> slope --> arcsec
        return dimg

    def sampleHF(self,image, step = 400):
        """
        This function samples a given frame according to the specification for the HF error and returns the 95-th percentile of the sample. The input frame is obtained after calibrating a M4-segment realization with the local OTT calibration
        Parameters
        ----------
        image     :   array
            the calibrated frame (surface error)

        Returns
        -------
        val    : float
            The 95-th percentile of the WF sample
        """
        #req, list_ima, result_vect =  raz.patches_analysis(image, hf_samplerad, 1/pixscaleott, step=step, n_patches=None)
        image = image.copy()*2
        print('Data converted to WF (2x)')
        _, _, resvect =  raz.patches_analysis(image, hf_samplerad, 1/pixscaleott, step=step, n_patches=None)
        return resvect

    def sampleCurf(self, img):
        """
        This function samples a given frame according to the specification for the curvature error and returns the curvature values of the image. The input frame is obtained after calibrating a M4-segment realization with the local OTT calibration
        Parameters
        ----------
        image     :   array
            the calibrated frame (surface error)

        Returns
        -------
        cval    : array
            The curvature values computed over the frame
        """

        cval = raz.patches_analysis_map(img, 0.04, 1/px_ott, cpatch = 60)
        return cval


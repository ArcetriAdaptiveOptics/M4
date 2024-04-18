#!/usr/bin/env python
__version__= "$Id$"


import os
import Queue
import time

import numpy as np
import pyfits


from argos.lgsw.device.frame_buffer_constant import FrameBufferConstant
from argos.util.executor import MultiThreadExecutor, FunctionTask,\
    TaskTimeoutError
from argos.util.snapshotable import Snapshotable
from argos.util.laser_beam_map import LaserBeamMap
from argos.laser_beam import LaserBeam


class SnapshotPrefix(object):

    AARB= "AARB"
    LAN= "LAN"
    LAS= "LAS"
    LGSW= "LGSW"
    TIPTILT= "TIPTILT"
    CAL= "CAL"


class SnapshotEntry(object):

    USER_COMMENT= "USER_COMMENT"
    TELESCOPE_SIDE= "TELESCOPE_SIDE"
    BEAM_BLUE_LASER= "BEAM_BLUE_LASER"
    BEAM_RED_LASER= "BEAM_RED_LASER"
    BEAM_YELLOW_LASER= "BEAM_YELLOW_LASER"


class Snapshooter(object):


    FILE_AO_GAIN= 'ao_gain.fits'
    FILE_APD_COUNTERS= 'apd_counters.fits'
    FILE_ASM_DIST_AVGS= 'asm_dist_average.fits'
    FILE_ASM_INT_MODES= 'asm_integrated_modes.fits'
    FILE_FFCOMMANDS= 'asm_ffcommands.fits'
    FILE_FIELD_CAM_FRAMES= "field_camera.fits"
    FILE_HVC_TT_JITTER= 'jitter.fits'
    FILE_LGSW_CCD_FRAMES= 'lgsw_ccdframes.fits'
    FILE_LGSW_CCD_FRAME_COUNTERS= 'lgsw_ccdframe_counters.fits'
    FILE_OVERVIEW= "overview.fits"
    FILE_PAT_CAM_BLUE= 'patrol_camera_blue.fits'
    FILE_PAT_CAM_RED= 'patrol_camera_red.fits'
    FILE_PAT_CAM_YELLOW= 'patrol_camera_yellow.fits'
    FILE_PIEZO_COMMANDS= 'jitter_piezo_commands.fits'

    FILE_PIEZO_POSITIONS= 'jitter_piezo_positions.fits'
    FILE_PUPIL_CAM_FRAMES = "pupil_camera.fits"
    FILE_OSCILLOSCOPE = 'oscilloscope.fits'
    FILE_SLOPES= 'slopes.fits'
    FILE_SLOPE_FRAME_COUNTERS= 'slopes_frame_counters.fits'
    FILE_SLOPES_OFFSET= 'slopes_offset.fits'
    FILE_CALIBRATION_DEVICES= 'calibration_devices.fits'


    def __init__(self,
                 indexFile,
                 telescopeSide,
                 lasClientBuilder,
                 lasDevicePool,
                 lanClientBuilder,
                 lanDevicePool,
                 lgswClientBuilder,
                 lgswDevicePool,
                 tiptiltClientBuilder,
                 tiptiltDevicePool,
                 aarbClientBuilder,
                 calibrationDevicePool,
                 logger):
        self._indexFile= indexFile
        self._telescopeSide= telescopeSide
        self._lasClientBuilder= lasClientBuilder
        self._lasDevicePool= lasDevicePool
        self._lanClientBuilder= lanClientBuilder
        self._lanDevicePool= lanDevicePool
        self._lgswClientBuilder= lgswClientBuilder
        self._lgswDevicePool= lgswDevicePool
        self._tiptiltClientBuilder= tiptiltClientBuilder
        self._tiptiltDevicePool= tiptiltDevicePool
        self._aarbClientBuilder= aarbClientBuilder
        self._calibrationDevicePool= calibrationDevicePool
        self._logger= logger

        self._currentConfig= None
        self._executor= None
        self._futureQueue= Queue.Queue()
        self._initialized= False
        self._overviewHdr= None
        self._snapshotDir= None


    def initialize(self):
        assert not self._initialized
        self._executor= MultiThreadExecutor(8, self._logger)
        self._initialized= True


    def _createSnapshotDir(self):
        self._snapshotDir= self._currentConfig.getSnapshotDirectory()
        self._logger.notice("Creating snapshot directory %s" % (
                                                       self._snapshotDir))
        os.makedirs(self._snapshotDir)


    def _waitForCompletionOfTasks(self):
        while not self._futureQueue.empty():
            future= self._futureQueue.get_nowait()
            maxTimeoutForAllTasksSec= self._currentConfig.getTimeout()
            future.waitForCompletion(timeoutSec= maxTimeoutForAllTasksSec)


    def _lasCtrlClient(self):
        return self._lasClientBuilder.build()


    def _lanCtrlClient(self):
        return self._lanClientBuilder.build()


    def _tiptiltCtrlClient(self):
        return self._tiptiltClientBuilder.build()


    def _aarbClient(self):
        return self._aarbClientBuilder.build()


    def _snapshotLGSWController(self, lgswOverview):
        lgswOverview.update(self._getLgswCtrlClient().getSnapshot(
            SnapshotPrefix.LGSW))


    def _snapshotLASController(self, lasOverview):
        lasOverview.update(self._lasCtrlClient().getSnapshot(
            SnapshotPrefix.LAS))


    def _snapshotLANController(self, lanOverview):
        lanOverview.update(self._lanCtrlClient().getSnapshot(
            SnapshotPrefix.LAN))


    def _snapshotTIPTILTController(self, tiptiltOverview):
        tiptiltOverview.update(self._tiptiltCtrlClient().getSnapshot(
            SnapshotPrefix.TIPTILT))


    def _snapshotAArbController(self, aarbOverview):
        aarbOverview.update(self._aarbClient().getSnapshot(
            SnapshotPrefix.AARB))


    def _snapshotCalibrationDevices(self, calibrationDevicesOverview):
        devs= \
            {"SWING_ARM": self._calibrationDevicePool.getCalibrationSwingArm(),
             "WHITE_LIGHT": self._calibrationDevicePool.getWhiteLightSource(),
             "DIODE_SOURCE": self._calibrationDevicePool.getDiodeSource()}
        for each in devs:
            calibrationDevicesOverview.update(
                devs[each].getSnapshot("%s.%s" %(SnapshotPrefix.CAL, each)))


    def _addLaserBeamMapTo(self, overview):
        L= LaserBeam
        m= {L.Blue: SnapshotEntry.BEAM_BLUE_LASER,
            L.Red: SnapshotEntry.BEAM_RED_LASER,
            L.Yellow: SnapshotEntry.BEAM_YELLOW_LASER}
        for each in L.ALL_BEAMS:
            overview[m[each]]= int(self._indexFile.lookUp(
                LaserBeamMap.beamLaserIndexFileKey(self._telescopeSide, each)))


    def _buildOverviewHeaderAndDumpIt(self):
        overview= {}
        overview[SnapshotEntry.TELESCOPE_SIDE]= str(self._telescopeSide)
        overview[SnapshotEntry.USER_COMMENT]= self._currentConfig.getComment()
        self._addLaserBeamMapTo(overview)

        tiptiltOverview= {}
        self._createTaskFor(self._snapshotTIPTILTController, tiptiltOverview)
        lgswOverview= {}
        self._createTaskFor(self._snapshotLGSWController, lgswOverview)
        lasOverview= {}
        self._createTaskFor(self._snapshotLASController, lasOverview)
        lanOverview= {}
        self._createTaskFor(self._snapshotLANController, lanOverview)
        aarbOverview= {}
        self._createTaskFor(self._snapshotAArbController, aarbOverview)
        calibrationDevicesOverview= {}
        self._createTaskFor(
            self._snapshotCalibrationDevices, calibrationDevicesOverview)
        try:
            self._waitForCompletionOfTasks()
        except(TaskTimeoutError):
            self._logger.warn("Overview is not complete")

        for each in [calibrationDevicesOverview,
                     lasOverview,
                     lanOverview,
                     lgswOverview,
                     tiptiltOverview,
                     aarbOverview]:
            overview.update(each)

        self._overviewHdr= Snapshotable.asFITSHeader(overview)
        pyfits.writeto(self._snapshotPath(self.FILE_OVERVIEW),
                       np.array([1, 2]), self._overviewHdr)


    def _getOverviewHeader(self):
        return self._overviewHdr.copy()


    def _generateLASSnapshot(self):
        self._createTaskFor(self._takeLASSnapshot)


    def _generateAarbAoLoopSnapshot(self):
        gainVector= self._aarbClient().getAoLoopGain()
        pyfits.writeto(self._snapshotPath(self.FILE_AO_GAIN),
                       gainVector,
                       self._getOverviewHeader())


    def _generateAARBSnapshot(self):
        self._createTaskFor(self._generateAarbAoLoopSnapshot)


    def takeSnapshot(self, config):
        self._logger.debug("Config: %s" % config)
        self._currentConfig= config
        self._createSnapshotDir()
        self._buildOverviewHeaderAndDumpIt()
        self._generateLGSWSnapshot()
        self._generateLASSnapshot()
        self._generateAARBSnapshot()
        self._waitForCompletionOfTasks()


    def _getLasDevicePool(self):
        return self._lasDevicePool


    def _takeFieldCameraSnapshot(self):
        images= self._getLasDevicePool().getFieldCamera().getFrames(
            self._currentConfig.getNumOfFieldCameraFrames())
        pyfits.writeto(self._snapshotPath(self.FILE_FIELD_CAM_FRAMES),
                       images,
                       self._getOverviewHeader())


    def _createTaskFor(self, f, *args, **kwds):
        self._futureQueue.put(self._executor.execute(
            FunctionTask(f, *args, **kwds)))


    def _takePupilCameraSnapshot(self):
        images= self._getLasDevicePool().getPupilCamera().getFrames(
            self._currentConfig.getNumOfPupilCameraFrames())
        pyfits.writeto(self._snapshotPath(self.FILE_PUPIL_CAM_FRAMES),
                       images,
                       self._getOverviewHeader())


    def _takeOsciSnapshot(self):
        channels= 4
        curves= self._getLasDevicePool().getOsci().getCurve(channels)
        pyfits.writeto(self._snapshotPath(self.FILE_OSCILLOSCOPE),
                       curves,
                       self._getOverviewHeader())


    def _lasController(self):
        return self._lasTerminal.controller


    def _takeLASSnapshot(self):
        self._createTaskFor(self._takeFieldCameraSnapshot)
        self._createTaskFor(self._takePupilCameraSnapshot)
        self._createTaskFor(self._takeOsciSnapshot)


    def _snapshotPath(self, fileName):
        return os.path.join(self._snapshotDir, fileName)


    def _getLgswCtrlClient(self):
        return self._lgswClientBuilder.build()


    def _takeFrameBufferSnapshot(self, config):
        howmany= config.getNumOfBCUDiagnosticRecords()
        t1=time.time()
        frameBuffer= self._lgswDevicePool.getWFSCCDFrameBuffer()
        if self._isSlopeTransmissionEnabled():
            loopdata= frameBuffer.getMatchingDiagnosticsFuture(howmany)
        else:
            loopdata= frameBuffer.getOverallDiagnosticsFuture(howmany)
        t2=time.time()
        hdr= self._getOverviewHeader()
        hdr.update('hierarch ACQUISITION_TIME_IN_SECOND', t2 - t1)
        hdr.update('hierarch ACQUISITION_DATE_START', t1)
        hdr.update('hierarch ACQUISITION_DATE_END', t2)
        FB= FrameBufferConstant
        pyfits.writeto(self._snapshotPath(self.FILE_SLOPES),
                       loopdata[FB.KEY_DIAG_ALL_SLOPES], hdr)
        pyfits.writeto(self._snapshotPath(self.FILE_SLOPE_FRAME_COUNTERS),
                       loopdata[FB.KEY_DIAG_SLOPE_FRAME_COUNTERS], hdr)
        pyfits.writeto(self._snapshotPath(self.FILE_HVC_TT_JITTER),
                       loopdata[FB.KEY_DIAG_HVC_TT_JITTER], hdr)
        pyfits.writeto(self._snapshotPath(self.FILE_PIEZO_COMMANDS),
                       loopdata[FB.KEY_DIAG_HVC_COMMANDS], hdr)
        pyfits.writeto(self._snapshotPath(self.FILE_PIEZO_POSITIONS),
                       loopdata[FB.KEY_DIAG_HVC_POSITIONS], hdr)
        pyfits.writeto(self._snapshotPath(self.FILE_APD_COUNTERS),
                       loopdata[FB.KEY_DIAG_APD_COUNTERS], hdr)
        pyfits.writeto(self._snapshotPath(self.FILE_ASM_INT_MODES),
                       loopdata[FB.KEY_DIAG_ASM_INT_MODES], hdr)
        pyfits.writeto(self._snapshotPath(self.FILE_FFCOMMANDS),
                       loopdata[FB.KEY_DIAG_ASM_FF_COMMANDS], hdr)
        pyfits.writeto(self._snapshotPath(self.FILE_ASM_DIST_AVGS),
                       loopdata[FB.KEY_DIAG_ASM_DIST_AVGS], hdr)

        pyfits.writeto(self._snapshotPath(self.FILE_LGSW_CCD_FRAMES),
                       loopdata[FB.KEY_DIAG_CCD_RAW_FRAMES], hdr)
        pyfits.writeto(self._snapshotPath(self.FILE_LGSW_CCD_FRAME_COUNTERS),
                       loopdata[FB.KEY_DIAG_CCD_FRAME_COUNTERS], hdr)


    def _takeCurrentSlopesOffsetVectorSnapshot(self, config):
        hdr= self._getOverviewHeader()
        slopesOffset= self._getLgswCtrlClient().getCurrentSlopeOffsetVector()
        pyfits.writeto(self._snapshotPath(self.FILE_SLOPES_OFFSET),
                       slopesOffset.toNumPyArray(), hdr)


    def _takePatrolCameraSnapshot(self, config):
        overviewHdr= self._getOverviewHeader()

        def takeFramesFromBlueCamera():
            patBlue= self._lgswDevicePool.getPatrolCameraBlue()
            pyfits.writeto(self._snapshotPath(self.FILE_PAT_CAM_BLUE),
                      patBlue.getFrames(config.getNumOfPatrolCameraFrames()),
                      overviewHdr)
        self._createTaskFor(takeFramesFromBlueCamera)


        def takeFramesFromRedCamera():
            patRed= self._lgswDevicePool.getPatrolCameraRed()
            pyfits.writeto(self._snapshotPath(self.FILE_PAT_CAM_RED),
                      patRed.getFrames(config.getNumOfPatrolCameraFrames()),
                      overviewHdr)
        self._createTaskFor(takeFramesFromRedCamera)

        def takeFramesFromYellowCamera():
            patYellow= self._lgswDevicePool.getPatrolCameraYellow()
            pyfits.writeto(self._snapshotPath(self.FILE_PAT_CAM_YELLOW),
                      patYellow.getFrames(config.getNumOfPatrolCameraFrames()),
                      overviewHdr)
        self._createTaskFor(takeFramesFromYellowCamera)


    def _generateLGSWSnapshot(self):
        self._futureQueue.put(
            self._executor.execute(FunctionTask(
                self._takeFrameBufferSnapshot, self._currentConfig)))

        self._futureQueue.put(
            self._executor.execute(FunctionTask(
                self._takeCurrentSlopesOffsetVectorSnapshot,
                self._currentConfig)))

        self._futureQueue.put(
            self._executor.execute(FunctionTask(
                self._takePatrolCameraSnapshot, self._currentConfig)))


    def _asFITSHeader(self, dictionary):
        return Snapshotable.asFITSHeader(dictionary)


    def _isSlopeTransmissionEnabled(self):
        return self._getLgswCtrlClient().isSlopeTransmissionEnabled()

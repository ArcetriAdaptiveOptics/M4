"""
Authors
  - C. Selmi: written in 2019
              modified in 2021
"""

import os
import logging
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration import config_folder_names as fold_name
from m4.utils.optical_calibration import OpticalCalibration
from m4.ground import zernike
from m4.configuration.ott_parameters import OttParameters, OtherParameters
from m4.ground.timestamp import Timestamp
from m4.utils import roi
from m4.configuration import ott_status
from m4.ground import geo
from arte.utils.zernike_projection_on_subaperture import (
    ZernikeProjectionOnSubaperture)


class OpticalAlignment:
    """
    Class for the optical alignment

    HOW TO USE IT::

        from m4.utils.optical_alignment import OpticalAlignment
        al = OpticalAlignment(tt, ott, interf)
        par_command, rm_command, dove = al.opt_aligner(n_images, delay,
        zernike_to_be_corrected, dof_command_id)
    """

    def __init__(self, tt_cal, ott, interf):
        """The constructor"""
        self._logger = logging.getLogger("OPT_ALIGN:")
        self._loggerRuna = logging.getLogger("ALIG_ZER:")
        self.tt_cal = tt_cal
        self.cal = OpticalCalibration.loadCalibrationObjectFromFits(tt_cal)
        self._who = self.cal.getWho()
        self._mask = self.cal.getMask()
        self._interf = interf
        self._ott = ott
        # post selection
        self._intMatModesVector = None
        self._commandId = None
        self._rec = None
        self._intMat = None
        self._cmat = None
        self._zernikeVectorSelected = None
        # output
        self.par_command = None
        self.rm_command = None

    @staticmethod
    def _storageFolder():
        """Creates the path where to save data"""
        return fold_name.ALIGNMENT_ROOT_FOLDER

    def opt_aligner(
        self,
        n_images,
        delay,
        zernike_to_be_corrected=None,
        dof_command_id=None,
        tnPar = None,
        subapOffsets=None
    ):
        """
        Parameters
        ----------
        n_images: int
            number of interferometers frames
        delay: int [s]
            delay between images

        Other Parameters
        ----------
            zernike_to_be_corrected: numpy array
                        None is equal to np.array([0,1,2,3,4])
                        for tip, tilt, fuoco, coma, coma
            commandId: numpy array
                    array containing the number of degrees of freedom 
                    to be commanded
            subapOffsets: tuple containing the subaperture offset in 
                    meters (axis from the center) and 
                    degrees (azimuthal angle) according Negro94

        Note: Zernike idx can be [0,1], [0,1,3,4], [0,1,2,3,4]
              and are coupled with dof [3,4], [1,2,3,4], [0,1,2,3,4] (??)

        Returns
        -------
                cmd: numpy array
                    final delta command for the optical alignment
                dove: string
                    file path containing measurements
        """
        par_position = self._ott.parabola.getPosition()
        rm_position = self._ott.referenceMirror.getPosition()
        # m4_position = self._ott.m4Exapode.getPosition()
        self._logger.info(
            "Calculation  of the alignment command using calibration "
            "measurements in tt = %s",
            self.tt_cal,
        )
        self._intMatModesVector = zernike_to_be_corrected
        self._commandId = dof_command_id
        tout = self.selectModesInIntMatAndRecConstruction(
            zernike_to_be_corrected, dof_command_id, subapOffsets=subapOffsets
        )
        self._intMat, self._rec, self._cmat = tout

        img = self._interf.acquire_phasemap(n_images, delay)
        # modRB20231027 to implement fullframe
        img = self._interf.intoFullFrame(img)
        # par = self._load_registeredPar(tnpar)
        # img = img -2*par
        # print('Par Removed')
        # endmod
        name = "StartImage.fits"
        calfilename = "CalibrationTracknum.txt"
        self.tt_al = Timestamp.now()
        dove = os.path.join(
            self._storageFolder(), self.tt_al
        )  # old folder naming convention, with CalibratioTN--ResultTN.
        os.makedirs(dove)
        # modRB 20231024 to save in text file the calibration tracknum
        calfile = open(dove + "/" + self.tt_cal, "w")
        calfile.write("")
        calfile.close()
        calfile = open(dove + "/" + calfilename, "w")
        calfile.write(self.tt_cal)
        calfile.close()
        # ott_status.save_positions(dove, self._ott) #saving the ott status
        # end of modRB
        self._interf.save_phasemap(dove, name, img)

        if self._who == "PAR + RM":
            cmd, self._zernikeVectorSelected, total_zernike_vector = (
                self._commandGenerator(img, tnPar)
            )
            self.par_command, self.rm_command = self._reorgCmdForParAndRm(
                cmd, dof_command_id
            )
            self._saveData(dove, par_position, rm_position)
            self._alignmentLog(total_zernike_vector, self.tt_al)
            self._loggerRuna.info("Calibration tt used = %s", self.tt_cal)
            self._loggerRuna.info(
                "Zernike calculate on image before alignment =  %s",
                str(total_zernike_vector),
            )
            self._loggerRuna.info(
                "Tracking number for alignment measurements = %s", self.tt_al
            )
            ott_status.save_positions(dove, self._ott)  # saving the ott status
            return self.par_command, self.rm_command, dove
        elif self._who == "M4":
            pass

            # cmd, zernike_vector = self._commandGenerator(img)
            # m4_command, zernike_vector_selected, total_zernike_vector =
            #  self._reorgCmdForM4(cmd)
            # self._saveAllDataM4(dove, m4_position, m4_command)
            # self._saveZernikeVector(dove, zernike_vector_selected)
            # return m4_command, dove

    def _alignmentLog(self, start_total_coef, tt):
        """The Log prints the Calibration tracknum and the initial Zernike 
        coeff [m]
        Prints also the alignment tracknum (i.e. the folder to be saved).
        Both is TRUE or FALSE are passed
        """
        fits_file_name = os.path.join(self._storageFolder(), 
                                      "AlignmentLog.txt")
        file = open(fits_file_name, "a+")
        file.write("Calib. Trackn & IniZern:  %s " % self.tt_cal)
        for i in range(start_total_coef.size):
            file.write("%7.3e " % start_total_coef[i])  # was 9.3
        file.write("\n")
        file.write("Result Trackn & EndZern:  %s " % tt)
        #         for i in range(total_coef.size):
        #             file.write('%9.3e ' %total_coef[i])
        #         file.write('\n')
        #         file.write('%s \n ************\n' %commandId)
        file.close()

    def selectModesInIntMatAndRecConstruction(
        self, zernike2control=None, commandId=None, subapOffsets=None
    ):
        """
        Other Parameters
        ----------
        zernike2control: numpy array
                    None is equal to np.array([0,1,2,3,4])
                    for tip, tilt, fuoco, coma, coma
        commandId: numpy array
                array containing the number of degrees of freedom 
                to be commanded

        nota: gli zernike possono essere [0,1], [0,1,3,4], [0,1,2,3,4]
             e vanno in coppia con i dof [3,4], [1,2,3,4], [0,1,2,3,4]

        if use_xmp_style is True, the function returns the interaction matrix 
        and the command matrix according Negro84 Zernike weighting
        """
        if subapOffsets is not None:
            subapOffAxisRadiusInMeter = subapOffsets[0]
            subapOffAxisAzimuthInDegrees = subapOffsets[1]
            if zernike2control is None:
                zernike2control = np.array([0, 1, 2, 3, 4])
            intMat = self.cal.getfullLocalInteractionMatrix()
            cmat = self.cal.getFullCommandMatrix()
            # Simulate M4 case sub-pupil
            pupilRadiusInMeter = 0.7  # to be dynamically set
            subapsRadiusInMeter = 0.3  # to be dynamically set
            # parpos = self._ott.parabola.getPosition()
            # rmpos  = (self._ott.referenceMirror.getPosition()-600)/1e3
            # subapOffAxisRadiusInMeter = 0.15 #to be dynamically set
            # subapOffAxisAzimuthInDegrees = 0 #to be dynamically set
            Zproj = ZernikeProjectionOnSubaperture(
                pupilRadiusInMeter,
                subapsRadiusInMeter,
                subapOffAxisRadiusInMeter,
                subapOffAxisAzimuthInDegrees,
            )
            # print("Negro84 sensitivity matrix small pupil coeffs (rows) and
            # large pupil coeffs (columns)")
            S = Zproj.get_projection_matrix()
            # print("Negro84 sensitivity matrix large pupil coeffs (rows) and
            # small pupil coeffs (columns)")
            S_inv = np.linalg.inv(S)
            SMat = (S_inv[:, zernike2control].T)[:, zernike2control].T
            W = SMat.T @ SMat
            new_intMat = intMat
            new_cmat = cmat
            rec = cmat @ np.linalg.inv(new_intMat.T @ W @ new_intMat
                                       ) @ new_intMat.T @ W

        else:
            intMat = self.cal.getInteractionMatrix()
            cmat = self.cal.getCommandMatrix()

            if zernike2control is None:
                new_intMat = intMat
                new_cmat = cmat
            else:
                new_intMat = intMat[zernike2control, :]

            if commandId is not None:
                new_intMat = new_intMat[:, commandId]
                new_cmat = cmat[commandId, :]
                new_cmat = new_cmat[:, commandId]

            rec = np.linalg.pinv(new_intMat)
            # self._plotIntMat()
        return new_intMat, rec, new_cmat

    def _plotIntMat(self):
        """ """
        # y = ['PAR_PIST', 'PAR_TIP', 'PAR_TILT', 'RM_TIP', 'RM_TILT']
        plt.clf()
        plt.imshow(self._intMat, origin="lower")
        plt.colorbar()
        plt.xlabel("Commands")
        plt.ylabel("Zernike Modes")
        return

    def _reorgCmdForParAndRm(self, cmd, commandId=None):
        """reorganizes the delta command in the
        right positions for par and rm"""
        dofIndex = np.append(OttParameters.PARABOLA_DOF, OttParameters.RM_DOF)
        par_command = np.zeros(6)
        rm_command = np.zeros(6)

        if commandId is not None:
            mycomm = np.zeros(5)
            mycomm[commandId] = cmd
            cmd = mycomm

        for i in range(cmd.size):
            if i < OttParameters.PARABOLA_DOF.size:
                par_command[dofIndex[i]] = cmd[i]
            else:
                rm_command[dofIndex[i]] = cmd[i]

        return par_command, rm_command

    def getReorganizatedCommandForParAndRm(self, cmd, commandId=None):
        """
        Function that returns the command reorganized according to
        the selected degrees of freedom
        """
        par_command, rm_command = self._reorgCmdForParAndRm(cmd, commandId)
        return par_command, rm_command

    def _commandGenerator(self, img, tnPar):
        """
        args:
            img: numpy array
                image

        returns:
            cmd: numpy array
                command not ready to apply
            zernike_vector_selected: numpy array
                vector of zernike to be corrected (selected by zernike2control)
            total_zernike_vector: numpy array
                all five zernike measured before the selection by 
                total_zernike_vector
        """
        tout = self._zernikeCoeffCalculator(
            img, tnPar
        )
        total_zernike_vector, zernike_vector_selected = tout
        print("zernike:")
        print(zernike_vector_selected)
        M = np.dot(self._cmat, self._rec)  # non serve la trasposta
        cmd = -np.dot(M, zernike_vector_selected)  # serve il meno
        print("mix command:")
        print(cmd)
        return cmd, zernike_vector_selected, total_zernike_vector

    def _zernikeCoeffCalculator(self, img, tnPar):
        """
        Returns:
                final_coef = zernike coeff on the image
                            (zernike modes 2,3,4,7,8)
                final_coef_selected = zernike selected using intMatModesVector 
                (zernike2control)
        """
        if fold_name.simulated_interf is True:
            mask_index = OtherParameters.MASK_INDEX_SIMULATORE
        else:
            mask_index = OtherParameters.MASK_INDEX_TOWER
        rois = roi.roiGenerator(img)
        mask = rois[mask_index]
        mm = np.ma.mask_or(img.mask, mask)

        new_image = np.ma.masked_array(img, mask=mm)
        # coef, mat = zernike.zernikeFit(new_image, np.arange(10)+1)
        # modRB20231026 to implement aux mask fitting. the following lines 
        # replace the previous one
        # from m4.utils import image_registration_lib as imgreg
        # from m4.analyzers import timehistory as th
        # img = th.frame2ottFrame(new_image,[580,20])
        if tnPar is not None:  #modRB 20240719 to implement the tnPAR parameter
            #tnpar = "20240521_161525"  
            par = self._load_registeredPar(tnPar)
            img = new_image - 2 * par
            print("Using global modes fitting, TNPar: " + tnPar)

            # par = imgreg.load_registeredPar(tnpar)
            cir = geo.qpupil(-1 * par.mask + 1)
            mm = geo.draw_mask(par.data * 0, cir[0], 
                               cir[1], 1.44 / 0.00076 / 2, out=0)
            # img = img-2*par !!!! removing the PAR twice!!noooooooo
            print("Removing the PAR shape")
            coef, mat = zernike.zernikeFitAuxmask(img, mm, np.arange(10) + 1)
        else:
            coef, mat = zernike.zernikeFit(img, np.arange(10) + 1)
        # end of modRB
        z = np.array([1, 2, 3, 6, 7])
        all_final_coef = coef[z]

        if self._intMatModesVector is None:
            final_coef_selected = all_final_coef
        else:
            final_coef_selected = np.zeros(self._intMatModesVector.size)
            for i in range(self._intMatModesVector.size):
                final_coef_selected[i] = all_final_coef[
                    self._intMatModesVector[i]]
        return all_final_coef, final_coef_selected

    def getZernikeWhitAlignerObjectOptions(self, img, tnPar):
        """
        Returns
        -------
        all_final_coef: numpy array
            zernike five coefficients on the image (zernike modes 2,3,4,7,8)
        final_coef_selected: numpy array
            zernike selected using intMatModesVector (zernike2control)
        """
        all_final_coef, final_coef_selected = self._zernikeCoeffCalculator(img, tnPar)
        return all_final_coef, final_coef_selected

    def _saveData(self, dove, par_position, rm_position):
        fits_file_name = os.path.join(dove, "AlignmentInfo.fits")
        header = pyfits.Header()
        header["WHO"] = self._who
        pyfits.writeto(fits_file_name, self._intMat, header)
        pyfits.append(fits_file_name, self._rec, header)
        pyfits.append(fits_file_name, self._cmat, header)
        if self._intMatModesVector is None:
            self._intMatModesVector = np.array([0, 1, 2, 3, 4])
        pyfits.append(fits_file_name, self._intMatModesVector, header)
        if self._commandId is None:
            self._commandId = np.array([0, 1, 2, 3, 4])
        pyfits.append(fits_file_name, self._commandId, header)
        pyfits.append(fits_file_name, self._zernikeVectorSelected, header)

        # salvataggi espliciti
        name = "PositionAndDeltaCommand.fits"
        vector = np.array(
            [par_position, rm_position, self.par_command, self.rm_command]
        )
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, vector)
        if self._intMatModesVector is not None:
            fits_file_name = os.path.join(dove, "intMatModesVector.fits")
            pyfits.writeto(fits_file_name, self._intMatModesVector)
        if self._commandId is not None:
            fits_file_name = os.path.join(dove, "commandId.fits")
            pyfits.writeto(fits_file_name, self._commandId)
        name = "Zernike.fits"
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, self._zernikeVectorSelected)

    @staticmethod
    def loadAlignmentObjectFromFits(tt):
        """Creates the object using information contained in alignment fits 
            file

        Parameters
        ----------
        tt: string
            tracking number

        Returns
        -------
        theObject: ibjecct
                 opt_alignment class object
        """
        ott = None
        interf = None
        tt_cal = tt.split("--")[0]
        theObject = OpticalAlignment(tt_cal, ott, interf)
        theObject.tt_cal = tt_cal
        theObject.tt_al = tt.split("--")[1]
        dove = os.path.join(theObject._storageFolder(), tt)
        file = os.path.join(dove, "AlignmentInfo.fits")
        header = pyfits.getheader(file)
        hduList = pyfits.open(file)
        theObject._who = header["WHO"]
        theObject._intMat = hduList[0].data
        theObject._rec = hduList[1].data
        theObject._cmat = hduList[2].data
        theObject._intMatModesVector = hduList[3].data
        theObject._commandId = hduList[4].data
        theObject._zernikeVectorSelected = hduList[5].data
        return theObject

    def _load_registeredPar(self, tn, zlist=[1, 2, 3, 4]):
        # fold=th.foldname.PARABOLA_REMAPPED_FOLDER+'/'+tn+'/'
        name = fold_name.PARABOLA_REMAPPED_FOLDER +\
                  "/" + tn + "/" + "par_remapped.fits"
        print("Loading registered Par " + name)
        hdu = pyfits.open(name)
        img = hdu[0].data
        mask = hdu[1].data
        imgout = np.ma.masked_array(img, mask)
        coeff, mat = zernike.zernikeFit(imgout, zlist)
        surf = zernike.zernikeSurface(imgout, coeff, mat)
        imgout = imgout - surf
        return imgout

    # M4 ###
    # va riscritta
    def _reorgCmdForM4(self, cmd):
        dofIndex = OttParameters.M4_DOF
        m4_command = np.zeros(6)
        for i in range(cmd.size):
            m4_command[dofIndex[i]] = cmd[i]
        return m4_command

    def _saveAllDataM4(self, dove, m4_position, m4_command):
        name = "PositionAndDeltaCommand.fits"
        vector = np.array([m4_position, m4_command])
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, vector)

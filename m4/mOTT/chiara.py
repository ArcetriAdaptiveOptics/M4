'''
file temporaneo
'''

import os
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration import start
from m4.configuration.ott_parameters import OpcUaParameters
from m4.alignment import Alignment
from m4.ground.opc_ua_controller import OpcUaController

class Measure_mOTT():

    def __init__(self):
        """The constructor """
        self._ott, self._interf  = start.create_ott()
        self._opc = OpcUaController()
        self._a = Alignment(self._ott)

    def deltaPositionVsSlide(self, n_object_to_move, shift, tt_for_align):
        n_frames_alignment = 10
        par0 = self._ott.parab()
        rm0 = self._ott.refflat()
        delta_par = []
        delta_rm = []
        delta_object = []

        if n_object_to_move == 0:
            n_iter = np.int((OpcUaParameters.max_angle-OpcUaParameters.min_angle)/shift)
            angle0 = self._ott.angle()
            for i in range(n_iter):
                angle = self._ott.angle(angle0+((i+1)*shift))
                par_cmd, rm_cmd = self._a.ott_alignment(n_frames_alignment, 1,
                                                        np.array([0,1]), np.array([3, 4]),
                                                        tt_for_align)
                par = self._ott.parab()
                rm = self._ott.refflat()
                delta_par.append(par - par0)
                delta_rm.append(rm - rm0)
                delta_object.append(angle - angle0)
                self._saveData(delta_par, delta_rm, delta_object, n_object_to_move)
        elif n_object_to_move == 1:
            n_iter = np.int((OpcUaParameters.max_r_slide-OpcUaParameters.min_r_slide)/shift)
            r_slide0 = self._ott.rslide()
            for i in range(n_iter):
                rslide = self._ott.rslide(r_slide0+((i+1)*shift))
                par_cmd, rm_cmd = self._a.ott_alignment(n_frames_alignment, 1,
                                                        np.array([0,1]), np.array([3, 4]),
                                                        tt_for_align)
                par = self._ott.parab()
                rm = self._ott.refflat()
                delta_par.append(par - par0)
                delta_rm.append(rm - rm0)
                delta_object.append(rslide - r_slide0)
                self._saveData(delta_par, delta_rm, delta_object, n_object_to_move)
                self._saveData(delta_par, delta_rm, delta_object, n_object_to_move)
        elif n_object_to_move == 2:
            n_iter = np.int((OpcUaParameters.max_slide-OpcUaParameters.min_slide)/shift)
            slide0 = self._ott.slide()
            for i in range(n_iter):
                slide = self._ott.slide(slide0+((i+1)*shift))
                par_cmd, rm_cmd = self._a.ott_alignment(n_frames_alignment, 1,
                                                        np.array([0,1]), np.array([3, 4]),
                                                        tt_for_align)
                par = self._ott.parab()
                rm = self._ott.refflat()
                delta_par.append(par - par0)
                delta_rm.append(rm - rm0)
                delta_object.append(slide - slide0)
                self._saveData(delta_par, delta_rm, delta_object, n_object_to_move)
        else:
            raise Exception('Incorrect number of object to move')
        return np.array(delta_object), np.array(delta_par), np.array(delta_rm)

    def _saveData(self, delta_par, delta_rm, delta_object, n_object_to_move):
        dove = '?'
        if n_object_to_move == 0:
            name = 'delta_angle.fits'
        elif n_object_to_move == 1:
            name = 'delta_rslide.fits'
        elif n_object_to_move == 2:
            name = 'delta_slide.fits'

        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, np.array(delta_object), overwrite=True)
        fits_file_name = os.path.join(dove, 'delta_PAR_positions.fits')
        pyfits.writeto(fits_file_name, np.array(delta_par), overwrite=True)
        fits_file_name = os.path.join(dove, 'delta_RM_positions.fits')
        pyfits.writeto(fits_file_name, np.array(delta_rm), overwrite=True)
        
        
        
        
        
        
ttp = '20210111_152430'
from m4.utils.optical_calibration import OpticalCalibration
cal = OpticalCalibration.loadCommandMatrixFromFits(ttp)
cal.createCube(ttp, 77)
cube = cal.getCube()
ima = cube[:,:,0]
from m4.utils.roi import ROI
r = ROI()
roi = r.roiGenerator(ima)
mask = roi[1]
intMat = cal.getInteractionMatrix(mask)

cmd0 = np.array([0,0,0,1,1])
z0 = np.dot(intMat, cmd0)/3
intMatModesVector = np.array([0,1])
commandId = np.array([3,4])

intMatModesVector = np.array([0,1,2,3,4])
commandId = np.array([0,1,2,3,4])
z_new = z0[intMatModesVector]

#standard
new_intMat = intMat[intMatModesVector, :]
new_intMat = new_intMat[:,commandId]
new_rec = np.linalg.pinv(new_intMat)
z_new = z0[intMatModesVector]
cmd_new = np.zeros(5)
cmd_new[commandId] = np.dot(new_rec, z_new)

#with cmat
cmat = np.diag(np.array([0.8, 1.5, 1.5, 3. , 3. ]))
new_intMat2 = intMat[intMatModesVector, :]
cmat2 = cmat[commandId, :]
new_rec2 = np.linalg.pinv(new_intMat2)
M = np.dot(cmat2, new_rec2)
cmd_new2 = np.zeros(5)
cmd_new2[commandId] = np.dot(M, z_new)


cmat = np.diag(np.array([0.8, 1.5, 1.5, 3. , 3. ]))
pp = intMat[intMatModesVector, :]
new_intMat3 = pp[:, commandId]
cc = cmat[commandId, :]
cmat3 = cc[:, commandId]
new_rec3 = np.linalg.pinv(new_intMat3)
M = np.dot(cmat3, new_rec3)
cmd_new3 = np.zeros(5)
cmd_new3[commandId] = np.dot(M, z_new)





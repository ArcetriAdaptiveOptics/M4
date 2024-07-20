#4D Configurations
phasecam_markerconfig =  'D:/config/20240608_negativeMarkersMask50mm.4Dini'
phasecam_alignmentconfig = 'D:/config/alignmentConf-noMask-dataFill-lowModThr.4Dini'
phasecam_noiseconfig = 'D:/config/20240608_negativeMarkersMask50mm.4Dini'
phasecam_monitorconfig = 'D:/config/20240608_negativeMarkersMask50mm.4Dini'
phasecam_baseconfig = 'D:/config/20240719_2000x2000_RMCenter.4Dini'

#Parabola registration
remappedpar_tn = '20240521_161525'


#alignment
alignmentCalibration_tn = '20240521_120211' #calibrated with par subtraction

#alignment Calibration
alignCal_parPist = 0.7
alignCal_parTip  = 100
alignCal_parTilt = 100
alignCal_rmTip   = 6
alignCal_rmTilt  = 6

#noise acquisition parameters
import numpy as np
noise_nframes = 2000
noise_tau_vector = np.arange(1,100,2)
noise_difftemplate = np.array([3,11,25,37,51])

#items translation
refMirror_maxstepBeforeAlignment = 100
truss_maxstepBeforeAlignment = 100

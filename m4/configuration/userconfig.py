# 4D Configurations
phasecam_markerconfig    =     "D:/config/20240608_negativeMarkersMask50mm.4Dini"
# phasecam_alignmentconfig = "D:/config/alignmentConf-noMask-dataFill-lowModThr.4Dini"
phasecam_OTTalignmentconfig =   "D:/config/20251106_RMmask-Alignment.4Dini" 
phasecam_M4alignmentconfig  =   "D:/config/20260312_203000-segm0-alignment.4Dini"
phasecam_dmconfig        =       "D:/config/20260317_151400-2000X2000-splremoved-OPTIMIZED-minSlope3x3.4Dini"  #20250910_AllSegments-nospl.4Dini"
phasecam_segmentconfig  =       "D:/config/20260317_151400-2000X2000-splremoved-OPTIMIZED-minSlope3x3.4Dini"
phasecam_noiseconfig     = "D:/config/20240608_negativeMarkersMask50mm.4Dini"
phasecam_monitorconfig   = "D:/config/20251001_monitoring.4Dini"
phasecam_baseconfig      = "D:/config/20250909_RMAlignment.4Dini"

# Parabola registration
remappedpar_tn = "20240521_161525"


# alignment
alignmentCalibration_tn = "20260207_192654"


# alignment Calibration
alignCal_parPist = 0.7
alignCal_parTip = 100
alignCal_parTilt = 100
alignCal_rmTip = 6
alignCal_rmTilt = 6

# noise acquisition parameters
import numpy as np

noise_nframes = 2000
noise_tau_vector = np.arange(1, 100, 2)
noise_difftemplate = np.array([3, 11, 25, 37, 51])

# items translation
refMirror_maxstepBeforeAlignment = 100
truss_maxstepBeforeAlignment = 100

# RM and PAR sliders configuration
parslider4segment = 0.757  # corresponds to current configuration with DP
rmslider4alignment = 0.4
rmsliderout = -0.1


### Configurations for the DM
dm_defaultFlatCmd = '20260316_204445'
dm_defaultIFF     = '20260318_225155' #calibrated with optimized mask #'20260313_102030'

### Configuration for IFF sampling
iff_zonal_amp     = 150e-9
iff_modal_ampTN   = '20260218_220911'


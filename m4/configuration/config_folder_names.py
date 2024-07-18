"""
Author(s)
    - Chiara Selmi: written in 2021

Description
-----------
This is a library module, containing the information about the devices, whether
they are simulated or not, and all the important filepaths stored in variables.

How to Use it
-------------
This module alone is useless. First of all, must be populated through the module
''m4.configuration.update_folder_paths'', then the variables are accessible directly
from the update_folder_paths module, or either way the can be accessed importing
this module.

Example
-------
>>> from m4.configuration import update_folder_paths
>>> from m4.configuration import config_folder_names as fn
>>> fn.BASE_PATH
'.../M4/m4/data'
"""
mirror_conf = '20170430'
optical_conf = '20150730'
# Devices settings
simulated_interf            = None
simulated_dm                = None
simulated_accelerometers    = None
simulated_angleRotator      = None
simulated_m4Exapode         = None
simulated_parSlider         = None
simulated_par               = None
simulated_rmSlider          = None
simulated_rm                = None
simulated_tempSensors       = None
# Base Paths
BASE_PATH                           = None
CONFIGURATION_ROOT_FOLDER           = None
ALL_CALIBRATION_DATA_ROOT_FOLDER    = None
OPT_DATA_FOLDER                     = None
OUT_FOLDER                          = None
MIRROR_FOLDER                       = None
OPTICAL_FOLDER                      = None
# Built paths on base paths
OTT_CALIB_CONF_FOLDER               = None
OPD_IMAGES_ROOT_FOLDER              = None
PHASECAM_ROOT_FOLDER                = None
LOG_ROOT_FOLDER                     = None
IFFUNCTIONS_ROOT_FOLDER             = None
FLAT_ROOT_FOLD                      = None
CALIBRATION_ROOT_FOLDER             = None
ALIGNMENT_ROOT_FOLDER               = None
ZERNIKECOMMANDTEST_ROOT_FOLDER      = None
NOISE_ROOT_FOLDER                   = None
SPL_ROOT_FOLDER                     = None
CALIBALL_ROOT_FOLDER                = None
MODESVECTOR_ROOT_FOLDER             = None
MODALBASE_ROOT_FOLDER               = None
MODALAMPLITUDE_ROOT_FOLDER          = None
COMMANDHISTORY_ROOT_FOLDER          = None
GEOTRANSFORM_ROOT_FOLDER            = None
ROT_OPT_ALIGN_ROOT_FOLDER           = None
PT_ROOT_FOLDER                      = None
OPD_SERIES_ROOT_FOLDER              = None
REPEATABILITY_ROOT_FOLDER           = None
PISTON_TEST_ROOT_FOLDER             = None
MAPPING_TEST_ROOT_FOLDER            = None
ACC_ROOT_FOLDER                     = None
SIMUL_DATA_CALIB_DM_FOLDER          = None
PARABOLA_CGH_FOLDER                 = None
PARABOLA_REMAPPED_FOLDER            = None
INTMAT_ROOT_FOLDER                  = None
DM_CONFIGURATION_ID                 = None
MARKERS_ROOT_FOLDER                 = None
MONITORING_ROOT_FOLDER              = None

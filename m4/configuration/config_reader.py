"""
Authors
    - Chiara Selmi: written in 2021
    - Pietro Ferraiuolo: updated in 2024

Description
-----------
Module containing the Class for reading data from the software's configuration 
file.
"""
import os
import yaml

#    simulated = 0  # 1 per il simulatore
#    BASE_PATH = '/home/labot/data/M4/Data'

class configuration_path():
    """
    This class contains all the methods to read and write the configuration file
    for the M4's software, both for the devices and the data path tree. If a file
    path is not fount, it will create it, so that at the end, the data tree structure
    will be present on the machine.
    """
    def __init__(self, confFile):
        """ The constructor"""
        with open(confFile) as file:
            self._conf = yaml.load(file, Loader=yaml.FullLoader)

    @property
    def simulated_interf(self):
        return self._conf['simulated_interf']

    @property
    def simulated_dm(self):
        return self._conf['simulated_dm']

    @property
    def simulated_accelerometers(self):
        return self._conf['simulated_accelerometers']

    @property
    def simulated_angleRotator(self):
        return self._conf['simulated_angleRotator']

    @property
    def simulated_m4Exapode(self):
        return self._conf['simulated_m4Exapode']
    @property
    def simulated_parSlider(self):
        return self._conf['simulated_parSlider']

    @property
    def simulated_par(self):
        return self._conf['simulated_par']

    @property
    def simulated_rmSlider(self):
        return self._conf['simulated_rmSlider']

    @property
    def simulated_rm(self):
        return self._conf['simulated_rm']

    @property
    def simulated_tempSensors(self):
        return self._conf['simulated_tempSensors']

### PATH_NAMES ###
    @property
    def BASE_PATH(self):
        if 'base_path' in self._conf.keys():
            return self._conf['base_path']
        else:
            return '/mnt/m4storage/Data'

    @property
    def CONFIGURATION_ROOT_FOLDER(self):
        if 'configuration_root_folder' in self._conf.keys():
            return self._conf['configuration_root_folder']
        else:
            path = os.path.join(self.BASE_PATH, 'SYSCONFData')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def ALL_CALIBRATION_DATA_ROOT_FOLDER(self):
        if 'all_calibration_data_root_folder' in self._conf.keys():
            return self._conf['all_calibration_data_root_folder']
        else:
            path = os.path.join(self.BASE_PATH, 'M4Data')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def OPT_DATA_FOLDER(self):
        if 'opt_data_folder' in self._conf.keys():
            return self._conf['opt_data_folder']
        else:
            path = os.path.join(self.ALL_CALIBRATION_DATA_ROOT_FOLDER, 'OPTData')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def OUT_FOLDER(self):
        if 'out_folder' in self._conf.keys():
            return self._conf['out_folder']
        else:
            path = os.path.join(self.BASE_PATH, 'Results')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def MIRROR_FOLDER(self):
        if 'mirror_folder' in self._conf.keys():
            return self._conf['mirror_folder']
        else:
            path = os.path.join(self.BASE_PATH, 'MIRROR_System')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def OPTICAL_FOLDER(self):
        if 'optical_folder' in self._conf.keys():
            return self._conf['optical_folder']
        else:
            path = os.path.join(self.BASE_PATH, 'OPTICAL_System')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

### FOLD_NAME ###
    @property
    def OTT_CALIB_CONF_FOLDER(self):
        if 'ott_calib_conf_folder' in self._conf.keys():
            return self._conf['ott_calib_conf_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'OTTCalibConf')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def OPD_IMAGES_ROOT_FOLDER(self):
        if 'opd_images_root_folder' in self._conf.keys():
            return self._conf['opd_images_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'OPDImages')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def LOG_ROOT_FOLDER(self):
        if 'log_root_folder' in self._conf.keys():
            return self._conf['log_root_folder']
        else:
            path1 = os.path.join(self.BASE_PATH, 'LOGData')
            path2 = os.path.join(path1, 'mylog')
            if os.path.exists(path1) is False:
                os.mkdir(path1)
                os.mkdir(path2)
            return path2

    @property
    def PHASECAM_ROOT_FOLDER(self):
        if 'phasecam_root_folder' in self._conf.keys():
            return self._conf['phasecam_root_folder']
        else:
            path = '/home/m4/4d/Zcopy/'
            # if os.path.exists(path) is False:
            #     os.mkdir(path)
            return path

    @property
    def IFFUNCTIONS_ROOT_FOLDER(self):
        if 'iffunctions_root_folder' in self._conf.keys():
            return self._conf['iffunctions_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'IFFunctions')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def FLAT_ROOT_FOLD(self):
        if 'flat_root_folder' in self._conf.keys():
            return self._conf['flat_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'Flattening')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def CALIBRATION_ROOT_FOLDER(self):
        if 'calibration_root_folder' in self._conf.keys():
            return self._conf['calibration_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'AlignmentCalibration')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def ALIGNMENT_ROOT_FOLDER(self):
        if 'alignment_root_folder' in self._conf.keys():
            return self._conf['alignment_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'Alignment')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def ZERNIKECOMMANDTEST_ROOT_FOLDER(self):
        if 'zernikecommandtest_root_folder' in self._conf.keys():
            return self._conf['zernikecommandtest_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'ZernikeCommandTest')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path
        
    @property
    def NOISE_ROOT_FOLDER(self):
        if 'noise_root_folder' in self._conf.keys():
            return self._conf['noise_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'Noise')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def SPL_ROOT_FOLDER(self):
        if 'spl_root_folder' in self._conf.keys():
            return self._conf['spl_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'SPL')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def CALIBALL_ROOT_FOLDER(self):
        if 'caliball_root_folder' in self._conf.keys():
            return self._conf['caliball_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'Caliball')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def MODESVECTOR_ROOT_FOLDER(self):
        if 'modesvector_root_folder' in self._conf.keys():
            return self._conf['modesvector_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'ModesVector')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def MODALBASE_ROOT_FOLDER(self):
        if 'modalbase_root_folder' in self._conf.keys():
            return self._conf['modalbase_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'ModalBase')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def MODALAMPLITUDE_ROOT_FOLDER(self):
        if 'modalamplitude_root_folder' in self._conf.keys():
            return self._conf['modalamplitude_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'ModalAmplitude')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def COMMANDHISTORY_ROOT_FOLDER(self):
        if 'commandhistory_root_folder' in self._conf.keys():
            return self._conf['commandhistory_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'CommandHistory')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def GEOTRANSFORM_ROOT_FOLDER(self):
        if 'geotransform_root_folder' in self._conf.keys():
            return self._conf['geotransform_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'GeomTransf')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def ROT_OPT_ALIGN_ROOT_FOLDER(self):
        if 'rot_opt_align_root_folder' in self._conf.keys():
            return self._conf['rot_opt_align_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'RotOptAlignment')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def PT_ROOT_FOLDER(self):
        if 'pt_root_folder' in self._conf.keys():
            return self._conf['pt_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'PTCalibration')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def OPD_SERIES_ROOT_FOLDER(self):
        if 'opd_series_root_folder' in self._conf.keys():
            return self._conf['opd_series_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'OPDSeries')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def REPEATABILITY_ROOT_FOLDER(self):
        if 'repeatability_root_folder' in self._conf.keys():
            return self._conf['repeatability_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'Repeatability')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def PISTON_TEST_ROOT_FOLDER(self):
        if 'piston_test_root_folder' in self._conf.keys():
            return self._conf['piston_test_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'PistonTest')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def MAPPING_TEST_ROOT_FOLDER(self):
        if 'mapping_test_root_folder' in self._conf.keys():
            return self._conf['mapping_test_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'Mapping')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def ACC_ROOT_FOLDER(self):
        if 'acc_root_folder' in self._conf.keys():
            return self._conf['acc_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'AccelerometersData')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

#     @property
#     def POINTER_ALIGN_ROOT_FOLDER(self):
#         if 'pointer_align_root_folder' in self._conf.keys():
#             return self._conf['pointer_align_root_folder']
#         else:
#             return os.path.join(self.OPT_DATA_FOLDER, 'PointerAlign')

    @property
    def SIMUL_DATA_CALIB_DM_FOLDER(self):
        if 'simul_data_calib_dm_folder' in self._conf.keys():
            return self._conf['simul_data_calib_dm_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'SimDataCalibDM')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def PARABOLA_CGH_FOLDER(self):
        if 'parabola_cgh_measurements_folder' in self._conf.keys():
            return self._conf['parabola_cgh_measurements_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'ParabolaCGHMeasurements')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def PARABOLA_REMAPPED_FOLDER(self):
        if 'parabola_remapped_folder' in self._conf.keys():
            return self._conf['parabola_remapped_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'ParabolaRemapped')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path
        
    @property
    def INTMAT_ROOT_FOLDER(self):
        if 'intmat_root_folder' in self._conf.keys():
            return self._conf['intmat_root_folder']
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'INTMatrices')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def DM_CONFIGURATION_ID(self):
        if 'dm_configuration' in self._conf.keys():
            return str(self._conf['dm_configuration'])
        else:
            path = os.path.join(self.MIRROR_FOLDER, '')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    def MONITORING_ROOT_FOLDER(self):
        if 'monitoring_root_folder' in self._conf.keys():
            return str(self._conf['monitoring_root_folder'])
        else:
            path = os.path.join(self.OPT_DATA_FOLDER, 'MonitoringData')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path

    @property
    def MARKERS_ROOT_FOLDER(self):
        if 'markers_root_folder' in self._conf.keys():
            return str(self._conf['markers_root_folder'])
        else:
            path = os.path.join(self.OPD_DATA_FOLDER, 'Markers')
            if os.path.exists(path) is False:
                os.mkdir(path)
            return path        

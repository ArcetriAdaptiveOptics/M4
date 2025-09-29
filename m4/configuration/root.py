from opticalib.core import root as _oroot
from os.path import join as _join


class Folders(_oroot._folds):

    def __init__(self):
        super().__init__()
        self.ACCELEROMETERS_ROOT_FOLDER = _join(
            self.BASE_DATA_PATH, "AccelerometersData"
        )
        self.CALIBALL_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "Caliball")
        self.CALIBRATION_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "Calibration")
        self.COMMAND_HISTORY_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "CommandHistory")
        self.GEOTRANSF_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "GeomTransf")
        self.MAPPING_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "Mapping")
        self.MARKERS_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "Markers")
        self.MODALAMP_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "ModalAmplitude")
        self.MODALBASE_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "ModalBase")
        self.MODEVEC_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "ModesVector")
        self.NOISE_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "Noise")
        self.OTTCALIB_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "OTTCalibConf")
        self.PARCGH_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "ParabolaCGHMeasurements")
        self.PARCALIB_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "ParabolaRemapped")
        self.PARPOS_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "PARPosition")
        self.SPL_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "SPL")
        self.MIRROR_FOLDER = _join(self.BASE_DATA_PATH, "MIRROR_System")
        self.OPTICAL_FOLDER = _join(self.BASE_DATA_PATH, "OPTICAL_System")
        self.mirror_conf = "20170430"
        self.optical_conf = "20150730"
        self.SYSDATA_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "SYSCONFData")
        self.SIMDATACALIB_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "SimDataCalibDM")
        self.MONITORING_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "Monitoring")
        # self.PARTEST_ROOT_FOLDER = _join(self.BASE_DATA_PATH, 'PARTest')
        # self.PISTONTEST_ROOT_FOLDER = _join(self.BASE_DATA_PATH, 'PistonTest')
        # self.PTCALIB_ROOT_FOLDER = _join(self.BASE_DATA_PATH, 'PTCalibration')
        # self.REFMIR_ROOT_FOLDER = _join(self.BASE_DATA_PATH, 'RefMirror')
        # self.REPEAT_ROOT_FOLDER = _join(self.BASE_DATA_PATH, 'Repeatability')
        # self.ROTOPTALIGN_ROOT_FOLDER = _join(self.BASE_DATA_PATH, 'RotOptAlignment')
        # self._ROOT_FOLDER = _join(self.BASE_DATA_PATH, 'Tilt_test')
        # self._ROOT_FOLDER = _join(self.BASE_DATA_PATH, 'ZernikeCommandTest')


folders = Folders()

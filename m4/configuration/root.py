from opticalib.core import root as _oroot
from os.path import join as _join


class Folders(_oroot._folds):

    def __init__(self):
        super().__init__()
        self.ACCELEROMETERS_ROOT_FOLDER = _join(
            self.BASE_DATA_PATH, "AccelerometersData"
        )
        self.CALIBALL_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "Caliball")
        self.CALIBRATION_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "Calibration")
        self.COMMAND_HISTORY_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "CommandHistory")
        self.GEOTRANSF_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "GeomTransf")
        self.MAPPING_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "Mapping")
        self.MARKERS_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "Markers")
        self.MODALAMP_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "ModalAmplitude")
        self.MODALBASE_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "ModalBase")
        self.MODEVEC_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "ModesVector")
        self.NOISE_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "Noise")
        self.OTTCALIB_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "OTTCalibConf")
        self.PARCGH_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "ParabolaCGHMeasurements")
        self.PARCALIB_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "ParabolaRemapped")
        self.PARPOS_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "PARPosition")
        self.SPL_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "SPL")
        self.SPL_DATA_ROOT_FOLDER = _join(self.SPL_ROOT_FOLDER, "Data")
        self.SPL_RESULTS_ROOT_FOLDER = _join(self.SPL_ROOT_FOLDER, "Results")
        self.SPL_FRINGES_ROOT_FOLDER = _join(self.SPL_ROOT_FOLDER, "Fringes")
        self.MIRROR_FOLDER = _join(self.BASE_DATA_PATH, "MIRROR_System")
        self.OPTICAL_FOLDER = _join(self.BASE_DATA_PATH, "OPTICAL_System")
        self.mirror_conf = "20170430"
        self.optical_conf = "20150730"
        self.SYSDATA_ROOT_FOLDER = _join(self.BASE_DATA_PATH, "SYSCONFData")
        self.SIMDATACALIB_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, "SimDataCalibDM")
        # self.PARTEST_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, 'PARTest')
        # self.PISTONTEST_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, 'PistonTest')
        # self.PTCALIB_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, 'PTCalibration')
        # self.REFMIR_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, 'RefMirror')
        # self.REPEAT_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, 'Repeatability')
        # self.ROTOPTALIGN_ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, 'RotOptAlignment')
        # self._ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, 'Tilt_test')
        # self._ROOT_FOLDER = _join(self.OPT_DATA_ROOT_FOLDER, 'ZernikeCommandTest')


folders = Folders()

from m4.configuration.config import *
from m4.configuration.ott_parameters import *
from m4.ground import object_from_fits_file_name as obj

class Mirror():
    def __init__(self):
        """The constructor """
        curr_conffolder = os.path.join(path_name.CONFIGURATION_ROOT_FOLDER, tnconf_mirror)
        self.vmat = obj.readFits_object(os.path.join(curr_conffolder,'vmat.fits'))
        self.ff = obj.readFits_object(os.path.join(curr_conffolder,'ff_matrix.fits'))

        self.m4od = OttParameters.m4od
        self.m4optod = OttParameters.m4optod
        self.m4id = OttParameters.m4id

class Parabola():
    def __init__(self):
        """The constructor """
        self.radius = OttParameters.parab_radius
        self.dof = OttParameters.PARABOLA_DOF

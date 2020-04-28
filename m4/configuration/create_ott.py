from m4.configuration.config import *
from m4.configuration.ott_parameters import *
from m4.ground import object_from_fits_file_name as obj

class Mirror():
    curr_conffolder = os.path.join(path_name.CONFIGURATION_ROOT_FOLDER, tnconf_mirror)
    vmat = obj.readFits_object(os.path.join(curr_conffolder,'vmat.fits'))
    ff = obj.readFits_object(os.path.join(curr_conffolder,'ff_matrix.fits'))

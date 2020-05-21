from m4.configuration.config import *
from m4.configuration.ott_parameters import *
from m4.ground import object_from_fits_file_name as obj
from m4.utils.roi import ROI

class OTT():
    def __init__(self):
        """The constructor """
        self._r = ROI()
        self._slide = 0
        self._rslide = 0
        self._angle = 0
        self.start_position = np.zeros(6)

# Elements position
    def slide(self, par_trans=None):
        ''' parabola translation: -0.9 m +0.9m
        '''
        if par_trans is None:
            self._slide = self._slide
        else:
            self._slide = par_trans
        return self._slide

    def rslide(self, ref_flat=None):
        ''' ref-flat: -0.05 m to 0.4 m
        '''
        if ref_flat is None:
            self._rslide = self._rslide
        else:
            self._rslide = ref_flat
        return self._rslide

    def angle(self, rot_ring_angle=None):
        ''' rotating ring angle: 0 to 360°
        '''
        if rot_ring_angle is None:
            self._angle = self._angle
        else:
            self._angle = rot_ring_angle
        return self._angle
# Elements alignment
    def parb(self, start_position=None):
        if start_position is None:
            parb = self.start_position
        else:
            parb = start_position
        return parb

    def refflat(self, start_position=None):
        if start_position is None:
            refflat = self.start_position
        else:
            refflat = start_position
        return refflat

    def m4(self, start_position=None):
        if start_position is None:
            m4 = self.start_position
        else:
            m4 = start_position
        return m4
### Sensitivity matrices
    def _readMatFromTxt(self, file_name):
        ''' 11Zernike x 6 displacements, m RMS, per 1 m displacement - or 1 radiant rotation
        '''
        file = open(file_name, 'r')
        triplets=file.read().split()
        x = np.array(triplets)
        mat = x.reshape(11, 6)
        return mat

    def zmx_parpos2z(self, file_name):
        file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_PAR_pos2z.txt'
        mat = self._readMatFromTxt(file_name)
        return mat

    def zmx_refflatpos2z(self, file_name):
        file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_FM_pos2z.txt'
        mat = self._readMatFromTxt(file_name)
        return mat

    def zmx_m4pos2z(self, file_name):
        file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_M4_pos2z.txt'
        mat = self._readMatFromTxt(file_name)
        return mat
# Zmat
    def zmat(self):
        smap = np.zeros((Interferometer.N_PIXEL[0], Interferometer.N_PIXEL[1]))
        mask_par    = self._r.create_circular_mask(Interferometer.N_PIXEL[0]/2,
                                               Interferometer.N_PIXEL[1]/2,
                                               OttParameters.parab_radius*OttParameters.PIXEL_SCALE-2)
        mask    = self._r.create_circular_mask(Interferometer.N_PIXEL[0]/2,
                                               Interferometer.N_PIXEL[1]/2,
                                               OttParameters.parab_radius*OttParameters.fold_radius)
        rmask   = self._r.create_circular_mask(OttParameters.rflat_radius*OttParameters.PIXEL_SCALE,
                                               OttParameters.rflat_radius*OttParameters.PIXEL_SCALE,
                                               OttParameters.rflat_radius*OttParameters.PIXEL_SCALE-1)

        return mask_par, mask, rmask

class DMmirror():
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

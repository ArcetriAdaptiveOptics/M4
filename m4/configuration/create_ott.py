'''
Autors
  - C. Selmi: written in 2020
'''

from m4.configuration.config import *
from m4.configuration.ott_parameters import *
from astropy.io import fits as pyfits
from m4.ground import object_from_fits_file_name as obj
from m4.utils.roi import ROI
from m4.ground.zernikeGenerator import ZernikeGenerator

class OTT():
    def __init__(self):
        """The constructor """
        self._r = ROI()
        self._zg = ZernikeGenerator(2*OttParameters.parab_radius*OttParameters.pscale)
        self._slide = 0
        self._rslide = 0
        self._angle = 0
        self.start_position = np.zeros(6)
        self.smap = np.zeros((Interferometer.N_PIXEL[0], Interferometer.N_PIXEL[1]))
        self.rmap = np.zeros(((2*OttParameters.rflat_radius*OttParameters.pscale).astype(int),
                              (2*OttParameters.rflat_radius*OttParameters.pscale).astype(int)))

# Elements position
    def slide(self, par_trans=None):
        ''' Function to set the parabola translation (range: -0.9 m +0.9 m)

        Other Parameters
        ----------
        par_trans: int, optional
                If par_trans is not set it's equal to zero

        Returns
        -------
            par_trans: int
                    parabola translation
        '''
        if par_trans is None:
            self._slide = self._slide
        else:
            self._slide = par_trans
        return self._slide

    def rslide(self, ref_flat=None):
        '''  Function to set the reference flat mirror (range: -0.05 m to 0.4 m)

        Other Parameters
        ----------
        ref_flat: int, optional
                If ref_flat is not set it's equal to zero

        Returns
        -------
        ref_flat: int
                reference flat mirror position
        '''
        if ref_flat is None:
            self._rslide = self._rslide
        else:
            self._rslide = ref_flat
        return self._rslide

    def angle(self, rot_ring_angle=None):
        ''' Function to set the rotating ring angle (range: 0 to 360°)

        Other Parameters
        ----------
            rot_ring_angle: int, optional
                If rot_ring_angle is not set it's equal to zero

        Returns
        -------
            rot_ring_angle: int
                            rot_ring_angle
        '''
        if rot_ring_angle is None:
            self._angle = self._angle
        else:
            self._angle = rot_ring_angle
        return self._angle
# Elements alignment
    def parab(self, start_position=None):
        '''Function to set the start position of the parable

        Other Parameters
        ----------
        start_position: numpy array, optional
                        vector of six position
                        If start_position is not set it's equal to zero vector

        Returns
        -------
            start_position: numpy array
                        start position of the parable
            
        '''
        if start_position is None:
            parab = self.start_position
        else:
            parab = start_position
        return parab

    def refflat(self, start_position=None):
        '''Function to set the start position of the reference flat

        Other Parameters
        ----------
        start_position: numpy array, optional
                        vector of six position
                        If start_position is not set it's equal to zero vector

        Returns
        -------
            start_position: numpy array
                        start position of the reference flat
        '''
        if start_position is None:
            refflat = self.start_position
        else:
            refflat = start_position
        return refflat

    def m4(self, start_position=None):
        '''Function to set the start position of the deformable mirror

        Other Parameters
        ----------
        start_position: numpy array, optional
                        vector of six position
                        If start_position is not set it's equal to zero vector

        Returns
        -------
            start_position: numpy array
                        start position of the deformable mirror
        '''
        if start_position is None:
            m4 = self.start_position
        else:
            m4 = start_position
        return m4
### Sensitivity matrices
    def _readMatFromTxt(self, file_name):
        ''' Function to read matrix of 11 Zernike x 6 displacements, m RMS, per 1 m displacement - or 1 radiant rotation

        Parameters
        ----------
        file_name: string
                    matrix file path

        Returns
        -------
                mat: numpy array [11,6]
                    matrix from txt file
        '''
        file = open(file_name, 'r')
        triplets=file.read().split()
        x = np.array(triplets)
        mat = x.reshape(11, 6)
        return mat.astype(float)

    def zmx_parpos2z(self):
        '''
        Returns
        -------
                mat: numpy array [11,6]
                    matrix parable positions to zernike
        '''
#         conffolder = os.path.join(path_name.CONFIGURATION_ROOT_FOLDER, tnconf)
#         file_name =  os.path.join(conffolder, 'ZST_PAR_pos2z.txt')
        file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_PAR_pos2z.txt'
        mat = self._readMatFromTxt(file_name)
        return mat

    def zmx_refflatpos2z(self):
        '''
        Returns
        -------
                mat: numpy array [11,6]
                    matrix reference flat positions to zernike
        '''
#         conffolder = os.path.join(path_name.CONFIGURATION_ROOT_FOLDER, tnconf)
#         file_name =  os.path.join(conffolder, 'ZST_FM_pos2z.txt')
        file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_FM_pos2z.txt'
        mat = self._readMatFromTxt(file_name)
        return mat

    def zmx_m4pos2z(self):
        '''
        Returns
        -------
                mat: numpy array [11,6]
                    matrix deformable mirror positions to zernike
        '''
#         conffolder = os.path.join(path_name.CONFIGURATION_ROOT_FOLDER, tnconf)
#         file_name =  os.path.join(conffolder, 'ZST_M4_pos2z.txt')
        file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_M4_pos2z.txt'
        mat = self._readMatFromTxt(file_name)
        return mat
# Zmat
    def zmat(self):
        '''
        Returns
        -------
            zmat: numpy array
            
        
        '''
        file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ott_mask.fits'
        hduList = pyfits.open(file_name)
        final_mask = np.invert(hduList[0].data.astype(bool))

        prova = np.ma.masked_array(np.ones(((2*OttParameters.parab_radius*OttParameters.pscale).astype(int),
                                            (2*OttParameters.parab_radius*OttParameters.pscale).astype(int))),
                                            mask=final_mask)
        zernike_mode = np.arange(2, 12)
        zmat_nopist = np.zeros((prova.compressed().shape[0], zernike_mode.size))
        for i in range(0, zernike_mode.size):
            z = self._zg.getZernike(zernike_mode[i])
            aa = np.ma.masked_array(z, mask=final_mask)
            zmat_nopist.T[i] = aa.compressed()
        zmat = np.ones((zmat_nopist.shape[0], zmat_nopist.shape[1]+1))
        for j in range(1, zmat_nopist.shape[1]):
            zmat[:, j] = zmat_nopist[:, j-1]
        return zmat

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

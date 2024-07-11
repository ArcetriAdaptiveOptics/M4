from astropy.io import fits
from arte.utils.snapshotable import Snapshotable
from m4.ground.read_data import read_phasemap
import datetime
from matplotlib import pyplot as plt
import functools
import pprint
from pathlib import Path, PurePath


class SurfaceMeasure():
    '''
    Class to handle the interferometer data in a time series
    '''

    def __init__(self, surface, overview, temperatures=None, zernikes=None):
        self._surface = surface
        self._zernikes = zernikes
        self._temperatures = temperatures
        self._overview = overview
        self._timestamp = datetime.datetime.strptime(
            overview['TN'], '%Y%m%d_%H%M%S')

    def __repr__(self):
        return pprint.PrettyPrinter().pformat(self.overview)

    @property
    def surface(self):
        return self._surface

    @property
    def temperatures(self):
        return self._temperatures

    @property
    def zernikes(self):
        return self._zernikes

    @property
    def shape(self):
        return self._surface.shape

    @property
    def overview(self):
        return Snapshotable.as_fits_header(self._overview)

    @property
    def tn(self):
        return self.overview['TN']

    @property
    def n_meas(self):
        return self.overview['N_MEAS']

    def save(self, filepath):
        # os.makedirs(os.path.dirname(filepath), exist_ok=True)
        # hdr = Snapshotable.as_fits_header(self._overview)
        # fits.writeto(filepath, self._interf, header=hdr, overwrite=True)
        pass

    @staticmethod
    def load(filepath_in):
        filename = PurePath(filepath_in)
        current_meas = str(filename.name)
        current_path = str(filename.parent)
        overview = Snapshotable.from_fits_header(fits.getheader(filename))
        if 'N_MEAS' not in list(overview.keys()):
            overview['N_MEAS'] = list(
                Path(current_path).iterdir()).index(filename)
        overview['TN'] = current_meas[0:15]

        imas = read_phasemap(str(filename))
        tmp = SurfaceMeasure._load_temperatures(current_path)
        if tmp is not None:
            temperatures = tmp[overview['N_MEAS'], :]
        else:
            temperatures = None
        tmp = SurfaceMeasure._load_zernikes(current_path)
        if tmp is not None:
            zernikes = tmp[overview['N_MEAS'], :]
        else:
            zernikes = None

        return SurfaceMeasure(imas, overview, temperatures, zernikes)

    @functools.cache
    @staticmethod
    def _load_temperatures(filepath):
        tfile = Path(filepath, 'temperature.fits')
        print(tfile)
        if Path.exists(tfile):
            with fits.open(tfile) as hdul:
                return hdul[0].data
        else:
            return None

    @functools.cache
    @staticmethod
    def _load_zernikes(filepath):
        zfile = Path(filepath, 'zernike.fits')
        if Path.exists(zfile):
            with fits.open(zfile) as hdul:
                return hdul[0].data
        else:
            return None

    def imshow(self, **kwargs):
        plt.imshow(self._surface, **kwargs)
        plt.colorbar()
        plt.show()

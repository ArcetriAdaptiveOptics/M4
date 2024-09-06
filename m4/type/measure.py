from astropy.io import fits
from arte.utils.snapshotable import Snapshotable
from m4.ground.read_data import read_phasemap
from m4.ground.timestamp import Timestamp
import datetime
from matplotlib import pyplot as plt
import functools
import pprint
from pathlib import Path, PurePath
import pandas as pd
import numpy as np
import copy
from datetime import datetime as dt
from m4.ground import zernike


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
        self._timestamp = Timestamp.fromString(overview['TN']).datetime

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
    def timestamp(self):
        return self._timestamp

    @property
    def n_meas(self):
        return int(self.overview['N_MEAS'])

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
        tmp = SurfaceMeasure.load_temperatures(current_path)
        if tmp is not None:
            temperatures = tmp[overview['N_MEAS'], :]
        else:
            temperatures = None
        tmp = SurfaceMeasure.load_zernikes(current_path)
        if tmp is not None:
            zernikes = tmp[overview['N_MEAS'], :]
        else:
            zernikes = None

        return SurfaceMeasure(imas, overview, temperatures, zernikes)

    @functools.lru_cache(maxsize=None)
    @staticmethod
    def load_temperatures(filepath):
        tfile = Path(filepath, 'temperature.fits')
        if Path.exists(tfile):
            with fits.open(tfile) as hdul:
                return hdul[0].data
        else:
            return None

    @functools.lru_cache(maxsize=None)
    @staticmethod
    def load_zernikes(filepath):
        zfile = Path(filepath, 'zernike.fits')
        if Path.exists(zfile):
            with fits.open(zfile) as hdul:
                return hdul[0].data
        else:
            return None


class SurfaceSequence():
    '''
    Class to handle a sequence of surface measurements
    Usage:


        from m4.utils.osutils import FileWalker
        from m4.type.measure import SurfaceSequence, SurfaceMeasure

        % init
        fw = FileWalker(path)
        tn = '20230715_150944'
        fullopath = fw.findTracknum(tn)
        seq = SurfaceSequence(fullopath)
        print(seq[1].surface.std())

        %use static methods to load temperatures and zernikes
        aa = SurfaceMeasure.load_temperatures(fullopath[0])
        print(aa)
        aa = SurfaceMeasure.load_zernikes(fullopath[0])

        %slicing
        seq2 = seq[0:50]

        %compute structure function bulk (fitting zernike polynomials each time)
        time, stf = seq2.compute_time_stf_bulk()
        seq2.compute_intersection_mask()


        %compute structure function bulk (fitting zernike polynomials on intersection mask)
        time, stf = seq2.compute_fast_time_stf()
        plt.figure()
        plt.plot(time, stf)


    '''

    def __init__(self, root_dirs):

        if not isinstance(root_dirs, list):
            root_dirs = [root_dirs]

        if root_dirs[0].is_file():
            print("Sequence from file not implemented yet")
            pass

        else:
            for root_dir in root_dirs:
                if root_dir == root_dirs[0]:
                    self._root_dir = root_dir
                    self._filenames = list(
                        Path(root_dir).rglob('????????_??????.fits'))
                    print("found  %d frames" % len(self._filenames))
                    tmp_tmp = SurfaceMeasure.load_temperatures(root_dir)
                    if tmp_tmp is not None:
                        self._temperature = pd.DataFrame(
                            tmp_tmp.byteswap().newbyteorder('='))
                    tmp_zn = SurfaceMeasure.load_zernikes(root_dir)
                    if tmp_zn is not None:
                        self._zernike = pd.DataFrame(
                            tmp_zn.byteswap().newbyteorder('='))

                    self._tau_vector = np.array(list(map(lambda l: dt.strptime(
                        l.name[0:15], '%Y%m%d_%H%M%S'), self._filenames)))
                else:
                    self._filenames += list(
                        Path(root_dir).rglob('????????_??????.fits'))
                    print("found  %d frames" % len(self._filenames))
                    tmp_tmp = SurfaceMeasure.load_temperatures(root_dir)
                    self._temperature = pd.concat(
                        [self._temperature, pd.DataFrame(tmp_tmp.byteswap().newbyteorder('='))])
                    tmp_zn = SurfaceMeasure.load_zernikes(root_dir)
                    self._zernike = pd.concat(
                        [self._zernike, pd.DataFrame(tmp_zn.byteswap().newbyteorder('='))])
                    self._tau_vector = np.concatenate(self._tau_vector, np.array(list(map(lambda l: dt.strptime(
                        str(l.name)[0:15], '%Y%m%d_%H%M%S'), self._filenames))))
            self._intersetion_mask = None
            print("Sequence loaded with %d frames" % len(self))

    def __len__(self):
        return len(self._filenames)

    def __getitem__(self, idx):

        if isinstance(idx, slice):
            new_seq = copy.deepcopy(self)
            new_seq._filenames = self._filenames[idx]
            new_seq._temperature = self._temperature.iloc[idx]
            new_seq._zernike = self._zernike.iloc[idx]
            new_seq._tau_vector = self._tau_vector[idx]
            return new_seq
        else:
            img_path = Path(self._filenames[idx])
            meas = SurfaceMeasure.load(img_path)
        return meas

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __repr__(self):
        return f"SurfaceSequence({self._root_dir}  [{len(self)} frames])"

    def __str__(self):
        return self.__repr__()

    @ property
    def temperture(self):
        return self._temperture

    @ property
    def zernike(self):
        return self._zernike

    @ property
    def tau_vector(self):
        return self._tau_vector

    @ property
    def tau_vector_as_timestamp(self):
        return np.array(list(map(lambda l: dt.timestamp(l), self._tau_vector)))

    def set_tau_vector(self, tau_vector):
        if not isinstance(tau_vector, np.ndarray):
            raise TypeError("tau_vector must be a numpy array")
        if not isinstance(tau_vector[0], dt):
            raise TypeError(
                "tau_vector must be a numpy array of datetime objects")
        if len(tau_vector) != len(self):
            raise ValueError(
                "tau_vector must have the same length as the sequence")

        self._tau_vector = tau_vector

    def set_intersection_mask(self, mask):
        self._intersetion_mask = mask

    @ property
    def intersection_mask(self):
        return self._intersetion_mask

    def compute_intersection_mask(self):
        for ff in self:

            if self._intersetion_mask is None:
                self._intersetion_mask = ff.surface.mask
            else:
                self._intersetion_mask = self._intersetion_mask | ff.surface.mask

    def compute_time_stf_bulk(self, nzern=3):
        '''Compute the structure function of the sequence of surfaces fitting the first nzern zernike polynomials'''

        if nzern is not None:
            zv = np.arange(nzern) + 1

        dtime = self.tau_vector_as_timestamp
        xx, yy = np.meshgrid(dtime, dtime)
        dd = np.sqrt((xx-xx.T)**2 + (yy-yy.T)**2)
        datalen = len(self)
        stf = np.zeros(datalen-1)
        time = np.zeros(datalen-1)

        stf_mat = np.zeros((datalen, datalen))
        for i in range(datalen):
            print(i, end=" ")

            if nzern is not None:
                coeff, mat = zernike.zernikeFit(self[i].surface, zv)
                sur = zernike.zernikeSurface(self[i].surface, coeff, mat)
                ima_i = self[i].surface - sur

            else:
                ima_i = self[i].surface

            for j in range(datalen-i):
                if nzern is not None:
                    coeff, mat = zernike.zernikeFit(self[i+j].surface, zv)
                    sur = zernike.zernikeSurface(self[i+j].surface, coeff, mat)
                    ima_ij = self[i+j].surface - sur

                else:
                    ima_ij = self[i+j].surface
                dd[i+j, j] = 0
                stf_mat[i, i+j] = np.ma.std(ima_i - ima_ij)**2
            # calculate the structure function as the average of the elements of the matrix which distance is the same
        for i in range(datalen-1):
            stf[i] = np.sqrt(np.mean(stf_mat[dd == dd[0, i]]))
            time[i] = dd[0, i]

        return time, stf

    def compute_fast_time_stf(self, nzern=3):
        '''Compute the structure function of the sequence of surfaces fitting the first nzern zernike polynomials on intersection mask'''

        if nzern is not None:
            zv = np.arange(nzern) + 1

        dtime = self.tau_vector_as_timestamp
        xx, yy = np.meshgrid(dtime, dtime)
        dd = np.sqrt((xx-xx.T)**2 + (yy-yy.T)**2)
        datalen = len(self)
        stf = np.zeros(datalen-1)
        time = np.zeros(datalen-1)

        stf_mat = np.zeros((datalen, datalen))
        if self._intersetion_mask is None and nzern is not None:
            print("Computing intersection mask", end="...")
            self.compute_intersection_mask()
            print("done")

        print("Computing structure function")
        print("Frame distance:", end=" ")
        for i in range(datalen):
            print(i, end=" ")

            if nzern is not None:

                data = np.ma.masked_array(
                    self[i].surface.data, mask=self._intersetion_mask)
                if i == 0:
                    coeff, mat = zernike.zernikeFit(data, zv)
                    mat_pinv = np.linalg.pinv(mat)
                else:
                    coeff = np.dot(mat_pinv, data[data.mask == False])
                sur = zernike.zernikeSurface(data, coeff, mat)
                ima_i = data - sur

            else:
                ima_i = data

            for j in range(datalen-i):
                if nzern is not None:

                    data = np.ma.masked_array(
                        self[i+j].surface.data, mask=self._intersetion_mask)
                    coeff = np.dot(mat_pinv, data[data.mask == False])

                    sur = zernike.zernikeSurface(data, coeff, mat)
                    ima_ij = data - sur

                else:
                    ima_ij = self[i+j].surface
                dd[i+j, j] = 0
                stf_mat[i, i+j] = np.ma.std(ima_i - ima_ij)**2
            # calculate the structure function as the average of the elements of the matrix which distance is the same
        for i in range(datalen-1):
            stf[i] = np.sqrt(np.mean(stf_mat[dd == dd[0, i]]))
            time[i] = dd[0, i]

        return time, stf

    def compute_time_stf_w_preloaded_data(self, nzern=3):
        '''Compute the structure function of the sequence of surfaces using intersection mask and zernike coefficients from the folder'''

        if nzern is not None:
            zv = np.arange(nzern) + 1

        dtime = self.tau_vector_as_timestamp
        xx, yy = np.meshgrid(dtime, dtime)
        dd = np.sqrt((xx-xx.T)**2 + (yy-yy.T)**2)
        datalen = len(self)
        stf = np.zeros(datalen-1)
        time = np.zeros(datalen-1)

        stf_mat = np.zeros((datalen, datalen))
        if self._intersetion_mask is None and nzern is not None:
            print("Computing intersection mask", end="...")
            self.compute_intersection_mask()
            print("done")

        zern_coeff = self._zernike.to_numpy()
        print("Computing structure function")
        print("Frame distance:", end=" ")
        for i in range(datalen):
            print(i, end=" ")

            if nzern is not None:

                data = np.ma.masked_array(
                    self[i].surface.data, mask=self._intersetion_mask)
                if i == 0:
                    coeff, mat = zernike.zernikeFit(data, zv)
                    # print(coeff)
                    # print(zern_coeff[i, 0:nzern])
                    # print(data.mean())
                    # print(self[i].surface.mean())

                    plt.imshow(data)
                    try:
                        np.testing.assert_allclose(
                            coeff, zern_coeff[i, 0:nzern], rtol=0.1)
                    except AssertionError:
                        print(
                            "WARNING: Zernike coefficients in folder are not consistent with the data")
                sur = zernike.zernikeSurface(data, zern_coeff[i, 0:nzern], mat)
                ima_i = data - sur

            else:
                ima_i = data

            for j in range(datalen-i):
                if nzern is not None:

                    data = np.ma.masked_array(
                        self[i+j].surface.data, mask=self._intersetion_mask)
                    sur = zernike.zernikeSurface(
                        data, zern_coeff[i+j, 0:nzern], mat)
                    ima_ij = data - sur

                else:
                    ima_ij = self[i+j].surface
                dd[i+j, j] = 0
                stf_mat[i, i+j] = np.ma.std(ima_i - ima_ij)**2
            # calculate the structure function as the average of the elements of the matrix which distance is the same
        for i in range(datalen-1):
            stf[i] = np.sqrt(np.mean(stf_mat[dd == dd[0, i]]))
            time[i] = dd[0, i]

        return time, stf

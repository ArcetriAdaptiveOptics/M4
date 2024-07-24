import numpy as np
from m4.type.measure import Measure
from minigrant.file_walker import FileWalker
from minigrant.package_data import data_root_dir
from appoppy.mask import sector_mask
from arte.types.mask import BaseMask
from arte.dataelab.cache_on_disk import cache_on_disk
from arte.dataelab.base_timeseries import BaseTimeSeries
from arte.dataelab.base_analyzer import BaseAnalyzer


def mask_sectors(ima, angs, centres, radii):
    mask_sectors = []
    for ang, cen, radius in zip(angs, centres, radii):
        smask = sector_mask(ima.shape,
                            (ang[1], ang[0]),
                            centre=cen,
                            radius=radius)
        mask_sectors.append(smask)
    return mask_sectors


def wrap_opd(opd_in_nm, wv):
    b = opd_in_nm / wv * 2 * np.pi
    c = np.cos(b)
    s = np.sin(b)
    return np.arctan2(s, c) / (2 * np.pi) * wv


def unwrap_opd(opd_nm, wv):
    phase = opd_nm * 2 * np.pi / wv
    return unwrap_phase(phase) * wv / (2 * np.pi)


def unwrap_opd_cedricco(opd_nm, wv):
    for i in np.arange(1, opd_nm.shape[0]):
        if opd_nm[i-1] - opd_nm[i] > wv / 2:
            opd_nm[i:] += wv
    return opd_nm


class InterferogramsTimeSeries(BaseTimeSeries):

    def __init__(self, data, time_vector=None, **kwargs):
        super().__init__(data, time_vector, **kwargs)


class OpdTimeSeries(BaseTimeSeries):

    def __init__(self, data, time_vector=None, **kwargs):
        super().__init__(data, time_vector, **kwargs)


class PistonTimeSeries(BaseTimeSeries):

    def __init__(self, data, time_vector=None, **kwargs):
        super().__init__(data, time_vector, **kwargs)


class Analyzer(BaseAnalyzer):
    '''
    '''

    def __init__(self, tn, recalc=False):

        super().__init__(tn, recalc)

        self._tn = tn

    @cache_on_disk
    def measure(self):
        return Measure.load(FileWalker(data_root_dir()).rsi_measure_path(self._tn))

    def _read_interferograms(self):
        return InterferogramsTimeSeries(self.measure().interferograms)

    @cache_on_disk
    def interferograms(self):
        return self._read_interferograms().get_data()

    @cache_on_disk
    def interferograms_time_average(self):
        return self._read_interferograms().time_average()

    def wavelength(self):
        return self.measure().wavelength

    def analysis_mask(self):
        mask = np.ones(self.interferograms_time_average().shape)
        for sector in self.mask_sectors():
            mask[sector] = 0
        return mask

    @cache_on_disk
    def mask_sectors(self):
        return mask_sectors(self.interferograms_time_average(),
                            self.measure().mask_sectors_angs,
                            self.measure().mask_sectors_centres,
                            self.measure().mask_sectors_radii)

    def _opd(self, interferogram):
        rec = OpdReconstructor(interferogram, self.wavelength())
        rec.set_mask(BaseMask(self.analysis_mask()))
        rec.compute_intensity_fourier_transform()
        rec.set_window(
            TopFlatGaussianCircularWindow(self.measure().window_x,
                                          self.measure().window_y,
                                          self.measure().window_sigma))
        rec.compute_intensity_filtering()
        rec.compute_phase()
        opd = rec.get_unwrapped_OPD()
        return opd

    @cache_on_disk
    def opd_on_time_average(self):
        return OpdTimeSeries(np.expand_dims(self._opd(self.interferograms_time_average()), 0))

    @cache_on_disk
    def opd_time_series(self):
        return OpdTimeSeries(
            np.ma.stack([self._opd(x) for x in self.interferograms()]))

    def piston_on_sector_on_time_average(self, sector_idx):
        ima = self.opd_on_time_average().get_data()[0]
        mean_opd = np.mean(ima[self.mask_sectors()[sector_idx]])
        return wrap_opd(mean_opd, self.wavelength())

    def piston_on_sector_list(self, sector_idx):
        ima_list = self.opd_time_series().get_data()
        p_list = []
        for ima in ima_list:
            mean_opd = np.mean(ima[self.mask_sectors()[sector_idx]])
            p_list.append(wrap_opd(mean_opd, self.wavelength()))
        return PistonTimeSeries(np.expand_dims(np.array(p_list), 1))

    # def _compute_zernike_after_ciaociao(self, n_modes, circ_mask):
    #     if circ_mask is None:
    #         circular_mask = self.mask
    #     else:
    #         circular_mask = circ_mask
    #     opd_map = self.get_unwrapped_OPD()
    #     md = ModalDecomposer(n_modes)
    #     zern_coeff = md.measureZernikeCoefficientsFromWavefront(
    #         wavefront=Wavefront.fromNumpyArray(opd_map),
    #         circular_mask=circular_mask,
    #         user_mask=self.mask
    #     )
    #     self._zernike_out = zern_coeff.toNumpyArray()

    # def get_zernike_at_ciaociao_entrance(self, zernike_j, circ_mask=None):
    #     '''
    #     Return the Zernike coefficients in [um] RMS representing the aberrations at the CiaoCiao entrance.

    #     Parameters
    #     ----------
    #     zernike_j: ~numpy.ndarray
    #         List of Zernike indices. It must start from j=2.
    #     '''
    #     czr = CiaoCiaoZernikeReconstructor(zernike_j, self._rot)
    #     rec_mat = czr.get_reconstructor()
    #     try:
    #         zernike_in = np.dot(
    #             rec_mat,
    #             self.get_zernike_at_ciaociao_exit(zernike_j, circ_mask))
    #     except ValueError:
    #         zernike_in = np.dot(
    #                 rec_mat,
    #                 self.get_zernike_at_ciaociao_exit(
    #                     np.append(zernike_j, zernike_j[-1]+1), circ_mask)
    #         )
    #     return zernike_in

    # def get_zernike_at_ciaociao_exit(self, zernike_j, circ_mask=None):
    #     '''
    #     Return the Zernike coefficients in [um] RMS representing the aberrations at the CiaoCiao exit.

    #     Parameters
    #     ----------
    #     zernike_j: ~numpy.ndarray
    #         List of Zernike indices. It must start from j=2.
    #     '''
    #     self._compute_zernike_after_ciaociao(n_modes=len(zernike_j), circ_mask=circ_mask)
    #     return self._zernike_out


class AnalyzerWithReference(Analyzer):

    def __init__(self, tn, tn_ref, recalc):
        super().__init__(tn, recalc)
        self._tn_ref = tn_ref

    @cache_on_disk
    def opd_on_time_average(self):
        opd = OpdTimeSeries(np.expand_dims(
            self._opd(self.interferograms_time_average()), 0))
        ref_opd = Analyzer(self._tn_ref).opd_on_time_average()
        return opd - ref_opd

    @cache_on_disk
    def opd_time_series(self):
        opd = OpdTimeSeries(np.ma.stack(
            [self._opd(x) for x in self.interferograms()]))
        ref_opd = Analyzer(self._tn_ref).opd_on_time_average()
        return opd - ref_opd

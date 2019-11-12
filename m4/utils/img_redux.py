'''
@author: cs
'''

import logging
import numpy as np
from m4.ground.configuration import Configuration
from m4.utils.roi import ROI
from m4.utils.zernike_on_m_4 import ZernikeOnM4


class TipTiltDetrend():

    def __init__(self):
        self._logger = logging.getLogger('TIP_TILT_DETREND:')
        self._pupilXYRadius = Configuration.PARABOLA_PUPIL_XYRADIUS
        self._zOnM4 = ZernikeOnM4()
        self._totalMatList = None

    def tipTiltRemover(self, image, roi, final_index, analysis_ind=None):
        roi_copy = np.copy(roi)
        self._logger.debug('Removal of tip-tilt from roi[%d]',
                           final_index)
        '''
            arg:
                image= immagine da analizzare (np.ma-masked_array)
                roi= roi dell'immagine (list)
                final_index= indice della roi finale (int)
                analysis_ind= indice delle roi da utilizzare per
                                l'analisi (np.array[...])
        '''
        self._totalMatList = []
        coefList = []
        for r in roi:
            imag = np.ma.masked_array(image.data, mask=r)
            ima = np.ma.MaskedArray.copy(imag)
            coef, mat = self._zOnM4.zernikeFit(ima, np.array([2, 3]))
            self._totalMatList.append(mat)
            coefList.append(coef)

        if analysis_ind is None:
            coef_list = coefList
            del coef_list[final_index]
        else:
            coef_list = []
            for i in range(len(analysis_ind)):
                coef_list.append(coefList[analysis_ind[i]])
        tip, tilt = np.average(coef_list, axis=0)

        surfcoef = np.array([tip, tilt])
        surface_map = \
                    self._zOnM4.zernikeSurface(surfcoef, roi_copy[final_index],
                                               self._totalMatList[final_index])

        cx = self._pupilXYRadius[0]
        cy = self._pupilXYRadius[1]
        r = self._pupilXYRadius[2]
        ima_cut = image[cy-r:cy+r, cx-r:cx+r]
        image_ttr = np.ma.masked_array(ima_cut.data - surface_map,
                                       mask=roi[final_index])

        return image_ttr






class PhaseSolve():

    def __init__(self):
        self._roi = ROI()
        self._lambda = Configuration.LAMBA
        self._n = None


    def n_calculator(self, spl_values):
        n = np.zeros(spl_values.shape[0])
        for i in range(spl_values.shape[0]):
            n[i] = (2.* spl_values[i]) / self._lambda
        self._n = n
        return self._n


    def m4PhaseSolver(self, m4_ima, spl_values):
        self.n_calculator(spl_values)
        roiList = self._roi._ROIonM4(m4_ima)
        m4_new_image = None

        media = []
        imgList = []
        for roi in roiList:
            imgg = np.ma.masked_array(m4_ima.data, mask=roi)
            m = imgg.mean()
            media.append(m)
            imgList.append(imgg)

        aa = np.arange(self._n.shape[0])
        zipped = zip(aa, imgList)
        img_phase_solve_list = []
        for i, imgg in zipped:
            img_phase_solve = np.ma.masked_array(imgg.data - self._n[i],
                                                mask=imgg.mask)
            img_phase_solve_list.append(img_phase_solve)

        img_phase_solve_list[len(img_phase_solve_list)-1] = \
                       np.ma.masked_array(imgList[len(imgList)-2].data, 
                                          mask= imgList[len(imgList)-2].mask)


        for j in range(1, len(img_phase_solve_list)):
            if m4_new_image is None:
                m4_new_image = np.ma.array(img_phase_solve_list[0].filled(1) * \
                                           img_phase_solve_list[j].filled(1),
                                           mask=(img_phase_solve_list[0].mask * \
                                                 img_phase_solve_list[j].mask))
            else:
                m4_new_image = np.ma.array(m4_new_image.filled(1) * \
                                           img_phase_solve_list[j].filled(1),
                                           mask=(m4_new_image.mask * \
                                                 img_phase_solve_list[j].mask))

        return m4_new_image, img_phase_solve_list, imgList



    def masterRoiPhaseSolver(self, seg_ima, spl_value):
        self.n_calculator(spl_value)
        roiList = self._roi.roiGenerator(seg_ima)

        imgg = np.ma.masked_array(seg_ima.data, mask=roiList[3])

        img_phase_solve = np.ma.masked_array(imgg.data - self._n,
                                            mask=imgg.mask)

        return img_phase_solve

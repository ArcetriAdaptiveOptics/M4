'''
Authors
  - C. Selmi: written in 2021
'''
import os
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration import config_folder_names as fold_name
from m4.ground import tracking_number_folder

IP = '192.168.1.18'
PORT = 7100 #ma anche 7110

class Pointer_align():
    '''
    https://github.com/ArcetriAdaptiveOptics/pysilico
    '''

    def __init__(self, pysilico):
        self._camera = pysilico(IP, PORT)
        self._exposureTimeInMilliSeconds = 0.01
        self._dl = 0.60
        self._dl0 = 0.5
        self._pixs = 7.2e-6
        self._rad2as = 206265 #[as/ras]

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.POINTER_ALIGN_ROOT_FOLDER

    def _saveInfo(self):
        self._dove, self._tt = tracking_number_folder.createFolderToStoreMeasurements(self._storageFolder())
        fits_file_name = os.path.join(self._dove, 'BenchLength-PixScale.fits')
        pyfits.writeto(fits_file_name, [self._dl0, self._dl, self._pixs])

    def _saveImage(self, image, name):
        fits_file_name = os.path.join(self._dove, name)
        pyfits.writeto(fits_file_name, image)

    def doAlign(self, dy_offset=None):
        if dy_offset is not None:
            print('che devo ruotare?') #box

        images = self.take_images(2)
        coord0 = self.centroid_calculator(images[0])
        coord1 = self.centroid_calculator(images[1])

        dy = coord1[1]-coord0[1]
        print('Y offset [um, as]')
        print(dy*self._pixs*1e6, dy*self._pixs/self._dl * self._rad2as)

        image = images[0] + 2*images[1]
        return dy, image

    def doDecenter(self):
        'da finire'
        images = self.take_images(2)
        coord0 = self.centroid_calculator(images[0])
        coord1 = self.centroid_calculator(images[1])

        image = images[0] + 2*images[1]
        d = np.max(image) #p?

        dc = coord1[1]-coord0[1]
        print('Y decenter [um]')
        print(dc*self._pixs*1e6)
        return image




    def main(self):
        dy0, image_a = self.doAlign()
        self._saveImage(image_a, 'xTilt.fits')

        dy1, image_b = self.doAlign(dy0)
        self._saveImage(image_b, 'yTilt.fits')

        image_c = self.doDecenter()
        self._saveImage(image_c, 'xDecenter.fits')






    def take_images(self, numberOfReturnedImages):
        '''
        '''
        self._camera.setExposureTime(self._exposureTimeInMilliSeconds)
        images = self._camera.getFutureFrames(numberOfReturnedImages,
                                     numberOfFramesToAverageForEachImage=1)
                                    #timeout
        return images

    def centroid_calculator(self, data):
        from photutils import centroids
        x, y = centroids.centroid_com(data, mask=None) #Calculates the object “center of mass” from 2D image moments
        x, y = centroids.centroid_1dg(data, error=None, mask=None) #Calculates the centroid by fitting 1D Gaussians to the marginal x and y distributions of the data.
        x, y = centroids.centroid_2dg(data, error=None, mask=None) #Calculates the centroid by fitting a 2D Gaussian to the 2D distribution of the data.
        #https://photutils.readthedocs.io/en/v0.3.1/photutils/centroids.html
        par = centroids.fit_2dgaussian(data)._parameters #uguale a 2dg
        return [x, y]


    def center_calculator(self, data):
        pass

'''
Authors
  - C. Selmi: written in 2021
'''
import os
from photutils import centroids
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
from m4.type.pointer_align import PointerAlign
from m4.ground import tracking_number_folder

IP = '192.168.1.18'
PORT = 7100 #ma anche 7110

#vecchio IP dal file conf di Lorenzo 192.168.29.159

class PointerAligner():
    '''
    https://github.com/ArcetriAdaptiveOptics/pysilico
    '''

    def __init__(self, pysilico, pointerId):
        '''
        Parameters
        ----------
        pysilico: object
            camera forframe acquisition
        poniterId: string
            'NGS' or 'LGS'
        '''
        self._camera = pysilico(IP, PORT)
        self._idData = PointerAlign(pointerId)
        self._tt = None
        self._dove = None

    def main(self):
        self._dove, self._tt = tracking_number_folder.createFolderToStoreMeasurements(self._idData.folder)
        self._idData.dove = self._dove
        self._idData.tt = self._tt
        y00, y01, image_a = self.dyCalculator()
        self._saveImage(image_a, 'xTilt.fits')
        print('Move camera away')
        self._pause() #spostare camera
        y11, y10, image_b = self.dyCalculator()
        self._saveImage(image_b, 'yTilt.fits')

        point0, point1 = self._idData.operationAndPrint(y00, y01, y11, y10)
        self._idData.saveInfo()

        self.target_center(point0) #1 vite
        self._pause()
        self.target_center(point1) #2 vite


    def dyCalculator(self):
        images = self.take_images(1)
        coord0 = self.centroid_calculator(images[0])
        print('Rotate ±180°')
        self._pause() #rotazione
        images = self.take_images(1)
        coord1 = self.centroid_calculator(images[0])

        dy = coord1[1]-coord0[1]
        print('Y offset [um, as]')
        print(dy*self._idData.pixs*1e6) #(dy*self._pixs/self._dl * self._rad2as)

        image = images[0] + 2*images[1]
        return coord0[1], coord1[1], image


    def target_center(self, point):
        '''
        Function that acquires an image and draws a cross on the input point

        Parameters
        ----------
        point: list
            point y coord
        '''
        image = self.take_images(1)
        plt.imshow(image, origin='lower')
        plt.colorbar()
        size = 250
        plt.scatter(point, point, s=size, c='red', marker='+')

    def _saveImage(self, image, name):
        fits_file_name = os.path.join(self._dove, name)
        pyfits.writeto(fits_file_name, image)

    def _pause(self):
        '''this function will pause the script with a default massage'''
        os.system("read -p 'Press Enter to continue...' var")
        return

    def take_images(self, numberOfReturnedImages):
        '''
        Parameters
        ----------
        numberOfReturnedImages: int
            numebers of sequential images

        Returns
        -------
        images: ?
        '''
        self._camera.setExposureTime(self._idData.exposureTimeInMilliSeconds)
        images = self._camera.getFutureFrames(numberOfReturnedImages,
                                     numberOfFramesToAverageForEachImage=1)
                                    #timeout
        return images

    def centroid_calculator(self, data):
        '''
        Parameters
        ----------
        data: numpy array
            image

        Returns
        -------
        coord: list
            centroid's coordinates
        '''
        x, y = centroids.centroid_com(data, mask=None) #Calculates the object “center of mass” from 2D image moments
        x, y = centroids.centroid_1dg(data, error=None, mask=None) #Calculates the centroid by fitting 1D Gaussians to the marginal x and y distributions of the data.
        x, y = centroids.centroid_2dg(data, error=None, mask=None) #Calculates the centroid by fitting a 2D Gaussian to the 2D distribution of the data.
        #https://photutils.readthedocs.io/en/v0.3.1/photutils/centroids.html
        par = centroids.fit_2dgaussian(data)._parameters #uguale a 2dg
        return [x, y]



# def test():
#     from m4.pointer_aligner import PointerAligner
#     from m4.configuration import start
#     conf = '/Users/rm/eclipse-workspace/M4/m4/configuration/myConfig.yaml'
#     ott, interf = start.create_ott(conf)
#     def camera(IP, PORT): 
#         IP = 'ciao' 
#         PORT = 'riciao'
#     p = PointerAligner(camera, 'NGS')

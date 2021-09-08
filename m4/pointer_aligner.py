'''
Authors
  - C. Selmi: written in 2021
'''
import os
import numpy as np
from photutils import centroids
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
from m4.configuration import config_folder_names as fold_name
from m4.ground import tracking_number_folder

IP = '192.168.1.18'
PORT = 7100 #ma anche 7110

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
        self._exposureTimeInMilliSeconds = None
        self._pixs = 7.2e-6
        self._rad2as = 206265 #[as/ras]
        self._l1 = None
        self._l2 = None
        self._p = None #distance between actuators
        self._defineParameters(pointerId)

    def _defineParameters(self, pointerId):
        '''
        Parameters
        ----------
        poniterId: string
            'NGS' or 'LGS'
        '''
        if pointerId == 'NGS':
            self._exposureTimeInMilliSeconds = 0.0002
            self._l1 = 0.495
            self._l2 = 1.163
            self._p = 0.038
            self._folder = os.path.join(fold_name.POINTER_ALIGN_ROOT_FOLDER, 'NGSPointer')
        if pointerId == 'LGS':
            self._exposureTimeInMilliSeconds = 0.001*5
            self._l1 = 0.485 #distance from II actuator and first measurement point
            self._l2 = 1.03 #distance from II actuator and second measurement point
            self._p = 0.033
            self._folder = os.path.join(fold_name.POINTER_ALIGN_ROOT_FOLDER, 'LGSPointer')


    def main(self):
        self._dove, self._tt = tracking_number_folder.createFolderToStoreMeasurements(self._folder)
        y00, y01, image_a = self.dyCalculator()
        self._saveImage(image_a, 'xTilt.fits')
        print('Move camera away')
        self._pause() #spostare camera
        y11, y10, image_b = self.dyCalculator()
        self._saveImage(image_b, 'yTilt.fits')

        point0, point1 = self.operationAndPrint(y00, y01, y11, y10)

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
        print(dy*self._pixs*1e6) #(dy*self._pixs/self._dl * self._rad2as)

        image = images[0] + 2*images[1]
        return coord0[1], coord1[1], image

    def operationAndPrint(self, y00, y01, y11, y10):
        m0  = np.mean([y00,y01]) #mid-point at the first distance
        m1  = np.mean([y10,y11]) #mid point at the second distance
        r0  = y00 - m0
        r1  = y10 - m1
        s0  = (m1-m0)/(self._l2-self._l1) #slope of the rotation axis

        s   = (y10-m1-(y00-m0))/(self._l2-self._l1) #slope of the beam wrt the rotation axis
        z0  = np.mean([y10,y11]) #== m1, expected beam position at second distance whit no tilt and no decenter
        h2 = (y10-m1)-s*self._l2 #error (displacement) of the actuator II
        h1 = (y10-m1)-s*(self._l2+self._p) #error (displacement) of the actuator I
        z1 = m1-h1/self._p*self._l2 #expected beam position 

        self._beamSlope = s*self._pixs*self._rad2as
        self._rotationAxsSl = (m1-m0)/(self._l2-self._l1)*self._pixs*self._rad2as
        self._radii = (np.abs(y10-y11)-np.abs(y00-y01))*self._pixs*1e6
        self._meanOffs = np.mean([h1*self._pixs*1e6, h2*self._pixs*1e6])
        self._offset1 = h1*self._pixs*1e6
        self._offset2 = h2*self._pixs*1e6
        print('Tracknum: %s' %self._tt)
        print('Beam slope [arcsec]: %f' %self._beamSlope)
        print('Radii difference [um]: %f' %self._radii)
        print('Rotation axis slope [arcsec]: %f' %self._rotationAxsSl)
        print('Pointer offset at actuators [um]: %f %f' %(self._offset1, self._offset2) )

        return m0, m1

    def target_center(self, point):
        '''
        Function that acquires an image and draws a cross on the input point

        Parameters
        ----------
        point: list
            point coordinates
        '''
        image = self.take_images(1)
        plt.imshow(image, origin='lower')
        plt.colorbar()
        size = 250
        plt.scatter(point[0], point[1], s=size, c='red', marker='+')

    def _saveInfo(self):
        txt_file_name = os.path.join(self._dove, 'BenchLength-Info.txt')
        file = open(txt_file_name, 'w+')
        file.write('tt = %s \n' %self._tt)
        file.write('Beam slope = %f \n' %self._beamSlope)
        file.write('Rot.axis sl. = %f \n' %self._rotationAxsSl)
        file.write('Mean offs. [um] = %f \n' %self._meanOffs)
        file.write('Offset I [um] = %f \n' %self._offset1)
        file.write('Offset II [um] = %f \n' %self._offset2)
        file.close()

        fits_file_name = os.path.join(self._dove, 'BenchLength-Info.fits')
        header = pyfits.Header()
        header['TT'] = self._tt
        header['BEAMSL'] = self._beamSlope
        header['ROTAXSL'] = self._rotationAxsSl
        header['MEANOFFS'] = self._meanOffs
        pyfits.writeto(fits_file_name, np.array([self._offset1, self._offset2]), header)

#     @staticmethod
#     def readInfo(self, tt):
#         theObject = PointerAligner()
#         theObject._dove = tracking_number_folder.findTrackingNumberPath(tt)
#         filename = os.path.join(self._dove, 'BenchLength-Info.fits')
#         header = pyfits.getheader(filename)
#         hduList = pyfits.open(filename)
# 
#         theObject._tt = header['TT']
#         theObject._beamSlope = header['BEAMSL']
#         theObject._rotationAxsSl = header['ROTAXSL']
#         theObject._meanOffs = header['MEANOFFS']
# 
#         offsets = hduList[0].data
#         theObject._offset1 = offsets[0]
#         theObject._offset2 = offsets[1]

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
        self._camera.setExposureTime(self._exposureTimeInMilliSeconds)
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
    



def test():
    from m4.pointer_aligner import PointerAligner
    from m4.configuration import start
    conf = '/Users/rm/eclipse-workspace/M4/m4/configuration/myConfig.yaml'
    ott, interf = start.create_ott(conf)
    def camera(IP, PORT): 
        IP = 'ciao' 
        PORT = 'riciao'
    p = PointerAligner(camera, 'NGS')

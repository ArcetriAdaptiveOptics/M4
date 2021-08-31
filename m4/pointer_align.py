'''
Authors
  - C. Selmi: written in 2021
'''

IP = '192.168.1.18'
PORT = 7100 #ma anche 7110

class Pointer_align():
    '''
    https://github.com/ArcetriAdaptiveOptics/pysilico
    '''

    def __init__(self, pysilico):
        self._camera = pysilico(IP, PORT)
        self._exposureTimeInMilliSeconds = 0.01

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


    def center_calculator(self, data):
        pass

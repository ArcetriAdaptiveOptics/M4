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

    def centroid_calculator(self):
        pass

    def center_calculator(self):
        pass

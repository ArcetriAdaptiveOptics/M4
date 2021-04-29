'''
Authors
  - C. Selmi: written in 2020
'''

class BaseAccelerometers():
    '''
    Abstract class for parabola
    '''

    def acquireData(self, recording_seconds):
        ''' Function for data acquisition

        Parameters
        ----------
        recording_seconds: int [s]
            recording seconds for data acquisition
        '''
        raise Exception('Implement me!')

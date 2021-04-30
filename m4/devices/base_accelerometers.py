'''
Authors
  - C. Selmi: written in 2020
'''
from abc import ABC, abstractmethod


class BaseAccelerometers(ABC):
    '''
    Abstract class for parabola
    '''

    @abstractmethod
    def acquireData(self, recording_seconds):
        ''' Function for data acquisition

        Parameters
        ----------
        recording_seconds: int [s]
            recording seconds for data acquisition
        '''
        raise Exception('Implement me!')

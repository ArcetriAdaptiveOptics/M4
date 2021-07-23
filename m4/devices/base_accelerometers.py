'''
Authors
  - C. Selmi: written in 2020
'''
import six
from abc import ABCMeta, abstractmethod

@six.add_metaclass(ABCMeta)
class BaseAccelerometers():
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

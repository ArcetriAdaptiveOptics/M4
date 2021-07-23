'''
Authors
  - C. Selmi: written in 2020
'''
import six
from abc import ABCMeta, abstractmethod

@six.add_metaclass(ABCMeta)
class BaseInterferometer():
    '''
    Abstract class for the interferometer
    '''

    @abstractmethod
    def acquire_phasemap(self, nframes, show):
        ''' Function for data acquisition '''
        raise Exception('Implement me!')

    @abstractmethod
    def save_phasemap(self, location, file_name, masked_image):
        ''' Function for saving data '''
        raise Exception('Implement me!')

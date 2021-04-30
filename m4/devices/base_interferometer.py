'''
Authors
  - C. Selmi: written in 2020
'''
from abc import ABC, abstractmethod


class BaseInterferometer(ABC):
    '''
    Abstract class for the interferometer
    '''

    @abstractmethod
    def acquire_phasemap(self, nframes_or_ott, show):
        ''' Function for data acquisition '''
        raise Exception('Implement me!')

    @abstractmethod
    def save_phasemap(self, location, file_name, masked_image):
        ''' Function for saving data '''
        raise Exception('Implement me!')

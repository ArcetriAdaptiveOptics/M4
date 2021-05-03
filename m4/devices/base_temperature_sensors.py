'''
Authors
  - C. Selmi: written in 2020
'''
from abc import ABC, abstractmethod


class BaseTemperatureSensors(ABC):
    '''
    Abstract class for PT
    '''

    @abstractmethod
    def getTemperature(self):
        ''' Function for getting PT temperature
        '''
        raise Exception('Implement me!')

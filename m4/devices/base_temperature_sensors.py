'''
Authors
  - C. Selmi: written in 2020
'''
import six
from abc import ABCMeta, abstractmethod

@six.add_metaclass(ABCMeta)
class BaseTemperatureSensors():
    '''
    Abstract class for PT
    '''

    @abstractmethod
    def getTemperature(self):
        ''' Function for getting PT temperature
        '''
        raise Exception('Implement me!')

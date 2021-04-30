'''
Authors
  - C. Selmi: written in 2020
'''
from abc import ABC, abstractmethod


class BaseM4(ABC):
    '''
    Abstract class for M4
    '''

    @abstractmethod
    def getPosition(self):
        ''' Function for getting object position
        '''
        raise Exception('Implement me!')

    @abstractmethod
    def setPosition(self, absolute_position_in_mm):
        ''' Function for setting object position

        Parameters
        ----------
        absolute_position_in_mm: int [mm]
        '''
        raise Exception('Implement me!')

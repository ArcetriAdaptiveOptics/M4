'''
Authors
  - C. Selmi: written in 2020
'''
from abc import ABC, abstractmethod


class BaseReferenceMirror(ABC):
    '''
    Abstract class for reference mirror
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

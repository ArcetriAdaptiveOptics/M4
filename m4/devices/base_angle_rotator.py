
'''
Authors
  - C. Selmi: written in 2020
'''
from abc import ABC, abstractmethod


class BaseAngleRotator(ABC):
    '''
    Abstract class for the ring angle rotation
    '''

    @abstractmethod
    def getPosition(self):
        ''' Function for getting object position
        '''
        raise Exception('Implement me!')

    @abstractmethod
    def setPosition(self, absolute_position_in_deg):
        ''' Function for setting object position

        Parameters
        ----------
        absolute_position_in_deg: int [deg]
            rotating ring position
        '''
        raise Exception('Implement me!')

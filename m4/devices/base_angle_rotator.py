
'''
Authors
  - C. Selmi: written in 2020
'''
import six
from abc import ABCMeta, abstractmethod

@six.add_metaclass(ABCMeta)
class BaseAngleRotator():
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

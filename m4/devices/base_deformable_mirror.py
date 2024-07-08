'''
Authors
  - C. Selmi: written in 2022
'''
import six
from abc import ABCMeta, abstractmethod

@six.add_metaclass(ABCMeta)
class BaseDeformableMirror():
    '''
    Abstract class for the deformable mirror
    '''

    def set_shape(self, command):
        ''' 
        Function for setting actuators command

        '''

    def get_shape(self):
        '''
        Function to read the actuators current position

        '''

    def nActuators(self):
        '''
        Function for reading the number of actuators

        '''



    """
    The following lines are the old implementation of the base deformable mirror
    @abstractmethod
    def setActsCommand(self, command):
        ''' Function for setting actuators command '''
        raise Exception('Implement me!')

    @abstractmethod
    def getActsCommand(self):
        ''' Function for getting actuators command '''
        raise Exception('Implement me!')

    @abstractmethod
    def getNActs(self):
        ''' Function for number of actuators '''
        raise Exception('Implement me!')
    """

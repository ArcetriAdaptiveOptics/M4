
'''
Authors
  - C. Selmi: written in 2020
'''

class BaseAngleRotator():
    '''
    Abstract class for the ring angle rotation
    '''

    def getPosition(self):
        raise Exception('Implement me!')

    def setPosition(self, absolute_position_in_mm):
        raise Exception('Implement me!')
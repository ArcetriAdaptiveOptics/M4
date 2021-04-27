'''
Authors
  - C. Selmi: written in 2020
'''

class BaseM4():
    '''
    Abstract class for M4
    '''

    def getPosition(self):
        raise Exception('Implement me!')

    def setPosition(self, absolute_position_in_mm):
        raise Exception('Implement me!')

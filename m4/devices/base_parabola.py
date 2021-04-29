'''
Authors
  - C. Selmi: written in 2020
'''

class BaseParabola():
    '''
    Abstract class for parabola
    '''

    def getPosition(self):
        ''' Function for getting object position
        '''
        raise Exception('Implement me!')

    def setPosition(self, absolute_position_in_mm):
        ''' Function for setting object position

        Parameters
        ----------
        absolute_position_in_mm: int [mm]
        '''
        raise Exception('Implement me!')

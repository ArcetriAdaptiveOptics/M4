
'''
Authors
  - C. Selmi: written in 2020
'''

class BaseAngleRotator():
    '''
    Abstract class for the ring angle rotation
    '''

    def getPosition(self):
        ''' Function for getting object position
        '''
        raise Exception('Implement me!')

    def setPosition(self, absolute_position_in_deg):
        ''' Function for setting object position

        Parameters
        ----------
        absolute_position_in_deg: int [deg]
            rotating ring position
        '''
        raise Exception('Implement me!')

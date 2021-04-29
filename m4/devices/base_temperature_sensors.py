'''
Authors
  - C. Selmi: written in 2020
'''

class BaseTemperatureSensors():
    '''
    Abstract class for PT
    '''

    def getTemperature(self):
        ''' Function for getting PT temperature
        '''
        raise Exception('Implement me!')

'''
Authors
  - C. Selmi: written in 2020
'''

class BaseTemperatureSensors():
    '''
    Abstract class for PT
    '''

    def getTemperature(self):
        raise Exception('Implement me!')

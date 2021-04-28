'''
Authors
  - C. Selmi: written in 2020
'''

class BaseAccelerometers():
    '''
    Abstract class for parabola
    '''

    def acquireData(self, recording_seconds):
        raise Exception('Implement me!')

'''
Authors
  - C. Selmi: written in 2020
'''

class BaseInterferometer():
    '''
    Abstract class for the interferometer
    '''

    def acquire_phasemap(self, nframes_or_ott, show):
        raise Exception('Implement me!')

    def save_phasemap(self, dove, name, image):
        raise Exception('Implement me!')

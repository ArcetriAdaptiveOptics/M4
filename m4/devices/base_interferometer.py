'''
Authors
  - C. Selmi: written in 2020
'''

class BaseInterferometer():
    '''
    Abstract class for the interferometer
    '''

    def acquire_phasemap(self, nframes_or_ott, show):
        ''' Function for data acquisition '''
        raise Exception('Implement me!')

    def save_phasemap(self, location, file_name, masked_image):
        ''' Function for saving data '''
        raise Exception('Implement me!')

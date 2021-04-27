'''
Authors
  - C. Selmi: written in 2020
'''

class BaseReferenceMirror():
    '''
    Abstract class for reference mirror
    '''

    def getPosition(self):
        raise Exception('Implement me!')

    def setPosition(self, absolute_position_in_mm):
        raise Exception('Implement me!')

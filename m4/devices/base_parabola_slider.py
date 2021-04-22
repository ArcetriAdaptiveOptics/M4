

class BaseParabolaSlider():
    '''
    Abstract class for the parabola slider
    '''

    def getPosition(self):
        raise Exception('Implement me!')

    def setPosition(self, absolute_position_in_mm):
        raise Exception('Implement me!')


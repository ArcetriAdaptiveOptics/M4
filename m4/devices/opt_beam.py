'''
Author:
    P. Ferraiuolo - Written in 2024

Description
--------------
High-level , end-user functions to move both the parabola and the reference mirror slider, with respect to the optical aligned centre. both the simulated and real case are handled passing by the configuration file

How To Use
--------------
with the configuration file defined and the ott created

    conf = 'path/to/config.yaml'
    ott, _, _ = start.create_ott(conf)

Stand alone use:

    from m4.devices.opt_beam import Parabola, ReferenceMirror
    par = Parabola(ott, conf)
    rm = ReferenceMirror(ott, conf)
'''
import yaml
from m4.configuration.ott_parameters import OttParameters

class Parabola:
    def __init__(self, ott, conf):
        '''The Constructor'''
        self._slider = ott.parabolaSlider
        self._pos = self._slider.getPosition()

        try:
            with open(conf, 'r') as file:
                self._config = yaml.safe_load(file)
        except Exception as e:
            print(f"Error loading configuration file: {e}")
            self._config = {}

    def _conversion(self, pos: float, get=False) -> float:
        '''
        Converts the given position with relative to M4 Center. Pos is in meters
        '''
        if get==False:
            return (pos + OttParameters.PAR_SLIDER_KIN_OFFSET*1000)
        elif get==True: 
            return (pos - OttParameters.PAR_SLIDER_KIN_OFFSET*1000)

    def trussGetPosition(self) -> float:
        '''
        Returns the current position of the parabola slider in meters

        Returns
        -------
        current_pos : float
            Current position of the parabola slider, in meters
        '''
        self._pos = self._slider.getPosition()
        current_pos = self._conversion(self._pos, get=True)
        return current_pos

    def moveTrussTo(self, pos_in_m: float) -> float:
        '''
        Moves the parabola slider to a given coordinate in meters

        Parameters
        ----------
        pos_in_m : float
            Position, in meters, where to move the parabola

        Returns
        -------
        current_pos : float
            The current position in meters of the parabola relative to M4's center
        '''
        pos_in_mm = pos_in_m * 1000

        if self._config['simulated_parSlider'] == False:
            opcua_pos = self._conversion(pos_in_mm, get=False)
        else: 
            opcua_pos = pos_in_mm

        self._slider.setPosition(opcua_pos)
        current_pos = self._pos = self.trussGetPosition()
        return current_pos

    def moveTrussBy(self, change_in_m: float) -> float:
        '''
        Moves the parabola slider by the specified amount in meters

        Parameters
        ----------
        change_in_m : float
            Change in meters of the position to apply

        Returns
        -------
        current_pos : float
            Current position, in meters, of the parabola slider
        '''
        old_pos = self._slider.getPosition()
        new_pos = old_pos + change_in_m*1000
        self._slider.setPosition(new_pos)
        current_pos = self._pos = self.trussGetPosition()
        return current_pos

    def angleRotatorGetPosition(self) -> float:
        '''
        Gets the current position of the angle rotator in degrees

        Parameters
        ---------------

        Returns
        ---------------
        current_pos : float
            Current position of the angle rotator, in degrees
        '''
        current_pos = self._rotator.getPosition()
        return current_pos

    def angleRotatorSetPositions(self, absolute_deg: float) -> float:
        '''
        Sets the angular position to a desired degree

        Parameters
        ---------------
        absolute_deg : float
            Angular position where to set the parabola

        Returns
        ---------------
        current_pos : float
            Current angular position of the parabola
        '''
        current_pos = self._rotator.setPosition(absolute_deg)
        return self.angleRotatorGetPosition()

    def rotateBy(self, rel_deg) -> float:
        '''
        Rotates the parabola, from the current angular position, by a desired amount

        Parameters
        ---------------
        rel_deg : float
            Change in degrees of the current angular position

        Returns
        ---------------
        current_pos : float
            Current angular position of the parabola
        '''
        old_pos = self.angleRotatorGetPosition()
        new_pos = old_pos + rel_deg
        self.angleRotatorSetPosition(new_pos)
        return self.angleRotatorGetPosition()

class ReferenceMirror:
    def __init__(self, ott, conf):
        '''The Constructor'''
        self._slider = ott.referenceMirrorSlider
        self._pos = self._slider.getPosition()

        try:
            with open(conf, 'r') as file:
                self._config = yaml.safe_load(file)
        except Exception as e:
            print(f"Error loading configuration file: {e}")
            self._config = {}

    def _conversion(self, pos: float, get=False) -> float:
        if get==False:
            return (pos + OttParameters.RM_SLIDER_KIN_OFFSET*1000)
        elif get==True: 
            return (pos - OttParameters.RM_SLIDER_KIN_OFFSET*1000)

    def rmGetPosition(self) -> float:
        '''
        Returns the current position in meters of the reference mirror slider

        Returns
        -------
        current_pos : float
            Current position of the reference mirror slider, in meters
        '''
        self._pos = self._slider.getPosition()
        current_pos = self._conversion(self._pos, get=True)
        return current_pos

    def moveRmTo(self, pos_in_m: float) -> float:
        '''
        Moves the reference mirror slider to a given coordinate in meters

        Parameters
        ----------
        pos_in_m : float
            Position, in meters, where to move the reference mirror

        Returns
        -------
        current_pos : float
            The current position in meters of the reference mirror relative to M4's center
        '''
        pos_in_mm = pos_in_m * 1000

        if self._config['simulated_rmSlider'] == False:
            opcua_pos = self._conversion(pos_in_mm, get=False)
        else:
            opcua_pos = pos_in_mm

        self._slider.setPosition(opcua_pos)
        current_pos = self._pos = self.rmGetPosition()
        return current_pos

    def moveRmBy(self, change_in_m: float) -> float:
        '''
        Moves the reference mirror slider by the specified amount in meters

        Parameters
        ----------
        change_in_m : float
            Change in meters of the position to apply

        Returns
        -------
        current_pos : float
            Current position, in meters, of the reference mirror slider
        '''
        old_pos = self._slider.getPosition()
        new_pos = old_pos + change_in_m*1000
        self._slider.setPosition(new_pos)
        current_pos = self._pos = self.rmGetPosition()
        return current_pos


"""
Author(s):
    - Pietro Ferraiuolo : written in 2024

Description
-----------
High-level, end-user functions to move both the parabola and the reference mirr
or slider, with respect to the optical aligned centre. both the simulated and r
eal case are handled passing by the configuration file

How to Use it
-------------
Once the ott object has been created:
    
    >>> from m4.configuration.start import create_ott
    >>> ott, interf, dm = create_ott()

    >>> from m4.devices.opt_beam import Parabola, ReferenceMirror, AngleRotator
    >>> truss   = Parabola(ott, conf)
    >>> flat    = ReferenceMirror(ott, conf)
    >>> angrot  = AngleRotator(ott, conf)
"""
import numpy as np
from m4.configuration.ott_parameters import OttParameters
from m4.configuration import update_folder_paths as ufp
config = ufp.folders

class Parabola:
    """
    
    """
    def __init__(self, ott):
        """The Constructor"""
        self._par=ott.parabola
        self._slider = ott.parabolaSlider
        self._pos = self._slider.getPosition()

        try:
            if hasattr(config, 'simulated_parSlider'):
                self._config = config.simulated_parSlider
            else:
                raise KeyError("Parameter not found")
        except Exception as e:
            raise e

    def _conversion(self, pos: float, get=False) -> float:
        """
        

        Parameters
        ----------
        pos : float
            DESCRIPTION.
        get : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        float
            DESCRIPTION.

        """
        if get is False:
            new_pos = pos + OttParameters.PAR_SLIDER_KIN_OFFSET*1000
        else:
            new_pos = pos - OttParameters.PAR_SLIDER_KIN_OFFSET*1000
        return new_pos

    def trussGetPosition(self) -> float:
        """
        

        Returns
        -------
        current_pos : float
            DESCRIPTION.

        """
        self._pos = self._slider.getPosition()
        if self._config is False:
            current_pos = self._conversion(self._pos, get=True)
        else:
            current_pos = self._pos
        current_pos /= 1000
        return current_pos

    def moveTrussTo(self, pos_in_m: float) -> float:
        """
        

        Parameters
        ----------
        pos_in_m : float
            DESCRIPTION.

        Returns
        -------
        float
            DESCRIPTION.

        """
        pos_in_mm = pos_in_m * 1000

        if self._config is False:
            opcua_pos = self._conversion(pos_in_mm, get=False)
        else:
            opcua_pos = pos_in_mm

        self._slider.setPosition(opcua_pos)
        current_pos = self._pos = self.trussGetPosition()
        return current_pos

    def moveTrussBy(self, change_in_m: float) -> float:
        """
        Moves the parabola slider by the specified amount in meters

        Parameters
        ----------
        change_in_m : float
            Change in meters of the position to apply

        Returns
        -------
        current_pos : float
            Current position, in meters, of the parabola slider
        """
        old_pos = self._slider.getPosition()
        new_pos = old_pos + change_in_m*1000
        self._slider.setPosition(new_pos)
        current_pos = self._pos = self.trussGetPosition()
        return current_pos

    def parabolaPiston(self, intensity):
        """
        Applies a relative piston command to the parabola. For absolute positio
        n movements, refer to ott.parabola.

        Parameters
        ---------------
        intensity : int
            Relative change in position to apply, in mllimeters, of the parabol
            a piston.

        Returns
        ---------------
        pos : int
            Current position of the piston, in millimeters
        """
        coords = self._par.getPosition()
        coords[2] += intensity
        self._par.setPosition(coords)
        return self._par.getPosition()[2]

    def parabolaTipTilt(self, tt):
        """
        Applies a relative tip/tilt command to the parabola. For asbolute posit
        ioning, refer to ott.parabola.

        Parameters
        ---------------
        intensity : int
            Relative change in position, in arcseconds, of the parabola tip and
            tilt.

        Returns
        ---------------
        pos : int
            Current position of the tip and tilt, in millimeters.
        """
        tt = np.array(tt)
        if tt.shape != (2,):
            raise ValueError(f"The input array has shape {tt.shape}, but shape\
                             (2,) is expected!")
        coords = self._par.getPosition()
        coords[3] += tt[0]
        coords[4] += tt[1]
        self._par.setPosition(coords)

class ReferenceMirror:
    """
    
    """
    def __init__(self, ott):
        """"The Constructor"""
        self._rm = ott.referenceMirror
        self._slider = ott.referenceMirrorSlider
        self._pos = self._slider.getPosition()
        try:
            if hasattr(config, 'simulated_rmSlider'):
                self._config = config.simulated_rmSlider
            else:
                raise KeyError("Parameter not found")
        except Exception as e:
            raise e

    def _conversion(self, pos: float, get=False) -> float:
        """
        

        Parameters
        ----------
        pos : float
            DESCRIPTION.
        get : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        float
            DESCRIPTION.
        """
        if get is False:
            new_pos = pos + OttParameters.RM_SLIDER_KIN_OFFSET*1000
        else:
            new_pos = pos - OttParameters.RM_SLIDER_KIN_OFFSET*1000
        return new_pos

    def rmGetPosition(self) -> float:
        """
        Returns the current position in meters of the reference mirror slider

        Returns
        -------
        current_pos : float
            Current position of the reference mirror slider, in meters
        """
        self._pos = self._slider.getPosition()
        if self._config is False:
            current_pos = self._conversion(self._pos, get=True)
        else:
            current_pos = self._pos
        current_pos /= 1000
        return current_pos

    def moveRmTo(self, pos_in_m: float) -> float:
        """
        Moves the reference mirror slider to a given coordinate in meters

        Parameters
        ----------
        pos_in_m : float
            Position, in meters, where to move the reference mirror

        Returns
        -------
        current_pos : float
            The current position in meters of the reference mirror relative to 
            M4's center
        """
        pos_in_mm = pos_in_m * 1000
        if self._config is False:
            opcua_pos = self._conversion(pos_in_mm, get=False)
        else:
            opcua_pos = pos_in_mm
        self._slider.setPosition(opcua_pos)
        current_pos = self._pos = self.rmGetPosition()
        return current_pos

    def moveRmBy(self, change_in_m: float) -> float:
        """
        Moves the reference mirror slider by the specified amount in meters

        Parameters
        ----------
        change_in_m : float
            Change in meters of the position to apply

        Returns
        -------
        current_pos : float
            Current position, in meters, of the reference mirror slider
        """
        old_pos = self._slider.getPosition()
        new_pos = old_pos + change_in_m*1000
        self._slider.setPosition(new_pos)
        current_pos = self._pos = self.rmGetPosition()
        return current_pos

    def rmPiston(self, intensity: int) -> int:
        """
        Applies a relative piston command to the reference mirror. For absolute
        positioning, refer to ott.referenceMirror.

        Parameters
        ---------------
        intensity : int
            Relative change in position to apply, in millimeters, to the refere
            nce mirror piston.

        Returns
        ---------------
        pos : int
            Current position of the piston, in millimeters
        """
        coords = self._rm.getPosition()
        coords[2] += intensity
        self._rm.setPosition(coords)
        return self._rm.getPosition()[2]

    def rmTipTilt(self, tt):
        """
        Applies a relative tip/tilt command to the reference mirror. For absolu
        te positioning, refer to ott.referenceMirror

        Parameters
        ---------------
        intensity : int
            Relative change in position to apply, in arcseconds, to the referen
            ce mirror tip and tilt.

        Returns
        ---------------
        pos : int
            Current position of the tip and tilt, in millimeters
        """
        tt = np.array(tt)
        if tt.shape != (2,):
            raise ValueError(f"The input array has shape {tt.shape}, but shape\
                             (2,) is expected!")
        coords = self._rm.getPosition()
        coords[3] += tt[0]
        coords[4] += tt[1]
        self._rm.setPosition(coords)
        return self._rm.getPosition()[3:5]

class AngleRotator:
    """
    
    """
    def __init__(self, ott):
        """"The Constructor"""
        self._rotator = ott.angleRotator
        self._pos = ott.angleRotator.getPosition()
        try:
            if hasattr(config, 'simulated_angleRotator'):
                self._config = config.simulated_angleRotator
            else:
                raise KeyError("Parameter not found")
        except Exception as e:
            raise e

    def getPosition(self) -> float:
        """
        Returns the current position of the angle rotator in degrees.

        Returns
        -------
        current_pos : float
            Current position of the angle rotator, in degrees.
        """
        current_pos = self._rotator.getPosition()
        return current_pos

    def setPosition(self, absolute_deg: float) -> float:
        """
        Sets the angular position to a desired degree

        Parameters
        ----------
        absolute_deg : float
            Angular position where to set the parabola

        Returns
        -------
        current_pos : float
            Current angular position of the parabola
        """
        self._rotator.setPosition(absolute_deg)
        return self.getPosition()

    def rotateBy(self, rel_deg) -> float:
        """
        Rotates the parabola, from the current angular position, by a desired 
        amount counter-clockwise

        Parameters
        ---------------
        rel_deg : float
            Relative change, in degrees, of the current angular position

        Returns
        ---------------
        current_pos : float
            Current angular position of the parabola
        """
        old_pos = self.getPosition()
        new_pos = old_pos + rel_deg
        self.setPosition(new_pos)
        return self.getPosition()

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
    >>> par     = Parabola(ott, conf)
    >>> flat    = ReferenceMirror(ott, conf)
    >>> angrot  = AngleRotator(ott, conf)
"""

import numpy as np
from m4.configuration.ott_parameters import OttParameters
from opticalib.core.read_config import load_yaml_config as _lyc

config = _lyc()["SYSTEM"]["simulated.devices"]


class Parabola:
    """
    Class for moving the parabola slider, i.e. the Truss, with respect to M4's
    optical center, and to command the parabola itself, with piston, tip and ti
    lt.

    Methods
    =======
    Truss Methods
    -------------
    trussGetPosition():
        Gets the current position, in meters, of the Truss, relative to M4'
        s optical center.
    moveTrussTo(pos_in_m):
        Sets the position of the Truss at the desired location.
    moveTrussBy(change_in_m):
        Move the Truss by a relative amount from the current position.
    Parabola Methods
    ----------------
    parabolaPiston(intensity):
        Applies a piston mode to the parabola.
    parabolaTipTilt(tt):
        Applies tip and tilt to the parabola.

    More information on the specific functions uses can be found in the functio
    ns documentation.
    """

    def __init__(self, ott: object):
        """The Constructor"""
        self._par = ott.parabola
        self._slider = ott.parabolaSlider
        self._pos = self._slider.getPosition()

        try:
            if "parSlider" in config.keys():
                self._config = config['parSlider']
            else:
                raise KeyError("Parameter not found")
        except Exception as e:
            raise e

    def _conversion(self, pos: float, get: bool = False) -> float:
        """
        Internal function which handles the offsets beetween M4's optical centr
        e and the OPCUA reference frame.

        Parameters
        ----------
        pos : float
            Input, end-user, position (in meters).
        get : boolean, optional
            Option which handles the cases of conversion for the ''getPosition(
            )'' function and the ''setPosition()'' one. The default is False, t
            hat is for the ''setPosition()''.

        Returns
        -------
        new_pos : float
            Position scaled for offsets and converted to millimiters, to be pas
            sed to the OPCUA.
        """
        if get is False:
            new_pos = pos*1000 + OttParameters.PAR_SLIDER_KIN_OFFSET * 1000
        else:
            new_pos = (pos - OttParameters.PAR_SLIDER_KIN_OFFSET * 1000)
        return new_pos

    def trussGetPosition(self) -> float:
        """
        Returns the current position, in metrs, of the Truss (parabola slider)
        relative to M4's centre.

        Returns
        -------
        current_pos : float
            Current Truss position, in meters.
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
        Moves the Truss to a specifide position, in meters, along the slider's
        range

        Parameters
        ----------
        pos_in_m : float
            Position to be reached along the parabola slider, in meters.

        Returns
        -------
        current_pos : float
            The current position of the Truss in meters.
        """
        pos_in_mm = pos_in_m #* 1000

        if self._config is False:
            opcua_pos = self._conversion(pos_in_mm, get=False)
        else:
            opcua_pos = pos_in_mm

        self._slider.setPosition(opcua_pos)
        current_pos = self._pos = self.trussGetPosition()
        return current_pos

    def moveTrussBy(self, change_in_m: float) -> float:
        """
        Moves the Truss along the parabola slider range, by the specified amoun
        t,from the current position.

        Parameters
        ----------
        change_in_m : float
            Truss position shift from current position, in meters.

        Returns
        -------
        current_pos : float
            The current position of the Truss in meters.
        """
        old_pos = self._slider.getPosition()
        new_pos = old_pos + change_in_m * 1000
        self._slider.setPosition(new_pos)
        current_pos = self._pos = self.trussGetPosition()
        return current_pos

    def parabolaPiston(self, intensity: int):
        """
        Applies a relative piston command to the parabola. For absolute positio
        n movements, refer to ''ott.parabola''.

        Parameters
        ---------------
        intensity : int
            Relative (differential) change in position to apply, in mllimeters,
            of the parabola piston.

        Returns
        ---------------
        pos : int
            Current position of the piston, in millimeters
        """
        coords = self._par.getPosition()
        coords[2] += intensity
        self._par.setPosition(coords)
        return self._par.getPosition()[2]

    def parabolaTipTilt(self, tt: int):
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
            raise ValueError(
                f"The input array has shape {tt.shape}, but shape\
                             (2,) is expected!"
            )
        coords = self._par.getPosition()
        coords[3] += tt[0]
        coords[4] += tt[1]
        self._par.setPosition(coords)

    def parabolaGetPosition(self):
        """
        Returns the current position of the parabola in meters.

        Returns
        -------
        current_pos : float
            Current position of the parabola, in meters.
        """
        return self._par.getPosition()


class ReferenceMirror:
    """
    Class for moving the reference mirror slider with respect to M4's optical
    center, and to command the reference mirror itself, with piston, tip and ti
    lt.

    Methods
    =======
    RM Slider Methods
    -----------------
    rmGetPosition():
        Gets the current position, in meters, of the reference mirror slider, r
        elative to M4's optical center.
    moveRmTo(pos_in_m):
        Sets the position of the reference mirror slider at the desired locatio
        n.
    moveRmBy(change_in_m):
        Move the reference mirror slider by a relative amount from the current
        position.
    Reference Mirror Methods
    ------------------------
    rmPiston(intensity):
        Applies a piston mode to the reference mirror.
    rmTipTilt(tt):
        Applies tip and tilt to the reference mirror.

    More information on the specific functions uses can be found in the functio
    ns documentation.
    """

    def __init__(self, ott: object):
        """ "The Constructor"""
        self._rm = ott.referenceMirror
        self._slider = ott.referenceMirrorSlider
        self._pos = self._slider.getPosition()
        try:
            if "rmSlider" in config.keys():
                self._config = config['rmSlider']
            else:
                raise KeyError("Parameter not found")
        except Exception as e:
            raise e

    def _conversion(self, pos: float, get: bool = False) -> float:
        """
        Internal function which handles the offsets beetween M4's optical centr
        e and the OPCUA reference frame.

        Parameters
        ----------
        pos : float
            Input, end-user, position (in meters).
        get : boolean, optional
            Option which handles the cases of conversion for the ''getPosition(
            )'' function and the ''setPosition()'' one. The default is False, t
            hat is for the ''setPosition()''.

        Returns
        -------
        new_pos : float
            Position scaled for offsets and converted to millimiters, to be pas
            sed to the OPCUA.
        """
        if get is False:
            new_pos = pos*1000 + OttParameters.RM_SLIDER_KIN_OFFSET * 1000
        else:
            new_pos = (pos - OttParameters.RM_SLIDER_KIN_OFFSET * 1000)
        return new_pos

    def rmsGetPosition(self) -> float:
        """
        Returns the current position, in meters, of the reference mirror slider.

        Returns
        -------
        current_pos : float
            Current position of the reference mirror slider, in meters.
        """
        self._pos = self._slider.getPosition()
        if self._config is False:
            current_pos = self._conversion(self._pos, get=True)
        else:
            current_pos = self._pos
        current_pos /= 1000
        return current_pos

    def moveRmsTo(self, pos_in_m: float) -> float:
        """
        Moves the reference mirror slider to a given position in the reference
        mirror slider's range.

        Parameters
        ----------
        pos_in_m : float
            Position to be reached along the reference mirror slider's range, i
            n meters.

        Returns
        -------
        current_pos : float
            The current position, in meters, of the reference mirror relative t
            o M4's center.
        """
        pos_in_mm = pos_in_m  # * 1000
        if self._config is False:
            opcua_pos = self._conversion(pos_in_mm, get=False)
        else:
            opcua_pos = pos_in_mm
        self._slider.setPosition(opcua_pos)
        current_pos = self._pos = self.rmGetPosition()
        return current_pos

    def moveRmsBy(self, change_in_m: float) -> float:
        """
        Moves the reference mirror slider from the current position by the spec
        ified amount, in meters.

        Parameters
        ----------
        change_in_m : float
            Reference mirror slider shift from current position, in meters.

        Returns
        -------
        current_pos : float
            Current position, in meters, of the reference mirror slider.
        """
        old_pos = self._slider.getPosition()
        new_pos = old_pos + change_in_m * 1000
        self._slider.setPosition(new_pos)
        current_pos = self._pos = self.rmGetPosition()
        return current_pos

    def rmTipTilt(self, tt: int):
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
            raise ValueError(
                f"The input array has shape {tt.shape}, but shape\
                             (2,) is expected!"
            )
        coords = self._rm.getPosition()
        coords[3] += tt[0]
        coords[4] += tt[1]
        self._rm.setPosition(coords)
        return self._rm.getPosition()[3:5]

    def rmGetPosition(self):
        """
        Returns the current position of the reference mirror in meters.

        Returns
        -------
        current_pos : float
            Current position of the reference mirror, in meters.
        """
        return self._rm.getPosition()


class AngleRotator:
    """
    Class for the control of the exapode angle rotator.

    Methods
    -------
    getPosition():
        Gets the current angular position relative to a set zero position where
        the segment 1 is bottom left # to be checked
    setPosition(absolute_deg):
        Rotate the exapode to reach the desired angular position.
    rotateBy(rel_deg):
        Rotates the exapode, in the counter-clockwise direction, by a specific
        amount

    More information on the specific functions uses can be found in the functio
    ns documentation.
    """

    def __init__(self, ott: object):
        """ "The Constructor"""
        self._rotator = ott.angleRotator
        self._pos = ott.angleRotator.getPosition()
        try:
            if "angleRotator" in config.keys():
                self._config = config['angleRotator']
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

    def rotateBy(self, rel_deg: float) -> float:
        """
        Rotates the exapode, from the current angular position, by a desired
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

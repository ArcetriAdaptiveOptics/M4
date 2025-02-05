"""
Author(s):
----------
Written by:
        - Marco Xompero
        - Chiara Selmi
Modified by
        - Pietro Ferraiuolo, 2024

Desctiption
===========
Initialization script that sets up the interactive python shell to work with th
e OTT, allowing the user for the control of the Parabola, the Reference Mirror 
and the M4 deformable mirror.
All the subsequent information can be accessed through help() in the initialize
d shell.

High-Level Use
--------------
--> par: 
    Parabola slider with respect to optical axis.

        trussGetPosition():
            Get current position of the parabola slider.

        moveTrussTo(pos_in_m):
            Move parabola slider to desired position.

        moveTrussBy(change_in_m):
            Change the current position of the parabola slider by the desired q
            uantity.

        parabolaPiston(intensity):
            Change piston position of the desired quantity, from the starting p
            osition. Takes the 3th coordinate of the full 6D coordinates array 
            as imput.
            Example 
                >>> coord = np.array([0,0,x,-,-,0] --> intensity = x

        parabolaTipTilt(tt):
            Change the tip and tilt of the parabola of the desired quantity, fr
            om the starting position. Takes the 4th and 5th coordinates of the 
            full 6D coordinate array as imput.
            Example 
                >>> coord = np.array([0,0,-,x,y,0]) --> tt = np.array([x,y])

--> angrot:
    Angle rotator class.

        getPosition():
            Returns the current angolar position of the parabola.

        setPosition(absolute_deg):
            Rotates the parabola of the desired angular quantity.

        rotateBy(rel_deg):
            Rotate the parabola from the current position of the desired angle.


-->  rm:
    Reference mirror slider with respect to optical axis.

        rmGetPosition():
            Get current position of the reference mirror slider.

        moveRmTo(pos_in_m):
            Move reference mirror slider to desired position.

        moveRmBy(change_in_m):
            Change the current position of the reference mirror slider by the d
            esired quanty.

        rmPiston(intensity):
            Change the current piston position of the desiredamount. Takes the 
            3th coordinate of the full 6D coordinate array as imput.
            Example
                >>> coord = np.array([0,0,x,-,-,0]) --> intensity = x

        rmTipTilt(tt):
            Change the current tip and tilt position of the parabola by the des
            ired quantity. Takes the 4th and 5th coordinates of the full 6D coo
            rdinate array as imput.
            Example
                >>> coord = np.array([0,0,-,x,y,0]) --> tt = np.array([x,y])


--> interf:
    Class controlling the interferometer, containing all the commands documente
    d by 4D.

        acquire_phasemap(n_frames=1, delay=0):
            Acquisition of n_frames images with a delay between one another.

        save_phasemapr(location, file_name, masked_ima):
            Save the acquired phasemap masked_ima into a fits file into the giv
            en location

        capture(numberOfFrames, folder_name=None):
            Function for burst acquisition of numberOfFrames images, to be save
            d into folder_name. If no folder name is given, a tracking number i
            s generated and added at the standard capture folder for the 4D.

        produce(folder_name):
            Function for burst conversion of raw files into .4D files. This inc
            ludes transferring the data from the interferometer WS to the tower
            WS in the OPDimages folder.

        getCameraSettings():
            Reads the camera setting from the 4D configuration file of the inte
            rferometer.

        getFrameRate():
            Gets the frame rate of the interferometer from the 4D configuration
            file.

        intoFullFrame(img):
            Takes the cropped acquisition and returns it into the full 2048x2048
            interferometer frame, after reading the cropping parameters.


--> dm:
    Class which contains the information and the functions to control the defor
    mable mirror

        setActsCommand(command, rel=True):
            Set the given command array to the dm actuators. If rel=True, then 
            relative commands are applied, else absolute.

        getActsCommand():
            Gets the current commands applied to the actuators.

        getNActs():
            Gets the number of actuators of the dm.

        setIncreaseToSegActs(inc, rel=False):
            Applies a piston command to all actuators.

        setZeroToSegActs():
            Resets the piston position for all actuators.

        setRandomCommandToSegActs(amp, rel=False):
            Generates a set of random distribution of piston between 0 and amp.

Low-Level Use
-------------
Functions that talk directly to the OPCUA, thus relative to its reference frame.

    ott.parabolaSlider :
        getPosition()
        
        setPosition(pos)

    ott.parabola:
        getPosition()
        
        setPosition(coords)

    ott.referenceMirrorSlider:
        getPosition()
        
        setPosition(pos)

    ott.referenceMirror:
        getPosition()
        
        setPosition(coords)

    ott.angleRotator:
        getPosition()
        
        setPosition(deg)
"""
from m4 import main, noise
import os, time, numpy as np
from m4.devices.i4d import I4D
from m4.utils import osutils as osu
from m4.mini_OTT import measurements
from matplotlib import pyplot as plt
from m4.userscripts import OTTScripts
from m4.analyzers import timehistory as th
from m4.configuration.start import create_ott
from m4.ground import read_data as rd, zernike as zern
from m4.configuration import update_folder_paths as ufp
from m4.configuration.ott_parameters import Interferometer
from m4.devices.opt_beam import Parabola, ReferenceMirror, AngleRotator

ott, interf, dm = create_ott()
fn = ufp.folders
par = Parabola(ott)
flat = ReferenceMirror(ott)
angrot = AngleRotator(ott)
meas = measurements.Measurements(ott, interf)
phcamfocus = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)
myott = OTTScripts(ott, interf, dm)

text = \
"""Using the IPython console for OTT operations.

           |X|           |X|
           |X|_____ _____|X|
           |X|           |X|
   O.T.T.  |X|_____--___<|X|
           |X|    |\  /| |X|
     __    |X|    |_\/_| |X|
    |  |   |X|    | /\ | |X|
    |  |   |X|    |/  \| |X|
 ___|__|___|X|____ ---- _|X|_____

Type help() for information on the available operations."""
print(text)

def docs(arg=None):
    """
    Prints the software documentation.

    Parameters
    ----------
    arg : str, optional
        The specified object you want info about, without the need to search
        it in the full documentation. If no argument is given, the full documentation
        will be printed.
        Args:
            - par
            - flat
            - angrot
            - interf
            - dm
            - ott

    Returns
    -------
    The documentation.
    """
    INIT = [
        "    OTT  DOCUMENTATION",
        "",
        "     M       M         44",
        "     M M   M M        4 4",
        "     M  M M  M       4  4",
        "     M   M   M      4   4",
        "     M       M     4 4 4 4 4",
        "     M       M          4",
        "     M       M          4",
        "     M       M          4",
    ]

    HIGH_LEVEL_USAGE = """
The OTT is initialized.

HIGH-LEVEL USAGE:

--> par : parabola slider with respect to optical axis.

    Available functions:

        trussGetPosition() : Get current position of the parabola slider.

        moveTrussTo(pos_in_m) : Move parabola slider to desired position.

        moveTrussBy(change_in_m) : Change the current position of the parabola 
        slider by the desired quantity.

        parabolaPiston(intensity) : Change piston position of the desired 
        quantity, from the starting position. Takes the 3rd coordinate of the 
        full 6D coordinates array as input.
            Example: coord = np.array([0,0,x,-,-,0]) --> intensity = x

        parabolaTipTilt(tt) : Change the tip and tilt of the parabola of the 
        desired quantity, from the starting position. Takes the 4th and 5th 
        coordinates of the full 6D coordinate array as input.
            Example: coord = np.array([0,0,-,x,y,0]) --> tt = np.array([x,y])
    -------------------

--> angrot : Angle rotator class.

    Available functions:

        getPosition() : Returns the current angular position of the parabola.

        setPosition(absolute_deg) : Rotates the parabola to the desired angular 
        position.

        rotateBy(rel_deg) : Rotate the parabola from the current position by 
        the desired angle.

    -------------------
--> flat : Reference mirror slider with respect to optical axis.

    Available functions:

        rmGetPosition() : Get current position of the reference mirror slider.

        moveRmTo(pos_in_m) : Move reference mirror slider to desired position.

        moveRmBy(change_in_m) : Change the current position of the reference 
        mirror slider by the desired quantity.

        rmPiston(intensity) : Change the current piston position by the desired 
        amount. Takes the 3rd coordinate of the full 6D coordinate array as 
        input.
            Example: coord = np.array([0,0,x,-,-,0]) --> intensity = x

        rmTipTilt(tt) : Change the current tip and tilt position of the 
        parabola by the desired quantity. Takes the 4th and 5th coordinates of 
        the full 6D coordinate array as input.
            Example: coord = np.array([0,0,-,x,y,0]) --> tt = np.array([x,y])

    -------------------

--> interf : Class controlling the interferometer, containing all the commands 
documented by 4D.

    Available functions:

        acquire_phasemap(n_frames=1, delay=0) : Acquisition of n_frames images 
        with a delay between one another.

        save_phasemap(location, file_name, masked_ima) : Save the acquired 
        phasemap masked_ima into a fits file at the given location.

        capture(numberOfFrames, folder_name=None) : Function for burst 
        acquisition of numberOfFrames images, to be saved into folder_name. If 
        no folder name is given, a tracking number is generated and added to 
        the standard capture folder for the 4D.

        produce(folder_name) : Function for burst conversion of raw files into 
        .4D files. This includes transferring the data from the interferometer 
        WS to the tower WS in the OPDimages folder.

        getCameraSettings() : Reads the camera settings from the 4D 
        configuration file of the interferometer.

        getFrameRate() : Gets the frame rate of the interferometer from the 4D 
        configuration file.

        intoFullFrame(img) : Takes the cropped acquisition and returns it into 
        the full 2048x2048 interferometer frame, after reading the cropping 
        parameters.

    -------------------

--> dm : Class controlling the actuators on the deformable mirror.

    Available functions:

        setActsCommand(command, rel=True) : Set the given command array to the 
        dm actuators. If rel=True, then relative commands are applied, else 
        absolute.

        getActsCommand() : Gets the current commands applied to the actuators.

        getNActs() : Gets the number of actuators of the dm.

        setIncreaseToSegActs(inc, rel=False) : Applies a piston command to all 
        actuators.

        setZeroToSegActs() : Resets the piston position for all actuators.

        setRandomCommandToSegActs(amp, rel=False) : Generates a set of random 
        distribution of piston between 0 and amp.

    -------------------
"""

    LOW_LEVEL_USAGE = """
LOW LEVEL USAGE
Classes that talk directly to the OPCUA, thus relative to its reference frame.

    ott.parabolaSlider :
        getPosition()
        setPosition(pos)

    ott.parabola:
        getPosition()
        setPosition(coords)

    ott.referenceMirrorSlider:
        getPosition()
        setPosition(pos)

    ott.referenceMirror:
        getPosition()
        setPosition(coords)

    ott.angleRotator:
        getPosition()
        setPosition(deg)
    """

    if arg == 'par':
        print(HIGH_LEVEL_USAGE.split('-------------------')[0])
    elif arg == 'angrot':
        print(HIGH_LEVEL_USAGE.split('-------------------')[1])
    elif arg == 'flat':
        print(HIGH_LEVEL_USAGE.split('-------------------')[2])
    elif arg == 'interf':
        print(HIGH_LEVEL_USAGE.split('-------------------')[3])
    elif arg == 'dm':
        print(HIGH_LEVEL_USAGE.split('-------------------')[4])
    elif arg == 'ott':
        print(LOW_LEVEL_USAGE)
    else:
        for line in INIT:
            print(line)
        print(HIGH_LEVEL_USAGE)
        print(LOW_LEVEL_USAGE)


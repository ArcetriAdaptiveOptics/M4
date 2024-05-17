import numpy as np
from matplotlib import pyplot as plt
#from m4.misc import par_meas as pp
from m4.ground import read_data as rr
from m4.mini_OTT import timehistory as th
from m4.ground import zernike as zern
from m4.mini_OTT.measurements import Measurements
from m4 import main
from astropy.io import fits as pyfits
from m4.configuration import start
from m4.devices.i4d import I4D
from m4.devices.opt_beam import Parabola, ReferenceMirror
from m4.configuration.ott_parameters import Interferometer
from m4 import noise
import time
#conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
conf = os.path.join(os.path.expanduser('~'), 'git/M4/m4/configuration','myConfig.yaml')
ott, interf, dm = start.create_ott(conf)
truss = Parabola(ott, conf) # could name 'optical_beam'
rm = ReferenceMirror(ott, conf)
meas = Measurements(ott,interf)
phcamfocus = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)

print("\n Type help() for information on the available operations")

def help():
    INIT = [
            "OTT INITIALIZATION DOCUMENTATION",
            "",
            "     M       M        4 4",
            "     M M   M M       4  4",
            "     M  M M  M      4   4",
            "     M   M   M     4 4 4 4 4",
            "     M       M          4",
            "     M       M          4",
            "     M       M          4",
            ]


TEXT = """
The (Simulated) OTT has been initialized.

HIGH-LEVEL USAGE:

    truss : parabola slider with respect to optical axis.

    Available functions:

        trussGetPosition() : Get current position of the parabola slider.

        moveTrussTo(pos_in_m) : Move parabola slider to desired position.

        moveTrussBy(change_in_m) : Change the current position of the parabola slider by the desired quantity.

        parabolaSetPosition(coord) : Change piston position of the desired quantity. Takes the 6-D coordinates array as imput, where the piston/focus coordinate must be the third.
            Example: coord = np.array([0,0,x,-,-,0]

        parabolaTilt(coords) : Change tilt in the parabola. Takes the 6-D coordinte array as imput, where the tilt coordinates must be the fourth and the fifth.
            Example: coord = np.array([0,0,-,x,x,0])

        angleRotatorGetPosition() : Returns the current angolar position of the parabola.

        angleRotatorSetPosition(absolute_deg) : Rotates the parabola of the desired angular quantity.

        rotateBy(rel_deg) : Rotate the parabola from the current position of the desired angle.

    -------------------
 rm : Reference mirror slider with respect to optical axis.

    Available functions:

        rmGetPosition() : Get current position of the reference mirror slider.

        moveRmTo(pos_in_m) :  Move reference mirror slider to desired position.

        moveRmBy(change_in_m) : Change the current position of the reference mirror slider by the desired quanty.

    -------------------

    interf : Class controlling the interferometer, containing all the commands documented by 4D.

    Available functions:

        acquire_phasemap(n_frames=1, delay=0) : Acquisition of n_frames images with a delay between one another

        save_phasemapr(location, file_name, masked_ima) : Save the acquired phasemap masked_ima into a fits file into the given location

        capture(numberOfFrames, folder_name=None) : Function for burst acquisition of numberOfFrames images, to be saved into folder_name. If no folder name is given, a tracking number is generated and added at the standard capture folder for the 4D.

        produce(folder_name) : Function for burst conversion of raw files into .4D files. This includes transferring the data from the interferometer WS to the tower WS in the OPDimages folder.

        getCameraSettings() : Reads the camera setting from the 4D configuration file of the interferometer.

        getFrameRate() : Gets the frame rate of the interferometer from the 4D configuration file.

        intoFullFrame(img) : Takes the cropped acquisition and returns it into the full 2048x2048 interferometer frame, after reading the cropping parameters.
         dm : Class containing all the Deformable Mirror function able to controll it.

    Available functions:

        setActsCommand(command, rel=True) : Set the given command array to the dm actuators. If rel=True, then relative commands are applied, else absolute.

        getActsCommand() : Gets the current commands applied to the actuators.

        getNActs() : Gets the number of actuators of the dm.

        setIncreaseTOSegActs(inc, rel=False) : Applies a piston command to all actuators.

        setZeroToSegActs() : Resets the piston position for all actuators.

        setRandomCommandToSegActs(amp, rel=False) : Generates a set of random distribution of piston between 0 and amp.
 -------------------

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

    for line in INIT:
        print(line)

    print(TEXT)


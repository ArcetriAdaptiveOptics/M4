"""
Alignment Software Configuration File
=====================================
Author(s)
---------
    - Pietro Ferraiuolo : written in 2024
        e-mail : pietro.ferraiuolo@inaf.it

Description
-----------
This python file is a configuration file to be passed to the main alignment code.
This file contains all the informations about the kinematic motors that move the
system and how to acces them through higher level calls.

Variables Explanation
---------------------
cmdDof : list of int | int
    This variable represents the total degrees of freedom each device has, so
    it will be the length of the accepted command vector by the devices.
    Example:
        Normally a device could move up to 6 degreed of freedom, 3 linear and 3
        angular, so cmdDof would be 6.
        If a device accepts a vector of different size for its actuation, then 
        this variable should be a list of integers, each would indicate the dimention
        of the vector corresponding to the respective device.

dof : list of lists
    This variable, initialized as an empty list, will contain the degrees of
    freedom each device can actually move.
    Example:
        If a device can move only 3 DoF, say piston, tip and tilt, then, if the
        cmdDof if 6 (i.e the device accept a vector of 6 elements), the 'dof' list
        should be the index at which these degrees of freedom are located in the
        accepted command vector. In this case, the list could be [2, 3, 4] (as 
        in the case for the M4's OTT).

slices : list of slices
    This variable is a list of slices that will be used to extract the right dof
    positions from the command matrix full vector. 
    Example:
        If the command matrix is a column vector (the full command) of 8 elements,
        in which the first 3 elements are for the 'device_1' dof, the second 2 
        for 'device_2' dof and the last 3 for 'device_3' dof, then the slices 
        list should be [slice(0,3), slice(3,5), slice(5,8)].

zernike_to_use : list of int
    This variable contains the zernike modes indices which will be used for the 
    Interaction Matrix Construction, extracted from the general zernike fit.

push_pull_template : list of int
    Template for the push-pull algorithm. Being the command differential, the 
    push and pulls are to be written with respect to zero.
    Example:
        A classic [1,-1,1] template, translates into [1,-2,1].

devices_move_calls : list of str
    This list must contain the callables for moving each device, in the same 
    order as all the othe lists.
    Example:
        following the previus examples, here the order for the methods would be 
        device_1, device_2 and then device_3 functions.

devices_read_calls : list of str
    This list must contain the callables used to read the device positions, in 
    the same order as all the othe lists.

names : list of str
    This list contains the names of the devices, in the same order as all the 
    other lists. This is used for logging purposes and print fancyness.

ccd_acquisition : list of str
    This list contains the callables for the acquisition CCD/device. The first
    entry must be the acquisition method.

base_read_data_path : str
    Path to the data directory where the data is read.

base_write_data_path : str
    Path to the data directory where the produced data is stored.

log_path : str
    Path to the log file, file name included.

logging_level : int
    Logging level for the logger.
    Warning = 30, Info = 20, Debug = 10, Notset = 0

commandMatrix : str
    Path to the actuation command matrix file.

calibration : str
    Path to the calibration fits file, whose mask will be use for zernike fitting
    purposes. If no mask is to be provided, leave it as an empty string.

Important Notes
---------------
The element order of the lists is very important, as the code will use the same index
for all the lists to access the right information. The order of the devices in the
lists must be the same as the order in the command matrix supplied, eg. if the device_y 
is the second in place going through the command vector, then it's position must be
the secondo for every list in this configuration file.
It's easier to see, as the following example will show.
"""
# M4's specific data path
from m4.configuration import update_folder_paths as ufp
_base_path = ufp.folders.BASE_PATH

# Variables
cmdDof = 6                # Total DoF per device
        # or 
        # cmdDof = []
        # cmdDof.append(6)
        # ...
dof = []                  # Available Degrees of Freedom (DoF)
dof.append([2, 3, 4])     # Parabola DoF
dof.append([3, 4])        # Reference Mirror DoF
dof.append([3, 4])        # M4 Exapode DoF
slices = []               # Full cmd vector devices indices
slices.append(slice(0,3)) # Parabola
slices.append(slice(3,5)) # Reference Mirror
slices.append(slice(5,7)) # M4 Exapode

zernike_to_use = [1,2,3,6,7]
push_pull_template = [+1,-2,+1]

# Devices calls
devices_move_calls = [
    'parabola.setPosition',
    'referenceMirror.setPosition',
    'm4Exapode.setPosition'
    ]

devices_read_calls = [
    'parabola.getPosition',
    'referenceMirror.getPosition',
    'm4Exapode.getPosition'
    ]

names = [
    'Parabola',
    'Reference Mirror',
    'M4 Exapode'
    ]

ccd_acquisition = [
    'acquire_phasemap'
    ]

# Data paths
base_read_data_path     = _base_path+'/M4Data/OPTData/AlignmentCalibration'
base_write_data_path    = _base_path+'/M4Data/OPTData/AlignmentCalibration'
log_path                = _base_path+'/M4Data/OPTData/AlignmentCalibration/alignment.log'
logging_level           =  20
commandMatrix           = _base_path.strip('/data')+'/configuration/cmdMat.fits'
calibrated_parabola     = ''

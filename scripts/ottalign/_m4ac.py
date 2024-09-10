"""
Author(s):
    - Pietro Ferraiuolo : written in 2024

M4's Alignment Configuration File
=================================
Description
-----------
This python file is a configuration file to be passed to the alignment procedure,
which contains all the informations about the kinematic motors that move the
system and how to acces them through the higher level library.

How to Use it
-------------
Check the 'opt_alignment.py' library to see how it's used
"""
import numpy as np
#_____________________________________________________________________________#
"""
Mechanical motors - available degrees of freedom
'dof' is a list containing the arrays of available degree of fredom for each 
      ordered as the full cmd vector is
'cmdDof' is the total degrees of fredom the device accept as a command
'slices' is a list in which each element is the slice relative to the corresponding
      device in the full cmd vector
"""
dof    = [np.array([2, 3, 4]), np.array([3, 4]), np.array([3, 4])]
cmdDof = 6
slices = [slice(0,3), slice(3,5), slice(5,7)]
#_____________________________________________________________________________#
"""Functions to actuate the motors
This variable should be a list of ordered functions, i.e the first function 
commands the first device in full command vector, the second one controls the
second device in the full command vector and so on..."""

dev_calls = [
    'parabola.setPosition',
    'referenceMirror.setPosition',
    'm4Exapode.setPosition'
         ]
ccd_calls = [
    'acquire_phasemap'
    ]
"""
Names of the devices, for print fancyness. Ordered as usual.
"""
names = [
    'Parabola',
    'Reference Mirror',
    'M4 Exapode'
    ]
#_____________________________________________________________________________#

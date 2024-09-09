"""
Author(s):
    - Pietro Ferraiuolo : written in 2024

Description
-----------

How to Use it
-------------
"""
import numpy as np

# Mechanical motors - available degrees of freedom
dof    = [np.array([2, 3, 4]), np.array([3, 4]), np.array([3, 4])]
cmdDof = 6
slices = [slice(0,3), slice(3,5), slice(5,7)]

# Functions to actuate the motors
ottpar = ott.parabola.setPosition
ottrm  = ott.referenceMirror.setPosition
ottm4  = ott.dm.setPosition

"""
Kernel startup script of the M4 GUI.

Run by the GUI in the kernel after opticalib's ``initCalpy.py`` (which has
imported opticalib with the experiment configuration), and again after every
kernel restart.  Keep it light: it must not connect to the hardware, which
is done from the procedure windows.
"""

import os as _os
import sys as _sys

import m4

# userscripts.py imports ott_measurements as a top-level module.
_m4_dir = _os.path.dirname(m4.__file__)
if _m4_dir not in _sys.path:
    _sys.path.append(_m4_dir)

del _os, _sys, _m4_dir
print("M4 environment loaded: use the 'OTT user scripts' procedure to set up the OTT.")

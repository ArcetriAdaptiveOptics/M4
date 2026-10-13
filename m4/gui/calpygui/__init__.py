"""
M4 CalpyGUI
===========

The opticalib GUI (:class:`opticalib.gui.CalpyGUI`) with the M4 procedures.

How it is built
---------------
opticalib's GUI is customised by subclassing it (see :mod:`.app`):

* :attr:`~opticalib.gui.CalpyGUI.PROCEDURES` lists the cards of the
  *Procedures* panel; each :class:`opticalib.gui.Procedure` opens a window.
* A window (see :mod:`.procedures`) is a
  :class:`~opticalib.gui.procedures.base.ProcedureWindow`: a list of
  :class:`~opticalib.gui.procedures.base.Step` objects, each with a form of
  parameters and a code template.  Running a step runs its code, visibly,
  in the GUI console, so everything stays reproducible from the shell.
* :attr:`~opticalib.gui.CalpyGUI.INIT_SCRIPTS` are run in the kernel after
  opticalib's ``initCalpy.py`` (see ``initM4GUI.py``): they prepare the
  namespace the step code relies on (here, ``m4``).

Launch it with ``m4gui -f <experiment>``.

Adding a function of ``userscripts.py``
---------------------------------------
Add a :class:`~opticalib.gui.procedures.base.Step` to
:meth:`.procedures.OttScriptsWindow.steps`: give it a key, a title, a
description, its parameters (``IntParam``, ``BoolParam``, ``ExprParam``,
...) and a ``lambda v: "..."`` returning the code, where ``v`` holds the
Python source of each parameter.  Add ``confirm="..."`` when it moves
hardware.  A new card is a new ``ProcedureWindow`` subclass listed in
:attr:`.app.M4GUI.PROCEDURES`.
"""

from .app import M4GUI, launch_gui

__all__ = ["M4GUI", "launch_gui"]

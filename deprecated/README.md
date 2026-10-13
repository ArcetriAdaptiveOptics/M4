# Deprecated modules

The modules listed below were retired from the `m4` package during the tidy-up.
They are kept here for reference only: they are not installed with the package
and are not importable as `m4` modules. Each file lives at `deprecated/<original path>`
and was moved with `git mv`, so its history is preserved.

| original path | reason | replacement (if any) |
|---|---|---|
| m4/caliball_average.py | Imports modules deleted in earlier purges (m4.ground.zernike/geo/read_data/timestamp/tracking_number_folder, m4.analyzers.timehistory, m4.main, m4.configuration.update_folder_paths) | none |
| m4/howto_OTT.py | Imports modules deleted in earlier purges (see above) | none |
| m4/howto_OTT2.py | Does not parse (SyntaxError) | none |
| m4/requirements_checker.py | Imports modules deleted in earlier purges (see above); also uses removed opticalib.ground.geo | none |
| m4/zernikeCommandTest.py | Imports modules deleted in earlier purges (see above) | none |
| m4/analyzers/executor.py | No users anywhere in the repo | none |
| m4/devices/m4_exapode.py | Superseded by devices/eelt_exapode | m4/devices/eelt_exapode.py |
| m4/devices/opc_ua_controller_old.py | Superseded by devices/opc_ua_controller | m4/devices/opc_ua_controller.py |
| m4/ground/ottInformation.py | Uses opticalib.ground.logger.txtLogger, removed from opticalib | none |
| m4/ground/ott_configurations.py | Unused; docstring refers to non-existent m4.configuration.start | none |
| m4/ground/roi.py | Replaced by opticalib.ground.roi | opticalib.ground.roi |
| m4/ground/zabbix.py | Unused; opens an OPC-UA connection at import time | none |
| m4/gui/accelerometers_GUI.py | guietta GUI, superseded by m4.gui.calpygui (opticalib CalpyGUI) | m4/gui/calpygui |
| m4/gui/geometry_GUI.py | guietta GUI, superseded by m4.gui.calpygui; fails to import (PYOTTCONF env var, m4.ground.package_data, opticalib.ground.zernike) | m4/gui/calpygui |
| m4/gui/geometry_GUI_control.py | guietta GUI (part of geometry GUI), superseded by m4.gui.calpygui; fails to import (see geometry_GUI.py) | m4/gui/calpygui |
| m4/gui/geometry_GUI_new.py | guietta GUI (part of geometry GUI), superseded by m4.gui.calpygui; fails to import (see geometry_GUI.py) | m4/gui/calpygui |
| m4/gui/geometry_GUI_view.py | guietta GUI (part of geometry GUI), superseded by m4.gui.calpygui; fails to import (see geometry_GUI.py) | m4/gui/calpygui |
| m4/gui/monitor_GUI.py | guietta GUI, superseded by m4.gui.calpygui; fails to import (PYOTTCONF env var, m4.ground.package_data, opticalib.ground.zernike) | m4/gui/calpygui |
| m4/gui/temperature_GUI.py | guietta GUI, superseded by m4.gui.calpygui; fails to import (PYOTTCONF env var, m4.ground.package_data, opticalib.ground.zernike) | m4/gui/calpygui |
| m4/mini_OTT/__init__.py | Part of the mini_OTT package; imports modules deleted in earlier purges | none |
| m4/mini_OTT/analysis.py | Part of the mini_OTT package; imports modules deleted in earlier purges | none |
| m4/mini_OTT/measurements.py | Superseded by m4.ott_measurements | m4/ott_measurements.py |
| m4/mini_OTT/reqanalysis.py | Part of the mini_OTT package; imports modules deleted in earlier purges | none |
| m4/mini_OTT/spider_test.py | Part of the mini_OTT package; imports modules deleted in earlier purges | none |
| m4/mini_OTT/timehistory.py | Superseded by opticalib.analyzer timeseries functions | opticalib.analyzer (timeseries functions) |
| m4/simulator/actuator_positions_generator.py | No users anywhere in the repo | none |
| m4/simulator/cap_sens_calib_simulator.py | Imports modules deleted in earlier purges (see above) | none |
| m4/type/measure.py | Imports modules deleted in earlier purges (see above) | none |
| m4/type/temperature_sensors.py | No users anywhere (temperature data handled by devices/temperature_sensors) | m4/devices/temperature_sensors.py |
| m4/utils/PistonMeasures.py | Imports modules deleted in earlier purges (see above) | none |
| m4/utils/old_image_registration_lib.py | Imports modules deleted in earlier purges (see above) | none |
| m4/utils/read_amos_xyz.py | Uses arte.utils.surface_fitter, which no longer exists | none |
| m4/utils/rotation_and_optical_axis_alignment.py | Imports modules deleted in earlier purges (see above) | none |
| m4/utils/shell_actuator_footprint_registration.py | Uses photutils.centroids.fit_2dgaussian, removed from photutils | none |
| m4/utils/snapshooter.py | Python-2 era code (Queue, pyfits, argos) | none |
| m4/utils/snapshot_tweaker.py | Python-2 era code (Queue, pyfits, argos) | none |
| m4/utils/snapshotable.py | Python-2 era code (Queue, pyfits, argos) | none |
| m4/utils/parallelAcq.py | No users anywhere in the repo | none |
| test/analyzers/test_compute_reconstructor.py | Tests m4.analyzers.compute_reconstructor, which does not exist in the repo | none |
| test/ground/test_timestamp.py | Tests m4.ground.timestamp, deleted in an earlier purge | none |
| test/type/test_measure.py | Tests m4.type.measure, moved to deprecated | none |
| test/ground/test_zabbix.py | Tests m4.ground.zabbix, moved to deprecated | none |
| test/ground/test_groundSimpleFunctions.py | Tests m4.ground.geo, smooth_function, logger_set_up, tracking_number_folder, zernike and read_data, none of which exist | none |
| test/utils/test_roi.py | Tests m4.ground.roi (body commented out; imports the removed m4.ground.read_data) | opticalib.ground.roi |
| test/devices/test_i4d.py | Tests m4.devices.i4d, which does not exist in the repo | none |
| test/ground/test_ott_configurations.py | Tests m4.ground.ott_configurations (deprecated) | none |
| test/utils/test_opticalAlignment.py | Tests m4.utils.optical_alignment, which does not exist in the repo | m4.alignment.OttAligner |
| __init_scripts__/pyott.py | primordial calpy launcher, superseded by opticalib `calpy` | `calpy -f path/to/m4` (opticalib) |
| __init_scripts__/initOTT.py | IPython start-up script for pyott; superseded by calpy (+ m4.gui.calpygui.initM4GUI for the GUI) | `calpy -f path/to/m4` (opticalib) |
| __init_scripts__/initAnalysis.py | IPython start-up script for pyott; superseded by calpy (+ m4.gui.calpygui.initM4GUI for the GUI) | `calpy -f path/to/m4` (opticalib) |
| m4/spl_meas.py | superseded by opticalib.procedures.phasing | opticalib.procedures.phasing |
| m4/spl_meas2.py | superseded by opticalib.procedures.phasing | opticalib.procedures.phasing |

To resurrect a module, `git mv` it back to its original path (e.g.
`git mv deprecated/m4/utils/snapshooter.py m4/utils/snapshooter.py`) and fix its imports.

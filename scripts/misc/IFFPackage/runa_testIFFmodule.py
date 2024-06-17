from importlib import reload


from m4.devices import deformable_mirror as dm
from scripts.misc.IFFPackage import iff_acquisition_preparation as ifa
m4u = dm.M4AU()

ifc=ifa.IFFCapturePreparation(m4u)

ifc._updateModalBase('zonal')
t = ifc._createTriggerPadding()
r = ifc._createRegistrationPattern()

mlist = arange(0,5)
mm = ifc._createCmdMatrix(mlist)
tm = ifc.createCmdMatrixHistory(mlist, 1e-7)

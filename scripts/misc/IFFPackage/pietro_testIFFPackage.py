# Test IFF Package Pietro
# run '/home/pietrof/git/M4/m4/initOTT.py'
run '/home/pietrof/git/M4/m4/analysis_imports.py'
from m4.devices import deformable_mirror as dm
from scripts.misc.IFFPackage import iff_acquisition_preparation as ifa
from scripts.misc.IFFPackage import iff_processing as ifp

# Matrices creation and visualization
ifa     = ifa.IFFCapturePreparation(dm.M4AU())
ifa._updateModalBase('zonal')
amp     = 1e-7
nmodes  = 5
mlist   = np.arange(0, nmodes, 1)

t       = ifa._createTriggerPadding()
r       = ifa._createRegistrationPattern()
cmd     = ifa._createCmdMatrix(mlist)
cmdh    = ifa.createCmdMatrixHistory(mlist)
tcmdh   = ifa.createTimedCmdHistory(mlist,amp)

imshow(t)
imshow(r)
imshow(cmd)
imshow(cmdh)
imshow(tcmdh)

# Processing verification
registrationList, registrationFrames = ifp.getTriggerAndRegistrationFrames()
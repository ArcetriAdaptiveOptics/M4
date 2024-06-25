# Test IFF Package Pietro
run '/home/pietrof/git/M4/m4/initOTT.py'
run '/home/pietrof/git/M4/m4/analysis_imports.py'
from m4.devices import deformable_mirror as dfm
from scripts.misc.IFFPackage import iff_acquisition_preparation as ifa
from scripts.misc.IFFPackage import iff_processing_pietro as ifp
from m4.configuration import config_folder_names as fn
from m4.ground import read_data as rd
from m4.configuration import read_iffconfig
from m4.mini_OTT import timehistory as th
m4 = dfm.M4AU()

# ciclo per creare dati simulati
# Si assume che si sia runnato pyott -i ( o run inittOTT.py)

# Matrices creation and visualization
ifa     = ifa.IFFCapturePreparation(dfm.M4AU())
ifa._updateModalBase('zonal')
amp     = 1e-7
nModes  = 10
mlist   = np.arange(0, nmodes, 1)
t       = ifa._createTriggerPadding();imshow(t);colorbar()
r       = ifa._createRegistrationPattern();imshow(r);colorbar()
cmd     = ifa._createCmdMatrix(mlist);imshow(cmd);colorbar()
cmdh    = ifa.createCmdMatrixHistory(mlist);imshow(cmdh);colorbar()
tcmdh   = ifa.createTimedCmdHistory(mlist,amp);imshow(tcmdh);colorbar()

# Processing verification
tn                      = '20160516_114916'
parameters              = ifp._getAcqPar(tn)                # v
infoT, infoR, infoIF    = ifp._getAcqInfo()                 # v
filelist                = ifp._getFileList(tn)              # v
trigF                   = ifp.getTriggerFrame(tn)           # v
regEnd, regMat          = ifp.getRegFileMatrix(tn)          # v
iffMat                  = ifp.getIffFileMatrix(tn)          # v
regImgList              = ifp.registrationRedux(regMat)     # /v\
ifp.iffRedux(tn, iffMat, ampVect, modesList, template)      # v
ifp.saveCube(tn)                                            # v
# check cube
cube = rd.read_phasemap(os.path.join(fn.INTMAT_ROOT_FOLDER, tn, 'IMCube.fits'))

##### !!! Plots the filelist data
#----------------------------------------
# for file in filelist:
#     img = rd.readFits_maskedImage(file)
#     plt.figure()
#     plt.imshow(img)
#----------------------------------------
##### For visualization only

# !!!
def rename4D(tn):
    fold = os.path.join(imgFold, tn)
    files = os.listdir(fold)
    for file in files:
        if file.endswith('.4D'):
            num_str = file.split('.')[0]
            if num_str.isdigit():
                num = int(num_str)
                new_name = f"{num:04d}.4D"
                old_file = os.path.join(fold, file)
                new_file = os.path.join(fold, new_name)
                os.rename(old_file, new_file)

def create_sample_files(tn, number):
    directory = os.path.join(imgFold, tn)
    for i in range(0, number):
        file_name = f"{i}.4D"
        file_path = os.path.join(directory, file_name)
        with open(file_path, 'w') as f:
            f.write(f"Sample content for file {file_name}")
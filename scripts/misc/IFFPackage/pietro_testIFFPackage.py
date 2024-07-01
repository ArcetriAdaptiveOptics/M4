# Test IFF Package Pietro
import os
from m4.devices import deformable_mirror as dfm
from astropy.io import fits
from m4.dmutils import iff_acquisition_preparation as ifa
from scripts.misc.IFFPackage import iff_processing as ifp
from m4.configuration import config_folder_names as fn
from m4.ground import read_data as rd
from m4.configuration import read_iffconfig
from m4.utils.roi import ROI
from m4.utils import image_registration_lib as iml
roi = ROI()
m4 = dfm.M4AU()
ifa= ifa.IFFCapturePreparation(m4)
opd_fold = fn.OPD_IMAGES_ROOT_FOLDER
iff_fold = fn.IFFUNCTIONS_ROOT_FOLDER
tn = '20160516_114916'
tn2= '20160516_114917'
#______Data_Simulation_________________________________________________________
# Position correctly the segment
truss.moveTrussBy(0.500)
angrot.rotateBy(-30)
dm.setZerosToSegActs()
figure(); imshow(interf.acquire_phasemap()); colorbar()
# Create and run the timed command matrix history
ifa._updateModalBase('mirror')
t           = ifa._createTriggerPadding()
tcmdh       = ifa.createTimedCmdHistory() #;imshow(tcmdh);colorbar()
iff_par     = ifa.getInfoToSave()
regActs     = read_iffconfig.getConfig(('REGISTRATION'))['modes']
images = []
for column in tcmdh.T:
    dm.setActsCommand(column, rel=False)
    images.append(interf.acquire_phasemap())
masked_images = []
for i in range(len(images)):
    mask = roi.roiGenerator(images[i])[9]
    new_masked_img = np.ma.masked_array(images[i], mask=mask)
    cropped_img = iml.crop_frame(new_masked_img)
    masked_images.append(cropped_img)
#### Plot Check
#-------------------------
    plt.figure()
    plt.imshow(masked_images[i])
    plt.colorbar()
    plt.title(f"image_{i:04d}")
#-------------------------
####
# Save. This should happen when the timed_cmd_history is sent and the acquisi-
# tion starts, so that a tn is created
if os.path.exists(opd_fold) is False:
    os.mkdir(opd_fold)
if os.path.exists(iff_fold) is False:
    os.mkdir(iff_fold)
for i in range(len(masked_images)):
    name = os.path.join(opd_fold, f"img_{i:04d}.fits")
    rd.save_phasemap(name, masked_images[i])
fits.writeto(os.path.join(iff_fold, 'cmdMatrix.fits'), iff_par['cmdMatrix'])
fits.writeto(os.path.join(iff_fold, 'ampVector.fits'), iff_par['ampVector'])
fits.writeto(os.path.join(iff_fold, 'indexList.fits'), iff_par['indexList'])
fits.writeto(os.path.join(iff_fold, 'modesVector.fits'), iff_par['modesList'])
fits.writeto(os.path.join(iff_fold, 'Template.fits'), iff_par['template'])
fits.writeto(os.path.join(iff_fold, 'registrationActs.fits'), regActs)
with open(os.path.join(iff_fold, 'shuffle.dat'), 'w') as file:
    file.write(str(iff_par['shuffle']))
#______________________________________________________________________________
# Processing tn
ampVect,modesList,template,_,_,_ = ifp._getAcqPar(tn)
trigF                   = ifp.getTriggerFrame(tn)
regEnd, regMat          = ifp.getRegFileMatrix(tn)
iffMat                  = ifp.getIffFileMatrix(tn)
regImgList              = ifp.registrationRedux(regMat)
ifp.iffRedux(tn, iffMat, ampVect, modesList, template)
ifp.saveCube(tn)
# Processing tn2
ampVect,modesList,template,_,_,_ = ifp._getAcqPar(tn2)
trigF                   = ifp.getTriggerFrame(tn2)
regEnd, regMat          = ifp.getRegFileMatrix(tn2)
iffMat                  = ifp.getIffFileMatrix(tn2)
regImgList              = ifp.registrationRedux(regMat)
ifp.iffRedux(tn2, iffMat, ampVect, modesList, template)
ifp.saveCube(tn2)

tnlist = [tn,tn2]
ifp.stackCubes(tnlist)

#_Data_Check_______________________________________________________________
fold = os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn)
filelist = ifp._getFileList(tn)
for file in filelist:
    plt.figure()
    plt.imshow(rd.read_phasemap(file))
    plt.colorbar()

#!!!___________________________________________________________________________
def rename4D(tn):
    fold = os.path.join(opd_fold, tn)
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
    directory = os.path.join(opd_fold, tn)
    for i in range(0, number):
        file_name = f"{i}.4D"
        file_path = os.path.join(directory, file_name)
        with open(file_path, 'w') as f:
            f.write(f"Sample content for file {file_name}")
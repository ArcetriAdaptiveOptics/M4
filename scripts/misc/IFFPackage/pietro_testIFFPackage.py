# Test IFF Package Pietro
#run '/home/pietrof/git/M4/m4/analysis_imports.py'
import os
from m4.devices import deformable_mirror as dfm
from astropy.io import fits
from scripts.misc.IFFPackage import iff_acquisition_preparation as ifa
from scripts.misc.IFFPackage import iff_processing_pietro as ifp
from m4.configuration import config_folder_names as fn
from m4.ground import read_data as rd
from m4.configuration import read_iffconfig
from m4.utils.roi import ROI
from m4.utils import image_registration_lib as iml
roi = ROI()
m4 = dfm.M4AU()
ifa= ifa.IFFCapturePreparation(dfm.M4AU())

# ciclo per creare dati simulati (pyott -i)
# Position the segment
truss.moveTrussBy(0.500)
angrot.rotateBy(-30)
dm.setZerosToSegActs()
figure(); imshow(interf.acquire_phasemap()); colorbar()
# Create and run the timed command matrix history
ifa._updateModalBase('hadamard')
t           = ifa._createTriggerPadding()
tcmdh       = ifa.createTimedCmdHistory(); imshow(tcmdh);colorbar()
cmdM        = ifa._cmdMatrix
ampVect     = ifa._modesAmp
indexList   = ifa._indexingList
modesVect   = ifa._modesList
template    = ifa._template
regActs     = read_iffconfig.getConfig(('REGISTRATION'))['modes']
images = []
for column in tcmdh.T:
    dm.setActsCommand(column)
    images.append(interf.acquire_phasemap())
    #dm.setZerosToSegActs()
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
# Saving the data simulated data 2 and stack cubes
tn  = '20160516_114916'
tn2 = '20160516_114916-2'
# Save
fold = os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn2)
ifFold = os.path.join(fn.IFFUNCTIONS_ROOT_FOLDER, tn2)
if os.path.exists(fold) is False:
    os.mkdir(fold)
if os.path.exists(ifFold) is False:
    os.mkdir(ifFold)
for i in range(len(masked_images)):
    name = os.path.join(fold, f"img_{i:04d}.fits")
    rd.save_phasemap(name, masked_images[i])
fits.writeto(os.path.join(ifFold, 'cmdMatrix.fits'), cmdM)
fits.writeto(os.path.join(ifFold, 'ampVector.fits'), ampVect)
fits.writeto(os.path.join(ifFold, 'indexList.fits'), indexList)
fits.writeto(os.path.join(ifFold, 'modesVector.fits'), modesVect)
fits.writeto(os.path.join(ifFold, 'Template.fits'), template)
fits.writeto(os.path.join(ifFold, 'registrationActs.fits'), regActs)
with open(os.path.join(ifFold, 'shuffle.dat'), 'w') as file:
    file.write('0')
# Processing
ampVect,modesList,template,_,_,_ = ifp._getAcqPar(tn2)
trigF                   = ifp.getTriggerFrame(tn2)
regEnd, regMat          = ifp.getRegFileMatrix(tn2)
iffMat                  = ifp.getIffFileMatrix(tn2)
regImgList              = ifp.registrationRedux(regMat)
ifp.iffRedux(tn2, iffMat, ampVect, modesList, template)
# cube save

# cubes stack
tn = '20160516_114916'
tn2 = '20160516_114916 (Copy)'
ifp.saveCube(tn)
ifp.saveCube(tn2)
ifp.stackCubes([tn,tn2])

#_Old_Data_Check_______________________________________________________________
tn  = '20160516_114916'
fold = os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn)
filelist = ifp._getFileList(tn)
for file in filelist:
    plt.figure()
    plt.imshow(rd.read_phasemap(file))
    plt.colorbar()

#!!!___________________________________________________________________________
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
# Test IFF Package Pietro
run '/home/pietrof/git/M4/m4/initOTT.py'
# run '/home/pietrof/git/M4/m4/analysis_imports.py'
from m4.devices import deformable_mirror as dfm
from scripts.misc.IFFPackage import iff_acquisition_preparation as ifa
from scripts.misc.IFFPackage import iff_processing as ifp
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
tn = '20160516_114916'
filelist = sorted([os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, (tn+'/'+image)) for image in os.listdir(os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn))])
trigImage, regFrames = ifp.getTriggerAndRegistrationFrames(tn)

##### !!! Plots the data
#----------------------------------------
# for file in filelist:
#     img = rd.readFits_maskedImage(file)
#     plt.figure()
#     plt.imshow(img)
#----------------------------------------
##### For visualization only

ifp.iffRedux(tn)
#############!!!
_,modesList,_,template,_  = read_iffconfig.getConfig('IFFUNC')
nPushPull = len(template)
indexingList = rd.readFits_data(os.path.join(fn.IFFUNCTIONS_ROOT_FOLDER, tn)+'/indexList.fits')    # to be implemented
amplitude = rd.readFits_data(os.path.join(fn.IFFUNCTIONS_ROOT_FOLDER, tn)+'/ampVector.fits')      # to be implemented
filelist = sorted([os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, (tn+'/'+image)) for image in os.listdir(os.path.join(fn.OPD_IMAGES_ROOT_FOLDER, tn))])
filelist = filelist[10:] #regEnd:10
for i in range(len(modesList)):
    print('Mode', i) 
    for k in range(nPushPull):
        p = nPushPull * i + k
        n = indexingList[i] # was 'p'
        mis_amp = k * indexingList.shape[0] + n
        mis = mis_amp*template.shape[0]
        # Crea il pacchetto di immagini del modo 'i', contenente nPushPull images
        file_name = filelist[p]
        image_list = []
        for l in range(0, template.shape[0]-1):
            file_name = filelist[p+l]
            ima = rd.readFits_maskedImage(file_name)
            image_list.append(ima)
        image = np.zeros((ima.shape[0], ima.shape[1]))
        # Algorimo differenziale
        for x in range(1, len(image_list)):
            opd2add = image_list[x]*template[x] + image_list[x-1]*template[x-1]
            master_mask2add = np.ma.mask_or(image_list[x].mask, image_list[x-1].mask)
            if x==1:
                master_mask = master_mask2add
            else:
                master_mask = np.ma.mask_or(master_mask, master_mask2add)
            image += opd2add
        image = np.ma.masked_array(image, mask=master_mask)
        norm_image = image / (2*amplitude[n] * (template.shape[0]-1))
    fold = os.path.join(fn.IFFUNCTIONS_ROOT_FOLDER, tn)
    img_name = os.path.join(fold, 'mode_{:04d}.fits'.format(i))
    rd.save_phasemap(img_name, norm_image)
##################
ifp.createCube(tn)


##!!!

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
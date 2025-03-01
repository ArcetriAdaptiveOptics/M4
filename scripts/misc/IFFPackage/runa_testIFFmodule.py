from importlib import reload
#!!!!!
#!!!!
#warning
#the management of data saving is missing
# who is in charge of saving data?
#where is the Tracknum created?
#maybe the save shall be operated somewhere else

from m4.devices import deformable_mirror as dm
from m4.iffutils import iff_acquisition_preparation as ifa
m4u = dm.M4AU()

ifc=ifa.IFFCapturePreparation(m4u)

#general comments:
#registration and trigger modes are defined via config file
#we don't want to pass parameters for the trigger/registration modes


#auziliary patterns: based on:
ifc._updateModalBase('zonal')  #è un comando interno, per aggiornare le modalBase prima di creare comandi. non serve mandarlo da terminare

t = ifc._createTriggerPadding()
r = ifc._createRegistrationPattern()
imshow(t)

#create the cmdMatrix for acquisition, with mirror modes (default)
nmodes = 5
mlist = arange(0,nmodes)
mm = ifc._createCmdMatrix(mlist)

#but in reality, the upper level script to be used is  createTimedCmdHistory
amp = 1e-7
ifc.createTimedCmdHistory(mlist,amp)

#to switch type of modal base use this (default is None, or 'mirror')
ifc.modalBaseId = 'zonal'
nmodes = 5
mlist = arange(0,nmodes)
mm = ifc._createCmdMatrix(mlist) #the default is as in the conf file
#or:
mm = ifc._createCmdMatrix(mlist,'zonal')
amp = 1e-7
tm = ifc.createCmdMatrixHistory(mlist,amp )
amp = np.zeros(nmodes)+1e-7
amp[2:3]=2e-7
tm = ifc.createCmdMatrixHistory(mlist,amp )

tm = ifc.createTimedCmdHistory(mlist,1e-7)


#### test processing
from scripts.misc.IFFPackage import iff_processing_runa as ifp
from m4.configuration import config_folder_names as fn
from importlib import reload
tn = '20160516_114916'

infoT, infoR, infoIF = ifp._getAcqInfo(tn)
ampVector, modesVector, template,indexList, registrationActs, shuffle=ifp._getAcqPar(tn)
filelist=ifp.getFileList(tn)
regMat, modesMat = ifp.findTriggerStartFrame(tn)
actImgList = ifp.registrationRedux(regMat, [1,-1])
ifp.iffRedux(modesMat, tn, ampVector, modesVector, template, shuffle)



#used to debug the trigger
    img0 = rd.read_phasemap(fileList[0])
    go = 1
    i = 0
    #repeat from here
    i = i+1
    img1 = rd.read_phasemap(fileList[i])
    rr2check = zern.removeZernike(img1-img0,[1,2,3]).std()
    imshow(zern.removeZernike(img1-img0,[1,2,3]))
    print(i)
    print(rr2check)
    print(rr2check > triggAmp/3)


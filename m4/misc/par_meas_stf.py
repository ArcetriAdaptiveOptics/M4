conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
#conf = '/home/marco/work/git/M4/m4/configuration/myConfig.yaml'
from m4.configuration import start
import numpy as np
import matplotlib.pyplot as plt
ott, interf, dm = start.create_ott(conf)
from m4.configuration.ott_parameters import Interferometer
from m4.mini_OTT import timehistory as th
from m4.analyzers.noise_data_analyzer import Noise

print("Frequency acquired: " , Interferometer.BURST_FREQ)

from m4 import noise

comp_stf=False
comp_ast=False
comp_ave_stf = True

#tnlist = ['20230126_125637', '20230126_130004', '20230126_131049', '20230126_132012', '20230126_160816']
tnlist = ['20230127_104054','20230127_115010','20230127_142613','20230127_143310','20230127_120016','20230127_120707','20230127_130826','20230127_125418','20230127_130203','20230127_121545','20230127_122227','20230127_123255','20230127_124244','20230127_145150','20230127_145945','20230127_151719','20230127_131524','20230127_131708','20230127_132159','20230127_145027','20230127_161043','20230127_163227','20230129_100739', '20230131_181616', '20230201_183722']
#tnlist = ['20230127_163227','20230129_100739']
path_series = '/mnt/data/M4/Data/M4Data/OPTData/OPD_series/'
path_images = '/mnt/data/M4/Data/M4Data/OPTData/OPDImages/'
path_optdata = '/mnt/data/M4/Data/M4Data/OPTData/'


tn = tnlist[0]
if comp_stf:
    #structure function
    for tn in tnlist:
        dfpath =  '/home/m4/4d/M4/Produced/'+tn
        print(dfpath)
        tau_vector = np.arange(1,100,2)
        noise.convection_noise(dfpath, tau_vector)


if comp_ast:
    zv = []
    q = []
    #astigmatism analysis
    for tn in tnlist:
        fl = th.fileList(tn)
        img = th.averageFrames(0,499,fl)
        zz, mat =th.zernike.zernikeFit(img, [1,2,3,4,5,6])
        img = th.removeZernike(img)
        q.append(img)
        print(tn, zz)
        zv.append(zz)

    zv = np.array(zv)

img_list =[]
be_list  =[]
stf_list =[]
if comp_ave_stf:
    #zv = []
    #img_list = []
    #astigmatism analysis
    nn = Noise()
    plt.figure()
    for tn in tnlist:
        fl = th.fileList(tn)
        print(tn+' files: %i ' %len(fl))
        if len(fl) > 1000:
            img = (th.averageFrames(0,int(len(fl)/2)-1,fl)+th.averageFrames(int(len(fl)/2),len(fl)-1,fl))/2
        else:
            img = th.averageFrames(0,len(fl)-1,fl)
            
        zz, mat =th.zernike.zernikeFit(img, [1,2,3,4,5,6,7,8])
        img_noz = th.removeZernike(img)
        plt.imshow(img_noz)
        plt.title(tn+' ave no8Z')
        plt.colorbar()
        plt.savefig("aveimg_"+tn+'.png')
        plt.clf()
        
        #zz, mat =th.zernike.zernikeFit(img, [1,2,3,4,5,6])
        #img = th.removeZernike(img)
        be, stf = nn.comp_spatial_stf(img, pixsc=3e-3)
        img_list.append(img)
        be_list.append(be)
        stf_list.append(np.sqrt(stf))
        #plt.figure()
        plt.plot(be, np.sqrt(stf))
        plt.title("stf "+tn)
        plt.xlabel('m')
        plt.ylabel('m RMS')
        plt.show()
        plt.savefig("stf_aveimg_"+tn+'.png')
        plt.clf()

        #zv.append(zz)

    #zv = np.array(zv)
    
        
fig=plt.figure(figsize=(10,15))
for ii in np.arange(len(tnlist)):
    pp=plt.plot(be_list[ii], stf_list[ii])
plt.title( "STF on average")
plt.legend(tnlist, loc='upper left')
plt.xlabel('m')
plt.ylabel('m RMS')
plt.show()
plt.savefig("stf_aveimg1.png")
#

fig=plt.figure()
for ii in np.arange(4)+len(tnlist)-4:
    print(tnlist[ii])
    pp=plt.plot(be_list[ii], stf_list[ii])
plt.title( "STF on average")
plt.legend(tnlist[-4:], loc='upper left')
plt.xlabel('m')
plt.ylabel('m RMS')
plt.show()
plt.savefig("stf_aveimg2.png")




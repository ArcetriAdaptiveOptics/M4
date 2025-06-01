conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
#conf = '/home/marco/work/git/M4/m4/configuration/myConfig.yaml'
from m4.configuration import ott
import numpy as np
import matplotlib.pyplot as plt
ott, interf, dm = ott.create_ott(conf)
from m4.configuration.ott_parameters import Interferometer
from m4.mini_OTT import timehistory as th
from m4.analyzers.noise_data_analyzer import Noise

print("Frequency acquired: " , Interferometer.BURST_FREQ)

from m4 import noise

comp_stf=True
comp_ast=False
comp_ave_stf = True

#tnlist = ['20230126_125637', '20230126_130004', '20230126_131049', '20230126_132012', '20230126_160816']
tnlist = ['20230127_104054','20230127_115010','20230127_142613','20230127_143310','20230127_120016','20230127_120707','20230127_130826','20230127_125418','20230127_130203','20230127_121545','20230127_122227','20230127_123255','20230127_124244','20230127_145150','20230127_145945','20230127_151719','20230127_131524','20230127_131708','20230127_132159','20230127_145027','20230127_161043','20230127_163227','20230129_100739', '20230131_181616', '20230201_183722','20230203_182356', '20230204_160455']
tnlist = ['20230127_163227','20230129_100739', '20230131_181616', '20230201_183722', '20230202_175230','20230203_182356', '20230204_160455' ]
path_series = '/mnt/m4storage/Data/M4Data/OPTData/OPD_series/'
path_images = '/mnt/data/M4/Data/M4Data/OPTData/OPDImages/'
path_optdata = '/mnt/data/M4/Data/M4Data/OPTData/'


tn = tnlist[0]
if comp_stf:
    #structure function
    for tn in tnlist:
        dfpath =  path_series+tn
        print(dfpath)
        tau_vector = np.arange(1,200,2)
        noise.convection_noise(dfpath, tau_vector)
    
        plt.savefig("timestf_img_"+tn+'.png')


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
freq_list=[]
psd_list=[]
if comp_ave_stf:
    #zv = []
    #img_list = []
    #astigmatism analysis
    nn = Noise()
    plt.figure()
    for tn in tnlist:
        cum_stf = None
        cum_psd = None
        fl = th.fileList(tn)
        print(tn+' files: %i ' %len(fl))
#        if len(fl) > 1000:
#            img = (th.averageFrames(0,int(len(fl)/2)-1,fl)+th.averageFrames(int(len(fl)/2),len(fl)-1,fl))/2
#        else:
            #img = th.averageFrames(0,len(fl)-1,fl)
        for ff in np.arange(len(fl[:-1])):
            ph1 = th.read_phasemap(fl[ff])
            ph2 = th.read_phasemap(fl[ff+1])
            img=ph2-ph1
            #zz, mat =th.zernike.zernikeFit(img, [1,2,3])
            img_noz = th.removeZernike(img)
        
            #zz, mat =th.zernike.zernikeFit(img, [1,2,3,4,5,6])
            #img = th.removeZernike(img)
            be, stf = nn.comp_spatial_stf(img_noz, pixsc=3e-3)
            freq,psd = th.comp_psd(img_noz)
            if cum_stf is None:
                cum_stf = stf/(len(fl)-1)
                cum_psd = stf/(len(fl)-1)
            else:
                cum_stf += stf/(len(fl)-1)
                cum_psd = psd/(len(fl)-1)
        plt.figure()
        plt.plot(be, np.sqrt(cum_stf), drawstyle='steps-post')
        plt.title("stf "+tn)
        plt.xlabel('m')
        plt.ylabel('m RMS')
        plt.show()
        plt.savefig("avestf_dimg_"+tn+'.png')
    
        plt.figure()
        plt.loglog(freq, np.sqrt(cum_psd),drawstyle='steps-post')
        plt.title("stf "+tn)
        plt.xlabel('m')
        plt.ylabel('m RMS')
        plt.show()
        plt.pause(1)
        #plt.clf()
        be_list.append(be)
        stf_list.append(np.sqrt(cum_stf))
        freq_list.append(freq)
        psd_list.append(np.sqrt(cum_psd))
        plt.savefig("avepsd_dimg_"+tn+'.png')
        #zv.append(zz)

    #zv = np.array(zv)
    
        
    fig=plt.figure()
    for ii in np.arange(len(tnlist)):
        pp=plt.plot(be_list[ii], stf_list[ii])
    plt.title( "average STF")
    plt.legend(tnlist, loc='lower right')
    plt.xlabel('m')
    plt.ylabel('m RMS')
    plt.show()
    plt.savefig("avestf_dimg1.png")
    #

    fig=plt.figure()
    for ii in np.arange(len(tnlist)):
        pp=plt.loglog(freq_list[ii], psd_list[ii], drawstyle='steps-post')
    plt.title( "average PSD")
    plt.legend(tnlist, loc='lower left')
    plt.xlabel('m')
    plt.ylabel('m RMS')
    plt.show()
    plt.savefig("avepsd_dimg1.png")
    #




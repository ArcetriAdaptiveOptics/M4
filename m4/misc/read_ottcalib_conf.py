import numpy as np
import configparser
import json
from m4.mini_OTT import timehistory as th
#config=configparser.ConfigParser()
basepath = th.foldname.OPT_DATA_FOLDER
fold = 'OTTCalibConf/'
parname = 'PAR'
ncgh_tn_marker   = 'cgh_tn_marker'
ncgh_tn_img      = 'cgh_tn_img'
ntnpar           = 'tnpar'
nmark_cgh_list   = 'mark_cgh_list'
nf0              = 'f0'
nf1              = 'f1'

ottname         = 'OTT'
nott_tn_marker   = 'ott_tn_marker'
nott_tn_img      = 'ott_tn_img'
nmark_ott_list   = 'mark_ott_list'
npx_ott          = 'px_ott'
def get_nparray(key1,key2):
    out   = np.array(json.loads(key1[key2]))
    return out
    
def gimmetheconf(tn):
        config=configparser.ConfigParser()
        fname = basepath+'/'+fold+tn+'.ini'
        print(fname)
        config.read(fname)
        #PAR
        pp = config[parname]
        cgh_tn_marker = pp[ncgh_tn_marker]
        cgh_tn_img      = pp[ncgh_tn_img]
        tnpar           = pp[ntnpar]
        mark_cgh_list   = get_nparray(pp, nmark_cgh_list)
        f0              = float(pp[nf0])
        f1              = float(pp[nf1])

        oo              = config[ottname]
        ott_tn_marker   = json.loads(oo[nott_tn_marker])
        ott_tn_img      = oo[nott_tn_img]
        mark_ott_list   = get_nparray(oo, nmark_ott_list)
        px_ott          = float(oo[npx_ott])
        '''
        if config.has_option(ottname,'ott_tn_marker1'):
            ott_tn_marker = [ott_tn_marker]
            ott_tn_marker1 = oo['ott_tn_marker1']
            ott_tn_marker.append( ott_tn_marker1)

        if config.has_option(ottname,'ott_tn_marker2'):
            ott_tn_marker2 = oo['ott_tn_marker2']
            ott_tn_marker.append( ott_tn_marker2)
        '''
        return cgh_tn_marker,cgh_tn_img,tnpar,mark_cgh_list,f0,f1,ott_tn_marker,ott_tn_img,mark_ott_list,px_ott






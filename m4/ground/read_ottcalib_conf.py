import configparser
import json
from m4.configuration.root import folders as foldname
import os.path as _op

basepath = foldname.OPT_DATA_ROOT_FOLDER
calibfold = basepath + '/OTTCalibConf'  #foldname.OTTCALIB_ROOT_FOLDER
parname = "PAR"
ncgh_tn_marker = "cgh_tn_marker"
ncgh_tn_img = "cgh_tn_img"
ntnpar = "tnpar"
nmark_cgh_list = "mark_cgh_list"
nf0 = "f0"
nf1 = "f1"
ottname = "OTT"
nott_tn_marker = "ott_tn_marker"
nott_tn_img = "ott_tn_img"
nmark_ott_list = "mark_ott_list"
npx_ott = "px_ott"


def _get_nparray(key1, key2):
    data = list(json.loads(key1[key2]))
    pp = []
    for ii in data:
        if ii.__class__ == int:
            pp.append(ii)
        else:
            for jj in ii:
                pp.append(jj)
    # return np.array(pp)
    return data


def gimmetheconf(tn):
    """
    Parameters
    ---
    tn:string
    tracking number of the calibration configuration

    Returns
    ---
    cgh_tn_marker:string
    tracking number of the CGH markers
    cgh_tn_img:string
    tnpar,mark_cgh_list,f0,f1,ott_tn_marker,ott_tn_img,mark_ott_list,px_ott

    """
    config = configparser.ConfigParser()
    fname = _op.join(calibfold, tn + ".ini")
    print(fname)
    config.read(fname)
    # PAR
    pp = config[parname]
    cgh_tn_marker = pp[ncgh_tn_marker]
    cgh_tn_img = pp[ncgh_tn_img]
    tnpar = pp[ntnpar]
    mark_cgh_list = _get_nparray(pp, nmark_cgh_list)
    f0 = float(pp[nf0])
    f1 = float(pp[nf1])

    oo = config[ottname]
    ott_tn_marker = json.loads(oo[nott_tn_marker])
    ott_tn_img = oo[nott_tn_img]
    mark_ott_list = _get_nparray(oo, nmark_ott_list)
    px_ott = float(oo[npx_ott])

    # if config.has_option(ottname,'ott_tn_marker1'):
    #    ott_tn_marker = [ott_tn_marker]
    #    ott_tn_marker1 = oo['ott_tn_marker1']
    #    ott_tn_marker.append( ott_tn_marker1)

    # if config.has_option(ottname,'ott_tn_marker2'):
    #    ott_tn_marker2 = oo['ott_tn_marker2']
    #    ott_tn_marker.append( ott_tn_marker2)

    return (
        cgh_tn_marker,
        cgh_tn_img,
        tnpar,
        mark_cgh_list,
        f0,
        f1,
        ott_tn_marker,
        ott_tn_img,
        mark_ott_list,
        px_ott,
    )

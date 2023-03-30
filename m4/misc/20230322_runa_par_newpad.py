import numpy as np
from m4.mini_OTT import timehistory as th
from m4.misc import pad_intmat as padim
im = padim.invmat

tn = '20230303_102203'
tn = '20230201_183722'

tn = ['20230321_210648','20230321_225557','20230322_072934']
st = []
cc = []
im = []
for i in tn:
    fl=th.fileList(i)
    img = th.averageFrames(0,499, fl)
    img = th.removeZernike(th.averageFrames(0,499, fl),[1,2,3,4])
    im.append(img)
    st.append(img.std())
    cx, mat=th.zernike.zernikeFit(img,[1,2,3,4,5,6,7,8,9,10,11])
    cc.append(cx)

wf = np.mean(cc,0)*2


tn0 = '20230201_183722'
fl=th.fileList(tn0)
img = th.averageFrames(0,499, fl)
cx, mat=th.zernike.zernikeFit(img,[1,2,3,4,5,6,7,8,9,10,11])
tt = np.array([cx[4],cx[5],cx[8],cx[10]]) # array([-2.33241245e-08, -3.59469434e-09, -5.45536396e-09,  9.59816912e-09])
cmd = -np.matmul(im,tt) #array([-0.06490412, -0.07268578,  0.00436887, -0.05221968, -0.01599924,   -0.0822043 ])

tn1 = '20230322_072934'
tn1 = '20230321_210648'
fl=th.fileList(tn1)
img = th.averageFrames(0,499, fl)
cx, mat=th.zernike.zernikeFit(img,[1,2,3,4,5,6,7,8,9,10,11])
tt = np.array([cx[4],cx[5],cx[8],cx[10]]) #array([-3.25386752e-08,  6.04059910e-09, -4.05043035e-09,  1.20988564e-08])
cmd = -np.matmul(im,tt) #array([-0.07678485, -0.11027378,  0.01507068, -0.07554086,  0.0047513 ,  -0.11821489])

tn2 = '20230328_182436'
fl=th.fileList(tn2)
img = th.averageFrames(0,499, fl)
cx, mat=th.zernike.zernikeFit(img,[1,2,3,4,5,6,7,8,9,10,11])
tt = np.array([cx[4],cx[5],cx[8],cx[10]]) #array([-3.25386752e-08,  6.04059910e-09, -4.05043035e-09,  1.20988564e-08])
cmd = -np.matmul(im,tt) #array([-0.07678485, -0.11027378,  0.01507068, -0.07554086,  0.0047513 ,  -0.11821489])



tn0 = '20230201_183722'
fl=th.fileList(tn0)
img0 = th.averageFrames(0,199, fl)
img1 = th.averageFrames(200,399, fl)
dd = th.removeZernike((img1-img0),[1,2,3,4,5,6,7,8,9,10])


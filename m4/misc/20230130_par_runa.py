from m4.mini_OTT import timehistory as th

#Id	Tracknum	Box	Fan	Net	 	                Astigm	Ast.std	notes
# 1	20230127_104054	1	0	0	            	        -1; -16	11;14
# 2	20230127_115010	0	0	0	 	                -20;-8
# 3	20230127_142613	0	0	0	 	                -27; -3	5;6	after lunch break
# 4	20230127_143310	0	0	0	 	                -23; -7	7; 8
# 5	20230127_120016	0	1	0	large fan E-W,1m,30°	-5;-10
# 6	20230127_120707	0	3	0	large fan E-W,1m,30°	-18;-19
# 7	20230127_130826	0	1	0	large fan E-W,1m,30°	-15; -17	9; 10	second run
# 8	20230127_125418	0	2	0	large fan E-W,1m,30°	-16; -8	13; 10	second run
# 9	20230127_130203	0	3	0	large fan E-W,1m,30°	-17; -12	11; 11	second run
# 10	20230127_121545	0	1	0	large fan N-S,0.5m,30°	-24; -12
# 11	20230127_122227	0	3	0	large fan N-S,0.5m,30°	-18; -18
# 12	20230127_123255	0	1	0	small fan E-W, 0.5m,0°	-23; -2
# 13	20230127_124244	0	1	0	smallFan, N-S, 0.5m,0°	-40; -1	 	maybe we were closer??
# 14	20230127_145150	0	3	2	large fan E-W,1m,30°	-1; -2
# 15	20230127_145945	0	3	2	large fan E-W,1m,30°	1; -11
# 16	20230127_151719	0	3	2	large fan E-W,1m,30°	-3; 84	 	doors open

import numpy as np
from matplotlib import pyplot as plt

def companalysis(tna,tnb):
    fla = th.fileList(tna)
    flb = th.fileList(tnb)
    qa= th.averageFrames(0,99, fla)
    qb= th.averageFrames(0,99, flb)
    dq = qb-qa
    dq = th.removeZernike(dq)
    print('StDev' ,round(dq.std()*1e9,1))
    c,m = th.zernike.zernikeFit(dq,[1,2,3,4,5,6,7,8,9,10])
    print('Ast: ',np.round(c[4:6]*1e9,1))
    print('Coma: ', np.round(c[6:8]*1e9,1))
    print('Tref: ', np.round(c[8:10]*1e9,1))
    plt.clf();    plt.imshow(dq); plt.colorbar(); plt.title(tna+' - '+tnb)

tn1 = '20230127_104054'
tn2 = '20230127_115010'
tn3 = '20230127_142613'
tn4 = '20230127_143310'
tn6 = '20230127_120707' #0       3       0       large fan E-W,1m,30°
tn7 = '20230127_130826' #0       1       0       large fan E-W,1m,30°
tn8 = '20230127_125418' #0       2       0       large fan E-W,1m,30°
tn9 = '20230127_130203' #0       3       0       large fan E-W,1m,30°
tn10= '20230127_121545' #0       1       0       large fan N-S,0.5m,30°
tn11= '20230127_122227' #0       3       0       large fan N-S,0.5m,30°
tn12= '20230127_123255' #0       1       0       small fan E-W, 0.5m,0°
tn13= '20230127_124244' #0       1       0       smallFan, N-S, 0.5m,0°
tn14= '20230127_145150' #0       3       2       large fan E-W,1m,30°
tn15= '20230127_145945' #0       3       2       large fan E-W,1m,30°

companalysis(tn3,tn4)

tn16= '20230127_151719'
fl2 = th.fileList(tn2)
fl2 = th.fileList(tn2)
fl3 = th.fileList(tn3)
fl4 = th.fileList(tn4)
fl9 = th.fileList(tn9)
fl11 = th.fileList(tn11)
fl15 = th.fileList(tn15)
fl16= th.fileList(tn16)


q2= th.averageFrames(0,99, fl2)
q3= th.averageFrames(0,99, fl3)
q4= th.averageFrames(0,99, fl4)

q9= th.averageFrames(0,99, fl9)
q11= th.averageFrames(0,99, fl11)
q15= th.averageFrames(0,99, fl15)

q16= th.averageFrames(0,99, fl16)

dq11 = th.removeZernike(q11-q2)
dq9 = th.removeZernike(q9-q2)
dq15 = th.removeZernike(q15-q2)
dq16 = th.removeZernike(q16-q2)

plt.imshow(dq16); plt.colorbar()



tn7 = '20230127_130826'
tn8 = '20230127_125418'
tn9 = '20230127_130203'
tn14 = '20230127_145150'
tn15 = '20230127_145945'

fl7 = th.fileList(tn7)
cc7 = th.cubeFromList(fl7)
z7  = th.zernikePlot(fl7,[1,2,3,4,5,6])
st7 = np.std(cc7,0)
q7  = th.removeZernike(th.averageFrames(0,99, fl7),[1,2,3])
plt.imshow(st7)

fl8 = th.fileList(tn8)
cc8 = th.cubeFromList(fl8)
z8  = th.zernikePlot(fl8,[1,2,3,4,5,6])
st8 = np.std(cc8,0)
q8  = th.removeZernike(th.averageFrames(0,99, fl8),[1,2,3])
plt.imshow(st8)

fl9 = th.fileList(tn9)
cc9 = th.cubeFromList(fl9)
z9  = th.zernikePlot(fl9,[1,2,3,4,5,6])
st9 = np.std(cc9,0)
q9  = th.removeZernike(th.averageFrames(0,99, fl9),[1,2,3])
plt.imshow(st9)

fl14 = th.fileList(tn14)
cc14 = th.cubeFromList(fl14)
z14  = th.zernikePlot(fl14,[1,2,3,4,5,6])
st14 = np.std(cc14,0)
q14  = th.removeZernike(th.averageFrames(0,99, fl14),[1,2,3])
plt.imshow(st14)

fl15 = th.fileList(tn15)
cc15 = th.cubeFromList(fl15)
z15  = th.zernikePlot(fl15,[1,2,3,4,5,6])
st15 = np.std(cc15,0)
q15  = th.removeZernike(th.averageFrames(0,99, fl15),[1,2,3])
plt.imshow(st15)


print(np.mean(z7[3,:]), np.std(z7[3,:]))
print(np.mean(z8[3,:]), np.std(z8[3,:]))
print(np.mean(z9[3,:]), np.std(z9[3,:])) 

m = '?'
a5=th.zernike.zernikeSurface(q9,1,m[:,4])
a6=th.zernike.zernikeSurface(q9,1,m[:,5])

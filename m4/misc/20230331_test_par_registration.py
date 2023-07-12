import numpy as np
from m4.mini_OTT import timehistory as th
from m4.utils.parabola_identification import ParabolaActivities
from matplotlib.patches import Circle
from matplotlib import pyplot as plt

'''
OTT
TN
20230707_142156


CGH
2000x2000, offs 28,0
TN          	DIM MARKERS [MM]	N_MEAS
20230328_175944	20	                20
20230328_180252	20	                20
20230328_181123	10	                20
20230328_181333	10	                20

1680x1680, offs 180,180
TN	            DIM MARKERS [MM]	N_MEAS
20230331_181519	20	                20
20230331_181658	20	                20
20230331_181946	10	                20
20230331_182118	10	                20

relative offset is 180-28, 180
'''
pa = ParabolaActivities()
offs = [180-28,180]
nmark = 25
offvec = np.zeros([2,nmark])
for i in range(nmark):
    offvec[:,i]=offs

circ0 = [8,11,12]
circ1 = [6,7,13,14,17,18]
circ2 = [3,4,5,15,16,20,21]
circ3 = [0,1,2,9,10,19,22,23,24]
radii = np.array([0.1, 0.2505, 0.452, 0.6911])

tnlist = ['20230328_175944','20230328_180252','20230328_181123','20230328_181333']
tnlist0= ['20230331_181519','20230331_181658','20230331_181946','20230331_182118']
fl = th.fileList(tnlist[0])
img=th.frame(0, fl)
vv=pa.rawMarkersPos(img)
pos = pa.filterMarkersPos(vv, 400,500)
#a=vv['area']

def coord2ima(cc):
    cc1=np.array([cc[1,:],cc[0,:]])
    #cc2=np.array([2000-cc[1,:],cc[0,:]])
    #cc3=np.array([cc[1,:],2000-cc[0,:]])
    #cc4=np.array([2000-cc[1,:],2000-cc[0,:]])

    return cc1

def getMarkers(tn, diam=24, thr=0.2):
    npix = 3.14*(diam/2)**2
    fl = th.fileList(tn)
    nf = len(fl)
    pos = np.zeros([2,25,nf])
    for j in range(nf):
        img = th.frame(j,fl)
        imaf = pa.rawMarkersPos(img)
        c0 = pa.filterMarkersPos(imaf, (1-thr)*npix, (1+thr)*npix)
        pos[:,:,j]=c0
    pos = np.average(pos,2)
    return pos
#--- First test
tn0 = '20230328_175944' #2000x2000
tn1 = '20230331_181519' # 1680x1680
#--- test before-after new pads
tn0 = '20230328_175944' #2000x2000
tn1 = '20230418_171241'
fl = th.fileList(tn0)
img=th.frame(0, fl)
p0 = getMarkers(tn0)
p1 = getMarkers(tn1)
p1 = p1+offvec
p0 = coord2ima(p0)
p1 = coord2ima(p1)
dd = p1-p0
rr = np.sqrt(dd[0,:]**2+dd[1,:]**2)

imshow(img)
scatter(p0[0,:],p0[1,:],c=rr, s=100); plt.axis('equal')
for i in range(nmark):
    text(p0[0,i],p0[1,i],i)


figure(2)

scatter(p0[0,:],-p0[1,:],c=rr, s=100); plt.axis('equal')
for i in range(nmark):
    text(p0[0,i]+10,-p0[1,i],i)
for i in range(nmark):
    rrs = ';'+str(int(rr[i]))+'px'
    text(p0[0,i]+70,-p0[1,i],rrs)

#check distortion
c0, axs0, r0 = fitEllipse(p0[0,circ0],p0[1,circ0]); c0 = np.real(c0)
c1, axs1, r1 = fitEllipse(p0[0,circ1],p0[1,circ1]); c1 = np.real(c1)
c2, axs2, r2 = fitEllipse(p0[0,circ2],p0[1,circ2]); c2 = np.real(c2)
c3, axs3, r3 = fitEllipse(p0[0,circ3],p0[1,circ3]); c3 = np.real(c3)
pradii = np.array([r0, r1, r2, r3])
ps = pradii/radii
plt.plot(radii,ps, 'x')
plt.plot(radii,pradii, 'x');xlabel('Markers mech radius [m]'); ylabel('Optical radius [pix]')
z=np.polyfit(radii,pradii,2)#x, y
xx = np.linspace(radii[0],radii[-1],100)
yy = z[0]*xx**2+z[1]*xx+z[2]
plt.plot(xx, yy)
radiin = radii/radii[-1]
pradiin = pradii/pradii[-1]
plt.plot(radiin,pradiin, 'x');xlabel('Markers mech radius [m]'); ylabel('Optical radius [pix]')
z=np.polyfit(radiin,pradiin,2)#x, y
xx = np.linspace(radiin[0],radiin[-1],100)
yy = z[0]*xx**2+z[1]*xx+z[2]
plt.plot(xx, yy)
z1 = np.array([0.115, 0.879, 0]) # dati da Zemax
yy1 = z1[0]*xx**2+z1[1]*xx+z1[2]
plt.plot(xx, yy1)




#---------------
pos = np.zeros([2,25,20,4])
for i in range(len(tnlist)):
    fl = th.fileList(tnlist[i])
    for j in range(len(fl)):
        img = th.frame(j,fl)
        imaf = pa.rawMarkersPos(img)
        c0 = pa.filterMarkersPos(imaf, 100,500)
        pos[:,:,j,i]=c0
       
v=[]
for i in range(25):
    v.append([np.std(pos[0,i,:,:]),np.std(pos[1,i,:,:])])


pos1 = np.zeros([2,25,20,4])
for i in range(len(tnlist0)):
    fl = th.fileList(tnlist0[i])
    for j in range(len(fl)):
        img = th.frame(j,fl)
        imaf = pa.rawMarkersPos(img)
        c0 = pa.filterMarkersPos(imaf, 100,500)
        if shape(c0)[1] == 25:
            pos1[:,:,j,i]=c0


i=0; j=0
plt.plot(pos[0,:,i,0],pos[1,:,i,0],'x'); plt.axis('equal')
plt.plot(pos1[0,:,i,0]+180-28,pos1[1,:,i,0]+180,'ro')

dpos = np.zeros([2,25,20,2])
for i in range(25):
    for j in range(20):
        for k in range(2):
            dpos[:,i,j,k] = pos[:,i,j,k]-pos1[:,i,j,k]-[180-28,180]

plt.plot(dpos.flatten())
dd = dpos[:,:,0,0]
rr = np.sqrt(dd[0,:]**2+dd[1,:]**2)
scatter(pos[0,:,0,0],pos[1,:,0,0],  c=rr, marker='o', s=200)

scatter(pos[0,:,0,0],pos[1,:,0,0],  c=dpos[0,:,0,0], marker='o', s=200)
scatter(pos[0,:,0,0],pos[1,:,0,0],  c=dpos[1,:,0,0], marker='o', s=200)

fl = th.fileList(tnlist[0])
img=th.frame(0, fl)
imshow(img)
plt.plot(pos[0,:,0,0],pos[1,:,0,0],'x')

fig,ax = plt.subplots(1)
ax.set_aspect('equal')
ax.imshow(img)
for xx, yy in zip(imaf[1,:],imaf[0,:]):
    circ = Circle((xx,yy),20, color='red')
    ax.add_patch(circ)


def fitEllipse(x,y):
    x = x[:, np.newaxis]
    y = y[:, np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T, D)
    C = np.zeros([6, 6])
    C[0, 2] = C[2, 0] = 2
    C[1, 1] = -1
    E, V =  eig(np.dot(inv(S), C))
    #         import pdb
    #         pdb.set_trace()
    n = np.argmax(np.abs(E))
    a_vect = V[:, n]
    #center
    b, c, d, f, g, a = a_vect[1]/2, a_vect[2], a_vect[3]/2, a_vect[4]/2, a_vect[5], a_vect[0]
    num = b*b - a*c
    x0 = (c*d - b*f)/num
    y0 = (a*f - b*d)/num
    centro = np.array([x0, y0])
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1 = (b*b-a*c)*((c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2 = (b*b-a*c)*((a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1 = np.sqrt(np.abs(up/down1)) #prima non c'era il valore assoluto
    res2 = np.sqrt(np.abs(up/down2)) #attenzione
    axs = np.array([res1, res2])
    raggio = axs.mean()
    return centro, axs, raggio
 

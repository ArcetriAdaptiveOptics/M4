from m4.mini_OTT import timehistory as th
from m4.ground import zernike as zern

tnlist= ['20231005_214542','20231006_123231','20231006_213044','20231006_234107','20231007_114524','20231019_124135','20231019_155257','20231019_224726']

trusstempid = 18
zlist=[1,2,3,4,5,6]
zz=[]
t=[]
for i in tnlist:
    q=th.openAverage(i)
    z,m =zern.zernikeFit(q,zlist)
    zz.append(z)
    temp=th.readTemperatures(i)
    t.append(np.mean(temp[:,trusstempid]))

zz =np.array(zz)
t=np.array(t)
figure()
plot(zz[:,4],'-o');title('Astigmatism variability across datasets');grid('on')
plot(zz[:,5],'-d')
legend(['Z5','Z6'])
figure()
plot(t,'-o'); title('Truss Temperatures °C')



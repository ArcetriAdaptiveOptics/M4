par=ott.parabola.getPosition()
rm=ott.referenceMirror.getPosition()

tnc = '20231029_200447' #rm600
tnc = '20231029_204423' #rm200
off = 30
idx = 4
par[idx]-=off
rm[idx]+=off*2.05
ott.parabola.setPosition(par)
ott.referenceMirror.setPosition(rm)
tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf, dofidf, n_frames=1, doit=True)



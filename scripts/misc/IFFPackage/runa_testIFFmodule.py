from importlib import reload


from scripts.misc.IFFPackage import iff_acquisition_preparation as ifa
ifc=ifa.IFFCapturePreparation()
#voglio usare i modi dello specchio:
#non faccio niente
#voglio usare i modi zonali o altri:
ifc.updateCommandMatrix('zonal')
ifc.updateCommandMatrix() #ritorna a modale
t=ifc._createTriggerPadding()

mlist = arange(1,5)
mm = ifc._createCmdMatrix(mlist)
tm = ifc.createCmdMatrixHistory(arange(1,5), 1e-7)

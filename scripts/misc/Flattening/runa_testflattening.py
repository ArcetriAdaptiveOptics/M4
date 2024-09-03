#in generale, guarda /home/runa/git/M4/m4/flattening.py per confronto
#non c'è import
tnIntmat = 'xxx'
ff = Flattening(tnIntmat) #servono solo i dati. i dati relativi al DM sono importato con i dati della misura, anche perchè la cmdmatrix deve essere quella salvata, non quella attuale.
nmodes = 100
fmodes = arange(0,nmodes)
nframes = 20
tnflat = ff.flatten(fmodes, nframes)


#come importare il Rec
from m4.analyzers import compute_reconstructor as cmprec
rec = ComputeReconstructor(tn)
recmat =rec.run(nmodes)
#chiedere a MX di modificare l'init di cmprec in modo che serva solo il TN
#oppure che capisca se gli hai passato il tn o i dati (utile se ripeti diverse volte??)

#chiedere a MX di implementare la modalità "seleziona quanti modi tagliare" e fa vedere lo spettro

#chiedere a MX di implementare il filtraggio di N Zernike su analysis mask e predisporre per ROI e operazioni per rimozione tilt globale e pistone

#il filtraggio degli zernike può avvenire all'init? forse no, meglio farlo come operazione separata e prevedere due variabili "interne": _intMat e _intMatFiltered, e la Rec si fa sempre sulla filtered



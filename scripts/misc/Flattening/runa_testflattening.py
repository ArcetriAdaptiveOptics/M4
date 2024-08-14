#non c'è import
tnIntmat = 'xxx'
ff = Flattening(tnIntmat) #servono solo i dati. i dati relativi al DM sono importato con i dati della misura, anche perchè la cmdmatrix deve essere quella salvata, non quella attuale.
nmodes = 100
fmodes = arange(0,nmodes)
nframes = 20
tnflat = ff.flatten(fmodes, nframes)

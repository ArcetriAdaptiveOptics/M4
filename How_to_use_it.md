## Calibrazione ed allineamento OTT ##
 1) __Calibrazione della parabola con lo specchio di riferimento__
    - Scegliere un vettore di comandi con cui calibrare i gradi di libertà della parabola (PAR) e quello dello specchio di riferimento (RM).
 		 ```
 		  command_amp_vector = np.array([par_piston, par_tip, par_tilt, rm_tip, rm_tilt])
 		 ```

	espressi in [mm, arcsec, arcsec, arcsec, arcsec]  
 	NOTA: tenere presente che il metodo di calibrazione usato prevede che all'applicazione del comando par_tip corrisponda un'applicazione del comando rm_tip=-2.05*par_tip. Stessa relazione sussiste tra par_tilt ed rm_tilt.
    - Possibiltà di visualizzazione della matrice dei comandi di calibrazione prima dell'effetiva applicazione tramite il comando
 		  ```
 		  main.showCommandMatrixBeforeCalibration(command_amp_vector)
 		  ```
 
    - Il comando per effettuare la calibrazione è
      ```
      main.calibrate_PARAndRM(ott, interf, n_frames, command_amp_vector, nPushPull)
      ```
      dove:  
      ott e interf sono gli oggetti precedentemente creati nello start up  
      n_frame corrisponde al numero di frame da acquisire con l'interferometro per ogni misura da effettuare  
      command_amp_vector è il vettore da utilizzare per la calibrazione  
      nPushPull corrisponde al numero di volte che verrà effettuata la misura per ognuno dei 5 comandi applicati
    - La funzione di calibrazione ritorna il tracking number (tt_cal) della misura effettuata. Questo contiene:  
      CalibrationInfo.fits dove sono scritti insieme tutti i dati utilizzati per la calibrazione  
      Gli stessi dati separati in più file (CMat.fits, CommandAmplitude.fits, InteractionMatrix.fits, Mask.fits)  
      Tutti gli interferogrammi acquisiti

 2) __Allineamento della parabola e specchio di riferimento__
 	- Prima di allineare è possibile calcolare i comandi di allineamento di parabola e specchio di riferimento tramite il comando
 		```
 		main.showCommandForParAndRmBeforeAlignement(ott, interf, tt_cal, n_images, zernike_to_be_corrected, dof_command_id)
 		```
 		dove:  
 		ott e interf sono gli oggetti precedentemente creati nello start up  
 		tt_cal è la stringa contenente il tracking number di calibrazione che si vuole utilizzare  
 		n_images corrisponde al numero di frame da acquisire con l'interferometro per l'immagine da allineare  
 		zernike_to_be_corrected è il vettore che indica quali modi di zernike correggere (np.array([0,1,2,3,4]) indicano rispettivamente tip, tilt, fuoco, coma, coma)  
 		dof_command_id è il vettore che indica con quali gradi di libertà si vuole fare la correzione  
 		NOTA: gli zernike possono essere [0,1], [0,1,3,4], [0,1,2,3,4]  
 			e vanno in coppia con i dof [3,4], [1,2,3,4], [0,1,2,3,4]
    - Il comando per effettuare l'allineamento è
    ```
    main.align_PARAndRM(ott, interf, tt_calib, n_images, zernike_to_be_corrected=None, dof_command_id=None)
    ```
    - La funzione stampa i comandi applicati a PAR e RM, gli zernike calcolati e restituisce il traching number della misura di allineamento appena eseguita
	   	che contiene:  
	   	AlignmentInfo.fits in cui sono scritte tutte le info usate per allineare  
	   	le stesse informazioni separte (commandId.fits, IntMatModeVector.fits che sono gli zernike_to_be_corrected, PositionAndDeltaCommand.fits 
	   	e Zernike.fits che contiene gli zernike selezionati prima della correzione)
	   	Immagini prima e dopo l'allineamento
	   	
3) __Calibrazione M4__

4) __Allineamento M4__
 
 
 ## Attività su mOTT ##
 
 ### Allineamento asse ottico e asse di rotazione meccanico ###
__Acquisizione dati__
- La funzione ha come obiettivo quello di ottenere interferegrammi a diversi valori di angolo di rotazione della torre
e, da questi, calcolare i valori di tip e tilt presenti. I vettori ottenuti (tip = [tipa at angle0, tip at angle1...],
tilt = [tilt at angle0, tilt at angle1 ...]) vengono utilizzati per fittare l'ellisse e restituire le coordinate del centro.
In questo modo è possibile conoscere come muovere le viti dello spider per riallineare la L3.
- I comandi per eseguire la misura e l'analisi sono:
	```
	from m4.utils.rotation_and_optical_axis_alignment import RotOptAlign
    ro = RotOptAlign(ott, interf)
    tt = ro.image_acquisition(start_point, end_point, n_points)
    centro, axs, raggio = ro.data_analyzer()
	```
- Questi restituiscono l'oggetto _ro_, con cui è possibile rianalizzare o rivisualizzare i dati, il tracking number della misura effettuata
  che contiene:  
	Cube.fits che è il cubo delle immagini acquisite  
	angle_list.txt che è lista degli angoli usati per la misura e
	theta.fits, cioè la stessa lista di angoli in versione vettore  
  centro, axs e raggio corrispondono agli output del fitting dell'ellisse sui vettori di tip e tilt ottenuti con le misure effettuate
- Per analizzare nuovamente una misura di cui si conosce il tracking number usare:
	```
	from m4.utils.rotation_and_optical_axis_alignment import RotOptAlign
	ro = RotOptAlign.reloadROObject(tt)
	centro, axs, raggio = ro.data_analyzer()
	```
NOTA: per maggiori informazioni fare riferimento anche alla seguente [pagina wiki](https://redmine.ict.inaf.it/projects/adopt_oaa/wiki/MOTT-20200915)

### Misure ###
Per ottenere la classe che permette di fare le acquisizioni dati sulla mini ott è necessario usare i comandi
```
from m4.mini_ott.measurements import Measurements
meas = Measurements(ott, interf)
```
- meas.opticalMonitoring(n_images, delay): in ingresso alla funzione si stabilisce il numero di misure da mediare durante l' acquisire e il ritardo
in secondi tra una acquisizione ed un altra. La funzione salva le immagini nella cartella OPD_SERIES e, per ogni acquisizione, viengono salvati i
primi 10 coefficienti di zernike e le temperature.
NOTA: la prima colonna della matrice con gli zernike contiene il dt dato dal tempo a cui avviene la misura ed il tempo zero di lancio del comando
di acquisizione

- meas.diffOpticalMonitoring(n_images, delayshort, delaylong): con questa funzione vengono acquisiste coppie di immagini a distanza di delayshort
tra loro due e delaylong dalla coppia di immagini successiva. La funzione salva le immagini nella cartella OPD_SERIES e, per ogni acquisizione, viengono salvati i
primi 10 coefficienti di zernike e le temperature.
NOTA: la prima colonna della matrice con gli zernike contiene il dt dato dal tempo a cui avviene la misura ed il tempo zero di lancio del comando
di acquisizione

- meas.actsRepeatability(nMeas, pistonValue, n_frames):

- meas.scanAstigmComa(stepamp, nstep, nframes=10):

- meas.parPistonTest(pistonValue, deltaposFilepath, amp, ttForAlign):

- meas.parTiltTest(act, val_vec):

- meas.mappingPar(shift, nIter, ttForAlign):

- meas.alignTest(tt, nImages, perturbationVec, pre=False):


### Analisi ###
Per ottenere la classe che permette di fare l'analisi dati della mini ott è necessario usare i comandi
```
from m4.mini_ott.analysis import Analysis
an = Analysis(tt)
```
To be continued

### Spider Test ###
```
from m4.mini_OTT.spider_test import SpiderTest()
sp = SpiderTest()
```

### Analisi dei requisiti ###
Con i comandi
```
from m4 import requirements_checker as rc
rc.analysis_req(data_file_path, zernike_vector_to_subtract, step=None, offset=None)
```
dove zernike_vector_to_subtract è il vettore degli zernike che si vuole sottrarre all'immagine robusta, step è la distanza tra le patches (necessaria 
per il calcolo del raggio di curvatura e dell'rms interactuator) e offset distingue il metodo di creazione dell'immagine robusta: se è None l'immagine viene
creata sottraendo tra loro i due cubi creati con la metà del numero di misure indicato, altrimenti al numero di misure indicato viene sottratta un'immagine di offset
precedentemente salvata nella cartella con le misure,  
NOTA: l'immagine di offset viene creata con i seguenti comandi
```
from m4.analyzers import requirement_analyzer as ra
ra.imageOpticOffset(data_file_path, start, stop)
```
Offset=None permette di utilizzate la procedura standard che prevede la creazione di tre immagini robuste (create utilizzando 50, 100 e 300 file presenti nella cartella
delle misure da analizzare). Le immagini robuste vengono analizzate e vengono automaticamente plottati e salvati i seguenti risultati:  
	- slop  
	- differential piston  
	- radius of curvature  
	- rms at the interactuator scale 31 mm  
	- rms at the interactuator scale 500 mm  

Nel caso in cui si voglia applicare l'analisi dei requisiti ad una serie di immagini utilizzare la funzione
```
slop_list, diff_piston_list, roc_list, rms31, rms500 = rc.fromImagesToReq(image_list, pscale=None, step=None, n_patches=None)
```
dove n_patches indica il numero di patches per il secondo taglio (se non indicato viene eseguito un solo taglio al centro dell'immagine)
Per visualizzare e salvare i valori ottenuti dell'analisi si può sfruttare la funzione rc.plotAndSaveForReqAnalysis(...): usare l'help
della funzione per maggiori dettagli.

...provare se la procedura sopra descritta funziona anche con una sola immagine (lista di un elemento). Altrimenti usare i comandi qui sotto...
Nel caso in cui si voglia applicare l'analisi dei requisiti avando a disposizione una sola immagine utilizzare
```
from m4.analyzers import requirement_analyzer as ra
slope = ra.test242(image, pscale)
diff_piston = ra.diffPiston(image)
roc = ra.test283(image, pscale, step)
rms31 = ra.test243(image, 0.015, pscale, step, n_patches)
rms500 = ra.est243(image, 0.1, pscale, step, n_patches)
```
NOTA: per maggiori informazioni fare riferimento anche alla seguente [pagina wiki](https://redmine.ict.inaf.it/projects/adopt_oaa/wiki/MOTT-20210408)

### Analisi del rumore ###
Le funzione per l'analisi del rumore si ottengono con il comando _from m4 import noise_ e comprendono:  
- noise.noise_vibrations(dataFilePath, numbersArray, tidyOrShuffle)  
	Questa funzione analizza i file presenti nella cartella dataFilePath utilizzando l'algoritmo scritto per l'analisi delle funzioni di influenza
	(con numero di push/pull pari a 6 e vettore di attuatori lungo nTotImagesInFolder / (template.size * nPushPull))
	e utilizzando l'input tidyOrShuffle per determinare l'ordine di accoppiamento delle immagini.  
	numbersArray è il vettore che indica tutte le lunghezze del vettore template che si vogliono utilizzare. 
	Esempio: numbersArray = np.array([3,5,7]) corrisponde all'utilizzo dei vettore template uguali a [1,-1,1], [1,-1,1,-1,1], [1,-1,1,-1,1,-1,1]  
	La funzione crea quindi una cartella con un cubo per ogni vettore di template utilizzato e per ognuno di essi calcola
	rms, somma in quadratura di tip e tilt, solo tilt e ptv delle immagini (Nota: a queste immagini viengono sottratti gli zernike
	1, 2 e 3). Restituisce poi un unico valore per ognuno di questi output mediandoli.  
	Salva e visualizza i plot in funzione del numero di template.
- noise.spectrumFromData(dataFilePath)
	La funzione calcola tip e tilt per ogni immagine presente nella cartella dataFilePath. Dei due vettori viene quindi calcolata
	la trasformata di fourier reale e vengono visualizzati e salvati i plot.
- noise.convection_noise(dataFilePath, tauVector)
	Questa funzione analizza i dati tramite funzione di struttura. Il vettore tau_vector determina la distanza tra le 
	immagini da utilizzare in modo accoppiato e 2*tau_vector l'intervallo tra le misure. Per ogni elemento del vettore tau
	vengono calcolati rms medio e somma quadratica media di tip e tilt: questi vengono salvati e visualizzati in un plot assieme
	al fit calcolato sul vettore di rms.
	Qualche esempio [qui](https://redmine.ict.inaf.it/projects/adopt_oaa/wiki/MOTT-20210401)
- noise.piston_noise(dataFilePath)
	Ad ogni immagine presente nella cartella viene sottratto tip e til e ne viene calcolata la media: del vettore ottenuto viene fatto
	e visualizzato lo spettro.

# Dispositivi #
Di seguito alcune indicazione su dove trovare/come usare alcuni dispositivi collegati al progetto.

## Accelerometri ##
 __Acquisizione dati__
```
ott.accelerometrs.acquireData(recording_seconds)
```
 Nota: se non specificato acquisisce per 5 secondi
 - L'acquisizione restituisce il tracking number della misura appena effettuata dove sono salvati:
 	i dati, dt, plc_id, directions, time, plc_range, plc_totcounts, sensitivity
 	
 
 __Analisi dati__
 ```
 from m4.analyzers.accelerometers_data_analyzer import AccelerometersDataAnalyzer
 an = AccelerometersDataAnalyzer(tt)
 ```
 - Dopodichè si hanno a disposizione diversi modi di accedere/vedere i dati
  	- spe, freq = an.getSpecAndFreq() 
  	- data = an.datah5
  	- an.plot_power_spectrum()
  	- an.readAndShow()  
Nota: l'analisi non salva nulla

## SPL ##
Il repositorio contenente il codice SPL si trova in ArcetriAdaptiveOptics al seguente [link](https://github.com/ArcetriAdaptiveOptics/SPL).

## Interferometro ##
### 4D PhaseCam 6110 ###
In case you want to use the full functionality of the device use:
```
from m4.devices.i4d import I4D
from m4.configuration.ott_parameters import Interferometer
interf = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)
```

the "interf" class contains all the commands documented by 4D. In particular:

- interf.burstFramesToSpecificDirectory(directory, numberOfFrames)
- interf.convertRawFramesInDirectoryToMeasurementsInDestinationDirectory(measurementsDirectory, rawFramesDirectory). 

the "burst" instruction requires the full absolute destination path, including the folder to be created  
the "convert" instruction requires the full absolute destination path, including the folder to be created and where to store the phasemaps

NOTA: in Arcetri i4d_IP = '193.206.155.193' e i4d_port = 8011

For the simple case of capturing and saving an image use:
```
from m4.devices.interferometer import I4d6110()
i4d6110 = I4d6110()
```

the "i4d6110" class contains the command for acquisition (i4d6110.acquire_phasemap()) and saving fits file (i4d6110.save_phasemap(location, file_name, masked_image)).

###### Read saved files
- .4D files:
```
from m4.ground.read_data import InterferometerConverter
ic = InterferometerConverter()
image = ic.fromPhaseCam6110(i4dfilename)
```
- .fits files:
```
from m4.ground import read_data
image = read_data.readFits_maskedImage(fits_file_path)
```

# Gestione delle immagini #

__Zernike__
```
from m4.ground import zernike
coeff, mat = zernike.zernikeFit(img, zernike_index_vector)
or
coeff, mat = zernikeFitAuxmask(img, auxmask, zernike_index_vector)
surf_image = zernike.zernikeSurface(img, coef, mat)
```

__ROI__
```
from m4.utils.roi import ROI
roi = ROI()
roi.automatical_roi_selection(image, segment_view, ref_mirror_in)
```
where segment_view and ref_mirror_in are boolean that decribes the input image.
The function return:
- roi_dx, roi_sx, roi_c, roi_rm in the case of segment_view=True and ref_mirror_in=True or False
- segRoiList, roi_rm in the case of segment_view=False and ref_mirror_in=True of only segRoiList for segment_view=False and ref_mirror_in=False

__Masks__
- Mask the image with a circular mask:
```
from m4.ground import geo
masked_image = draw_mask(start_image, centre_x, centre_y, radius)
```
- Extracting the parabola mask using markers:
```
from m4.utils.parabola_identification import ParabolaCirculaPupil
pz = ParabolaCirculaPupil()
circle_mask = pz.par_mask_on_ott(image)
```
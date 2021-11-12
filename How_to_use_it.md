vecchio file in CartellaBella

## Calibrazione ed allineamento OTT ##
 1) __Calibrazione della parabola con lo specchio di riferimento__
    - Scegliere un vettore di comandi con cui calibrare i gradi di libertà della parabola (PAR) e quello dello specchio di riferimento (RM).
 		  command_amp_vector = [par_piston, par_tip, par_tilt, rm_tip, rm_tilt] espressi in millimetri
 	    Nota: tenere presente che il metodo di calibrazione usato prevede che all'applicazione del comando par_tip corrisponda un'applicazione del comando rm_tip=-2.05*par_ti.
 		    stessa relazione sussiste tra par_tip ed rm_tip
    - Possibiltà di visualizzazione della matrice dei comandi di calibrazione prima dell'effetiva applicazione tramite il comando
 		  _main.showCommandMatrixBeforeCalibration(command_amp_vector)_
    - Il comando per effettuare la calibrazione è _main.calibrate_PARAndRM(ott, interf, n_frames, command_amp_vector, nPushPull)_ dove:  
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
 		_main.showCommandForParAndRmBeforeAlignement(ott, interf, tt_cal, n_images, zernike_to_be_corrected, dof_command_id)_ dove:  
 		ott e interf sono gli oggetti precedentemente creati nello start up  
 		tt_cal è la stringa contenente il tracking number di calibrazione che si vuole utilizzare  
 		n_images corrisponde al numero di frame da acquisire con l'interferometro per l'immagine da allineare  
 		zernike_to_be_corrected è il vettore che indica quali modi di zernike correggere (np.array([0,1,2,3,4]) indicano rispettivamente tip, tilt, fuoco, coma, coma)  
 		dof_command_id è il vettore che indica con quali gradi di libertà si vuole fare la correzione  
 		NOTA: gli zernike possono essere [0,1], [0,1,3,4], [0,1,2,3,4]  
 			e vanno in coppia con i dof [3,4], [1,2,3,4], [0,1,2,3,4]
    - Il comando per effettuare l'allineamento è _main.align_PARAndRM(ott, interf, tt_calib, n_images, zernike_to_be_corrected=None, dof_command_id=None)  
    - La funzione stampa i comandi applicati a PAR e RM, gli zernike calcolati e restituisce il traching number della misura di allineamento appena eseguita
	   	che contiene:  
	   	AlignmentInfo.fits in cui sono scritte tutte le info usate per allineare  
	   	le stesse informazioni separte (commandId.fits, IntMatModeVector.fits che sono gli zernike_to_be_corrected, PositionAndDeltaCommand.fits 
	   	e Zernike.fits che contiene gli zernike selezionati prima della correzione)
	   	Immagini prima e dopo l'allinemanto

 ## Accelerometri ##
 __Acquisizione dati__
 - ott.accelerometrs.acquireData(recording_seconds)
 	Nota: se non specificato acquisisce per 5 secondi
 - L'acquisizione restituisce il tracking number della misura appena effettuata dove sono salvati:
 	i dati, dt, plc_id, directions, time, plc_range, plc_totcounts, sensitivity
 	
 
 __Analisi dati__
 - from m4.analyzers.accelerometers_data_analyzer import AccelerometersDataAnalyzer
 - an = AccelerometersDataAnalyzer(tt)
 - Dopodichè si hanno a disposizione diversi modi di accedere/vedere i dati
  	- spe, freq = an.getSpecAndFreq() 
  	- data = an.datah5
  	- an.plot_power_spectrum()
  	- an.readAndShow()  
Nota: l'analisi non salva nulla
 
 ## Attività su mOTT ##

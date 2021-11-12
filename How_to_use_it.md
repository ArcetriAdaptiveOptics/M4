vecchio file in CartellaBella
  
## Calibrazione ed allineamento OTT ##
 * Calibrazione della parabola con lo specchio di riferimento:
 - Scegliere un vettore di comandi con cui calibrare i gradi di libertà della parabola (PAR) e quello dello specchio di riferimento (RM).
 		command_amp_vector = [par_piston, par_tip, par_tilt, rm_tip, rm_tilt] espressi in millimetri
 	  Nota: tenere presente che il metodo di calibrazione usato prevede che all'applicazione del comando par_tip corrisponda un'applicazione del comando rm_tip=-2.05*par_ti.
 		    stessa relazione sussiste tra par_tip ed rm_tip
 - Possibiltà di visualizzazione della matrice dei comandi di calibrazione prima dell'effetiva applicazione tramite il comando
 		main.showCommandMatrixBeforeCalibration(command_amp_vector)
 - Il comando per effettuare la calibrazione è main.calibrate_PARAndRM(ott, interf, n_frames, command_amp_vector, nPushPull) dove
 		ott e interf sono gli oggetti precedentemente creati nello start up
 		n_frame corrisponde al numero di frame da acquisire con l'interferometro per ogni misura da effttuare
 		command_amp_vector è il vettore da utilizzare per la calibrazione
 		nPushPull corrisponde al numero di volte che verrà effettuata la misura per ognuno dei 5 comandi applicati
 - La funzione di calibrazione ritorna il tracking number (tt_cal) della misura effettuata. Questo contiene:
 		CalibrationInfo.fits dove sono scritti insieme tutti i dati utilizzati per la calibrazione
 		Gli stessi dati separati in più file (CMat.fits, CommandAmplitude.fits, InteractionMatrix.fits, Mask.fits)
 		Tutti gli interferogrammi acquisiti
 * Allineamento della parabola e specchio di riferimento:
 - Il comando per effettuare l'allineamento è main.
 	
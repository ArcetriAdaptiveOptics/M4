# M4
***
Software per la calibrazione dello specchio deformabile M4

## Installare il pacchetto M4 da github
Per installarlo potendo modificarne il contenuto:
- git clone https://github.com/ArcetriAdaptiveOptics/M4.git
- cd M4
- pip install -e .

Nota: installare in modalita' editable (`-e`). I certificati OPC-UA (`m4/devices/*.der`, `*.pem`) non fanno parte del pacchetto installato: vengono letti dalla cartella dei sorgenti, e in loro assenza `import m4` termina con errore.

Dipendenze opzionali (extras): `pip install -e ".[gui,sound,test,docs]"`.

Il pacchetto `opticalib` non e' su PyPI e va installato da github:
- pip install git+https://github.com/pietroferraiuolo/labott.git

Per scaricare i cambiamenti presenti su github:
- git pull origin master

Il codice non piu' in uso si trova in `deprecated/` (vedi [deprecated/README.md](deprecated/README.md)); non viene installato.

## Start up
Il file di configurazione di `opticalib` viene letto dalla variabile d'ambiente `AOCONF`.
- Da terminale:
  - calpy -f <cartella esperimento> --create : crea cartelle e file di configurazione (`SysConfig/configuration.yaml`) e apre ipython
  - calpy -f <cartella esperimento> : apre ipython con la configurazione caricata
  - calpy -h : mostra tutte le opzioni
  - m4gui -f <cartella esperimento> : apre la GUI
- Nella shell aperta con `calpy`, eseguire:
  - import m4
  - from m4 import create_ott

    NOTA: all'import, `m4` aggiunge al file di configurazione la sezione `SYSTEM: simulated.devices` (solo i flag mancanti, impostati a True), con un flag per ciascun dispositivo:
    - dm, interferometer, accelerometers, angleRotator, m4Exapode, dp, parSlider, par, rmSlider, rm, tempSensors: True (simulato) or False

Leggere [How_to_use_it.md](https://github.com/ChiaraSelmi/M4/blob/master/How_to_use_it.md) per l'utilizzo di tutte le funzioni specifiche.

## COVERAGE
[![Build Status](https://travis-ci.org/codecov/sourcegraph-codecov.svg?branch=master)](https://codecov.io/gh/ChiaraSelmi/M4/)
[![codecov](https://codecov.io/gh/ChiaraSelmi/M4/branch/master/graph/badge.svg?token=DjC9LhW9FP)](https://codecov.io/gh/ChiaraSelmi/M4)
[![Documentation Status](https://readthedocs.org/projects/m4/badge/?version=latest)](https://m4.readthedocs.io/en/latest/?badge=latest)


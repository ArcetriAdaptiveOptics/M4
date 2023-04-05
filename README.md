# M4
***
Software per la calibrazione dello specchio deformabile M4

## Installare il pacchetto M4 da github
Per installarlo potendo modificarne il contenuto:
- git clone https://github.com/ChiaraSelmi/M4.git
- cd M4
- python setup.py install

Per installarlo senza accesso al contenuto:
- pip install git+https://github.com/ChiaraSelmi/M4.git

Per scaricare i cambiamenti presenti su github:
- git pull origin master

## Start up
- Aprire python (ipython --pylab) ed eseguire:
  - from m4.configuration import start
  - ott, interf, dm = start.create_ott('/....../M4/Data/SYSCONFData/Config.yaml')

    NOTA: file di configurazione da modificare all'occorrenza inserendo:
    - base_path: '/....../M4/Data'
    - simulated_interf: True or False
    - simulated_dm: True or False
    - simulated_accelerometers: True or False
    - simulated_angleRotator: True or False
    - simulated_m4Exapode: True or False
    - simulated_parSlider: True or False
    - simulated_par: True or False
    - simulated_rmSlider: True or False
    - simulated_rm: True or False
    - simulated_tempSensors: True or False
  - from m4 import main
  - main.start_log(logging_level) NOTA: ritorna il cammino del file di log

Leggere [How_to_use_it.md](https://github.com/ChiaraSelmi/M4/blob/master/How_to_use_it.md) per l'utilizzo di tutte le funzioni specifiche.

## COVERAGE
[![Build Status](https://travis-ci.org/codecov/sourcegraph-codecov.svg?branch=master)](https://codecov.io/gh/ChiaraSelmi/M4/)
[![codecov](https://codecov.io/gh/ChiaraSelmi/M4/branch/master/graph/badge.svg?token=DjC9LhW9FP)](https://codecov.io/gh/ChiaraSelmi/M4)
[![Documentation Status](https://readthedocs.org/projects/m4/badge/?version=latest)](https://m4.readthedocs.io/en/latest/?badge=latest)


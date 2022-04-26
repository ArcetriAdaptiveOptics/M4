from m4.configuration import start
import time
import os
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters
from m4.ground.opc_ua_controller import OpcUaController
from m4.configuration.config import fold_name
from m4.alignment import Alignment
from m4.ground import tracking_number_folder
from m4.ground import zernike
from m4.ground.timestamp import Timestamp
#from m4.utils import accelerometer_data as acc
import threading
import logging


from m4.configuration import start as startOTT


class parallelAcq():

    def __init__(self,time2acq=3):

        logging.basicConfig(level=logging.DEBUG,format='[%(levelname)s] (%(threadName)-10s) %(message)s')
        self._threads = []
        self._time2acq = time2acq
        #self._acc = acc.Accelerometers() //to update with new class
        self._ott, self._interf =  startOTT.create_ott()
        self._accTN = ""
        self._interfTN = ""
        self._zabbixTN = ""
        
        self.resetThreads()

    def resetThreads(self):
        self._threads=[]
        t1 = threading.Thread(name='worker_acc', target=self.worker_acc)
        self._threads.append(t1)
        t2 = threading.Thread(name='worker_i4d', target=self.worker_i4d)
        self._threads.append(t2)
        t3 = threading.Thread(name='worker_zbx', target=self.worker_zbx)
        self._threads.append(t3)

    def startMeas(self):
        self._accTN = ""
        self._interfTN = ""
        self._zabbixTN = ""
        for tt in self._threads:
            tt.start()
        self.resetThreads()

    def worker_acc(self):
        logging.debug("start acc acquisition")
        self._accTN = self._ott.accelerometer.acquireData(self._time2acq)
        logging.debug("stop acc acquisition")
    

    def worker_i4d(self):
        logging.debug("start i4d")
        nframes = self._time2acq/self._nframespersec
        self._interfTN=self._interf.acquire_phasemap(nframes)
        logging.debug("stop i4d")
    
    def worker_zbx(self):
        logging.debug("start zbx")
        #self._ott.temperatures.getData()  (time2acq)
        time.sleep(self._time2acq+3)
        logging.debug("stop zbx")



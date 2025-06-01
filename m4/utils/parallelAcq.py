""" Da dove viene? A cosa serve?
RIVEDERE
"""

import time
import threading
import logging


class ParallelAcq:

    def __init__(self, ott, interf, time2acq=3):

        logging.basicConfig(
            level=logging.DEBUG,
            format="[%(levelname)s] (%(threadName)-10s) %(message)s",
        )
        self._threads = []
        self._time2acq = time2acq
        self._nframespersec = "?"
        # self._acc = acc.Accelerometers() //to update with new class
        self._ott = ott
        self._interf = interf
        self._accTN = ""
        self._interfTN = ""
        self._zabbixTN = ""

        self.resetThreads()

    def resetThreads(self):
        self._threads = []
        t1 = threading.Thread(name="worker_acc", target=self.worker_acc)
        self._threads.append(t1)
        t2 = threading.Thread(name="worker_i4d", target=self.worker_i4d)
        self._threads.append(t2)
        t3 = threading.Thread(name="worker_zbx", target=self.worker_zbx)
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
        nframes = self._time2acq / self._nframespersec
        self._interfTN = self._interf.acquire_phasemap(nframes)
        logging.debug("stop i4d")

    def worker_zbx(self):
        logging.debug("start zbx")
        # self._ott.temperatures.getData()  (time2acq)
        time.sleep(self._time2acq + 3)
        logging.debug("stop zbx")

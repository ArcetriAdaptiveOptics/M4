'''
Authors
  - C. Selmi: written in 2021
'''
import os
import yaml
import serial
import time

class SerialTimeoutException(Exception):
    def __init__(self, value=-1):
        print ("Missing response from serial after %i iterrations" % value)

class TunaFilt:

    def __init__(self, confLctfFile):
        self.currfilt = CurrentFilterReader(confLctfFile)
        self.ser = None

    def _pollSerial(self):
        nw = 0
        nw0 = 0
        it = 0
        while True:
            nw = self.ser.inWaiting()
            it = it + 1
            time.sleep(0.01)
            if (nw >0) and (nw0==nw) or (it==10000):
                break
            nw0 = nw
        if nw == 0:
            raise SerialTimeoutException(it)
        else:
            return nw

    def connect(self, port='/dev/ttyUSB1', speed=115200):
        if self.ser is None:
            self.ser = serial.Serial(port, speed)
            out = self.get_status()
            return out
        else:
            print ("Already connected")

    def disconnect(self):
        if self.ser is not None:
            self.ser.close()
            self.ser = None

    def get_wl(self):
        cmd = bytes(self.currfilt.read_wl, 'utf-8')
        tmp = self.ser.write(cmd)
        nw = self._pollSerial()
        out = self.ser.read(self.ser.inWaiting()) 
        return out

    def set_wl(self, wl):
        if wl < 650 or wl > 1100:
            raise BaseException()
        cmd = bytes(self.currfilt.write_wl % wl, 'utf-8')
        tmp = self.ser.write(cmd)
        nw = self._pollSerial()
        out = self.ser.read(self.ser.inWaiting()) 
        return out

    def reset(self):
        cmd = bytes(self.currfilt.reset, 'utf-8')
        tmp = self.ser.write(cmd) #modo di esprire la @ con python3 (con b ritorna quello di python2)
        nw = self._pollSerial()
        out = self.ser.read(self.ser.inWaiting()) 
        return out

    def get_status(self):
        cmd = bytes(self.currfilt.read_status, 'utf-8')
        tmp = self.ser.write(cmd)
        nw = self._pollSerial()
        out = self.ser.read(self.ser.inWaiting()) 
        return out

    def isbusy(self):
        cmd = bytes(self.currfilt.busy_check, 'utf-8')
        tmp = self.ser.write(cmd)
        nw = self._pollSerial()
        out = self.ser.read(self.ser.inWaiting()) 
        return out

    def cancel(self):
        cmd = bytes(self.currfilt.escape, 'utf-8')
        tmp = self.ser.write(cmd)
        nw = self._pollSerial()
        out = self.ser.read(self.ser.inWaiting()) 
        return out

class CurrentFilterReader():
    ''' class for reading data from yaml file
    '''

    def __init__(self, confLctfFile):
        ''' The constructor'''
        with open(confLctfFile) as file:
            self._currFilt = yaml.load(file, Loader=yaml.FullLoader)

    @property
    def write_wl(self):
        return self._currFilt['WRITE_WL']

    @property
    def read_wl(self):
        return self._currFilt['READ_WL']

    @property
    def read_status(self):
        return self._currFilt['READ_STATUS']

    @property
    def busy_check(self):
        return self._currFilt['BUSY_CHECK']

    @property
    def escape(self):
        return self._currFilt['ESCAPE']

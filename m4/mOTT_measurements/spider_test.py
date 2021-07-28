'''
Authors
  - C. Selmi: written in 2021
'''
import os
from m4.ground import tracking_number_folder
from m4.devices.i4d import I4D
from m4.configuration.ott_parameters import Interferometer

class SpiderTest():

    def __init__(self):
        """The constructor """
        self._interf = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)
        self._rawFramesDirectory = 'D:/M4/SpiderTest/Raw'
        self._configDirectory = 'D:/config/'
        self._measurementsDirectory = 'D:/M4/SpiderTest/h5'

    def rawAcquisition(self, numberOfFrames):
        dove, tt = tracking_number_folder.createFolderToStoreMeasurements(self._rawFramesDirectory)
        self._interf.burstFramesToSpecificDirectory(dove, numberOfFrames)
        return tt

    def convertRawWithSelectedConfig(self, tt_rawFrame, configurationName):
        confLabel = configurationName.split('.')[0]
        configurationPath = os.path.join(self._configDirectory, configurationName)
        self._interf.loadConfiguration(configurationPath)
        rawFramesFolder = os.path.join(self._rawFramesDirectory, tt_rawFrame)
        tt_conf = tt_rawFrame + '_' + confLabel
        measurementsFolder = os.path.join(self._measurementsDirectory, tt_conf)
        self._interf.convertRawFramesInDirectoryToMeasurementsInDestinationDirectory(measurementsFolder,
                                                                                rawFramesFolder)
        return tt_conf

    def seriesConvertion(self, tt_rawFrame, configurationNameList):
        for name in configurationNameList:
            tt_conf = self.convertRawWithSelectedConfig(tt_rawFrame, name)
            print(tt_conf)



###
def main0713TestSpider(numberOfFrames, measurementsDirectory):
    from m4.devices.i4d import I4D
    from m4.configuration.ott_parameters import Interferometer
    interf = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)
    rawFramesDirectory = "C:/Users/PhaseCam/Desktop/prova"
    interf.burstFramesToSpecificDirectory(rawFramesDirectory, numberOfFrames)

    configurationPath = 'D:/config/20210708-spider.4Dini'
    configurationPath1 = 'D:/config/20210708-spider.4Dini'
    configurationPath2 = 'D:/config/20210708-spider2.4Dini'
    interf.loadConfiguration(configurationPath)

    interf.convertRawFramesInDirectoryToMeasurementsInDestinationDirectory(measurementsDirectory,
                                                                            rawFramesDirectory)
    return measurementsDirectory

'''
Authors
  - C. Selmi: written in 2021
'''
import os
from m4.ground import tracking_number_folder
from m4.devices.i4d import I4D
from m4.configuration.ott_parameters import Interferometer
import h5py
from astropy.io import fits as pyfits
import glob
from m4.ground import read_data

class SpiderTest():

    def __init__(self):
        """The constructor """
        self._interf = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)
        self._rawFramesDirectory = 'D:/M4/SpiderTest/Raw'
        self._configDirectory = 'D:/config/'
        self._measurementsDirectory = 'D:/M4/SpiderTest/h5'
        self._measurementsDirectoryInM4wsPc = '/mnt/data/M4/Data/SpiderTest/h5'
        self._fitsDirectory = 'D:/M4/SpiderTest/fits'
        self._mask_file_name = '/mnt/data/M4/Data/SpiderTest/h5/islandsmask.fits'

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

    def _convertAndSaveSingle4DDataToSlimFits(self, filePath4D, filePathfits):
        '''
        Parameters
        ----------
        filePath4D: string
            complete data file path
        filePathfits: string
            complete fits file path
        '''
        file = h5py.File(filePath4D, 'r')
        data = file.get('/Measurement/SurfaceInWaves/Data')
        meas = data[()]  * 632.8e-9
        pyfits.writeto(filePathfits, meas)

    def convert4DDataToSlimFits(self, folderName4D):
        '''
        Parameters
        ----------
        folderName4D: string
            folder name (tracking number) for .4D data
        '''
        folderPath4D = os.path.join(self._measurementsDirectoryInM4wsPc, folderName4D)
        list = glob.glob(os.path.join(folderPath4D,'*.4D'))
        list.sort()
        n_images = len(list)

        for i in range(n_images):
            number_fileName4D = '%d.4D' %i
            filePath4D = os.path.join(folderPath4D, number_fileName4D)
            fitsFileName = '%04d.fits' %i
            filePathFits = os.path.join(self._fitsDirectory, folderName4D, fitsFileName)
            self._convertAndSaveSingle4DDataToSlimFits(filePath4D, filePathFits)

    def analysis(self, dataFitsFolfer):
        list = glob.glob(os.path.join(dataFitsFolfer,'*.fits'))
        list.sort()
        n_images = len(list)
        mask = read_data.readFits_data(self._mask_file_name) #shape [1000, 1000, 3]

        for i in range(n_images):
            if i % 2 == 0:
                fits_file_path = os.path.join(dataFitsFolfer, '%04d.fits' %i)
                ima0 = read_data.readFitsSlimImage(fits_file_path)
                fits_file_path = os.path.join(dataFitsFolfer, '%04d.fits' %(i+1))
                ima1 = read_data.readFitsSlimImage(fits_file_path)
                diff = ima1 -ima0


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

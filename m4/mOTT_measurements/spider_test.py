'''
Authors
  - C. Selmi: written in 2021
'''
import os
import numpy as np
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
        self._fitsDirectory = '/mnt/data/M4/Data/SpiderTest/fits'
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
        pyfits.writeto(filePathfits, meas, overwrite=True)

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
            directory = os.path.join(self._fitsDirectory, folderName4D)
            if not os.path.exists(directory):
                os.makedirs(directory)
            filePathFits = os.path.join(directory, fitsFileName)
            self._convertAndSaveSingle4DDataToSlimFits(filePath4D, filePathFits)

    def analysis(self, dataFitsFolfer):
        list = glob.glob(os.path.join(self._fitsDirectory, dataFitsFolfer,'*.fits'))
        list.sort()
        n_images = len(list)
        mask = read_data.readFits_data(self._mask_file_name) #shape [1000, 1000, 3]

        cp_mat = np.zeros((np.int(n_images/2), mask.shape[2]))
        cr_mat = np.zeros((np.int(n_images/2), mask.shape[2]))
        diff_list = []
        for i in range(n_images):
            if i % 2 == 0:
                fits_file_path = os.path.join(self._fitsDirectory, dataFitsFolfer, '%04d.fits' %i)
                ima0 = read_data.readFitsSlimImage(fits_file_path)
                fits_file_path = os.path.join(self._fitsDirectory, dataFitsFolfer, '%04d.fits' %(i+1))
                ima1 = read_data.readFitsSlimImage(fits_file_path)
                diff = ima1 -ima0
                diff_list.append(diff)
                
                for j in range(mask.shape[2]):
                    k = np.int(i/2)
                    prod = np.ma.masked_array(diff, mask=np.invert(mask[:,:,j].astype(bool)))
                    cp = prod.mean()
                    if cp > 200e-09 or cp < -200e-09:
                        cp_mat[k, j]=1
                    cr = prod.std()
                    if cr > 20e-09 or cr < -20e-09:
                        cr_mat[k, j]=1
        return diff_list, cp_mat, cr_mat

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

def main0908ConvertionAndAnalysis():
    sp = SpiderTest()
    tt_list = ['20210721_103804_20210720-spider6',
               '20210721_104344_20210720-spider6',
               '20210721_104858_20210720-spider6']
    for tt in tt_list:
        print('Converting tt = %s' %tt)
        sp.convert4DDataToSlimFits(tt)
        print('Analysing tt = %s' %tt)
        diff, cp, cr = sp.analysis(tt)
        print(np.sum(cp, 0).max())



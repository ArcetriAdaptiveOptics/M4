'''
Authors
  - C. Selmi: written in 2024
'''
import unittest
import os
import mock
import numpy as np
import logging
from m4.devices.i4d import I4D

class Testi4d(unittest.TestCase):

    def setUp(self):
        self._setUpLogging()
        self._ip = '123.345.456.123'
        self._port = 1111
        self._interferometer = I4D(
            self._ip, self._port)

        self._interferometer._readJsonData = self._my_readJsonData

        self._calls = []

    def _setUpLogging(self):
        FORMAT = '%(asctime)s %(levelname)s %(message)s'
        logging.basicConfig(level=logging.DEBUG, format=FORMAT)
        self._logger = logging.getLogger('i4d library')

    ### DATA PROXY ###
    def _expected_url_GetFeatureAnalysisResults(self):
        return 'http://%s:%i/DataService/GetFeatureAnalysisResults' % (
            self._ip, self._port)

    def _expected_url_GetFirstNineZernikeTerms(self):
        return 'http://%s:%i/DataService/GetFirstNineZernikeTerms' % (
            self._ip, self._port)

    def _expected_url_GetFringeAmplitudeData(self):
        return 'http://%s:%i/DataService/GetFringeAmplitudeData' % (
            self._ip, self._port)

    def _expected_url_GetIntensityData(self):
        return 'http://%s:%i/DataService/GetIntensityData' % (
            self._ip, self._port)

    def _expected_url_GetInterferogram(self):
        return 'http://%s:%i/DataService/GetInterferogram' % (
            self._ip, self._port)

    def _expected_url_GetMeasurementInfo(self):
        return 'http://%s:%i/DataService/GetMeasurementInfo' % (
            self._ip, self._port)

    def _expected_url_GetModulationData(self):
        return 'http://%s:%i/DataService/GetModulationData/' % (
            self._ip, self._port)

    def _expected_url_GetPhaseStepCalculatorResults(self):
        return 'http://%s:%i/DataService/GetPhaseStepCalculatorResults' % (
            self._ip, self._port)

    def _expected_url_GetSurfaceData(self):
        return 'http://%s:%i/DataService/GetSurfaceData' % (
            self._ip, self._port)

    def _expected_url_GetUnprocessedSurfaceData(self):
        return 'http://%s:%i/DataService/GetUnprocessedSurfaceData' % (
            self._ip, self._port)

    def _expected_url_SaveDataToDisk(self):
        return 'http://%s:%i/DataService/SaveDataToDisk/' % (
            self._ip, self._port)

### SYSTEM PROXY ###
    def _expected_url_TakeSingleMeasurement(self):
        return 'http://%s:%i/SystemService/TakeSingleMeasurement/' % (
            self._ip, self._port)
        
    def _expected_url_SetDetectorMask(self):
        return 'http://%s:%i/SystemService/SetDetectorMask' % (
            self._ip, self._port)

    def _expected_url_GetSystemInfo(self):
        return 'http://%s:%i/SystemService/GetSystemInfo' % (
            self._ip, self._port)

    def _expected_url_ConvertRawFramesInDirectoryToMeasurementsInDestinationDirectory(self):
        return 'http://%s:%i/SystemService/ConvertRawFramesInDirectoryToMeasurementsInDestinationDirectory' % (
            self._ip, self._port)

    def _expected_url_SetTriggerMode(self):
        return 'http://%s:%i/SystemService/SetTriggerMode' % (
            self._ip, self._port)

    def _expected_url_TakeAveragedMeasurement(self):
        return 'http://%s:%i/SystemService/TakeAveragedMeasurement' % (
            self._ip, self._port)

    def _expected_url_LoadConfiguration(self):
        return 'http://%s:%i/SystemService/LoadConfiguration' % (
            self._ip, self._port)

### FRAME BURST PROXY ###
    def _expected_url_BurstFramesToSpecificDirectory(self):
        return 'http://%s:%i/FrameBurstService/BurstFramesToSpecificDirectory' % (
            self._ip, self._port)

    def _expected_url_BurstFramesToDisk(self):
        return 'http://%s:%i/FrameBurstService/BurstFramesToDisk' % (
            self._ip, self._port)

    def _my_readJsonData(self, url, data=None, timeout=None):
        '''
        This one register the calls to the WFC and return the proper
        structure according to the requested url
        '''
        this_call = {'url': url, 'data': data, 'timeout': timeout}
        self._logger.info("Got call %s" % this_call)
        self._calls.append(this_call)
        # parse url to mimick proper return values
        ret = {}
        # DATA PROXY #
        if url == self._expected_url_GetFeatureAnalysisResults():
            ret['Features'] = (1, 2, 3)
            ret['FractionalFeatureArea'] = 4
            ret['NumberOfFeatures'] = 5
            ret['TotalFeatureAreaInSquareMicrons'] = 6
            return ret
        elif url == self._expected_url_GetFirstNineZernikeTerms():
            return np.array(np.arange(10))
        elif url == self._expected_url_GetFringeAmplitudeData()\
            or url == self._expected_url_GetIntensityData():
            ret['Width'] = 2
            ret['Height'] = 3
            ret['PixelSizeInMicrons'] = 13.2
            ret['Data'] = (12, 13, 14, 15, 16, 17)
            return ret
        elif url == self._expected_url_GetInterferogram():
            pass
        elif url == self._expected_url_GetMeasurementInfo():
            ret['AverageFringeAmplitude'] = 2
            ret['AverageIntensity'] = 3
            ret['AverageModulation'] = 4
            ret['FringeAmpThresholdPercentage'] = 5
            ret['IntensityThresholdPercentage'] = 6
            ret['ModulationThresholdPercentage'] = 7
            ret['NumberOfSamples'] = 8
            ret['NumberOfValidPixels'] = 9
            ret['PathMatchPositionInMM'] = '../path'
            ret['RMSInNM'] = 10
            ret['UserSettingsFilePath'] = '.../path2'
            ret['WavelengthInNM'] = 11
            ret['Wedge'] = 12
            return ret
        elif url == self._expected_url_GetModulationData()\
            or url == self._expected_url_GetSurfaceData()\
            or url == self._expected_url_GetUnprocessedSurfaceData():
            ret['Width'] = 2
            ret['Height'] = 3
            ret['PixelSizeInMicrons'] = 13.2
            ret['Data'] = (12, 13, 14, 15, 16, 17)
            return ret
        elif url == self._expected_url_GetPhaseStepCalculatorResults():
            ret['AveragePhaseStepInDegrees'] = 1
            ret['Height'] = 2
            ret['PhaseStepsInDegrees'] = 3
            ret['Width'] = 4
            return ret
        # SYSTEM PROXY #
        elif url == self._expected_url_TakeSingleMeasurement():
            ret['Width'] = 2
            ret['Height'] = 3
            ret['PixelSizeInMicrons'] = 12.3
            ret['Data'] = (12, 13, 14, 15, 16, 17)
            return ret
        elif url == self._expected_url_GetSystemInfo():
            ret['SystemSerialNumber'] = 'ciao'
            return ret
        elif url == self._expected_url_TakeAveragedMeasurement():
            ret['Width'] = 2
            ret['Height'] = 3
            ret['PixelSizeInMicrons'] = 12.3
            ret['Data'] = (12, 13, 14, 15, 16, 17)
            return ret
        # FRAME BURST PROXY #
        else:
            self._logger.warning("unknown url %s" % this_call)
            return

    def tearDown(self):
        pass

### DATA PROXY ###
    def test_get_feature_analysis_results(self):
        ret = self._interferometer.getFeatureAnalysisResults()
        wanted_url = self._expected_url_GetFeatureAnalysisResults()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_first_nine_zernike_terms(self):
        ret = self._interferometer.getFirstNineZernikeTerms()
        wanted_url = self._expected_url_GetFirstNineZernikeTerms()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_fringe_amplitude_data(self):
        ret = self._interferometer.getFringeAmplitudeData()
        wanted_url = self._expected_url_GetFringeAmplitudeData()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_intensity_data(self):
        ret = self._interferometer.getIntensityData()
        wanted_url = self._expected_url_GetIntensityData()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_interferogram(self):
        index = 4
        ret = self._interferometer.getInterferogram(index)
        wanted_url = self._expected_url_GetInterferogram()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_measurement_info(self):
        ret = self._interferometer.getMeasurementInfo()
        wanted_url = self._expected_url_GetMeasurementInfo()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_modulation_data(self):
        ret = self._interferometer.dataServiceGetModulationData()
        wanted_url = self._expected_url_GetModulationData()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_phase_step_calculator_results(self):
        ret = self._interferometer.getPhaseStepCalculatorResults()
        wanted_url = self._expected_url_GetPhaseStepCalculatorResults()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_surface_data(self):
        ret = self._interferometer.getSurfaceData()
        wanted_url = self._expected_url_GetSurfaceData()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_unprocessed_surface_data(self):
        ret = self._interferometer.getUnprocessedSurfaceData()
        wanted_url = self._expected_url_GetUnprocessedSurfaceData()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_save_data_to_disk(self):
        path = '.../path'
        ret = self._interferometer.saveDataToDisk(path)
        wanted_url = self._expected_url_SaveDataToDisk()

### SYSTEM PROXY ###
    def test_take_single_measurement(self):
        wv = self._interferometer.takeSingleMeasurement()
        wanted_url = self._expected_url_TakeSingleMeasurement()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_set_detector_mask(self):
        mask = np.zeros((5, 5))
        self._interferometer.setDetectorMask(mask)
        wanted_url = self._expected_url_SetDetectorMask()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_get_system_info(self):
        ret = self._interferometer.getSystemInfo()
        wanted_url = self._expected_url_GetSystemInfo()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_convert_raw_frames_in_directory_to_measurements_in_destination_directory(self):
        measurementsDirectory = '.../path1'
        rawFramesDirectory = '.../path2'
        self._interferometer.convertRawFramesInDirectoryToMeasurementsInDestinationDirectory(measurementsDirectory,
                                                                                            rawFramesDirectory)
        wanted_url = self._expected_url_ConvertRawFramesInDirectoryToMeasurementsInDestinationDirectory()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_set_trigger_mode(self):
        trigger = True
        self._interferometer.setTriggerMode(trigger)
        wanted_url = self._expected_url_SetTriggerMode()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_take_averaged_measurement(self):
        ret = self._interferometer.takeAveragedMeasurement(3)
        wanted_url = self._expected_url_TakeAveragedMeasurement()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_load_configuration(self):
        configurationPath = '.../configurationFile'
        self._interferometer.loadConfiguration(configurationPath)
        wanted_url = self._expected_url_LoadConfiguration()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

### FRAME BURST PROXY ###
    def test_burst_frames_to_specific_directory(self):
        directory = '.../dir'
        self._interferometer.burstFramesToSpecificDirectory(directory, 10)
        wanted_url = self._expected_url_BurstFramesToSpecificDirectory()
        self.assertEqual(self._calls[-1]['url'], wanted_url)

    def test_burst_frames_to_disk(self):
        self._interferometer.burstFramesToDisk(10)
        wanted_url = self._expected_url_BurstFramesToDisk()
        self.assertEqual(self._calls[-1]['url'], wanted_url)


if __name__ == "__main__":
    unittest.main()

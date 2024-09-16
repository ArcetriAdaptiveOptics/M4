import unittest
import numpy as np
from unittest.mock import MagicMock, patch
from m4.opt_alignment import Alignment

class TestAlignment(unittest.TestCase):

	def setUp(self):
		# Mock mechanical and acquisition devices
		self.mock_mechanical_devices = MagicMock()
		self.mock_acquisition_devices = MagicMock()

		# Mock the readFits_data function
		self.mock_readFits_data = patch('m4.opt_alignment.rd.readFits_data').start()
		self.mock_readFits_data.return_value = np.eye(3)

		# Mock the saveFits_data function
		self.mock_saveFits_data = patch('m4.opt_alignment.rd.saveFits_data').start()

		# Mock the get_callables function
		self.mock_get_callables = patch('m4.opt_alignment.Alignment._get_callables').start()
		self.mock_get_callables.return_value = [MagicMock()]

		# Mock the mac module attributes
		self.mock_mac = patch('m4.opt_alignment.mac').start()
		self.mock_mac.devices_move_calls = ['move']
		self.mock_mac.devices_read_calls = ['read']
		self.mock_mac.ccd_acquisition = ['acquire']
		self.mock_mac.names = ['device1']
		self.mock_mac.dof = [0, 1, 2]
		self.mock_mac.cmdDof = 3
		self.mock_mac.slices = [slice(0, 3)]
		self.mock_mac.base_read_data_path = '/mock/path/to/read'
		self.mock_mac.base_write_data_path = '/mock/path/to/write'

		# Initialize the Alignment class
		self.alignment = Alignment(self.mock_mechanical_devices, self.mock_acquisition_devices)

	def tearDown(self):
		patch.stopall()

	def test_correct_alignment(self):
		# Test the correct_alignment method
		self.assertIsNone(self.alignment.correct_alignment())

	@patch('m4.opt_alignment.rd.saveFits_data')
	def test_calibrate_alignment(self, mock_saveFits_data):
		# Test the calibrate_alignment method
		cmdAmp = 1.0
		template = [1, -1, 1]
		n_repetitions = 1
		intMat = self.alignment.calibrate_alignment(cmdAmp, template, n_repetitions)
		self.assertIsInstance(intMat, np.ndarray)
		mock_saveFits_data.assert_called_once_with('/mock/path/to/write/intMat.fits', intMat)

	def test_images_production(self):
		# Test the images_production method
		template = [1, -1, 1]
		n_repetitions = 1
		results = self.alignment.images_production(template, n_repetitions)
		self.assertIsInstance(results, list)

	def test_create_intMat(self):
		# Test the create_intMat method
		imglist = [np.zeros((3, 3))]
		intMat = self.alignment.create_intMat(imglist)
		self.assertIsInstance(intMat, np.ndarray)

	@patch('m4.opt_alignment.Alignment._get_callables')
	def test_apply_command(self, mock_get_callables):
		# Test the apply_command method
		fullCmd = np.zeros(3)
		self.alignment.apply_command(fullCmd)
		mock_get_callables.assert_called()

	def test_read_positions(self):
		# Test the read_positions method
		positions = self.alignment.read_positions()
		self.assertIsInstance(positions, list)

	@patch('m4.opt_alignment.rd.readFits_maskedImage')
	def test_reload_parabola_tn(self, mock_readFits_maskedImage):
		# Test the reload_parabola_tn method
		mock_readFits_maskedImage.return_value = MagicMock(mask=np.zeros((3, 3)))
		filepath = '/mock/path/to/parabola.fits'
		message = self.alignment.reload_parabola_tn(filepath)
		self.assertIsInstance(message, str)

	def test_extract_cmds_to_apply(self):
		# Test the _extract_cmds_to_apply method
		fullCmd = np.zeros(3)
		device_commands = self.alignment._extract_cmds_to_apply(fullCmd)
		self.assertIsInstance(device_commands, list)

	def test_img_acquisition(self):
		# Test the _img_acquisition method
		k = 0
		template = [1, -1, 1]
		imglist = self.alignment._img_acquisition(k, template)
		self.assertIsInstance(imglist, list)

	def test_push_pull_redux(self):
		# Test the _push_pull_redux method
		imglist = [np.zeros((3, 3)), np.zeros((3, 3))]
		template = [1, -1, 1]
		image = self.alignment._push_pull_redux(imglist, template)
		self.assertIsInstance(image, np.ndarray)

	def test_get_callables(self):
		# Test the _get_callables method
		device = MagicMock()
		callables = ['method1', 'method2']
		functions = self.alignment._get_callables(device, callables)
		self.assertIsInstance(functions, list)

if __name__ == '__main__':
	unittest.main()
'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import os
import numpy as np
from skimage.draw import circle
from m4.ground import geo
from m4.ground import smooth_function
from m4.ground import logger_set_up
from m4.ground import tracking_number_folder
from m4.ground import zernike
from m4.ground import read_data

TESTDATA_DIR = os.path.dirname(__file__)


class Test(unittest.TestCase):

    @unittest.skip('Mettere dei file')
    def testFindDirectory(self):
        tt = '2021...'
        tracking_number_folder.findTrackingNumberPath(tt)

    def testGeometry(self):
        img = np.random.rand(500, 500)
        masked_ima = geo.draw_mask(img, 250, 250, 50)
        geo.qpupil(masked_ima)
        geo.rotate(img, 30)

    def testGUI(self):
        pass

    @unittest.skip('Dove mettere il file')
    def testLogger(self):
        import tempfile
        file = tempfile.NamedTemporaryFile(mode='w', delete=False)
        file_path = os.path.join(file, 'mylog')
        logging_level = 0
        logger_set_up.set_up_logger(file_path.name, logging_level)

    def testReadData(self):
        ic = read_data.InterferometerConverter()
        file_path = os.path.join(TESTDATA_DIR, 'img_0046.h5')
        ima = ic.from4D(file_path)
        file_path = os.path.join(TESTDATA_DIR, '0.4D')
        ima2 = ic.fromNew4D(file_path)

    def testSmoothFunction(self):
        data = np.arange(100)
        smooth_function.smooth(data, 4)

    @unittest.skip('Dove mettere il file')
    def testTtFolder(self):
        store_in_folder = '?'
        dove, tt = tracking_number_folder.createFolderToStoreMeasurements(store_in_folder)

    def testZernike(self):
        img = np.random.rand(500, 500)
        mask = np.ones((500, 500), dtype=bool)
        rr, cc = circle(250, 250, 100)
        mask[rr, cc] = 0
        masked_ima = np.ma.masked_array(img, mask=mask)

        coef, mat = zernike.zernikeFit(masked_ima, np.arange(10) + 1)
        zernike.zernikeSurface(masked_ima, coef, mat)

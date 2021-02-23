'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import mock
from m4.ground.interface_4D import comm4d


class Test(unittest.TestCase):


    def setUp(self):
        self.interf = comm4d()


    def tearDown(self):
        del self.interf

    @unittest.skip('Non trova il modulo e ott non ok')
    @mock.patch('oaautils.i4d', autospect=True)
    def testComm4D(self):
        ott = mock.MagicMock()
        image = self.interf.acq4d(1, 0, ott)
        image2 = self.interf.acq4d(7, 1)

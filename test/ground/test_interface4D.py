'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import mock


class Test(unittest.TestCase):


    @mock.patch('oaautils.i4d', autospect=True)
    def setUp(self, mock_lib):
        from m4.ground.interface_4D import comm4d
        self.interf = comm4d()


    def tearDown(self):
        del self.interf

    @unittest.skip('Non trova il modulo oaautils e ott non ok')
    def testComm4D(self):
        ott = mock.MagicMock()
        image = self.interf.acq4d(1, 0, ott)
        image2 = self.interf.acq4d(7, 1)

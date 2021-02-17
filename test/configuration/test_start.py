'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
from m4.configuration import start


class Test(unittest.TestCase):

    @unittest.skip('Mancano i file per costruire il simulatore')
    def testStart(self):
        ott, interf = start.create_ott()

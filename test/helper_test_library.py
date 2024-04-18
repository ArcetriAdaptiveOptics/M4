'''
Created on 19 mag 2021

@author: lbusoni
'''
import os

def testDataRootDir():
    return os.path.join(os.path.dirname(__file__), 'data')

def testTmpDir():
    return os.path.join(os.path.dirname('/'), 'tmp')

'''
Created on 19 mag 2021

@author: lbusoni
'''
import os
import pathlib3x as pathlib
import tempfile


def testDataRootDir():
    return os.path.join(os.path.dirname(__file__), 'data')


def testTmpDir():
    return os.path.join(os.path.dirname('/'), 'tmp')


def testDataRootPath():
    return pathlib.Path(pathlib.Path(__file__).parent, 'data')


def testTmpPath():
    return pathlib.Path(pathlib.Path(tempfile.gettempdir()), 'tmp')

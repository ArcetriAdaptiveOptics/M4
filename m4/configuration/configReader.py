'''
Authors
  - C. Selmi: written in 2021
'''
import os
import yaml

class Configuration():

    def __init__(self, confFile):
        with open(confFile) as file:
            self._conf = yaml.load(file, Loader=yaml.FullLoader)
            self._basepath = None

    @staticmethod
    def _basePath():
        return Configuration._basepath

    @property
    def BASE_PATH(self):
        if 'base_path' in self._conf.keys():
            self._basepath = self._conf['base_path']
        else:
            self._basepath = '/mnt/m4storage/Data'
        return self._basepath

    @BASE_PATH.getter
    def BASE_PATH(self):
        return self._basepath

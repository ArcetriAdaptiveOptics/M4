__version__= "$Id$"

from argos.util.timestamp import Timestamp


class Tag(object):

    def __init__(self, tagString):
        assert tagString.count("_") > 0
        self._tagString= tagString


    def getDayAsString(self):
        return self._tagString.split("_")[0]


    def __str__(self):
        return self._tagString


    @staticmethod
    def createTag():
        return Tag(Timestamp().now())

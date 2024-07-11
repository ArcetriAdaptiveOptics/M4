__version__= "$Id$"

from argos.util.fits_header import FitsHeader

import abc


def joinEntries(aList):
    return '.'.join(aList)


class Snapshotable():

    __metaclass__= abc.ABCMeta

    @abc.abstractmethod
    def getSnapshot(self, prefix):
        assert False


    @staticmethod
    def _isSuitableForFITSHeader(value):
        return value is not None


    @staticmethod
    def _truncateStringIfAny(key, value):
        maxValueLenInChars= 67 - len(key)
        if len(value) > maxValueLenInChars:
            return value[0:maxValueLenInChars]
        return value


    @staticmethod
    def _updateHeader(hdr, key, value):
        MAX_KEY_LEN_CHARS= 59
        assert len(key) <= MAX_KEY_LEN_CHARS
        if isinstance(value, str):
            value= Snapshotable._truncateStringIfAny(key, value)
        hdr.update('hierarch '+ key, value)


    @staticmethod
    def asFITSHeader(snapshotDictionary):
        # TODO 2014-05-14 mk:  Replace long header lines
        #                      with a level of indirection!
        hdr= FitsHeader()
        for k in sorted(snapshotDictionary.iterkeys()):
            value= snapshotDictionary[k]
            if Snapshotable._isSuitableForFITSHeader(value):
                Snapshotable._updateHeader(hdr, k, value)
        return hdr


    @staticmethod
    def prepend(prefix, snapshotDict):
        assert len(prefix) > 0, "Prefix length must be greater than zero"
        for each in snapshotDict.keys():
            value= snapshotDict[each]
            del snapshotDict[each]
            newKey= prefix + "." + each
            snapshotDict[newKey]= value

        return snapshotDict


    @staticmethod
    def fromFITSHeader(hdr):
        # TODO 2014-05-14 mk:  Restore snapshots with long values!
        snapshot= {}
        for each in hdr:
            snapshot[each]= hdr[each]

        return snapshot


    @staticmethod
    def removeEntriesWithValueNone(snapshotDict):
        for each in snapshotDict.keys():
            if snapshotDict[each] is None:
                del snapshotDict[each]

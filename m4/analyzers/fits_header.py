__version__ = "$Id$"

import pyfits


class FitsHeader(pyfits.Header):

    def __init__(self):
        pyfits.Header.__init__(self)

    def update(self, key, value):
        if self._isOldPyfits(pyfits.__version__):
            pyfits.Header.update(self, key, value)
        else:
            pyfits.Header.update(self, {key: value})

    @staticmethod
    def _isOldPyfits(version):
        firstDotPos = version.find(".")
        majorVersion = int(version[0:firstDotPos])

        minorVersionStr = ""
        for each in version[firstDotPos + 1 :]:
            if each.isdigit():
                minorVersionStr += each
            else:
                break

        minorVersion = int(minorVersionStr)

        return majorVersion < 3 or (majorVersion == 3 and minorVersion < 1)

from opticalib.dmutils import iff_module as ifm, iff_processing as ipf # type: ignore
from opticalib import typings as ot# type: ignore


class OttIffAcquisition():
    
    def __init__(self, dm: ot.DeformableMirrorDevice, interf: ot.InterferometerDevice):
        """The constructor"""
        if ot.isinstance_(dm, 'DeformableMirrorDevice'):
            self._dm = dm
        else:
            if all(hasattr(dm, method) for method in ['set_shape', 'get_shape', 'uploadCmdHistory', 'runCmdHistory']):
                self._dm = dm
            else:
                raise TypeError("""
`dm` must be a DeformableMirrorDevice instance, or at least implements the methods:
 - set_shape
 - get_shape
 - uploadCmdHistory
 - runCmdHistory""")
        if ot.isinstance_(interf, 'InterferometerDevice'):
            self._interf = interf
        else:
            if all(hasattr(interf, method) for method in ['acquire_map', 'capture', 'produce']):
                self._interf = interf
            else:
                raise TypeError("""
`interf` must be an InterferometerDevice instance, or at least implements the methods:
 - acquire_map
 - capture
 - produce""")



    
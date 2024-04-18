'''
Authors
  - C. Selmi: written in 2021
'''
from  unittest import TestCase
import numpy as np

from unittest.mock import call, MagicMock, Mock
from numpy import testing
from m4.devices.reference_mirror import OpcUaReferenceMirror
from m4.configuration.ott_parameters import OpcUaParameters


class TestOpcUaReferenceMirror(TestCase):

    def setUp(self):
        self._opc = MagicMock()
        self._rm = OpcUaReferenceMirror(self._opc)

    def testGetPositionReturnSixTuple(self):
        pos = self._rm.getPosition()
        message = "Get position doesn't return a vector of shape 6."
        self.assertEqual(pos.shape, (6,), message)
#         if pos.shape == (6,):
#             test_boolean = True
#         else:
#             test_boolean = False
#         self.assertTrue(test_boolean, message)

    def testCallOpcUaOnGetPosition(self):
        piston = 1e-16
        tip = 1e-6
        tilt = 0.007

        self._values = {OpcUaParameters.RM_PISTON: piston,
                  OpcUaParameters.RM_TIP: tip,
                  OpcUaParameters.RM_TILT: tilt}

        def side_effect(arg):
            return self._values[arg]
        self._opc.get_position = Mock(side_effect=side_effect)

        pos = self._rm.getPosition()
        calls = [call(OpcUaParameters.RM_PISTON),
                 call(OpcUaParameters.RM_TIP),
                 call(OpcUaParameters.RM_TILT)]
        self._opc.get_position.assert_has_calls(calls, any_order=True)
        testing.assert_almost_equal(pos,
                                    np.array([0, 0, piston, tip, tilt, 0]))

    def testCallOpcUaOnSetPosition(self):
        piston = -7
        tip = 30.23
        tilt = -12e-3
        pos = np.array([0, 0, piston, tip, tilt, 0])

        self._values = {}

        def side_effect_set(dof, value):
            self._values[dof] = value

        def side_effect_get(arg):
            return self._values[arg]

        self._opc.get_position = Mock(side_effect=side_effect_get)
        self._opc.set_target_position = Mock(side_effect=side_effect_set)

        newpos = self._rm.setPosition(pos)

        calls = [call(OpcUaParameters.RM_PISTON, piston),
                 call(OpcUaParameters.RM_TILT, tilt),
                 call(OpcUaParameters.RM_TIP, tip)]
        self._opc.set_target_position.assert_has_calls(calls, any_order=True)
        self._opc.move_object.assert_called_once_with(OpcUaParameters.RM_KIN)
        self._opc.wait_for_stop.assert_called_once_with(OpcUaParameters.RM_KIN)

        testing.assert_almost_equal(pos, newpos)

if __name__ == '__main__':
    from m4.utils.roi import TestOpcUaReferenceMirror
    TestOpcUaReferenceMirror().testGetPositionReturnSixTuple()
    TestOpcUaReferenceMirror().testCallOpcUaOnGetPosition()
    TestOpcUaReferenceMirror().testCallOpcUaOnSetPosition()
    
    print('All test passed')

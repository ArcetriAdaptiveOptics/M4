'''
@author: cselmi
'''

import numpy as np
from opcua import Client
from opcua import ua

server = "opc.tcp://192.168.22.100:48050"

class OpcUaController():

    def __init__(self):
        """The constructor """
        self._client = Client(url=server)

    def _test(self):
        self._client.connect()
        var = self._client.get_node("ns=7;s=MAIN.i_Temperature_Sensor[16]")
        value = var.get_value()
        type = var.get_data_type_as_variant_type()
        new_value = ua.DataValue(ua.Variant(int(1), type))
        var.set_value(new_value)


    def stop(self):
        stop_node = self._client.get_node("ns=7;s=MAIN.b_StopCmd")
        stop_type = stop_node.get_data_type_as_variant_type()
        stop_node.set_value(ua.DataValue(ua.Variant(True, stop_type)))

    def get_temperature_vector(self):
        temperature_node = self._client.get_node("ns=7;s=MAIN.i_Temperature_Sensor")
        temperature_vector = np.array(temperature_node.get_value())
        return temperature_vector

    def get_slide(self):
        slide_node = self._client.get_node("ns=7;s=MAIN.Drivers_input.f_PosAct[2]")
        slide = slide_node.get_value()
        return slide

    def get_rslide(self):
        rslide_node = self._client.get_node("ns=7;s=MAIN.Drivers_input.f_PosAct[1]")
        rslide = rslide_node.get_value()
        return rslide

    def get_rotation_angle(self):
        rot_angle_node = self._client.get_node("ns=7;s=MAIN.Drivers_input.f_PosAct[0]")
        rot_angle = rot_angle_node.get_value()
        return rot_angle

    def set_slide(self, par_trans):
        slide_node = self._client.get_node("ns=7;s=MAIN.f_TargetPosition_input[2]")
        type_slide_node = slide_node.get_data_type_as_variant_type()
        slide_node.set_value(ua.DataValue(ua.Variant(par_trans, type_slide_node)))
        slide = slide_node.get_value()
        return slide

    def set_rslide(self, ref_flat):
        rslide_node = self._client.get_node("ns=7;s=MAIN.f_TargetPosition_input[1]")
        type_rslide_node = rslide_node.get_data_type_as_variant_type()
        rslide_node.set_value(ua.DataValue(ua.Variant(ref_flat, type_rslide_node)))
        rslide = rslide_node.get_value()
        return rslide

    def set_rotation_angle(self, rot_ring_angle):
        rot_angle_node = self._client.get_node("ns=7;s=MAIN.f_TargetPosition_input[0]")
        type_rot_angle = rot_angle_node.get_data_type_as_variant_type()
        rot_angle_node.set_value(ua.DataValue(ua.Variant(rot_ring_angle, type_rot_angle)))
        rot_angle = rot_angle_node.get_value()
        return rot_angle

    def move_slide(self):
        move_slide = self._client.get_node("ns=7;s=MAIN.b_MoveCmd[2]")
        move_slide_type = move_slide.get_data_type_as_variant_type()
        move_slide.set_value(ua.DataValue(ua.Variant(True, move_slide_type)))

    def move_rslide(self):
        move_rslide = self._client.get_node("ns=7;s=MAIN.b_MoveCmd[1]")
        move_rslide_type = move_rslide.get_data_type_as_variant_type()
        move_rslide.set_value(ua.DataValue(ua.Variant(True, move_rslide_type)))

    def move_ring_angle(self):
        move_rot = self._client.get_node("ns=7;s=MAIN.b_MoveCmd[0]")
        move_rot_type = move_rot.get_data_type_as_variant_type()
        move_rot.set_value(ua.DataValue(ua.Variant(True, move_rot_type)))

    def _get_command_state(self, int_number):
        self._client.connect()
        node = self._client.get_node("ns=7;s=MAIN.b_MoveCmd[%d]" %int_number)
        value = node.get_value()
        self._client.disconnect()
        return value

    def wait_for_stop(self, int_number):
        value = self._opcUa._get_command_state(0)
        while value == True:
            time.sleep(0.1)
            value = self._opcUa._get_command_state(0)

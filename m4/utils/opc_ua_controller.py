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

    def test(self):
        self._client.connect()
        var = self._client.get_node("ns=7;s=MAIN.i_Temperature_Sensor[16]")
        value = var.get_value()
        type = var.get_data_type_as_variant_type()
        new_value = ua.DataValue(ua.Variant(int(1), type))
        var.set_value(new_value)

    def stop(self):
        stop_node = self._client.get_node("ns=7;s=MAIN.b_StopCmd")
        stop_type = stop_node.get_data_type_as_variant_type()
        stop_node.set_value(True, stop_type)

    def get_temperature_vector(self):
        temperature_node = self._client.get_node("ns=7;s=MAIN.i_Temperature_Sensor")
        temperature_vector = np.array(temperature_node.get_value())
        return temperature_vector

    def get_rotation_angle(self):
        rot_angle_node = self._client.get_node("ns=7;s=MAIN.f_targetPosition_input[0]")
        rot_angle = rot_angle_node.get_value()
        return rot_angle

    def set_rotation_angle(self, rot_ring_angle):
        rot_angle_node = self._client.get_node("ns=7;s=MAIN.f_targetPosition_input[0]")
        type_rot_angle = rot_angle_node.get_data_type_as_variant_type()
        rot_angle_node.set_value(rot_ring_angle, type_rot_angle)
        rot_angle = rot_angle_node.get_value()
        return rot_angle

    def move_ring_angle(self):
        ''' Function to move the rotating ring angle (range: 0 to 180) '''
        move_rot = self._client.get_node("ns=7;s=MAIN.b_MoveCmd[0]")
        move_rot_type = move_rot.get_data_type_as_variant_type()
        move_rot.set_value(True, move_rot_type)


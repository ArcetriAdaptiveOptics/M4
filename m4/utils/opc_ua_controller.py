'''
Autors
  - C. Selmi:  written in September 2020
'''

import logging
import numpy as np
import time
from opcua import Client
from opcua import ua
from m4.configuration.ott_parameters import OttParameters, OpcUaParameters

server = "opc.tcp://192.168.22.100:48050"

class OpcUaController():
    """
    Function for test tower management via OpcUa::
    
        from m4.utils.opc_ua_controller import OpcUaController
        opc = OpcUaController()
    """
    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('OPCUA:')
        self._client = Client(url=server)

    def stop(self):
        """
        Stop all commands
        """
        self._client.connect()
        stop_node = self._client.get_node("ns=7;s=MAIN.b_StopCmd")
        stop_type = stop_node.get_data_type_as_variant_type()
        stop_node.set_value(ua.DataValue(ua.Variant(True, stop_type)))
        self._client.disconnect()

    def get_temperature_vector(self):
        """
        Returns
        -------
            temperature_vector: numpy array
                                values obtained from PT
        """
        self._client.connect()
        temperature_node = self._client.get_node("ns=7;s=MAIN.i_Temperature_Sensor")
        temperature_vector = np.array(temperature_node.get_value())/100.
        self._client.disconnect()
        return temperature_vector

    def get_variables_positions(self):
        """
        Returns
        -------
            variables: numpy array
                    all variables value
        """
        self._client.connect()
        var_list = []
        for i in range(len(OpcUaParameters.zabbix_variables_name)):
            node = self._client.get_node("ns=7;s=MAIN.Drivers_input.f_PosAct[%d]" %i)
            var = node.get_value()
            var_list.append(var)
        self._client.disconnect()
        return np.array(var_list)


### Command for object ###
    def get_position(self, int_number):
        """
        Parameters
        ----------
            int_number: int
                    number of the chosen object

        Returns
        -------
            position: float
                    position of the requested object
        """
        self._client.connect()
        node = self._client.get_node("ns=7;s=MAIN.Drivers_input.f_PosAct[%d]" %int_number)
        position = node.get_value()
        self._client.disconnect()
        self._logger.debug('Position = %f', position)
        return position

    def set_target_position(self, int_number, value):
        """
        Parameters
        ----------
            int_number: int
                    number of the chosen object
            value: float
                value to assign to the chosen object

        Returns
        -------
            target_position: float
                    value assigned to the chosen object
                    (not applied)
        """
        self._client.connect()
        node = self._client.get_node("ns=7;s=MAIN.f_TargetPosition_input[%d]" %int_number)
        type_node = node.get_data_type_as_variant_type()
        node.set_value(ua.DataValue(ua.Variant(value, type_node)))
        target_position = node.get_value()
        self._client.disconnect()
        self._logger.debug('Target position = %f', target_position)
        return target_position

    def move_object(self, int_number):
        """
        Function that applies command and moves object

        Parameters
        ----------
            int_number: int
                    number of the chosen object
        """
        self._client.connect()
        node = self._client.get_node("ns=7;s=MAIN.b_MoveCmd[%d]" %int_number)
        node_type = node.get_data_type_as_variant_type()
        node.set_value(ua.DataValue(ua.Variant(True, node_type)))
        self._client.disconnect()
        self._logger.debug('Object moved successfully')

    def _get_command_state(self, int_number):
        """
        Parameters
        ----------
            int_number: int
                    number of the chosen object

        Returns
        -------
            value: boolean 
                    position of the requested object
        """
        self._client.connect()
        node = self._client.get_node("ns=7;s=MAIN.b_MoveCmd[%d]" %int_number)
        value = node.get_value()
        self._client.disconnect()
        return value

    def wait_for_stop(self, int_number):
        """
        Function to wait for the movement to be completed
        
        Parameters
        ----------
            int_number: int
                    number of the chosen object
        """
        value = self._get_command_state(int_number)
        while value == True:
            time.sleep(0.1)
            value = self._get_command_state(int_number)

###
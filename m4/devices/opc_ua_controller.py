'''
Authors
  - C. Selmi:  written in September 2020
'''

import logging
import numpy as np
import time
from opcua import Client
from opcua import ua
from m4.configuration.ott_parameters import OpcUaParameters

server = OpcUaParameters.server


class OpcUaController():
    """
    Function for test tower management via OpcUa

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
    """
    STOP_NODE = "ns=7;s=MAIN.b_StopCmd"
    TEMPERATURE_NODE = "ns=7;s=MAIN.i_Temperature_Sensor"

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('OPCUA:')
        self._client = Client(url=server)

    def stop(self):
        """
        Stop all commands
        """
        self._client.connect()
        stop_node = self._client.get_node(self.STOP_NODE)
        stop_type = stop_node.get_data_type_as_variant_type()
#         import pdb
#         pdb.set_trace()
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
        temperature_node = self._client.get_node(self.TEMPERATURE_NODE)
        temperature_vector = np.array(temperature_node.get_value()) / 100.
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
            node = self._client.get_node("ns=7;s=MAIN.Drivers_input.f_PosAct[%d]" % i)
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
        node = self._client.get_node("ns=7;s=MAIN.Drivers_input.f_PosAct[%d]" % int_number)
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
        node = self._client.get_node("ns=7;s=MAIN.f_TargetPosition_input[%d]" % int_number)
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
        node = self._client.get_node("ns=7;s=MAIN.b_MoveCmd[%d]" % int_number)
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
        node = self._client.get_node("ns=7;s=MAIN.b_MoveCmd[%d]" % int_number)
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

    def readActsPositions(self, n1, n2, n3):
        '''
        Function to read actuators positions

        Parameters
        ----------
            n1, n2, n3: int
                    number of the chosen object

        Returns
        -------
            acts: numpy array
                vector of actuators position
        '''
        act1 = self.get_position(n1)
        act2 = self.get_position(n2)
        act3 = self.get_position(n3)
        return np.array([act1, act2, act3])

    def setActsPositions(self, n1, n2, n3, v1, v2, v3):
        '''
        Function to set actuators positions

        Parameters
        ----------
            n1, n2, n3: int
                    number of the chosen object
            v1, v2, v3: int, float
                    value to pass to actuators

        Returns
        -------
            acts: numpy array
                vector of actuators position
        '''
        act1 = self._setAct(n1, v1)
        act2 = self._setAct(n2, v2)
        act3 = self._setAct(n3, v3)
        return np.array([act1, act2, act3])

    def _setAct(self, number, value):
        ''' specific function for actuators because on these
        does not work the wait for stop (not set the transition
        from true to false by ads)'''
        self.set_target_position(number, value)
        self.move_object(number)
        time.sleep(10)
        # self.wait_for_stop(number)
        act = self.get_position(number)
        return act


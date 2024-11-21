"""
Author(s)
---------
    - Pietro Ferraiuolo: written in 2024
    - Runa Briguglio: written in 2024
"""

import zmq
import time
import struct
import numpy as np
from m4.ground import logger_set_up as logger
from m4.configuration.ott_parameters import OttParameters
from m4.devices.base_m4_exapode import BaseM4Exapode

# The kinematics is written in the form (normalization-less):

#           C     ^Y
#                 |
#        A     B  -->X
#   Pist  Tx  Ty
# A    1  -1  -1
# B    1   1  -1
# C    1   0   1


class DpMotors(BaseM4Exapode):
    """Class for M4 control via opc ua??

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.m4_exapode import OpcUaM4Exapode
    """

    DPinterface = 0  # dummy, to replicate the sintax

    def __init__(self, DPinterface):
        """The constructor"""
        self._dpinterface = DPinterface
        self._kinematrix = np.array([[1, 1, 1], [-1, 1, 0], [-1, -1, 1]])
        self._invkinematrix = np.linalg.inv(self._kinematrix)
        self._context = None
        self._socket = None
        self.remote_ip = "192.168.22.51"  # final?
        self.remote_port = 6660  # final?
        self._m4dof = slice(OttParameters.M4_DOF[0], OttParameters.M4_DOF[1] + 1)
        self._minVel = 30 # um/s
        self._timeout = 180 # seconds
        # logger.set_up_logger('path/dpMotors.log', 10)

    def getPosition(self):
        """
        Function to get the current position of the actuators.

        Returns
        -------
        pos : list or ArrayLike
            The current position of the actuators.
        """
        pp = self._getMotorPosition()
        kk = self._motor2kinematics(pp)
        pos = np.array([0, 0, 0] + [k for k in kk[1:]] + [0])
        print(pos)
        return pos

    def setPosition(self, absolute_position_in_mm: list):
        """
        Function to set the position of the DP actuators.

        Parameters
        ----------
        absolute_position_in_mm : list or ArrayLike
            The 6-dimentional array of the absolute positions, in mm, of all
            the degrees of freedom of the DP. Since only Tip/Til can be commanded,
            the array should have the form [0,0,0,x,y,0].

        Returns
        -------
        pos : list or ArrayLike
            The current position of the actuators.
        """
        v = np.array([0] + absolute_position_in_mm[self._m4dof])
        vm = self._kinematics2motor(v)
        self._setMotorPosition(vm)
        return self.getPosition()

    def _getMotorPosition(self):
        """
        Middle level function to get the current position of the actuators.

        Returns
        -------
        positions : list
            A list containing, in order, the encoder position of the actuators
            read from the bus box, in millimeters.
        """
        connected = self._connectBusBox()
        if connected is not True:
            self._recconnectBusBox()
        reading = self._send_read_message()
        positions = self._extract_motor_position(reading)
        positions =  [p*1e-3 for p in positions]  # conversion in mm
        self._disconnectBusBox()
        return positions

    def _setMotorPosition(self, motorcmd):
        """
        Middle level function to set the position of the actuators.

        Parameters
        ----------
        motorcmd : list
            A list containing, in order, the absolute position in mm of the
            actuators to be set.

        Returns
        -------
        response : dict
            The response from the BusBox.
        """
        pos = motorcmd * 1e3  # conversion in um
        print(pos)
        connected = self._connectBusBox()
        if connected is not True:
            self._recconnectBusBox()
        self._send_message(pos)
        self._disconnectBusBox()
        return

    def _motor2kinematics(self, pos):
        """
        Function to convert the absolute positions of the actuators from motor
        to zernike positions.

        Parameters
        ----------
        pos: list or ArrayLike
            The 3-component only vector containing the motor encoder positions.

        Returns
        -------
        kin: list
            The 3-component only vector containing the kinematic positions of the
            Motors.
        """
        kin = np.dot(pos, self._kinematrix)
        return kin

    def _kinematics2motor(self, pos):
        """
        Function to convert the absolute positions of the actuators from zernike
        to motor positions.

        Parameters
        ----------
        pos: list or ArrayLike
            The 3-component only vector containing the kinematic positions of the
            Motors.

        Returns
        -------
        motor_pos: list
            The 3-component only vector containing the motor encoder positions.
        """
        motor_pos = np.dot(pos, self._invkinematrix)
        return motor_pos

    def _extract_motor_position(self, response):
        """
        Function to extract the motor position from the response of the BusBox.

        Parameters
        ----------
        response : dict
            The response from the BusBox.

        Returns
        -------
        positions : list
            A list containing, in order, the encoder position of the actuators
            read from the bus box.
        """
        position = [response["actuators"][act]["encoder_position"] for act in range(3)]
        return position

    def _send_message(self, acts_pos_in_um: list):
        """
        Function wich comunicates with the BusBox for moving the actuators.

        Parameters
        ----------
        acts_pos_in_um: int or list
            A List containing the absolute positions, in um, to actuate of the
            DP actuators.

        Returns
        -------
        response : bytearray
            The 116 bytes response from the BusBox, to be decoded.
        """
        act_positions = np.array([pos for pos in acts_pos_in_um], dtype=int)
        cmd = struct.pack(
            "<BBBhhh", 17, 7, 7, act_positions[0], act_positions[1], act_positions[2]
        )
        cmd += bytearray(10)
        out = self._send(cmd)
        time.sleep(3)
        response = self._decode_message(out)
        print(f"{response['actuators'][0]['target_position']}")
        check = self._check_actuation_success(act_positions, cmd, response)
        if check is False:
            print("Actuation failed")
            return
        else:
            return self._send_read_message()

    def _check_actuation_success(self, target_pos, cmd, response):
        """
        Function to check if the actuation was successful.

        Parameters
        ----------
        trget_pos: list of int
            The target position to reach of the actuators.
        cmd : bytearray
            The message sent to the BusBox.
        response : dict
            The response from the BusBox.

        Returns
        -------
        success : bool
            True if the actuation was successful, False otherwise.
        """
        act_pos = [b['actual_position'] for b in [c for c in response['actuators']]][:3]
        act_enc_pos = [b['encoder_position'] for b in [c for c in response['actuators']]][:3]
        tot_time = 0
        pos_err = np.max(np.abs(target_pos - act_enc_pos))
        timeout = self._timeout
        while pos_err > 1:
            waittime = np.min([pos_err/self._minVel, 3])
            print(waittime)
            tot_time += waittime
            self._timeout_check(tot_time, timeout)
            out = self._send(cmd)
            time.sleep(waittime)
            response = self._decode_message(out)
            act_pos = [b['actual_position'] for b in [c for c in response['actuators']]][:3]
            act_enc_pos = [b['encoder_position'] for b in [c for c in response['actuators']]][:3] 
            pos_err = np.max(np.abs(target_pos - act_enc_pos))
            print(f"ActPos: {act_pos} ; EncPos: {act_enc_pos}\n")
        return True

    def _send_read_message(self):
        """
        Function wich comunicates with the BusBox for reading the actuators
        position, done through a null byte command.

        Returns
        -------
        response : bytearray
            The 116 bytes response from the BusBox, to be decoded.
        """
        read_cmd = bytearray(hex(17), "utf-8") + bytearray(18)
        out = self._send(read_cmd)
        response = self._decode_message(out)
        return response

    def _send(self, bytestream):
        """
        Function wich comunicates with the BusBox.

        Parameters
        ----------
        bytestream : bytearray
            The message to be sent to the BusBox.
        """
        try:
            self._socket.send(bytestream)
            out = self._socket.recv()
            return out
        except zmq.ZMQError as ze:
            print(ze)
            return False

    def _decode_message(self, message):
        """
        Function which decodes the message from the BusBox into human readable format.

        Parameters
        ----------
        message : bytearray
            The 116 bytes response from the BusBox.

        Returns
        -------
        response : dict
            The decoded message.
        """
        act_struct = "hhib3bH"  # short*2, int, char, 3-char array, unsigned short
        struct_size = struct.calcsize("Bbbb") # 1 unsigned char, 3 char
        outer_data = struct.unpack("Bbbb", message[:struct_size])
        heart_beat_counter, voltage, temperature, status = outer_data
        actuators = []
        act_size = struct.calcsize(act_struct)
        # Unpack each inner struct
        for i in range(3):
            start = struct_size + i * act_size
            end = start + act_size
            status_bytes = np.zeros(3, dtype=int)
            (
                actual_position,
                target_position,
                encoder_position,
                actual_velocity,
                *status_bytes,
                bus_voltage,
            ) = struct.unpack(act_struct, message[start:end])
            status_hex = status_bytes # bytearray([b for b in status_bytes])
            actuators.append(
                {
                    "actual_position": actual_position,
                    "target_position": target_position,
                    "encoder_position": encoder_position,
                    "actual_velocity": actual_velocity,
                    "bus_voltage": bus_voltage,
                    "status": status_hex,
                }
            )
        decodified_out_message = {
            "heart_beat_counter": heart_beat_counter,
            "voltage": voltage,
            "temperature": temperature,
            "status": status,
            "actuators": actuators,
        }
        return decodified_out_message

    def _timeout_check(self, total_time, timeout):
        """
        Function to check if the timeout is reached.

        Parameters
        ----------
        total_time : float
            The total time spent in the loop.
        timeout : float
            The maximum time allowed for the loop.

        Returns
        -------
        timeout_reached : bool
            True if the timeout is reached, False otherwise.
        """
        if total_time >= timeout:
            print("Timeout reached, probably BusBox freezed. Trying reconnecting...")
            reconnected = self._recconnectBusBox()
            if reconnected:
                print("Reconnection successful, re-applying command")
                pass
            else:
                raise ConnectionError("Reconnection failed, exiting...")
        if total_time > 2*timeout:
            raise ConnectionError("Something's wrong. Exiting...")
        

    def _connectBusBox(self):
        """
        Connection to the BusBox controlling the DP motors actuators, via
        ZMQ Pair protocol.
        """
        try:
            self._context = zmq.Context()
            self._socket = self._context.socket(zmq.PAIR)
            self._socket.connect(f"tcp://{self.remote_ip}:{self.remote_port}")
            return True
        except Exception as e:
            # logger.log(f"ConnectZMQ: {e}")
            return False
        
    def _recconnectBusBox(self):
        """
        Reconnection to the BusBox controlling the DP motors actuators.
        """
        self._disconnectBusBox()
        connected = self._connectBusBox()
        if connected is not True:
            for ntry in range(5):
                print(f"Connection failed, retrying... {ntry}")
                connected = self._connectBusBox()
                if connected is True:
                    break
                time.sleep(1)
        if connected is False:
            print("Connection failed")
            return False
        return True

    def _disconnectBusBox(self):
        """
        Disconnection from the BusBox controlling the DP motors actuators.
        """
        try:
            self._socket.disconnect(f"tcp://{self.remote_ip}:{self.remote_port}")
            self._socket.close()
            self._context.term()
            return True
        except Exception as e:
            # logger.log(f"DisconnectZMQ: {e}")
            return False


# _______________________
# Dictionary
# TCS Interface KAMAL (C#)
#
# struct status_actuator {
#    short actual_position;         2by - h
#    short target_position;         2by - h
#    int encoder_position;          4by - i
#    char actual_velocity;          1by - b
#    unsigned short bus_voltage;    2by - H
#  //unsigned short run_current;
#    char status[3];                1by x 3 - 3b
# }__attribute((packed)); -> 14 bytes - hhibH3b

# struct status_packet
#    {
#    unsigned char heart_beat_counter;      1by - B
#    char voltage;                          1by - b
#    char temperature;                      1by - b
#    char status;                           1by - b
#    struct status_actuator actuators[8];   14by x 8 -> 112
# }__attribute((packed)); -> 116by - (Bbbb + hhibH3b * 8)

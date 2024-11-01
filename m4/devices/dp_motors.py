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
    ''' Class for M4 control via opc ua??

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.m4_exapode import OpcUaM4Exapode
    '''
    DPinterface = 0 #dummy, to replicate the sintax 
    def __init__(self, DPinterface):
        """The constructor """
        self._dpinterface = DPinterface
        self._kinematrix = np.array([[1,1,1],[-1,1,0],[-1,-1,1]])
        self._invkinematrix = np.linalg.inv(self._kinematrix)
        self._context = None
        self._socket = None
        self.remote_ip = "192.168.22.51" # final?
        self.remote_port = 6660 # final?
        logger.set_up_logger('path/dpMotors.log', 10)

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
        pos = np.zeros(6)
        pos[OttParameters.M4_DOF] = kk[1:]
        print(pos)
        return pos

    def setPosition(self, absolute_position_in_mm):
        ''' to be implemented
        '''
        v = np.array([0,absolute_position_in_mm[OttParameters.M4_DOF]])
        vm = np.dot(v, self._kinematrix)
        print(vm)
        self._setMotorPosition(vm)
        pass

    def _getMotorPosition(self):
        """
        Function to get the current position of the actuators.
        
        Returns
        -------
        positions : list
            A list containing, in order, the encoder position of the actuators
            read from the bus box, in millimeters.
        """
        self._connectBusBox()
        reading = self._decode_message(self._send_read_message())
        positions = self._extract_motor_position(reading)
        positions = positions*10e-3 # conversion in mm
        self._disconnectBusBox()
        return positions

    def _setMotorPosition(self, motorcmd):
        ''' to be implemented, at low level
        '''
        pos = motorcmd*10e3 # conversion in um
        self._connectBusBox()
        # Si manda il comando
        # si riceve risposta 
        # si fanno le conversioni del caso
        self._disconnectBusBox()
        pass

    def _motor2kinematics(self, pos):
        '''
        '''
        ptt = np.dot(pos, self._kinematrix)
        #return ptt
        pass

    def _kinematics2motor(self,pos):
        '''
        '''
        abc = np.dot(pos, self._invkinematrix)
        #return abc
        pass

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
        position = [response['actuators'][act]['encoder_position'] for act in range(3)]
        return position
    
    def _send_message(self, abs_pos_in_um, act_n:int=3):
        """
        Function wich comunicates with the BusBox for moving the actuators.

        Parameters
        ----------
        abs_pos: int or list
            The absolute position in um.
        act_n: int
            The actuator to move. default is 3, which means all actuators are 
            moved, in which case abs_pos should be a list or array, whose order
            will be the same as the actuators.

        Returns
        -------
        response : bytearray
            The 116 bytes response from the BusBox, to be decoded.
        """
        act_positions = np.zeros(3, dtype=int)
        act_positions[act_n] = abs_pos_in_um
        cmd = struct.pack('<BBBhhh', 
                        17,
                        7,
                        7, 
                        act_positions[0],
                        act_positions[1],
                        act_positions[2])
        cmd += bytearray(12)
        out = self._send(cmd)
        response = self._decode_message(out)
        act_pos = response['actuators'][act_n]['actual_position']
        act_enc_pos = response['actuators'][act_n]['encoder_position']
        act_targ_pos = abs_pos_in_um
        tot_time = 0
        while np.abs(act_targ_pos-act_enc_pos) > 4:
            waittime = 3
            tot_time += waittime
            print("waiting ", waittime)
            time.sleep(waittime)
            out = self._send(cmd)
            response = self._decode_message(out)
            act_pos = response['actuators'][act_n]['actual_position']
            act_enc_pos = response['actuators'][act_n]['encoder_position']
            print(f"actual position: {act_pos}\nencoder_position: {act_enc_pos}")
        check = self._send_read_message()
        response = self._decode_message(check)
        response['execution_time'] = tot_time
        return response

    def _send_read_message(self):
        """
        Function wich comunicates with the BusBox for reading the actuators
        position, done through a null byte command.

        Returns
        -------
        response : bytearray
            The 116 bytes response from the BusBox, to be decoded.
        """
        read_cmd = bytearray(hex(17), 'utf-8') + bytearray(18)
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
        act_struct = 'hhibH3b' # short*2, int, char, unsigned short, 3-char array
        struct_size = struct.calcsize('Bbbb')
        outer_data = struct.unpack('Bbbb', message[:struct_size])
        heart_beat_counter, voltage, temperature, status = outer_data
        actuators = []
        act_size = struct.calcsize(act_struct)
        # Unpack each inner struct
        for i in range(3):
            start = struct_size + i * act_size
            end = start + act_size
            status_bytes = np.zeros(3, dtype=int)
            actual_position, target_position, encoder_position, actual_velocity, bus_voltage, *status_bytes = struct.unpack(act_struct, message[start:end])
            # Decode the status field from bytes to string and strip null bytes
            status_str = ''.join(chr(b if b >= 0 else 256 + b) for b in status_bytes if b != 0)
            actuators.append({
                'actual_position': actual_position,
                'target_position': target_position,
                'encoder_position': encoder_position,
                'actual_velocity': actual_velocity,
                'bus_voltage': bus_voltage,
                'status': status_str
            })
        decodified_out_message = {
            'heart_beat_counter': heart_beat_counter,
            'voltage': voltage,
            'temperature': temperature,
            'status': status,
            'actuators': actuators
        }
        return decodified_out_message

    def _connectBusBox(self):
        try:
            self._context = zmq.Context()
            self._socket = self._context.socket(zmq.PAIR)
            self._socket.connect(f"tcp://{self.remote_ip}:{self.remote_port}")
            return True
        except Exception as e:
            logger.log(f"ConnectZMQ: {e}")
            return False

    def _disconnectBusBox(self):
        try:
            self._socket.disconnect(f"tcp://{self.remote_ip}:{self.remote_port}")
            self._socket.close()
            self._context.term()
            return True
        except Exception as e:
            logger.log(f"DisconnectZMQ: {e}")
            return False
        pass


#_______________________
# Dictionary

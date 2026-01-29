"""
Authors
  - C. Selmi:  written in September 2020
"""

import logging
import numpy as np
import time
from opcua import Client
from opcua import ua
from m4.configuration.ott_parameters import OpcUaParameters as _opcpar


##################### CERTIFICATE CHECKING AND GENERATION #####################

def generate_certificate(cert_path, key_path):
    from cryptography import x509
    from cryptography.x509.oid import NameOID
    from datetime import datetime, timedelta, timezone
    from cryptography.hazmat.primitives import hashes, serialization
    from cryptography.hazmat.primitives.asymmetric import rsa

    print("Generating client certificate and key...")
    private_key = rsa.generate_private_key(public_exponent=65537, key_size=2048)

    subject = issuer = x509.Name(
        [
            x509.NameAttribute(NameOID.COUNTRY_NAME, "EN"),
            x509.NameAttribute(NameOID.ORGANIZATION_NAME, "PythonOpcUaClient"),
        ]
    )

    cert = (
        x509.CertificateBuilder()
        .subject_name(subject)
        .issuer_name(issuer)
        .public_key(private_key.public_key())
        .serial_number(x509.random_serial_number())
        .not_valid_before(datetime.now(timezone.utc) - timedelta(days=1))
        .not_valid_after(datetime.now(timezone.utc) + timedelta(days=3650))
        .sign(private_key, hashes.SHA256())
    )

    with open(key_path, "wb") as f:
        f.write(
            private_key.private_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PrivateFormat.TraditionalOpenSSL,
                encryption_algorithm=serialization.NoEncryption(),
            )
        )

    with open(cert_path, "wb") as f:
        f.write(cert.public_bytes(serialization.Encoding.DER))

    print("Certificate and key created.")

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CLIENT_CERT_DER = os.path.join(SCRIPT_DIR, "client_cert.der")
CLIENT_KEY_PEM = os.path.join(SCRIPT_DIR, "client_key.pem")
SERVER_CERT_DER = os.path.join(SCRIPT_DIR, "server_cert.der")

if not (os.path.exists(CLIENT_CERT_DER) and os.path.exists(CLIENT_KEY_PEM)):
    generate_certificate(CLIENT_CERT_DER, CLIENT_KEY_PEM)

if not os.path.exists(SERVER_CERT_DER):
    print(f"Error: Server certificate not found at '{SERVER_CERT_DER}'")
    sys.exit(1)

del (
    SCRIPT_DIR,
    generate_certificate,
    os,
    sys
)

##################### CERTIFICATE CHECKING AND GENERATION #####################


class OpcUaController:
    """
    Function for test tower management via OpcUa

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
    """

    STOP_NODE = "ns=7;s=MAIN.b_StopCmd"
    TEMPERATURE_NODE = "ns=7;s=MAIN.i_Temperature_Sensor"

    def __init__(self):
        """The constructor"""
        self._logger = logging.getLogger("OPCUA:")
        self._client = Client(url=_opcpar.server)
        self._client.set_user(_opcpar.username)
        self._client.set_password(_opcpar.password)
        self._client.set_security_string(
            "Basic256Sha256,SignAndEncrypt,%s,%s,%s"
            % (_opcpar.CLIENT_CERT_DER, _opcpar.CLIENT_KEY_PEM, _opcpar.SERVER_CERT_DER)
        )

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
        temperature_vector = np.array(temperature_node.get_value()) / 100.0
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
        for i in range(len(_opcpar.zabbix_variables_name)):
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
        node = self._client.get_node(
            "ns=7;s=MAIN.Drivers_input.f_PosAct[%d]" % int_number
        )
        position = node.get_value()
        self._client.disconnect()
        self._logger.debug("Position = %f", position)
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
        node = self._client.get_node(
            "ns=7;s=MAIN.f_TargetPosition_input[%d]" % int_number
        )
        type_node = node.get_data_type_as_variant_type()
        node.set_value(ua.DataValue(ua.Variant(value, type_node)))
        target_position = node.get_value()
        self._client.disconnect()
        self._logger.debug("Target position = %f", target_position)
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
        self._logger.debug("Object moved successfully")

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
        """
        Function to read actuators positions

        Parameters
        ----------
            n1, n2, n3: int
                    number of the chosen object

        Returns
        -------
            acts: numpy array
                vector of actuators position
        """
        act1 = self.get_position(n1)
        act2 = self.get_position(n2)
        act3 = self.get_position(n3)
        return np.array([act1, act2, act3])

    def setActsPositions(self, n1, n2, n3, v1, v2, v3):
        """
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
        """
        act1 = self._setAct(n1, v1)
        act2 = self._setAct(n2, v2)
        act3 = self._setAct(n3, v3)
        return np.array([act1, act2, act3])

    def _setAct(self, number, value):
        """specific function for actuators because on these
        does not work the wait for stop (not set the transition
        from true to false by ads)"""
        self.set_target_position(number, value)
        self.move_object(number)
        time.sleep(10)
        # self.wait_for_stop(number)
        act = self.get_position(number)
        return act


    def _check_plc_is_running(self):
        try:
            st = self._client.get_node(_opcpar.NODE_ID_PLC_STATE).get_value()
            if st == _opcpar.PLC_RUN_STATE_VALUE:
                return True
            print(f"PLC is not RUN. State={st}")
            return False
        except Exception as e:
            print(f"PLC state read error (ignored): {e}")
            return True


    def _read_telemetry(self):
        """
        Access the PLC telemetry data structure, for reading.
        
        Returns
        -------
            tel: telemetry object
        """
        try:
            node = self._client.get_node(_opcpar.NODE_ID_TELEMETRY)
            tel = node.get_value()
            
            ## FOR DEBUGGING PURPOSES ##
            
            def pretty_print_obj(obj, indent=0, max_list=30):
                """Helper function to print telemetry readings"""
                sp = " " * indent

                if isinstance(obj, list):
                    for i, item in enumerate(obj[:max_list]):
                        print(f"{sp}[{i}]:")
                        pretty_print_obj(item, indent + 4, max_list=max_list)
                    if len(obj) > max_list:
                        print(f"{sp}... ({len(obj) - max_list} more)")
                    return

                if not hasattr(obj, "__dict__"):
                    print(f"{sp}{obj}")
                    return

                for k, v in obj.__dict__.items():
                    if k.startswith("_"):
                        continue
                    print(f"{sp}{k}:")
                    pretty_print_obj(v, indent + 4, max_list=max_list)
                    
            print("\n--- Telemetry ---")
            pretty_print_obj(tel)

            ###############################

            return tel
        except Exception as e:
            print(f"Telemetry read error: {e}")
            return None
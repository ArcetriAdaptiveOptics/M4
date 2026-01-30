"""
Authors
-------
Chiara Selmi:  written in September 2020
Pietro Ferraiuolo: Modified in 2026 w/ ADS-International
"""

import numpy as np
from contextlib import contextmanager
from opcua import Client
import opcua as ua
from opticalib.ground.logger import SystemLogger
from m4.configuration.ott_parameters import OpcUaParameters as _opcpar

## Necessary patch to fix python-opcua bug ##
try:
    from m4.configuration import opcua_patch
except Exception:
    pass


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

del (SCRIPT_DIR, generate_certificate, os, sys)

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
        self._logger = SystemLogger(__class__)
        self._client = None
        self._connect_client()
        self._cmd_keywords = {
            "MOVE": ["FULL_MOVE"],
            "HOME": ["HOME"],
            "STOP": ["STOP"],
        }

    # LEGACY COMMAND TO (maybe) IMPLEMENT
    # stop (all devices)
    # get var positions
    # get position
    # set target position   For CAR, RA, STW
    # move object           For CAR, RA, STW
    # wait for stop
    # read acts positions   For PAR and RM
    # set acts positions    For PAR and RM

    def get_temperature_vector(self):
        """
        Returns
        -------
            temperature_vector: numpy array
                                values obtained from PT
        """
        with self.connected():
            telemetry = self._read_telemetry()
            try:
                if hasattr(telemetry.System, "Temperatures"):
                    tempvec = telemetry.System.Temperatures
                    return np.array(tempvec)
                else:
                    self._logger.error(
                        "Telemetry has no 'System.Temperatures' attribute"
                    )
                    return np.array([])
            except Exception as e:
                self._logger.error(f"Error accessing temperatures from telemetry: {e}")
                print(f"Error accessing temperatures from telemetry: {e}")

    def send_command(self, mode: str, device: str, target_values: float | list[float]):
        """ """
        plc_is_running = self._check_plc_is_running()
        if not plc_is_running is True:
            raise RuntimeError("PLC is not in RUN state.") from plc_is_running
        
        ## Getting the mode
        mode = mode.upper().strip()
        if not any([mode == m for m in ["MOVE", "HOME", "STOP"]]):
            raise ValueError(
                f"Invalid mode: {mode}. Must be one of 'MOVE', 'HOME', 'STOP'."
            )

        component = self.component_selector(device)
        is_axis = component.value in (1, 2, 3)  # RA/CAR/STW
        is_tripod = component.value in (4, 5)  # RM/PAR

        if not any([is_axis, is_tripod]):
            raise ValueError(
                f"Component {component.name}={component.value} not supported for MOVE/HOME/STOP by your rules."
            )

        node, cmd = self._get_cmd_struct()
        AxisCmdEnum = type(cmd.axisCmd)
        TripodCmdEnum = type(cmd.tripodCmd)
        AxisIdEnum = type(cmd.axisID)
        SystemEnum = type(cmd.systemCmd) if hasattr(cmd, "systemCmd") else None

        # 1) Apply component
        cmd.componentID = component

        # 2) Force axisID = 0 (hidden)
        cmd.axisID = self.__zero_or_first(AxisIdEnum)

        # 3) Clear numeric fields (avoid stale values)
        cmd.targetPos = 0.0
        cmd.tip = 0.0
        cmd.tilt = 0.0
        cmd.piston = 0.0

        # 4) Clear enums to neutral (0 if exists)
        if SystemEnum is not None:
            cmd.systemCmd = self.__zero_or_first(SystemEnum)
        cmd.axisCmd = self.__zero_or_first(AxisCmdEnum)
        cmd.tripodCmd = self.__zero_or_first(TripodCmdEnum)

        # 5) Set proper command enum + ask needed params
        match mode:

            case "MOVE":

                if is_axis:
                    cmd.axisCmd = self._choose_enum_by_mode(AxisCmdEnum, mode)
                    cmd.targetPos = float(target_values)

                elif is_tripod:
                    cmd.tripodCmd = self._choose_enum_by_mode(TripodCmdEnum, mode)
                    if not isinstance(target_values, list) or len(target_values) != 3:
                        raise ValueError(
                            "For TRIPOD MOVE, target_values must be a list of three floats: [tip, tilt, piston]."
                        )
                    cmd.tip = float(target_values[0])
                    cmd.tilt = float(target_values[1])
                    cmd.piston = float(target_values[2])

            case "STOP":
                pass

            case "HOME":
                pass

            case _:
                raise ValueError(
                    f"Invalid mode: {mode}. Must be one of 'MOVE', 'HOME', 'STOP'."
                )
        # 6) Write back
        self._write_struct(node, cmd)


    def component_selector(self, device: str) -> object:
        """
        Select the active component in the PLC.

        Parameters
        ----------
        device : str
            The component to select.

        Returns
        -------
        component: object
            The selected component ID, to be used to write the command to send.
        """
        command_struct = self._get_cmd_struct()
        members = list(type(command_struct.componentID))
        component = next((f for f in members if f.name == device), None)
        if component is None:
            raise ValueError(f"Invalid component ID: {component}")
        return component


    @contextmanager
    def connected(self):
        """
        Context manager to connect/disconnect the OPC UA client.

        Usage
        -----
        with opcua_controller.connected():
            # perform operations with self._client
        """
        try:
            plc_is_running = self._check_plc_is_running()
            if plc_is_running is True:
                self._client.connect()
                self._logger.info("Connected to OPC UA server.")
                self._client.load_type_definitions()
                self._logger.info("Type definitions loaded.")
            else:
                raise (plc_is_running)
            yield
        finally:
            self._client.disconnect()

    def _choose_enum_by_mode(self, enum_cls: object, mode: str):
        """
        mode: MOVE/HOME/STOP
        keywords_map example:
            {"MOVE":["FULL_MOVE"], "HOME":["HOME"], "STOP":["STOP"]}
        """
        mode = mode.upper().strip()
        mem = self.__enum_members(enum_cls)
        if not mem:
            raise RuntimeError(f"Enum {enum_cls} has no members")

        keywords = self._cmd_keywords.get(mode, [mode])

        matches = []
        for m in mem:
            name = getattr(m, "name", str(m)).upper()
            if any(k.upper() in name for k in keywords):
                matches.append(m)

        if len(matches) == 1:
            return matches[0]
        else:
            raise ValueError(
                f"Could not uniquely identify command for mode '{mode}'. Matches found: {[m.name for m in matches]}"
            )

    def _get_cmd_struct(self):
        """
        Access the PLC command data structure, for reading/writing.

        Returns
        -------
        command_struct: command object

        Raises
        -------
        RuntimeError
            If the command struct is NULL or has no 'componentID' field.
        """
        with self.connected():
            try:
                command_node = self._client.get_node(_opcpar.NODE_ID_COMMAND_VARIABLE)
                command_struct = command_node.get_value()

                if command_struct is None:
                    raise RuntimeError(
                        "Command struct is NULL on server (PLC not initialized)."
                    )

                # Required Fields for CMD struct
                required = [
                    "axisCmd",
                    "tripodCmd",
                    "componentID",
                    "axisID",
                    "targetPos",
                    "tip",
                    "tilt",
                    "piston",
                ]
                missing = [f for f in required if not hasattr(command_struct, f)]
                if missing:
                    self._logger.error(
                        "Decoded ST_Command is missing fields: %s", missing
                    )
                    raise RuntimeError("Decoded ST_Command is missing fields:", missing)

                return command_node, command_struct
            except Exception as e:
                self._logger.error(f"Command struct read error: {e}")
                raise e

    def _check_plc_is_running(self) -> bool | Exception:
        """
        Check if the PLC is in RUN state.

        Returns
        -------
            is_running: bool
                True if PLC is running, False otherwise
        """
        try:
            self._client.connect()
            st = self._client.get_node(_opcpar.NODE_ID_PLC_STATE).get_value()
            self._client.disconnect()
            return st == _opcpar.PLC_RUN_STATE_VALUE
        except Exception as e:
            print(f"PLC state read error (ignored): {e}")
            return e

    def _read_telemetry(self, print: bool = False):
        """
        Access the PLC telemetry data structure, for reading.

        Returns
        -------
            tel: telemetry object
        """
        with self.connected():
            try:
                node = self._client.get_node(_opcpar.NODE_ID_TELEMETRY)
                tel = node.get_value()

                ## FOR DEBUGGING PURPOSES ##

                def pretty_print_obj(obj: list, indent: int = 0, max_list: int = 30):
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

                if print is True:
                    print("\n--- Telemetry ---")
                    pretty_print_obj(tel)

                ###############################

                return tel

            except Exception as e:
                print(f"Telemetry read error: {e}")
                self._logger.error(f"Telemetry read error: {e}")
                return None

    def _write_struct(self, node: ua.Node, struct_instance: object):
        """
        Write a structured data instance to the given OPC UA node.
        """
        dv = ua.DataValue(ua.Variant(struct_instance, ua.VariantType.ExtensionObject))
        dv.SourceTimestamp = None
        dv.ServerTimestamp = None
        node.set_value(dv)

    def _connect_client(self):
        """
        Try to connect to the OPC UA server.
        """
        import time

        self._client = Client(url=_opcpar.server)
        self._client.set_user(_opcpar.username)
        self._client.set_password(_opcpar.password)
        self._client.set_security_string(
            "Basic256Sha256,SignAndEncrypt,%s,%s,%s"
            % (CLIENT_CERT_DER, CLIENT_KEY_PEM, SERVER_CERT_DER)
        )

        ## Try the connection to get errors early ##
        try:
            self._client.connect()
            time.sleep(0.5)
            self._client.disconnect()
        except Exception as e:
            self._logger.error(f"Failed to connect to OPC UA server: {e}")
            raise ConnectionError(f"Failed to connect to OPC UA server: {e}")

    ### HELPER FUNCTIONS FOR ENUMERATING THE OBECTS ###
    def __enum_members(self, enum_cls: object):
        # IntEnum-like
        try:
            return list(enum_cls)
        except Exception:
            members = []
            for k in dir(enum_cls):
                if k.startswith("_"):
                    continue
                v = getattr(enum_cls, k)
                try:
                    if isinstance(v, enum_cls):
                        members.append(v)
                except Exception:
                    pass
            if not members:
                raise RuntimeError(f"Enum {enum_cls} has no members")
            return members

    def __zero_or_first(self, enum_cls: object):
        """Return the enum member with value zero, or the first member if none has value zero."""
        mem = self.__enum_members(enum_cls)
        for m in mem:
            if int(m) == 0:
                return m
        return mem[0]

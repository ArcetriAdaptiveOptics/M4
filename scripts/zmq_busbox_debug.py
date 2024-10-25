"""
Author(s) 
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import zmq
import time
import struct
import numpy as np
import pprint
ip="tcp://192.168.22.51:6660"
read = bytearray(19)
read[0] = 17

# move_act_1 = bytearray(19)
# move_act_1[0] = 17 # byte del comando
# move_act_1[1] = 255 # byte dell'enable/disable x act (tutti its fine)
# move_act_1[2] = 255 # byte del move x act (tutti its fine, ma anche 7 x primi 3)
# move_act_1[3] = 1   # the following two bites should form a 500, which is the
# move_act_1[4] = 244 # commanded position to actuator 1 (takes a short, 2 bytes)
# [3:] sono, a due a due, tutti i parametri da dare agli attuatori, 
# che sono 16 bit (un int)
command_byte = 17
enable_disable_byte = 255
move_byte = 255
commanded_position = 500  # This is a short (2 bytes)

# velocità attuatori ~ 30(req) - 50 (meas) um/s
act_vel = 40 

# Pack the data into a byte array
move_act_1 = struct.pack('<BBBh', command_byte, enable_disable_byte, move_byte, commanded_position)
move_act_1 += bytearray(14)

def connect():
    context = zmq.Context()
    socket = context.socket(zmq.PAIR)
    try:
        socket.connect(ip)
    except zmq.ZMQError as e:
        print(f"Connection failed: {e}")
    return socket

sock = connect()

def iterative_actuators_test():
    cmd = [5, -5]*25 + [10,-10]*25 + [15,-15]*25 + [20,-20]*25
    actuators = np.zeros((3, len(cmd)), dtype=dict)
    for n,c in enumerate(cmd):
        for act in range(0,3):
            print(f"Actuator n.{act}")
            actuators[act][n] = command(c, act)
    return actuators

def sequential_actuators_test():
    cmd = [5, -5]*25 + [10,-10]*25 + [15,-15]*25 + [20,-20]*25
    actuators = np.zeros((3, len(cmd)), dtype=dict)
    for act in range(0,3):
        print(f"Actuator n.{act}")
        for n,c in enumerate(cmd):
            actuators[act][n] = command(c, act)
    return actuators
            
def command(position_um:int, act:int=0):
    act_positions = np.zeros(3, dtype=int)
    act_positions[act] = position_um
    cmd = struct.pack('<BBBhhh', 
                      command_byte,
                      enable_disable_byte,
                      move_byte, 
                      act_positions[0],
                      act_positions[1],
                      act_positions[2])
    cmd += bytearray(12)
    out = send(cmd)
    response = decode(out)
    act_pos = response['actuators'][act]['actual_position']
    act_enc_pos = response['actuators'][act]['encoder_position']
    act_targ_pos = position_um
    tot_time = 0
    while np.abs(act_targ_pos-act_enc_pos) > 4:
        waittime = 3
        #waittime = np.abs(act_pos-act_targ_pos)/act_vel
        tot_time += waittime
        print("waiting ", waittime)
        time.sleep(waittime)
        out = send(cmd)
        response = decode(out)
        act_pos = response['actuators'][act]['actual_position']
        act_enc_pos = response['actuators'][act]['encoder_position']
        print(f"actual position: {act_pos}\nencoder_position: {act_enc_pos}")
    check = send(read)
    check = decode(check)
    check['execution_time'] = tot_time
    #pprint.pprint(decode(check))
    return check

def send(bb):
    try:
        #print("Sending message\n", bb)
        sock.send(bb)
        out = sock.recv()
        #print("Response\n", out)
        return out
    except zmq.ZMQError as ze:
        print(ze)
        return "status 0"
    
def disconnect(socket):
    socket.disconnect(ip)
    socket.close()
    return
    
#  TCS Interface KAMAL
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

def decode(outMsg):
    act_struct = 'hhibH3b' # short*2, int, char, unsigned short, 3-char array
    struct_size = struct.calcsize('Bbbb')
    outer_data = struct.unpack('Bbbb', outMsg[:struct_size])
    heart_beat_counter, voltage, temperature, status = outer_data
    # Initialize a list to hold the actuator data
    actuators = []
    # Calculate the size of the inner struct
    act_size = struct.calcsize(act_struct)
    # Unpack each inner struct
    for i in range(3):
        start = struct_size + i * act_size
        end = start + act_size
        status_bytes = np.zeros(3, dtype=int)
        actual_position, target_position, encoder_position, actual_velocity, bus_voltage, status_bytes[0], status_bytes[1], status_bytes[2] = struct.unpack(act_struct, outMsg[start:end])
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
    

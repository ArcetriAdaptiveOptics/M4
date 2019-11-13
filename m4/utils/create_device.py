'''
@author: cs
'''

from m4.type import deformable_mirror

def myDevice(name_str):
    '''
    arg:
        device: optical element to be calibrated.
             ex. m4 = all the mirror, segment = one petal
    '''
    if name_str == "m4":
        device = deformable_mirror.m4()

    elif name_str == "segment":
        print('Segment number: ')
        input_number = int(input())
        device = deformable_mirror.segment(input_number)

    else:
        raise OSError('Device %s doesnt exists' % device)

    return device

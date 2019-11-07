'''
@author: cs
'''

from m4.type import deformableMirror

def myDevice(name_str):
    '''
    arg:
        device: elemento ottico da calibrare.
             ex. m4= tutto lo speccio, segment= un petalo
    '''
    if name_str == "m4":
        device = deformableMirror.m4()

    elif name_str == "segment":
        print('Segment number: ')
        input_number = int(input())
        device = deformableMirror.segment(input_number)

    else:
        raise OSError('Device %s doesnt exists' % device)

    return device

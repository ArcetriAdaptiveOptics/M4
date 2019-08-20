'''
@author: cs
'''

import m4 as _m4

def myDevice(nameStr):
    '''
    arg:
        device: elemento ottico da calibrare.
             ex. m4= tutto lo speccio, segment= un petalo
    '''
    if nameStr == "m4":
        device = _m4.allDevices.m4()
        deviceInfo= 'All segments'
        device= device, deviceInfo
    elif nameStr == "segment":
        print('Segment number: ')
        x = int(input())
        device = _m4.allDevices.segment(x)
        deviceInfo= 'Segment number %d' %x 
        device= device, deviceInfo 
    else:
        raise OSError('Device %s doesnt exists' % device)
    
    return device
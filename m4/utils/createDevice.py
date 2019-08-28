'''
@author: cs
'''

from m4 import deformableMirror

def myDevice(nameStr):
    '''
    arg:
        device: elemento ottico da calibrare.
             ex. m4= tutto lo speccio, segment= un petalo
    '''
    if nameStr == "m4":
        device = deformableMirror.m4()

    elif nameStr == "segment":
        print('Segment number: ')
        x = int(input())
        device = deformableMirror.segment(x)
         
    else:
        raise OSError('Device %s doesnt exists' % device)
    
    return device
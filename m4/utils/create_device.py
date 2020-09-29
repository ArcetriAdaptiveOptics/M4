'''
Autors
  - C. Selmi: written in 2019
'''

from m4.type import deformable_mirror

def myDevice(name_str):
    """
    Parameters
    ----------
        name_str: string
                name of optical element to be calibrated
                'm4' for all the mirror
                'segment' for one petal

    Returns
    -------
        device: object
                optical element to be calibrated
    """
    if name_str == "m4":
        device = deformable_mirror.m4()

    elif name_str == "segment":
        print('Segment number: ')
        input_number = int(input())
        device = deformable_mirror.segment(input_number)

    else:
        raise OSError('Device %s doesnt exists' % device)

    return device

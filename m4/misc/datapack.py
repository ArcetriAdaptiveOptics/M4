import numpy as np
def pack(img):
    '''
        img shall be a masked array
    '''
    m1 = np.max(img)
    m0 = np.min(img)
    dm = m1-m0
    pimg = (img-m0)/dm * 32767*2


    return pimg, exval

# test routine, in IDL!!!
#mm = minmax(img)
#dm = mm[1]-mm[0]
#q1 = img-mm[0]
#q1 = q1/dm
#q2 = q1*32767*2
#q3 = fix(q2, type=12)
#idx = where(mask ne 1)
#q3[idx] = 65535


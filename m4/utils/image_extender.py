'''
Autors
  - C. Selmi: written in 2019

Function used to extend the dimensions of the captured images to
those of the pupil on which I build the Zernike::

    from m4.utils import imageExterder
    image_ex = imageExterder.imageExterder(image)
'''
import numpy as np
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.configuration.ott_parameters import OttParameters

def imageExtender(cube_element):
    ''' 
    Parameters
    ----------
        cube_element: numpy masked array
                    an image from cube image

    Returns
    -------
        image: numpy masked array
            cube_element extended to the size of the Zernike pupil
    '''
    zg = ZernikeGenerator(2*OttParameters.PARABOLA_PUPIL_XYRADIUS[2])
    dim_y = (2 * zg.getRadius() - cube_element.shape[0]).astype(int) #512-500
    vv = np.ma.masked_array(np.zeros((dim_y, cube_element.shape[1])),
                            mask=np.ones((dim_y, cube_element.shape[1])).astype(bool)) #496
    dim_x = (2 * zg.getRadius() - cube_element.shape[1]).astype(int)   #512-496
    vv2 = np.ma.masked_array(np.zeros(((2 * zg.getRadius()).astype(int), dim_x)),
                            mask=np.ones(((2 * zg.getRadius()).astype(int), dim_x)).astype(bool))
    pp = np.ma.append(cube_element, vv, axis=0)
    image = np.ma.append(pp, vv2, axis=1)
    return image

'''
Authors
  - C. Selmi: written in 2021
'''
import os
import numpy as np
from astropy.io import fits as pyfits
from m4.ground import tracking_number_folder
from m4.configuration import config_folder_names as fold_name

class PointerAlign():
    '''
    '''

    def __init__(self, pointerId):
        self.id = pointerId
        self.tt = None
        self.dove = None

        self.beamSlope = None
        self.rotationAxsSl = None
        self.meanOffs = None
        self.offset1 = None
        self.offset2 = None
        self.exposureTimeInMilliSeconds = None
        self.l1 = None
        self.l2 = None
        self.p = None #distance between actuators
        self.defineParameters(pointerId)
        self.pixs = 7.2e-6
        self.rad2as = 206265 #[as/ras]

    def defineParameters(self, pointerId):
        '''
        Parameters
        ----------
        poniterId: string
            'NGS' or 'LGS'
        '''
        if pointerId == 'NGS':
            self.exposureTimeInMilliSeconds = 0.0002
            self.l1 = 0.495
            self.l2 = 1.163
            self.p = 0.038
            self.folder = os.path.join(fold_name.POINTER_ALIGN_ROOT_FOLDER, 'NGSPointer')
        if pointerId == 'LGS':
            self.exposureTimeInMilliSeconds = 0.001*5
            self.l1 = 0.485 #distance from II actuator and first measurement point
            self.l2 = 1.03 #distance from II actuator and second measurement point
            self.p = 0.033
            self.folder = os.path.join(fold_name.POINTER_ALIGN_ROOT_FOLDER, 'LGSPointer')

    def operationAndPrint(self, y00, y01, y11, y10):
        m0  = np.mean([y00,y01]) #mid-point at the first distance
        m1  = np.mean([y10,y11]) #mid point at the second distance
        r0  = y00 - m0
        r1  = y10 - m1
        s0  = (m1-m0)/(self.l2-self.l1) #slope of the rotation axis

        s   = (y10-m1-(y00-m0))/(self.l2-self.l1) #slope of the beam wrt the rotation axis
        z0  = np.mean([y10,y11]) #== m1, expected beam position at second distance whit no tilt and no decenter
        h2 = (y10-m1)-s*self.l2 #error (displacement) of the actuator II
        h1 = (y10-m1)-s*(self.l2+self.p) #error (displacement) of the actuator I
        z1 = m1-h1/self.p*self.l2 #expected beam position 

        self.beamSlope = s*self.pixs*self.rad2as
        self.rotationAxsSl = (m1-m0)/(self.l2-self.l1)*self.pixs*self.rad2as
        self.radii = (np.abs(y10-y11)-np.abs(y00-y01))*self.pixs*1e6
        self.meanOffs = np.mean([h1*self.pixs*1e6, h2*self.pixs*1e6])
        self.offset1 = h1*self.pixs*1e6
        self.offset2 = h2*self.pixs*1e6
        print('Tracknum: %s' %self.tt)
        print('Beam slope [arcsec]: %f' %self.beamSlope)
        print('Radii difference [um]: %f' %self.radii)
        print('Rotation axis slope [arcsec]: %f' %self.rotationAxsSl)
        print('Pointer offset at actuators [um]: %f %f' %(self.offset1, self.offset2) )

        return m0, m1

    @staticmethod
    def readInfo(tt):
        dove = tracking_number_folder.findTrackingNumberPath(tt)
        filename = os.path.join(dove, 'BenchLength-Info.fits')
        header = pyfits.getheader(filename)
        hduList = pyfits.open(filename)

        pointerId = header['ID']
        theObject = PointerAlign(pointerId)
        theObject.dove = dove
        theObject.tt = header['TT']
        theObject.beamSlope = header['BEAMSL']
        theObject.rotationAxsSl = header['ROTAXSL']
        theObject.meanOffs = header['MEANOFFS']

        offsets = hduList[0].data
        theObject.offset1 = offsets[0]
        theObject.offset2 = offsets[1]
        return theObject

    def saveInfo(self):
        txt_file_name = os.path.join(self.dove, 'BenchLength-Info.txt')
        file = open(txt_file_name, 'w+')
        file.write('tt = %s \n' %self.tt)
        file.write('Beam slope = %f \n' %self.beamSlope)
        file.write('Rot.axis sl. = %f \n' %self.rotationAxsSl)
        file.write('Mean offs. [um] = %f \n' %self.meanOffs)
        file.write('Offset I [um] = %f \n' %self.offset1)
        file.write('Offset II [um] = %f \n' %self.offset2)
        file.close()

        fits_file_name = os.path.join(self.dove, 'BenchLength-Info.fits')
        header = pyfits.Header()
        header['ID'] = self.id
        header['TT'] = self.tt
        header['BEAMSL'] = self.beamSlope
        header['ROTAXSL'] = self.rotationAxsSl
        header['MEANOFFS'] = self.meanOffs
        pyfits.writeto(fits_file_name, np.array([self.offset1, self.offset2]), header)

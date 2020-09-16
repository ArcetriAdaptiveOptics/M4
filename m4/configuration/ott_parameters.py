# This is a configuration file for 8s, the OTT Simulator
'''
Tracking numbers::

    tn_conf = {mirror:'20170203',optical:'20150730',mechanical:'20150730'} ; last mirror configuration with no slave acts, 892 acts
    tn_conf = {mirror:'20150730',optical:'20150730',mechanical:'20150730'} ; initial mirror configuration
    tn_conf = {mirror:'20170430',optical:'20150730',mechanical:'20150730'} ; mirror configuration with slave acts
    ott_folder = {mirror:'MIRROR_System/',optical:'OPTICAL_System/',mechanical:'MECHANICAL_System/', fea:'FEA/',zemax:'ZST/',config:tn_conf}
    interf: configuration parameters for the interferometer
'''
import numpy as np

tnconf_mirror = ''
tnconf = ''

class Interferometer():
    ''' Interferometer parameters
    '''
    N_PIXEL = np.array([512,512]) #2048
    BIN_PIX = 1
    HORIZ_CROP = 100
    VERT_CROP = 100
    WEDGE = 0.5
    WAVEL = 632.8e-9
    QUANTIZATION = 1

class OttParameters():
    ''' Optical tower parameters
    '''
    parab_radius    = 1.420/2 #was 1.440
    parab_dist      = 5.4
    rflat_dist      = 4.24
    rflat_radius    = 0.3
    fold_radius     = 0.025  #was 0.025
    frame2m4center  = 0.887
    segm_gap        = 0.002
    m4od            = 2.540  #was 2.387
    m4optod         = 2.387
    m4id            = 0.54
    outarea         = 4
    fullrslide      = 0.85
    segment_angle   = 60
    rflat_cell      = 0.01
    pscale  = Interferometer.N_PIXEL[0]/2/parab_radius
    parab_max_displacement = np.array([0,0,1,1,1,0]) #range of maximum allowed displacement
    rm_max_displacement = np.array([0,0,0,1,1,0]) #range of maximum allowed displacement
    m4_max_displacement = np.array([0,0,0,1,1,0]) #range of maximum allowed displacement


    PARABOLA_PUPIL_XYRADIUS = np.array([256, 256, 256])
    PARABOLA_DOF = np.array([2,3,4])
    RM_DOF = np.array([3,4])
    PIXEL_SCALE = 360.5 #PIXEL/METRI
    RADIUS_FIDUCIAL_POINT = 0.3

    M4_MECHANICAL_PUPIL_XYRADIUS = np.array([458, 458, 458]) #np.array([512, 512, 512])
    M4_OPTICAL_DIAMETER = 858
    N_SEG = 6
    N_ACTS_TOT = 5352
    N_ACT_SEG = 892
    M4_DOF = np.array([3,4])
    REFERENCE_ANGLE_RAD = np.pi / 3
    REFERENCE_ANGLE_DEGREES = 60
    SEGMENT_DISTANCE_FROM_CENTRE = 320
    DIAMETER_IN_PIXEL_FOR_SEGMENT_IMAGES = 512
    BIG_IMAGE_DIAMETER = 1236

    #SPL
    TN_FRINGES = '20181108_1'

class OpcUaParameters():
    ''' Numbers for opc ua parameters '''
    RA = 0 #angolo di rotazione
    CAR = 1 #carrello dell RM
    ST = 2 #slitta della parabola
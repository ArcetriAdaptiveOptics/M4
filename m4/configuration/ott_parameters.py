# This is a configuration file for 8s, the OTT Simulator
'''
Tracking numbers::

    tn_conf = {mirror:'20170203',optical:'20150730',mechanical:'20150730'} ;
    last mirror configuration with no slave acts, 892 acts
    tn_conf = {mirror:'20150730',optical:'20150730',mechanical:'20150730'} ;
    initial mirror configuration
    tn_conf = {mirror:'20170430',optical:'20150730',mechanical:'20150730'} ;
    mirror configuration with slave acts
    ott_folder = {mirror:'MIRROR_System/',optical:'OPTICAL_System/',mechanical:'MECHANICAL_System/',
                  fea:'FEA/',zemax:'ZST/',config:tn_conf}
    interf: configuration parameters for the interferometer
'''
import numpy as np

tnconf_mirror = ''
tnconf = ''

class Sound():
    AUDIO_FILE_PATH = '/mnt/m4storage/Data/Audio'
    PLAY = True

class Interferometer():
    ''' Interferometer parameters
    '''
    i4d_IP =  '192.168.22.78' #'193.206.155.193'
    i4d_port = 8011
    N_PIXEL = np.array([512, 512]) #2048
    BIN_PIX = 1
    HORIZ_CROP = 100
    VERT_CROP = 100
    WEDGE = 0.5
    WAVEL = 632.8e-9
    QUANTIZATION = 1
    BURST_FREQ = 20.0 #28.57 #PhaseCam4020
    CAPTURE_FOLDER_NAME_4D_PC = 'D:/M4/Capture'
    PRODUCE_FOLDER_NAME_4D_PC = 'D:/M4/Produced'
    PRODUCE_FOLDER_NAME_M4OTT_PC = '/home/m4/4d/M4/Produced'
    SETTINGS_CONF_FILE_M4OTT_PC = '/home/m4/4dConfig/AppSettings.ini'

class OttParameters():
    ''' Optical tower parameters
    '''
    par_rm_coef_for_coma_measuremets = -2.05

    parab_radius = 1.420/2 #was 1.440
    parab_dist = 5.4
    rflat_dist = 4.24
    rflat_radius = 0.3
    fold_radius = 0.025  #was 0.025
    frame2m4center = 0.887
    segm_gap = 0.002
    m4od = 2.540  #was 2.387
    m4optod = 2.387
    m4id = 0.54
    outarea = 4
    fullrslide = 0.85
    segment_angle = 60
    rflat_cell = 0.01
    pscale = Interferometer.N_PIXEL[0]/2/parab_radius
    parab_max_displacement = np.array([0, 0, 3, 10, 10, 0]) #range of maximum allowed displacement
    rm_max_displacement = np.array([0, 0, 0, 10, 10, 0]) #range of maximum allowed displacement
    m4_max_displacement = np.array([0, 0, 0, 1, 1, 0]) #range of maximum allowed displacement


    PARABOLA_DOF = np.array([2, 3, 4])
    RM_DOF_PISTON = np.array([2, 3, 4])
    RM_DOF = np.array([3, 4])
    PIXEL_SCALE = 360.5 #PIXEL/METRI
    RADIUS_FIDUCIAL_POINT = 0.5
    INNER_MARKERS_REJECTION_RADIUS = 100


    M4_MECHANICAL_PUPIL_XYRADIUS = np.array([458, 458, 458]) #np.array([512, 512, 512])
    M4_OPTICAL_DIAMETER = 858
    M4_DOF = np.array([3, 4])
    REFERENCE_ANGLE_RAD = np.pi / 3
    REFERENCE_ANGLE_DEGREES = 60
    SEGMENT_DISTANCE_FROM_CENTRE = 320
    DIAMETER_IN_PIXEL_FOR_SEGMENT_IMAGES = 512
    BIG_IMAGE_DIAMETER = 1236

    #SPL
    TN_FRINGES = '20181108_1'

class M4Parameters():
    ''' '''
    N_SEG = 6
    N_ACTS_TOT = 5352
    N_ACT_SEG = 892


    V_MATRIX_FOR_SEGMENT_ROOT_811 = \
        '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions/20170630_105105/modeMatrix.fits'

    M4COORDINATE_ROOT_FOLDER = \
        '/Users/rm/Desktop/Arcetri/M4/ActuatorCoordinates.fits'


class OpcUaParameters():
    ''' Numbers for opc ua parameters '''
    server = "opc.tcp://192.168.22.100:48050"
    num_PT_sensor = 24
    min_angle, max_angle = np.array([-171, 181])
    min_r_slide, max_r_slide = np.array([-9000, 9000])  #modified RB 20210423: was -50; 9000
    min_slide, max_slide = np.array([-10000, 10000])

    RA = 0 #angolo di rotazione
    CAR = 1 #carrello dell RM
    ST = 2 #slitta della parabola
    RM1 = 3
    RM2 = 4
    RM3 = 5
    PAR1 = 6
    PAR2 = 7
    PAR3 = 8

    RM_PISTON = 11
    RM_TIP = 9
    RM_TILT = 10
    PAR_PISTON = 14
    PAR_TIP = 12
    PAR_TILT = 13

    RM_KIN = 9
    PAR_KIN = 10

    zabbix_variables_name = ['RA', 'CAR', 'ST',
                             'RM1', 'RM2', 'RM3',
                             'PAR1', 'PAR2', 'PAR3',
                             'RM_TIP', 'RM_TILT', 'RM_PISTON',
                             'PAR_TIP', 'PAR_TILT', 'PAR_PISTON']

    zabbix_hostname = 'M4OTT'
    zabbix_server = '192.168.22.22'
    zabbix_port = 10051
    accelerometers_server = 'tcp://192.168.22.100:6660'
    accelerometers_data_folder = '/mnt/acc_data'
    accelerometers_dt_plc = 2.5e-4
    accelerometers_dt = 5e-3
    accelerometers_sn = ['', '', '', '', 'a', 'a', 'a', 'b']
    accelerometers_plc_id = np.array([5,6,7,8]) #final vector is time, accelerom, --> first accelerom is 1
    accelerometrs_directions = ['', '', '', '', 'X', 'Z', 'Y', 'Z']
    accelerometers_sensitivity = np.array([0.,0,0,0,0.1,0.1,0.1,10])   #accel sentivity [V/g]
    accelerometers_plc_range = np.array([0.,0,0,0,0.32,0.32,0.32,1.28]) #this is the full range of measurement. i.e. the total span, where the signal is sampled with 2**24 counts
    accelerometers_plc_totcounts = 2**24 #tot counts is 2**24, for positive and negative range
 
class OtherParameters():
    ''' '''
    MASK_INDEX_SIMULATORE = 3
    MASK_INDEX_TOWER = 0

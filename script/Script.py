'''
Authors
  - L. Oggioni: written in 2022
'''

from m4.configuration import start
import numpy as np
import os
from matplotlib import pyplot as plt
from m4.ground import read_data
from m4.configuration import config_folder_names as fold_name
from m4.ground import tracking_number_folder
from m4 import noise


conf='G:\Il mio Drive\Lavoro_brera\M4\LucaConfig.yaml'
ott, interf, dm = start.create_ott(conf)

    
def main_test():
    
    '''
    script di test per verificare che il software faccia quello che vogliamo
    
    se ind = False non aggiunge indeterminazione di lambda nella lettura dell interferometro
    
    '''
    
    ind=True 
    plot=False
    
    move2petalo(6,0)
     
    ampiezza=10e-9; rel=False 
    dm.act_random(ampiezza,rel)
    immagine(1,True,ind)
    
    inc=100e-9;rel=True
    dm.act_incr(inc,rel)
    immagine(2,True,ind)
    
    rel=False
    #dm.act_zero(rel)
    dm.act_zero()
    immagine(3,True,ind)

def test_noise():
    
    data_file_path = fold_name.SIMUL_DATA_CALIB_DM_FOLDER+"\\PassettiWithNoise\\20220801_110511"
    numbers_array=np.array([2,3,5,7])
    tidy_or_shuffle=0    
    noise.noise_vibrations(data_file_path, numbers_array, tidy_or_shuffle)

def main_passetti_and_Unwrap(Nstep,Nact):
    
    '''  
    simulo l'acquisizione di immagini con 
    l'interferometro ed estrae la posizione 
    di un determinato attuatore. Poi unwrappo

 
    Parameters
        ----------
        Nstep: integer
            number of images
            
        Nact: integer
            actuator number
             
    '''    
  
    D=0;
    cartella_imm = passetti_con_Rumore(Nstep,D,freq=25,ind=True)
    cartella_act='PositionActuators\\20220519_130221'

    Z, mask, mask2=passetti_postporocess(cartella_imm,cartella_act,Nact)

    M=np.size(Z)
    x=np.linspace(1,M,M)

    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.imshow(mask+mask2)
    ax1.set_title("position actuator "+str(Nact))

    U=unwrapp(Z, 632.8e-9, 150e-9)
    ax2.plot(x,U,'-x')
    ax2.plot(x,Z,'-x')
    plt.xlabel('step')
    plt.ylabel('position (m)')
    ax2.set_title('unwrap result')
    
    plt.savefig(cartella_imm+"/results.png")
    
    mat=np.array([Z,U])
    mat=mat.transpose()
    np.savetxt(cartella_imm+'/position_Actuator'+str(Nact)+'.txt',mat,header='original    unwrap')
   
def main_Unwrap(Nact): 
    '''  
    Faccio l'unwrap delle immagini generate con "passetti" o "passettiwithnoise"

 
    Parameters
        ----------
            
        Nact: integer
            actuator number
             
    '''    
  
    D=0;
    dir0 = fold_name.SIMUL_DATA_CALIB_DM_FOLDER+"\\PassettiWithNoise"
    cartella_imm=os.path.join(dir0,"20220801_110511")
    cartella_act='PositionActuators\\20220519_130221'

    Z, mask, mask2=passetti_postporocess(cartella_imm,cartella_act,Nact)

    M=np.size(Z)
    x=np.linspace(1,M,M)

    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.imshow(mask+mask2)
    ax1.set_title("position actuator "+str(Nact))

    U=unwrapp(Z, 632.8e-9, 150e-9)
    ax2.plot(x,U,'-x')
    ax2.plot(x,Z,'-x')
    plt.xlabel('step')
    plt.ylabel('position (m)')
    ax2.set_title('unwrap result')
    
    plt.savefig(cartella_imm+"/results.png")
    
    mat=np.array([Z,U])
    mat=mat.transpose()
    np.savetxt(cartella_imm+'/position_Actuator'+str(Nact)+'.txt',mat,header='original    unwrap')
   
    
def passetti(N=1,D=20e-9,ind=True):
    
    '''
    
    vado avanti a passetti di D su tutti gli attuatori
    e salvo le immagini in una cartella 
    
    Parameters
        ----------
        N: integer 
            number of steps
        D: float
            ampiezza passo sugli attuatori (nm)
        ind: logic
            ad or not indetermination on the measure
    '''

    move2petalo(num=6,RM=0)

    dir0 = fold_name.SIMUL_DATA_CALIB_DM_FOLDER+"\\Passetti"
    dir , tt =tracking_number_folder.createFolderToStoreMeasurements(dir0)



    inc=D; rel=True; plot=True

    dm.act_incr(inc,rel) 
    ima=interf.acquire_phasemap()

    inc=D; rel=True; plot=False

    dm.act_incr(inc,rel) 
    ima=interf.acquire_phasemap()
    interf.save_phasemap(dir,'0000.fits',ima)

    np.save(os.path.join(dir, '0mask'),ima.mask)


    for x in range(0, N):
        print("step "+str(x+1)+"/"+str(N))
        dm.act_incr(inc,rel) 
        ima=interf.acquire_phasemap()
        #imname=f"{xx+1:04}"+".fits"
        imname="{xx+1:04}"+".fits"
        interf.save_phasemap(dir,imname,ima)
        
    return dir    

        
#    for x in range(1, 6):
#        imname=str(x)+".fits"
#        data=read_data.readFits_data(os.path.join(dir, imname))
#        mask=np.load(os.path.join(dir,'00mas               
#        plt.figure()
#        plt.clf()
#        plt.imshow(data)
#        plt.imshow(mask)
#        plt.colorbar()



def passetti_con_Rumore(N=1,D=20e-9,freq=25,ind=True):

    '''
    vado avanti a passetti di 20 nm su tutti gli attuatori
    e salvo le immagini in una cartella definita da dir
    Ad ogni step aggiungo poi un rumore oscillante dovuto al 
    movimento dei vari componenti (per ora solo la parabola)

    Parameters
        ----------
        N: integer 
            number of steps
        D: float
            ampiezza passo sugli attuatori (nm)
        freq: integer
            sampling frequency (Hz)
        ind: logic
            ad or not indetermination on the measure
    '''    
      
    
    move2petalo(num=6,RM=0)
    dir0 = fold_name.SIMUL_DATA_CALIB_DM_FOLDER+"\\PassettiWithNoise"
    dir , tt =tracking_number_folder.createFolderToStoreMeasurements(dir0)
    
        
    inc=D; rel=True; plot=False
                
    dm.act_incr(inc,rel) 
    ima=interf.acquire_phasemap()
    interf.save_phasemap(dir,'img_0000.fits',ima)


    np.save(os.path.join(dir, '0mask'),ima.mask)

    
    '''
    compiono a frequenza freq 
    '''
    t=np.linspace(0,N/freq,N)
    
    #definisco i parametri dei rumori sulla posizione della parabola
    
    # parametri ricavati da simulazioni ADS, europa no dampers
    #x=4e-5*1e3; Frx=14; start_x=2*np.pi*np.random.rand(); # A in mm, Fr in Hz 
    #y=3e-5*1e3; Fry=62; start_y=2*np.pi*np.random.rand(); # A in mm, Fr in Hz 
    #z=2e-5*1e3; Frz=58; start_z=2*np.pi*np.random.rand(); # A in mm, Fr in Hz 
    #tx=3e-6*1e3; Frtx=59; start_tx=2*np.pi*np.random.rand(); # A in mm, Fr in Hz 
    #ty=6e-6*1e3; Frty=62; start_ty=2*np.pi*np.random.rand(); # A in mm, Fr in Hz 
    
    
    Ax=2.1e-8*1e3; Frx=14; start_x=2*np.pi*np.random.rand(); # A in mm, Fr in Hz 
    Ay=2.2e-8*1e3; Fry=62; start_y=2*np.pi*np.random.rand(); # A in mm, Fr in Hz 
    Az=2.1e-8*1e3; Frz=158; start_z=2*np.pi*np.random.rand(); # A in mm, Fr in Hz 
    Atx=2.2e-8*1e3; Frtx=10; start_tx=2*np.pi*np.random.rand(); # A in mm, Fr in Hz 
    Aty=2.4e-8*1e3; Frty=90; start_ty=2*np.pi*np.random.rand(); # A in mm, Fr in Hz  
 
    plt.figure()

    l1,=plt.plot(t,1e-3*Ax*np.sin(2*np.pi*Frx*t+start_x))
    l2,=plt.plot(t,1e-3*Ay*np.sin(2*np.pi*Fry*t+start_y))
    l3,=plt.plot(t,1e-3*Az*np.sin(2*np.pi*Frz*t+start_z))
    l4,=plt.plot(t,1e-3*Atx*np.sin(2*np.pi*Frtx*t+start_tx))
    l5,=plt.plot(t,1e-3*Aty*np.sin(2*np.pi*Frty*t+start_ty))
    plt.xlabel('time(s)')
    plt.ylabel('noise')
    plt.legend((l1, l2, l3, l4, l5),('x(m)', 'y(m)', 'z(m)', 'tx(rad)', 'ty(rad)'))
    
    plt.savefig(dir+"/noise.png")
    
    f = open(dir+"/log.txt", "a")
    f.write("Simulation parameters\n")
    f.write("\nAcquisition frequency:    "+str(freq)+"Hz")
    f.write("\nDelta actuators position:    "+str(D)+"m")
    f.write("\nNumber of acquisitions:    "+str(N))
    f.write("\nAmplitues noise (mm,mrad) (x,y,z,tx,ty):    "+str(Ax)+"    "+str(Ay)+"    "+str(Az)+"    "+str(Atx)+"    "+str(Aty))
    f.write("\nFrequencies noise (Hz) (x,y,z,tx,ty):    "+str(Frx)+"    "+str(Fry)+"    "+str(Frz)+"    "+str(Frtx)+"    "+str(Frty))
    f.close()
    
    for xx in range(0, N):
        print("step "+str(xx+1)+"/"+str(N))
        
        x=Ax*np.sin(2*np.pi*Frx*t[xx]+start_x)
        y=Ay*np.sin(2*np.pi*Fry*t[xx]+start_y)
        z=Az*np.sin(2*np.pi*Frz*t[xx]+start_z)
        ty=Atx*np.sin(2*np.pi*Frty*t[xx]+start_tx)
        tx=Aty*np.sin(2*np.pi*Frtx*t[xx]+start_ty)      
        
        ott.parabola.setPosition(np.array([x,y,z,tx,ty,0]))
        
        dm.act_incr(inc,rel) 
        ima=interf.acquire_phasemap()
        imname="img_"+f"{xx+1:04}"+".fits"
        interf.save_phasemap(dir,imname,ima)
    return dir
        
def passetti_postporocess(cartella_imm,cartella_act,N): 
    
    '''
    segue la posizione media definita dalla maschera dell'attuatore N
    trovata nella cartella_act, nelle immagini trovate nella cartella_imm 
    
    ATTENZIONE: NON E' ASSICURATO CHE MASCHERA E IMMAGINI SIANO ALLINEATE.
    DIPENDE DA COME SONO STATE GENERATE
    '''
    
    dir0=fold_name.SIMUL_DATA_CALIB_DM_FOLDER
    dir1=os.path.join(dir0, cartella_act, 'PositionActuator_mask.fits' )   
    act_pos_matrix=read_data.readFits_data(dir1)
    mask=act_pos_matrix[:,N-1].reshape([512,512])
    
        
    dir2=cartella_imm   
    D= os.listdir(dir2)          
    mask2=np.load(os.path.join(dir2,'0mask.npy'))    
    
    
    z=np.zeros(int(D[-1][4:-5]))
    for x in range(0, int(D[-1][4:-5])):
        print("step "+str(x+1)+"/"+str(int(D[-1][4:-5])))
        dir3=os.path.join(dir2,"img_"+f"{x+1:04}"+".fits")
        im=read_data.readFits_data(dir3)
        
        imm=np.ma.masked_where(mask==1, im)
        
        z[x]=imm.mean()
        
    return z, mask, mask2
        
        
        

def move2petalo(num=1,RM=0):
    
    '''move the interferometer to a specific section of the DM
    
    Parameters
        ----------
        num: integer 
            number of the target section
        RM: integer
            Reference mirror in(1) or not(0) 
    '''
    
    ott.parabolaSlider.setPosition(844)
    if RM==1:
        ott.referenceMirrorSlider.setPosition(844)
    if RM==0:
        ott.referenceMirrorSlider.setPosition(0)

    ott.angleRotator.setPosition(30+60*(num-1))
    
    return
    


def immagine(numFig=0,plotta=True,indet=True):
    
    ''' plot immagine
    '''
    
    ima=interf.acquire_phasemap(indet)
    
    if plotta==True:
        
        if numFig==0:
            plt.figure()
        else:
            plt.figure(numFig)
        plt.clf()
        plt.imshow(ima)
        plt.colorbar()

    return ima

def unwrapp(x, phase, threshold):
    
    '''  Apply to the positions vector x the unwrapping algorithm. 
    
    Parameters
    ----------
    x: float
        vector of positions extracted from the interferometric measures (m)
    phase: float
        wavelength of the interferometer (m)
    threshold: float
        greater then the increments between the x elements. lower then the phase (m)
        
    '''

    a=np.insert(np.diff(x),0,x[0])
    
    dx=a-phase*np.fix(a/phase)
    x_un = x-phase*np.fix(x/phase) 
    
    n = x.size
    
    if np.absolute(x_un[0])>threshold:
        x_un[0] = x_un[0] - np.sign(dx[0])*phase
      
      
    for ii in range (1,n):
        if np.absolute(dx[ii])<threshold:
            x_un[ii] = x_un[ii-1] + dx[ii];
        else:
            x_un[ii] = x_un[ii-1] - np.sign(dx[ii])*(phase-np.absolute(dx[ii]));
    
    
    return x_un
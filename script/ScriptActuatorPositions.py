from m4.configuration import start
import numpy as np
import os
from matplotlib import pyplot as plt
from m4.ground import read_data
from m4.configuration import config_folder_names as fold_name
from m4.ground import tracking_number_folder

conf='G:\Il mio Drive\Lavoro_brera\M4\LucaConfig.yaml'
ott, interf, dm = start.create_ott(conf)


def DefineActuatorMask_interf_PointOfView():
    
    '''
    applica le IF zonali (un attuatore alla volta) 
    e crea una maschera attorno ad ognuno applicando una treshold
    
    poi salva in una matrice complessiva 512*512 x 892 sia le immagini
    che le maschere 
    
    CONTROLLARE CHE FUNZIONI SUGLI ATTUATORI AL BORDO!!! 
    '''
        
    move2petalo(num=6,RM=0)
    
    date=datetime.datetime.now()
    name=str(date.year)+str(date.month)+str(date.day) 
    dir0 = fold_name.SIMUL_DATA_CALIB_DM_FOLDER+"\\PositionActuators"
    dir , tt =tracking_number_folder.createFolderToStoreMeasurements(dir0)
    
    act_pos_matrix=np.zeros([512*512,892])
    act_zonal=np.zeros([512*512,892])
    
    for x in range(0, 892):
        print("step "+str(x+1)+"/892")
        ampiezza=(10*1e-9)*1e6
        comm=np.zeros(892)
        comm[x]=ampiezza
        rel=False 
        dm.setActsCommand(comm,rel)
        indet=False
        ima=interf.acquire_phasemap(indet)
        
        act_zonal[:,x]=ima.data.flatten()
                
        im=np.ma.masked_where(ima < 5e-9 , ima)
        act_mask=im.mask.flatten() 
        act_pos_matrix[:,x]=act_mask
        
    hdu=fits.PrimaryHDU(act_pos_matrix)
    hdu.writeto(os.path.join(dir, 'PositionActuator_mask.fits'))
    hdu=fits.PrimaryHDU(act_zonal)
    hdu.writeto(os.path.join(dir, 'PositionActuator_data.fits'))
    
#    np.save(os.path.join(dir, 'PositionActuator_mask'),act_pos_matrix)
#    np.save(os.path.join(dir, 'PositionActuator_data'),act_zonal)
            
    return act_pos_matrix, act_zonal


    
        
def Check_actposition(cartella,N,arg=1):
    
    name=cartella
    dir0='G:/Il mio Drive/Lavoro_brera/M4/'
    if arg==1:
        dir=os.path.join(dir0, name, 'PositionActuator_data.fits' )    
    if arg==0:
        dir=os.path.join(dir0, name, 'PositionActuator_mask.fits' )
    
    act_pos_matrix=read_data.readFits_data(dir)
    
    im=act_pos_matrix[:,N-1].reshape([512,512])
          
    return im        
    
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
"""
Author(s):
    - Pietro Ferraiuolo : written in 2024

Description
-----------

How to Use it
-------------
"""
import logging
import numpy as np
from scripts.ottalign import _m4ac as mac # to change
from m4.ground import zernike as zern, read_data as rd

class Alignment():
    """
    Class for the alignment procedure: calibration and correction.
    """
    def __init__(self, mechanical_devices, acquisition_devices):
        """The Constructor"""
        self.ott        = mechanical_devices
        self.ccd        = acquisition_devices
        self.cmdMat     = rd.readFits_data('/home/pietrof/git/M4/scripts/ottalign/cmdMat.fits') #temporary
        self._moveFnc   = self._get_callables(self.ott, mac.devices_move_calls)
        self._readFnc   = self._get_callables(self.ott, mac.devices_read_calls)
        self._acquire   = self._get_callables(self.ccd, mac.ccd_acquisition)
        self._devName   = mac.names
        self._dof       = mac.dof
        self._dofTot    = mac.cmdDof
        self._idx       = mac.slices
        self._readPath  = mac.base_read_data_path
        self._writePath = mac.base_write_data_path
        self._zvec2fit  = np.arange(1,11)
        self._zvec2use  = [1,2,3,6,7]
        self._template  = [+1,-2,+1]
        self._auxMask   = None#rd.readFits_data(mac.calibrated_parabola)
        self._cmdAmp    = None        
        logging.basicConfig(filename=(self._writePath+'/alignment.log'), 
                                      level=logging.DEBUG,
                                      format='%(asctime)s - %(levelname)s - %(message)s')
        
    def calibrate_alignment(self, cmdAmp, template:list=None, n_repetitions:int=1):
        """


        Parameters
        ----------
        cmdAmp : float or ndarray
            DESCRIPTION.
        template : list, optional
            DESCRIPTION. The default is None.
        n_repetitions : int, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        self._cmdAmp = cmdAmp
        template = template if template is not None else self._template
        results = self.images_production(template, n_repetitions)
        return results, np.array([img.std() for img in results]) #!!! to remove
    
    def images_production(self, template, n_repetitions):
        """
        

        Parameters
        ----------
        template : TYPE
            DESCRIPTION.
        n_repetitions : TYPE
            DESCRIPTION.

        Returns
        -------
        n_results : TYPE
            DESCRIPTION.

        """
        results = []
        n_results = []
        for i in range(n_repetitions):
            logMsg = ''
            logMsg += f"Repetition n.{i}\n"
            for k in range(self.cmdMat.shape[1]):
                logMsg += f"Matrix Column {k+1} : {self.cmdMat.T[k]}"
                logging.info(logMsg)
                imglist = self._img_acquisition(k, template)
                image = self._push_pull_redux(imglist, template)
                results.append(image)
            if n_repetitions!=1:
                n_results.append(results)
            else:
                n_results = results
        return n_results

    def zern_calculator(self, imglist):
        """
        

        Parameters
        ----------
        imglist : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        coefflist = []
        for i in imglist:
            if self._auxMask is None:
                coeff, _ = zern.zernikeFit(i, self._zvec2fit)
            else:
                coeff, _ = zern.zernikeFitAuxmask(i, self._auxMask, self._zvec2fit)
            coefflist.append(coeff[self._zvec2use])
            return np.array(coefflist)  #intmat

    def apply_command(self, fullCmd):
        """
        Function which takes

        Parameters
        ----------
        fullCmd : list or ndarray
            Full command of the interaction matrix which command all device's
            available motors..
        """
        logMsg = '' #!!!
        device_commands = self._extract_cmds_to_apply(fullCmd)
        for cmd,fnc,dev in zip(device_commands,self._moveFnc,self._devName):
            if cmd.to_ignore:
                logMsg += f'Skipping null command for {dev}\n' # debug
                logging.info(f'Skipping null command for {dev}\n')
            else:
                try:
                    logMsg += f"Commanding {cmd} to {dev}\n"# debug
                    logging.info(f"Commanding {cmd} to {dev}\n")
                    fnc(cmd.vect)
                except Exception as e:
                    logging.warning(f"Someting went wrong with {dev}: ", e)
        logMsg += '\n'+'-'*30 # debug
        print(logMsg) #!!! debug only

    def read_positions(self):
        """


        Returns
        -------
        pos : TYPE
            DESCRIPTION.

        """
        logMsg = ''
        pos = []
        logMsg += "Current Positions\n"
        for fnc,dev_name in zip(self._readFnc, self._devName):
            temp = fnc()
            pos.append(_Command(temp))
            logMsg += f"{dev_name}"+' '*(16-len(dev_name))+f" : {temp}\n"
        logMsg += '-'*30
        logging.info(logMsg)
        print(logMsg) #!!! debug only
        return pos

    def _extract_cmds_to_apply(self, fullCmd):
        """


        Parameters
        ----------
        fullCmd : TYPE
            DESCRIPTION.

        Returns
        -------
        device_commands : TYPE
            DESCRIPTION.

        """
        commands = []
        for d,dof in enumerate(self._dof):
            dev_cmd = np.zeros(self._dofTot)
            dev_idx = fullCmd[self._idx[d]]
            for i,idx in enumerate(dev_idx):
                dev_cmd[dof[i]] = idx
            commands.append(_Command(dev_cmd))
        positions= self.read_positions()
        device_commands = []
        for pos,cmd in zip(positions, commands):
            res_cmd = pos+cmd 
            if np.array_equal(res_cmd.vect, pos.vect): 
                res_cmd.is_reset=False 
            device_commands.append(res_cmd)
        return device_commands

    def _img_acquisition(self, k, template):
        """
        

        Parameters
        ----------
        k : TYPE
            DESCRIPTION.
        template : TYPE
            DESCRIPTION.

        Returns
        -------
        imglist : TYPE
            DESCRIPTION.

        """
        imglist = [self._acquire[0](15)]
        for t in template:
            logMsg = ''
            logMsg += f"t = {t}"
            cmd = self.cmdMat.T[k] * self._cmdAmp * t
            logMsg += f" - Full Command : {cmd}"
            logging.info(logMsg)
            self.apply_command(cmd)
            imglist.append(self._acquire[0](15))
        print(logMsg) #!!! debug only
        return imglist

    def _push_pull_redux(self, imglist, template):
        """
        

        Parameters
        ----------
        imglist : TYPE
            DESCRIPTION.
        template : TYPE
            DESCRIPTION.

        Returns
        -------
        image : TYPE
            DESCRIPTION.

        """
        template.insert(0,-1)
        image = np.zeros((imglist[0].shape[0], imglist[0].shape[1]))
        for x in range(1, len(imglist)):
            opd2add = imglist[x] * template[x] + imglist[x-1] * template[x-1]
            mask2add = np.ma.mask_or(imglist[x].mask, imglist[x-1].mask)
            if x==1: master_mask = mask2add
            else: master_mask = np.ma.mask_or(master_mask, mask2add)
            image += opd2add
        image = np.ma.masked_array(image, mask=master_mask) / self._cmdAmp
        template.pop(0)
        return image
    
    @staticmethod
    def _get_callables(device, callables):
        """
        Returns a list of callables for the instanced object, taken from the
        configuration .py file

        Returns
        -------
        functions : list
            List of callables, which interacts with the input object of the class.
        """
        functions = []
        for dev_call in callables:
            obj, *methods = dev_call.split('.')
            call = getattr(device, obj)
            for method in methods:
                call = getattr(call, method)
            functions.append(call)
        return functions        

class _Command:
    """
    """
    def __init__(self, vector=None, is_reset=False):
        """The constructor"""
        self.update(vector, is_reset)

    def __repr__(self):
        return f"Command({self.vect}, is_reset={self.is_reset})"
    
    def __str__(self):
        return self.vect.__str__()
    
    def __add__(self, other):
        if not isinstance(other, _Command):
            return NotImplemented
        elif not isinstance(self.vect, np.ndarray) and not isinstance(other.vect, np.ndarray):
            raise NotImplementedError(f"Operation not supported for operands types {type(self.vect)} and {type(other)}")
        else:
            combined_vect = self.vect + other.vect
        is_reset = True if np.all(combined_vect==0) else False
        return _Command(combined_vect, is_reset)
    
    @property
    def to_ignore(self):
        """
        Returns True if the command is flagged to be ignored, False if it is null
        but should not be ignored.
        """
        return self.is_null and not self.is_reset
    
    @property
    def is_null(self):
        """
        Determins whether the command is null, i.e a sequence of zeros.
        """
        return np.all(self.vect == 0)

    def update(self, vector, is_reset=None):
        """
        Recycle method to reuse the class, in order to not overload the memory.
        """
        self.vect = np.array(vector) if isinstance(vector, list) else vector
        self.is_reset = is_reset

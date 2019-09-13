'''
@author: cs
'''


# TEST
        folder= os.path.join(Configuration.CALIBRATION_ROOT_FOLDER, "IFFunctions", tt)
        who, tt_cmdH, actsVector, cmdMatrix, cmdAmpl, nPushPull, indexingList= IF.loadInfo(folder)
        from m4.type.commandHistory import CmdHistory
        cmdH= CmdHistory.load(device, tt)
        commandHistory= cmdH._cmdHToApply
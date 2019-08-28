from m4.utils.configuration import Configuration




class ModalAmplitude():
    
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                   "ModalAmplitude")


    def save(self, tag):
        storeInFolder= ModalAmplitude._storageFolder()
        pass


    def load(filename):
        pass
        
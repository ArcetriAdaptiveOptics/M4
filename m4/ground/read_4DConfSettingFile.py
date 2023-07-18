import configparser
import numpy as np

class ConfSettingReader():
    '''
    How to use it:
        from m4.ground.read_4DConfSettingFile import ConfSettingReader
        cr = ConfSettingReader(file_path)
    '''

    def __init__(self, file_path):
        self.config = configparser.ConfigParser()
        self.config.read(file_path)
        self.camera_section = 'ACA2440'
        self.path_section = 'Paths'
    
    
    #CAMERA
    def getFrameRate(self):
        frame_rate = self.config.get(self.camera_section, 'FrameRate')
        return np.float(frame_rate)
    
    def getImageWidhtInPixels(self):
        image_width_in_pixels = self.config.get(self.camera_section, 'ImageWidthInPixels')
        return np.int(image_width_in_pixels)
    
    def getImageHeightInPixels(self):
        image_height_in_pixels = self.config.get(self.camera_section, 'ImageHeightInPixels')
        return np.int(image_height_in_pixels)
    
    def getOffsetX(self):
        offset_x = self.config.get(self.camera_section, 'OffsetX')
        return np.int(offset_x)
    
    def getOffsetY(self):
        offset_y = self.config.get(self.camera_section, 'OffsetY')
        return np.int(offset_y)

    def getPixelFormat(self):
        pixel_format = self.config.get(self.camera_section, 'PixelFormat')
        return pixel_format
    
    #PATH
    def getUserSettingFilePath(self):
        user_setting_file_path = self.config.get(self.path_section, 'UserSettingsFilePath')
        return user_setting_file_path



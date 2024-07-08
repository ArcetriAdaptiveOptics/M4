"""
Author(s)
    - Chiara Selmi
    
Description
-----------
Module containing the class for reading setting configuration files for the int
erferometers.

Hoe to Use it
-------------
Import the module and initialize the class with the complete filepath of the camera
settings file (this can be retrieved through a tracking number and the use of 
''m4.utils.osutils.getConf4DSettingsPath''.

>>> from m4.utils.osutils import getConf4DSettingsPath
>>> from m4.ground.read_4DConfSettingFile import ConfSettingReader
>>> tn = '20160516_114916' # example tracking number
>>> file_path = getConf4DSettingsPath(tn)
>>> cr = ConfSettingReader(file_path)

Then call methods on 'cr'.
"""
import configparser

class ConfSettingReader():
    """
    Class which reads an interferometer configuration settings file '4DSettings.ini'
    
    Methods
    -------
    getFrameRate() : 
        Gets the camera frame rate in Hz.
    
    getImageWidthInPixels() : 
        Get the width of the frame in pixel units.
    
    getImageHeightInPixels() : 
        Get the height of the frame in pixel units.
    
    getOffsetX() : 
        Get the frame offset in x-axis.
    
    getOffsetY() : 
        Get the frame offset in y-axis.
    
    getPixelFormat() : 
        Get the format of the pixels.
    
    getUserSettingFilePath() : 
        Get the path of the configuration file.
    
    How to Use it
    -------------
    After initializing the class with a file path, just call methods on the defined
    object
    
    >>> cr = ConfSettingReader(file_path)
    >>> cr.getImageWidhtInPixels()
    2000
    >>> cr.getImageHeightInPixels()
    2000
    
    Notes
    -----
    Note that there is no need to directly use this module, as the settings information
    retrievement is handled by m4.urils.osutils, with its functions
    ''getConf4DSettingsPath'' and ''getCameraSettings''.
    """
    def __init__(self, file_path):
        self.config = configparser.ConfigParser()
        self.config.read(file_path)
        self.camera_section = 'ACA2440'
        self.path_section = 'Paths'
    #CAMERA
    def getFrameRate(self):
        """
        Returns the acquisition frame rate of the interferometer in Hz

        Returns
        -------
        frame_rate : float
            The frame rate.
        """
        frame_rate = self.config.get(self.camera_section, 'FrameRate')
        return float(frame_rate)

    def getImageWidhtInPixels(self):
        """
        Returns the image widht in pixel scale

        Returns
        -------
        image_wight_in_pixels : int
            Image pixel width.
        """
        image_width_in_pixels = self.config.get(self.camera_section, 'ImageWidthInPixels')
        return int(image_width_in_pixels)

    def getImageHeightInPixels(self):
        """
       Returns the image height in pixel scale

       Returns
       -------
       image_height_in_pixels : int
           Image pixel height.
        """
        image_height_in_pixels = self.config.get(self.camera_section, 'ImageHeightInPixels')
        return int(image_height_in_pixels)

    def getOffsetX(self):
        """
        Returns the camera offset, in pixels, along the x-axis.

        Returns
        -------
        offset_x : int
            Pixel offset in the x-axis.
        """
        offset_x = self.config.get(self.camera_section, 'OffsetX')
        return int(offset_x)

    def getOffsetY(self):
        """
        Returns the camera offset, in pixels, along the y-axis.

        Returns
        -------
        offset_y : int
            Pixel offset in the y-axis.
        """
        offset_y = self.config.get(self.camera_section, 'OffsetY')
        return int(offset_y)

    def getPixelFormat(self):
        """
        Returns the format of the pixel.

        Returns
        -------
        pixel_format : str
            Pixel format.
        """
        pixel_format = self.config.get(self.camera_section, 'PixelFormat')
        return pixel_format
    #PATH
    def getUserSettingFilePath(self):
        """
        Returns the complete filepath of the settings configuration file.

        Returns
        -------
        user_setting_file_path : str
            Settings file path.
        """
        user_setting_file_path = self.config.get(self.path_section, 'UserSettingsFilePath')
        return user_setting_file_path
    
"""
Author(s) 
---------
    - P. Ferraiuolo : written in 2024

Description
-----------
Script which updates the folder tree on the base of the configuration file, acc
essed through the environment variable 'PYOTTCONF'.

How to Use it
-------------
Simply import the module, and the ''m4.configuration.config_folder_names'' will
be populated.

    >>> from m4.configuration import update_folder_paths as ufp
    >>> fn = ufp.folders
    >>> fn.BASE_PATH
    '.../M4/m4/data'
"""
import os
from m4.configuration.config_reader import configuration_path
from m4.configuration.config_uploader import config_rewriter
try:
    config_file_name = os.environ['PYOTTCONF']
except KeyError:
    raise KeyError("Environment variable not found! Please set up your PYOTTCONF env variable that points the configuration file")
folders        = configuration_path(config_file_name)
cr             = config_rewriter(folders)
cr.upload()
print(f"Folder tree updated.\nBase path is '{folders.BASE_PATH}'")
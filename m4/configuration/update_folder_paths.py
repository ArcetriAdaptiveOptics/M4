"""
Author(s) 
---------
    - P. Ferraiuolo : written in 2024

Description
-----------
Script which updates the folder tree on the base of the configuration file, acc
essed through the environment variable 'PYOTTCONF'.
"""
import os
from m4.configuration.config_reader import configuration_path
from m4.configuration.config_uploader import config_rewriter
try:
    config_file_name = os.environ['PYOTTCONF']
except KeyError:
    raise KeyError("Environment variable not found! Please set up your PYOTTCONF env variable that points the configuration file")
conf_obj         = configuration_path(config_file_name)
cr               = config_rewriter(conf_obj)
cr.upload()
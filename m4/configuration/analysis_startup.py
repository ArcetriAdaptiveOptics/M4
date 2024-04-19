import os
from m4.configuration import config_folder_names as foldname
from m4.configuration.config_reader import configuration_path
from m4.configuration.config_uploader import config_rewriter
#config_file_name = "/home/labot/git/M4/m4/configuration/labotConfig.yaml"

def start(config_file_name=None):

    if config_file_name is None:
        try:
            config_file_name= os.environ["PYOTTCONF"]
        except KeyError:
            raise ValueError("Environment variable PYOTTCONF is not set")
        if ~os.path.exists(config_file_name):
            raise ValueError("Environment variable PYOTTCONF is not set")
    
    conf_obj         = configuration_path(config_file_name)
    cr               = config_rewriter(conf_obj)
    cr.upload()
    print(r"Using the Python console for data analysis.")
    print(r"Base data path is: "+foldname.BASE_PATH)
    print(r"")
    print(r"           |X|_____ _____|X|")
    print(r"           |X|           |X|")
    print(r"           |X|_____--___<|X|")
    print(r"           |X|    |\  /| |X|")
    print(r"     __    |X|    |_\/_| |X|")
    print(r"    |  |   |X|    | /\ | |X|")
    print(r"    |  |   |X|    |/  \| |X|")
    print(r" ___|__|___|X|____ ---- _|X|_____")
    return

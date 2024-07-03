from m4.configuration.config_reader import configuration_path
from m4.configuration.config_uploader import config_rewriter

def create_conf_paths(config_file_name):
    
    conf_obj         = configuration_path(config_file_name)
    cr               = config_rewriter(conf_obj)
    cr.upload()
    from m4.configuration import config_folder_names as foldname
    # print('Using the Python console for data analysis.')
    # print('Base data path is: '+foldname.BASE_PATH)
    # print('')
    # print('           |X|_____ _____|X|')
    # print('           |X|           |X|')
    # print('           |X|_____--___<|X|')
    # print('           |X|    |\  /| |X|')
    # print('     __    |X|    |_\/_| |X|')
    # print('    |  |   |X|    | /\ | |X|')
    # print('    |  |   |X|    |/  \| |X|')
    # print(' ___|__|___|X|____ ---- _|X|_____')
    return

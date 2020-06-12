'''
'''

    #per l'applicazione
def set_up_logger(file_path, logging_level):
    """ Set up logger for the application

    Parameters
    ----------
    file_path: string
            file path and log name (txt) where to create the log
    logging_level: int
                    Warning = 30, Info = 20, Debug = 10, Notset = 0
    
    """
    import logging
    import logging.handlers
    FORMAT = '%(asctime)s %(levelname)s %(name)s %(message)s'
    f = logging.Formatter(fmt = FORMAT)
    handler = logging.handlers.RotatingFileHandler(file_path,
                                                    encoding='utf8',
                                                    maxBytes=10000,
                                                    backupCount=3)
    root_logger = logging.getLogger()
    root_logger.setLevel(logging_level)
    handler.setFormatter(f)
    handler.setLevel(logging_level)
    root_logger.addHandler(handler)
    handler.doRollover()

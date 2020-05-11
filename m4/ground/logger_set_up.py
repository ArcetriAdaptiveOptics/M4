'''
'''

    #per l'applicazione
def set_up_logger(file_path, logging_level):
    """ Set up logger for the application """
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

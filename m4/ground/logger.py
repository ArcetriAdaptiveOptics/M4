'''
@author: cs
'''
import os
import logging
from m4.ground.configuration import Configuration

def set_up_logger(file_path, logging_level):
    import logging.handlers
    FORMAT = '%(asctime)s %(levelname)s %(name)s %(message)s'
    f = logging.Formatter(fmt = FORMAT)
    handler = logging.handlers.RotatingFileHandler(file_path, encoding='utf8',
                                                   maxBytes=1000, backupCount=3)
    root_logger = logging.getLogger()
    root_logger.setLevel(logging_level)
    handler.setFormatter(f)
    handler.setLevel(logging_level)
    root_logger.addHandler(handler)
    handler.doRollover()


def pippo(a):
    logger= logging.getLogger('IlLoggerDiPippo')
    logger.debug('debug mi hai passato %d' % a)
    if a < 0:
        logger.error('error a<0 (a=%g)'% a)



class Pluto():

    def __init__(self):
        self._logger= logging.getLogger('PLUTO')

    def funz1(self, a):
        self._logger.debug('debug mi hai passato %d' % a)
        if a < 0:
            self._logger.error('error a<0 (a=%g)'% a)


def log(info1, info2, info3, info4= None ):
    if info4 is None:
        string=" ".join([info1, info2, info3])
    else:
        string=" ".join([info1, info2, info3, info4])
    dove= Configuration.LOG_ROOT_FOLDER
    filename= os.path.join(dove, 'example.log')

    logging.basicConfig(filename= filename,
                        format='%(asctime)s - %(message)s')
    logging.warning(string)
        
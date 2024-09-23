"""
Author(s)
---------
    - Chiara Selmi : written in 2020
    - Pietro Ferraiuolo : modified in 2024

Description
-----------
Sets up the logger for the application.
"""
import logging
import logging.handlers

def set_up_logger(file_path, logging_level):
    """
    Set up a rotating file logger.
    This function configures a logger to write log messages to a file with
    rotation. The log file will be encoded in UTF-8 and will rotate when it
    reaches a specified size, keeping a specified number of backup files.

    Parameters
    ----------
    file_path : str
        The path to the log file where log messages will be written.
    logging_level : int
        The logging level to set for the logger. This should be one of the
        logging level constants defined in the `logging` module:
            Warning = 30, Info = 20, Debug = 10, Notset = 0

    Notes
    -----
    - The log file will rotate when it reaches 10,000,000 bytes (10 MB).
    - Up to 3 backup log files will be kept.
    - The log format includes the timestamp, log level, logger name, and message.
    - The logger is configured at the root level, affecting all loggers in the application.
    - The handler will perform an initial rollover when set up.

    Examples
    --------
    >>> set_up_logger('/path/to/logfile.log', logging.DEBUG)
    """
    FORMAT = "%(asctime)s %(levelname)s %(name)s %(message)s"
    formato = logging.Formatter(fmt=FORMAT)
    handler = logging.handlers.RotatingFileHandler(
        file_path, encoding="utf8", maxBytes=10000000, backupCount=3
    )
    root_logger = logging.getLogger()
    root_logger.setLevel(logging_level)
    handler.setFormatter(formato)
    handler.setLevel(logging_level)
    root_logger.addHandler(handler)
    handler.doRollover()


def log(message, level: str = "INFO"):
    """
    Log a message at the specified level.

    Parameters
    ----------
    message : str
        The message to log.
    level : str, optional
        The logging level to use for the message. This should be one of the
        following strings: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'. (can
        use lowercase too).
        The default is 'DEBUG'.

    Notes
    -----
    - The message will be logged using the logger configured by `set_up_logger`.
    - The message will be logged with the specified level.
    - If the specified level is not recognized, the message will be logged at the
      'DEBUG' level.
    """
    level = level.upper()
    if level == "DEBUG":
        logging.debug(message)
    elif level == "INFO":
        logging.info(message)
    elif level == "WARNING":
        logging.warning(message)
    elif level == "ERROR":
        logging.error(message)
    elif level == "CRITICAL":
        logging.critical(message)
    else:
        logging.debug(message)
        logging.warning(f"Invalid log level '{level}'. Defaulting to 'DEBUG'.")

class txtLogger:
    """
    A simple logger class that writes log messages to a text file.

    Attributes:
        file_path (str): The path to the log file, name included.

    Methods:
        __init__(file_path):
            Initializes the logger with the specified file path.
        
        log(message):
            Writes a log message to the file.
    """

    def __init__(self, file_path):
        """
        Initializes the txtLogger with the specified file path.

        Args:
            file_path (str): The path to the log file.
        """
        self.file_path = file_path

    def log(self, message):
        """
        Writes a log message to the file.

        Args:
            message (str): The log message to be written to the file.
        """
        with open(self.file_path, 'a') as f:
            f.write(message + '\n')
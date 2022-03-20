import logging
import os
import time


def setup_logging(config, prefix):
    # Setup logging
    log_filename = f"{prefix}_{time.strftime('%d%m%Y-%H%M%S')}.log"

    log_dir = os.path.join(config['data_dir'], "logs")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    logging.basicConfig(level=logging.ERROR,
                        format='%(asctime)s %(levelname)s %(message)s',
                        filename=os.path.join(log_dir, log_filename),
                        filemode='w')

    logging.info("Config: " + str(config))


def plog(message, verbose=True):
    # Print and log the message
    if verbose:
        print(message)
    logging.debug(message)


def elog(message, verbose=True):
    # Log an error
    logging.error(message)
    print(message)


def wlog(message, verbose=True):
    # Print and log warning
    if verbose:
        print(message)
    logging.warning(message)

import logging
from copy import copy
from threading import Thread
from queue import Empty
from time import sleep
from traceback import print_exc
from sys import stderr, exit as sys_exit
from multiprocessing import Queue, Lock
from logging.handlers import QueueHandler


__copyright__ = 'Copyright (C) 2017-2018 ICTP'
__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__credits__ = ["Stefano Piani"]


# The Queue where the processes will put their logs to be printed by the
# main process. It must be initialized  with the function "collect_logs()"
LOG_QUEUE = None


# See bug https://bugs.python.org/issue6721 to understand the reason of this
# Lock. Essentially, we need to avoid to log something and, at the very same
# moment, to execute a Fork because this would replicate the Lock of the log
# and create a deadlock. It must be initialized  with the function
# "collect_logs()"
FORK_LOCK = None


class OnelinerFormatter(logging.Formatter):
    """
    Alter the text of the logs replacing "\n" with " | ".
    """
    def __init__(self, fmt=None, datefmt=None):
        super(OnelinerFormatter, self).__init__(fmt, datefmt)

    def format(self, record):
        new_record = copy(record)
        new_record.msg = new_record.msg.strip('\n').replace('\n', ' | ')
        return super(OnelinerFormatter, self).format(new_record)


class LogCollector(Thread):
    """
    A LogCollector is a thread that checks every "wait" seconds for
    log records inside the LOGS_QUEUE. If records are found, they
    are handled by the appropriate logger inside this thread.
    """

    def __init__(self, logs_queue=LOG_QUEUE, lock=FORK_LOCK, wait=0.01):
        super(LogCollector, self).__init__()
        self.daemon = True
        self.wait = wait
        self.logs_queue = logs_queue
        self._lock = lock

    def run(self):
        while True:
            try:
                record = self.logs_queue.get(block=False)
                logger = logging.getLogger(record.name)
                with self._lock:
                    logger.handle(record)

            except Empty:
                sleep(self.wait)
            except Exception:
                print(
                    'WARNING - Logging service anomaly terminated!',
                    file=stderr
                )
                print_exc(file=stderr)
                sys_exit(1)


class WaitQueueHandler(QueueHandler):
    def enqueue(self, record):
        self.queue.put(record, block=True, timeout=None)


def collect_logs():
    global LOG_QUEUE
    LOG_QUEUE = Queue(1000)

    global FORK_LOCK
    FORK_LOCK = Lock()

    log_collector = LogCollector(LOG_QUEUE, FORK_LOCK)
    log_collector.start()

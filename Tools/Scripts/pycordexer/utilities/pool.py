from multiprocessing import Process
from time import sleep
from threading import Thread, Lock as TLock
from logging import getLogger
import utilities.log_utilities


__copyright__ = 'Copyright (C) 2017-2018 ICTP'
__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__credits__ = ["Stefano Piani"]


LOGGER = getLogger(__name__)


class Pool(object):
    """
    A Pool is an object similar to multiprocessing.Pool.

    The main differences are that Pool spawns a new process for each function
    that must be computed and that, if there are already max_num_or_process
    processes spawned, trying to spawn a new process will block until one of the
    other processes terminates.

    The main reason to use this object is that you may want the new processes
    spawned to inherit some global variables from the main process that could
    have changed their values during the run.

    In this case, you can pass these variables in a list as obj_list. This
    ensures that the garbage collector does not remove them until the execution
    ends
    """
    def __init__(self, max_num_of_process, name='process'):
        self.max_num_of_process = max(1, max_num_of_process)
        self._spawned = []
        self._name = name
        self._total_spawned = 0
        self._stop = False

        # _objects is a dictionary that links some global objects to a list of
        # processes that depends on that object. Therefore, this ensures that
        # the garbage collector does not remove such objects until the processes
        # that are using that processes end their execution
        self._objects = {}
        self._objects_lock = TLock()

        self._process_cleaner = ProcessCleaner(
            self._spawned,
            self._objects,
            self._objects_lock,
        )
        self._process_cleaner.start()

    def stop_now(self):
        """
        Inhibit the possibility of this pool to launch new processes. Moreover,
        if another thread is waiting for executing the execute method, it will
        return immediately without launching any process (the method will return
        None)
        """
        self._stop = True

    def set_max_number_of_processes(self, max_num):
        max_num = max(1, max_num)
        LOGGER.debug(
            'Changing max time of processes from %s to %s',
            self.max_num_of_process,
            max_num
        )
        self.max_num_of_process = max_num

    def execute(self, call_obj, args, obj_list=None):
        """
        Execute a callable_object in a different process. If you have reached
        the maximum number of processes, it will wait until one of them exit

        :param call_obj: a callable object that must be executed by a new
               process
        :param args: the arguments that must be passed to call_obj
        :param obj_list: A list of global objects that the process will use

        :return: The process that is executing the callable_object. If this
        method is interrupted by someone calling the stop_now method, it will
        return None
        """
        if obj_list is None:
            obj_list = []

        while len(self._spawned) >= self.max_num_of_process and not self._stop:
            LOGGER.log(
                1,
                'There are %s processes running (%s is the max). Waiting in '
                'queue',
                len(self._spawned),
                self.max_num_of_process
            )
            sleep(0.25)

        if self._stop:
            return None

        self._total_spawned += 1
        process_name = '{} {}'.format(self._name, self._total_spawned)
        args = args[:] + (utilities.log_utilities.LOG_QUEUE,)
        with utilities.log_utilities.FORK_LOCK:
            LOGGER.debug('Spawning process %s', process_name)
            p = Process(target=call_obj, args=args, name=process_name)
            LOGGER.debug('Starting process %s', process_name)
            p.start()
            self._spawned.append(p)

        with self._objects_lock:
            for obj in obj_list:
                if obj in self._objects:
                    self._objects[obj].append(p)
                else:
                    self._objects[obj] = [p, ]
        return p

    def close(self):
        LOGGER.debug('Waiting all the processes to terminate')
        while len(self._spawned) > 0:
            sleep(0.5)
        LOGGER.debug('All processes terminated')
        return True


class ProcessCleaner(Thread):
    """
    A ProcessCleaner is a thread that takes a list of processes and continuously
    checks if a process is not running anymore. In that case, it removes it from
    the list
    """
    def __init__(self, process_list, object_dict, object_dict_lock):
        self._process_list = process_list
        self._objects = object_dict
        self._objects_lock = object_dict_lock

        super(ProcessCleaner, self).__init__()
        self.daemon = True

    def run(self):
        while True:
            for p in self._process_list[:]:
                if not p.is_alive():
                    with utilities.log_utilities.FORK_LOCK:
                        LOGGER.debug(
                            'Process %s is stopped. Removed from the list',
                            p.name
                        )
                    self.remove_p_obj(p)
                    self._process_list.remove(p)
            sleep(0.5)

    def remove_p_obj(self, p):
        with self._objects_lock:
            to_be_removed = []
            for obj in self._objects:
                if p in self._objects[obj]:
                    with utilities.log_utilities.FORK_LOCK:
                        LOGGER.debug(
                            'Removing %s from dependencies of %s',
                            p.name,
                            obj
                        )
                    self._objects[obj].remove(p)
                    if len(self._objects[obj]) == 0:
                        to_be_removed.append(obj)

            for obj in to_be_removed:
                with utilities.log_utilities.FORK_LOCK:
                    LOGGER.debug(
                        'No more processes running that refer to %s',
                        obj
                    )
                del self._objects[obj]

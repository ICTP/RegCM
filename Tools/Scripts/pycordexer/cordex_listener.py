#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Alberto Chiusole, alberto.chiusole@exact-lab.it
# Date: 2017-09-29


# This script provides a "listener" on a directory that executes the
# command requested by the user when two or more files are found.
# It is implemented using a busy-waiting polling technique with a long
# delay.
#
# When a file named "COMPLETED" is found in the directory, it means that
# all the other files have been written on disk and can be processed:
# after they have been processed, the script is terminated.
#
# The script can also be executed manually by using the `--manual-mode`
# option, which does not use the "listener" feature and searches the
# given directory for files matching the patterns provided, executing the
# command passed as input.
#
# A namelist file argument is required, in order to import
# configuration parameters (the python package f90nml is required).

import argparse
import os
from os.path import join
from traceback import format_exc
import time
from datetime import datetime
from sys import exit as sys_exit, argv
from threading import Thread
import glob
import f90nml
import logging
import signal
import re
from multiprocessing import cpu_count


# The following module is used to extract the current timezone
# If it is missing, the date will be saved without timezone
try:
    from tzlocal import get_localzone
except ImportError:
    get_localzone = None
    tzlocal_import_error = format_exc()

# These modules are used only in daemon mode
try:
    import daemon
    import daemon.pidfile
except ImportError:
    daemon = None
    daemon_import_error = format_exc()

from pycordexer import save_vars as save_cordex_vars
import utilities.log_utilities


__copyright__ = 'Copyright (C) 2017-2018 ICTP'
__author__ = 'Alberto Chiusole <alberto.chiusole@exact-lab.it>'
__credits__ = ["Alberto Chiusole", "Stefano Piani"]


_GLOB_MASKS = ["*_ATM.??????????.nc", "*_SRF.??????????.nc",
               "*_STS.??????????.nc", "*_RAD.??????????.nc",
               "*_SHF.??????????.nc", "*.??????????.txt"]
_SLEEP_TIMER = 10

# A flag to stop the execution of this software as soon as possible
EXIT_NOW = False

# A flag to stop the execution of this software when all the data have been
# processed
EXIT_WHEN_DONE = False

# The string that will be appended to the txt files to show that the data have
# been processed
APPEND_TXT_STRING = 'Elaborated by Pycordexer at %(now)'

# The string that will be appended to the txt files to show that the elaboration
# failed
APPEND_FAILED_TXT_STRING = 'Failed elaboration attempt by Pycordexer at %(now)'

# The set of the files that are already been elaborated by this software
ALREADY_ELABORATED = set()

# The pools of processes that must elaborate the data
WORKERS_POOL = None

# A list with all the threads that are running to delete the original files
# produced by RegCM
FILE_DELETERS = []

# A list with all the threads that are running to manage the txt files produced
# by RegCM
TXT_HANDLERS = []

# Save the default handlers for the SIGINT and SIGTERM
DEFAULT_SIGINT_HANDLER = signal.getsignal(signal.SIGINT)
DEFAULT_SIGTERM_HANDLER = signal.getsignal(signal.SIGTERM)


LOGGER = logging.getLogger()


class MissingInputDir(ValueError):
    pass


class FileDeleter(Thread):
    """
    A FileDeleter deletes a particular file as soon as the processes that are
    elaborating it stop their execution
    """

    def __init__(self, file_path, process_list):
        super(FileDeleter, self).__init__()
        self.file_path = file_path
        self.process_list = list(process_list)
        self.daemon = False

    def run_file_deleter(self):
        if None in self.process_list:
            with utilities.log_utilities.FORK_LOCK:
                LOGGER.debug(
                    'File deleter for %s exited without deleting the file '
                    'because the computation was interrupted before it was '
                    'completely elaborated', self.file_path
                )
            return
        for p in self.process_list:
            p.join()
        with utilities.log_utilities.FORK_LOCK:
            LOGGER.debug(
                'All processes that were elaborating the file %s have '
                'completed their execution', self.file_path
            )

        exit_status = [p.exitcode for p in self.process_list]

        all_ok = True
        for exit_sts in exit_status:
            if exit_sts != 0:
                all_ok = False

        if all_ok:
            with utilities.log_utilities.FORK_LOCK:
                LOGGER.debug('Deleting file %s', self.file_path)
                try:
                    os.remove(self.file_path)
                    LOGGER.info('File %s removed', self.file_path)
                except Exception:
                    LOGGER.warning(
                        'Error deleting file {}:\n{}'.format(self.file_path,
                                                             format_exc())
                    )

        else:
            with utilities.log_utilities.FORK_LOCK:
                LOGGER.warning(
                    'File %s will not be eliminated because some errors '
                    'occurred while processing the variables it contains',
                    self.file_path
                )

    def run(self):
        return self.run_file_deleter()


class TxtHandler(Thread):
    """
    A TxtHandler takes care of updating a txt file as soon as all the processes
    that are elaborating data it refers to end their execution
    """

    def __init__(self, txt_path, process_list, must_be_deleted=False):
        super(TxtHandler, self).__init__()
        self.txt_path = txt_path
        self.process_list = list(process_list)
        self.delete = must_be_deleted
        self.daemon = False

    def run_txt_handler(self):
        if None in self.process_list:
            with utilities.log_utilities.FORK_LOCK:
                LOGGER.debug(
                    'Txt handler for %s exited without editing the file '
                    'because the computation was interrupted before the data '
                    'it refers to was completely elaborated', self.txt_path
                )
            return
        for p in self.process_list:
            p.join()
        with utilities.log_utilities.FORK_LOCK:
            LOGGER.debug(
                'All processes that were elaborating data which the file %s '
                'refers to have terminated their execution', self.txt_path
            )

        exit_status = [p.exitcode for p in self.process_list]

        all_ok = True
        for exit_sts in exit_status:
            if exit_sts != 0:
                all_ok = False

        if self.delete and all_ok:
            deleted = False
            try:
                os.remove(self.txt_path)
                deleted = True
            except Exception:
                with utilities.log_utilities.FORK_LOCK:
                    LOGGER.error(
                        'Error deleting file {}:\n{}'.format(
                            self.txt_path,
                            format_exc()
                        )
                    )
            if deleted:
                with utilities.log_utilities.FORK_LOCK:
                    LOGGER.info('File %s removed', self.txt_path)
        else:
            append_string = _get_append_string(all_ok)
            with utilities.log_utilities.FORK_LOCK:
                LOGGER.debug(
                    'Appending string {} to file {}'.format(
                        append_string,
                        self.txt_path
                    )
                )
            try:
                with open(self.txt_path, 'a', encoding='ASCII') as f:
                    f.write(append_string)
            except Exception:
                with utilities.log_utilities.FORK_LOCK:
                    LOGGER.error(
                        'Error appending string "{}" to {}:\n{}'.format(
                            append_string,
                            self.txt_path,
                            format_exc()
                        )
                    )

    def run(self):
        return self.run_txt_handler()


def _clean_file_deleters():
    to_be_removed = []
    for fd in FILE_DELETERS:
        if not fd.is_alive():
            to_be_removed.append(fd)
    for fd in to_be_removed:
        FILE_DELETERS.remove(fd)


def _clean_txt_handler():
    to_be_removed = []
    for fd in TXT_HANDLERS:
        if not fd.is_alive():
            to_be_removed.append(fd)
    for fd in to_be_removed:
        TXT_HANDLERS.remove(fd)


def _get_append_string(success=True):
    """
    Produce the string that will be appended to the txt files
    """
    now_naive = datetime.now()
    if get_localzone is not None:
        now = get_localzone().localize(now_naive)
    else:
        now = now_naive

    now = str(now)

    if success:
        return '\n' + APPEND_TXT_STRING.replace('%(now)', now) + '\n'
    else:
        return '\n' + APPEND_FAILED_TXT_STRING.replace('%(now)', now) + '\n'


def _collect_files_by_masks(base_directory):
    collected_files = []

    LOGGER.debug('Examining dir %s', base_directory)
    for mask in _GLOB_MASKS:
        collected_files.extend(glob.glob(os.path.join(base_directory, mask)))

    return collected_files


def _file_name_2_check_mask(file_name):
    """
    Starting from a netcdf output file of RegCM (like ATM or SRF), produce the
    name of the txt file that is produced when the other files are closed and
    ready to be processed
    """
    file_name_split = file_name.split('.')
    file_date = file_name_split[-2]
    file_name_main = '.'.join(file_name_split[:-2])[:-4]
    return file_name_main + '.' + file_date + '.txt'


def _check_txt_file_is_elaborated(file_path):
    """
    Check if the txt_file contains a line that could have been appended by
    this script. Return True if such a line is found
    """
    append_string_match = '^' + APPEND_TXT_STRING.replace('%(now)', '.*') + '$'
    LOGGER.debug(
        'Looking for a line that matches this regular expression: "%s"',
        append_string_match
    )
    try:
        with open(file_path, 'r') as f:
            lines = [l.strip('\n') for l in f.readlines()]
    except Exception:
        LOGGER.warning(
            'Error opening file {}:\n{}'.format(file_path, format_exc())
        )
        raise

    for line in lines:
        if re.match(append_string_match, line) is not None:
            LOGGER.debug(
                'Line "%s" honours regexp "%s", file elaborated!',
                line,
                append_string_match,
            )
            return True

    LOGGER.debug(
        'No line found that honours "%s". File not elaborated',
        append_string_match
    )
    return False


def _process_regcm_output_file(file_path, pycordexer_opts, remove_file, pool):
    """
    Process a single RegCM output file. Return a list of all the process spawned
    to perfor the elaboration
    """
    if EXIT_NOW:
        LOGGER.debug(
            'Conversion of file %s will be skipped because this script must '
            'terminate', file_path
        )
        return [None]

    pycordexer_opts_copy = dict(pycordexer_opts)
    pycordexer_opts_copy['sigterm_handler'] = DEFAULT_SIGTERM_HANDLER
    pycordexer_opts_copy['sigint_handler'] = DEFAULT_SIGINT_HANDLER

    LOGGER.info('Converting file %s...', file_path)
    try:
        procs = save_cordex_vars(file_path, 'ALL', pool, **pycordexer_opts_copy)
        LOGGER.debug('All process for conversion of file %s spawned', file_path)
    except Exception:
        LOGGER.warning(
            'Error converting file {}:\n{}'.format(file_path, format_exc())
        )
        return [None]

    if remove_file:
        file_deleter = FileDeleter(file_path, procs)
        file_deleter.start()
        FILE_DELETERS.append(file_deleter)

    return procs


def _process_file_group_in_manual_mode(group_of_files, pycordexer_options,
                                       remove_files, pool):
    """
    A group of files is a list of RegCM output files. All the entries but the
    last are NetCDF files that refers to the same period of simulation.
    The last entry could be the txt file produced by RegCM related of the same
    period of simulation or None (if that file is missing)
    This function will elaborate all the netcdf files using pycordex and, if
    the remove_files flag is set, will also remove them.
    Afterwards, if the last entry of the group is a txt fie and the
    remove_files flag is set, it will be removed. If the flag is not set, a
    string containing the elaboration time will be appended to the txt file.
    """
    regcm_netcdf_files = group_of_files[:-1]
    txt_file = group_of_files[-1]

    spawned_processes = []

    log_message = 'Elaborating a group with the following files:\n'
    log_message += '\n'.join(['    ' + i for i in regcm_netcdf_files])
    if txt_file is None:
        log_message += '\nand no txt file'
    else:
        log_message += '\nand this txt file:\n    {}'.format(txt_file)
    LOGGER.debug(log_message)

    for fl in regcm_netcdf_files:
        prc = _process_regcm_output_file(
            fl,
            pycordexer_options,
            remove_files,
            pool
        )
        spawned_processes.extend(prc)

    if txt_file is not None:
        LOGGER.debug('Launching a txt_handler for %s', txt_file)
        txt_handler = TxtHandler(txt_file, spawned_processes, remove_files)
        txt_handler.start()
        TXT_HANDLERS.append(txt_handler)


def _process_file_group(group_of_files, pycordexer_options, remove_files, pool):
    """
    A group of files is a list of RegCM output files. All the entries but the
    last are NetCDF files that refers to the same period of simulation.
    The last entry could be the txt file produced by RegCM related of the same
    period of simulation or None (if that file is missing)
    """
    regcm_netcdf_files = group_of_files[:-1]
    txt_file = group_of_files[-1]

    log_message = 'Elaborating a group with the following files:\n'
    log_message += '\n'.join(['    ' + i for i in regcm_netcdf_files])
    if txt_file is None:
        log_message += '\nand no txt file'
    else:
        log_message += '\nand this txt file:\n    {}'.format(txt_file)
    LOGGER.debug(log_message)

    if txt_file is None:
        LOGGER.debug('Discarding the current group because it has not txt file')
        return

    if txt_file in ALREADY_ELABORATED:
        LOGGER.debug(
            'Discarding the current group because it has already '
            'been elaborated'
        )
        return

    # We need to check if the current txt file has already been elaborated by
    # a previous process
    LOGGER.debug(
        'Reading file "%s" to check if its group was already being elaborated',
        txt_file
    )
    try:
        txt_file_elaborated = _check_txt_file_is_elaborated(txt_file)
    except:
        LOGGER.warning(
            'Error reading file "%s". Its group will be elaborated anyway',
            txt_file
        )
        txt_file_elaborated = False

    if txt_file_elaborated:
        LOGGER.debug(
            'Accordling to the content of the file "%s", its group has '
            'already been elaborated. Stopping elaboration of this group',
            txt_file,
        )
        ALREADY_ELABORATED.add(txt_file)
        return

    spawned_processes = []

    LOGGER.debug(
        'It seems that this group has not been elaborated before. Proceding...'
    )
    for fl in regcm_netcdf_files:
        prc = _process_regcm_output_file(
            fl,
            pycordexer_options,
            remove_files,
            pool
        )
        spawned_processes.extend(prc)

    ALREADY_ELABORATED.add(txt_file)
    if txt_file is not None:
        LOGGER.debug('Launching a txt_handler for %s', txt_file)
        txt_handler = TxtHandler(txt_file, spawned_processes, remove_files)
        txt_handler.start()
        TXT_HANDLERS.append(txt_handler)


def _process_file_list(list_of_files, pycordexer_options, manual_mode,
                       remove_files, pool):
    """
    Given a list of files to process, execute pycordexer on all of them
    starting from the oldest and move to the newer.

    If not in manual mode, the files will not be processed until a txt file
    with the same name appears (to be sure that RegCM has closed the files).

    If remove file is True, the files will be deleted after the extraction
    """

    # Sort the list of files from the oldest to the most recent.
    list_of_files = sorted(list_of_files, key=os.path.getctime)[::-1]
    LOGGER.debug('The following files will be processed: %s', list_of_files)

    regcm_output_files = [f for f in list_of_files if f.lower().endswith('.nc')]
    regcm_txt_files = [f for f in list_of_files if f.lower().endswith('.txt')]

    txt_associated = {}
    for f in regcm_output_files:
        txt_associated[f] = _file_name_2_check_mask(f)
        LOGGER.debug(
            'The txt file related to the NetCDF file %s should be %s',
            f,
            txt_associated[f]
        )

    # Pack the different files together in groups with the same txt file
    LOGGER.debug('Associating files in groups')
    file_groups = []
    while len(regcm_output_files):
        grp_txt = txt_associated[regcm_output_files[0]]
        group = [f for f in regcm_output_files if txt_associated[f] == grp_txt]
        for f in group:
            regcm_output_files.remove(f)
        if grp_txt in regcm_txt_files:
            group.append(grp_txt)
        else:
            group.append(None)
        file_groups.append(group)

    for group in file_groups:
        if manual_mode:
            _process_file_group_in_manual_mode(
                group,
                pycordexer_options,
                remove_files,
                pool
            )
        else:
            _process_file_group(
                group,
                pycordexer_options,
                remove_files,
                pool
            )
        if EXIT_NOW:
            break


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("{}: no such file or directory.".format(arg))
    return arg


def is_valid_dir(parser, arg):
    if not os.path.isdir(arg):
        parser.error("{}: no such file or directory.".format(arg))
    return arg


def _prepare_logger(verbosity, file_verbosity, silent_mode, log_file):
    """
    Prepare the logging system and return the two handlers (the one that writes
    on the standard output and the one that writes on the logfile).

    For each handler that has not been created, return None as a placeholder

    """
    verbosity_level = getattr(logging, verbosity.upper())
    file_verbosity_level = getattr(logging, file_verbosity.upper())
    LOGGER.setLevel(min(verbosity_level, file_verbosity_level))

    if not silent_mode:
        streamformatter = logging.Formatter(
            '%(levelname)s - %(processName)s: %(message)s'
        )
        streamhandler = logging.StreamHandler()
        streamhandler.setLevel(verbosity_level)
        streamhandler.setFormatter(streamformatter)
        LOGGER.addHandler(streamhandler)
    else:
        streamhandler = None

    if log_file is not None:
        fileformatter = utilities.log_utilities.OnelinerFormatter(
            '%(asctime)s - %(levelname)s - %(processName)s - %(funcName)s: '
            '%(message)s',
            datefmt='%Y/%m/%d %H:%M:%S'
        )
        filehandler = logging.FileHandler(log_file)
        filehandler.setLevel(file_verbosity_level)
        filehandler.setFormatter(fileformatter)
        LOGGER.addHandler(filehandler)
    else:
        filehandler = None

    # This avoid problems if no other handler is applied on the LOGGER
    LOGGER.addHandler(logging.NullHandler())

    return streamhandler, filehandler


def _create_input_parser():
    parser = argparse.ArgumentParser()

    v_levs = ['debug', 'info', 'warning']

    parser.add_argument("-f", "--namelist_file", help="Path to a namelist file,"
                        " containing options stored as Fortran instructions",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-d", "--directory", help="Path to a directory to "
                        "keep listening on (overrides the `dirout' parameter "
                        "contained in the namelist file)", default=None,
                        type=lambda x: is_valid_dir(parser, x))
    parser.add_argument("-m", "--manual-mode", help="Does not make the process"
                        " 'listen' on the directory, but processes every file "
                        "inside it and then exit",
                        action="store_true", default=False)
    parser.add_argument("-t", "--timer", help="Time to wait between each check"
                        " on the directory", type=int, default=_SLEEP_TIMER)
    parser.add_argument("-s", "--silent", help="Do not print logs on "
                        "current terminal", action="store_true", default=False)
    parser.add_argument("-r", "--remove-files", help="Automatically remove "
                        "RegCM files after having extracted the variables",
                        action="store_true", default=False)
    parser.add_argument("--daemon", help="Using this flag, this script acts as "
                        "a daemon, working in background", action="store_true",
                        default=False)
    parser.add_argument("-p", "--pidfile", help='Create a file with the pid of '
                        'the daemon created by this script. It is ignored if '
                        '--daemon is not submitted', type=str, default=None)
    parser.add_argument("--processes", type=int, default=1, help='The max '
                        'number of concurrency workers that this script is '
                        'allowed to run')
    parser.add_argument("-l", "--logfile", help="Save the logs on this file",
                        type=str, default=None)
    parser.add_argument('-v', '--verbosity', choices=v_levs, default='info',
                        help='The level of verbosity of this script on screen '
                        '(by default is "info")')
    parser.add_argument('-w', '--file-verbosity', choices=v_levs, default=None,
                        help='The level of verbosity of this script on the log '
                        'file. If this is not specified, it has the same of '
                        '"verbosity"')
    parser.add_argument("--mail", help="The mail of the user that is "
                        "generating the CORDEX files", type=str, default=None)
    parser.add_argument("--global-model", help="The global model of this "
                        "simulation", type=str, default=None)
    parser.add_argument("--experiment", help="The experiment of this "
                        "simulation", type=str, default=None)
    parser.add_argument("--ensemble", help="The ensemble of this "
                        "simulation", type=str, default=None)
    parser.add_argument("--domain", help="The domain of the simulation the "
                        "files belong to", type=str, default=None)
    parser.add_argument("--notes", help="Some notes about the simulation the "
                        "files belong to", type=str, default=None)
    parser.add_argument('--regcm-version', help='The version of the current '
                        'model of RegCM. If this is not specified, this script '
                        'will try to guess the version from the input file. '
                        'Otherwise, submit a string like "5.1" or '
                        '"rev6446". This is the string that will be appended '
                        'after "RegCM-"', type=str, default=None)
    parser.add_argument('--regcm-version-id', help='The version id of the '
                        'current model of RegCM. Usually, this is a number '
                        'like 0 or 1 and will be appended to the output file '
                        'names as "v0" or "v1". This argument is ignored if '
                        '--regcm-version is not submitted and it is mandatory '
                        'otherwise', type=int, default=None)
    parser.add_argument('--regcm_nest_tag', help='For the FPS where the model '
                        'is nested into a previous RegCM run, we can provide '
                        'here the information on the nesting. Usually this is '
                        'a string like fpsconv-x2yn2-v1. Read the experiment '
                        'protocol documentation to find correct format. This '
                        'argument is ignored if --regcm-version is not '
                        'submitted and it is mandatory otherwise',
                        type=str, default=None)
    return parser


def _read_pycordexer_options(args):
    """
    Read the Fortran configuration file specified in the command line
    arguments and return the options needed to run the pycordexer script.

    Returns:
        - *pycordexer_options*: a dictionary with the options for the pycordexer
          script
        - *input_dir*: the path of the directory that must be monitored by this
          script
    """
    pycordexer_options = {}

    if args['namelist_file'] is not None:
        namelist_path = args['namelist_file']
        # Try to interpret the namelist file passed as input and extract a
        # dict-like object
        LOGGER.debug('Reading Fortran configuration file "%s"', namelist_path)
        try:
            args_from_text_file = f90nml.read(namelist_path)
        except:
            LOGGER.error(
                'Error reading file "{}":{}\nPlease, submit a '
                'valid fortran namelist file'.format(
                    namelist_path,
                    format_exc(),
                )
            )
            raise
    else:
        LOGGER.debug(
            'No namelist file specified. Pycordexer options will be '
            'read from command line parameters only'
        )
        namelist_path = None
        args_from_text_file = None

    if args['directory'] is not None:
        LOGGER.debug(
            'Using input dir submitted by command line: "%s"',
            args['directory']
        )
        input_dir = args['directory']
    elif args_from_text_file is None:
        LOGGER.error(
            'Neither "directory" not "namelist file" have been specified in '
            'through the command line interface. At least one among these '
            'parameters must be submitted!'
        )
        raise MissingInputDir('No input dir and no namelist file submitted')
    else:
        outparam_options = args_from_text_file.get("outparam")
        if outparam_options is None:
            LOGGER.error(
                'No section "outparam" found in file "%s". Therefore, monitor '
                'dir can not be read from the namelist file. Please, run '
                'again this script specifying the directory with the -d flag '
                '(or edit your namelist file)',
                namelist_path
            )
            raise MissingInputDir(
                'No section "outparam" found in the namelist file'
            )

        input_dir = outparam_options.get("dirout")
        if input_dir is None:
            LOGGER.error(
                'No field "dirout" found in the outparam section of the '
                'file "%s". Therefore, monitor dir can not be read from '
                'the namelist file. Please, run again this script '
                'specifying the directory with the -d flag (or edit your '
                'namelist file)',
                namelist_path
            )
            raise MissingInputDir(
                '"dirout" parameter not found in section "outparam"'
            )

        LOGGER.debug(
            'Read directory to monitor from Fortran file: "%s"',
            input_dir
        )

    # Add the output directory to the configuration of pycordexer
    pycordexer_options['cordex_root_dir'] = join(input_dir, 'CORDEX')

    if args_from_text_file is None:
        postprocess_options = None
    else:
        postprocess_options = args_from_text_file.get("postprocess")

    if postprocess_options is None:
        if args_from_text_file is not None:
            LOGGER.warning(
                'No "postprocess" section found in the namefile. Pycordexer '
                'will use the default values for the parameter not submitted '
                'through command line'
            )
        else:
            LOGGER.debug(
                'No Fortran file specified: skipping read of "postprocess" '
                'section'
            )
    else:
        LOGGER.debug('"Postprocess" section of file %s read', namelist_path)

    if postprocess_options is None:
        def read_parameter(p_name, p_in_postprocess):
            if args[p_name] is not None:
                LOGGER.debug(
                    'Using the value from the command line for field "%s": %s',
                    p_name,
                    args[p_name],
                )
                pycordexer_options[p_name] = args[p_name]
            else:
                LOGGER.debug(
                    'No value for field "%s" found. Pycordexer will use '
                    'the default one',
                    p_name
                )
    else:
        def read_parameter(p_name, p_in_postprocess):
            if args[p_name] is not None:
                LOGGER.debug(
                    'Using the value from the command line for field "%s": %s',
                    p_name,
                    args[p_name],
                )
                pycordexer_options[p_name] = args[p_name]
            else:
                try:
                    namelist_option = postprocess_options.get(p_in_postprocess)
                    if namelist_option is None:
                        raise ValueError('{} not found in args'.format(p_name))
                    pycordexer_options[p_name] = namelist_option
                    LOGGER.debug(
                        '"%s" field read from the namelist file: %s',
                        p_name,
                        pycordexer_options[p_name]
                    )
                except:
                    LOGGER.debug(
                        'No value for field "%s" found. Pycordexer will use '
                        'the default one',
                        p_name
                    )

    read_parameter('mail', 'user')
    read_parameter('global_model', 'global_model')
    read_parameter('experiment', 'experiment')
    read_parameter('ensemble', 'ensemble')
    read_parameter('notes', 'notes')

    if args['domain'] is not None:
        LOGGER.debug(
            'Using the value from the command line for field "domain": %s',
            args['domain']
        )
        pycordexer_options['domain'] = args['domain']
    else:
        # The domain can be extracted from 'terrainparam' => 'domname'.
        if args_from_text_file is None:
            terrainparam_options = None
        else:
            terrainparam_options = args_from_text_file.get("terrainparam")

        if terrainparam_options is None:
            if args_from_text_file is not None:
                LOGGER.warning(
                    'No "terrainparam" section found in the namefile. '
                    'Pycordexer will use the default value'
                )
            else:
                LOGGER.debug(
                    'No Fortran file specified. Skipping read "terrainparam" '
                    'section'
                )
        else:
            LOGGER.debug('"terraparam" section read')
            domain = terrainparam_options.get("domname", None)
            if domain is None:
                LOGGER.warning(
                    'No value for field "domain" found in the namelist '
                    'file because field "domname" is missing in the '
                    '"terraparam" section. Pycordexer will use the default '
                    'one.'
                )
            else:
                LOGGER.debug(
                    'Domain field read from the namelist file: %s',
                    domain
                )
                pycordexer_options['domain'] = domain

    regcm_version = args['regcm_version']
    if regcm_version is None:
        regcm_version_id = None
        regcm_nest_tag = None
    else:
        regcm_version_id = args['regcm_version_id']
        regcm_nest_tag = args['regcm_nest_tag']
        if regcm_version_id is None and regcm_nest_tag is None:
            raise ValueError(
                'If --regcm-version is submitted, then one in the couple '
                '--regcm-version-id,--regcm_nest_tag is mandatory and must '
                'be specified.'
            )

    if regcm_version is not None:
        pycordexer_options['regcm_version'] = regcm_version
        pycordexer_options['regcm_version_id'] = regcm_version_id
        pycordexer_options['regcm_nest_tag'] = regcm_nest_tag

    return pycordexer_options, input_dir


def stop_execution_on_signal(signum, frame):
    """
    Stop the execution as soon as all the processes that are running right now
    terminate their execution. Avoid to spawn other processes.
    """
    global EXIT_NOW
    global WORKERS_POOL

    LOGGER.info(
        'Received a signal number %s. Execution will be stopped as soon as '
        'all operations terminate',
        signum
    )

    EXIT_NOW = True
    if WORKERS_POOL is not None:
        WORKERS_POOL.stop_now()


def stop_execution_when_done(signum, frame):
    """
    Stop the execution as soon as all the files that have been produced at the
    time being have been elaborated.
    """
    global EXIT_WHEN_DONE

    LOGGER.info(
        'Received a signal number %s. Execution will be stopped when all the '
        'current files will be converted',
        signum
    )

    EXIT_WHEN_DONE = True


def stop_execution_when_done_max_procs(signum, frame):
    """
    Stop the execution as soon as all the files that have been produced at the
    time being have been elaborated. Use the maximum number of processes until
    the execution ends.
    """
    global EXIT_WHEN_DONE

    LOGGER.info(
        'Received a signal number %s. Execution will be stopped when all the '
        'current files will be converted',
        signum
    )

    EXIT_WHEN_DONE = True
    if WORKERS_POOL is not None:
        WORKERS_POOL.set_max_number_of_processes(cpu_count())


class Listener(object):

    def __init__(self, input_dir, pycordexer_options, manual_mode, remove_files,
                 timer=_SLEEP_TIMER):
        self.input_dir = input_dir
        self.pycordexer_options = pycordexer_options
        self.manual_mode = manual_mode
        self.remove_files = remove_files
        self.timer = timer

    def run(self, processes):
        global WORKERS_POOL
        # This import must be here to ensure that the global variables of the
        # module are already initialized
        from utilities.pool import Pool
        WORKERS_POOL = Pool(max_num_of_process=processes, name='Subordinate')
        files_checked_while_exit_when_done = False

        while not EXIT_NOW:
            files_found = _collect_files_by_masks(self.input_dir)

            if EXIT_WHEN_DONE:
                files_checked_while_exit_when_done = True

            if len(files_found) > 0:
                LOGGER.debug('Some files have been found')
                _process_file_list(
                    files_found,
                    self.pycordexer_options,
                    self.manual_mode,
                    self.remove_files,
                    WORKERS_POOL
                )

            if self.manual_mode:
                break
            elif files_checked_while_exit_when_done:
                break
            elif not (EXIT_NOW or EXIT_WHEN_DONE):
                LOGGER.debug('Sleeping for %s seconds', self.timer)
                sleep_time = 0
                timer_step = 0.5
                stop_waiting = EXIT_NOW or EXIT_WHEN_DONE
                while sleep_time < self.timer and not stop_waiting:
                    time.sleep(timer_step)
                    stop_waiting = EXIT_NOW or EXIT_WHEN_DONE
                    sleep_time += timer_step
                if stop_waiting:
                    LOGGER.debug(
                        'Sleep interrupted because this script is terminating '
                        'its execution'
                    )

            LOGGER.debug(
                'Removing the FileDeleters that have already performed their '
                'execution'
            )
            _clean_file_deleters()
            _clean_txt_handler()

        LOGGER.info('Waiting while the spawned processes terminate')
        WORKERS_POOL.close()

        if self.remove_files:
            LOGGER.debug('Waiting while deleting the RegCM files')
            while len(FILE_DELETERS) > 0:
                _clean_file_deleters()

        LOGGER.debug('Waiting to manage the txt files')
        while len(TXT_HANDLERS) > 0:
            _clean_txt_handler()

        LOGGER.info('Execution complete!')
        time.sleep(0.5)


def main():
    # Read the options from command line
    parser = _create_input_parser()

    args = parser.parse_args()
    args = vars(args)

    if args['file_verbosity'] is None:
        args['file_verbosity'] = args['verbosity']

    # If no arguments are submitted, print help
    if len(argv) <= 1:
        parser.print_help()
        sys_exit(1)

    # Set the logging system
    streamhandler, filehandler = _prepare_logger(
        args['verbosity'],
        args['file_verbosity'],
        args['silent'],
        args['logfile']
    )

    LOGGER.info(
        'Starting execution. Executed with the following parameters: ' +
        ' '.join(argv)
    )

    LOGGER.debug('This execution will use %s processes', args['processes'])

    if get_localzone is None:
        LOGGER.warning(
            'Packet tzlocal is missing. Dates and times will be saved without '
            'timezones. To install this packet, use:\n"pip install tzlocal"'
        )
        LOGGER.debug(
            'Import of tzlocal package failed with the following error:\n{}'
            .format(tzlocal_import_error)
        )

    try:
        pycordexer_options, input_dir = _read_pycordexer_options(args)
    except Exception as e:
        LOGGER.error(
            'Error reading command line input. Execution will be terminated!\n'
            '{}'.format(e)
        )
        sys_exit(1)

    LOGGER.info(
        'Pycordexer will run with the following options: %s',
        pycordexer_options
    )

    LOGGER.debug(
        'Manual mode is: ' + ('ON' if args['manual_mode'] else 'OFF')
    )

    if args['remove_files']:
        LOGGER.warning(
            'This script will remove the RegCM files during its execution!'
        )

    listener = Listener(
        input_dir,
        pycordexer_options,
        args['manual_mode'],
        args['remove_files'],
        args['timer']
    )

    if args['daemon'] is False:
        utilities.log_utilities.collect_logs()
        try:
            listener.run(args['processes'])
        except:
            LOGGER.error(format_exc())
    else:
        LOGGER.debug('Preparing to go in daemon mode')

        if daemon is None:
            LOGGER.error(
                'No daemon module found. Please install it with:\n'
                '    pip install python-daemon\n'
                'Execution aborted'
            )
            LOGGER.debug(
                'Import of daemon package failed with the following error:\n{}'
                .format(daemon_import_error)
            )
            sys_exit(2)

        LOGGER.debug('Creating a context for the daemon')
        context_options = {
            'working_directory': os.getcwd(),
            'umask': 0o007,
        }
        if filehandler is not None:
            context_options['files_preserve'] = [filehandler.stream, ]

        if args['pidfile'] is not None:
            try:
                pidfile = daemon.pidfile.PIDLockFile(args['pidfile'])
            except:
                LOGGER.error('Error creating the pid file:\n{}'.format(
                    format_exc()
                ))
                raise
            context_options['pidfile'] = pidfile
        else:
            LOGGER.warning(
                'No PID file specified. This script will run anyway, but there '
                'is no way to call the cordex_listener_stop script to stop '
                'this run (beside using the kill command)'
            )

        LOGGER.debug(
            'Context will be created with the following options: %s',
            context_options
        )
        context = daemon.DaemonContext(**context_options)

        LOGGER.debug('Preparing the daemon to handling signals')
        context.signal_map = {
            signal.SIGTERM: stop_execution_on_signal,
            signal.SIGINT: stop_execution_on_signal,
            signal.SIGUSR1: stop_execution_when_done,
        }

        if streamhandler is not None:
            LOGGER.debug('Detaching from standard output')
            streamhandler.flush()
            LOGGER.removeHandler(streamhandler)

        LOGGER.info('Starting daemon')
        try:
            with context:
                utilities.log_utilities.collect_logs()
                listener.run(args['processes'])
        except:
            LOGGER.error(format_exc())
            raise


signal.signal(
    signal.SIGTERM,
    stop_execution_on_signal,
)

signal.signal(
    signal.SIGINT,
    stop_execution_on_signal,
)

signal.signal(
    signal.SIGUSR1,
    stop_execution_when_done,
)

signal.signal(
    signal.SIGUSR2,
    stop_execution_when_done_max_procs,
)

if __name__ == "__main__":
    main()

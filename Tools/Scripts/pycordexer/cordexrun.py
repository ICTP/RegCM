#!/usr/bin/env python3

"""
..Author: Stefano Piani, stefano.piani@exact-lab.it
..Date: 2017-10-23
"""

from __future__ import print_function

import argparse
from getpass import getuser
from os import path, kill, makedirs
from datetime import datetime
from subprocess import Popen, PIPE
from sys import executable, argv, exit as sys_exit, stderr, stdout
from threading import Lock
import signal


__copyright__ = 'Copyright (C) 2018 ICTP'
__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__credits__ = ["Stefano Piani"]


USERDIR = path.expanduser('~')
CORDEXRUN_DIR = path.join(USERDIR, '.cordexrun')
PIDDIR = CORDEXRUN_DIR

CORDEX_LISTENER_PATH = '/marconi/home/userexternal/ggiulian/data/exact-lab/' \
                       'regcm-cordex-repo/pycordexer/cordex_listener.py'
CORDEX_STOP_PATH = '/marconi/home/userexternal/ggiulian/data/exact-lab/' \
                   'regcm-cordex-repo/pycordexer/cordex_listener_stop.py'

cordex_listener_pid_file = None

# This flags communicate if it is possible (or not) to start other processes.
# When we receive a termination signal, we set this flag to False so, while we
# wait for all the other processes to terminate, no other processes can be
# started
start_other_processes = True

# Stopping execution is a lock that ensure that a process can not be started
# while we try to send signals to all the activated processes
editing_processes = Lock()

children = []

_SLEEP_TIMER = 10


def print_title():
    print(
        '''
        _____________________________________________
          __   __   __   __   ___      __
         /  ` /  \ |__) |  \ |__  \_/ |__) |  | |\ |
         \__, \__/ |  \ |__/ |___ / \ |  \ \__/ | \|
        _____________________________________________

        '''
    )
    stdout.flush()


def print_message(msg):
    print('\n' + '-' * 31 + ' CORDEXRUN ' + '-' * 32)
    print(msg)
    print('-'*74 + '\n')
    stdout.flush()


def stop_execution_spreading_signal(signum, frame):
    global start_other_processes

    print_message('Received signal {}; stopping execution!'.format(signum))

    with editing_processes:
        start_other_processes = False

        for process in children:
            process.send_signal(signum)

    if cordex_listener_pid_file is not None:
        try:
            with open(cordex_listener_pid_file, 'r') as f:
                raw_pid = f.read().strip('\n')
            listener_pid = int(raw_pid)
        except:
            listener_pid = None

        if listener_pid is not None:
            try:
                kill(listener_pid, signal.SIGTERM)
            except Exception:
                pass

def _create_input_parser():
    parser = argparse.ArgumentParser(
        usage=argv[0] + ' [-h] [-t TIMER] [-r] '
            '[-l LOG_FILE_PATH] [-v {debug,info,warning}] '
            '[--mail MAIL] [--global-model GLOBAL_MODEL] '
            '[--experiment EXPERIMENT] [--ensemble ENSEMBLE] '
            '[--domain DOMAIN] [--notes NOTES] [--regcm-version REGCM_VERSION] '
            '[--regcm-version-id REGCM_VERSION_ID] [MPIRUN OPTIONS] '
            'regcm_executable namelist',
        epilog='Beside these options, any other parameter accepted by mpirun is'
               ' allowed'
    )

    v_levels = ['debug', 'info', 'warning']

    parser.add_argument("-t", "--timer", help="Time to wait between each check"
                        " on the directory", type=int, default=_SLEEP_TIMER)
    parser.add_argument("-r", "--remove-files", help="Automatically remove "
                        "RegCM files after having extracted the variables into "
                        "the CORDEX files", action="store_true", default=False)
    parser.add_argument("-l", "--logfile", help="Save the logs of cordex_run "
                        "on this file", type=str, default=None)
    parser.add_argument('-v', '--verbosity', choices=v_levels, default='info',
                        help='The level of verbosity of the logs of the CORDEX '
                        'extractor (ignored if -l is not submitted)')
    parser.add_argument("--mail", help="The mail of the user that is "
                        "generating the CORDEX files", type=str, default=None)
    parser.add_argument("--global-model", help="The global model of the "
                        "simulation", type=str, default=None)
    parser.add_argument("--experiment", help="The experiment of the simulation "
                        "the files belong to", type=str, default=None)
    parser.add_argument("--ensemble", help="The ensemble of the simulation the "
                        "files belong to", type=str, default=None)
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

    return parser


signal.signal(
    signal.SIGTERM,
    stop_execution_spreading_signal,
)


signal.signal(
    signal.SIGINT,
    stop_execution_spreading_signal,
)


if __name__ == '__main__':
    # Create a parser for the command line
    parser = _create_input_parser()

    makedirs(CORDEXRUN_DIR, exist_ok=True)

    if len(argv) < 3:
        parser.print_help()
        sys_exit(1)

    print_title()

    username = getuser()
    current_time = datetime.now().strftime('%Y%m%d%H%M%S')
    args, mpirun_options = parser.parse_known_args(argv[:-2])
    namelist = argv[-1]
    regcm_executable = argv[-2]

    # Remove the first value of the mpirun_options list (which is the name of
    # this script)
    mpirun_options = mpirun_options[1:]

    cordex_listener_pid_file = path.join(
        PIDDIR,
        username + '_' + current_time + '.pid'
    )

    # Create a command line to run the cordex listener script
    listener_args = [
        executable, CORDEX_LISTENER_PATH, '-f', namelist,
        '-p', cordex_listener_pid_file, '--daemon', '-v', 'info'
    ]

    # Copy the command line of this script into the one that will be used to
    # load the listener.
    # If the command is a flag (and therefore its value is True or False), copy
    # the flag only if True. If the argument is not set (and therefore is None),
    # do not copy the flag
    for par_name, par_value in vars(args).items():
        if par_value is not None:
            if par_value is True:
                listener_args.append('--' + par_name.replace('_', '-'))
            elif par_value is False or par_value is None:
                continue
            else:
                if par_name == 'verbosity':
                    par_name = 'file-verbosity'
                listener_args.append('--' + par_name.replace('_', '-'))
                listener_args.append(str(par_value))

    listener_process = None
    with editing_processes:
        if start_other_processes:
            print_message(
                'Starting cordex_listener with the following command line:\n' +
                ' '.join(listener_args)
            )
            listener_process = Popen(listener_args, stderr=PIPE)
            children.append(listener_process)

            # Print listener_process standard error on the current standard
            # output
            while listener_process.poll() is None:
                stderr_line = listener_process.stderr.readline()
                print(stderr_line.decode('utf-8'), end='')
            last_line = listener_process.stderr.read().decode('utf-8')
            if last_line.strip('') != '':
                print(last_line, end='')

    # Wait listener to demonize itself
    if listener_process is not None:
        listener_status = listener_process.wait()
        with editing_processes:
            children.remove(listener_process)
        if listener_status != 0:
            print(
                'Error launching the cordex_listener script. RegCM will not be '
                'executed!',
                file=stderr
            )
            stderr.flush()
            sys_exit(listener_status)

    if not start_other_processes:
        sys_exit(1)

    # Run regcm with mpirun
    regcm_args = ['mpirun'] + mpirun_options + [regcm_executable, namelist]

    regcm_process = None
    with editing_processes:
        if start_other_processes:
            print_message(
                'Executing RegCM using the following command line:\n' +
                ' '.join(regcm_args)
            )
            regcm_process = Popen(regcm_args)
            children.append(regcm_process)

    if regcm_process is not None:
        regcm_status = regcm_process.wait()

        with editing_processes:
            children.remove(regcm_process)

    # Run cordex_listener_stop
    listener_stop_arguments = [
        executable, CORDEX_STOP_PATH, '-p', cordex_listener_pid_file, '-f'
    ]
    with editing_processes:
        if start_other_processes:
            print_message(
                'Stopping cordex listener with the following command:\n' +
                ' '.join(listener_stop_arguments)
            )
            listener_stop_process = Popen(listener_stop_arguments, stderr=PIPE)

            # The same that we already did for listener_process
            while listener_stop_process.poll() is None:
                stderr_line = listener_stop_process.stderr.readline()
                print(stderr_line.decode('utf-8'), end='')
            last_line = listener_stop_process.stderr.read().decode('utf-8')
            if last_line != '':
                print(last_line, end='')

    print_message('Execution complete!')

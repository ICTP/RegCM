#!/usr/bin/env python3

import argparse
import logging
import os
from netCDF4 import Dataset
from traceback import format_exc
from sys import exit as sys_exit
import signal

from utilities.cordex_utils import RegcmOutputFile, RegcmSimulation
from utilities.globals import OUTPUTDIR, CORDEX_VARS, ASSOCIATIONS
from utilities.pycordexer_methods import Action, get_methods
from utilities.log_utilities import collect_logs, WaitQueueHandler
from utilities.pool import Pool


__copyright__ = 'Copyright (C) 2018 ICTP'
__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__credits__ = ["Stefano Piani", "Graziano Giuliani", "Adriano Fantini"]


if __name__ == '__main__':
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


# This must be a global variable to be inherit by the new processes spawned
REGCM_FILE = None

# Read all the methods available in the pycordexer_methods file
METHODS_AVAIL = get_methods()


class InvalidVar(Exception):
    pass


class VariableFailed(Exception):
    pass


def _prepare_logger(verbosity):
    """
    Prepare the logging system
    """
    verbosity_level = getattr(logging, verbosity.upper())
    LOGGER.setLevel(verbosity_level)

    streamformatter = logging.Formatter(
        '%(levelname)s - %(processName)s: %(message)s'
    )
    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(verbosity_level)
    streamhandler.setFormatter(streamformatter)
    LOGGER.addHandler(streamhandler)


def _read_vars_arg(vars_arg):
    if vars_arg == 'ALL':
        return 'ALL'

    vars_split = vars_arg.split(',')

    for var_name in (v.lower() for v in vars_split):
        if var_name not in CORDEX_VARS:
            raise InvalidVar('Invalid cordex var: {}'.format(var_name))

    if len(vars_split) > len(list(set(vars_split))):
        raise ValueError(
            'The same variable can not be specified multiple times'
        )

    return [v.lower() for v in vars_split]


def get_vars_in_file(filename):
    for data_file_type in ASSOCIATIONS:
        if data_file_type['mask'].search(filename):
            ftype = data_file_type['name']
            LOGGER.debug(
                'The filename %s honours the mask of %s %s file',
                filename,
                'an' if ftype.lower()[0] in ('a', 'e', 'i' 'o', 'u') else 'a',
                ftype
            )
            return [t.lower() for t in data_file_type['vars']]
    else:
        raise ValueError(
            '"{}" name does not match any typical file mask. This script can '
            'not identify the variables it contains.'.format(filename)
        )


def parse_input():
    """
    Parse the input using argparse, in order to retrieve the following
    parameters:

    - datafile
    - variable(s)
    - email
    - domain model
    - experiment
    - ensemble
    - notes
    - disable-corrflag
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        'datafile',
        type=str,
        help='The path of the file'
    )
    parser.add_argument(
        'vars',
        type=str,
        metavar='var1,var2,var3,...',
        help='the variables to extract (if more than one, use comma to '
             'separate different name). Use "ALL" to read all variables from '
             'the datafile'
    )
    parser.add_argument(
        '-m',
        '--mail',
        default='esp@ictp.it',
        type=str,
    )
    parser.add_argument(
        '-d',
        '--domain',
        default='NONE',
        type=str,
    )
    parser.add_argument(
        '-g',
        '--global-model',
        default='ERA5',
        type=str,
    )
    parser.add_argument(
        '-e',
        '--experiment',
        default='none',
        type=str,
    )
    parser.add_argument(
        '-b',
        '--ensemble',
        default='NN',
        type=str,
    )
    parser.add_argument(
        '-n',
        '--notes',
        default='none',
        type=str,
    )
    parser.add_argument(
        '-o',
        '--output-dir',
        default=OUTPUTDIR,
        type=str,
    )
    parser.add_argument(
        '-p',
        '--processes',
        default=1,
        type=int,
    )
    parser.add_argument(
        '-c',
        '--disable-corrflag',
        action='store_true',
        help='Disable the time correction on the variables'
    )
    parser.add_argument(
        '--regcm-model-name',
        type=str,
        default='RegCM',
        help='The name of the current model. If this is not '
             'specified, this script will use the name RegCM'
    )
    parser.add_argument(
        '-r',
        '--regcm-version',
        type=str,
        default=None,
        help='The version of the current model of RegCM. If this is not '
             'specified, this script will try to guess the version from the '
             'input file. Otherwise, submit a string like "5.1" or '
             '"rev6446". This is the string that will be appended after '
             '"RegCM-"'
    )
    parser.add_argument(
        '--regcm-version-id',
        type=str,
        default=None,
        help='The version id of the current model of RegCM. Usually, this is a '
             'string like "v1-r1", "v2-r1", "v1-r2" and will be appended to '
             'the output file names. Defaults to v1-r1'
    )
    parser.add_argument(
        '--regcm-nest-tag',
        type=str,
        default=None,
        help='The nesting tag of the current model of RegCM. Usually, this is '
             'a complex string like x2yz2v1. Ignored if not present.'
    )
    parser.add_argument(
        '-v',
        '--verbosity',
        choices=['debug', 'info', 'warning'],
        default='info',
        help='The level of verbosity of this script'
    )

    args = parser.parse_args()
    return args


def save_vars(datafile, requested_vars, worker_pool,
              mail='esp@ictp.it', domain='NONE', global_model='NONE',
              experiment='none', ensemble='NN', notes='none', corrflag=True,
              regcm_model_name='RegCM', regcm_version=None,
              regcm_version_id=None, regcm_nest_tag=None,
              cordex_root_dir=OUTPUTDIR, sigterm_handler=None,
              sigint_handler=None):

    global REGCM_FILE

    # A list of all the processes that will be spawned by this function
    processes = []

    # If requested_vars is ALL, then we must check how many variables can be
    # saved from the file
    if requested_vars == 'ALL':
        LOGGER.debug('"ALL" flag recognized for the variables')
        requested_vars = get_vars_in_file(os.path.basename(datafile))

    LOGGER.debug(
        'The following requested_vars will be saved: %s',
        ', '.join([v for v in requested_vars])
    )

    # Create an object that stores all the metadata for the current simulation
    simulation = RegcmSimulation(
        mail,
        domain,
        global_model,
        experiment,
        ensemble,
        notes
    )

    LOGGER.debug('Reading file %s', datafile)
    with Dataset(datafile, 'r') as ncf:
        REGCM_FILE = RegcmOutputFile(
            ncf,
            datafile,
            regcm_model_name,
            regcm_version,
            regcm_version_id,
            regcm_nest_tag
        )

    for variable in requested_vars:
        LOGGER.debug('Spawning process for variable %s', variable)
        p = worker_pool.execute(
            execute_var_chain_wrapper,
            (
                variable,
                datafile,
                simulation,
                corrflag,
                cordex_root_dir,
                sigterm_handler,
                sigint_handler
            ),
            [REGCM_FILE, ]
        )
        if p is not None:
            LOGGER.debug('Process for variable %s spawned', variable)
        else:
            LOGGER.debug('Process spawn aborted')

        if p is not None:
            processes.append(p)

    return processes



def execute_var_chain_wrapper(variable, regcm_file_path, simulation, corrflag,
                              cordex_root_dir, sigterm_handler, sigint_handler,
                              logs_queue):
    """
    A wrapper that is useful if execute_var_chain runs in another process. It
    configures the logging system to send all the logs to a queue and ensures
    that all the exception are correctly reported on the logs.
    """
    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        root_logger.removeHandler(handler)

    newHandler = WaitQueueHandler(logs_queue)
    root_logger.addHandler(newHandler)

    if sigterm_handler is not None:
        LOGGER.debug('Setting sigterm handler for this process')
        signal.signal(
            signal.SIGTERM,
            sigterm_handler,
        )

    if sigint_handler is not None:
        LOGGER.debug('Setting sigint handler for this process')
        signal.signal(
            signal.SIGINT,
            sigint_handler,
        )

    try:
        LOGGER.debug('Calling method "execute_var_chain"')
        execute_var_chain(
            variable,
            regcm_file_path,
            simulation,
            corrflag,
            cordex_root_dir
        )
    except:
        LOGGER.error(
            'Error elaborating variable {}:\n{}'.format(variable, format_exc())
        )
        sys_exit(1)


def execute_var_chain(variable, regcm_file_path, simulation, corrflag,
                      cordex_root_dir):
    LOGGER.info('Saving variable %s', variable)

    if variable not in CORDEX_VARS:
        raise ValueError(
            'Variable {} not found (it must be added as a json file in the '
            'variables/cordex_vars directory)'.format(variable))

    var_definition = CORDEX_VARS[variable]
    for action_num, action_definition in enumerate(var_definition):
        LOGGER.debug(
            'Trying to save variable %s with action %s',
            variable,
            action_num + 1
        )

        try:
            LOGGER.debug('Building Action from definition')
            methods_list = []

            for method_name, method_options in action_definition:
                LOGGER.debug(
                    'Trying to associate a method to string "%s"',
                    method_name
                )
                if method_name not in METHODS_AVAIL:
                    raise Exception(
                        'No method called "{}" found'
                        .format(method_name)
                    )
                MethodClass = METHODS_AVAIL[method_name]

                LOGGER.debug(
                    'Trying to create an object from the class "%s" using '
                    'the submitted options',
                    method_name
                )
                methods_list.append(MethodClass(**method_options))

            LOGGER.debug('Gluing the methods to create an action')
            current_action = Action(
                methods_list,
                REGCM_FILE,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_root_dir
            )

            LOGGER.debug('Executing action')
            current_action.execute()

            LOGGER.info(
                'Variable %s saved using action number %s',
                variable,
                action_num
            )
            break
        except Exception:
            LOGGER.debug(
                'Building of Action {} failed with the following '
                'exception:\n{}'.format(action_num + 1, format_exc())
            )
            continue
    else:
        raise VariableFailed(
            'This software has not found a viable action to elaborate the '
            'variable "{}" from file {}.'.format(variable, regcm_file_path)
        )


def main():
    args = parse_input()

    _prepare_logger(args.verbosity)

    LOGGER.info('Execution started')

    # Read the value from the args
    datafile = args.datafile
    vars = _read_vars_arg(args.vars)
    mail = args.mail
    domain = args.domain
    global_model = args.global_model
    experiment = args.experiment
    ensemble = args.ensemble
    notes = args.notes

    regcm_model_name = args.regcm_model_name
    regcm_version = args.regcm_version
    regcm_version_id = args.regcm_version_id
    regcm_nest_tag = args.regcm_nest_tag

    LOGGER.debug(
        'Starting a thread to collect the logs generated by others processes'
    )
    collect_logs()

    LOGGER.debug('Building a pool of subordinate processes to perform elaborations')
    worker_pool = Pool(max_num_of_process=args.processes, name='Subordinate')

    if args.disable_corrflag:
        corrflag = False
    else:
        corrflag = True
    output_dir = args.output_dir

    procs = save_vars(
        datafile,
        vars,
        worker_pool,
        mail,
        domain,
        global_model,
        experiment,
        ensemble,
        notes,
        corrflag,
        regcm_model_name,
        regcm_version,
        regcm_version_id,
        regcm_nest_tag,
        output_dir
    )

    worker_pool.close()

    # Get the exit status of the processes and, if one is different from 0,
    # return it as the exit code of this process
    exit_status = [p.exitcode for p in procs]
    for exitst in exit_status:
        if exitst != 0:
            LOGGER.error('Some variables have not been correctly saved!')
            return exitst
    return 0


if __name__ == '__main__':
    sys_exit(main())

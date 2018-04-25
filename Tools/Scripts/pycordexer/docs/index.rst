.. PyCordexer documentation master file, created by
   sphinx-quickstart on Mon Feb 19 17:42:51 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========================
PyCordexer documentation
========================

.. contents:: Table of Contents
    :depth: 2
    :backlinks: none


Introduction
==============

The *PyCordexer* scripts have been developed to ease the RegCM Model User
in converting variables from RegCM model output files into the netCDF format
respecting the CORDEX experiment conventions.


Requirements
~~~~~~~~~~~~

To run the python script, you need a working Python3 interpreter. Moreover,
the following python libraries are required:

- *numpy*
- *netCDF4* (from https://unidata.github.io/netcdf4-python)

You might also consider to install other libraries to enable additional features:

- *daemon* if you want to use the daemon mode of cordex_listener
- *tzlocal* if you want to save dates with a reference to te current time zone


Installation
~~~~~~~~~~~~

Part of the interpolation and derived variable calculation has been taken
from RegCM Model source code, and you need to compile the Fortran source
codes in the *PyCordexer* directory using the f2py program of the Python
numpy package, which is a requirement of the netcdf4-python module.
Just type `make` in the *PyCordexer* directory to compile the Fortran
code to be used by Python

.. raw:: latex
    
        \begin{tcolorbox}
         cd pycordexer\\
         make        
        \end{tcolorbox}

If the name of your f2py executable is not "f2py3", please change the
first line of `Makefile` in the main directory with the name
of your f2py executable.

The scripts can be run from the *PyCordexer* directory afterwards.

PyCordexer Scripts
==================

The *PyCordexer* is split into four separate scripts:

1. *pycordexer.py*
   This script extracts variables from RegCM output files and creates a
   new file in the CORDEX netCDF format.
2. *cordex_listener.py*
   This script allows the extraction of CORDEX variables from RegCM output while the
   execution of RegCM is running
3. *cordex_listener_stop.py*
   This script stops the execution of cordex_listener 
4. *means.py*
   This script computes daily and monthly averages of CORDEX files at higher
   temporal resolution and saves them in CORDEX netCDF format files.


PyCordexer
============
The pycordexer.py script is the base tool of the PyCordexer suite.
It reads a RegCM model output file and extracts one or more variables
creating for each of them a CORDEX file with the details
specified from the command line arguments.

The syntax to invoke the program is:

.. raw:: latex
    
        \begin{tcolorbox}
         ./pycordexer.py datafile variable
        \end{tcolorbox}

The *datafile* is a RegCM file which has the variable specified by *variable*
inside.

If you want to extract more than one variable, you can supply more than one name,
separating them using a comma (without spaces)

.. raw:: latex
    
        \begin{tcolorbox}
         ./pycordexer.py datafile var1,var2
        \end{tcolorbox}


List of Implemented Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are the names of the currently implemented variables:

.. raw:: latex

        \begin{center}


.. include:: table.inc


.. raw:: latex

        \end{center}

Some of these variables generate more than one file: for example, the `ua` variable generates
several files for different heights. To have a clear insight of what operations are
performed when a specific variables is saved, please read the section :ref:`variable-definition`.


Command line options
~~~~~~~~~~~~~~~~~~~~~~

Beside the first two (mandatory) arguments, other options can be submitted by command line.
The syntax for these options is 

.. raw:: latex
    
        \begin{tcolorbox}
         ./pycordexer.py datafile var1,var2 --option value
        \end{tcolorbox}

Here you can find an exhaustive list of these options:

.. raw:: latex

        \begin{center}
        
+------------------+------------------+-----------------------------------------------+
| Argument         | Default          | Meaning                                       |
+==================+==================+===============================================+
| mail             | esp@ictp.it      | E-Mail of the Person responsible of this run  |
+------------------+------------------+-----------------------------------------------+
| domain           | NONE             | CORDEX domain (EUR-44,AFR-44,etc)             |
+------------------+------------------+-----------------------------------------------+
| global-model     | NONE             | CMIP5 driving model (MOHC-HadGEM2-ES,etc)     |
+------------------+------------------+-----------------------------------------------+
| experiment       | none             | CMIP5 experiment (historical,rcp85,etc)       |
+------------------+------------------+-----------------------------------------------+
| ensemble         | NN               | CMIP5 model experiment (r1i1p1,r2i1p1,etc)    |
+------------------+------------------+-----------------------------------------------+
| notes            | none             | Eventual notes to be added in the file        |
+------------------+------------------+-----------------------------------------------+
| regcm-version    | Read from input  | The version of the current model of RegCM     |
+------------------+------------------+-----------------------------------------------+
| regcm-version-id | Read from input  | The version id of the current model of RegCM  |
+------------------+------------------+-----------------------------------------------+
| output-dir       | .                | The main dir where the output will be saved   |
+------------------+------------------+-----------------------------------------------+
| processes        | 1                | The number of parallel processes              |
+------------------+------------------+-----------------------------------------------+
| verbosity        | info             | The level of verbosity of the script          |
+------------------+------------------+-----------------------------------------------+

.. raw:: latex

        \end{center}

The first 8 options in the previous table are metadata that will be used to generate
the names of the directories where the output files will be generated.

Beside the options in the previous table, there is a flag called --disable-corrflag.
If submitted, it disables the time correction performed by the method "CorrectTime".
Actually, there is hardly a good reason to submit such a flag.

.. _all-variable:

The "ALL" variable
~~~~~~~~~~~~~~~~~~

If you want, rather than specifying several different names for the variables, you can simply
use "ALL". This means that pycordexer will try to automatically recognize the file by its
contents, and will extract all the variables accordingly to the following table

.. raw:: latex

        \begin{center}


.. include:: associations_table.inc


.. raw:: latex

        \end{center}


Cordex Listener
================
`cordex_listener.py` is a script that automatically calls pycordexer.py on several
RegCM output files using the meta-variable "ALL" (see section :ref:`all-variable`).

The cordex listener continuously monitors a directory. Whenever a new file is
created inside the directory, it checks if the name of the file matches a particular
mask that recognizes if the file has been generated by RegCM. 

As soon as new files are generated, the cordex listener tries to put together these
files in one or more groups. A "group" is a bunch of files that describe the same
time interval. So, for example, in a group there can be an ATM file and a STS file
that contain the simulation of the same month.
A group will not be processed until it is completed. A group is completed as soon as
a txt file appears with the same name of the other files in the group (beside the part
of the name that describes the type: ATM, SRF, STS, ...). Please, ensure that your
version of RegCM produces the txt file.
At this point, the cordex listener reads the txt file: if there is a line that reports
that the group has been already elaborated, the group is discarded. Otherwise, pycordexer
is launched on each file of the group. When pycordex ends its execution, cordex_listener
will report in the txt file that the group has been elaborated.


Command line options
~~~~~~~~~~~~~~~~~~~~~~

Cordex listener supports several options from command line.

The most relevant one are:

+------------------+-----------------+------------------------------------------------+
| Argument         | Default         | Meaning                                        |
+==================+=================+================================================+
| directory        |                 | The path of the directory with the RegCM files |
+------------------+-----------------+------------------------------------------------+
| namelist         |                 | The namelist you used for the run of RegCM     |
+------------------+-----------------+------------------------------------------------+

At list one of these two options must be submitted. If cordex_listen can read the
namelist, then it will infer from it the location of the RegCM output files (and,
therefore, the right value for the "directory" option) and all the other metadata that
you can submit to pycordexer.py by command line.

Beside these options, cordex_listener supports all the options that you can submit to
pycordexer except the output-dir:

+------------------+------------------+-----------------------------------------------+
| Argument         | Default          | Meaning                                       |
+==================+==================+===============================================+
| mail             | esp@ictp.it      | E-Mail of the Person responsible of this run  |
+------------------+------------------+-----------------------------------------------+
| domain           | NONE             | CORDEX domain (EUR-44,AFR-44,etc)             |
+------------------+------------------+-----------------------------------------------+
| global-model     | NONE             | CMIP5 driving model (MOHC-HadGEM2-ES,etc)     |
+------------------+------------------+-----------------------------------------------+
| experiment       | none             | CMIP5 experiment (historical,rcp85,etc)       |
+------------------+------------------+-----------------------------------------------+
| ensemble         | NN               | CMIP5 model experiment (r1i1p1,r2i1p1,etc)    |
+------------------+------------------+-----------------------------------------------+
| notes            | none             | Eventual notes to be added in the file        |
+------------------+------------------+-----------------------------------------------+
| regcm-version    | Read from input  | The version of the current model of RegCM     |
+------------------+------------------+-----------------------------------------------+
| regcm-version-id | Read from input  | The version id of the current model of RegCM  |
+------------------+------------------+-----------------------------------------------+
| processes        | 1                | The number of parallel processes              |
+------------------+------------------+-----------------------------------------------+
| verbosity        | info             | The level of verbosity of the script          |
+------------------+------------------+-----------------------------------------------+

These options have the same meaning that they have for pycordexer. If you specify a 
namelist and also some other options like "experiment" or "ensemble" there is an
overlapping of information. In this case, the rule is the following: cordex_listener.py
will use the value you have submit by the command line with the specific option; if you
do not have submit any value, the pycordexer will try to use the value from the namelist
and, if the value is missing from the namelist or you do not have submitted any namelist,
it will fall back to the default values of the options.

The same rule is also valid if you submit a namelist and a directory at the same time.

The output-dir flag has been suppressed because, by design, cordex_listener saves its
output in a directory named CORDEX inside the directory that contains the RegCM files.

Moreover, cordex_listener.py has the following flags (options without a value):

+------------------+------------------------------------------------------------------+
| Argument         | Meaning                                                          |
+==================+==================================================================+
| manual-mode      | Enable the manual mode (see section :ref:`manual-mode`)          |
+------------------+------------------+-----------------------------------------------+
| daemon           | Enable the daemon mode (see section :ref:`daemon-mode`)          |
+------------------+------------------+-----------------------------------------------+
| silent           | Do not print anything on the screen                              |
+------------------+------------------+-----------------------------------------------+
| remove-files     | Delete the original RegCM files after having elaborated them     |
+------------------+------------------+-----------------------------------------------+

If you use the "silent" flag, then maybe you want to save the logs in a file. For this
reason, you can use the option "logfile" and specify a path. cordex_listener will save
its logs on that file. Obviously, you do not need to use the silent flag to specify a
logfile. In this case, the logs will be displayed on the terminal and also saved on the
file. The lines of the logs written on the files will have more information: the datetime,
the name of the process and of the function which generated that line will be prepended.
Finally, you can have a different level of verbosity between the display and the log file:
use the "file-verbosity" option to specify a different verbosity for the file (otherwise,
the same verbosity of the display will be used, the one that you can set with the
"verbosity" flag).

Finally, there are two more options:

+------------------+---------+--------------------------------------------------------+
| Argument         | Default | Meaning                                                |
+==================+=========+========================================================+
| timer            | 10      | The amount of time between two checks on the directory |
+------------------+---------+--------------------------------------------------------+
| pidfile          |         | See section :ref:`daemon-mode`                         |
+------------------+---------+--------------------------------------------------------+


.. _daemon-mode:

The daemon mode
---------------

In "daemon" mode, the cordex_listener.py script detaches itself from the current terminal
and behaves as a system service (a daemon).

The "daemon" python library must be installed to use this functionality.

The main problem you may address when using the daemon mode is that you no longer have any
idea of what operation the script is performing nor you have a way to stop its execution.
Indeed, by design, cordex listener will run forever and, in the daemon mode, you can not
simply press Ctrl+C to stop its execution.

That is the reason way we provide a script called "cordex_listener_stop.py" and the option
"--pidfile" for the cordex_listener. If you want to run Cordex Listener in daemon mode,
please put a valid path for the pidfile option. If you do, cordex_listener will create a
file where the PID of its process will be saved. In this way, you can send a TERM or a 
KILL signal to the process.

Cordex Listener can be triggered using signals to perform several different shutdown
procedures:

- *raw shutdown* (triggered by SIGINT): cordex_listener will not spawn any other
  process. The running ones will terminate before executing any other line of
  python code
- *strong shutdown* (triggered by SIGTERM): cordex_listener will not spawn any other
  process. The processes that are extracting variables will terminate their
  execution as soon as the CORDEX variable is saved. The main process will
  terminate as soon as the other processes stop their execution
- *clean shutdown* (triggered by SIGUSR1): cordex_listener will search for new
  files and, if they are present, they will be elaborated. As soon as all the
  elaborations end, Cordex Listener will stop its execution
- *fast and clean shutdown* (triggered by SIGUSR2): cordex_listener will search
  for new files and, if they are present, they will be elaborated using a number
  of processes equal to the total number of CPU cores of the machine. As soon as
  all the elaborations end, Cordex Listener will stop its execution

You can use the cordex_listener_stop.py script to send these signals to Cordex Listener. 


.. _manual-mode:

The manual mode
----------------
In principle, Cordex Listener is designed to elaborate data produced by RegCM during
its execution. It is possible, however, that you already have some files produced by
RegCM that you want to process with this script.

For this reason, there is a "manual mode". In manual mode, the cordex listener checks
once for the RegCM output file in its directory and, after that, it processes any of
these files, whether or not the txt file for their group exists and independently of
its content.

When all the files have been processed, the execution will end. Therefore, Cordex
Listener will run indefinitely only when the manual mode is not active.

Manual mode and daemon mode can be activated together. Obviously, "clean shutdown" has
no effect in this situation because, when in manual mode, Cordex Listener would perform
a "clean shutdown" anyway.


Cordex Listener Stop
====================
cordex_listener_stop.py is a small script that stops the execution of cordex_listener.py. To
do that, cordex_listener_stop.py needs a pidfile created by cordex_listener.py.

The syntax to invoke the program is:

.. raw:: latex
    
        \begin{tcolorbox}
         ./cordex\_listener\_stop.py -p pidfile
        \end{tcolorbox}

By default, this script invokes a "clean shutdown" (check section :ref:`daemon-mode`). Using the
"-f" flag, you can trigger instead a "fast and clean shutdown".

If the execution does not terminate after a fixed amount of seconds (by default 600, but you can
change this with the -t option), a "strong shutdown" will be invoked. You can disable this
feature setting 0 as value for the "-t" option.

Moreover, you can set another timer after the previous one using the "-k" flag (by default, this
is disable). When K seconds are passed after the beginning of the procedure of strong shutdown,
if cordex_listener.py is still running, it will be killed.

Please, use the -k flag with care. During its execution, cordex_listener.py spawns several processes.
The -k option kills the main process, but not its children. Therefore, it is possible that you
will have pending processes if you kill cordex_listener.py.



Advanced
========

In this chapter you will find some details that you do not strictly need to run the code. But, if you want
to understand how the code works or what operations are performed on each variable, you will find
this section extremely useful.


.. _variable-definition:

The definition of the variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The scripts rely on a series of external JSON-encoded files (stored in the
``variables/cordex_vars`` directory), containing instructions to apply to
variables in the process of extraction.

If you look in that directory, you will find that, for each variable, there
is a JSON file that describes the operations that must be accomplished to
save that variable.

To understand the syntax of these JSON files, you need to know what is an
action: an action is a chain of methods that will be executed one after the
other. The output of the first method will be the input of the second one, and
so on. For example, an action can be something like::

    read data --> change unit of measurement --> write output on disk

For each variable, you can specify one or more actions. Pycordexer will try
to execute the first action. If this action fails, it will fall back to the
second one and so on...

In this way, you can implement different procedures to extract data accordingly
to the different version of RegCM.

The root element of a JSON file is a dictionary that associates to the variable
name a list of actions::

    {
        variable: [
            action1,
            action2,
            ...
        ]
    }

But how can an action be described? With a list of methods! ::

    {
        variable: [
            [
                method1,
                method2,
                ...
            ],
            [
                method1,
                method2,
                ...
            ],
            ...
        ]
    }


Unfortunately, a method is again a list. In particular, it is a couple. The first
element is the name of the method, the second one is a dictionary with the options
of the method. Therefore, our JSON becomes a little bit cumbersome::

    {
        variable: [
            [
                [
                    read_data,
                    {from: regcm_variable_name}
                ],
                [
                    change_unit_of_measurement,
                    {from: km, to: m}
                ],
                ...
            ],
            [
                [
                    read_data,
                    {from: another_regcm_variable}
                ],
                [
                    change_unit_of_measurement,
                    {from: mm, to: m}
                ],
                ...
            ]
        ]
    }

A few words about the methods. There are two kinds of methods: ActionStarters and
Filters. An ActionStarter is a method that does not need the output of another
method to be executed. For example, reading a RegCM output file is an ActionStarter.
On the other side, Filters uses a previous output to perform their operations.
All the methods of an action must be Filters beside the first one, that must be
an ActionStarter.

If you need a complete list of all the methods that you can use, see the
section :ref:`methods`.


Variable associations
~~~~~~~~~~~~~~~~~~~~~

In this section, it is explained how to edit the table that appears at the end
of section :ref:`all-variable`, the one that associates to each RegCM output file
the variables that must be extracted from it.

This information is stored in the `file_associations.json` file inside the
`variables` directory. This files describe a list of dictionary. Each dictionary
represent a particular type of output file of RegCM. For each file type (ATM,
SRF, ...) there is its name, a list of variables that can be extracted from that
file and a mask. The mask is the regular expression that a file need to match to
be recognized as of that type.

Therefore, when you use the ALL flag with pycordexer.py or when cordex_listener.py
elaborates a file, the software tries to match the filename with one of the masks.
If one mask matches, all the variables that are in the list "`vars`"
will be extracted from that file.



.. _methods:

PyCordexer Methods
~~~~~~~~~~~~~~~~~~~

The following are the Methods that can be used to manipulate RegCM variables.

Inside ``variables/cordex_vars/`` there is a JSON file for each RegCM variable:
each JSON contains a sequence of Methods to be applied to the original variable
to produce it.

The Methods are presented in alphabetical order; when not specified, each method
returns a ``Variable`` object.

.. automodule:: utilities.pycordexer_methods
    :members:
    :undoc-members:
    :show-inheritance:
    :exclude-members: Action, ActionStarter, Method, Filter, get_methods



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

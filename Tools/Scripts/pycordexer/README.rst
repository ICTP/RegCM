####
PyCordexer
####


The *PyCordexer* scripts have been developed to ease the RegCM Model User
in converting variables in RegCM model output files into the netCDF format
respecting the CORDEX experiment conventions.

It has been developed to satisfy the requirements in the 05/03/2014 release
of the CORDEX Specification Document available from the DMI site at:

   http://cordex.dmi.dk/joomla

The *PyCordexer* script have been developed by ICTP staff to post-process
the CREMA experiment output dataset to be contributed to the CORDEX archive,
and have been adapted afterwards to post-process the latest RegCM model output
formats.

Any output file of the RegCM model from version 4.3 onward should be compatible
with the *PyCordexer* scripts.

Requirements
------------

To run the python script, you need a working Python3 interpreter and the
numpy and netCDF4 Python libraries:

   https://unidata.github.io/netcdf4-python

The numba library is needed only to compute CAPE and CIN from the
model vertical profiles (the *capecin* variable).

Installation
------------

No compilation is needed: the interpolation and derived variable
calculations, taken from the RegCM Model source code, are written in
Python. The scripts can be run from the *PyCordexer* directory.

.. raw:: latex

        \newpage

PyCordexer Scripts
------------------

The *PyCordexer* is split into two separate scripts:

1. **pycordexer.py**
   This script extracts variables from RegCM output files and creates a
   new file in the CORDEX netCDF format.
2. **means.py**
   This script computes day and monthly averages of CORDEX files at higher
   temporal resolution saving them in CORDEX netCDF format files.
3. **cordex_finalize.py**
   This script assembles, at the end of an experiment, the archive-size files
   that pycordexer.py could not complete (see "Archive-size files" below).

The pycordexer.py script
--------------------

The **pycordexer.py** script reads a RegCM model output file and extracts
the variable creating a CORDEX file out of the details specified from the
command line arguments.

The syntax to invoke the program is:

.. raw:: latex

        \begin{tcolorbox}
         pycordexer.py [-h] [-m MAIL] [-d DOMAIN] [-g GLOBAL_MODEL] [-e EXPERIMENT] \\
         [-b ENSEMBLE] [-n NOTES] [-o OUTPUT_DIR] [-c] [--no-monthly] \\
         [--no-chunks] [--first-year YEAR] [--last-year YEAR] \\
         datafile var1,var2,var3,...
        \end{tcolorbox}

The parameter *var1,var2,var3,...* is a list of CORDEX variables.
The *datafile* is a RegCM file which has the variables specified by *var1,var2,var3,...*
inside.
If var1,var2,var3,... is "ALL" (without quotes), the script will try to save all the
values that are saved inside the RegCM file.

Unless the --no-monthly flag is submitted, every daily variable that is also
in the CORDEX-CMIP6 monthly table is saved as monthly means too. The mean of
each month is computed from the daily values, and the cell_methods attribute
is taken from the monthly table. A warning is printed if the RegCM file does
not cover the whole month.

Archive-size files
~~~~~~~~~~~~~~~~~~

The CORDEX-CMIP6 archiving specifications want files of one year for the
sub-daily data (1hr, 3hr, 6hr), five years for the daily data (starting with
the years ending in 1 or 6) and ten years for the monthly data (starting with
the years ending in 1); the first and last files of an experiment may be
shorter. pycordexer.py produces them while the RegCM output is processed one
month at a time:

* each RegCM file gives, for every variable and frequency, a segment: a
  complete CORDEX file with the records of that month, saved in the hidden
  directory *.segments* of the variable directory;
* when the segments cover a whole chunk, they are merged into the chunk file,
  which appears in the variable directory under its final name; the segments
  are then removed.

No file is modified in place: segments and chunks are written under a
temporary name and renamed when complete, so an interrupted job leaves the
segments and the files already assembled intact, and rerunning it repeats the
work. The months can be processed in any order and in parallel jobs; a lock
file in *.segments* lets only one process assemble a chunk (a lock older than
two hours is considered left by an interrupted job). A record found in more
than one file is taken once, and a month processed again replaces its records,
also in a chunk already assembled. A chunk with missing records is not
assembled.

With *--first-year* and *--last-year* the chunks are clipped to the years of
the experiment, so that its first and last chunks are completed as soon as
their data are there; the first chunk may then start up to one day late
(point values start one step after the initial time), and records outside
these years are not written. Otherwise, once the whole experiment has been processed,
run

.. raw:: latex

        \begin{tcolorbox}
         cordex\_finalize.py [--first-year YEAR] [--last-year YEAR] [--dry-run] OUTPUT\_DIR
        \end{tcolorbox}

to assemble the chunks left at the edges of the experiment (and any other
chunk still in segments). It reports the segments it could not assemble
because of missing records, and exits with an error if any is left.

With *--no-chunks* pycordexer.py writes one file per RegCM file, as before.
The time independent (fx) variables are always written directly.

The tests of the assembly are in *tests/test_chunks.py*.

Example: post-processing a whole experiment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The script *example/regcmcordex.sh* is a Slurm job that post-processes one
month of a RegCM run: it runs pycordexer.py on the SRF, RAD, STS and ATM
files of the month at the same time, and removes the RegCM files when all
four have succeeded. Edit at its top the account, the paths, the simulation
metadata (domain, driving model, experiment, ensemble, model name and
version) and the years of the experiment:

.. raw:: latex

        \begin{tcolorbox}
         first\_year=1980 \\
         last\_year=2024
        \end{tcolorbox}

Then submit one job per month, in any order and as many at the same time as
the queue allows, for example for 1980 to 2024:

.. raw:: latex

        \begin{tcolorbox}
         for y in \$(seq 1980 2024); do \\
         \ \ for m in 01 02 03 04 05 06 07 08 09 10 11 12; do \\
         \ \ \ \ sbatch regcmcordex.sh /path/to/regcm/output /path/to/cordex \${y}\${m}0100 \\
         \ \ done \\
         done
        \end{tcolorbox}

Each job adds its month to the segments and completes the archive files
whose months are all there; for this evaluation run the files are

.. raw:: latex

        \begin{tcolorbox}
         .../1hr/pr/pr\_MED-12\_ERA5\_evaluation\_r1i1p1f1\_ICTP\_RegCM-ES1-1\_v1-r1\_1hr\_198001010030-198012312330.nc \\
         .../6hr/ta850/ta850\_MED-12\_..._6hr\_198001010600-198012311800.nc \\
         .../6hr/ta850/ta850\_MED-12\_..._6hr\_198101010000-198112311800.nc \\
         .../day/pr/pr\_MED-12\_..._day\_19800101-19801231.nc \\
         .../day/pr/pr\_MED-12\_..._day\_19810101-19851231.nc \\
         .../mon/pr/pr\_MED-12\_..._mon\_198001-198012.nc \\
         .../mon/pr/pr\_MED-12\_..._mon\_198101-199012.nc
        \end{tcolorbox}

and so on: one file per year at 1hr and 6hr, five years per daily file
(starting in 1981, 1986, ...) and ten years per monthly file (starting in
1981, 1991, ...), with shorter files at the start and at the end of the
experiment. The time range of the name comes from the first and last records:
the hourly means are at the middle of the hour, and the 6-hourly point values
of a RegCM month go from the 1st at 06:00 to the 1st of the next month at
00:00, so the first year of the experiment starts at 06:00 (a record at the
first time of a chunk is taken from the previous month, and the record after
the last year of the experiment is not written). The segments of the chunks not yet complete stay in the hidden
*.segments* directories; a failed job can be submitted again.

When all the months have been processed, check that nothing is left behind:

.. raw:: latex

        \begin{tcolorbox}
         regcmcordex.sh /path/to/regcm/output /path/to/cordex finalize
        \end{tcolorbox}

This runs cordex\_finalize.py on the output tree. With the years of the
experiment set it only has something to do if months are missing: it then
lists the segments it could not assemble and exits with an error, so that the
missing months can be found and processed. If first\_year and last\_year are
left empty in the script, this step also assembles the first and last
(shorter) files of the experiment.

Some variables can be extracted from different kind of files: for example, "tas" variable
can be extracted both from SRF and STS files. To allow users to convert all files
generated by RegCM without converting the same varibles multiple times, the ALL flag is
configured to avoid repetitions. Therefore, using ALL on a STS file will not save the
"tas" variable that is instead saved from the SRF files.

If you want to extract the tas variable from a STS file, you must not use the ALL flag and,
instead, submit a complete list of variables.

List of Implemented Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The variable which can be extracted from the RegCM file are:

.. raw:: latex

        \begin{center}

+------------+------------+-----------------------------------------------+
| variable   | RegCM file | Description                                   |
+------------+------------+-----------------------------------------------+
| tas        | SRF,STS    | Near-Surface Air Temperature                  |
+------------+------------+-----------------------------------------------+
| pr         | SRF,STS    | Precipitation                                 |
+------------+------------+-----------------------------------------------+
| prc        | SRF        | Convective Precipitation                      |
+------------+------------+-----------------------------------------------+
| huss       | SRF        | Near-Surface Specific Humidity                |
+------------+------------+-----------------------------------------------+
| hurs       | SRF        | Near-Surface Relative Humidity                |
+------------+------------+-----------------------------------------------+
| evspsbl    | SRF        | Evaporation                                   |
+------------+------------+-----------------------------------------------+
| mrros      | SRF        | Surface Runoff                                |
+------------+------------+-----------------------------------------------+
| ps         | SRF        | Surface Air Pressure                          |
+------------+------------+-----------------------------------------------+
| psl        | ATM        | Sea Level Pressure                            |
+------------+------------+-----------------------------------------------+
| tasmax     | STS        | Daily Maximum Near-Surface Air Temperature    |
+------------+------------+-----------------------------------------------+
| tasmin     | STS        | Daily Minimum Near-Surface Air Temperature    |
+------------+------------+-----------------------------------------------+
| sfcWindmax | STS        | Daily Maximum Near-Surface Wind Speed         |
+------------+------------+-----------------------------------------------+
| mrro       | SRF        | Total Runoff                                  |
+------------+------------+-----------------------------------------------+
| sfcWind    | SRF        | Near-Surface Wind Speed                       |
+------------+------------+-----------------------------------------------+
| ua850      | ATM,ATMp   | Eastward Wind (at 850 hPa)                    |
+------------+------------+-----------------------------------------------+
| va850      | ATM,ATMp   | Northward Wind (at 850 hPa)                   |
+------------+------------+-----------------------------------------------+
| ta850      | ATM,ATMp   | Air Temperature (at 850 hPa)                  |
+------------+------------+-----------------------------------------------+
| hus850     | ATM,ATMp   | Specific Humidity (at 850 hPa)                |
+------------+------------+-----------------------------------------------+
| ua500      | ATM,ATMp   | Eastward Wind (at 500 hPa)                    |
+------------+------------+-----------------------------------------------+
| va500      | ATM,ATMp   | Northward Wind (at 500 hPa)                   |
+------------+------------+-----------------------------------------------+
| ta500      | ATM,ATMp   | Air Temperature (at 500 hPa)                  |
+------------+------------+-----------------------------------------------+
| zg500      | ATM,ATMp   | Geopotential Height (at 500 hPa)              |
+------------+------------+-----------------------------------------------+
| ua200      | ATM,ATMp   | Eastward Wind (at 200 hPa)                    |
+------------+------------+-----------------------------------------------+
| va200      | ATM,ATMp   | Northward Wind (at 200 hPa)                   |
+------------+------------+-----------------------------------------------+
| ta200      | ATM,ATMp   | Air Temperature (at 200 hPa)                  |
+------------+------------+-----------------------------------------------+
| zg200      | ATM,ATMp   | Geopotential Height (at 200 hPa)              |
+------------+------------+-----------------------------------------------+

.. raw:: latex

        \end{center}

Using the ALL flag on a particular file type (for example, STS), this script
will extract only the variables that, in the previous table, have that file
type as the first one that appears in the "RegCM file" type.


Arguments
~~~~~~~~~

The first two arguments are required, while the others have basic defaults
which do not comply to CORDEX standard but allow the script to be used as
a generic post-processing tool for RegCM

.. raw:: latex

        \vfill
        \begin{center}

+--------------+-------------+-----------------------------------------------+
| Argument     | Default     | Meaning                                       |
+--------------+-------------+-----------------------------------------------+
| mail         | esp@ictp.it | E-Mail of the Person responsible of this run  |
+--------------+-------------+-----------------------------------------------+
| domain       | NONE        | CORDEX domain (EUR-44,AFR-44,etc)             |
+--------------+-------------+-----------------------------------------------+
| global_model | NONE        | CMIP5 driving model (MOHC-HadGEM2-ES,etc)     |
+--------------+-------------+-----------------------------------------------+
| experiment   | none        | CMIP5 experiment (historical,rcp85,etc)       |
+--------------+-------------+-----------------------------------------------+
| ensemble     | NN          | CMIP5 model experiment (r1i1p1,r2i1p1,etc)    |
+--------------+-------------+-----------------------------------------------+
| notes        | none        | Eventual notes to be added in the file        |
+--------------+-------------+-----------------------------------------------+
| corrflag     | True        | Apply a correction time to dates.             |
+--------------+-------------+-----------------------------------------------+

.. raw:: latex

        \end{center}

Usage
~~~~~

Example usage is the following. Suppose you have in a directory RegCM model
output file, and want to extract in CORDEX format the variable tas from SRF
file for the AFR-44 domain run on ECMWF ERA15 dataset:

.. raw:: latex

        \begin{tcolorbox}
          pycordexer.py /data/output/path/Africa\_SRF.2002030100.nc tas -m me@here\\
           -d AFR-44 -g ECMWF-ERAINT -e evaluation -b r1i1p1 -n 'Some text'
        \end{tcolorbox}

This will create in the current directory the file:

.. raw:: latex

        \begin{tcolorbox}
         tas\_AFR-44\_ECMWF-ERAINT\_evaluation\_r1i1p1\_ICTP-RegCM4-3\_v4\_3hr\_200203010300-200204010000.nc
        \end{tcolorbox}

with all the mandatory attributes and naming conventions required by CORDEX
experiment convention document.

If just the name of the file and the variable name are specified, the output
file name will be:

.. raw:: latex

        \begin{tcolorbox}
         tas\_NONE\_NONE\_none\_NN\_ICTP-RegCM4-3\_v4\_3hr\_200203010300-200204010000.nc
        \end{tcolorbox}


The means.py script
-------------------

The **means.py** script reads a CORDEX file created by the **cordex.py**
script and computes day or monthly average of the variable in the file as
instructed by command line argument, saving the result in a CORDEX
conforming netCDF file.

The syntax to invoke the program is:

.. raw:: latex

        \begin{tcolorbox}
          means.py datafile [mon/day]
        \end{tcolorbox}

The *datafile* is a file produced by cordex.py, the secund argument if not
specified defaults to *mon*, i.e. creates a monthly average file.

Usage
~~~~~

Example usage is the following. Suppose you have in a directory the file

.. raw:: latex

        \begin{tcolorbox}
         tas\_NONE\_NONE\_none\_NN\_ICTP-RegCM4-3\_v4\_3hr\_200203010300-200204010000.nc
        \end{tcolorbox}


and want to compute daily averages in a CORDEX conforming netCDF file.

.. raw:: latex

        \begin{tcolorbox}
          means.py tas\_NONE\_NONE\_none\_NN\_ICTP-RegCM4-3\_v4\_3hr\_200203010300-200204010000.nc day
        \end{tcolorbox}


This will create in the current directory the file:

.. raw:: latex

        \begin{tcolorbox}
         tas\_NONE\_NONE\_none\_NN\_ICTP-RegCM4-3\_v4\_day\_2002030112-2002033112.nc
        \end{tcolorbox}


which contains daily mean values with all the mandatory attributes and naming
conventions required by CORDEX experiment convention document.

If you want to compute monthly average, you can issue:

.. raw:: latex

        \begin{tcolorbox}
          means.py tas\_NONE\_NONE\_none\_NN\_ICTP-RegCM4-3\_v4\_day\_2002030112-2002033112.nc mon
        \end{tcolorbox}


and obtain:

.. raw:: latex

        \begin{tcolorbox}
         tas\_NONE\_NONE\_none\_NN\_ICTP-RegCM4-3\_v4\_mon\_20020301-20020331.nc
        \end{tcolorbox}


Contacts
--------

The *pycordexer* scripts are now kept in the *Scripts/Tools* directory
of the RegCM model code package, and are maintained by ICTP as part
of the RegCM codebase. Please address any problem/suggestion to the RegCNET
mailing list:

    https://lists.ictp.it/mailman/listinfo.cgi/regcnet

* file: lats4d.gs 
*
*    A minimum fuss gs script for writing NetCDF, HDF-SDS or GRIB
*    files from GrADS using the PCMDI (http://www-pcmdi.llnl.gov) LATS
*    interface. This script requires GrADS Version 1.7.beta.9 or later.
*
*    See usage() for documentation; from GrADS enter "lats4d -h"
*
* REVISION HISTORY:
*
* 21dec1998   da Silva   First crack.
* 24dec1998   da Silva   Further development: options -lat -lon -cal -mean,
*                        better error trapping.
* 31dec1998   da Silva   First beta version; automatic generation of
*                        LATS parameter table; now variables with 1 vertical
*                        level are identified as sfc variables.
* 04jan1999   da Silva   Added -ftype option, automatic detection of
*                        file types, removed -sdf/-xdf options.
* 05jan1999   da Silva   Cosmetic changes in man page.
* 06jan1999   da Silva   Added -precision option.
* 08jan1999   da Silva   Added -precision as valid option in
*                        parsecmd(), small other changes
* 11jan1999   da Silva   now file name can have 3 tokens, so SDFopen
*                        can use templates.
* 17jan1999   da Silva   Added -func and -freq options, fixed bug when
*                        there was only 1 time step on file; -func is
*                        currently undocumented.
* 22Jan1999   da Silva   Fixed bug with -time option: now 0Z15dec1991
*                        works.
* 09Feb1999   da Silva   Added -grid option; added detection of invalid
*                        time frequency in minutes (1.0.Beta.8)
* 17Mar1999   da Silva   Fixed parsing bug in -freq; fixed bug in the
*                        evaluation of -func arguments; increased
*                        long name column size in parm table to 70
*                        to accomodate files going to ECS; removed FUNCTION
*                        from variable long name - this information is
*                        now part of the global attribute "comments"
*                        (1.0.Beta.9)
* 18Mar1999   da Silva   Fixed another -func related bug in getgid().
* 19Mar1999   da Silva   Added -func to options to avoid parsing error; 
*                        revised documentation + small mod in -freq.
* 29Mar1999   da Silva   Added -zrev option
* 28Apr1999   da Silva   Added "-de" option (1.0.beta.13); now we
*                        can regrid; also added undocumented options
*                        -geos1x1, -geos1x1a, -geos1x1s, and same for
*                        2x2.5, 4x5; these are DAO specific.
* 04May1999   da Silva   Fixed small bug in -levs; now -levs 0.4 is
*                        possible. 
* 30Jun1999   da Silva   Added CPTEC to list of centers, fixed -vars
*                        to allow variables starting with 'e'.
* 13Jul1999   da Silva   Fixed "-i" bug: now file names can start with "e".
* 16Jul1999   da Silva   Changed order of upper-air var/level. This
*                        way the data is written in the same order they
*                        come in a standard IEEE grads file. Implemented
*                        stream (fwrite) output; no automatic ctl generation
*                        yet.
* 23Jul1999   da Silva   Added -title option.
* 22Nov1999   da Silva   Made v1.0.beta.18 the first stable version: v1.0.
* 23Feb2000   da Silva   Special handling for grib output and p<1; now it
*                        uses hybrid level number for vertical coordinate
*                        in such cases.
* 01Jun2000   da Silva   Added -xvars option.
* 16Jun2000   da Silva   Added -xsfc/xupper options.
* 10Oct2001   Owens      Added -tmean option to allow the user to specify
*                        a date for mean file other than tmid
* 17dec2002   da Silva   Bug fix in getgid for pre-projected data; change
*                        in openf() to detect URLs.
*
*......................................................................
*
* main(args)  Main driver checking config, handling error conditions, etc.
*
function main ( args )

      _myname = 'lats4d: '
      _version = '1.4 of 19 Dec 2002'


*     Parse command line
*     ------------------
      rc = parsecmd(args)
      if ( rc != 0 )
           if ( _quit = 1 ) 
                say _myname 'exiting from GrADS...'
                'quit'
           endif
           return rc
      endif


*     Check GrADS configuration
*     -------------------------
      rc = chkcfg()
      if ( rc != 0 )
           if ( _quit = 1 ) 
                say _myname 'exiting from GrADS...'
                'quit'
           endif
           return rc
      endif


*     Output file message
*     -------------------
      if ( _format = 'stream' )
           _fmsg = _format ' file ' _ofname '.bin'
      else
         if ( _format = 'coards' )
           _fmsg = _format ' file ' _ofname '.nc'
         else
           _fmsg = _format ' file ' _ofname '.grb'
         endif
      endif


*     OK, ready for the real work
*     ---------------------------
      rc = lats4d ( args )
      if ( rc != 0 )
         say _myname 'error creating ' _fmsg
      else
         say _myname 'created ' _fmsg
      endif

      if ( _quit = 1 ) 
         say _myname 'exiting from GrADS...'
         'quit'
      endif


return rc


*............................................................................

*
* lats4d(args)  The main lats4d entry point.
*

function lats4d ( args )

  if ( _verb = 1 )
      say _myname 'Version ' _version 
  endif

* Open input file and reset dimension environment if desired...
* Assumes global grid, i.e., longitudinal wrapping
* -------------------------------------------------------------
  if ( _ifname != '' )
      'reinit'
      rc=openf(_ifname,_ftype)
      if ( rc != 0 )
         say _myname 'cannot open file "'_ifname'"'
         return rc
      endif
      rc = xyrange()
      'set x ' _xmin ' ' _xmax
      'set y ' _ymin ' ' _ymax
  endif

* Get file number
* ---------------
  'q file'
   if ( rc != 0 )
       say _myname 'no files open; specify "-i" option or open a file'
       return rc
   endif
   _datf = subwrd(result,2)


*  Open dimension environment file if so desired,
*  otherwise let dim env file be the same as the current data file
*  ---------------------------------------------------------------
   rc = setdimf()
   if ( rc != 0 )
       say _myname 'could not open dimension env file: ' _dimenv
       return rc
   endif


* File information
* ----------------
   if ( _verb = 1 )
        'q file ' _datf
        say _myname 'Data file is ' 
        say result
        if ( _datf != _dimf )
           'q file ' _dimf
           say _myname 'Dimension environment file is '
           say result
        else
           say _myname 'Dimension environment file same as data file' 
        endif
   endif


* Metadata describing origin of the dataset
*------------------------------------------
  'set lats model       ' _model
  'set lats center      ' _center

* Parameter table
* ---------------
  rc = getvars()
  ch = substr(_table,1,1)
*                                     generate parameter table
  if ( ch = '@' )
       _table = substr(_table,2,256)
       rc = mkptab()
       if ( rc != 0 ); return 1; endif
  endif
*                                     use internal table
  if ( ch = '=' )
       _table = ''
  endif
  if ( _table != '' )
      'set lats parmtab     ' _table
      id = subwrd(result,5)
      if ( id = 0 ); return 1; endif
  endif

* LATS comments (Title for GRIB files)
* ------------------------------------
  if ( _title = '' )
     if ( _func = '' )
       comment = 'File created by GrADS using LATS4D available from '
     else
       comment = 'File created by GrADS using LATS4D '
       comment = comment '(-func ' _func ') available from '
     endif
     url = 'http://dao.gsfc.nasa.gov/software/grads/lats4d/'
     'set lats comment ' comment url
  else
     'set lats comment ' _title
  endif


* Metadata describing conventions/format
* --------------------------------------
  if ( _format != 'stream' )
    'set lats convention ' _format
    if ( rc != 0 )
       say _myname 'not really, we stop right here'
       return 1
    endif
  endif

* Get (z,t) range on file
* -----------------------
  rc = ztrange()

* Metadata describing time frequency of grids
* -------------------------------------------
  if ( _freq = '' )
       rc = getfreq()
       if ( rc != 0 )
          say _myname 'could not determine time frequency; specify -freq option'
          return rc
       endif
  endif
  'set lats calendar  '  _cal
  if ( rc != 0 ); return rc; endif
  'set lats frequency '  _freq
  if ( rc != 0 ); return rc; endif
  if ( _tinc='' ); _tinc=1; endif
  odeltat = _tinc * _deltat
  'set lats deltat    '  odeltat
  if ( rc != 0 ); return rc; endif

* Redefine time range (_tmin,_tmax) if user wants to
* --------------------------------------------------
  rc = gettim()
  if ( rc = 1 )
       say _myname 'invalid time range ' _time
       return 
  endif
  if ( _verb = 1 )
      'q time'
      t1 = subwrd(result,3)
      t2 = subwrd(result,5)
      say _myname 'time range: ' t1 ' ' t2 ' by ' _tinc ', delta t: ' _deltat ' ' _freq
  endif
        
* Define the vertical dimension
* -----------------------------
  if ( _levels = '' )
      rc = getlevs()
  endif
  if ( _verb = 1 )
      say _myname 'vertical levels: ' _levels
  endif
  rc = setlevs()
  if ( rc = 1 ); return 1; endif

* Set LATS grid type
* ------------------
  'set lats gridtype ' _gridtype
  if ( rc != 0 ); return rc; endif

* Instruct LATS to use the dimension environment to set the LATS time
* -------------------------------------------------------------------
  'set lats timeoption dim_env'
  if ( rc != 0 ); return rc; endif

* Define the horizontal grid (based on dim env file)
* --------------------------------------------------
  'set dfile ' _dimf
  rc = setgrid()
  if ( rc != 0 ); return rc; endif
  if ( _verb = 1 )
      say _myname 'latitudinal  range: ' _lat
      say _myname 'longitudinal range: ' _lon
  endif
  _gid = getgid()
   if ( _gid < 1 ); return 1; endif
  'set dfile ' _datf


* Define the output file name
* ---------------------------
  if ( _format = 'stream' )
      'set fwrite '_ofname'.bin'
  else
      'set lats create ' _ofname
      _fid = subwrd(result,5)
       if ( _fid < 1 ); return 1; endif
  endif


* Define all the variables to be written to the LATS file
* -------------------------------------------------------
  rc = defvars()
  if ( rc != 0 ) 
       say _myname 'fatal error defining variables'
       return rc 
  endif
  if ( _verb = 1 )
      if ( _func  !='' ); say _myname 'Function expression: ' _func; endif
      if ( _svars !='' ); say _myname 'surface   variables: ' _svars ; endif
      if ( _uvars !='' ); say _myname 'upper air variables: ' _uvars ; endif
  endif
  if ( _svars='' & _uvars='' )
	say _myname 'All invalid sfc/upper variables: ' _vars
	return 1
  endif

*  Put GrADS in LATS output mode
*  -----------------------------
   if ( _format = 'stream' )
      'set gxout fwrite'
   else
      'set gxout latsdata'
      if ( rc != 0 ); return rc; endif
   endif

*  Now, if user wants time mean then redefine time environment
*  -----------------------------------------------------------
   if ( _mean = 1 )
      _tbeg = _tmin
      _tend = _tmax
      if ( _tinc = '' ); _tinc = 1; endif
      _tmid = ( _tmin + _tmax ) / 2
      if ( _tmean != '' ); _tmid = _tmean; endif
      _tmin = _tmid
      _tmax = _tmid
      if ( _verb = 1 )
         'set t ' _tmid
         'q time'
         t1 = subwrd(result,3)
         say _myname 'time average grid valid at ' t1
         say _myname 'averaging increment: ' _tinc ', delta t: ' _deltat ' ' _freq
      endif
   endif

* The LATS interface only allows one horizontal slice of data
* to be written at a time (z and t dimensions must be fixed),
* so we setup 2 loops
* --------------------------------------------------------------


* Loop over time...
* -----------------
  _t = _tmin
  while ( _t <= _tmax )

     'set t ' _t
    if ( rc != 0 ); return rc; endif
    'set z ' _zmin
    if ( rc != 0 ); return rc; endif

      if ( _verb = 1 )
        'q time'
        time = subwrd(result,3)
        say _myname 'writing to ' _fmsg ' on ' time
      endif

*    Write surface variables
*    -----------------------
     n = 1
     while ( n <= _nsvar )
        var = subwrd(_svars,n)
        if ( _dimenv != '' )
           var = var '.' _datf '(t=' _t ')'
        endif 
        var = subst(var,_func)
        if ( _mean = 1 )
           var = 'ave('var',t='_tbeg',t='_tend','_tinc')'
        endif
        vid = subwrd(_svids,n)
        if ( _format != 'stream' )
            'set lats write  ' _fid ' ' vid
            if ( rc != 0 ); return rc; endif
        endif
        'set dfile ' _dimf
        'display ' var
        if ( rc != 0 ); return rc; endif
        'set dfile ' _datf
        n = n + 1
     endwhile


*
*    Loop over upper air variables...
*    --------------------------------
     n = 1
     while ( n <= _nuvar )

         var = subwrd(_uvars,n)
         if ( _dimenv           != '' )
             var = var '.' _datf '(t=' _t ')'
         endif 
         var = subst(var,_func)
         if ( _mean = 1 )
              var = 'ave('var',t='_tbeg',t='_tend','_tinc')'
         endif
         vid = subwrd(_uvids,n)


*        Loop over levels...
*        -------------------
         k = 1
         _level = subwrd(_levels,k)
         while ( _level != '' )

              _latlevel = subwrd(_latlevels,k)

              'set lev ' _level

              if ( _format != 'stream' )
                   'set lats write  ' _fid ' ' vid ' ' _latlevel
                   if ( rc != 0 ); return rc; endif
              endif

              'set dfile ' _dimf
              'display ' var
              if ( rc != 0 ); return rc; endif
              'set dfile ' _datf

              k = k + 1
              _level = subwrd(_levels,k)

         endwhile

         n = n + 1

     endwhile

     _t = _t + _tinc

   endwhile

   rc = restdim()

*  Close LATS file
*  ---------------
  if ( _format = 'stream' )
     'disable fwrite'
  else
     'set lats close ' _fid
     if ( rc != 0 ); return rc; endif
  endif

* All done
* --------
  return 0

*.........................................................................

function usage(flag)

      say ''
      say 'NAME'
      say ''
      say '     lats4d - LATS for Dummies (Version ' _version ')'
      say ''
      say 'SYNOPSIS'
      say ''
      say '     lats4d  [-i fn] [-o fn] [-cal calendar] [-center ctr] [-de fn]'
      say '             [-format fmt] [-ftype ctl|sdf|xdf] [-freq ...] '
      say '             [-func expr] [-h] [-grid type]'
      say '             [-lat y1 y2] [-levs ...] [-lon x1 x2] '
      say '             [-model mod] [-mean] [-precision nbits] [-table tab] '
      say '             [-time t1 t2 [tincr]] [-title ...]'
      say '             [-v] [-vars ...] [-xvars] [-zrev] [-q] '
      say ''

      if ( flag != 1 )
                say '     Enter "lats4d -h" for additional information'
                say ''
                return
      endif

      say 'DESCRIPTION'
      say ''
      say '     A minimum fuss gs script for writing NetCDF, HDF-SDS or '
      say '     GRIB files from GrADS using the PCMDI LATS interface '
      say '     (http://www-pcmdi.llnl.gov).  This script can serve as a'
      say '     general purpose file conversion and subsetting utility.'
      say '     Any GrADS readable file (GrADS IEEE, GSFC Phoenix, GRIB, '
      say '     NetCDF or HDF-SDS) can be subset and converted to GRIB, '
      say '     NetCDF, HDF-SDS or flat binary (direct access) using a'
      say '     single command line.'
      say '     '
      say '     When invoked without arguments this script will create a' 
      say '     COARDS compliant NetCDF or HDF-SDS file named '
      say '     "grads.lats.nc", with all the contents of the default '
      say '     file (all variables, levels, times). The file name and '
      say '     several other attributes can be customized at the command' 
      say '     line, see OPTIONS below.'
      say '     '
      say '     NetCDF files are obtained by running this script under the'
      say '     executable "gradsnc".  HDF-SDS files can be produced with'
      say '     the "gradshdf" executable. Notice that the classic version'
      say '     of grads, "gradsc", does not include support for LATS and'
      say '     therefore cannot be used with lats4d. This script requires'
      say '     GrADS Version 1.7.beta.9 or later.' 
      say ''
      say 'OPTIONS'
      say ''
      say '     -i      fn              input file name; it can be any'
      say '                             of the following:'
      say '                             - an ASCII control (ctl) file'
      say '                               used for GRIB and IEEE files'
      say '                             - a binary NetCDF file/template' 
      say '                             - a binary HDF-SDS file/template'
      say '                             - an ASCII data descriptor file (ddf)'
      say '                               used for non-COARDS compliant'
      say '                               NetCDF/HDF-SDS files through'
      say '                               the "xdfopen" command'
      say '                             If the option "-ftype" is not'
      say '                             specified lats4d attempts to'
      say '                             determine the file type using'
      say '                             a heuristic algorithm.'
      say '                             NOTE:  When the option "-i" is '
      say '                             specified a GrADS "reinit" is '
      say '                             issued before the file is opened.'
      say '                             For NetCDF/HDF-SDS templates in'
      say '                             GrADS consult the SDFopen home'
      say '                             page listed under SEE ALSO'
      say ''
      say '     -o      fn              output (base) file name; default: '
      say '                             "grads.lats"'
      say ''
      say '     -cal    calendar        calendar type: "standard", "noleap", '
      say '                             "clim", or "climleap"; default: '
      say '                             "standard"'
      say ''
      say '     -center ctr             center, e.g., PCMDI, GSFC, NCEP, etc'
      say ''
      say '     -de     fn              Dimension environment file name;'
      say '                             defaut: same as "-i" argument.'
      say '                             This option is useful for using'
      say '                             lats4d with the user defined function'
      say '                             (udf) regrid2. See REGRIDDING below'
      say '                             for more information.'
      say ''
      say '     -format fmt             LATS file format: coards, grib,'
      say '                             grads_grib or stream; specify '
      say '                             "grads_grib" instead of "grib" for'
      say '                             getting ctl and gribmap files as well.'
      say '                             NOTE: The option "stream" creates'
      say '                             a flat binary file using the GrADS'
      say '                             command "set gxout fwrite" which'
      say '                             is not part of LATS.'                  
      say '                             '
      say ''
      say '    -ftype ctl|sdf|xdf       Specifies the input file type:'
      say '                             ctl  standard GrADS control (ctl)'
      say '                                  file used for IEEE and GRIB '
      say '                                  files'
      say '                             sdf  COARDS compliant NetCDF/HDF-SDS'
      say '                                  binary data file'
      say '                             xdf  data descriptor file (ddf)'
      say '                                  used for non-COARDS compliant'
      say '                                  NetCDF/HDF-SDS files through'
      say '                                  the "xdfopen" command'
      say '                             By default lats4d attempts to '
      say '                             determine the file type using a'
      say '                             heuristic algorithm; use this'
      say '                             option if lats4d fails to properly'
      say '                             detect the input file type'
      say '                             '
      say '     -freq  [n] unit         Time frequency of the input file.'
      say '                             LATS4D usually detects this from'
      say '                             the GrADS metadata, but sometimes'
      say '                             it fails with an error message.'
      say '                             In such cases use this option.'
      say '                             Example: -freq 6 hourly '
      say '                             NOTE: unlike GrADS, LATS does not'
      say '                             support time frequency in minutes'
      say '                             Default: n=1, e.g., -freq daily'
      say '                             '
      say '    -func expr               Evaluates the expression "expr"'
      say '                             before writing  to the output'
      say '                             file. The character "@" is used'
      say '                             to denote the variable name in'
      say '                             "expr". Example:'
      say '                             '
      say '                               -func ave(@,t-1,t+1)'
      say '                             '
      say '                             will replace "@" with each '
      say '                             variable name and produce a file'
      say '                             with running means. Default:'
      say '                             expr = @'
      say '                             '
      say '     -grid type              Grid type: linear, gaussian or'
      say '                             generic; default: linear'
      say '                             '
      say '     -h                      displays this man page'
      say ''
      say '     -lat    y1 y2           latitude range, e.g., "-30 30" for '
      say '                             30S thru 30N;  default: latitude '
      say '                             dimension environment'
      say ''
      say '     -levs   lev1 ... levN   list of levels; default: all levels'
      say ''
      say '     -lon    x1 x2           longitude range, e.g., "-50 20" for '
      say '                             50W thru 20E; default: longitude '
      say '                             dimension environment'
      say ''
      say '     -mean                   saves time mean to file; the actual'
      say '                             averaging period is specified with'
      say '                             the "-time" option; the "tincr" '
      say '                             parameter is the time increment'
      say '                             for the average (see GrADS ave()'
      say '                             function)'
      say ''
      say '     -model  mod             model name, e.g., GEOS/DAS'
      say ''
      say '     -precision nbits        specify the number of bits of'
      say '                             precision when storing in GRIB.'
      say '                             This option is only used when'
      say '                             lats4d automatically generates'
      say '                             a parameter table file (see option'
      say '                             -table below), and the output'
      say '                             format is "grib" or "grads_grib".'
      say '                             Default: nbits = 16'
      say ''
      say '     -table  tab             LATS parameter table file, e.g., '
      say '                             "dao.lats.table". If the table name'
      say '                             starts with "@" (e.g., @my.table)'
      say '                             then lats4d automatically generates'
      say '                             a LATS parameter table appropriate'
      say '                             for the current file and saves it '
      say '                             to a file; the file name in this'
      say '                             case is the same as "tab" with the'
      say '                             @ character removed (e.g., my.table).'
      say '                             Specify tab as "=" for using the'
      say '                             internal LATS parameter table.'
      say '                             See below for additional info on'
      say '                             parameter tables.'
      say '                             Default: @.grads.lats.table'
      say ''
      say ''
      say '     -time   t1 t2 [tincr]   time range and time increment in '
      say '                             units of the "delta t" in the'
      say '                             input file; "tincr" is optional;'
      say '                             Example: "0z1dec91 18z31dec91 2"'
      say '                                       to write every other '
      say '                                       time step' 
      say '                             Defaults: (t1,t2) is taken from '
      say '                             the time dimension environment,' 
      say '                             and tincr=1. Note: enter "= ="'
      say '                             for keeping the default values'
      say '                             for (t1,t2) while specifying tincr'
      say ''
      say '     -title text             output dataset TITLE for GRIB files,'
      say '                             COMMENTS for NetCDF/HDF files'
      say ''
      say '     -v                      verbose mode'
      say ''
      say '     -vars   var1 ... varN   list of variables; default: all '
      say '                             variables on the current file will'
      say '                             be written to the output file'
      say ''
      say '     -xsfc                   exclude all surface variables'
      say ''
      say '     -xvars  var1 ... varN   list of variables to exclude;'
      say '                             default: none'
      say ''
      say '     -xupper                 exclude all upper air variables'
      say ''
      say '     -zrev                   reverse order of vertical levels'
      say ''
      say '     -q                      quits GrADS upon return'
      say '     '
      say 'LATS PARAMETER TABLES'
      say '     '
      say '     LATS maintains an internal parameter table that prescribes'
      say '     variable names, description, units, datatype, basic'
      say '     structure (e.g., upper air or surface), and compression'
      say '     (GRIB options). These descriptors are inferred from the'
      say '     parameter name only, and thus most of the metadata needed'
      say '     to write GRIB and/or netCDF data are located in the'
      say '     parameter table and need not be specified at the command'
      say '     line. The option "-table" is provided to override the '
      say '     internal table with an external parameter file. For'
      say '     additional information on LATS parameter tables'
      say '     consult http://www-pcmdi.llnl.gov/software/lats/.'
      say ''
      say '     The only inconvenience of this approach is that variables'
      say '     names being written to file must match those defined in '
      say '     this internal parameter table (which is essentially the '
      say '     same as the "AMIPS2" LATS table, see URL above).'
      say '     To circumvent this problem lats4d can automatically'
      say '     generate a parameter table based on the current file'
      say '     metadata. Since GrADS keeps no units or GRIB packing'
      say '     information, this parameter file sets the units entry'
      say '     to blank and uses defaults for the GRIB packing parameters.'
      say '     The default GRIB packing algorithm is "16-bit fixed width '
      say '     compression" and produces GRIB files which are about half'
      say '     the size of NetCDF/HDF-SDS files. The option "-precision"'
      say '     allows the user to define the number of bits of precision'
      say '     at the command line; see EXAMPLES ex2a,b,c below.' 
      say '     If you care about having proper metadata written to'
      say '     your file or need more efficient GRIB packing then you can '
      say '     either change your variable names to match those in the '
      say '     internal LATS table, or customize an existing LATS parameter'
      say '     table; see URL above for sample parameter tables.'
      say '     '
      say 'LATS QUALITY CONTROL WARNINGS'
      say '     '
      say '     Quality control (QC) information is included in some '
      say '     LATS parameter tables to help the user ensure that their'
      say '     data is being written properly. In such cases, if LATS'
      say '     detects suspect data it writes a warning message to the'
      say '     screen and saves additional information in a log file.'
      say '     Consult the LATS home page for additional information.'
      say '     '
      say 'REGRIDDING'
      say '     '
      say '     This script can be used with Mike Fiorino s user'
      say '     defined function (udf) regrid2(). This combination'
      say '     allows you to convert any GrADS redable file to any'
      say '     other horizontal resolution/domain of your choice. '
      say '     Here is a quick roadmap:'
      say '     1. Start by installing regrid2() available from'
      say '        ftp://grads.iges.org/grads/sprite/udf/regrid2beta.tar'
      say '     2. If you already have a sample file at the desired new'
      say '        resolution, great! Otherwise you can get one by creating '
      say '        a fake GrADS control file. There are a few samples'
      say '        on the last4d home page: geos1x1.ctl, geos4x5.ctl and'
      say '        geos2x25.ctl. This file is used to define the dimension'
      say '        environment at the new desired resolution through the'
      say '        "-de" option.'
      say '     3. Here is an example which converts the sample model.???'
      say '        data file from 4x5 (latxlon) resolution to 1x1:'
      say '     '
      say '     lats4d -i model -de geos1x1 -func regrid2(@,1,1,bs_p1,-180,-90)'
      say '     '
      say '       The resulting "grads.lats.nc" file is at 1x1 degree'
      say '       resolution.'
      say '     '
      say 'EXAMPLES'
      say ''
      say '     Download files "model.ctl", "model.gmp" and "model.grb"'
      say '     from http://dao.gsfc.nasa.gov/software/grads/lats4d/.'
      say '     Then start "gradsnc" or "gradshdf" and try these,'
      say '     carefully examining the files produced:'
      say '     '
      say '     lats4d -h'
      say '     lats4d -v -q -i model -o ex1 '
      say '     lats4d -v -q -i model -o ex2a -format grads_grib'
      say '     lats4d -v -q -i model -o ex2b -format grads_grib -precision 8'
*     say '     lats4d -v -q -i model -o ex2c -format grads_grib -precision 32'
      say '     lats4d -v -q -i model -o ex3 -levs 700 500 -vars ua va'
      say '     lats4d -v -q -i model -o ex4 -time 1jan1987 3jan1987'
      say '     lats4d -v -q -i model -o ex5 -time = = 2'
      say '     lats4d -v -q -i model -o ex6 -mean'
      say '     lats4d -v -q -i model -o ex7 -mean -time = = 2'
      say '     lats4d -v -q -i model -o ex8 -lat 20 70 -lon -140 -60'
      say ''
      say '     Note: the "-q" option above is to make sure you'
      say '           restart GrADS; see BUGS below. You may want to'
      say '           enter these from your OS shell, e.g.,'
      say ''
      say '     % gradsnc -blc "lats4d -v -q -i model -o ex1"'
      say ''
      say '    The sh(1) script "lats4d" allows you to enter lats4d'
      say '    options directly from the Unix command line, e.g.,'
      say ''
      say '    % lats4d -v -i model -o ex1 '
      say ''
      say 'BUGS'
      say ''
      say '     Sometimes lats4d will only work if you exit and'
      say '     restart GrADS.'
      say ''
      say '     The option "-precision 32" does not quite work. This'
      say '     appears to be a LATS bug.'
      say ''
      say '     Because of a limitation in the GRIB format, "grib" or '
      say '     "grads_grib" output cannot have levels where p<1.'
      say '     To circumvent this problem, a hybrid level number is'
      say '     is used in such cases.'
      say ''
      say 'SEE ALSO    '
      say ''
      say '     GrADS   http://grads.iges.org/grads/   '
      say '     LATS    http://www-pcmdi.llnl.gov/software/lats'
      say '     LATS4D  http://dao.gsfc.nasa.gov/software/grads/lats4d'
      say '     SDFopen http://www.cdc.noaa.gov/~hoop/grads.html'
      say '     XDFopen http://www.cdc.noaa.gov/~hoop/xdfopen.shtml'
      say '     NetCDF  http://www.unidata.ucar.edu/packages/netcdf/'
      say '     HDF     http://hdf.ncsa.uiuc.edu/'
      say '     GRIB    ftp://ncardata.ucar.edu/docs/grib/prev-vers/guide.txt'
      say '             http://www.wmo.ch/web/www/reports/Guide-binary-2.html'
      say '     '
      say 'COPYRIGHT'
      say ''
      say '     Copyright (c) 1998-1999 A. da Silva'
      say '     Permission is granted to copy, modify and distribute this'
      say '     software provided it is not sold for profit, and provided '
      say '     this notice is included. '
      say ''
      say 'NO WARRANTY'
      say '     '
      say '     Because lats4d is provided free of charge, it is provided'
      say '     "as is" WITHOUT WARRANTY OF ANY KIND, either expressed or'
      say '     implied, including, but not limited to, the implied'
      say '     warranties of merchantability and fitness for a particular'
      say '     purpose. USE AT YOUR OWN RISK.'
      say '     '
      say 'CREDITS'
      say ''
      say '    Arlindo da Silva (NASA/GSFC) wrote the lats4d.gs script. '
      say '    Mike Fiorino (PCMDI/LLNL) wrote the LATS interface to'
      say '    GrADS. Robert Drach, Mike Fiorino and Peter Gleckler'
      say '    (PCMDI/LLNL) wrote the LATS library.'
      say ''

return 1

*.........................................................................

*
* parsecmd() Parse command line arguments
*

function parsecmd(args)

*
*     Note: customize defaults for your site
*
      _ifname = ''
      _ofname = 'grads.lats'
      _format = 'coards'
      _model  = 'geos/das'
      _center = 'gsfc'
      _table  = '@.grads.lats.table'
      _precision = 16
      _vars = ''
      _xvars = ''
      _xsfc = ''
      _xupper = ''
      _func = ''
      _dimenv = ''
      _title = ''

      _lat = ''
      _lon = ''
      _levels = ''
      _time = ''
      _tmean = ''
      _tinc = ''
      _freq = ''
      _ftype  = ''

      _gridtype = 'linear'

      _help = 0
      _verb = 0
      _quit = 0
      _mean = 0
      _zrev = 0

      _cal  = 'standard'

      options = '-model -center -table -v -q -i -o -format -vars -levs -time'
      options = options ' -h -help -lat -lon -cal -mean -precision -ftype -tmean'
      options = options ' -grid -func -zrev -de -title -xvars -xsfc -xupper'
      options = options ' -geos1x1 -geos4x5 -geos2x25'
      options = options ' -geos1x1a -geos4x5a -geos2x25a'
      options = options ' -geos1x1s -geos4x5s -geos2x25s'

      i = 1
      token = subwrd(args,i)
      while ( token != '' )

*        Handle each option separately ...

         if ( token = '-help' | token = '-h' )
              _help = 1
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-v' )
              _verb = 1
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-q' )
              _quit = 1
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-mean' )
              _mean = 1
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-xsfc' )
              _xsfc = 1
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-xupper' )
              _xupper = 1
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-zrev' )
              _zrev = 1
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-i' )
              i = i + 1
              token = subwrd(args,i)
              first = x '' substr(token,1,1)
              while ( first != 'x-' & first != 'x' )
                _ifname= _ifname ' ' token
                 i = i + 1
                 token = subwrd(args,i)
                 first = x '' substr(token,1,1)
              endwhile
*                                 allows sdfopen templates
              if ( critique(1,3,'-i',_ifname) = 1 ); return 1; endif
         endif

         if ( token = '-o' )
              _ofname = subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-o',_ofname) = 1 ); return 1; endif
         endif

         if ( token = '-cal' )
              _cal = subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-cal',_cal) = 1 ); return 1; endif
         endif

         if ( token = '-format' )
              _format = subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-format',_format) = 1 ); return 1; endif
              if ( _format = 'fwrite' ); _format = 'stream'; endif
         endif

         if ( token = '-grid' )
              _gridtype = subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-grid',_gridtype) = 1 ); return 1; endif
         endif

         if ( token = '-ftype' )
              _ftype= subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-ftype',_format) = 1 ); return 1; endif
              if ( _ftype='nc'|_ftype='netcdf'|_ftype='hdf'|_ftype='hdf-sds'|_ftype='gfio' )
                   _ftype = 'sdf'
              endif
              if ( _ftype='grib'|_ftype='grads_grib'|_ftype='ieee'|_ftype='sequential')
                   _ftype = 'ctl'
              endif
         endif

         if ( token = '-model' )
              _model = subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-model',_model) = 1 ); return 1; endif
         endif

         if ( token = '-center' )
              _center = subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-center',_center) = 1 ); return 1; endif
         endif

         if ( token = '-de' )
              _dimenv = subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-de',_dimenv) = 1 ); return 1; endif
         endif

         if ( token = '-table' )
              _table = subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-table',_table) = 1 ); return 1; endif
         endif

         if ( token = '-precision' )
              _precision = subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-precision',_precision) = 1 ); return 1; endif
         endif

         if ( token = '-func' )
              _func= subwrd(args,i+1) 
              i = i + 2
              token = subwrd(args,i)
              if ( critique(1,1,'-func',_func) = 1 ); return 1; endif
         endif

         if ( token = '-vars' )
              i = i + 1
              token = subwrd(args,i)
              first = x '' substr(token,1,1)
              while ( first != 'x-' & first != 'x' )
                _vars = _vars ' ' token
                 i = i + 1
                 token = subwrd(args,i)
                 first = x '' substr(token,1,1)
              endwhile
              if ( critique(1,999,'-vars',_vars) = 1 ); return 1; endif
         endif

         if ( token = '-xvars' )
              i = i + 1
              token = subwrd(args,i)
              first = x '' substr(token,1,1)
              while ( first != 'x-' & first != 'x' )
                _xvars = _xvars ' ' token
                 i = i + 1
                 token = subwrd(args,i)
                 first = x '' substr(token,1,1)
              endwhile
              if ( critique(1,999,'-xvars',_xvars) = 1 ); return 1; endif
         endif

         if ( token = '-title' )
              i = i + 1
              token = subwrd(args,i)
              first = x '' substr(token,1,1)
              while ( first != 'x-' & first != 'x' )
                _title = _title ' ' token
                 i = i + 1
                 token = subwrd(args,i)
                 first = x '' substr(token,1,1)
              endwhile
              if ( critique(1,999,'-title',_title) = 1 ); return 1; endif
         endif

*
*        Note: we force -lon/-lat to have 2 args because lats do
*              not allow single point grids. 
*
         if ( token = '-lat' )
              tok1 = subwrd(args,i+1)
              tok2 = subwrd(args,i+2)
              _lat = tok1 ' ' tok2
              i = i + 3
              token = subwrd(args,i)
              if ( critique(2,2,'-lat',_lat) = 1 ); return 1; endif
         endif

         if ( token = '-lon' )
              tok1 = subwrd(args,i+1)
              tok2 = subwrd(args,i+2)
              _lon = tok1 ' ' tok2
              i = i + 3
              token = subwrd(args,i)
              if ( critique(2,2,'-lon',_lon) = 1 ); return 1; endif
         endif

         if ( token = '-levs' )
              i = i + 1
              token = subwrd(args,i)
              first = x '' substr(token,1,1)
              while ( first != 'x-' & first != 'x' )
                _levels = _levels ' ' token
                 i = i + 1
                 token = subwrd(args,i)
                 first = x '' substr(token,1,1)
              endwhile
         endif

         if ( token = '-time' )
              i = i + 1
              token = subwrd(args,i)
              first = x '' substr(token,1,1)
              while ( first != 'x-' & first != 'x' )
                _time = _time ' ' token
                 i = i + 1
                 token = subwrd(args,i)
                 first = x '' substr(token,1,1)
              endwhile
              if ( critique(1,3,'-time',_time) = 1 ); return 1; endif
              _tinc = subwrd(_time,3)
              t1 = subwrd(_time,1)
              t2 = subwrd(_time,2)
              if ( _tinc != '' )
                   _time = t1 ' ' t2
              else
                   _tinc = 1
              endif
              if ( t1 = '=' | t2 = '=' ); _time = ''; endif
              if ( _tinc < 1 )
                   rc = usage()
                   say _myname 'invalid time increment ' _tinc
                   return 1
              endif
         endif

         if ( token = '-freq' )
              _deltat = subwrd(args,i+1)
              _freq   = subwrd(args,i+2)
              token = _freq ' ' _deltat
              if ( critique(1,2,'-freq',token) = 1 ); return 1; endif
              first = substr(_freq,1,1)
              if ( first = '-' | first = '' )
                   _freq = _deltat
                   _deltat = 1
                   i = i + 2
              else 
                   i = i + 3
              endif
              token = subwrd(args,i)
         endif

*        Undocumented options (DAO specific)
*        -----------------------------------
         if ( token = '-geos1x1' )
              rc = mkgeosf('geos1x1.ctl',1,1,360,181,'bl_p1')
              if ( rc != 0 ); return 1; endif
              i = i + 1
              token = subwrd(args,i)
         endif
         if ( token = '-geos1x1a' )
              rc = mkgeosf('geos1x1.ctl',1,1,360,181,'ba_p1')
              if ( rc != 0 ); return 1; endif
              i = i + 1
              token = subwrd(args,i)
         endif
         if ( token = '-geos1x1s' )
              rc = mkgeosf('geos1x1.ctl',1,1,360,181,'bs_p1')
              if ( rc != 0 ); return 1; endif
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-geos2x2.5' )
              rc = mkgeosf('geos2x25.ctl',2.5,2,144,91,'bl_p1')
              if ( rc != 0 ); return 1; endif
              i = i + 1
              token = subwrd(args,i)
         endif
         if ( token = '-geos2x2.5a' )
              rc = mkgeosf('geos2x25.ctl',2.5,2,144,91,'ba_p1')
              if ( rc != 0 ); return 1; endif
              i = i + 1
              token = subwrd(args,i)
         endif
         if ( token = '-geos2x2.5s' )
              rc = mkgeosf('geos2x25.ctl',2.5,2,144,91,'bs_p1')
              if ( rc != 0 ); return 1; endif
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-geos4x5' )
              rc = mkgeosf('geos4x5.ctl',5,4,72,46,'bl_p1')
              if ( rc != 0 ); return 1; endif
              i = i + 1
              token = subwrd(args,i)
         endif
         if ( token = '-geos4x5a' )
              rc = mkgeosf('geos4x5.ctl',5,4,72,46,'ba_p1')
              if ( rc != 0 ); return 1; endif
              i = i + 1
              token = subwrd(args,i)
         endif
         if ( token = '-geos4x5s' )
              rc = mkgeosf('geos4x5.ctl',5,4,72,46,'bs_p1')
              if ( rc != 0 ); return 1; endif
              i = i + 1
              token = subwrd(args,i)
         endif

         if ( token = '-tmean' )
              _tmean =  subwrd(args,i+1)
              i = i + 2
              token = subwrd(args,i)
         endif

*        Validate option
*        ---------------
         j = 1
         opt = subwrd(options,j)
         valid = 0
         while ( opt != '' )
            if ( token = opt ); valid = 1; endif
            j = j + 1
            opt = subwrd(options,j)
         endwhile

         if ( valid = 0 & token != '' )
            rc = usage()
            say _myname 'invalid option "'token'"'
            return 1
         endif

      endwhile

      if ( _help = 1 ) 
           rc = usage(1)
           return 1
      endif


return 0

function printopt()

  say '   Input  fname: ' _ifname
  say '   Output fname: ' _ofname
  say '         Format: ' _format
  say '      Variables: ' _vars
  say '     Time Range: ' _time
  say '         Levels: ' _levels

return

*.......................................................................

*
* openf(fname,ftype)  Opens a file according to the file type (ftype).
*                     If ftype is not specified it attempts to determine it
*                     by a heuristic algorithm;
*                     ftype can be 'ctl', 'sdf' or 'xdf'
*

function openf(fname,ftype)

*  Determine file type
*  -------------------
   if ( ftype = '' | ftype ='ftype' )
*                                       fname may be a template...
*                                       filen is always a file name
    filen = subwrd(fname,1)
    http = substr(filen,1,7)
     if ( http = 'http://' )
        ftype = 'sdf'
 
    else

      buf = read(filen)
      rc = sublin(buf,1)
      if ( rc != 0 )
           buf = read(filen'.ctl')
           rc = sublin(buf,1)
      endif
      if ( rc != 0 )
           say _myname 'cannot read file ' filen ' or ' filen '.ctl'
           return rc
      endif
      rec = sublin(buf,2)
      tok = subwrd(rec,1)
      if ( tok = 'dset' | tok='DSET' )
         is_xdf = 0
         i = 1
         tok = substr(filen,i,4)
         while ( tok != '' )
           if ( tok='.ddf' | tok='.DDF' ); is_xdr = 1; endif
           i = i + 1
           tok = substr(filen,i,4)
         endwhile
         if ( is_xdr = 1 )
              ftype = 'xdf'
         else
              ftype = 'ctl'
         endif
      else
         ftype = 'sdf'
      endif
   
    endif

   endif   

*  Open according to file type
*  ---------------------------
   if ( ftype = 'ctl' )
        'open ' fname
        _result = result
        return rc
   endif
   if ( ftype = 'sdf' ) 
        'sdfopen ' fname
        _result = result
        return rc
   endif
   if ( ftype = 'xdf' ) 
        'xdfopen ' fname
        _result = result
        return rc
   endif

   say _myname 'cannot handle file type "' ftype '"'

return 1

*
* setdimf - Open dimension env file 
* 

function setdimf()

   _dimf = 0
   if ( _dimenv != '' )
      rc=openf(_dimenv,'')
      if ( rc != 0 )
         say _myname 'cannot open file "'_dimenv'"'
         return rc
      endif
      result = _result
      i = 0
      while ( token != '' )
          i = i + 1
          lin = sublin(result,i)
          token = subwrd(lin,1)
          word  = subwrd(lin,2)
          if ( token = 'Data' & word = 'file' )
               _dimf = subwrd(lin,8)
          endif
      endwhile
      if ( _dimf = 0 ); return 1; endif
      'set dfile ' _dimf
      rc = xyrange()
      'set x ' _xmin ' ' _xmax
      'set y ' _ymin ' ' _ymax
      'set dfile ' _datf
   else
     _dimf = _datf
   endif

return 0


*
* gettim()  Redefines time range according to user defined _time
*
function gettim()

      if ( _time = '' ); return 0; endif

      tbeg = subwrd(_time,1)
      tend = subwrd(_time,2)

      'set time ' tbeg ' ' tend
      if ( rc = 1 ); return 1; endif

      rc = savedim()

      if ( _t1save < _tmin | _t2save > _tmax )
           return 1
      else
           _tmin = _t1save
           _tmax = _t2save
      endif

return 0
  



*
* getfreq() Examines two subsequent times and determine
*           time increment. This is heuristic and not
*           guaranteed to always work.
*

function getfreq()

*     Special case: only 1 time step on file

      if ( _tmax = 1 )
           _freq   = hourly
           _deltat = 6
           say _myname" only 1 time step on file, assuming frequency is " _deltat ' ' _freq
           return 0
      endif         

*
*     Get year-month-day-hour for 2 time intervals
*
      
      'q time'
      tb = subwrd(result,3)
      te = subwrd(result,5)
      'set t 1 2'
      'q time'
      t1 = subwrd(result,3)
      t2 = subwrd(result,5)
      'set time ' tb ' ' te
      rc = split( t1 )
      y1=_y; m1=_m; d1=_d; h1=_h 
      rc = split( t2 )
      y2=_y; m2=_m; d2=_d; h2=_h 
*
*     Minutes: not supported
*
      if ( y1=y2 & m1=m2 & d1=d2 & h1=h2 )
           say _myname 'LATS does not support time frequency in minutes'
           return 1
      endif
*
*     Hourly
*      
      if ( y1=y2 & m1=m2 & d1=d2 )
           _freq   = hourly
           _deltat = h2 - h1
           if ( _deltat <= 0 ); return 1; endif
           return 0
      endif
*
*     Daily
*      
      if ( y1=y2 & m1=m2 )
           _freq   = daily
           _deltat = d2 - d1
           if ( _deltat <= 0 ); return 1; endif
           return 0
      endif
*
*     Monthly
*      
      if ( y1=y2 )
           _freq   = monthly
           _deltat = m2 - m1
           if ( _deltat <= 0 ); return 1; endif
           return 0
      endif
*
*     Yearly
*
      if ( y1 < y2 )
           _freq   = yearly
           _deltat = y2 - y1
           return 0
      endif

return 1

*
* split()  Returns year, month, day & hour from date of the form
* 00Z02JAN1987
*
function split ( t )

   ch  = substr(t,3,1)
   if ( ch = ':' )
        off = 3
   else
        off = 0
   endif
   _h  = substr(t,1,2)
   _d  = substr(t,4+off,2)
   m3  = substr(t,6+off,3)
   _y  = substr(t,9+off,4)
   mons = 'JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC'
   i = 1
   while(i<=12)
      mon = subwrd(mons,i)
      if ( m3 = mon ) 
           _m = i
           return
      endif
      i = i + 1
   endwhile

return

*
* ztrange() Get (z,t) range on file
*
function ztrange()

      'q file'
      tmp = sublin ( result, 5 )
      _zmin = 1
      _zmax = subwrd(tmp,9)
      _tmin = 1
      _tmax = subwrd(tmp,12)

return

function xyrange()

      'q file'
      tmp = sublin ( result, 5 )
      _xmin = 1
      _xmax = subwrd(tmp,3)
      _ymin = 1
      _ymax = subwrd(tmp,6)

return


*
* getlevs() Get all levels from dimension environment
*
function getlevs()

      _levels = ''
      if ( _zrev = 1 )

         z = _zmax
         while ( z >= _zmin )
            'set z ' z
            lev = subwrd(result,4)
            _levels = _levels ' ' lev
             z = z - 1
         endwhile

      else

         z = _zmin
         while ( z <= _zmax )
            'set z ' z
            lev = subwrd(result,4)
            _levels = _levels ' ' lev
             z = z + 1
         endwhile

     endif

return

*
* setlevs() Set LATS vertical dimension. 
*
function setlevs()


* GRIB does not allow pressure < 1 hPa
* ------------------------------------
  plev = 1
  zlevs = ''
  if ( _format = 'grads_grib' | _format = 'grib' )

     k = 1
     lev = subwrd(_levels,k)
     while ( lev != '' )

         if  ( lev < 1 )
               say _myname 'invalid plev ' lev ' for GRIB output'
               plev = 0
         endif

         'set lev ' lev
         'q dims'
         tmp   = sublin(result,4)
         zlevs = zlevs ' ' subwrd(tmp,9)

         k = k + 1
        lev = subwrd(_levels,k)

     endwhile

   endif


*  Set pressure levels
*  -------------------
   if ( plev = 1 )

        if ( _verb = 1 )
             say _myname 'using PRESSURE for vertical coordinate'
        endif
        _latlevels = _levels
       'set lats vertdim plev ' _latlevels
       _lid = subwrd(result,5)
       if ( _lid < 1 ); return 1; endif
       _sid = 0

*  Set hybrid levels
*  -----------------
   else

        if ( _verb = 1 )
             say _myname 'using HYBRID level number for vertical coordinate'
        endif
        _latlevels = zlevs
       'set lats vertdim hybrid ' _latlevels
       _lid = subwrd(result,5)
       if ( _lid < 1 ); return 1; endif
       _sid = 0

   endif       


return 0


*
* setgrid()   Sets the horizontal grid if user specified 
*              -lat or -lon at the command line
*

function setgrid()

   rc = savedim()
   if ( _lon != '' )
      'set lon ' _lon
      if ( rc!=0 ); return rc; endif
   else
      _lon = _lonmin ' ' _lonmax
   endif

   if ( _lat != '' )
      'set lat ' _lat
      if ( rc!=0 ); return rc; endif
   else
      _lat = _latmin ' ' _latmax
   endif

return 0

*
* getgid()  Get horizontal grid id
*
function getgid()

  'q file'
  tmp = sublin(result,7)
  var = subwrd(tmp,1)

  rc = savedim()
  'set z 1'
  'set t 1'
  'set gxout latsgrid' 
  'd ' var
  bufr = sublin(result,1)
  grid = subwrd(bufr,2)
  if ( grid != 'GRID' )
       bufr = sublin(result,2)
       grid = subwrd(bufr,2)
       if ( grid != 'GRID' ); return -1; endif
  endif
  gid = subwrd(bufr,5)
  rc = restdim()

return gid

*
* savedim()  Save dimension environment
*
function savedim()

   'q dims'

   tmp = sublin(result,2)
   type = subwrd(tmp,3)
   if ( type = 'varying' )
      _x1save = subwrd(tmp,11)
      _x2save = subwrd(tmp,13)
      _lonmin = subwrd(tmp,6)
      _lonmax = subwrd(tmp,8)
   else
      _x1save = subwrd(tmp,9)
      _x2save = _x1save
      _lonmin = subwrd(tmp,6)
      _lonmax = _lonmin
   endif

   tmp = sublin(result,3)
   type = subwrd(tmp,3)
   if ( type = 'varying' )
      _y1save = subwrd(tmp,11)
      _y2save = subwrd(tmp,13)
      _latmin = subwrd(tmp,6)
      _latmax = subwrd(tmp,8)
   else
      _y1save = subwrd(tmp,9)
      _y2save = _y1save
      _latmin = subwrd(tmp,6)
      _latmax = _latmin
   endif

   tmp = sublin(result,4)
   type = subwrd(tmp,3)
   if ( type = 'varying' )
      _z1save = subwrd(tmp,11)
      _z2save = subwrd(tmp,13)
   else
      _z1save = subwrd(tmp,9)
      _z2save = _z1save
   endif

   tmp = sublin(result,5)
   type = subwrd(tmp,3)
   if ( type = 'varying' )
      _t1save = subwrd(tmp,11)
      _t2save = subwrd(tmp,13)
   else
      _t1save = subwrd(tmp,9)
      _t2save = _t1save
   endif

return

*
* restdim()  Restore saved dimension environment
*
function restdim()

     if ( _x1save = '_x1save' )
         say 'restdim: dimensions not saved'
         return 1
     endif

     'set x  ' _x1save ' ' _x2save
     'set y  ' _y1save ' ' _y2save
     'set z  ' _z1save ' ' _z2save
     'set t  ' _t1save ' ' _t2save

return

*
* getvars()  Get all variables from current file
*
function getvars()

  'q file'
  tmp = sublin(result,6)
  nvars = subwrd(tmp,5)

  n = 1
  _svars = ''; _slong=''; _nsvar = 0
  _uvars = ''; _ulong=''; _nuvar = 0
  while ( n <= nvars )

     tmp = sublin(result,6+n)
     var = subwrd(tmp,1)
     long = ''
     i = 4
     token = subwrd(tmp,i)
     while ( token != '' )
         long = long ' ' token
         i = i + 1
         token = subwrd(tmp,i)
     endwhile
     nlevs = subwrd(tmp,2)
     rc = validate(var,_vars)
     if ( rc = 1 )
          rc = xvalidat(var,_xvars)
          if ( rc = 0 )
             if ( _verb = 1 )
                  say _myname 'excluding variable ' var
             endif
          endif
     endif
     if ( rc = 1 )
      if ( nlevs < 2 )
        _nsvar = _nsvar + 1
        _svars = _svars ' ' var
        j = _nsvar; _slong.j = long
        if ( _verb = 1 )
***        say _myname 'selecting sfc var ' var ': ' _slong.j
        endif
      else
        _nuvar = _nuvar + 1
        _uvars = _uvars ' ' var
        j = _nuvar; _ulong.j = long
        if ( _verb = 1 )
***        say _myname 'selecting up   var ' var ': ' _ulong.j
        endif
      endif
     endif

     n = n + 1

  endwhile

* Exclude all surface variables if user so chooses
* ------------------------------------------------
  if ( _xsfc = 1 )
      _nsvar = 0
      _svars = ''
  endif

* Exclude all upper air variables if user so chooses
* --------------------------------------------------
  if ( _xupper = 1 )
      _nuvar = 0
      _uvars = ''
  endif

return

*
* mkptab()  Creates a LATS parameter table from variable list.
*

function mkptab()

*
* Temporary file
*

fname = _table

if ( _verb = 1 )
      say _myname 'creating LATS PARAMETER TABLE file ' fname
endif

*
* Identify ourselves...
*

rc=write(fname,'# LATS PARAMETER TABLE automatically generated by lats4d')
rc=subwrd(rc,1)
if ( rc > 0 )
   say rc
   say _myname 'cannot create LATS table file ' fname ' on current directory'
   return rc
endif

rc=write(fname,'# Note: This table lacks the UNITS and GRIB decimal_scale_factor columns.', append)
rc=write(fname,'#       It also does not include QC information. For additional', append)
rc=write(fname,'#       information on LATS parameter tables consult:'      , append)
rc=write(fname,'#           http://www-pcmdi.llnl.gov/software/lats/'       , append)

*
* General info
*

rc=write(fname,'#---------------------------------------------------------------------------------------------------',append)
rc=write(fname,'#',append)
#rc=write(fname,'# A parameter file is divided into sections, indicated by #! comments.',append) 
rc=write(fname,'# The sections may appear in any order. The 'center' section is only required for GRIB output.',append)
rc=write(fname,'#',append)
rc=write(fname,'# #!variable',append)
rc=write(fname,'#',append)
rc=write(fname,'#   Variable table: defines variable-specific parameters',append)
rc=write(fname,'#',append)
rc=write(fname,'# #!vert',append)
rc=write(fname,'#',append)
rc=write(fname,'#   Vertical dimension type table: defines categories of vertical dimensions',append)
rc=write(fname,'#',append)
rc=write(fname,'# #!center',append)
rc=write(fname,'#',append)
rc=write(fname,'#   Center table: defines GRIB parameters which identify the originating process, center, and subcenter.',append)
rc=write(fname,'#',append)
rc=write(fname,'# #!qc',append)
rc=write(fname,'#',append)
rc=write(fname,'#   Quality control marks table: defines the values controlling the quality control routines.',append)
rc=write(fname,'# ',append)
rc=write(fname,'#---------------------------------------------------------------------------------------------------',append)

*
* Variable list
* 

rc=write(fname,'#!variable',append)
rc=write(fname,'#',append)
rc=write(fname,'# Variables',append)
rc=write(fname,'#   (max number of entries = LATS_MAX_PARMS in lats.h)',append)
rc=write(fname,'#',append)
rc=write(fname,'# The format of each record is:',append)
rc=write(fname,'#   name | id | title | units | datatype | surface | decimal_scale_factor | precision | comments_1 | comments_2',append)
rc=write(fname,'#',append)
rc=write(fname,'# name = variable name (no blanks)',append)
rc=write(fname,'# id = GRIB parameter number (>127 => AMIP-2 specific)',append)
rc=write(fname,'# title = long name (description)',append)
rc=write(fname,'# units = variable units',append)
rc=write(fname,'# datatype = 'float' or 'int'',append)
rc=write(fname,'# level_type = level_type in vertical dimension table, or blank if values must be defined via lats_vert_dim',append)
rc=write(fname,'# decimal_scale_factor = GRIB decimal scale factor, or -999 if no decimal scaling',append)
rc=write(fname,'# precision = number of bits of precision if stored in GRIB,',append)
rc=write(fname,'#             or -999 for level-dependent bit length (ignored if decimal_scale_factor is set)',append)
rc=write(fname,'# comments_1 = comments, ignored by LATS',append)
rc=write(fname,'# comments_2 = comments, ignored by LATS',append)
rc=write(fname,'#',append)
rc=write(fname,'#---------------------------------------------------------------------------------------------------',append)


* Fake these for now
units = '    '
scale_factor = -999
precision = _precision
dattype = 'float'

* Sfc variables
id = 1
i = 1
levtype = 'sfc'
while(i <= _nsvar )
      name   = subwrd(_svars,i) '        '
      title  = _slong.i '                                                                     '
      name   = substr(name,1,8)
      title  = substr(title,1,70) 
      idc = substr(id' ',1,2)
      record = name '| ' idc ' |' title ' | ' units ' | ' dattype ' | ' levtype ' | ' scale_factor ' | ' precision ' |   |   |' 
      rc = write(fname,record,append)
      i = i + 1
      id = id + 1
endwhile

* Upper air variables
i = 1
levtype = '   '
while(i <= _nuvar )
      name   = subwrd(_uvars,i) '        '
      title  = _ulong.i '                                                                     '
      name   = substr(name,1,8)
      title  = substr(title,1,70) 
      idc = substr(id' ',1,2)
      record = name '| ' idc ' |' title ' | ' units ' | ' dattype ' | ' levtype ' | ' scale_factor ' | ' precision ' |   |   |' 
      rc = write(fname,record,append)
      i = i + 1
      id = id + 1
endwhile

*
* Vertical dimension
*
rc=write(fname,'#---------------------------------------------------------------------------------------------------',append)
rc=write(fname,'#!  vert',append)
rc=write(fname,'# Vertical dimension types',append)
rc=write(fname,'#   (max number of entries = LATS_MAX_VERT_TYPES in lats.h)',append)
rc=write(fname,'#',append)
rc=write(fname,'# The format of each record is:',append)
rc=write(fname,'#   level_type | description | units | verticality | positive | default | GRIB_id | GRIB_p1 | GRIB_p2 | GRIB_p3',append)
rc=write(fname,'#',append)
rc=write(fname,'# level_type = level type',append)
rc=write(fname,'# description = level description',append)
rc=write(fname,'# units = units for this level type',append)
rc=write(fname,'# verticality = 'single' (single surface) or 'multi' (variable can have multiple levels of this type)',append)
rc=write(fname,'# positive = 'up' (increasing values point up) or 'down' (increasing values point down)',append)
rc=write(fname,'# GRIB_id = GRIB level type indicator (PDS octet 10)',append)
rc=write(fname,'# GRIB_p1 = GRIB PDS octet 11',append)
rc=write(fname,'# GRIB_p2 = GRIB PDS octet 12',append)
rc=write(fname,'# GRIB_p3 = combined GRIB octets 11, 12 - overrides values of GRIB_p1, GRIB_p2 if specified',append)
rc=write(fname,'',append)
rc=write(fname,'0degiso	 | 0 deg isotherm    	     | hPa	| single |   up	|    4 | 0 |  0 | 0',append)
rc=write(fname,'atm	 | Atmosphere (entire)	     |          | single |   up |  200 | 0 |  0 | 0 ',append)
rc=write(fname,'ocn	 | Ocean (entire depth)	     |          | single |   up |  201 | 0 |  0 | 0 ',append)
rc=write(fname,'clhbot	 | High Cloud Bottom Level   | hPa      | single |   up	|  232 | 0 |  0 | 0',append)
rc=write(fname,'clhlay	 | High Cloud Top Layer      |          | single |   up	|  234 | 0 |  0 | 0',append)
rc=write(fname,'clhtop	 | High Cloud Top Level      | hPa      | single |   up	|  233 | 0 |  0 | 0',append)
rc=write(fname,'cllbot	 | Low Cloud Bottom Level    | hPa      | single |   up	|  212 | 0 |  0 | 0',append)
rc=write(fname,'clllay	 | Low Cloud Top Layer       |          | single |   up	|  214 | 0 |  0 | 0',append)
rc=write(fname,'clltop	 | Low Cloud Top Level       | hPa      | single |   up	|  213 | 0 |  0 | 0',append)
rc=write(fname,'clmbot	 | Mid Cloud Bottom Level    | hPa      | single |   up	|  222 | 0 |  0 | 0',append)
rc=write(fname,'clmlay	 | Mid Cloud Top Layer       |          | single |   up	|  224 | 0 |  0 | 0',append)
rc=write(fname,'clmtop	 | Mid Cloud Top Level       | hPa      | single |   up	|  223 | 0 |  0 | 0',append)
rc=write(fname,'cltbot	 | Cloud base level 	     | hPa	| single |   up	|    2 | 0 |  0 | 0',append)
rc=write(fname,'cltlay	 | Total Cloud layer 	     |		| single |   up	|    3 | 0 |  0 | 0',append)
rc=write(fname,'cltmax	 | Highest Cloud height      | m        | single |   up	|  105 | 0 |  0 | 0',append)
rc=write(fname,'landd	 | Below ground, 10 to 200 cm|		| single |   up |  112 |10 |200 | 0',append)
rc=write(fname,'lands	 | Below ground, 0 to 10 cm  |		| single |   up |  112 | 0 | 10 | 0',append)
rc=write(fname,'landt	 | Below ground, 0  to 200 cm|		| single |   up |  112 | 0 |200 | 0',append)
rc=write(fname,'lcl      | Adiabatic cond level      | hPa	| single |   up	|    5 | 0 |  0 | 0',append)
rc=write(fname,'maxwnd   | Maximum wind speed level  | hPa 	| single |   up	|    6 | 0 |  0 | 0',append)
rc=write(fname,'msl	 | Mean Sea Level 	     |		| single |   up	|  102 | 0 |  0 | 0',append)
rc=write(fname,'ocnbot	 | Ocean bottom      	     |		| single |   up	|    9 | 0 |  0 | 0',append)
rc=write(fname,'plev	 | Pressure level	     | hPa	| multi  | down |  100 | 0 |  0 | 0',append)
rc=write(fname,'pbltop	 | Top of PBL       	     |		| single |   up	|   21 | 0 |  0 | 0',append)
rc=write(fname,'sfc      | Earth surface             |          | single |   up |    1 | 0 |  0 | 0',append)
rc=write(fname,'sfclo    | Sfc Layer Ocean           |          | single |   up |  112 | 0 |300 | 0',append)
rc=write(fname,'sfc10m	 | 10 meters above earth surface| m	| single |   up	|  105 | 0 |  0 | 10',append)
rc=write(fname,'sfc2m	 | 2 meters above earth surface| m	| single |   up	|  105 | 0 |  0 | 2',append)
rc=write(fname,'toa	 | Top of atmosphere	     |		| single |   up	|    8 | 0 |  0 | 0',append)
rc=write(fname,'modtop	 | Top of Model     	     |		| single |   up	|   20 | 0 |  0 | 0',append)
rc=write(fname,'toasat   | TOA satellite             |     	| single |   up	|   22 | 0 |  0 | 0',append)
rc=write(fname,'troplev  | Tropopause level          | hPa 	| single |   up	|    7 | 0 |  0 | 0',append)
rc=write(fname,'theta    | Isentropic Level          | K        | multi  |   up	|  113 | 0 |  0 | 0',append)
rc=write(fname,'sigma	 | Sigma level               |          | multi  | down	|  107 | 0 |  0 | 0',append)
rc=write(fname,'hybrid   | Hybrid Model level number |          | multi  |   up	|  109 | 0 |  0 | 0',append)
rc=write(fname,'zocean   | Depth below sea level     | m        | multi  | down	|  160 | 0 |  0 | 0',append)
rc=write(fname,'',append)
rc=write(fname,'#---------------------------------------------------------------------------------------------------',append)
rc=write(fname,'#!	Center',append)
rc=write(fname,'# Modeling centers (GRIB only)',append)
rc=write(fname,'#   (max number of entries = LATS_MAX_CENTERS in lats.h)',append)
rc=write(fname,'#',append)
rc=write(fname,'# The format of each record is:',append)
rc=write(fname,'#   center | GRIB_id | GRIB_center | GRIB_subcenter',append)
rc=write(fname,'#',append)
rc=write(fname,'# center = mnemonic for the center',append)
rc=write(fname,'# GRIB_id = GRIB generating process id (PDS octet 6)',append)
rc=write(fname,'# GRIB_center = the id of center managing the data (for AMIP II this is PCMDI) - see GRIB Table 0',append)
rc=write(fname,'# GRIB_subcenter = the id of the subcenter',append)
rc=write(fname,'# ',append)
rc=write(fname,'#',append)
rc=write(fname,'#  Acronym           AMIP Group                                                    Location',append)
rc=write(fname,'#  -------           ----------                                                    --------',append)
rc=write(fname,'#',append)
rc=write(fname,'#  bmrc              Bureau of Meteorology Research Centre                         Melbourne, Australia',append)
rc=write(fname,'#  ccc               Canadian Centre for Climate Modelling and Analysis            Victoria, Canada',append)
rc=write(fname,'#  ccsr              Center for Climate System Research                            Tokyo, Japan',append)
rc=write(fname,'#  cnrm              Centre National de Recherches Meteorologiques                 Toulouse, France',append)
rc=write(fname,'#  cola              Center for Ocean-Land-Atmosphere Studies                      Calverton, Maryland',append)
rc=write(fname,'#  csiro             Commonwealth Scientific & Industrial Research Organization    Mordialloc, Australia',append)
rc=write(fname,'#  csu               Colorado State University                                     Fort Collins, Colorado',append)
rc=write(fname,'#  derf              Dynamical Extended Range Forecasting (at GFDL)                Princeton, New Jersey',append)
rc=write(fname,'#  dnm               Department of Numerical Mathematics                           Moscow, Russia',append)
rc=write(fname,'#  ecmwf             European Centre for Medium-Range Weather Forecasts            Reading, England',append)
rc=write(fname,'#  gfdl              Geophysical Fluid Dynamics Laboratory                         Princeton, New Jersey',append)
rc=write(fname,'#  giss              Goddard Institute for Space Studies                           New York, New York',append)
rc=write(fname,'#  gla               Goddard Laboratory for Atmospheres                            Greenbelt, Maryland',append)
rc=write(fname,'#  gsfc              Goddard Space Flight Center                                   Greenbelt, Maryland',append)
rc=write(fname,'#  iap               Institute of Atmospheric Physics                              Beijing, China',append)
rc=write(fname,'#  jma               Japan Meteorological Agency                                   Tokyo, Japan',append)
rc=write(fname,'#  llnl              Lawrence Livermore National Laboratory                        Livermore, California',append)
rc=write(fname,'#  lmd               Laboratoire de Meteorologie Dynamique                         Paris, France',append)
rc=write(fname,'#  mgo               Main Geophysical Observatory                                  St. Petersburg, Russia',append)
rc=write(fname,'#  mpi               Max-Planck-Institut fur Meteorologie                          Hamburg, Germany',append)
rc=write(fname,'#  mri               Meteorological Research Institute                             Ibaraki-ken, Japan',append)
rc=write(fname,'#  ncar              National Center for Atmospheric Research                      Boulder, Colorado',append)
rc=write(fname,'#  nmc               National Meteorological Center                                Suitland, Maryland',append)
rc=write(fname,'#  nrl               Naval Research Laboratory                                     Monterey, California',append)
rc=write(fname,'#  ntu               National Taiwan University                                    Taipei, Taiwan',append)
rc=write(fname,'#  pcmdi             Program for Climate Model Diagnosis and Intercomparison, LLNL Livermore, California',append)
rc=write(fname,'#  rpn               Recherche en Privision Numerique                              Dorval, Canada',append)
rc=write(fname,'#  sunya             State University of New York at Albany                        Albany, New York',append)
rc=write(fname,'#  sunya/ncar        State University of New York at Albany/NCAR                   Albany, New York/Boulder, Colorado',append)
rc=write(fname,'#  ucla              University of California at Los Angeles                       Los Angeles, California',append)
rc=write(fname,'#  ugamp             The UK Universities Global Atmospheric Modelling Programme   Reading, England',append)
rc=write(fname,'#  uiuc              University of Illinois at Urbana-Champaign                    Urbana, Illinois',append)
rc=write(fname,'#  ukmo              United Kingdom Meteorological Office                          Bracknell, UK',append)
rc=write(fname,'#  yonu              Yonsei University                                             Seoul, Korea',append)
rc=write(fname,'#',append)
rc=write(fname,'',append)
rc=write(fname,'bmrc	  |  1  | 100 | 2',append)
rc=write(fname,'ccc	  |  2  | 100 | 2',append)
rc=write(fname,'cnrm	  |  3  | 100 | 2',append)
rc=write(fname,'cola	  |  4  | 100 | 2',append)
rc=write(fname,'csiro	  |  5  | 100 | 2',append)
rc=write(fname,'csu	  |  6  | 100 | 2',append)
rc=write(fname,'dnm	  |  7  | 100 | 2',append)
rc=write(fname,'ecmwf	  |  8  | 100 | 2',append)
rc=write(fname,'gfdl	  |  9  | 100 | 2',append)
rc=write(fname,'derf      | 10  | 100 | 2',append)
rc=write(fname,'giss	  | 11  | 100 | 2',append)
rc=write(fname,'gla	  | 12  | 100 | 2',append)
rc=write(fname,'gsfc	  | 13  | 100 | 2',append)
rc=write(fname,'iap	  | 14  | 100 | 2',append)
rc=write(fname,'jma	  | 15  | 100 | 2',append)
rc=write(fname,'lmd	  | 16  | 100 | 2',append)
rc=write(fname,'mgo	  | 17  | 100 | 2',append)
rc=write(fname,'mpi	  | 18  | 100 | 2',append)
rc=write(fname,'mri	  | 19  | 100 | 2',append)
rc=write(fname,'ncar	  | 20  | 100 | 2',append)
rc=write(fname,'ncep	  | 21  | 100 | 2',append)
rc=write(fname,'nrl	  | 22  | 100 | 2',append)
rc=write(fname,'rpn	  | 23  | 100 | 2',append)
rc=write(fname,'sunya	  | 24  | 100 | 2',append)
rc=write(fname,'sunya/ncar| 25  | 100 | 2',append)
rc=write(fname,'ucla	  | 26  | 100 | 2',append)
rc=write(fname,'ugamp	  | 27  | 100 | 2',append)
rc=write(fname,'uiuc	  | 28  | 100 | 2',append)
rc=write(fname,'ukmo	  | 29  | 100 | 2',append)
rc=write(fname,'yonu	  | 30  | 100 | 2',append)
rc=write(fname,'ccsr      | 31  | 100 | 2',append)
rc=write(fname,'llnl      | 32  | 100 | 2',append)
rc=write(fname,'ntu       | 33  | 100 | 2',append)
rc=write(fname,'cptec     | 46  | 100 | 2',append)
rc=write(fname,'pcmdi	  | 100 | 100 | 2',append)
rc=write(fname,'#---------------------------------------------------------------------------------------------------',append)
rc=write(fname,'#!qc',append)
rc=write(fname,'# Quality control marks',append)
rc=write(fname,'#   (no limit on number of entries)',append)
rc=write(fname,'#',append)
rc=write(fname,'# The format of each record is:',append)
rc=write(fname,'#   variable | level_type | level | mean | std | tolerance | range | rangetol',append)
rc=write(fname,'#',append)
rc=write(fname,'# variable = variable name',append)
rc=write(fname,'# level_type = type of level, as defined in the leveltypes section, or blank if no associated level',append)
rc=write(fname,'# level = level value, or blank if no specified level',append)
rc=write(fname,'# mean = observed mean at specified level',append)
rc=write(fname,'# std = observed standard deviation at specified level',append)
rc=write(fname,'# tolerance = number of standard deviations about mean',append)
rc=write(fname,'#     (if abs(calculated_mean - mean) > (std * tolerance), flag the value as 'mean out of range')',append)
rc=write(fname,'# range = observed (maximum - minimum)',append)
rc=write(fname,'# rangetol = range tolerance:',append)
rc=write(fname,'#     (if calculated_range > (rangetol * range), flag as 'range is too large')',append)
rc=write(fname,'',append)
rc=write(fname,'# NOTE: no QC table yet',append)
rc=write(fname,'#',append)
rc=subwrd(rc,1)
if(rc!=0)
   say _myname 'problems creating LATS table file ' fname 
   return rc
endif

rc=close(fname)
rc=subwrd(rc,1)
if(rc>0)
   say _myname 'problems closing LATS table file ' fname 
   return rc
endif

return 0


*
* mkgeosf()  Creates a template GEOS ctl for regridding;
*            also sets _dimenv, _func
*

function mkgeosf(fname,dlon,dlat,nlon,nlat,algo)

*
* Temporary file
*

if ( _verb = 1 )
      say _myname 'creating GEOS template file ' fname
endif

rc=write(fname,'DSET nofile')
rc=subwrd(rc,1)
if ( rc > 0 )
   say rc
   say _myname 'cannot create GEOS template file ' fname ' on current directory'
   return rc
endif

xdef = 'xdef ' nlon ' linear -180 ' dlon 
ydef = 'ydef ' nlat ' linear -90  ' dlat
rc=write(fname,'title Template GrADS for regridding',append)
rc=write(fname,'options template',append)
rc=write(fname,'undef 1e+20',append)
rc=write(fname,xdef,append)
rc=write(fname,ydef,append)
rc=write(fname,'zdef   1 levels  1000 ',append)
rc=write(fname,'tdef 1 linear 0Z1jan1900 1dy',append)
rc=write(fname,'vars 1',append)
rc=write(fname,'var       0 0        generic sfc variable',append)
rc=write(fname,'endvars',append)

rc=close(fname)
rc=subwrd(rc,1)
if(rc>0)
   say _myname 'problems closing GEOS template file ' fname 
   return rc
endif

_dimenv = fname
_func = 'regrid2(@,' dlon ',' dlat ',' algo ',-180,-90)'


return 0


*
* validate(var,vars)  See whether token var is in the list vars
*
function validate(var,vars)
      if ( vars = '' ); return 1; endif
      valid = 0
      i = 1
      v = subwrd(vars,i)
      while ( v != '' )
         if ( var = v ); valid = 1; endif
         i = i + 1
         v = subwrd(vars,i)
      endwhile
return valid

*
* xvalidate(var,vars)  See whether token var is NOT in the list vars
*
function xvalidat(var,vars)
      if ( vars = '' ); return 1; endif
      valid = 1
      i = 1
      v = subwrd(vars,i)
      while ( v != '' )
         if ( var = v ); valid = 0; endif
         i = i + 1
         v = subwrd(vars,i)
      endwhile
return valid

*
* defvars()  Define variables in output LATS file
*

function defvars()

* BUG: for now, cannot handle "accum" variables, oh well...
* ---------------------------------------------------------
  if ( _mean = 1 )
       timestat = 'average'
  else
       timestat = 'instant'
  endif

* Surface variables
* ------------------
  n = 1; _svids = ''
  while ( n <= _nsvar )
     var = subwrd(_svars,n)
     if ( _format != stream )
       'set lats var  ' _fid ' ' var ' ' timestat ' ' _gid ' '  _sid
     endif
     vid = subwrd(result,5)
     if ( vid = 0 ); return 1; endif
     _svids = _svids ' ' vid
     n = n + 1
  endwhile


* Upper air variables
* -------------------
  n = 1; _uvids = ''
  while ( n <= _nuvar )
     var = subwrd(_uvars,n)
     if ( _format != stream )
        'set lats var  ' _fid ' ' var ' ' timestat ' ' _gid ' '  _lid
     endif
     vid = subwrd(result,5)
     _uvids = _uvids ' ' vid
     n = n + 1
  endwhile

return 0

*............................................................................

* subst(var,expr) Substitute all the occurences of '@' in expr with the 
*                 contents of var. Example:
*
*                 expr = cos(lat*pi/180)*@*hgrad(@)
*                 var  = ua
*
*                 results in cos(lat*pi/180)*ua*hgrad(ua).
*
*                 If expr='' it returns var.
*
function subst(var,expr)

      if ( expr ='' ); return var; endif

      str = ''
      i = 1
      ch = substr(expr,i,1)
      while( ch != '' )

           if ( ch = '@' )
                str = str '' var
           else
                str = str '' ch
           endif

           i = i + 1
           ch = substr(expr,i,1)

      endwhile

return str


*............................................................................

*
* chkcfg()  Check whether grads version is powerful enough to run lats4d
*

function chkcfg()

      'q config'
      if ( rc != 0 )
         say _myname 'GrADS version appears earlier than 1.7Beta'
         return 1
      endif

*     Check for LATS
*     --------------
      lats = 0
      cfg = sublin(result,1)
      i = 1
      version = subwrd(cfg,2)
      token = subwrd(cfg,i)
      while ( token != '' )
         if ( token = 'lats' ); lats = 1; endif
         i = i + 1
         token = subwrd(cfg,i)
      endwhile

      if ( lats = 0 )
         say _myname 'this build of GrADS ' version ' does not include the LATS option'
         say _myname 'please try "gradsnc" or "gradshdf", version 1.7beta9 or later'
         return 1
      endif

return 0

*............................................................................

*
* crtique(nmin,nmax,opt,string)  Returns 0 if the number of words in "string" 
*                            is in between nmin and nmax, returns 1 otherwise.
*

function critique ( nmin, nmax, opt, string)

      n = 1
      while ( subwrd(string,n) != '' ); n = n + 1; endwhile
      n = n - 1

      if ( n < nmin )
          rc = usage()
          say _myname 'you specified "' opt ' ' string '"'
          say _myname 'option "' opt '" requires at least ' nmin ' argument(s)'
          return 1
      endif

      if ( n > nmax )
          rc = usage()
          say _myname 'you specified "' opt ' ' string '"'
          say _myname 'option "' opt '" requires at most ' nmax ' argument(s)'
          return 1
      endif
   

return 0



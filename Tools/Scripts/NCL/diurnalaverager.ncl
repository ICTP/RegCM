;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;
;    This file is part of ICTP RegCM.
;
;    ICTP RegCM is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    ICTP RegCM is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
;
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin

;*******************************************************************************
;************************ Input Arguments **************************************
;*******************************************************************************
;    infilename = "./WNA_UW_OPT.2008060100.nc"
;    outfilename = "tmp.nc"
;    starthour = 0

;*******************************************************************************
;********************* Initalization *******************************************
;*******************************************************************************

  NDIMMAX = 6

  nmaxvars = 0
  nminvars = 0
  nstdvars = 0


  

;*******************************************************************************
;********************* Open the Input File and Read ****************************
;*******************************************************************************
  fin = addfile(infilename,"r")

  ;Read in the first two steps of the time variable
  time01 = fin->time(0:1)
  time01_parsed = ut_calendar(time01,0)
  ;Determine the time spacing
  dt = time01_parsed(1,3) - time01_parsed(0,3)
  ;Determine the starting hour
  starthour = time01_parsed(0,3)
  if(dt.lt.1.or.dt.gt.23)then
    print("Warning: " + infilename + " does not appear to be an hourly output file")
    exit
  end if
  ;Set the output data index offset; this allows the 
  ;starting hour to be different from 0
  kr0 = round(starthour/dt,3)
  ;Check if the start hour is divisible by dt
  if(mod(starthour,dt).ne.0)then
    print("Error: the variable starthour does not appear to be consistent with the output timestep.  Such a file does not work with this analysis algorithm.")
  end if
          
  ;Determine the number of steps per day
  stepsperday = round(24/dt,3)
  if(mod(24,dt).ne.0)then
    ;Give an error message
    print("Error: dt does not divide in to 24 nicely")
  end if

    

;Get all the variable and dimension names
;Go through each variable; get its dimensionality  and determine if it is a
;min/max variable

  varnames = getfilevarnames(fin)
  nvartot = dimsizes(varnames)

  dimlengths = getfiledimsizes(fin)
  numdims = dimsizes(dimlengths)
  

  maxvarindstmp = new(nvartot,integer)
  minvarindstmp = new(nvartot,integer)
  stdvarindstmp = new(nvartot,integer)

  dimnames = new((/nvartot,numdims/),string)
  varranks = new(nvartot,integer)
  istimevar = new(nvartot,logical)
  timedimind = new(nvartot,integer)
  ntime = new(nvartot,integer)

  istimevar = False

  do i = 0, nvartot-1
    isminmax = False
    if(.not.ismissing(str_index_of_substr(varnames(i),"max",-1)))then
      isminmax = True
      maxvarindstmp(nmaxvars) = i
      nmaxvars = nmaxvars + 1
    end if
    if(.not.ismissing(str_index_of_substr(varnames(i),"min",-1)))then
      isminmax = True
      minvarindstmp(nminvars) = i
      nminvars = nminvars + 1
    end if
    if(.not.isminmax)then
      stdvarindstmp(nstdvars) = i
      nstdvars = nstdvars + 1
    end if

    vardimlengths = getfilevardimsizes(fin,varnames(i))
    varranks(i) = dimsizes(vardimlengths)
    do j = 0,varranks(i)-1
      dimnames(i,j) = fin->$varnames(i)$!j
      if(dimnames(i,j).eq."time".and.varranks(i).gt.1)then
        istimevar(i) = True
        timedimind(i) = j
        ntime(i) = vardimlengths(j)
      end if
    end do

    delete(vardimlengths)
  end do


  if(nmaxvars.ne.0)then
    maxvarinds = new(nmaxvars,integer)
    maxvarinds = maxvarindstmp(0:nmaxvars-1)
  end if
  if(nminvars.ne.0)then
    minvarinds = new(nminvars,integer)
    minvarinds = minvarindstmp(0:nminvars-1)
  end if
  stdvarinds = new(nstdvars,integer)
  stdvarinds = stdvarindstmp(0:nstdvars-1)

  delete(maxvarindstmp)
  delete(minvarindstmp)
  delete(stdvarindstmp)

;*******************************************************************************
;********************* Open the Output File ************************************
;*******************************************************************************
  system("rm " + outfilename)
  fout = addfile(outfilename,"c")

;*******************************************************************************
;********************* Do the Averaging ****************************************
;*******************************************************************************

;Go through each normal variable and average over each hour-increment and copy
;all attributes

  do m = 0,nstdvars-1
    varind = stdvarinds(m)
    vardata = fin->$varnames(varind)$


    ;If it is a time variable, then do averaging
    if(istimevar(varind))then
      ;Set the dimension lengths
      dsize = getfilevardimsizes(fin,varnames(varind))
      ;Overwrite the time dimension
      dsize(timedimind(varind)) = stepsperday
      
      ;Declare a variable of the proper dimensionality
      dumavg = new(dsize,typeof(vardata))

      ;Copy the dimensions
      do j = 0, varranks(varind)-1
        dname = dimnames(varind,j)
        dumavg!j = dname
        if(dname.ne."time")then
          if(isfilevarcoord(fin,varnames(varind),dname))then
            dumavg&$dname$ = vardata&$dname$
          end if
        else
          dumavg&time = vardata&time(0:(stepsperday-1))
          ;Set the time variable in such a way that the output time starts
          ;at 00Z
          if(kr0.ne.0)then
            dumavg&time(kr0:(stepsperday-1)) = vardata&time(0:(stepsperday-1-kr0))
            dumavg&time(0:(kr0-1)) = vardata&time((stepsperday-kr0):(stepsperday-1))
          end if
        end if
      end do
      ;Rename the time dimension
      dumavg!0 = "hour"
      ;Copy Attributes
      copy_VarAtts(vardata,dumavg)
  

      print("Averaging '" + varnames(varind) + "'")
      ;Average over the subset of the data that corresponds to the current 
      ;timestep
      do k = 0,stepsperday -1

        kvar = k+kr0
        if(kvar.gt.(stepsperday-1))then
          kvar = (k+kr0) - stepsperday
        end if

        if(varranks(varind).eq.2)then
          vardatasub = vardata(k:(ntime(varind)-1):stepsperday,:)
          dumavg(kvar,:) = dim_avg_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.3)then
          vardatasub = vardata(k:(ntime(varind)-1):stepsperday,:,:)
          dumavg(kvar,:,:) = dim_avg_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.4)then
          vardatasub = vardata(k:(ntime(varind)-1):stepsperday,:,:,:)
          dumavg(kvar,:,:,:) = dim_avg_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.5)then
          vardatasub = vardata(k:(ntime(varind)-1):stepsperday,:,:,:,:)
          dumavg(kvar,:,:,:,:) = dim_avg_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.6)then
          vardatasub = vardata(k:(ntime(varind)-1):stepsperday,:,:,:,:,:)
          dumavg(kvar,:,:,:,:,:) = dim_avg_n(vardatasub,timedimind(varind))
        end if
        delete(vardatasub)
      end do


      ;Write the data
      fout->$varnames(varind)$ = dumavg

      ;Cleanup
      delete(dsize)
      delete(dumavg)
    ;Otherwise just copy the variable
    else
      print("Copying '" + varnames(varind) + "'")
      fout->$varnames(varind)$ = vardata
    end if

    ;Cleanup
    delete(vardata)
  end do

    
;Go through each min/max variable and get the maximum/minimum value for each day

;Max vars
if(nmaxvars.ne.0)then
  do m = 0,nmaxvars-1
    varind = maxvarinds(m)
    vardata = fin->$varnames(varind)$


    ;If it is a time variable, then do averaging
    if(istimevar(varind))then
      ndays = (ntime(varind)+kr0)/stepsperday
      ;Set the dimension lengths
      dsize = getfilevardimsizes(fin,varnames(varind))
      ;Overwrite the time dimension
      dsize(timedimind(varind)) = ndays
      
      ;Declare a variable of the proper dimensionality
      dummax = new(dsize,typeof(vardata))

      ;Copy the dimensions
      do j = 0, varranks(varind)-1
        dname = dimnames(varind,j)
        dummax!j = dname
        if(dname.ne."time")then
          if(isfilevarcoord(fin,varnames(varind),dname))then
            dummax&$dname$ = vardata&$dname$
          end if
        else
          dummax&time = vardata&time(0:(ntime(varind)-1):stepsperday)
        end if
      end do
      ;Rename the time dimension
      dummax!0 = "hour"
      ;Copy Attributes
      copy_VarAtts(vardata,dummax)
  

      print("Averaging '" + varnames(varind) + "'")
      ;Average over the subset of the data that corresponds to the current 
      ;timestep
      do k = 0,ndays-1
        if(varranks(varind).eq.2)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:)
          dummax(k,:) = dim_max_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.3)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:,:)
          dummax(k,:,:) = dim_max_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.4)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:,:,:)
          dummax(k,:,:,:) = dim_max_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.5)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:,:,:,:)
          dummax(k,:,:,:,:) = dim_max_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.6)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:,:,:,:,:)
          dummax(k,:,:,:,:,:) = dim_max_n(vardatasub,timedimind(varind))
        end if
        delete(vardatasub)
      end do

      dummaxavg = dim_avg_n_Wrap(dummax,0)

      dummaxavg@method = "Average of daily maximum"

      ;Write the data
      fout->$varnames(varind)$ = dummaxavg

      ;Cleanup
      delete(dsize)
      delete(dummax)
      delete(dummaxavg)
    ;Otherwise just copy the variable
    else
      print("Copying '" + varnames(varind) + "'")
      fout->$varnames(varind)$ = vardata
    end if

    ;Cleanup
    delete(vardata)
  end do
end if

;Min vars
if(nminvars.ne.0)then
  do m = 0,nminvars-1
    varind = minvarinds(m)
    vardata = fin->$varnames(varind)$


    ;If it is a time variable, then do averaging
    if(istimevar(varind))then
      ndays = (ntime(varind) + kr0)/stepsperday
      ;Set the dimension lengths
      dsize = getfilevardimsizes(fin,varnames(varind))
      ;Overwrite the time dimension
      dsize(timedimind(varind)) = ndays
      
      ;Declare a variable of the proper dimensionality
      dummin = new(dsize,typeof(vardata))

      ;Copy the dimensions
      do j = 0, varranks(varind)-1
        dname = dimnames(varind,j)
        dummin!j = dname
        if(dname.ne."time")then
          if(isfilevarcoord(fin,varnames(varind),dname))then
            dummin&$dname$ = vardata&$dname$
          end if
        else
          dummin&time = vardata&time(0:(ntime(varind)-1):stepsperday)
        end if
      end do
      ;Rename the time dimension
      dummin!0 = "hour"
      ;Copy Attributes
      copy_VarAtts(vardata,dummin)
  

      print("Averaging '" + varnames(varind) + "'")
      ;Average over the subset of the data that corresponds to the current 
      ;timestep
      do k = 0,ndays-1
        if(varranks(varind).eq.2)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:)
          dummin(k,:) = dim_min_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.3)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:,:)
          dummin(k,:,:) = dim_min_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.4)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:,:,:)
          dummin(k,:,:,:) = dim_min_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.5)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:,:,:,:)
          dummin(k,:,:,:,:) = dim_min_n(vardatasub,timedimind(varind))
        end if
        if(varranks(varind).eq.6)then
          vardatasub = vardata((k*stepsperday):((k+1)*stepsperday-1-kr0),:,:,:,:,:)
          dummin(k,:,:,:,:,:) = dim_min_n(vardatasub,timedimind(varind))
        end if
        delete(vardatasub)
      end do

      dumminavg = dim_avg_n_Wrap(dummin,0)

      dumminavg@method = "Average of daily minimum"

      ;Write the data
      fout->$varnames(varind)$ = dumminavg

      ;Cleanup
      delete(dsize)
      delete(dummin)
      delete(dumminavg)
    ;Otherwise just copy the variable
    else
      print("Copying '" + varnames(varind) + "'")
      fout->$varnames(varind)$ = vardata
    end if

    ;Cleanup
    delete(vardata)
  end do
end if

end

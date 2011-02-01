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

;This script takes yearly averages of monthly average data.  It is made to be
;general for any type of data, with a maximum rank of 5
;
;Created 01/2011 by Travis A. O'Brien for RegCM v4.1

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
  ;Input arguments
;  startyear = 1984
;  startmonth = 12
;  fileinname = ftype+".all.nc"
;  fileoutname = ftype+".yearly.avg.nc"

  frcm = addfile(fileinname,"r")

  system("rm "+fileoutname) ;Remove any prexisting version of the file
  ;Add the output file
  fileout = addfile(fileoutname,"c")


  ;Copy file attributes
  ;Get the list of variables from the file
  attlist = getvaratts(frcm)
  ;Get the number of variables
  natts = dimsizes(attlist) 
  ;Go through the list of attributes and copy them directly from the input to
  ; the output file
  if(natts.gt.0)then
      do iatt = 0,natts - 1
        ;Check that the attribute itself isn't missing (this will happen if
        ; there are no global attributes in the file).
        if(.not.ismissing(attlist(iatt)))then
          fileout@$attlist(iatt)$ = frcm@$attlist(iatt)$
        end if
      end do
  end if

  ;Get the list of variables from the file
  varlist = getfilevarnames(frcm)
  ;Get the number of variables
  nvars = dimsizes(varlist) 

  ;Loop through each variable, and take the average of each file
  ; this will create a file with identical variables and dimensions as the
  ; original, except the time dimension will be about 12 times smaller
  do ivar = 0,nvars-1
    ;Pre-delete the variable 'myvar' if necessary
    if(isdefined("myvar"))then
      delete(myvar)
    end if
    ;Get the variable's data from the file
    myvar = frcm->$varlist(ivar)$

    ;Only attempt to process if we are dealing with a 'time' variable
    if(isfilevardim(frcm,varlist(ivar),"time"))then
      ;Pre-delete several more variables
      if(isdefined("dsize"))then
        delete(dsize)
      end if
      if(isdefined("dimnames"))then
        delete(dimnames)
      end if
      ;Get the rank, dimension sizes, and dimension names of the variable
      dsize = dimsizes(myvar) ;Dimension sizes
      drank = dimsizes(dsize) ;Variable rank
      dimnames = getvardims(myvar) ;Dimension names

      ;Only process variables whose first dimension is five and that have rank
      ; 2,3,4 or 5.  This is because there is code in here that cannot be 
      ; generalized for an arbitrary-rank variable.
      if(dimnames(0).eq."time".and.drank.le.5.and.drank.gt.1)then
        ;Pre-delete several more variables
        if(isdefined("varyear"))then
          delete(varyear)
        end if
        if(isdefined("dsizeyear"))then
          delete(dsizeyear)
        end if
        ;Get the number of time steps
        ntime = dsize(0)
        ;Make a new set of dimension sizes, where we'll overwrite the size of
        ; the time dimension with the number of years to average
        dsizeyear = dsize
        ;Calculate the number of years that we'll get from the data
        ;nyears = round(floor(ntime/12),3) +1
        nyears = round(floor(ntime/12),3)

        ;Set the time dimension to be the number of years
        dsizeyear(0) = nyears

        ;Define a variable that is dimensioned identically to myvar, except that
        ; the time dimension is the number of years
        ;varyear = new(dsizeyear,float)
        varyear = new(dsizeyear,typeof(myvar))

        ;Define the time dimension, since it won't be copied like the other
        ; dimensions
        varyear!0 = "time"
        varyear&time = ispan(0,nyears-1,1)
        varyear&time@units = "years since "+startyear+"-"+startmonth+"-01"

        ;Copy the variable's attributes
        copy_VarAtts(myvar,varyear)

        ;Go through each dimension and copy the dimension's name and coordinate
        ; variable (if applicable)
        do irank = 1,drank-1
          ;Copy the dim name
          varyear!irank = myvar!irank
          ;Check if the coordinate variable is defined
          if(isfilevarcoord(frcm,varlist(ivar),dimnames(irank)))then
            ;copy it if so
            varyear&$dimnames(irank)$ = myvar&$dimnames(irank)$ 
          end if
        end do

        ;Take the average of each year
        do i = 0,nyears-1
          ;Get the index of the first month of the current year
          startind = 12*i
          ;Get the index of the last month of the current year
          endind = startind+11


          ;Do the averaging as long as we haven't gone beyond the end of the
          ; data set
          if(startind.ne.ntime)then
            ;If the end index is beyond the end of the file then truncate the
            ; last year
            if(endind.gt.(ntime-1))then
              endind = ntime-1
            end if
            ;Depending on the rank of the variable, take the time average between
            ; the start and end months of this year
            if(drank.eq.2)then
              varyear(i,:) = dim_avg(myvar( \
                                          $dimnames(1)$|:,  \
                                          time|startind:endind)) 
            end if
            if(drank.eq.3)then
              varyear(i,:,:) = dim_avg(myvar( \
                                          $dimnames(1)$|:,  \
                                          $dimnames(2)$|:,  \
                                          time|startind:endind)) 
            end if
            if(drank.eq.4)then
              varyear(i,:,:,:) = dim_avg(myvar( \
                                          $dimnames(1)$|:,  \
                                          $dimnames(2)$|:,  \
                                          $dimnames(3)$|:,  \
                                          time|startind:endind)) 
            end if
            if(drank.eq.5)then
              varyear(i,:,:,:,:) = dim_avg(myvar(\
                                          $dimnames(1)$|:,  \
                                          $dimnames(2)$|:,  \
                                          $dimnames(3)$|:,  \
                                          $dimnames(4)$|:,  \
                                          time|startind:endind)) 
            end if
          end if
        end do 

        ;Output the variable
        fileout->$varlist(ivar)$ = varyear
      else
        ;As long as we aren't dealing with the time variable, output
        ; a warning to the screen if we come across a non-standard variable
        if(varlist(ivar).ne."time")then
          print("(yearly_average.ncl): Warning, variable " + varlist(ivar) \
                + " was skipped because its rank is greater than 5, or because" \
                + " its first dimension was not time.")
        end if
      end if
    else
      ;If we aren't dealing with a time-dimensioned variable, then just copy
      ;it out to the output file
      fileout->$varlist(ivar)$ = myvar
    end if

  end do


end

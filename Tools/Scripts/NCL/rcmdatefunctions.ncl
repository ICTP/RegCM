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

;This file contains a set of date-related fucntions to assist gfdl_icbc.ncl
;
;Created 04/2011 by Travis A. O'Brien for RegCM v4.1

;Parses an integer date of YYYYMMDDHH format
undef("parse_date")
function parse_date(idate)
local  year, month, day, hour, yyyymmddhh, mmddhh, ddhh, hh
begin
  
  yyyymmddhh = int2dble(idate)

  year = floor(yyyymmddhh/1e6) 
  mmddhh = yyyymmddhh - 1e6*year
  month = floor(mmddhh/1e4)
  ddhh = mmddhh - 1e4*month
  day = floor(ddhh/1e2)
  hh = ddhh - 1e2*day
  hour = floor(hh)

  yyyymmddhh@year = year
  yyyymmddhh@month = month
  yyyymmddhh@day = day
  yyyymmddhh@hour = hour

  return(yyyymmddhh)
end

undef("increment_date")
function increment_date(dbledate,bcdt,bLeap)
local newdate,imonth
begin
  daysinmonths = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  ;TODO: deal with leap years
;  if(isleap(dbledate@year))then
;    daysinmonths(1) = 29 
;  else
;    daysinmonths(1) = 28 
;  end if

  newdate = dbledate

  newdate@hour = newdate@hour + bcdt

  if(newdate@hour.ge.24.0)then
    ;Wrap the hour
    newdate@hour = newdate@hour - 24.0
    ;Increment the day
    newdate@day = newdate@day + 1
    imonth = round(newdate@month,3)

    if(newdate@day.gt.daysinmonths(imonth-1))
      ;Wrap the day
      newdate@day = 1
      ;Increment the month
      newdate@month = newdate@month + 1
      
      if(newdate@month.gt.12)then
        ;Wrap the month
        newdate@month = 1
        ;Increment the year
        newdate@year = newdate@year + 1
      end if
    end if
  end if

  return(newdate)
   
end

undef("equiv_date")
function equiv_date(dbledate1,dbledate2)
local datesareequivalent
begin
  if( \
            (dbledate1@year.eq.dbledate2@year)  \
      .and. (dbledate1@month.eq.dbledate2@month)  \
      .and. (dbledate1@day.eq.dbledate2@day)  \
      .and. (dbledate1@hour.eq.dbledate2@hour)  \
    )then
      datesareequivalent = True
  else 
      datesareequivalent = False
  end if

  return(datesareequivalent)
end

undef("date2str")
function date2str(dbledate)
local datestr,iyear,imonth,iday,ihour
begin
  iyear = round(dbledate@year,3)
  imonth = round(dbledate@month,3)
  iday = round(dbledate@day,3)
  ihour = round(dbledate@hour,3)

  syear = sprinti("%0.4i",iyear)
  smonth = sprinti("%0.2i",imonth)
  sday = sprinti("%0.2i",iday)
  shour = sprinti("%0.2i",ihour)
  datestr = syear + "-" + smonth + "-" + sday + " " + shour + ":00:00 UTC"

  return(datestr)
end

undef("idatestr")
function idatestr(dbledate)
local datestr,iyear,imonth,iday,ihour
begin
  iyear = round(dbledate@year,3)
  imonth = round(dbledate@month,3)
  iday = round(dbledate@day,3)
  ihour = round(dbledate@hour,3)

  syear = sprinti("%0.4i",iyear)
  smonth = sprinti("%0.2i",imonth)
  sday = sprinti("%0.2i",iday)
  shour = sprinti("%0.2i",ihour)
  datestr = syear + smonth + sday + shour

  return(datestr)
end

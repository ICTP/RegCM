; 
;    Copyright 2010 Travis O'Brien, University of California 
; 
;    This file is part of ucsc_ncl 
; 
;    ucsc_ncl is free software: you can redistribute it under the  
;    terms of the GNU General Public License as published by 
;    the Free Software Foundation, either version 3 of the License, or 
;    (at your option) any later version. 
; 
;    This program is distributed in the hope that it will be useful, 
;    but WITHOUT ANY WARRANTY; without even the implied warranty of 
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
;    GNU General Public License for more details. 
; 
;    You should have received a copy of the GNU General Public License 
;    along with this program.  If not, see <http://www.gnu.org/licenses/>. 
;
;
;Takes the cloud fraction and geopotential height fields from raw RCM output,
;and calculates grid cells where fog occurs.  The dimensionality of the
;variables does not matter, as long as the level dimension is the outer-most
;dimesnion.  The returned fog variable will have one less dimension than the
;input data (no level).  Fog is defined following Johnstone and Dawson (2010),
;where FC >= 5/8 and HGT < 400 m.  This routine calculates the cloud fraction
;below 400 m (using the product method) and sets the fog variable to 1 if the
;calculated fraction is greater than 5/8.  TODO: consider also passing LWC for a
;threshold to avoid counting high-coverage, optically invisible clouds.
undef("calculateFog")
function calculateFog(cld:numeric, geohgt:numeric,res)
local fog,cfthresh,dsizecld,dsizehgt,rankcld,rankhgt,lt400,invcf400,totcf,dimname
begin

    cfthresh = 5.0/8.0

    bUseLWC = False
    fLWCthreshold = 1.e-5 ;Only count clouds that would have a LWP of about
                          ; 1 g/m2 or more in the lowest layer, if LWC
                          ; is being used as the cloud fraction variable

    ;****Check if options were set in the res variables (not used currently)****
    if(res.eq.True)then
        ;Check for attribute UseLWC, which indicates whether the cld field is 
        ;LWC (instead of cloud fraction), meaning that we should use a threshold
        ;to determine cloud fraction.
        if(isatt(res,"UseLWC")) then
          if(islogical(res@UseLWC)) then 
            bUseLWC = res@UseLWC
          else
            print("Warning (calculateFog): resource UseLWC is not a logical. Ignoring.")
            bUseLWC = False
          end if
        end if

        ;Check for a user-defined threshold for LWC
        if(isatt(res,"LWCthreshold"))then
          if(isnumeric(res@LWCthreshold)) then 
            fLWCthreshold = res@LWCthreshold
          else
            print("Warning (calculateFog): resource LWCthreshold is non-numeric. Using a default value instead.")
          end if
        end if
    end if

    
    ;****Check that cld and geohgt have the same dimensionality****
    dsizecld = dimsizes(cld)
    dsizehgt = dimsizes(geohgt)
    rankcld = dimsizes(dsizecld)
    rankhgt = dimsizes(dsizehgt)
    if(rankcld.ne.rankhgt)then
      print("Error (calculateFog): arguments 1 and 2 do not have the same rank.  They should.  Returning -1.")
      return -1
    end if
    do i = 0,rankcld-1
      if(dsizecld(i).ne.dsizehgt(i))then
        print("Error (calculateFog): Dimension " + i + " does not have the same size for arguments 1 and 2.  They should have identical dimensionality.  Returning -1")
        return -1
      end if
    end do

    ;****Find all gridpoints that are less than 400m; set them to one****
    lt400 = where(geohgt.le.400,1.0,0.0)
    
    ;****Make a new cld, where only cells that are less than 400m in height
    ;have a non-zero (or missing?) value****
    if(bUseLWC)then
      ;Do a simple sanity check to make sure that the user passed in LWC
      ; instead of cloud fraction
      maxcld = max(cld)
      if(maxcld.ge.0.1)then
        print("Warning (calculateFog): the user-defined options indicate that the first argument is LWC, but the maximum value is greater than 0.1, which suggests that this is not so.")
      end if
      ;Make cloud fraction be 1 where LWC exceeds the threshold
      ; and 0 elsewhere.  Assumes that when there is a cloud, there is
      ; full grid-cell coverage
      cldtmp = where(cld.gt.fLWCthreshold,1.0,0.0)
    else
      cldtmp = cld
    end if
    invcf400 = 1.0 - cldtmp*lt400 ; Subrtact from one for the next step
    delete(lt400)

    ;****Calculate the total cloud fraction from this field****
    totcf = 1.0 - dim_product(invcf400)
    delete(invcf400)

    ;****Use the where function to determine all lat/lons where the total
    ;cloud fraction is greater than 5/8; set them to 1****
    fog = where(totcf.ge.cfthresh,1.0,0.0)
    delete(totcf)

    ;****Set the proper dimension names/attributes for the new fog variable****
    do i = 0,rankcld-2
      dimname = cld!i
      fog!i = dimname
;      if(dimname.eq."time")then
      fog&$dimname$ = cld&$dimname$
;      end if
    end do
    fog@long_name = "Fog occurrence"
    fog@units = "unitless"
    if(bUseLWC)then
      fog@method = "LWC Threshold"
      fog@lwcthresh = fLWCthreshold
    else
      fog@method = "Cloud Fraction Threshold"
    end if
  
  ;****Return the fog variable****
  return fog
end

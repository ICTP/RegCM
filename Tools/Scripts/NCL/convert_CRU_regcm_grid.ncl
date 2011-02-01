;**************************************************
;************ Load necessary ncl code *************
;**************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;User-defined vars
	;Open the regcm header, where we'll grab lat2d and lon2d
	fho	=	addfile(domfilename,"r")
	;Open the CRU Temperature data file
	f1	=	addfile(crutempfilename,"r")
	;Open the CRU Precip data file
	f2	=	addfile(cruprefilename,"r")

	;Grab lat2d and lon2d
	lat2df		=	fho->xlat(:,:)
	lon2df		=	fho->xlon(:,:)

	;Get the number of latitudes and longitudes
	dvar		=	dimsizes(lat2df)
	nlatf		=	dvar(0)
	nlonf		=	dvar(1)

        ;remove the border points, since the RCM output files
        ;do the same
        lat2d = lat2df(1:nlatf-2,1:nlonf-2)
        lon2d = lon2df(1:nlatf-2,1:nlonf-2)

	nlat		=	nlatf-2
	nlon		=	nlonf-2

	;read the CRU lat, lon, and time
	lat_cru 	= f1->lat
	lon_cru		= f1->lon
	time_cru	= f1->time

        ;Find the start and end times
        uttime = ut_calendar(time_cru,-1)
        istart = ind(uttime.eq.istartdate)
        iend = ind(uttime.eq.ienddate)
        if(ismissing(istart).or.ismissing(iend))then
          print("Error: couldn't find corresponding start and end dates in CRU file.  Abort.")
          return
        end if

        ntime = dimsizes(time_cru(istart:iend))
        

        print("Reading t2m (TMP) from " + time_cru(istart) + " to " + time_cru(iend))
	;Read in the CRU data and give it proper units
	v1_ta		=	short2flt(f1->TMP(istart:iend,:,:))
	v1_ta		=	v1_ta + 273.15
	v1_ta@units	= "K"

	;Remove the CRU file sitting in this directory
        sCRUFileName ="cru." + istartdate + "_" + ienddate + ".nc"
	system("/bin/rm -f "+ sCRUFileName)
	;Open the CRU file for outputting the re-gridded CRU data
	fout	=	addfile(sCRUFileName,"c")

	;Regrid the CRU temperature to the RCM grid
        print("Re-gridding t2m")
	ta_regrid	= rgrid2rcm(lat_cru,lon_cru,v1_ta,lat2d,lon2d,1)
        ta = new((/ntime,1,nlat,nlon/),float)
        ta(:,0,:,:) = ta_regrid
        
	;Give ta the RCM coordinates
	ta!0 = "time"
	ta!1 = "m2"
	ta!2 = "iy"
	ta!3 = "jx"
        ta&time = v1_ta&time
        ta&m2 = 2
        ta&iy = lat2d&iy
        ta&jx = lat2d&jx
        ta@_FillValue = v1_ta@_FillValue
        ta@_FillValue = 1e-30
        ta@missing_value = ta@_FillValue

	;Write out ta
        print("Writing t2m")
	fout->t2m = ta

	;Read in the CRU precip data
        print("Reading TMP from " + time_cru(istart) + " to " + time_cru(iend))
	v2_PRE		=	short2flt(f2->PRE(istart:iend,:,:))

	;Regrid the CRU precip data to the RCM grid
        print("Re-gridding tpr")
	tpr_regrid	= rgrid2rcm(lat_cru,lon_cru,v2_PRE,lat2d,lon2d,1)
        tpr = new((/ntime,1,nlat,nlon/),float)
        tpr(:,0,:,:) = tpr_regrid 
	;Give it the RCM coordinates
	tpr!0 = "time"
	tpr!1 = "m2"
	tpr!2 = "iy"
	tpr!3 = "jx"
        tpr&time = v2_PRE&time
        tpr&m2 = 2
        tpr&iy = lat2d&iy
        tpr&jx = lat2d&jx
        tpr@_FillValue = v2_PRE@_FillValue
        tpr@_FillValue = 1e-30
        tpr@missing_value = tpr@_FillValue


	;Output the precip data
        print("Writing tpr")
	fout->tpr = tpr 
	

end

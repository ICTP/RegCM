;>
;; @par Copyright
;;
;; @brief program to generate netcdf file containing parameters specifying the plumes calculated in mo_simple_plumes
;;
;;
;; @author Bjorn Stevens, MPI-M, Hamburg (2015-11-22)
;;
;; @tuning
;;    2015-12-03, Karsten Peters, added effective radius factor
;;    2015-12-03, Bjorn Stevens, simplified plumes, revised emission scaling, removed beta_mx
;;    2016-06-20, Bjorn Stevens, revised emission scaling following CEDS beta release, and finalized v1.
;;    2016-07-09, Bjorn Stevens, revised to incorporate CEDS June 30 2016 Emission release.
;;
;; $ID: n/a$
;;
;;
;; @par Copyright
;; 
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
; 
begin
  of_name     = "MACv2.0-SP_v1.nc"
  fo          = addfile(of_name,"c")

  np          = 9
  nf          = 2
  deg2rad     = asin(1.)/90.
  twoPI       = 4.*asin(1.)
  unit_vector = new((/np/),float)
  unit_vector = 1.
;
; define data 
;
  np@long_name = "number of plumes"
  plume_number = ispan(1,np,1)
  plume_number@long_name = "plume number"
  plume_number!0 = "plume_number"
  plume_number&plume_number = plume_number

  nf@long_name = "number of features per plume"
  plume_feature = ispan(1,nf,1)
  plume_feature@long_name = "plume feature"
  plume_feature!0 = "plume_feature"
  plume_feature&plume_feature = plume_feature

  years = ispan(1850,2100,1)
  years@long_name = "years of scale factor data"
  years@unit      = "Gregorian Year"
  years!0         = "years"
  years&years     = years
  
  ny              = dimsizes(years)
  ny@long_name    = "number of years of plume scaling data"
  ny@unit         = "Year"

  year_fr = fspan(0.,1.,52)
  year_fr@long_name = "fractional year (for seasonal cycle)"
  year_fr@unit      = "Year"
  year_fr!0         = "year_fr"
  year_fr&year_fr   = year_fr

  nt                = dimsizes(year_fr)
  nt@long_name      = "Resolution of annual cycle"
  nt@unit           = "Year"

  template = new( (/np/), float)
  template!0 = "plume_number"
  template&plume_number = plume_number

  template2 = new( (/np,nf/), float)
  template2!0 = "plume_number"
  template2!1 = "plume_feature"
  template2&plume_number = plume_number
  template2&plume_feature = plume_feature

  template3 = new( (/np,ny/), float)
  template3!0 = "plume_number"
  template3!1 = "years"
  template3&plume_number = plume_number
  template3&years = years

  template4 = new( (/np,nt,nf/), float)
  template4!0 = "plume_number"
  template4!1 = "year_fr"
  template4!2 = "plume_feature"
  template4&plume_number  = plume_number
  template4&year_fr       = year_fr
  template4&plume_feature = plume_feature

  is_biomass = (/0,0,0,0,1,1,1,1,0/)
  is_biomass@long_name = "Biomass Dominated"
  is_biomass@units     = "Binary (1-yes; 0-no)"
  is_biomass!0 = "plume_number"
  is_biomass&plume_number = plume_number

  biomass_inds = ind(is_biomass.eq.1)
  indstry_inds = ind(is_biomass.ne.1)
  
  plume_lat = template
  plume_lat@long_name   = "latitude of plume center"
  plume_lat@units       = "degrees North"
  plume_lat             =   (/49.4,  40.1,  30.,  23.3,  3.5,  -10.3,   -1.,   -3.5,  -20./)

  plume_lon = template
  plume_lon@long_name   = "longitude of plume center"
  plume_lon@units       = "degrees East"
  plume_lon             =   (/20.6, 277.5, 114.,   88., 22.5,   298., 106.0,   16.0,  135./)
  
  sig_lat_W = template2
  sig_lat_W@long_name   = "Latitude extension of feature West of centre (unrotated)"
  sig_lat_W@units       = "degrees North"
  sig_lat_W(:,0)        =  (/ 6.,    7.,    6.,    9.,   6.,     6.,    8.,     9.,    6./)
  sig_lat_W(:,1)        =  (/10.,   25.,   13.,   15.,   1.,     6.,    4.,     5.,   12./)
  
  sig_lat_E = template2
  sig_lat_E@long_name   = "Latitude extension of feature East of centre (unrotated)"
  sig_lat_E@units       = "degrees South"
  sig_lat_E(:,0)        =  (/ 6.,    7.,    6.,    8.,   6.,     6.,    8.,     9.,    6./)
  sig_lat_E(:,1)        =  (/10.,    8.,   13.,   17.,   1.,     6.,    4.,     5.,   12./)
  
  sig_lon_W = template2
  sig_lon_W@long_name   = "Longitude extension of feature West of centre (unrotated)"
  sig_lon_W@units       = "degrees East"
  sig_lon_W(:,0)        =  (/  7.,  20.,    9.,    15.,   32.,   10.,   12.,   23.,   10./)
  sig_lon_W(:,1)        =  (/ 35.,   8.,    15.,   40.,    3.,    8.,    4.,    7.,   20./)

  sig_lon_E = template2
  sig_lon_E@long_name   = "Longitude extension of feature East of centre (unrotated)"
  sig_lon_E@units       = "degrees West"
  sig_lon_E(:,0)        =  (/ 13.,   35.,   8.,    10.,   6.,   10.,    10.,   14.,    4./)
  sig_lon_E(:,1)        =  (/ 80.,   11.,   40.,   15.,   3.,    8.,    6.,     6.,   20./)

  theta = template2
  theta@long_name       = "Clockwise rotation of feature central axis"
  theta@units           = "radians"
  theta(:,0)            =  (/  0.,   15.,  40.,    0.,   0.,    -30.,    10.,  -15.,    0./) * deg2rad
  theta(:,1)            =  (/  10., 130.,  15.,   15.,   0.,    -30. ,    0.,  -15.,  -30./) * deg2rad

  fwght = template2
  fwght@long_name       = "Feature Weight"
  fwght@units           = " "
  fwght(:,0)            =  (/2./5. ,3./5., 6./7., 3./5., 4./5., 1./8., 2./5., 7./10., 4./5./)
  fwght(:,1)            =  1.0 - fwght(:,0)
  
  beta_a = template
  beta_a@long_name       = "Parameter a of beta function for vertical distribution"
  beta_a@units           = " "
  beta_a                 =  (/1.5, 1.7, 1.3, 1.3, 7.0, 1.2, 2.3, 2.4, 1.4/)

  beta_b = template
  beta_b@long_name       = "Parameter b of beta function for vertical distribution"
  beta_b@units           = " "
  beta_b                 =  (/ 17., 17., 13.,  8., 35.,  9., 23., 14., 11./)
  
  aod_spmx = template   
  aod_spmx@long_name     = "Optical depth (AOD) at source"
  aod_spmx@units         = " "
  aod_spmx               = (/ 0.148, 0.094, 0.636, 0.259, 0.211, 0.351, 0.257, 0.372, 0.075/)

  aod_fmbg = template   
  aod_fmbg@long_name     = "Background fine mode optical depth (AOD) at 550 nm at source"
  aod_fmbg@units         = " "
  aod_fmbg(indstry_inds) = 0.1
  aod_fmbg(biomass_inds) = 0.6
  
  ssa550 = template   
  ssa550@long_name       = "Single Scattering Albedo at 550 nm"
  ssa550@units           = " "
  ssa550(indstry_inds)   = 0.93
  ssa550(biomass_inds)   = 0.87
   
  asy550 = template   
  asy550@long_name       = "Asymmetry Factor at 550 nm"
  asy550@units           = " "
  asy550                 = 0.63
         
  angstrom = template   
  angstrom@long_name     = "Angstrom Exponent"
  angstrom@units         = " "
  angstrom               = 2.00
 
  ann_cyc                = template4
  ann_cyc@long_name      = "Annual cycle for each of a plume's features"
  ann_cyc@units          = " "
  
  ann_cyc(0,:,0) = 1.00 - 0.25 * cos(twoPI * (year_fr - 10.5/12.))
  ann_cyc(1,:,0) = 1.00 - 0.30 * cos(twoPI * (year_fr -  0.5/12.))
  ann_cyc(2,:,0) = 1.00 + 0.15 * cos(2*(twoPI * (year_fr - 9.0/12.)))
  ann_cyc(3,:,0) = 1.00 - 0.10 * cos(   twoPI * (year_fr -  0.5/12.))
  ann_cyc(4,:,0) = linint1 ( (/0.0, 1.5, 2.5, 9.5, 10.5, 12.0/)/12.                   , (/1.00, 0.40, 0.05, 0.05, 0.50, 1.0/)                      ,False,year_fr,0)
  ann_cyc(5,:,0) = linint1 ( (/0.0, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.0/)/12.   , (/0.10, 0.10, 0.15, 0.7, 1.00, 0.7, 0.10, 0.10, 0.1/)      ,False,year_fr,0)
  ann_cyc(6,:,0) = linint1 ( (/0.0, 4.5, 5.5, 6.5, 7.5, 8.5,  9.5, 10.5, 12.0/)/12.   , (/0.10, 0.10, 0.10, 0.5, 0.8, 1.00, 0.80, 0.50, 0.1/)      ,False,year_fr,0)
  ann_cyc(7,:,0) = linint1 ( (/0.0, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12./)/12., (/0.05, 0.05, 0.1, 0.4, 0.86, 1.00, 0.8, 0.20, 0.07, 0.05/),False,year_fr,0)
  ann_cyc(8,:,0) = 1.00 + 0.60 * cos(twoPI * (year_fr - 17./24.))

  ann_cyc(0,:,1) = ann_cyc(0,:,0)
  ann_cyc(1,:,1) = ann_cyc(1,:,0)
  ann_cyc(2,:,1) = ann_cyc(2,:,0)
  ann_cyc(3,:,1) = ann_cyc(3,:,0)
  ann_cyc(4,:,1) = ann_cyc(4,:,0)
  ann_cyc(5,:,1) = linint1 ( (/0.0, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.0/)/12.    , (/0.00, 0.00, 0.00, 0.3, 1.00, 0.4, 0.00, 0.00, 0.00/)         ,False,year_fr,0)
  ann_cyc(6,:,1) = ann_cyc(6,:,0)
  ann_cyc(7,:,1) = linint1 ( (/0.0, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.0/)/12., (/0.00, 0.00, 0.00, 0.00, 0.35, 1.00, 1.00, 0.77, 0.20, 0.05/) ,False,year_fr,0)
  ann_cyc(8,:,1) = ann_cyc(8,:,0)

  year_weight = template3
  year_weight@long_name = "Annual Scaling factor for plume amplitude"
  year_weight@unit      = "Amplitude relative to 2005, scaled by CEDS emissions"

  years_in = new( (/20/), "float")
  years_in(:)      =                  (/ 1850, 1860, 1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2005, 2010, 2014, 2016/)

  year_weight(0,:) = linint1(years_in,(/ 0.00, 0.05, 0.13, 0.22, 0.32, 0.47, 0.62, 0.57, 0.65, 0.81, 1.01, 1.75, 2.50, 2.94, 2.35, 1.23, 1.00, 0.86, 0.71, 0.71/), False, years, 0)
  year_weight(1,:) = linint1(years_in,(/ 0.00, 0.01, 0.04, 0.10, 0.21, 0.42, 0.78, 0.98, 0.94, 1.12, 1.25, 1.57, 1.85, 1.62, 1.42, 1.16, 1.00, 0.70, 0.53, 0.53/), False, years, 0)
  year_weight(2,:) = linint1(years_in,(/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.01, 0.02, 0.03, 0.04, 0.06, 0.07, 0.20, 0.33, 0.46, 0.62, 0.79, 1.00, 1.02, 0.99, 0.99/), False, years, 0)
  year_weight(3,:) = linint1(years_in,(/ 0.00, 0.00, 0.00, 0.01, 0.01, 0.02, 0.02, 0.03, 0.04, 0.07, 0.10, 0.20, 0.37, 0.55, 0.72, 0.87, 1.00, 1.17, 1.28, 1.28/), False, years, 0)
  year_weight(4,:) = linint1(years_in,(/ 0.00, 0.00, 0.01, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.09, 0.14, 0.22, 0.33, 0.53, 0.71, 0.88, 1.00, 1.14, 1.23, 1.23/), False, years, 0)
  year_weight(5,:) = linint1(years_in,(/ 0.00, 0.00, 0.01, 0.01, 0.02, 0.03, 0.06, 0.10, 0.16, 0.24, 0.28, 0.44, 0.66, 0.95, 1.04, 1.05, 1.00, 1.03, 1.09, 1.09/), False, years, 0)
  year_weight(6,:) = linint1(years_in,(/ 0.00, 0.00, 0.00, 0.01, 0.01, 0.02, 0.03, 0.04, 0.06, 0.07, 0.09, 0.17, 0.34, 0.54, 0.92, 1.04, 1.00, 1.06, 1.10, 1.10/), False, years, 0)
  year_weight(7,:) = linint1(years_in,(/ 0.00, 0.00, 0.01, 0.01, 0.02, 0.03, 0.06, 0.09, 0.14, 0.26, 0.36, 0.62, 0.89, 1.02, 0.94, 0.93, 1.00, 1.16, 1.29, 1.29/), False, years, 0)
  year_weight(8,:) = linint1(years_in,(/ 0.00, 0.01, 0.02, 0.04, 0.07, 0.14, 0.22, 0.21, 0.26, 0.33, 0.39, 0.60, 0.81, 0.89, 0.99, 0.95, 1.00, 0.98, 0.94, 0.94/), False, years, 0)

  fo@History       = "File created by sp_make-data_v1.ncl on "+systemfunc("date")
  fo@Author        = "Bjorn Stevens"
  fo@Reference     = "Stevens et al., 2016: Simple Plumes ... GMD"
  fo@plume1_region = "Europe"
  fo@plume2_region = "North America"
  fo@plume3_region = "East Asia"
  fo@plume4_region = "South Asia"
  fo@plume5_region = "North Africa (biomass)"
  fo@plume6_region = "South America (biomass)"
  fo@plume7_region = "Maritime Continent (biomass)"
  fo@plume8_region = "South Central Africa (biomass)"
  fo@plume9_region = "Australia"

  fo->plume_number  = plume_number
  fo->plume_feature = plume_feature
  fo->years         = years
  fo->year_fr       = year_fr
  
  fo->is_biomass = is_biomass
  fo->plume_lat  = plume_lat
  fo->plume_lon  = plume_lon
  fo->sig_lat_W  = sig_lat_W
  fo->sig_lat_E  = sig_lat_E
  fo->sig_lon_E  = sig_lon_E
  fo->sig_lon_W  = sig_lon_W
  fo->theta      = theta
  fo->ftr_weight = fwght
  fo->beta_a     = beta_a
  fo->beta_b     = beta_b
  fo->aod_spmx   = aod_spmx
  fo->aod_fmbg   = aod_fmbg
  fo->ssa550     = ssa550
  fo->asy550     = asy550
  fo->angstrom   = angstrom
  fo->ann_cycle  = ann_cyc
  fo->year_weight= year_weight
    
end
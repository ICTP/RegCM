
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine myabort
  implicit none
  stop ' Execution terminated because of runtime error'
end subroutine myabort

program interpinic

#ifndef CLM45
  use mod_stdio
  write(stdout,*) 'RegCM built without CLM45 support.'
  write(stdout,*) 'Please recompile it using --enable-clm45 flag'
#else

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_memutil
  use mod_message
  use mod_kdinterp
  use netcdf

  implicit none

  character (len=256) :: prgname
  character (len=256) :: inputfile
  character (len=256) :: outputfile

  integer(ik4) :: istat
  integer(ik4) :: ncin , ncout

  integer(ik4) :: hostnm
  integer(ik4) :: ihost , idir
  integer(ik4) :: getcwd
  integer(ik4) :: ingc , inlu , inco , inpft , ingrd
  integer(ik4) :: ongc , onlu , onco , onpft , ongrd

  integer(ik4) , parameter :: itypeveg = 1
  real(rkx) , parameter :: missl = -9999.0_rkx

  real(rkx) , pointer , dimension(:) :: i_grid1d_lon
  real(rkx) , pointer , dimension(:) :: i_grid1d_lat
  real(rkx) , pointer , dimension(:) :: i_pft1d_lon
  real(rkx) , pointer , dimension(:) :: i_pft1d_lat
  real(rkx) , pointer , dimension(:) :: i_col1d_lon
  real(rkx) , pointer , dimension(:) :: i_col1d_lat
  real(rkx) , pointer , dimension(:) :: o_grid1d_lon
  real(rkx) , pointer , dimension(:) :: o_grid1d_lat
  real(rkx) , pointer , dimension(:) :: o_pft1d_lon
  real(rkx) , pointer , dimension(:) :: o_pft1d_lat
  real(rkx) , pointer , dimension(:) :: o_col1d_lon
  real(rkx) , pointer , dimension(:) :: o_col1d_lat
  real(rkx) , pointer , dimension(:) :: o_pfts1d_wtxy
  real(rkx) , pointer , dimension(:) :: o_cols1d_wtxy
  real(rkx) , pointer , dimension(:) :: i_gval
  real(rkx) , pointer , dimension(:) :: o_gval
  real(rkx) , pointer , dimension(:) :: i_pval
  real(rkx) , pointer , dimension(:) :: o_pval
  real(rkx) , pointer , dimension(:) :: i_cval
  real(rkx) , pointer , dimension(:,:) :: i_cval_ng
  real(rkx) , pointer , dimension(:,:) :: i_cval_ng_t
  real(rkx) , pointer , dimension(:) :: o_cval
  real(rkx) , pointer , dimension(:,:) :: o_cval_ng
  real(rkx) , pointer , dimension(:,:) :: o_cval_ng_t
  integer(ik4) , pointer , dimension(:) :: i_ltype
  integer(ik4) , pointer , dimension(:) :: o_ltype
  integer(ik4) , pointer , dimension(:) :: i_vtype
  integer(ik4) , pointer , dimension(:) :: o_vtype
  integer(ik4) , pointer , dimension(:) :: i_ctype
  integer(ik4) , pointer , dimension(:) :: o_ctype
  integer(ik4) , pointer , dimension(:,:) :: i_mapf
  integer(ik4) , pointer , dimension(:,:) :: o_mapf
  integer(ik4) , pointer , dimension(:,:) :: i_mapc
  integer(ik4) , pointer , dimension(:,:) :: o_mapc
  integer(ik4) , pointer , dimension(:) :: colval

  type(h_interpolator) :: hint
  integer(ik4) :: imaxpft , omaxpft
  integer(ik4) :: imaxlun , omaxlun
  integer(ik4) :: imaxcol , omaxcol , maxcol

  integer(ik4) , dimension(8) :: tval
  character (len=32) :: cdata='?'
  character (len=5) :: czone='?'
  character (len=32) :: hostname='?'
  character (len=32) :: user='?'
  character (len=128) :: directory='?'
  character (len=*) , parameter :: f99001 = &
          '(2x," GIT Revision: ",a," compiled at: data : ",a,"  time: ",a,/)'

  write (stdout,  &
     "(/,2x,'This is interpinic part of RegCM package version 4')")
  write (stdout,f99001)  GIT_VER, __DATE__ , __TIME__

#ifdef IBM
  hostname='ibm platform '
  user= 'Unknown'
#else
  ihost = hostnm(hostname)
  call getlog(user)
#endif
  call date_and_time(zone=czone,values=tval)
  idir = getcwd(directory)

  write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
     tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
  write(stdout,*) ": this run start at  : ",trim(cdata)
  write(stdout,*) ": it is submitted by : ",trim(user)
  write(stdout,*) ": it is running on   : ",trim(hostname)
  write(stdout,*) ": in directory       : ",trim(directory)
  write(stdout,*) "                     "

  call get_command_argument(0,value=prgname)
  istat = command_argument_count()
  if ( istat < 2 ) then
    write(stderr, *) 'Usage : ', trim(prgname) , &
              ' inputfile.r.nc outputfile.r.nc'
    call die(__FILE__,'Input file name missing',__LINE__)
  end if
  call get_command_argument(1,value=inputfile,status=istat)
  if ( istat < 0 ) then
    write(stderr, *) 'Usage : ', trim(prgname) , &
              ' inputfile.r.nc outputfile.r.nc'
    call die(__FILE__,'Input file name missing',__LINE__)
  end if
  call get_command_argument(2,value=outputfile,status=istat)
  if ( istat < 0 ) then
    write(stderr, *) 'Usage : ', trim(prgname) , &
              ' inputfile.r.nc outputfile.r.nc'
    call die(__FILE__,'Output file name missing',__LINE__)
  end if

  istat = nf90_open(inputfile, nf90_nowrite, ncin)
  if ( istat /= nf90_noerr ) then
    write(stderr, *) 'Usage : ', trim(prgname) , &
              ' inputfile.r.nc outputfile.r.nc'
    write (stderr,*) trim(inputfile)//' : '//nf90_strerror(istat)
    call die(__FILE__,'Input file open error',__LINE__)
  end if
  istat = nf90_open(outputfile, nf90_write, ncout)
  if ( istat /= nf90_noerr ) then
    write(stderr, *) 'Usage : ', trim(prgname) , &
              ' inputfile.r.nc outputfile.r.nc'
    write (stderr,*) trim(outputfile)//' : '//nf90_strerror(istat)
    call die(__FILE__,'Output file open error',__LINE__)
  end if

  call memory_init

  ingc = dlen(ncin,'gridcell')
  inlu = dlen(ncin,'landunit')
  inco = dlen(ncin,'column')
  inpft = dlen(ncin,'pft')
  ingrd = dlen(ncin,'levgrnd')

  write(stdout,*) 'Input DIMENSIONS'
  write(stdout,*) 'Input gridcell  : ',ingc
  write(stdout,*) 'Input landunit  : ',inlu
  write(stdout,*) 'Input column    : ',inco
  write(stdout,*) 'Input pft       : ',inpft
  write(stdout,*) 'Input levgrnd   : ',ingrd
  write(stdout,*) ''

  ongc = dlen(ncout,'gridcell')
  onlu = dlen(ncout,'landunit')
  onco = dlen(ncout,'column')
  onpft = dlen(ncout,'pft')
  ongrd = dlen(ncout,'levgrnd')

  write(stdout,*) 'Output DIMENSIONS'
  write(stdout,*) 'Output gridcell : ',ongc
  write(stdout,*) 'Output landunit : ',onlu
  write(stdout,*) 'Output column   : ',onco
  write(stdout,*) 'Output pft      : ',onpft
  write(stdout,*) 'Output levgrnd  : ',ongrd
  write(stdout,*) ''

  if ( ingrd /= ongrd ) then
    write (stderr,*) 'ERROR: levgrnd dimension do not match'
    call die(__FILE__,'File read error',__LINE__)
  end if

  call getmem1d(i_grid1d_lon,1,ingc,'interpinic:i_grid1d_lon')
  call getmem1d(i_grid1d_lat,1,ingc,'interpinic:i_grid1d_lat')
  call getmem1d(i_gval,1,ingc,'interpinic:i_gval')
  call getmem1d(i_pft1d_lon,1,inpft,'interpinic:i_pft1d_lon')
  call getmem1d(i_pft1d_lat,1,inpft,'interpinic:i_pft1d_lat')
  call getmem1d(i_col1d_lon,1,inco,'interpinic:i_col1d_lon')
  call getmem1d(i_col1d_lat,1,inco,'interpinic:i_col1d_lat')
  call getmem1d(i_pval,1,inpft,'interpinic:i_pval')
  call getmem1d(i_cval,1,inco,'interpinic:i_cval')
  call getmem2d(i_cval_ng,1,ingrd,1,inco,'interpinic:i_cval_ng')
  call getmem2d(i_cval_ng_t,1,inco,1,ingrd,'interpinic:i_cval_ng_t')
  call getmem1d(i_ctype,1,inco,'interpinic:i_ctype')
  call getmem1d(i_ltype,1,inpft,'interpinic:i_ltype')
  call getmem1d(i_vtype,1,inpft,'interpinic:i_vtype')
  call getmem1d(o_grid1d_lon,1,ongc,'interpinic:o_grid1d_lon')
  call getmem1d(o_grid1d_lat,1,ongc,'interpinic:o_grid1d_lat')
  call getmem1d(o_gval,1,ongc,'interpinic:o_gval')
  call getmem1d(o_pft1d_lon,1,onpft,'interpinic:o_pft1d_lon')
  call getmem1d(o_pft1d_lat,1,onpft,'interpinic:o_pft1d_lat')
  call getmem1d(o_col1d_lon,1,onco,'interpinic:o_col1d_lon')
  call getmem1d(o_col1d_lat,1,onco,'interpinic:o_col1d_lat')
  call getmem1d(o_pfts1d_wtxy,1,onpft,'interpinic:o_pfts1d_wtxy')
  call getmem1d(o_cols1d_wtxy,1,onco,'interpinic:o_cols1d_wtxy')
  call getmem1d(o_pval,1,onpft,'interpinic:o_pval')
  call getmem1d(o_cval,1,onco,'interpinic:o_cval')
  call getmem2d(o_cval_ng,1,ongrd,1,onco,'interpinic:o_cval_ng')
  call getmem2d(o_cval_ng_t,1,onco,1,ongrd,'interpinic:o_cval_ng_t')
  call getmem1d(o_ctype,1,onco,'interpinic:o_ctype')
  call getmem1d(o_ltype,1,onpft,'interpinic:o_ltype')
  call getmem1d(o_vtype,1,onpft,'interpinic:o_vtype')

  call rval(ncin,'grid1d_lon',i_grid1d_lon)
  call rval(ncin,'grid1d_lat',i_grid1d_lat)
  call rval(ncin,'pfts1d_lon',i_pft1d_lon)
  call rval(ncin,'pfts1d_lat',i_pft1d_lat)
  call rval(ncin,'cols1d_lon',i_col1d_lon)
  call rval(ncin,'cols1d_lat',i_col1d_lat)
  call ival(ncin,'cols1d_ityp',i_ctype)
  call ival(ncin,'pfts1d_ityplun',i_ltype)
  call ival(ncin,'pfts1d_itypveg',i_vtype)
  call rval(ncout,'grid1d_lon',o_grid1d_lon)
  call rval(ncout,'grid1d_lat',o_grid1d_lat)
  call rval(ncout,'pfts1d_lon',o_pft1d_lon)
  call rval(ncout,'pfts1d_lat',o_pft1d_lat)
  call rval(ncout,'cols1d_lon',o_col1d_lon)
  call rval(ncout,'cols1d_lat',o_col1d_lat)
  call ival(ncout,'cols1d_ityp',o_ctype)
  call ival(ncout,'pfts1d_ityplun',o_ltype)
  call ival(ncout,'pfts1d_itypveg',o_vtype)
  call rval(ncout,'pfts1d_wtxy',o_pfts1d_wtxy)
  call rval(ncout,'cols1d_wtxy',o_cols1d_wtxy)

  imaxpft = maxval(i_vtype) + 1 ! The bare ground class
  imaxlun = maxval(i_ltype)
  imaxcol = countcols(i_ctype)
  omaxpft = maxval(o_vtype) + 1 ! The bare ground class
  omaxlun = maxval(o_ltype)
  omaxcol = countcols(o_ctype)

  call set_colvals(i_ctype,o_ctype,colval)
  maxcol = size(colval)

  write(stdout,*) 'Input  PFT number : ', imaxpft
  write(stdout,*) 'Output PFT number : ', omaxpft
  write(stdout,*) 'Input  LUN number : ', imaxlun
  write(stdout,*) 'Output LUN number : ', omaxlun
  write(stdout,*) 'Input  COL number : ', imaxcol
  write(stdout,*) 'Output COL number : ', omaxcol

  if ( imaxpft /= omaxpft ) then
    write(stderr,*) 'Number of PFT do not match !'
    call die(__FILE__,'CAN WORK ONLY IF THE PFT MATCH.',__LINE__)
  end if

  call getmem2d(i_mapf,1,ingc,1,imaxpft,'interpinic:i_mapf')
  call getmem2d(o_mapf,1,ongc,1,omaxpft,'interpinic:o_mapf')
  call getmem2d(i_mapc,1,ingc,1,maxcol,'interpinic:i_mapc')
  call getmem2d(o_mapc,1,ongc,1,maxcol,'interpinic:o_mapc')

  call h_interpolator_create(hint, &
                       ingc, i_grid1d_lat, i_grid1d_lon, &
                       ongc, o_grid1d_lat, o_grid1d_lon)

  write(stdout,*) 'Input/output grids successfully read.'

  call makemap_pft(inpft,ingc,imaxpft,i_ltype,i_vtype, &
                   i_grid1d_lat,i_grid1d_lon,i_pft1d_lat,i_pft1d_lon,i_mapf)
  call makemap_col(inco,ingc,maxcol,i_ctype, &
                   i_grid1d_lat,i_grid1d_lon,i_col1d_lat,i_col1d_lon,i_mapc)
  write(stdout,*) 'Input Map created.'

  call makemap_pft(onpft,ongc,omaxpft,o_ltype,o_vtype, &
                   o_grid1d_lat,o_grid1d_lon,o_pft1d_lat,o_pft1d_lon,o_mapf)
  call makemap_col(onco,ongc,maxcol,o_ctype, &
                   o_grid1d_lat,o_grid1d_lon,o_col1d_lat,o_col1d_lon,o_mapc)
  write(stdout,*) 'Output Map created.'

  call col_interpolate('fpg')
  call col_interpolate('fpi')
  call col_interpolate('irrig_rate')
  call col_interpolate('cannsum_npp')
  call col_interpolate('col_lag_npp')
  call col_interpolate('lfc')
  call col_interpolate('wf')
  call col_interpolate('farea_burned')
  call col_interpolate('baf_crop')
  call col_interpolate('baf_peatf')
  call col_interpolate('fbac')
  call col_interpolate('fbac1')
  call col_interpolate('altmax')
  call col_interpolate('altmax_lastyear')
  call col_interpolate('litr1c')
  call col_interpolate('litr2c')
  call col_interpolate('litr3c')
  call col_interpolate('soil1c')
  call col_interpolate('soil2c')
  call col_interpolate('soil4c')
  call col_interpolate('soil4c')
  call col_interpolate('cwdc')
  call col_interpolate('col_ctrunc')
  call col_interpolate('seedc')
  call col_interpolate('totlitc')
  call col_interpolate('totcolc')
  call col_interpolate('prod10c')
  call col_interpolate('prod100c')
  call col_interpolate('totsomc')
  call col_interpolate('litr1c_13')
  call col_interpolate('litr2c_13')
  call col_interpolate('litr3c_13')
  call col_interpolate('soil1c_13')
  call col_interpolate('soil2c_13')
  call col_interpolate('soil4c_13')
  call col_interpolate('soil4c_13')
  call col_interpolate('seedc_13')
  call col_interpolate('totlitc_13')
  call col_interpolate('totcolc_13')
  call col_interpolate('prod10c_13')
  call col_interpolate('prod100c_13')
  call col_interpolate('litr1c_14')
  call col_interpolate('litr2c_14')
  call col_interpolate('litr3c_14')
  call col_interpolate('soil1c_14')
  call col_interpolate('soil2c_14')
  call col_interpolate('soil4c_14')
  call col_interpolate('soil4c_14')
  call col_interpolate('seedc_14')
  call col_interpolate('totlitc_14')
  call col_interpolate('totcolc_14')
  call col_interpolate('prod10c_14')
  call col_interpolate('prod100c_14')
  call col_interpolate('litr1n')
  call col_interpolate('litr2n')
  call col_interpolate('litr3n')
  call col_interpolate('soil1n')
  call col_interpolate('soil2n')
  call col_interpolate('soil4n')
  call col_interpolate('soil4n')
  call col_interpolate('cwdn')
  call col_interpolate('col_ntrunc')
  call col_interpolate('totcoln')
  call col_interpolate('seedn')
  call col_interpolate('sminn')
  call col_interpolate('prod10n')
  call col_interpolate('prod100n')
  call col_interpolate('QFLX_SURF_LAG')
  call col_interpolate('FINUNDATED_LAG')
  call col_interpolate('FINUNDATED')
  call col_interpolate('Z0MG')
  call col_interpolate('annavg_somhr')
  call col_interpolate('annavg_finrw')

  call col_interpolate_grnd('fpi_vr')
  call col_interpolate_grnd('som_adv_coef_vr')
  call col_interpolate_grnd('som_diffus_coef_vr')
  call col_interpolate_grnd('litr1c_vr')
  call col_interpolate_grnd('litr2c_vr')
  call col_interpolate_grnd('litr3c_vr')
  call col_interpolate_grnd('cwdc_vr')
  call col_interpolate_grnd('soil1c_vr')
  call col_interpolate_grnd('soil2c_vr')
  call col_interpolate_grnd('soil3c_vr')
  call col_interpolate_grnd('col_ctrunc_vr')
  call col_interpolate_grnd('sminn_vr')
  call col_interpolate_grnd('litr1n_vr')
  call col_interpolate_grnd('litr2n_vr')
  call col_interpolate_grnd('litr3n_vr')
  call col_interpolate_grnd('cwdn_vr')
  call col_interpolate_grnd('soil1n_vr')
  call col_interpolate_grnd('soil2n_vr')
  call col_interpolate_grnd('soil3n_vr')
  call col_interpolate_grnd('col_ntrunc_vr')
  call col_interpolate_grnd('f_nit_vr_vr')
  call col_interpolate_grnd('pot_f_nit_vr_vr')
  call col_interpolate_grnd('smin_no3_vr')
  call col_interpolate_grnd('smin_nh4_vr')
  call col_interpolate_grnd('CONC_CH4_SAT')
  call col_interpolate_grnd('CONC_CH4_UNSAT')
  call col_interpolate_grnd('CONC_O2_SAT')
  call col_interpolate_grnd('CONC_O2_UNSAT')
  call col_interpolate_grnd('O2STRESS_SAT')
  call col_interpolate_grnd('O2STRESS_UNSAT')
  call col_interpolate_grnd('LAYER_SAT_LAG')
  call col_interpolate_grnd('O2_DECOMP_DEPTH_SAT')
  call col_interpolate_grnd('O2_DECOMP_DEPTH_UNSAT')
  call col_interpolate_grnd('LAKE_SOILC')

#ifdef CNDV
  call pft_interpolate('elai')
  call pft_interpolate('esai')
  call pft_interpolate('tlai')
  call pft_interpolate('tsai')
  call pft_interpolate('mlaidiff')
  call pft_interpolate('htop')
  call pft_interpolate('hbot')
#endif
  call pft_interpolate('annavg_t2m')
  call pft_interpolate('dormant_flag')
  call pft_interpolate('days_active')
  call pft_interpolate('onset_flag')
  call pft_interpolate('onset_counter')
  call pft_interpolate('onset_gddflag')
  call pft_interpolate('onset_fdd')
  call pft_interpolate('onset_gdd')
  call pft_interpolate('onset_swi')
  call pft_interpolate('offset_flag')
  call pft_interpolate('offset_counter')
  call pft_interpolate('offset_fdd')
  call pft_interpolate('offset_swi')
  call pft_interpolate('fert_counter')
  call pft_interpolate('fert')
  call pft_interpolate('lgsf')
  call pft_interpolate('bglfr')
  call pft_interpolate('bgtr')
  call pft_interpolate('gpp_pepv')
  call pft_interpolate('availc')
  call pft_interpolate('xsmrpool_recover')
  call pft_interpolate('xsmrpool_c13ratio')
  call pft_interpolate('alloc_pnow')
  call pft_interpolate('c_allometry')
  call pft_interpolate('n_allometry')
  call pft_interpolate('plant_ndemand')
  call pft_interpolate('annsum_potential_gpp')
  call pft_interpolate('annmax_retransn')
  call pft_interpolate('avail_retransn')
  call pft_interpolate('plant_nalloc')
  call pft_interpolate('plant_calloc')
  call pft_interpolate('excess_cflux')
  call pft_interpolate('downreg')
  call pft_interpolate('prev_leafc_to_litter')
  call pft_interpolate('prev_frootc_to_litter')
  call pft_interpolate('annsum_npp')
  call pft_interpolate('rc13_canair')
  call pft_interpolate('rc13_psnsun')
  call pft_interpolate('rc13_psnsha')
  call pft_interpolate('grain_flag')
  call pft_interpolate('leafc')
  call pft_interpolate('leafc_storage')
  call pft_interpolate('leafc_xfer')
  call pft_interpolate('frootc')
  call pft_interpolate('frootc_storage')
  call pft_interpolate('frootc_xfer')
  call pft_interpolate('livestemc')
  call pft_interpolate('livestemc_storage')
  call pft_interpolate('livestemc_xfer')
  call pft_interpolate('deadstemc')
  call pft_interpolate('deadstemc_storage')
  call pft_interpolate('deadstemc_xfer')
  call pft_interpolate('livecrootc')
  call pft_interpolate('livecrootc_storage')
  call pft_interpolate('livecrootc_xfer')
  call pft_interpolate('deadcrootc')
  call pft_interpolate('deadcrootc_storage')
  call pft_interpolate('deadcrootc_xfer')
  call pft_interpolate('gresp_storage')
  call pft_interpolate('gresp_xfer')
  call pft_interpolate('cpool')
  call pft_interpolate('xsmrpool')
  call pft_interpolate('pft_ctrunc')
  call pft_interpolate('totvegc')
  call pft_interpolate('leafc_13')
  call pft_interpolate('leafc_storage_13')
  call pft_interpolate('leafc_xfer_13')
  call pft_interpolate('frootc_13')
  call pft_interpolate('frootc_storage_13')
  call pft_interpolate('frootc_xfer_13')
  call pft_interpolate('livestemc_13')
  call pft_interpolate('livestemc_storage_13')
  call pft_interpolate('livestemc_xfer_13')
  call pft_interpolate('deadstemc_13')
  call pft_interpolate('deadstemc_storage_13')
  call pft_interpolate('deadstemc_xfer_13')
  call pft_interpolate('livecrootc_13')
  call pft_interpolate('livecrootc_storage_13')
  call pft_interpolate('livecrootc_xfer_13')
  call pft_interpolate('deadcrootc_13')
  call pft_interpolate('deadcrootc_storage_13')
  call pft_interpolate('deadcrootc_xfer_13')
  call pft_interpolate('gresp_storage_13')
  call pft_interpolate('gresp_xfer_13')
  call pft_interpolate('cpool_13')
  call pft_interpolate('xsmrpool_13')
  call pft_interpolate('pft_ctrunc_13')
  call pft_interpolate('totvegc_13')
  call pft_interpolate('leafc_14')
  call pft_interpolate('leafc_storage_14')
  call pft_interpolate('leafc_xfer_14')
  call pft_interpolate('frootc_14')
  call pft_interpolate('frootc_storage_14')
  call pft_interpolate('frootc_xfer_14')
  call pft_interpolate('livestemc_14')
  call pft_interpolate('livestemc_storage_14')
  call pft_interpolate('livestemc_xfer_14')
  call pft_interpolate('deadstemc_14')
  call pft_interpolate('deadstemc_storage_14')
  call pft_interpolate('deadstemc_xfer_14')
  call pft_interpolate('livecrootc_14')
  call pft_interpolate('livecrootc_storage_14')
  call pft_interpolate('livecrootc_xfer_14')
  call pft_interpolate('deadcrootc_14')
  call pft_interpolate('deadcrootc_storage_14')
  call pft_interpolate('deadcrootc_xfer_14')
  call pft_interpolate('gresp_storage_14')
  call pft_interpolate('gresp_xfer_14')
  call pft_interpolate('cpool_14')
  call pft_interpolate('xsmrpool_14')
  call pft_interpolate('pft_ctrunc_14')
  call pft_interpolate('totvegc_14')
  call pft_interpolate('leafn')
  call pft_interpolate('leafn_storage')
  call pft_interpolate('leafn_xfer')
  call pft_interpolate('frootn')
  call pft_interpolate('frootn_storage')
  call pft_interpolate('frootn_xfer')
  call pft_interpolate('livestemn')
  call pft_interpolate('livestemn_storage')
  call pft_interpolate('livestemn_xfer')
  call pft_interpolate('deadstemn')
  call pft_interpolate('deadstemn_storage')
  call pft_interpolate('deadstemn_xfer')
  call pft_interpolate('livecrootn')
  call pft_interpolate('livecrootn_storage')
  call pft_interpolate('livecrootn_xfer')
  call pft_interpolate('deadcrootn')
  call pft_interpolate('deadcrootn_storage')
  call pft_interpolate('deadcrootn_xfer')
  call pft_interpolate('retransn')
  call pft_interpolate('npool')
  call pft_interpolate('pft_ntrunc')
  call pft_interpolate('btran2')
  call pft_interpolate('aleaf')
  call pft_interpolate('aleafi')
  call pft_interpolate('astem')
  call pft_interpolate('astemi')
  call pft_interpolate('htmx')
  call pft_interpolate('hdidx')
  call pft_interpolate('vf')
  call pft_interpolate('cumvd')
  call pft_interpolate('gdd1020')
  call pft_interpolate('gdd820')
  call pft_interpolate('gdd020')
  call pft_interpolate('gddmaturity')
  call pft_interpolate('huileaf')
  call pft_interpolate('huigrain')
  call pft_interpolate('grainc')
  call pft_interpolate('grainc_storage')
  call pft_interpolate('grainc_xfer')
  call pft_interpolate('grainn')
  call pft_interpolate('grainn_storage')
  call pft_interpolate('grainn_xfer')
  call pft_interpolate('grainc_xfer_to_grainc')
  call pft_interpolate('livestemc_to_litter')
  call pft_interpolate('grainc_to_food')
  call pft_interpolate('grainn_xfer_to_grainn')
  call pft_interpolate('livestemn_to_litter')
  call pft_interpolate('grainn_to_food')
  call pft_interpolate('cpool_to_grainc')
  call pft_interpolate('cpool_to_grainc_storage')
  call pft_interpolate('npool_to_grainn')
  call pft_interpolate('npool_to_grainn_storage')
  call pft_interpolate('cpool_grain_gr')
  call pft_interpolate('cpool_grain_storage_gr')
  call pft_interpolate('transfer_grain_gr')
  call pft_interpolate('grainc_storage_to_xfer')
  call pft_interpolate('grainn_storage_to_xfer')

  call pft_interpolate('CROWNAREA')
  call pft_interpolate('annsum_litfall')
  call pft_interpolate('T_MO_MIN')
  call pft_interpolate('leafcmax')

  call pft_interpolate('T_VEG240_VALUE')
  call pft_interpolate('FSD240_VALUE')
  call pft_interpolate('FSI240_VALUE')
  call pft_interpolate('FSUN240_VALUE')
  call pft_interpolate('PREC365_VALUE')
  call pft_interpolate('AGDDTW_VALUE')
  call pft_interpolate('AGDD_VALUE')
  call pft_interpolate('GDD0_VALUE')
  call pft_interpolate('GDD8_VALUE')
  call pft_interpolate('GDD10_VALUE')
  call pft_interpolate('GDDPLANT_VALUE')
  call pft_interpolate('GDDTSOI_VALUE')

  call pft_interpolate('annavg_agnpp')
  call pft_interpolate('annavg_bgnpp')

  call h_interpolator_destroy(hint)

  call memory_destroy

  istat = nf90_close(ncin)
  if ( istat /= nf90_noerr ) then
    write (stderr,*) nf90_strerror(istat)
    call die(__FILE__,'Input file close error',__LINE__)
  end if
  istat = nf90_close(ncout)
  if ( istat /= nf90_noerr ) then
    write (stderr,*) nf90_strerror(istat)
    call die(__FILE__,'Output file close error',__LINE__)
  end if

  write(stdout,*) 'Successfully interpolated CLM initial condition on restart.'

  contains

  integer(ik4) function dlen(ncid,dname)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: dname
    integer(ik4) :: idimid
    istat = nf90_inq_dimid(ncid, dname, idimid)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Inquire dimension ',trim(dname), &
                 ' : ',nf90_strerror(istat)
      call die(__FILE__,'File read error',__LINE__)
    end if
    istat = nf90_inquire_dimension(ncid, idimid, len=dlen)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Read dimension ',trim(dname),' : ',nf90_strerror(istat)
      call die(__FILE__,'File read error',__LINE__)
    end if
  end function dlen

  logical function hasval(nc1,nc2,vname)
    implicit none
    integer(ik4) , intent(in) :: nc1 , nc2
    character(len=*) , intent(in) :: vname
    integer(ik4) :: ivarid
    hasval = .true.
    istat = nf90_inq_varid(nc1, vname, ivarid)
    if ( istat /= nf90_noerr ) then
      hasval = .false.
    end if
    istat = nf90_inq_varid(nc2, vname, ivarid)
    if ( istat /= nf90_noerr ) then
      hasval = .false.
    end if
  end function hasval

  subroutine rval(ncid,vname,var)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rkx) , dimension(:) , intent(inout) :: var
    integer(ik4) :: ivarid
    istat = nf90_inq_varid(ncid, vname, ivarid)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Inquire variable ',trim(vname), &
                 ' : ',nf90_strerror(istat)
      call die(__FILE__,'File read error',__LINE__)
    end if
    istat = nf90_get_var(ncid, ivarid, var)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Read variable ',trim(vname),' : ',nf90_strerror(istat)
      call die(__FILE__,'File read error',__LINE__)
    end if
  end subroutine rval

  subroutine rval2(ncid,vname,var)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rkx) , dimension(:,:) , intent(inout) :: var
    integer(ik4) :: ivarid
    istat = nf90_inq_varid(ncid, vname, ivarid)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Inquire variable ',trim(vname), &
                 ' : ',nf90_strerror(istat)
      call die(__FILE__,'File read error',__LINE__)
    end if
    istat = nf90_get_var(ncid, ivarid, var)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Read variable ',trim(vname),' : ',nf90_strerror(istat)
      call die(__FILE__,'File read error',__LINE__)
    end if
  end subroutine rval2

  subroutine ival(ncid,vname,var)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:) , intent(inout) :: var
    integer(ik4) :: ivarid
    istat = nf90_inq_varid(ncid, vname, ivarid)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Inquire variable ',trim(vname), &
                 ' : ',nf90_strerror(istat)
      call die(__FILE__,'File read error',__LINE__)
    end if
    istat = nf90_get_var(ncid, ivarid, var)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Read variable ',trim(vname),' : ',nf90_strerror(istat)
      call die(__FILE__,'File read error',__LINE__)
    end if
  end subroutine ival

  subroutine wval(ncid,vname,var)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rkx) , dimension(:) , intent(in) :: var
    integer(ik4) :: ivarid
    istat = nf90_inq_varid(ncid, vname, ivarid)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Inquire variable ',trim(vname), &
                 ' : ',nf90_strerror(istat)
      call die(__FILE__,'File write error',__LINE__)
    end if
    istat = nf90_put_var(ncid, ivarid, var)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Write variable ',trim(vname),' : ',nf90_strerror(istat)
      call die(__FILE__,'File write error',__LINE__)
    end if
  end subroutine wval

  subroutine wval2(ncid,vname,var)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rkx) , dimension(:,:) , intent(in) :: var
    integer(ik4) :: ivarid
    istat = nf90_inq_varid(ncid, vname, ivarid)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Inquire variable ',trim(vname), &
                 ' : ',nf90_strerror(istat)
      call die(__FILE__,'File write error',__LINE__)
    end if
    istat = nf90_put_var(ncid, ivarid, var)
    if ( istat /= nf90_noerr ) then
      write (stderr,*) 'Write variable ',trim(vname),' : ',nf90_strerror(istat)
      call die(__FILE__,'File write error',__LINE__)
    end if
  end subroutine wval2

  subroutine makemap_pft(np,ng,mp,ltype,ptype,glat,glon,plat,plon,map)
    implicit none
    integer(ik4) , intent(in) :: np , ng , mp
    integer(ik4) , intent(in) , dimension(np) :: ltype
    integer(ik4) , intent(in) , dimension(np) :: ptype
    real(rkx) , intent(in) , dimension(ng) :: glat
    real(rkx) , intent(in) , dimension(ng) :: glon
    real(rkx) , intent(in) , dimension(np) :: plat
    real(rkx) , intent(in) , dimension(np) :: plon
    integer(ik4) , intent(out) , dimension(ng,mp) :: map
    integer(ik4) :: ig , ip

    map(:,:) = -1
    ig = 1
    ip = 1
    ploop: &
    do while ( ip <= np )
      if ( abs(glat(ig)-plat(ip)) > epsilon(1.0) .and. &
           abs(glon(ig)-plon(ip)) > epsilon(1.0) ) then
        ig = ig + 1
      end if
      if ( ltype(ip) == itypeveg ) then
        map(ig,ptype(ip)+1) = ip
      end if
      ip = ip + 1
    end do ploop
  end subroutine makemap_pft

  subroutine makemap_col(nc,ng,mc,ctype,glat,glon,clat,clon,map)
    implicit none
    integer(ik4) , intent(in) :: nc , ng , mc
    integer(ik4) , intent(in) , dimension(nc) :: ctype
    real(rkx) , intent(in) , dimension(ng) :: glat
    real(rkx) , intent(in) , dimension(ng) :: glon
    real(rkx) , intent(in) , dimension(nc) :: clat
    real(rkx) , intent(in) , dimension(nc) :: clon
    integer(ik4) , intent(out) , dimension(ng,mc) :: map
    integer(ik4) :: ig , ic

    map(:,:) = -1
    ig = 1
    ic = 1
    cloop: &
    do while ( ic <= nc )
      if ( abs(glat(ig)-clat(ic)) > epsilon(1.0) .and. &
           abs(glon(ig)-clon(ic)) > epsilon(1.0) ) then
        ig = ig + 1
      end if
      map(ig,myfindloc(colval,ctype(ic))) = ic
      ic = ic + 1
    end do cloop
  end subroutine makemap_col

  subroutine pft_interpolate(vname)
    implicit none
    character(len=*) , intent(in) :: vname
    integer(ik4) :: ip , ig

    if ( .not. hasval(ncin,ncout,vname) ) return

    write(stdout,*) 'Interpolating ',trim(vname)

    call rval(ncin,vname,i_pval)
    call rval(ncout,vname,o_pval)

    do ip = 1 , imaxpft
      ! Extract values for pft and put on grid values
      do ig = 1 , ingc
        if ( i_mapf(ig,ip) > 0 ) then
          if ( i_pval(i_mapf(ig,ip)) >= -1.0e-20 .and. &
               i_pval(i_mapf(ig,ip)) <= 1.0e+20 ) then
            i_gval(ig) = i_pval(i_mapf(ig,ip))
          else
            i_gval(ig) = missl
          end if
        else
          i_gval(ig) = missl
        end if
      end do
      ! Interpolate on grid.
      call h_interpolate_cont(hint,i_gval,o_gval)
      ! Put grid values back on output pft values
      do ig = 1 , ongc
        if ( o_mapf(ig,ip) > 0 .and. &
             abs(o_gval(ig)-missl) > epsilon(1.0) .and. &
             o_pfts1d_wtxy(ip) > 0.0_rkx ) then
          o_pval(o_mapf(ig,ip)) = o_gval(ig)
        end if
      end do
    end do
    call wval(ncout,vname,o_pval)
  end subroutine pft_interpolate

  subroutine col_interpolate_grnd(vname)
    implicit none
    character(len=*) , intent(in) :: vname
    integer(ik4) :: ic , ig , igr

    if ( .not. hasval(ncin,ncout,vname) ) return

    write(stdout,*) 'Interpolating ',trim(vname)

    call rval2(ncin,vname,i_cval_ng)
    call rval2(ncout,vname,o_cval_ng)

    i_cval_ng_t = transpose(i_cval_ng)
    o_cval_ng_t = transpose(o_cval_ng)

    do igr = 1 , ingrd
      do ic = 1 , maxcol
        ! Extract values for pft and put on grid values
        do ig = 1 , ingc
          if ( i_mapc(ig,ic) > 0 ) then
            if ( .not. is_nan(i_cval_ng_t(i_mapc(ig,ic),igr)) ) then
              if ( i_cval_ng_t(i_mapc(ig,ic),igr) >= -1.0e-20 .and. &
                   i_cval_ng_t(i_mapc(ig,ic),igr) <= 1.0e+20 ) then
                i_gval(ig) = i_cval_ng_t(i_mapc(ig,ic),igr)
              else
                i_gval(ig) = missl
              end if
            else
              i_gval(ig) = missl
            end if
          else
            i_gval(ig) = missl
          end if
        end do
        o_gval(:) = missl
        ! Interpolate on grid.
        call h_interpolate_cont(hint,i_gval,o_gval)
        ! Put grid values back on output pft values
        do ig = 1 , ongc
          if ( o_mapc(ig,ic) > 0 .and. &
               abs(o_gval(ig)-missl) > epsilon(1.0) .and. &
               o_cols1d_wtxy(ic) > 0.0_rkx ) then
            o_cval_ng_t(o_mapc(ig,ic),igr) = o_gval(ig)
          end if
        end do
      end do
    end do
    o_cval_ng = transpose(o_cval_ng_t)
    call wval2(ncout,vname,o_cval_ng)
  end subroutine col_interpolate_grnd

  subroutine col_interpolate(vname)
    implicit none
    character(len=*) , intent(in) :: vname
    integer(ik4) :: ic , ig

    if ( .not. hasval(ncin,ncout,vname) ) return

    write(stdout,*) 'Interpolating ',trim(vname)

    call rval(ncin,vname,i_cval)
    call rval(ncout,vname,o_cval)

    do ic = 1 , maxcol
      ! Extract values for pft and put on grid values
      do ig = 1 , ingc
        if ( i_mapc(ig,ic) > 0 ) then
          if ( .not. is_nan(i_cval(i_mapc(ig,ic))) ) then
            if ( i_cval(i_mapc(ig,ic)) >= -1.0e-20 .and. &
                 i_cval(i_mapc(ig,ic)) <= 1.0e+20 ) then
              i_gval(ig) = i_cval(i_mapc(ig,ic))
            else
              i_gval(ig) = missl
            end if
          else
            i_gval(ig) = missl
          end if
        else
          i_gval(ig) = missl
        end if
      end do
      o_gval(:) = missl
      ! Interpolate on grid.
      call h_interpolate_cont(hint,i_gval,o_gval)
      ! Put grid values back on output pft values
      do ig = 1 , ongc
        if ( o_mapc(ig,ic) > 0 .and. &
             abs(o_gval(ig)-missl) > epsilon(1.0) .and. &
             o_cols1d_wtxy(ic) > 0.0_rkx ) then
          o_cval(o_mapc(ig,ic)) = o_gval(ig)
        end if
      end do
    end do
    call wval(ncout,vname,o_cval)
  end subroutine col_interpolate

  integer(ik4) function countcols(ctype) result(res)
    implicit none
    integer(ik4) , dimension(:) , intent(in) :: ctype
    integer(ik4) :: i , mcol
    res = 0
    mcol = 0
    do i = 1 , size(ctype)
      if ( ctype(i) > mcol ) then
        res = res + 1
        mcol = ctype(i)
      end if
    end do
  end function countcols

  subroutine set_colvals(c1,c2,cv)
    implicit none
    integer(ik4) , dimension(:) , intent(in) :: c1 , c2
    integer(ik4) , pointer , dimension(:) , intent(out) :: cv
    integer(ik4) :: i , mcol , ic , ic1 , ic2 , imax , imin
    integer(ik4) , dimension(:) , allocatable :: iv1 , iv2 , itemp , ival

    ic = 0
    ic1 = countcols(c1)
    allocate(iv1(ic1))
    mcol = 0
    do i = 1 , size(c1)
      if ( c1(i) > mcol ) then
        ic = ic + 1
        iv1(ic) = c1(i)
        mcol = c1(i)
      end if
    end do
    ic = 0
    ic2 = countcols(c2)
    allocate(iv2(ic2))
    mcol = 0
    do i = 1 , size(c2)
      if ( c2(i) > mcol ) then
        ic = ic + 1
        iv2(ic) = c2(i)
        mcol = c2(i)
      end if
    end do
    allocate(itemp(ic1+ic2))
    allocate(ival(ic1+ic2))
    itemp(1:ic1) = iv1
    itemp(ic1+1:ic1+ic2) = iv2
    ival(:) = 0
    ic = 0
    imin = minval(itemp)-1
    imax = maxval(itemp)
    do while ( imin < imax )
      ic = ic+1
      imin = minval(itemp, mask=itemp>imin)
      ival(ic) = imin
    end do
    call getmem1d(cv,1,ic,'set_colvals:cv')
    cv(:) = ival(1:ic)
    deallocate(itemp,ival,iv1,iv2)
  end subroutine set_colvals

  pure integer(ik4) function myfindloc(ac,c) result(ip)
    implicit none
    integer(ik4) , dimension(:) , intent(in) :: ac
    integer(ik4) , intent(in) :: c
    integer(ik4) :: i
    ip = 0
    do i = 1 , size(ac)
      if ( ac(i) == c ) then
        ip = i
        return
      end if
    end do
  end function myfindloc

#endif

end program interpinic
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

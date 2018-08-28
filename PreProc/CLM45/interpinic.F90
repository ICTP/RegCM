
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
  use mod_memutil
  use mod_message
  use mod_kdinterp
  use netcdf

  implicit none

  character (len=256) :: prgname
  character (len=256) :: inputfile
  character (len=256) :: outputfile
  character (len=64) :: sres

  integer(ik4) :: istat
  integer(ik4) :: ncin , ncout

  integer(ik4) :: hostnm
  integer(ik4) :: ihost , idir
  integer(ik4) :: getcwd
  integer(ik4) :: ingc , inlu , inco , inpft
  integer(ik4) :: ongc , onlu , onco , onpft

  integer(ik4) , parameter :: itypeveg = 1

  real(rk8) :: ires
  real(rk8) , pointer , dimension(:) :: i_grid1d_lon
  real(rk8) , pointer , dimension(:) :: i_grid1d_lat
  real(rk8) , pointer , dimension(:) :: i_pft1d_lon
  real(rk8) , pointer , dimension(:) :: i_pft1d_lat
  real(rk8) , pointer , dimension(:) :: o_grid1d_lon
  real(rk8) , pointer , dimension(:) :: o_grid1d_lat
  real(rk8) , pointer , dimension(:) :: o_pft1d_lon
  real(rk8) , pointer , dimension(:) :: o_pft1d_lat
  real(rk8) , pointer , dimension(:) :: i_gval
  real(rk8) , pointer , dimension(:) :: o_gval
  real(rk8) , pointer , dimension(:) :: i_pval
  real(rk8) , pointer , dimension(:) :: o_pval
  integer(ik4) , pointer , dimension(:) :: i_ltype
  integer(ik4) , pointer , dimension(:) :: o_ltype
  integer(ik4) , pointer , dimension(:) :: i_vtype
  integer(ik4) , pointer , dimension(:) :: o_vtype
  integer(ik4) , pointer , dimension(:,:) :: i_map
  integer(ik4) , pointer , dimension(:,:) :: o_map

  type(h_interpolator) :: hint
  integer(ik4) :: imaxpft , omaxpft
  integer(ik4) :: imaxlun , omaxlun

  integer(ik4) , dimension(8) :: tval
  character (len=32) :: cdata='?'
  character (len=5) :: czone='?'
  character (len=32) :: hostname='?'
  character (len=32) :: user='?'
  character (len=128) :: directory='?'
  character (len=*) , parameter :: f99001 = &
          '(2x," SVN Revision: ",a," compiled at: data : ",a,"  time: ",a,/)'

  write (stdout,  &
     "(/,2x,'This is interpinic part of RegCM package version 4')")
  write (stdout,f99001)  SVN_REV, __DATE__ , __TIME__

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
  call get_command_argument(1,value=inputfile,status=istat)
  if ( istat < 0 ) then
    write(stderr, *) 'Usage : ', trim(prgname) , &
              ' inputfile.r.nc outputfile.r.nc inreskm'
    call die(__FILE__,'Input file name missing',__LINE__)
  end if
  call get_command_argument(2,value=outputfile,status=istat)
  if ( istat < 0 ) then
    write(stderr, *) 'Usage : ', trim(prgname) , &
              ' inputfile.r.nc outputfile.r.nc inreskm'
    call die(__FILE__,'Output file name missing',__LINE__)
  end if
  call get_command_argument(3,value=sres,status=istat)
  if ( istat < 0 ) then
    write(stderr, *) 'Usage : ', trim(prgname) , &
              ' inputfile.r.nc outputfile.r.nc inreskm'
    call die(__FILE__,'inreskm missing',__LINE__)
  else
    read(sres,fmt=*,iostat=istat, err=100) ires
    if ( istat < 0 ) then
      write(stderr, *) 'Usage : ', trim(prgname) , &
                ' inputfile.r.nc outputfile.r.nc inreskm'
      call die(__FILE__,'inreskm not parsable',__LINE__)
    end if
  end if
  if ( ires < 0.0_rk8 ) then
    write(stderr, *) 'inreskm : ',ires,' ????'
    call die(__FILE__,'inreskm not parsable',__LINE__)
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

  write(stdout,*) 'Input DIMENSIONS'
  write(stdout,*) 'Input gridcell  : ',ingc
  write(stdout,*) 'Input landunit  : ',inlu
  write(stdout,*) 'Input column    : ',inco
  write(stdout,*) 'Input pft       : ',inpft
  write(stdout,*) ''

  ongc = dlen(ncout,'gridcell')
  onlu = dlen(ncout,'landunit')
  onco = dlen(ncout,'column')
  onpft = dlen(ncout,'pft')

  write(stdout,*) 'Output DIMENSIONS'
  write(stdout,*) 'Output gridcell : ',ongc
  write(stdout,*) 'Output landunit : ',onlu
  write(stdout,*) 'Output column   : ',onco
  write(stdout,*) 'Output pft      : ',onpft
  write(stdout,*) ''

  call getmem1d(i_grid1d_lon,1,ingc,'interpinic:i_grid1d_lon')
  call getmem1d(i_grid1d_lat,1,ingc,'interpinic:i_grid1d_lat')
  call getmem1d(i_gval,1,ingc,'interpinic:i_gval')
  call getmem1d(i_pft1d_lon,1,inpft,'interpinic:i_pft1d_lon')
  call getmem1d(i_pft1d_lat,1,inpft,'interpinic:i_pft1d_lat')
  call getmem1d(i_pval,1,inpft,'interpinic:i_pval')
  call getmem1d(i_ltype,1,inpft,'interpinic:i_ltype')
  call getmem1d(i_vtype,1,inpft,'interpinic:i_vtype')
  call getmem1d(o_grid1d_lon,1,ongc,'interpinic:o_grid1d_lon')
  call getmem1d(o_grid1d_lat,1,ongc,'interpinic:o_grid1d_lat')
  call getmem1d(o_gval,1,ongc,'interpinic:o_gval')
  call getmem1d(o_pft1d_lon,1,onpft,'interpinic:o_pft1d_lon')
  call getmem1d(o_pft1d_lat,1,onpft,'interpinic:o_pft1d_lat')
  call getmem1d(o_pval,1,onpft,'interpinic:o_pval')
  call getmem1d(o_ltype,1,onpft,'interpinic:o_ltype')
  call getmem1d(o_vtype,1,onpft,'interpinic:o_vtype')

  call rval(ncin,'grid1d_lon',i_grid1d_lon)
  call rval(ncin,'grid1d_lat',i_grid1d_lat)
  call rval(ncin,'pfts1d_lon',i_pft1d_lon)
  call rval(ncin,'pfts1d_lat',i_pft1d_lat)
  call ival(ncin,'pfts1d_ityplun',i_ltype)
  call ival(ncin,'pfts1d_itypveg',i_vtype)
  call rval(ncout,'grid1d_lon',o_grid1d_lon)
  call rval(ncout,'grid1d_lat',o_grid1d_lat)
  call rval(ncout,'pfts1d_lon',o_pft1d_lon)
  call rval(ncout,'pfts1d_lat',o_pft1d_lat)
  call ival(ncout,'pfts1d_ityplun',o_ltype)
  call ival(ncout,'pfts1d_itypveg',o_vtype)

  imaxpft = maxval(i_vtype) + 1 ! The bare ground class
  imaxlun = maxval(i_ltype)
  omaxpft = maxval(o_vtype) + 1 ! The bare ground class
  omaxlun = maxval(o_ltype)

  write(stdout,*) 'Input  PFT number : ', imaxpft
  write(stdout,*) 'Output PFT number : ', omaxpft
  write(stdout,*) 'Input  LUN number : ', imaxlun
  write(stdout,*) 'Output LUN number : ', omaxlun

  if ( imaxpft /= omaxpft ) then
    write(stderr,*) 'Number of PFT do not match !'
    call die(__FILE__,'CAN WORK ONLY IF THE PFT MATCH.',__LINE__)
  end if

  call getmem2d(i_map,1,ingc,1,imaxpft,'interpinic:i_map')
  call getmem2d(o_map,1,ongc,1,omaxpft,'interpinic:o_map')

  call h_interpolator_create(hint, &
                       ingc, i_grid1d_lat, i_grid1d_lon, &
                       ongc, o_grid1d_lat, o_grid1d_lon, &
                       ires)

  write(stdout,*) 'Input/output grids successfully read.'

  call makemap(inpft,ingc,imaxpft,i_ltype,i_vtype, &
               i_grid1d_lat,i_grid1d_lon,i_pft1d_lat,i_pft1d_lon,i_map)
  write(stdout,*) 'Input Map created.'

  call makemap(onpft,ongc,omaxpft,o_ltype,o_vtype, &
               o_grid1d_lat,o_grid1d_lon,o_pft1d_lat,o_pft1d_lon,o_map)
  write(stdout,*) 'Output Map created.'

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
  call pft_interpolate('lgsf')
  call pft_interpolate('bglfr')
  call pft_interpolate('bgtr')
  call pft_interpolate('gpp_pepv')
  call pft_interpolate('availc')
  call pft_interpolate('xsmrpool_recover')
  call pft_interpolate('alloc_pnow')
  call pft_interpolate('c_allometry')
  call pft_interpolate('n_allometry')
  call pft_interpolate('plant_ndemand')
  call pft_interpolate('tempsum_potential_gpp')
  call pft_interpolate('annsum_potential_gpp')
  call pft_interpolate('tempmax_retransn')
  call pft_interpolate('annmax_retransn')
  call pft_interpolate('avail_retransn')
  call pft_interpolate('plant_nalloc')
  call pft_interpolate('plant_calloc')
  call pft_interpolate('excess_cflux')
  call pft_interpolate('downreg')
  call pft_interpolate('prev_leafc_to_litter')
  call pft_interpolate('prev_frootc_to_litter')
  call pft_interpolate('tempsum_npp')
  call pft_interpolate('annsum_npp')
  call pft_interpolate('annsum_npp')
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
  call pft_interpolate('CROWNAREA')
  call pft_interpolate('tempsum_litfall')
  call pft_interpolate('annsum_litfall')
  !call pft_interpolate('nind')
  !call pft_interpolate('leafcmax')
  !call pft_interpolate('PREC365_VALUE')
  !call pft_interpolate('AGDDTW_VALUE')
  !call pft_interpolate('AGDD_VALUE')

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

  write(stdout,*) 'Successfully completed CLM interpolation.'

  call exit(0)

 100 call die(__FILE__,'inreskm not parsable',__LINE__)

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
    real(rk8) , dimension(:) , intent(inout) :: var
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
    real(rk8) , dimension(:) , intent(in) :: var
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

  subroutine makemap(np,ng,mp,ltype,ptype,glat,glon,plat,plon,map)
    implicit none
    integer(ik4) , intent(in) :: np , ng , mp
    integer(ik4) , intent(in) , dimension(np) :: ltype
    integer(ik4) , intent(in) , dimension(np) :: ptype
    real(rk8) , intent(in) , dimension(ng) :: glat
    real(rk8) , intent(in) , dimension(ng) :: glon
    real(rk8) , intent(in) , dimension(np) :: plat
    real(rk8) , intent(in) , dimension(np) :: plon
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
  end subroutine makemap

  subroutine pft_interpolate(vname)
    implicit none
    character(len=*) , intent(in) :: vname
    integer(ik4) :: ip , ig

    write(stdout,*) 'Interpolating ',trim(vname)

    if ( .not. hasval(ncin,ncout,vname) ) return

    call rval(ncin,vname,i_pval)
    call rval(ncout,vname,o_pval)

    do ip = 1 , imaxpft
      ! Extract values for pft and put on grid values
      do ig = 1 , ingc
        if ( i_map(ig,ip) > 0 ) then
          i_gval(ig) = i_pval(i_map(ig,ip))
        else
          i_gval(ig) = -9999.0_rk8
        end if
      end do
      ! Interpolate on grid.
      call h_interpolate_cont(hint,i_gval,o_gval)
      ! Put grid values back on output pft values
      do ig = 1 , ongc
        if ( o_map(ig,ip) > 0 .and. &
             abs(o_gval(ig)+9999.0_rk8) > epsilon(1.0) ) then
          o_pval(o_map(ig,ip)) = o_gval(ig)
        end if
      end do
    end do
    call wval(ncout,vname,o_pval)
  end subroutine pft_interpolate

#endif

end program interpinic
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_sst_gnmnc

  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_sst_grid
  use mod_kdinterp
  use mod_message
  use mod_nchelper
  use mod_hadgem_helper
  use mod_csiro_helper
  use mod_canesm_helper
  use mod_miroc_helper
  use mod_ipsl_helper
  use mod_gfdl_helper
  use mod_cnrm_helper
  use mod_ecearth_helper
  use mod_ccsm4_helper
  use netcdf

  private

  integer(ik4) :: ilon, jlat

  integer(ik4) :: inet1
  integer(ik4), dimension(2) :: ivar2
  integer(ik4) :: timlen
  integer(ik4) :: timid
  integer(ik4) :: istatus
  integer(ik4), dimension(3) :: istart, icount
  real(rkx), pointer, contiguous, dimension(:) ::  work1
  real(rkx), pointer, contiguous, dimension (:,:) :: work2, work3
  real(rkx), pointer, contiguous, dimension(:,:) :: sst
  type(rcm_time_and_date), save :: fidate1
  character(len=64) :: cunit, ccal
  character(len=256) :: inpfile
  character(len=8), dimension(2) :: varname

  data varname/'time', 'TOBESET'/

  real(rkx), pointer, contiguous, dimension(:) :: glat
  real(rkx), pointer, contiguous, dimension(:) :: glon
  real(rkx), pointer, contiguous, dimension(:,:) :: glat2
  real(rkx), pointer, contiguous, dimension(:,:) :: glon2

  type(h_interpolator) :: hint

  public :: sst_gnmnc

  contains
  !
  !****************************************************************************
  !
  ! This is a package of subroutines to read SST data on regular latlon grid
  ! in NETCDF format and interpolate on RCM grid.
  ! Need to configure datapath and variable names.
  !
  ! Written By Moetasim Ashfaq Dec-2005 @ PURDUE.EDU for POP dataset prepared
  ! on the regular 1x1 degree resolution using NCL+NCO programs.
  !
  ! Modified to be generic reader for monthly netCDF by Graziano Giuliani 2011
  !
  !****************************************************************************
  !
  subroutine sst_gnmnc
    implicit none
    type(rcm_time_and_date) :: idate, idatef, idateo, imm1
    integer(ik4) :: i, j, k, nsteps
    real(rkx) :: ufac

    imm1 = prevmon(globidate1)

    ufac = 0.0
    if ( ssttyp == 'CAM4N' ) then
      inpfile = trim(inpglob)// &
         '/SST/sst_HadOIBl_bc_0.9x1.25_1870_2008_c091020.nc'
      varname(2) = 'SST_cpl'
      ufac = 273.15
    else if ( ssttyp == 'JRA55' ) then
      inpfile = trim(inpglob)//'/SST/sstCOBE.monmean.nc'
      varname(2) = 'sst'
      ufac = 273.15
    else if ( ssttyp == 'CCSST' ) then
      inpfile = trim(inpglob)//'/SST/ccsm_mn.sst.nc'
      varname(2) = 'SST'
      ufac = 273.15
    else if ( ssttyp(1:3) == 'CA_' ) then
      call find_canesm_sst(inpfile,imm1)
      varname(2) = 'ts'
    else if ( ssttyp(1:3) == 'HA_' ) then
      call find_hadgem_sst(inpfile,imm1)
      varname(2) = 'tos'
    else if ( ssttyp(1:3) == 'CS_' ) then
      call find_csiro_sst(inpfile,imm1)
      varname(2) = 'tos'
    else if ( ssttyp(1:3) == 'MI_' ) then
      call find_miroc_sst(inpfile,imm1)
      varname(2) = 'tos'
    else if ( ssttyp(1:3) == 'EC_' ) then
      call find_ecearth_sst(inpfile,imm1,.false.)
      varname(2) = 'sst'
    else if ( ssttyp(1:3) == 'IP_' ) then
      call find_ipsl_sst(inpfile,imm1)
      varname(2) = 'tos'
    else if ( ssttyp(1:3) == 'GF_' ) then
      call find_gfdl_sst(inpfile,imm1)
      varname(2) = 'ts'
    else if ( ssttyp(1:3) == 'CN_' ) then
      call find_cnrm_sst(inpfile,imm1)
      varname(2) = 'tos'
    else if ( ssttyp(1:3) == 'CC_' ) then
      call find_ccsm4_sst(inpfile,imm1)
      varname(2) = 'tos'
    else
      call die('gnmnc_sst','Unknown ssttyp: '//ssttyp,1)
    end if

    call open_input(inpfile)

    call getmem2d(work2,1,ilon,1,jlat,'mod_gnmnc_sst:work2')
    call getmem2d(work3,1,ilon,1,jlat,'mod_gnmnc_sst:work3')
    call getmem2d(sst,1,ilon,1,jlat,'mod_gnmnc_sst:sst')

    idateo = monfirst(globidate1)
    idatef = monfirst(globidate2)
    if (idatef < globidate2) then
      idatef = nextmon(idatef)
    end if
    nsteps = imondiff(idatef,idateo) + 1

    call open_sstfile(idateo)

    idate = idateo
    do k = 1, nsteps
      call gnmnc_sst(idate)
      call h_interpolate_cont(hint,sst,sstmm)
      do i = 1, iy
        do j = 1, jx
          if ( sstmm(j,i) > -999.0 ) sstmm(j,i) = sstmm(j,i) + ufac
        end do
      end do

      call writerec(idate)

      write (stdout,*) 'WRITEN OUT SST DATA : ', tochar(idate)
      idate = nextmon(idate)

    end do

    call h_interpolator_destroy(hint)

  end subroutine sst_gnmnc
  !
  ! Subroutine to read required records from SST data file
  !
  subroutine gnmnc_sst(idate)
    implicit none
    type(rcm_time_and_date), intent (in) :: idate
    real(rkx) :: wt1, wt2
    type(rcm_time_and_date) :: prev, next

    integer(ik4) :: it, i, j
    type(rcm_time_interval) :: tdiff1, tdiff2

    icount(1) = ilon
    icount(2) = jlat
    icount(3) = 1
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1

    it = imondiff(idate,fidate1)

    call read_month(idate,it,work2)
    prev = timeval2date(work1(it),cunit,ccal)

    it = it + 1
    call read_month(idate,it,work3)
    next = timeval2date(work1(it),cunit,ccal)

    tdiff1 = next-idate
    tdiff2 = next-prev
    wt1 = real(tohours(tdiff1)/tohours(tdiff2))
    wt2 = 1.0 - wt1
    do j = 1, jlat
      do i = 1, ilon
        if (work2(i,j) < 0.9E+20 .and. work3(i,j) < 0.9E+20 ) then
          sst(i,j) = work2(i,j)*wt2+work3(i,j)*wt1
        else
          sst(i,j) = -9999.0
        end if
      end do
    end do
  end subroutine gnmnc_sst

  subroutine open_input(fname)
    implicit none
    character(len=*), intent(in) :: fname
    integer(ik4) :: istatus
    integer(ik4) :: latid, lonid, ndims
    logical, save :: firstpass = .true.

    istatus = nf90_open(fname,nf90_nowrite,inet1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error opening '//trim(fname))
    write (stdout,*) 'Opened ', trim(fname)

    if ( firstpass ) then
      istatus = nf90_inq_dimid(inet1,'lat',latid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(inet1,'j',latid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(inet1,'rlat',latid)
        end if
      end if
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find dim lat/j/rlat')
      istatus = nf90_inq_dimid(inet1,'lon',lonid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(inet1,'i',lonid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(inet1,'rlon',lonid)
        end if
      end if
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find dim lon/i/rlon')
      istatus = nf90_inquire_dimension(inet1,latid,len=jlat)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dim lat')
      istatus = nf90_inquire_dimension(inet1,lonid,len=ilon)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dim lon')

      istatus = nf90_inq_varid(inet1,'lat',latid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var lat')
      istatus = nf90_inq_varid(inet1,'lon',lonid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var lon')

      istatus = nf90_inquire_variable(inet1, latid, ndims=ndims)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'INQUIRE LAT ERROR')
      if ( ndims == 1 ) then
        call getmem1d(glat,1,jlat,'mod_gnmnc_sst:glat')
        call getmem1d(glon,1,ilon,'mod_gnmnc_sst:glon')
        istatus = nf90_get_var(inet1,latid,glat)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var lat')
        istatus = nf90_get_var(inet1,lonid,glon)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var lon')
        call h_interpolator_create(hint,glat,glon,xlat,xlon)
      else
        if ( ssttyp(1:3) == 'IP_' ) then
          ilon = ilon - 2
          jlat = jlat - 2
        end if
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        call getmem2d(glat2,1,ilon,1,jlat,'mod_gnmnc_sst:glat2')
        call getmem2d(glon2,1,ilon,1,jlat,'mod_gnmnc_sst:glon2')
        istatus = nf90_get_var(inet1,latid,glat2,istart(1:2),icount(1:2))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var lat')
        istatus = nf90_get_var(inet1,lonid,glon2,istart(1:2),icount(1:2))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var lon')
        call h_interpolator_create(hint,glat2,glon2,xlat,xlon)
      end if

      firstpass = .false.

    end if

    istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var '//varname(1))
    istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var '//varname(2))

    istatus = nf90_inq_dimid(inet1,'time',timid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim time')
    istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim time')
    call getmem1d(work1,1,timlen,'mod_gnmnc_sst:work1')

    istatus = nf90_get_var(inet1,ivar2(1),work1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(1))
    istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(1)//' units')
    if ( ssttyp == 'JRA55' ) then
      ccal = 'gregorian'
    else
      istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1)//' calendar')
    end if
    fidate1 = timeval2date(work1(1),cunit,ccal)

  end subroutine open_input

  subroutine read_month(idate,it,vv)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4), intent(inout) :: it
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: vv
    logical :: lswitch
    icount(1) = ilon
    icount(2) = jlat
    icount(3) = 1
    istart(1) = 1
    istart(2) = 1
    if ( it < 1 ) then
      call die('gnmnc_sst','Timestep not in file.',1)
    end if
    if ( it > timlen ) then
      lswitch = .false.
      if ( ssttyp(1:3) == 'CA_' ) then
        call find_canesm_sst(inpfile,idate)
        lswitch = .true.
      else if ( ssttyp(1:3) == 'HA_' ) then
        call find_hadgem_sst(inpfile,idate)
        lswitch = .true.
      else if ( ssttyp(1:3) == 'CS_' ) then
        call find_csiro_sst(inpfile,idate)
        lswitch = .true.
      else if ( ssttyp(1:3) == 'EC_' ) then
        call find_ecearth_sst(inpfile,idate,.false.)
        lswitch = .true.
      else if ( ssttyp(1:3) == 'MI_' ) then
        call find_miroc_sst(inpfile,idate)
        lswitch = .true.
      else if ( ssttyp(1:3) == 'IP_' ) then
        call find_ipsl_sst(inpfile,idate)
        lswitch = .true.
      else if ( ssttyp(1:3) == 'GF_' ) then
        call find_gfdl_sst(inpfile,idate)
        lswitch = .true.
      else if ( ssttyp(1:3) == 'CN_' ) then
        call find_cnrm_sst(inpfile,idate)
        lswitch = .true.
      else if ( ssttyp(1:3) == 'CC_' ) then
        call find_ccsm4_sst(inpfile,idate)
        lswitch = .true.
      end if
      if ( lswitch ) then
        write(stdout,*) 'Switching file...'
        istatus = nf90_close(inet1)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error Close file')
        call open_input(inpfile)
        it = imondiff(idate,fidate1) + 1
      end if
    end if
    istart(3) = it
    istatus = nf90_get_var(inet1,ivar2(2),vv,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(2))
    if ( ssttyp(1:3) == 'MI_' ) then
      where ( vv < 100.0_rkx )
        vv = 1e+20_rkx
      end where
    end if
  end subroutine read_month

end module mod_sst_gnmnc

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

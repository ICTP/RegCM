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

module mod_sst_gnhnc

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_sst_grid
  use mod_earth
  use mod_kdinterp
  use mod_message
  use mod_nchelper
  use mod_ccsm3_helper
  use mod_ecearth_helper
  use mod_mpiesm_helper
  use mod_lgm_helper
  use mod_date
  use netcdf

  private

  integer(ik4) :: ilon, jlat

  integer(ik4) :: inet1
  integer(ik4), dimension(2) :: ivar2
  integer(ik4) :: timlen
  integer(ik4) :: timid
  integer(ik4) :: istatus
  integer(ik4) :: itcfs
  integer(ik4) :: lyear
  integer(ik4), dimension(3) :: istart, icount
  integer(ik2), pointer, contiguous, dimension(:,:) :: work
  integer(ik2) :: fillvalue
  real(rkx), pointer, contiguous, dimension(:) ::  work1
  real(rkx), pointer, contiguous, dimension(:,:) :: work2
  real(rkx), pointer, contiguous, dimension(:,:) :: sst
  real(rkx) :: add_offset, scale_factor
  logical :: cds_beta = .false.
  type(rcm_time_and_date), save :: fidate1
  character(len=64) :: cunit, ccal
  character(len=256) :: inpfile
  character(len=8), dimension(2) :: varname
  type(global_domain) :: gdomain
  type(h_interpolator) :: hint

  data varname/'time', 'TOBESET'/

  public :: sst_gnhnc

  contains
  !
  !**************************************************************************
  !
  ! This is a package of subroutines to read SST data on regular latlon grid
  ! in NETCDF format and interpolate on RCM grid.
  ! Need to configure datapath and variable names.
  !
  !**************************************************************************
  subroutine sst_gnhnc
    implicit none
    real(rkx), pointer, contiguous, dimension(:) :: glat
    real(rkx), pointer, contiguous, dimension(:) :: glon
    real(rkx), pointer, contiguous, dimension(:) :: grev
    type(rcm_time_and_date) :: idate, idatef, idateo
    type(rcm_time_interval) :: tdif
    integer(ik4) :: k, nsteps, latid, lonid
    integer(ik4) :: year, month, day, hour

    call split_idate(globidate1, year, month, day, hour)

    if ( ssttyp(1:3) == 'MP_' ) then
      call find_mpiesm_sst(inpfile,globidate1,'M')
      varname(2) = 'tos'
    else if ( ssttyp(1:3) == 'MPL' ) then
      call find_mpiesm_sst(inpfile,globidate1,'L')
      varname(2) = 'tos'
    else if ( ssttyp(1:3) == 'ECC' ) then
      call find_ecearth_sst(inpfile,globidate1,.true.)
      varname(2) = 'sst'
    else if ( ssttyp == 'EIXXX' ) then
      write(inpfile,'(a)') trim(inpglob)//'/ERAIN_MEAN/SST/sst_xxxx_xxxx.nc'
      varname(2) = 'sst'
    else if ( ssttyp(1:3) == 'EIN' ) then
      write (inpfile,'(a,i0.4,a)') &
        trim(inpglob)//pthsep//ssttyp//pthsep//'SST'//pthsep//'sst.',year, '.nc'
      varname(2) = 'sst'
    else if ( ssttyp(1:4) == 'ERA5' ) then
      write (inpfile,'(a,i0.4,a,i0.2,a)') &
        trim(inpglob)//pthsep//ssttyp(1:4)//pthsep//'SST'//pthsep// &
         'sst_',year,'_',month,'.nc'
      varname(2) = 'sst'
    else if ( ssttyp(1:3) == 'LGM' ) then
      lyear = year
      if ( ssttyp(4:4) == 'P' ) then
        call find_lgm_sst(inpfile,globidate1,'P')
      else
        call find_lgm_sst(inpfile,globidate1,'C')
      end if
      varname(2) = 'tsurf'
    else if ( ssttyp(1:3) == 'CFS' ) then
      write(inpfile,'(a,i0.4,i0.2,i0.2,i0.2,a,i0.4,i0.2,i0.2,i0.2,a)') &
        trim(inpglob)//'/CFS/',year,month,day,hour, &
        '/'//ssttyp(4:5)//'/SST/sst.',year,month,day,hour,'.nc'
      varname(2) = 'sst'
      itcfs = 1
    else if ( ssttyp(1:2) == 'E5' ) then
      write(inpfile,'(a,i4,a,i4,a)') &
        trim(inpglob)//'/ECHAM5/SST/EH5_OM'//ssttyp(3:5)//'_1_TSW_', &
        year,'010100-',year+1,'010100.nc'
      varname(2) = 'tos'
    else if ( ssttyp == 'CCSM3' ) then
      call find_ccsm3_file(inpfile,year,month,day,hour)
      varname(2) = 'SST'
    else
      call die('gnhnc_sst','Unknown ssttyp: '//ssttyp,1)
    end if

    istatus = nf90_open(inpfile,nf90_nowrite,inet1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error opening '//trim(inpfile))
    write (stdout,*) 'Opened ', trim(inpfile)

    istatus = nf90_inq_dimid(inet1,'lat',latid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet1,'latitude',latid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find dim lat')
    end if
    istatus = nf90_inq_dimid(inet1,'lon',lonid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet1,'longitude',lonid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find dim lon')
    end if
    istatus = nf90_inquire_dimension(inet1,latid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim lat')
    istatus = nf90_inquire_dimension(inet1,lonid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim lon')

    istatus = nf90_inq_dimid(inet1,'time',timid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet1,'valid_time',timid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find dim time')
    end if
    istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim time')

    call getmem1d(work1,1,timlen,'mod_gnhnc_sst:work1')

    istatus = nf90_inq_varid(inet1,'lat',latid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(inet1,'latitude',latid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var lat')
    end if
    istatus = nf90_inq_varid(inet1,'lon',lonid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(inet1,'longitude',lonid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var lon')
    end if
    istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(inet1,'valid_'//varname(1),ivar2(1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//varname(1))
    end if
    istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var '//varname(2))

    call getmem1d(glat,1,jlat,'mod_gnhnc_sst:glat')
    call getmem1d(glon,1,ilon,'mod_gnhnc_sst:glon')
    call getmem1d(grev,1,max(ilon,jlat),'mod_gnhnc_sst:glon')

    istatus = nf90_get_var(inet1,latid,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lat')
    istatus = nf90_get_var(inet1,lonid,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lon')
    !
    ! Find window to read
    !
    call get_window(glat,glon,xlat,xlon,i_band,gdomain)
    grev(1:jlat) = glat
    jlat = gdomain%nj
    call getmem1d(glat,1,jlat,'mod_gnhnc_sst:glat')
    glat = grev(gdomain%jgstart:gdomain%jgstop)
    grev(1:ilon) = glon
    ilon = sum(gdomain%ni)
    call getmem1d(glon,1,ilon,'mod_gnhnc_sst:glon')
    glon(1:gdomain%ni(1)) = grev(gdomain%igstart(1):gdomain%igstop(1))
    if ( gdomain%ntiles == 2 ) then
      glon(gdomain%ni(1)+1:ilon) = grev(gdomain%igstart(2):gdomain%igstop(2))
    end if

    call h_interpolator_create(hint,glat,glon,xlat,xlon)

    call getmem2d(work2,1,ilon,1,jlat,'mod_gnhnc_sst:work2')
    call getmem2d(sst,1,ilon,1,jlat,'mod_gnhnc_sst:sst')

    if (ssttyp(1:3) == 'LGM' ) then
      ccal = 'gregorian'
    else if ( ssttyp /= 'EIXXX' .and. ssttyp(1:3) /= 'CFS' .and. &
              ssttyp(1:3) /= 'EIN' .and. ssttyp(1:4) /= 'ERA5' ) then
      istart(1) = 1
      icount(1) = timlen
      istatus = nf90_get_var(inet1,ivar2(1),work1,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1))

      istatus = nf90_myget_att_text(inet1,ivar2(1),'units',cunit)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read attribute units of '//varname(1))
      istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
      if ( istatus /= nf90_noerr ) ccal = 'gregorian'
      fidate1 = timeval2date(work1(1),cunit,ccal)
    else
      istatus = nf90_get_att(inet1,ivar2(2),'add_offset',add_offset)
      if ( istatus /= nf90_noerr ) then
        add_offset = 0.0_rkx
        scale_factor = 1.0_rkx
        cds_beta = .true.
      else
        istatus = nf90_get_att(inet1,ivar2(2),'_FillValue',fillvalue)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(2)//' _FillValue')
        istatus = nf90_get_att(inet1,ivar2(2),'scale_factor',scale_factor)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(2)//' scale_factor')
        write(stdout,*) 'Add offset   = ',add_offset
        write(stdout,*) 'Scale factor = ',scale_factor
        write(stdout,*) 'Fill Value   = ',fillvalue
        call getmem2d(work,1,ilon,1,jlat,'mod_gnhnc_sst:work')
      end if
    end if

    if ( ssttyp(1:3) == 'CFS' .or. &
         ssttyp(1:3) == 'EIN' .or. &
         ssttyp(1:4) == 'ERA5' .or. &
         ssttyp(1:3) == 'LGM' ) then
      idateo = globidate1
    else
      idateo = monfirst(globidate1)
    end if
    idatef = globidate2
    tdif = idatef-idateo
    nsteps = int(tohours(tdif))/ibdyfrq + 1

    write (stdout,*) 'GLOBIDATE1 : ', tochar(globidate1)
    write (stdout,*) 'GLOBIDATE2 : ', tochar(globidate2)
    write (stdout,*) 'NSTEPS = ', nsteps

    call open_sstfile(idateo)

    idate = idateo
    tdif = ibdyfrq*3600
    do k = 1, nsteps
      call gnhnc_sst(idate)
      call h_interpolate_cont(hint,sst,sstmm)
      call writerec(idate)
      write (stdout,*) 'WRITEN OUT SST DATA : ', tochar(idate)
      idate = idate + tdif
    end do

    call h_interpolator_destroy(hint)

  end subroutine sst_gnhnc
  !
  !     Subroutine to read required records from SST data file
  !
  subroutine gnhnc_sst(idate)
    implicit none
    type(rcm_time_and_date), intent (in) :: idate
    integer(ik4) :: it, i, j
    integer(ik4) :: year, month, day, hour
    type(rcm_time_interval) :: tdif
    integer(ik4), dimension(12) :: isteps

    data isteps /1,125,237,361,481,605,725,849,973,1093,1217,1337/

    istart(3) = 1

    if ( ssttyp == 'EIXXX' ) then
      call split_idate(idate, year, month, day, hour)
      it = isteps(month) + (day-1)*4 + hour/ibdyfrq
    else if ( ssttyp(1:3) == 'LGM' ) then
      call split_idate(idate, year, month, day, hour)
      if ( lyear /= year ) then
        lyear = year
        it = timlen + 1
      else
        it = (dayofyear(idate)-1)*4 + hour/ibdyfrq + 1
      end if
    else if ( ssttyp(1:3) == 'CFS' ) then
      it = itcfs
      itcfs = itcfs + 1
    else
      tdif = idate-fidate1
      it = int(tohours(tdif))/ibdyfrq + 1
    end if

    if ( it > timlen ) then
      ! Try switching to next file
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error Close file')
      call split_idate(idate, year, month, day, hour)
      if ( ssttyp(1:3) == 'MP_' ) then
        call find_mpiesm_sst(inpfile,idate,'M')
      else if ( ssttyp(1:3) == 'MPL' ) then
        call find_mpiesm_sst(inpfile,idate,'L')
      else if ( ssttyp(1:3) == 'ECC' ) then
        call find_ecearth_sst(inpfile,idate,.true.)
      else if ( ssttyp == 'CCSM3' ) then
        call find_ccsm3_file(inpfile,year,month,day,hour)
      else if ( ssttyp(1:3) == 'LGM' ) then
        if ( ssttyp(4:4) == 'P' ) then
          call find_lgm_sst(inpfile,idate,'P')
        else
          call find_lgm_sst(inpfile,idate,'C')
        end if
      else if ( ssttyp(1:3) == 'EIN' ) then
        write (inpfile,'(a,i0.4,a)') &
          trim(inpglob)//pthsep//ssttyp//pthsep//'SST'// &
          pthsep//'sst.',year, '.nc'
      else if ( ssttyp(1:4) == 'ERA5' ) then
        write (inpfile,'(a,i0.4,a,i0.2,a)') &
          trim(inpglob)//pthsep//ssttyp(1:4)//pthsep//'SST'//pthsep// &
           'sst_',year,'_',month,'.nc'
      end if
      istatus = nf90_open(inpfile,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error opening '//trim(inpfile))
      write (stdout,*) 'Opened ', trim(inpfile)
      istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(inet1,'valid_'//varname(1),ivar2(1))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var '//varname(1))
      end if
      istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//varname(2))
      if ( ssttyp(1:3) == 'CFS' .or. &
           ssttyp(1:3) == 'EIN' .or. &
           ssttyp(1:4) == 'ERA5' ) then
        istatus = nf90_get_att(inet1,ivar2(2),'add_offset',add_offset)
        if ( istatus /= nf90_noerr ) then
          add_offset = 0.0_rkx
          scale_factor = 1.0_rkx
          cds_beta = .true.
        else
          istatus = nf90_get_att(inet1,ivar2(2),'_FillValue',fillvalue)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(2)//' _FillValue')
          istatus = nf90_get_att(inet1,ivar2(2),'scale_factor',scale_factor)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(2)//' scale_factor')
          write(stdout,*) 'Add offset   = ',add_offset
          write(stdout,*) 'Scale factor = ',scale_factor
          write(stdout,*) 'Fill Value   = ',fillvalue
        end if
      end if
      if ( ssttyp(1:3) == 'LGM' ) then
        istatus = nf90_inq_dimid(inet1,'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim time')
        istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error inquire dim time')
        call split_idate(idate, year, month, day, hour)
        it = (dayofyear(idate)-1)*4 + hour/ibdyfrq + 1
      else
        istatus = nf90_inq_dimid(inet1,'time',timid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(inet1,'valid_time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find dim time')
        end if
        istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error inquire dim time')
        call getmem1d(work1,1,timlen,'mod_gnhnc_sst:work1')
        istart(1) = 1
        icount(1) = timlen
        istatus = nf90_get_var(inet1,ivar2(1),work1,istart(1:1),icount(1:1))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(1))
        istatus = nf90_myget_att_text(inet1,ivar2(1),'units',cunit)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read attribute units of '//varname(1))
        istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
        if ( istatus /= nf90_noerr ) ccal = 'gregorian'
        fidate1 = timeval2date(work1(1),cunit,ccal)
        tdif = idate-fidate1
        it = int(tohours(tdif))/ibdyfrq + 1
      end if
    end if
    icount(3) = 1
    istart(3) = it
    if ( ssttyp == 'EIXXX' .or. ssttyp(1:3) == 'CFS' .or. &
         ssttyp(1:3) == 'EIN' .or. ssttyp(1:4) == 'ERA5' ) then
      if ( cds_beta ) then
        call getworkf(2,work2)
        where ( is_nan(work2) )
          work2 = 1E+20_rkx
        end where
      else
        call getworki(2,work)
        work2 = 1E+20
        where ( work /= fillvalue )
          work2 = work * scale_factor + add_offset
        end where
      end if
    else
      call getworkf(2,work2)
      if ( ssttyp(1:2) == 'E5' .or. ssttyp(1:2) == 'MP' ) then
        where ( abs(work2-273.15) < 0.001 )
          work2 = 1E+20
        end where
      else if ( ssttyp == 'CCSM3' ) then
        where ( work2 < 1.0 )
          work2 = 1E+20
        end where
      end if
    end if
    do j = 1, jlat
      do i = 1, ilon
        if (work2(i,j) < 0.9E+20 ) then
          sst(i,j) = work2(i,j)
        else
          sst(i,j) = -9999.0
        end if
      end do
    end do

    contains

      subroutine getworkf(irec,wk)
        implicit none
        integer(ik4), intent(in) :: irec
        real(rkx), pointer, contiguous, dimension(:,:) :: wk
        integer(ik4) :: itile, iti, itf
        iti = 1
        do itile = 1, gdomain%ntiles
          istart(1) = gdomain%igstart(itile)
          icount(1) = gdomain%ni(itile)
          istart(2) = gdomain%jgstart
          icount(2) = gdomain%nj
          itf = iti + gdomain%ni(itile) - 1
          istatus = nf90_get_var(inet1,ivar2(irec),wk(iti:itf,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine getworkf

      subroutine getworki(irec,wk)
        implicit none
        integer(ik4), intent(in) :: irec
        integer(ik2), pointer, contiguous, dimension(:,:) :: wk
        integer(ik4) :: itile, iti, itf
        iti = 1
        do itile = 1, gdomain%ntiles
          istart(1) = gdomain%igstart(itile)
          icount(1) = gdomain%ni(itile)
          istart(2) = gdomain%jgstart
          icount(2) = gdomain%nj
          itf = iti + gdomain%ni(itile) - 1
          istatus = nf90_get_var(inet1,ivar2(irec),wk(iti:itf,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine getworki

  end subroutine gnhnc_sst

end module mod_sst_gnhnc

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

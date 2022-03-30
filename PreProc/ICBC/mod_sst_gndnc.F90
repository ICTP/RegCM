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

module mod_sst_gndnc

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_sst_grid
  use mod_earth
  use mod_date
  use mod_kdinterp
  use mod_message
  use mod_nchelper
  use mod_noresm_helper
  use netcdf

  private

  integer(ik4) :: ilon , jlat

  integer(ik4) :: inet1
  integer(ik4) , dimension(2) :: ivar2
  integer(ik4) :: timlen
  integer(ik4) :: timid
  integer(ik4) :: istatus
  integer(ik4) , dimension(3) :: istart , icount
  real(rkx) , pointer ::  work(:)
  real(rkx) , pointer , dimension(:,:) :: workf
  integer(ik2) , pointer , dimension(:,:) :: worki
  real(rkx) , pointer , dimension(:,:) :: sst
  type(rcm_time_and_date) , save :: fidate1
  character(len=64) :: cunit , ccal
  character(len=256) :: inpfile
  character(len=8), dimension(2) :: varname
  type(h_interpolator) :: hint
  integer(ik2) :: fillvalue
  real(rkx) :: add_offset , scale_factor

  data varname/'time', 'TOBESET'/

  public :: sst_gndnc

  contains
  !
  !**************************************************************************
  !
  ! This is a package of subroutines to read SST data on regular latlon grid
  ! in NETCDF format and interpolate on RCM grid.
  ! Need to configure datapath and variable names.
  !
  !**************************************************************************
  subroutine sst_gndnc
    implicit none
    real(rkx) , pointer , dimension(:,:) :: glat
    real(rkx) , pointer , dimension(:,:) :: glon
    real(rkx) , pointer , dimension(:) :: glat1
    real(rkx) , pointer , dimension(:) :: glon1
    type(rcm_time_and_date) :: idate , idatef , idateo
    type(rcm_time_interval) :: tdif
    integer(ik4) :: k , nsteps , latid , lonid
    integer(ik4) :: year , month , day , hour

    call split_idate(globidate1, year, month, day, hour)

    if ( ssttyp(1:3) == 'NO_' ) then
      call find_noresm_sst(inpfile,globidate1)
      varname(2) = 'tos'
    else if ( ssttyp == 'TMIST' ) then
      write (inpfile,'(a,i0.4,a)') &
         trim(inpglob)//pthsep//'SST'//pthsep//'TMI'//pthsep// &
            'tmisst', year, '.nc'
      varname(1) = 'time'
      varname(2) = 'sst'
    else if ( ssttyp(1:3) == 'EID' ) then
      write (inpfile,'(a,i0.4,a)') trim(inpglob)//pthsep// &
        dattyp//pthsep//'SSTD'//pthsep//'sst.',year, '.nc'
      varname(1) = 'time'
      varname(2) = 'sst'
    else if ( ssttyp(1:4) == 'ERA5' ) then
      write (inpfile,'(a,i0.4,a,i0.2,a)') &
        trim(inpglob)//pthsep//ssttyp(1:4)//pthsep//'SSTD'//pthsep// &
         'sst_',year,'_',month,'.nc'
      varname(1) = 'time'
      varname(2) = 'sst'
    else
      call die('gndnc_sst','Unknown ssttyp: '//ssttyp,1)
    end if

    istatus = nf90_open(inpfile,nf90_nowrite,inet1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error opening '//trim(inpfile))
    write (stdout,*) 'Opened ', trim(inpfile)

    istatus = nf90_inq_dimid(inet1,'lat',latid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet1,'latitude',latid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(inet1,'j',latid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(inet1,'LAT',latid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim lat')
        end if
      end if
    end if
    istatus = nf90_inq_dimid(inet1,'lon',lonid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet1,'longitude',lonid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(inet1,'i',lonid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(inet1,'LON',lonid)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_dimid(inet1,'LONN719_720',lonid)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error find dim lon')
          end if
        end if
      end if
    end if
    istatus = nf90_inquire_dimension(inet1,latid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim lat')
    istatus = nf90_inquire_dimension(inet1,lonid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim lon')

    istatus = nf90_inq_dimid(inet1,'time',timid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet1,'TIME',timid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find dim time')
    end if
    istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim time')

    call getmem1d(work,1,timlen,'mod_gndnc_sst:work')

    istatus = nf90_inq_varid(inet1,'lat',latid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(inet1,'latitude',latid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(inet1,'LAT',latid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var lat')
      end if
    end if
    istatus = nf90_inq_varid(inet1,'lon',lonid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(inet1,'longitude',lonid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(inet1,'LON',lonid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_varid(inet1,'LONN719_720',lonid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var lon')
        end if
      end if
    end if
    istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var '//varname(1))
    istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var '//varname(2))

    if ( ssttyp(1:3) == 'NO_' ) then
      call getmem2d(glat,1,ilon,1,jlat,'mod_gndnc_sst:glat')
      call getmem2d(glon,1,ilon,1,jlat,'mod_gndnc_sst:glon')
      istatus = nf90_get_var(inet1,latid,glat)
      call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lat')
      istatus = nf90_get_var(inet1,lonid,glon)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var lon')
      call h_interpolator_create(hint,glat,glon,xlat,xlon)
    else if ( ssttyp(1:4) == 'ERA5' .or. ssttyp(1:3) == 'EID' ) then
      istatus = nf90_get_att(inet1,ivar2(2),'_FillValue',fillvalue)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(2)//' _FillValue')
      istatus = nf90_get_att(inet1,ivar2(2),'add_offset',add_offset)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(2)//' add_offset')
      istatus = nf90_get_att(inet1,ivar2(2),'scale_factor',scale_factor)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(2)//' scale_factor')
      write(stdout,*) 'Add offset   = ',add_offset
      write(stdout,*) 'Scale factor = ',scale_factor
      write(stdout,*) 'Fill Value   = ',fillvalue
      call getmem2d(worki,1,ilon,1,jlat,'mod_gndnc_sst:worki')
      call getmem1d(glat1,1,jlat,'mod_gndnc_sst:glat')
      call getmem1d(glon1,1,ilon,'mod_gndnc_sst:glon')
      istatus = nf90_get_var(inet1,latid,glat1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lat')
      istatus = nf90_get_var(inet1,lonid,glon1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var lon')
      call h_interpolator_create(hint,glat1,glon1,xlat,xlon)
    else
      call getmem1d(glat1,1,jlat,'mod_gndnc_sst:glat')
      call getmem1d(glon1,1,ilon,'mod_gndnc_sst:glon')
      istatus = nf90_get_var(inet1,latid,glat1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lat')
      istatus = nf90_get_var(inet1,lonid,glon1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var lon')
      call h_interpolator_create(hint,glat1,glon1,xlat,xlon)
    end if

    call getmem2d(workf,1,ilon,1,jlat,'mod_gndnc_sst:workf')
    call getmem2d(sst,1,ilon,1,jlat,'mod_gndnc_sst:sst')

    istart(1) = 1
    icount(1) = timlen
    istatus = nf90_get_var(inet1,ivar2(1),work,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(1))

    istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(1)//' units')
    istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
    if ( istatus /= nf90_noerr ) then
      ccal = 'gregorian'
    end if
    fidate1 = timeval2date(work(1),cunit,ccal)
    idateo = globidate1
    idatef = globidate2
    tdif = idatef-idateo
    nsteps = int(tohours(tdif))/24 + 1

    write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
    write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
    write (stdout,*) 'NSTEPS = ', nsteps

    call open_sstfile(idateo)

    idate = idateo
    tdif = 86400
    do k = 1 , nsteps
      call gndnc_sst(idate)
      call h_interpolate_cont(hint,sst,sstmm)
      call writerec(idate)
      write (stdout,*) 'WRITEN OUT SST DATA : ' , tochar(idate)
      idate = idate + tdif
    end do

    call h_interpolator_destroy(hint)

  end subroutine sst_gndnc
  !
  !     Subroutine to read required records from SST data file
  !
  subroutine gndnc_sst(idate)
    implicit none
    type(rcm_time_and_date) , intent (in) :: idate
    integer(ik4) :: it , i , j
    integer(ik4) :: year , month , day , hour
    type(rcm_time_interval) :: tdif

    istart(3) = 1

    tdif = idate-fidate1
    it = int(tohours(tdif))/24 + 1

    if ( it > timlen ) then
      ! Try switching to next file
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error Close file')
      call split_idate(idate, year, month, day, hour)
      if ( ssttyp(1:3) == 'NO_' ) then
        call find_noresm_sst(inpfile,idate)
      else if ( ssttyp(1:3) == 'EID' ) then
        write (inpfile,'(a,i0.4,a)') &
          trim(inpglob)//pthsep//dattyp//pthsep//'SSTD'// &
          pthsep//'sst.',year, '.nc'
      else if ( ssttyp(1:4) == 'ERA5' ) then
        write (inpfile,'(a,i0.4,a,i0.2,a)') &
          trim(inpglob)//pthsep//ssttyp(1:4)//pthsep//'SSTD'//pthsep// &
           'sst_',year,'_',month,'.nc'
      else if ( ssttyp == 'TMIST' ) then
        write (inpfile,'(a,i0.4,a)') &
           trim(inpglob)//pthsep//'SST'//pthsep//'TMI'//pthsep// &
              'tmisst', year, '.nc'
      end if
      istatus = nf90_open(inpfile,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error opening '//trim(inpfile))
      write (stdout,*) 'Opened ', trim(inpfile)
      istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//varname(1))
      istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//varname(2))
      istatus = nf90_inq_dimid(inet1,'time',timid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(inet1,'TIME',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim time')
      end if
      istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dim time')
      call getmem1d(work,1,timlen,'mod_gndnc_sst:work')
      if ( ssttyp(1:4) == 'ERA5' .or. ssttyp(1:3) == 'EID' ) then
        istatus = nf90_get_att(inet1,ivar2(2),'_FillValue',fillvalue)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(2)//' _FillValue')
        istatus = nf90_get_att(inet1,ivar2(2),'add_offset',add_offset)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(2)//' add_offset')
        istatus = nf90_get_att(inet1,ivar2(2),'scale_factor',scale_factor)
        call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(2)//' scale_factor')
        write(stdout,*) 'Add offset   = ',add_offset
        write(stdout,*) 'Scale factor = ',scale_factor
        write(stdout,*) 'Fill Value   = ',fillvalue
      end if
      istart(1) = 1
      icount(1) = timlen
      istatus = nf90_get_var(inet1,ivar2(1),work,istart(1:1),icount(1:1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1))
      istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1)//' units')
      istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
      if ( istatus /= nf90_noerr ) ccal = 'gregorian'
      fidate1 = timeval2date(work(1),cunit,ccal)
      tdif = idate-fidate1
      it = int(tohours(tdif))/24 + 1
    end if
    istart(1) = 1
    icount(1) = ilon
    istart(2) = 1
    icount(2) = jlat
    icount(3) = 1
    istart(3) = it
    if ( ssttyp(1:4) == 'ERA5' .or. ssttyp(1:3) == 'EID' ) then
      call getworki(2,workf,worki)
    else
      call getworkf(2,workf)
    end if
    if ( ssttyp == 'TMIST' ) then
      do j = 1 , jlat
        do i = 1 , ilon
          if (workf(i,j) > 0.0_rkx ) then
            sst(i,j) = workf(i,j) + 273.15_rkx
          else
            sst(i,j) = -9999.0_rkx
          end if
        end do
      end do
    else
      do j = 1 , jlat
        do i = 1 , ilon
          if (workf(i,j) < 0.9E+20_rkx ) then
            sst(i,j) = workf(i,j)
          else
            sst(i,j) = -9999.0_rkx
          end if
        end do
      end do
    end if

    contains

      subroutine getworki(irec,wk,wi)
        implicit none
        integer(ik4) , intent(in) :: irec
        real(rkx) , pointer , dimension(:,:) :: wk
        integer(ik2) , pointer , dimension(:,:) :: wi
        istatus = nf90_get_var(inet1,ivar2(irec),wi(:,:),istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error read var '//varname(irec))
        wk = 1E+20
        where ( wi /= fillvalue )
          wk = wi * scale_factor + add_offset
        end where
      end subroutine getworki

      subroutine getworkf(irec,wk)
        implicit none
        integer(ik4) , intent(in) :: irec
        real(rkx) , pointer , dimension(:,:) :: wk
        istatus = nf90_get_var(inet1,ivar2(irec),wk(:,:),istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error read var '//varname(irec))
      end subroutine getworkf

  end subroutine gndnc_sst

end module mod_sst_gndnc

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

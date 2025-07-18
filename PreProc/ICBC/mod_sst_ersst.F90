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

module mod_sst_ersst

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_date
  use mod_dynparam
  use mod_sst_grid
  use mod_kdinterp
  use mod_message
  use mod_memutil
  use mod_nchelper
  use netcdf

  private

  public :: sst_ersst

  integer(ik4) :: ilon, jlat
  real(rkx), pointer, contiguous, dimension(:,:) :: sst
  real(rkx), pointer, contiguous, dimension(:) :: lati
  real(rkx), pointer, contiguous, dimension(:) :: loni
  integer(2), pointer, contiguous, dimension(:,:,:) :: work
  real(rkx), dimension(2) :: xadd, xscale, xmiss

  integer(ik4) :: inet

  character(len=256) :: inpfile

  type(h_interpolator) :: hint

  contains
  !
  ! Comments on dataset sources and location:
  !
  ! ERSST    Extended Reconstructed Sea Surface Temperature (ERSST)
  !          Monthly from January 1854 to the present
  !
  subroutine sst_ersst
    implicit none
    integer(ik4) :: it
    integer(ik4) :: istatus
    integer(ik4) :: year, month, day, hour
    integer(ik4) :: dimi, vari
    integer(ik4) :: nsteps
    type(rcm_time_and_date) :: idate, idatef, idateo
    logical :: lfirst

    data lfirst /.true./

    idateo = monfirst(globidate1)
    idatef = monfirst(globidate2)
    if (idatef < globidate2) then
      idatef = nextmon(idatef)
    end if
    nsteps = imondiff(idatef,idateo) + 1

    write (stdout,*) 'GLOBIDATE1 : ', tochar(globidate1)
    write (stdout,*) 'GLOBIDATE2 : ', tochar(globidate2)
    write (stdout,*) 'NSTEPS     : ', nsteps

    ! SET UP LONGITUDES AND LATITUDES FOR SST DATA

    idate = prevmon(globidate1)
    call split_idate(idate,year,month,day,hour)
    write (inpfile,'(a,i0.4,i0.2,a)') &
      trim(inpglob)//pthsep//'SST'//pthsep//ssttyp// &
#ifdef USE_ERSST_V6
      pthsep//'V6'//pthsep//'ersst.v6.',year, month, '.nc'
#else
      pthsep//'ersst.v5.',year, month, '.nc'
#endif
    istatus = nf90_open(inpfile,nf90_nowrite,inet)
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Cannot open file '//trim(inpfile))
    lfirst = .true.
    write (stdout,*) 'Opened ', trim(inpfile)
    istatus = nf90_inq_dimid(inet,'lon',dimi)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lon dim')
    istatus = nf90_inquire_dimension(inet,dimi,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lon dim')
    istatus = nf90_inq_dimid(inet,'lat',dimi)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lat dim')
    istatus = nf90_inquire_dimension(inet,dimi,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lat dim')

    call getmem1d(loni,1,ilon,'sst_ersst:loni')
    call getmem1d(lati,1,jlat,'sst_ersst:lati')
    call getmem2d(sst,1,ilon,1,jlat,'sst_ersst:sst')
    call getmem3d(work,1,ilon,1,jlat,1,2,'sst_ersst:work')

    istatus = nf90_inq_varid(inet,'lon',vari)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lon var')
    istatus = nf90_get_var(inet,vari,loni)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read lon var')
    istatus = nf90_inq_varid(inet,'lat',vari)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lat var')
    istatus = nf90_get_var(inet,vari,lati)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read lat var')

    call h_interpolator_create(hint,lati,loni,xlat,xlon)

    call open_sstfile(idateo)

    idate = idateo
    do it = 1, nsteps
      call sst_readersst(idate,lfirst)
      call h_interpolate_cont(hint,sst,sstmm)
      call writerec(idate)
      write(stdout,*) 'WRITING OUT SST DATA:', tochar(idate)
      idate = nextmon(idate)
    end do

    call h_interpolator_destroy(hint)

  end subroutine sst_ersst

  subroutine sst_readersst(idate,lfirst)
    use netcdf
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    logical, intent(inout) :: lfirst
    integer(ik4) :: i, j
    integer(ik4) :: year, month, day, hour
    integer(ik4) :: istatus, ivar
    !
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers. The array 'sst'
    ! will contain the unpacked data.
    !
    if ( lfirst ) then
      istatus = nf90_inq_varid(inet,'sst',ivar)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot find variable sst')
      istatus = nf90_get_att(inet,ivar,'scale_factor',xscale(1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute scale_factor')
      istatus = nf90_get_att(inet,ivar,'add_offset',xadd(1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute add_offset')
      istatus = nf90_get_att(inet,ivar,'_FillValue',xmiss(1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute _FillValue')
      istatus = nf90_get_var(inet,ivar,work(:,:,1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                     'Cannot read sst')
      istatus = nf90_close(inet)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot close file')
      lfirst = .false.
    else
      xscale(1) = xscale(2)
      xadd(1) = xadd(2)
      xmiss(1) = xmiss(2)
      work(:,:,1) = work(:,:,2)
    end if

    call split_idate(idate,year,month,day,hour)
    write (inpfile,'(a,i0.4,i0.2,a)') &
      trim(inpglob)//pthsep//'SST'//pthsep//ssttyp// &
#ifdef USE_ERSST_V6
      pthsep//'V6'//pthsep//'ersst.v6.',year, month, '.nc'
#else
      pthsep//'ersst.v5.',year, month, '.nc'
#endif
    istatus = nf90_open(inpfile,nf90_nowrite,inet)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Cannot open file '//trim(inpfile))
    write (stdout,*) 'Opened ', trim(inpfile)
    istatus = nf90_inq_varid(inet,'sst',ivar)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Cannot find variable sst')
    istatus = nf90_get_att(inet,ivar,'scale_factor',xscale(2))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Cannot get attribute scale_factor')
    istatus = nf90_get_att(inet,ivar,'add_offset',xadd(2))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Cannot get attribute add_offset')
    istatus = nf90_get_att(inet,ivar,'_FillValue',xmiss(2))
    call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute _FillValue')
    istatus = nf90_get_var(inet,ivar,work(:,:,2))
    call checkncerr(istatus,__FILE__,__LINE__, &
                   'Cannot read sst')
    istatus = nf90_close(inet)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Cannot close file')
    do j = 1, jlat
      do i = 1, ilon
        if ( work(i,j,1) /= xmiss(1) .and. work(i,j,2) /= xmiss(2) ) then
          sst(i,j) = 0.5_rkx * &
                ( real(work(i,j,1),rkx)*xscale(1) + xadd(1) + 273.15_rkx + &
                  real(work(i,j,2),rkx)*xscale(2) + xadd(2) + 273.15_rkx )
          ! Respect convention for ice
          if ( sst(i,j) < 271.34_rkx ) sst(i,j) = 271.46_rkx
        else
          sst(i,j) = -9999.0_rkx
        end if
      end do
    end do
  end subroutine sst_readersst

end module mod_sst_ersst
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

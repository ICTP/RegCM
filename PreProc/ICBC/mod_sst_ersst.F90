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

module mod_sst_ersst

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_sst_grid
  use mod_interp
  use mod_message
  use mod_memutil
  use mod_nchelper
  use netcdf

  private

  public :: sst_ersst

  integer(ik4) :: ilon , jlat
  integer(ik4) , parameter :: idtbc = 6
  real(rk8) , pointer , dimension(:,:) :: sst
  real(rk8) , pointer , dimension(:) :: lati
  real(rk8) , pointer , dimension(:) :: loni
  integer(2) , pointer , dimension(:,:) :: work

  integer(ik4) :: inet
  contains
  !
  ! Comments on dataset sources and location:
  !
  ! ERAIN    ERAIN_SST is provided by ERA-Interim project
  !          6 hourly frequncy, 1.5x1.5 degree resolution
  !          from 1989010100 to present.
  !          'ERSST' for using the sea surface temperature;
  !          'ERSKT' for using the skin temperature.
  !
  ! ML = 1 is   0.0; ML = 2 is   1.5; => ML = 240 is 358.5E
  ! NL = 1 is  90.0; ML = 2 is  88.5; => ML = 121 is -90.
  !
  subroutine sst_ersst
    implicit none
    integer(ik4) :: it
    integer(ik4) :: istatus
    integer(ik4) :: year , month , day , hour , isyear = -1
    integer(ik4) :: dimi , vari
    integer(ik4) :: ierrec , nsteps
    type(rcm_time_and_date) :: idate , ierastart
    type(rcm_time_interval) :: tdiff , itbc
    character(len=256) :: inpfile
    logical :: lfirst

    data lfirst /.true./

    itbc = rcm_time_interval(idtbc,uhrs)
    tdiff = globidate2-globidate1
    nsteps = idnint(tohours(tdiff))/idtbc + 1

    write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
    write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
    write (stdout,*) 'NSTEPS     : ' , nsteps

    call open_sstfile(globidate1)

    ! SET UP LONGITUDES AND LATITUDES FOR SST DATA

    call split_idate(globidate1,year,month,day,hour)
    write (inpfile,'(a,i0.4,a)') &
      trim(inpglob)//pthsep//ssttyp//pthsep//'SST'//pthsep//'sst.',year, '.nc'
    istatus = nf90_open(inpfile,nf90_nowrite,inet)
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Cannot open file '//trim(inpfile))
    isyear = year
    ierastart = year*1000000 + 10100
    lfirst = .true.
    write (stdout,*) 'Opened ', trim(inpfile)
    istatus = nf90_inq_dimid(inet,'longitude',dimi)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lon dim')
    istatus = nf90_inquire_dimension(inet,dimi,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lon dim')
    istatus = nf90_inq_dimid(inet,'latitude',dimi)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lat dim')
    istatus = nf90_inquire_dimension(inet,dimi,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lat dim')

    call getmem1d(loni,1,ilon,'sst_ersst:loni')
    call getmem1d(lati,1,jlat,'sst_ersst:lati')
    call getmem2d(sst,1,ilon,1,jlat,'sst_ersst:sst')
    call getmem2d(work,1,ilon,1,jlat,'sst_ersst:work')

    istatus = nf90_inq_varid(inet,'longitude',vari)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lon var')
    istatus = nf90_get_var(inet,vari,loni)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read lon var')
    istatus = nf90_inq_varid(inet,'latitude',vari)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lat var')
    istatus = nf90_get_var(inet,vari,lati)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read lat var')

    idate = globidate1
    do it = 1 , nsteps
      call split_idate(idate,year,month,day,hour)
      if ( year /= isyear ) then
        istatus = nf90_close(inet)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Cannot close file')
        write (inpfile,'(a,i0.4,a)') &
               trim(inpglob)//pthsep//ssttyp//pthsep// &
               'SST'//pthsep//'sst.',year, '.nc'
        istatus = nf90_open(inpfile,nf90_nowrite,inet)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Cannot open file '//trim(inpfile))
        write (stdout,*) 'Opened ', trim(inpfile)
        isyear = year
        ierastart = year*1000000 + 10100
        lfirst = .true.
      end if

      tdiff = idate-ierastart
      ierrec = idnint(tohours(tdiff))/idtbc+1

      if ( ssttyp == 'ERSKT' ) then
        call sst_erain(ierrec,lfirst,2)
      else
        call sst_erain(ierrec,lfirst,1)
      end if

      call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,jx,iy,1)
      call writerec(idate)
      write(stdout,*) 'WRITING OUT SST DATA:' , tochar(idate)
      idate = idate + itbc
    end do
  end subroutine sst_ersst

  subroutine sst_erain(it,lfirst,itype)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: it , itype
    logical , intent(inout) :: lfirst
    integer(ik4) :: i , j
    character(len=4) , dimension(2) :: varname
    integer(ik4) :: istatus
    integer(ik4) , save :: ivar
    real(rk8) , save :: xadd , xscale , xmiss
    integer(ik4) , dimension(3) , save :: icount , istart
    !
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers. The array 'sst'
    ! will contain the unpacked data.
    !
    data varname/'sst','skt'/
    if ( lfirst ) then
      istatus = nf90_inq_varid(inet,varname(itype),ivar)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot find variable '//trim(varname(itype)))
      istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute scale_factor')
      istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute add_offset')
      istatus = nf90_get_att(inet,ivar,'_FillValue',xmiss)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute _FillValue')
      istart(1) = 1
      istart(2) = 1
      icount(1) = ilon
      icount(2) = jlat
      lfirst = .false.
    end if

    istart(3) = it
    icount(3) = 1
    istatus = nf90_get_var(inet,ivar,work,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                   'Cannot read '//trim(varname(itype)))

    do j = 1 , jlat
      do i = 1 , ilon
        if (work(i,j) /= xmiss) then
          sst(i,j) = dble(work(i,j))*xscale + xadd
          ! Respect convention for ice
          if ( sst(i,j) < 271.5D0 ) sst(i,j) = 271.0D0
        else
          sst(i,j) = -9999.0D0
        end if
        ! Respect convention on seaice
      end do
    end do
  end subroutine sst_erain

end module mod_sst_ersst
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

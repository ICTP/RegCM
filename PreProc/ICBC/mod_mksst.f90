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

module mod_mksst

  use m_realkinds
  use m_die
  use m_stdio
  use m_zeit
  use netcdf
  use mod_memutil
  use mod_constants
  use mod_dynparam

  private

  logical :: lopen , lhasice
  integer :: ncst , ntime
  integer , dimension(3) :: ivar
  type(rcm_time_and_date) , dimension(:) , pointer :: itime
  real(sp) , dimension(:,:) , pointer :: xlandu
  real(sp) , dimension(:,:) , pointer :: work1 , work2
  real(sp) , dimension(:,:) , pointer :: work3 , work4
  real(dp) , dimension(:) , pointer :: xtime

  data lopen/.false./
  character(256) :: sstfile

  public :: readsst , closesst

  contains
!
!-----------------------------------------------------------------------
!
  subroutine readsst(tsccm, idate)
    implicit none
    real(sp) , dimension(jx,iy) , intent(inout) :: tsccm
    type(rcm_time_and_date) , intent(in) :: idate
    integer :: istatus , idimid , itvar
    integer , dimension(3) :: istart , icount
    character(64) :: timeunits , timecal
    integer :: i , j , irec
    type(rcm_time_interval) :: ks1 , ks2
    real(sp) :: wt

    call zeit_ci('readsst')
    if (.not. lopen) then
      sstfile=trim(dirglob)//pthsep//trim(domname)//'_SST.nc'
      istatus = nf90_open(sstfile, nf90_nowrite, ncst)
      call check_ok(istatus,'Error Opening SST file'//trim(sstfile))
      istatus = nf90_inq_dimid(ncst, 'time', idimid)
      call check_ok(istatus,'Error time dimension SST file'//trim(sstfile))
      istatus = nf90_inquire_dimension(ncst, idimid, len=ntime)
      call check_ok(istatus,'Error time dimension SST file'//trim(sstfile))
      istatus = nf90_inq_varid(ncst, "time", itvar)
      call check_ok(istatus,'Error time variable SST file'//trim(sstfile))
      istatus = nf90_inq_varid(ncst, "landuse", ivar(1))
      call check_ok(istatus,'Error landuse variable SST file'//trim(sstfile))
      istatus = nf90_inq_varid(ncst, "sst", ivar(2))
      call check_ok(istatus,'Error sst variable SST file'//trim(sstfile))
      lhasice = .true.
      istatus = nf90_inq_varid(ncst, "ice", ivar(3))
      if ( istatus /= nf90_noerr) then
        lhasice = .false.
      end if
      istatus = nf90_get_att(ncst, itvar, "units", timeunits)
      call check_ok(istatus,'Error time var units SST file'//trim(sstfile))
      istatus = nf90_get_att(ncst, itvar, "calendar", timecal)
      call check_ok(istatus,'Error time var units SST file'//trim(sstfile))

      call getmem2d(xlandu,1,jx,1,iy,'mod_mksst:xlandu')
      call getmem2d(work1,1,jx,1,iy,'mod_mksst:work1')
      call getmem2d(work2,1,jx,1,iy,'mod_mksst:work2')
      if (lhasice) then
        call getmem2d(work3,1,jx,1,iy,'mod_mksst:work3')
        call getmem2d(work4,1,jx,1,iy,'mod_mksst:work4')
      end if
      call getmem1d(xtime,1,ntime,'mod_mksst:xtime')
      call getmem1d(itime,1,ntime,'mod_mksst:itime')

      istatus = nf90_get_var(ncst, itvar, xtime)
      call check_ok(istatus,'Error time var read SST file'//trim(sstfile))
      do i = 1 , ntime
        itime(i) = timeval2date(xtime(i), timeunits, timecal)
      end do
      lopen = .true.
    end if

    if (idate > itime(ntime) .or. idate < itime(1)) then
      write (stderr,*) 'Cannot find ', idate%tostring(), ' in SST file'
      write (stderr,*) 'Range is : ', itime(1)%tostring() , ' to ', &
                       itime(ntime)%tostring()
      call die('readsst')
    end if

    irec = 1
    do i = 1 , ntime
      if ( idate < itime(i) ) then
        exit
      end if
      irec = irec + 1
    end do
    irec = irec - 1

    istart(3) = irec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy
    icount(1) = jx
    istatus = nf90_get_var(ncst, ivar(1), xlandu, istart, icount)
    call check_ok(istatus,'Error landuse var read SST file'//trim(sstfile))
    istatus = nf90_get_var(ncst, ivar(2), work1, istart, icount)
    call check_ok(istatus,'Error sst var read SST file'//trim(sstfile))
    if (lhasice) then
      istatus = nf90_get_var(ncst, ivar(3), work3, istart, icount)
      call check_ok(istatus,'Error ice var read SST file'//trim(sstfile))
    end if
    if (idate == itime(irec)) then
      do i = 1 , jx
        do j = 1 , iy
          if ( (xlandu(i,j) > 13.9 .and. xlandu(i,j) < 15.1) .and. &
               (work1(i,j) > -900.0) ) then
            tsccm(i,j) = work1(i,j)
            if (lhasice) then
              if ( work3(i,j) > -900.0) then
                if ( work3(i,j) > 35. ) tsccm(i,j) = 271.0
              end if
            end if
          end if
        end do
      end do
    else
      istart(3) = irec+1
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = iy
      icount(1) = jx
      istatus = nf90_get_var(ncst, ivar(2), work2, istart, icount)
      call check_ok(istatus,'Error sst var read SST file'//trim(sstfile))
      if (lhasice) then
        istatus = nf90_get_var(ncst, ivar(3), work4, istart, icount)
        call check_ok(istatus,'Error ice var read SST file'//trim(sstfile))
      end if
      ks1 = itime(irec+1)-idate
      ks2 = itime(irec+1)-itime(irec)
      wt = real(ks1%hours()/ks2%hours())
      do i = 1 , jx
        do j = 1 , iy
          if ( (xlandu(i,j) > 13.9 .and. xlandu(i,j) < 15.1) .and.  &
               (work1(i,j) > -900.0 .and. work2(i,j) > -900.0) ) then
            tsccm(i,j) = (1.0-wt)*work1(i,j) + wt*work2(i,j)
            if (lhasice) then
              if ( work3(i,j) > -900.0 .and. work4(i,j) > -900.0 ) then
                if ( (1.-wt)*work3(i,j)+wt*work4(i,j) > 35. ) then
                  tsccm(i,j) = 271.0
                endif
              end if
            end if
          end if
        end do
      end do
    end if
    call zeit_co('readsst')
  end subroutine readsst

  subroutine closesst
    implicit none
    integer :: istatus
    istatus = nf90_close(ncst)
  end subroutine closesst
!
  subroutine check_ok(ierr,message)
    implicit none
    integer , intent(in) :: ierr
    character(*) :: message
    if (ierr /= nf90_noerr) then
      call die('mod_mksst',message,1,nf90_strerror(ierr),ierr)
    end if
  end subroutine check_ok
!
end module mod_mksst

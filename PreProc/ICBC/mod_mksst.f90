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

  use mod_constants
  use m_realkinds
  use m_die
  use m_stdio
  use m_zeit
  use m_mall

  private

  logical , private :: lopen , lhasice
  integer , private :: ncst , ntime
  integer , dimension(3) , private :: ivar
  integer , dimension(:) , allocatable , private :: itime
  real(sp) , dimension(:,:) , allocatable , private :: xlandu
  real(sp) , dimension(:,:) , allocatable , private :: work1 , work2
  real(sp) , dimension(:,:) , allocatable , private :: work3 , work4

  data lopen/.false./

  public :: readsst , closesst

  contains
!
!-----------------------------------------------------------------------
!
  subroutine readsst(tsccm, idate)
    use netcdf
    use mod_dynparam        
    use mod_date
    implicit none
    real(sp) , dimension(jx,iy) , intent(inout) :: tsccm
    integer , intent(in) :: idate
    real(dp) , dimension(:) , allocatable :: xtime
    integer :: istatus , idimid , itvar
    integer , dimension(3) :: istart , icount
    character(256) :: sstfile
    character(64) :: timeunits
    integer :: i , j , irec , ks1 , ks2
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

      allocate(xlandu(jx,iy), stat=istatus)
      if (istatus /= 0) call die('readsst','allocate xlandu',istatus)
      call mall_mci(xlandu,'readsst')
      allocate(work1(jx,iy), stat=istatus)
      if (istatus /= 0) call die('readsst','allocate work1',istatus)
      call mall_mci(work1,'readsst')
      allocate(work2(jx,iy), stat=istatus)
      if (istatus /= 0) call die('readsst','allocate work2',istatus)
      call mall_mci(work2,'readsst')
      if (lhasice) then
        allocate(work3(jx,iy), stat=istatus)
        if (istatus /= 0) call die('readsst','allocate work3',istatus)
        call mall_mci(work3,'readsst')
        allocate(work4(jx,iy), stat=istatus)
        if (istatus /= 0) call die('readsst','allocate work4',istatus)
        call mall_mci(work4,'readsst')
      end if
      allocate(xtime(ntime), stat=istatus)
      if (istatus /= 0) call die('readsst','allocate xtime',istatus)
      call mall_mci(xtime,'readsst')
      allocate(itime(ntime), stat=istatus)
      if (istatus /= 0) call die('readsst','allocate itime',istatus)
      call mall_mci(itime,'readsst')

      istatus = nf90_get_var(ncst, itvar, xtime)
      call check_ok(istatus,'Error time var read SST file'//trim(sstfile))
      do i = 1 , ntime
        itime(i) = timeval2idate(xtime(i), timeunits)
      end do
      deallocate(xtime)
      call mall_mco(xtime,'readsst')
      lopen = .true.
    end if

    if (idate > itime(ntime) .or. idate < itime(1)) then
      write (stderr,*) 'Cannot find ', idate, ' in SST file'
      write (stderr,*) 'Range is : ', itime(1) , '-', itime(ntime)
      call die('readsst')
    end if

    irec = 0
    do i = 1 , ntime
      if (idate <= itime(i)) then
        irec = i
        exit
      end if
    end do

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
      istart(3) = irec-1
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
      ks1 = idatediff(itime(irec),idate)
      ks2 = idatediff(itime(irec),itime(irec-1))
      wt = float(ks1)/float(ks2)
      do i = 1 , jx
        do j = 1 , iy
          if ( (xlandu(i,j) > 13.9 .and. xlandu(i,j) < 15.1) .and.  &
               (work1(i,j) > -900.0 .and. work2(i,j) > -900.0) ) then
            tsccm(i,j) = (1.-wt)*work1(i,j) + wt*work2(i,j)
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
    use netcdf
    implicit none
    integer :: istatus
    istatus = nf90_close(ncst)
    if (allocated(itime)) then
      call mall_mco(itime,'readsst')
      deallocate(itime)
    end if
    if (allocated(work1)) then
      call mall_mco(work1,'readsst')
      deallocate(work1)
    end if
    if (allocated(work2)) then
      call mall_mco(work2,'readsst')
      deallocate(work2)
    end if
    if (allocated(work3)) then
      call mall_mco(work3,'readsst')
      deallocate(work3)
    end if
    if (allocated(work4))then
      call mall_mco(work4,'readsst')
      deallocate(work4)
    end if
    if (allocated(xlandu))then
      call mall_mco(xlandu,'readsst')
      deallocate(xlandu)
    end if
  end subroutine closesst
!
  subroutine check_ok(ierr,message)
    use netcdf
    implicit none
    integer , intent(in) :: ierr
    character(*) :: message
    if (ierr /= nf90_noerr) then
      call die('mod_mksst',message,1,nf90_strerror(ierr),ierr)
    end if
  end subroutine check_ok
!
end module mod_mksst

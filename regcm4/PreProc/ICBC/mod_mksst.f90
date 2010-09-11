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

      logical , private :: lopen , lhasice
      integer , private :: ncid , ntime
      integer , dimension(3) , private :: ivar
      integer , dimension(:) , allocatable , private :: itime
      real(4) , dimension(:,:) , allocatable , private :: xlandu
      real(4) , dimension(:,:) , allocatable , private :: work1 , work2
      real(4) , dimension(:,:) , allocatable , private :: work3 , work4

      data lopen/.false./

      contains
!
!-----------------------------------------------------------------------
!
      subroutine readsst(tsccm, topogm, idate)
        use netcdf
        use mod_dynparam        
        use mod_date
        implicit none
        real(4) , dimension(jx,iy) , intent(in) :: topogm
        real(4) , dimension(jx,iy) , intent(inout) :: tsccm
        integer , intent(in) :: idate
        real(8) , dimension(:) , allocatable :: xtime
        integer :: istatus , idimid , itvar
        integer , dimension(3) :: istart , icount
        character(256) :: sstfile
        character(64) :: timeunits
        integer :: i , j , irec , ks1 , ks2
        real(4) :: wt
        if (.not. lopen) then
          sstfile = trim(dirglob)//pthsep//trim(domname)//'_SST.nc'
          istatus = nf90_open(sstfile, nf90_nowrite, ncid)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error Opening SST file ', trim(sstfile)
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_inq_dimid(ncid, 'time', idimid)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error time dimension not present'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_inquire_dimension(ncid, idimid, len=ntime)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error time dimension not present'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_inq_varid(ncid, "time", itvar)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error variable time'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_inq_varid(ncid, "landuse", ivar(1))
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error variable landuse'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_inq_varid(ncid, "sst", ivar(2))
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error variable landuse'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_inq_varid(ncid, "ice", ivar(3))
          if ( istatus /= nf90_noerr) then
            lhasice = .false.
          end if
          istatus = nf90_get_att(ncid, itvar, "units", timeunits)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error variable time attribute units'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          allocate(xlandu(jx,iy))
          allocate(work1(jx,iy))
          allocate(work2(jx,iy))
          if (lhasice) then
            allocate(work3(jx,iy))
            allocate(work4(jx,iy))
          end if
          allocate(xtime(ntime))
          allocate(itime(ntime))
          istatus = nf90_get_var(ncid, itvar, xtime)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error reading time variable'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          do i = 1 , ntime
            itime(i) = timeval2idate(xtime(i), timeunits)
          end do
          deallocate(xtime)
          lopen = .true.
        end if
          print *, 'Range is : ', itime(1) , ' - ', itime(ntime)

        if (idate > itime(ntime) .or. idate < itime(1)) then
          print *, 'Cannot find ', idate, ' in SST file'
          stop
          print *, 'Range is : ', itime(1) , '-', itime(ntime)
        end if

        irec = 0
        do i = 1 , ntime
          if (idate <= itime(i)) then
            irec = i
            exit
           end if
        end do
!!!! why this below ? 

!        if (idate == itime(irec-1)) then
!          irec = irec-1
!        end if

        istart(3) = irec
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = iy
        icount(1) = jx
        istatus = nf90_get_var(ncid, ivar(1), xlandu, istart, icount)
        if ( istatus /= nf90_noerr) then
          write (6,*) 'Error reading landuse variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(ncid, ivar(2), work1, istart, icount)
        if ( istatus /= nf90_noerr) then
          write (6,*) 'Error reading sst variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lhasice) then
          istatus = nf90_get_var(ncid, ivar(3), work3, istart, icount)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error reading ice variable'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        if (idate == itime(irec)) then
          do i = 1 , jx
            do j = 1 , iy
              if ( (topogm(i,j)<=1.) .and.                              &
             &     (xlandu(i,j)>13.9 .and. xlandu(i,j)<15.1) .and.      &
             &     (work1(i,j)>-900.0) ) tsccm(i,j) = work1(i,j)
            end do
          end do
          if (lhasice) then
            if ( work3(i,j)>-900.0) then
              if ( work3(i,j)>35. ) tsccm(i,j) = 273.15 - 2.15
            end if
          end if
        else
          istart(3) = irec-1
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = iy
          icount(1) = jx
          istatus = nf90_get_var(ncid, ivar(2), work2, istart, icount)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error reading sst variable'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          if (lhasice) then
            istatus = nf90_get_var(ncid, ivar(3), work4, istart, icount)
            if ( istatus /= nf90_noerr) then
              write (6,*) 'Error reading ice variable'
              write (6,*) nf90_strerror(istatus)
              stop
            end if
          end if
          ks1 = idatediff(itime(irec),idate)
          ks2 = idatediff(itime(irec),itime(irec-1))
          wt = float(ks1)/float(ks2)
          do i = 1 , jx
            do j = 1 , iy
              if ( (topogm(i,j)<=1.) .and.                              &
                 & (xlandu(i,j)>13.9 .and. xlandu(i,j)<15.1) .and.      &
                 & (work1(i,j)>-900.0 .and. work2(i,j)>-900.0) )        &
                tsccm(i,j) = (1.-wt)*work1(i,j) + wt*work2(i,j)
            end do
          end do
          if (lhasice) then
            if ( work3(i,j)>-900.0 .and. work4(i,j)>-900.0 ) then
              if ( (1.-wt)*work3(i,j)+wt*work4(i,j)>35. )               &
                tsccm(i,j) = 273.15 - 2.15
            end if
          end if
        end if
      end subroutine readsst

      subroutine closesst
        use netcdf
        implicit none
        integer :: istatus
        istatus = nf90_close(ncid)
        if (allocated(itime)) deallocate(itime)
        if (allocated(work1)) deallocate(work1)
        if (allocated(work2)) deallocate(work2)
        if (allocated(work3)) deallocate(work3)
        if (allocated(work4)) deallocate(work4)
        if (allocated(xlandu)) deallocate(xlandu)
      end subroutine closesst

      end module mod_mksst

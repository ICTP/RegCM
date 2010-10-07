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

      module mod_rdldtr

      use mod_block

      contains
!
!
! Read a netcdf file with a global dataset on a regular latitude/longitude grid
!
!  cfile = name of the input file
!  cvar = name of the var in the file
!  iires = file resolution in arc seconds
!  iores = output wanted resolution in minutes
!  lreg = flag for grid vertex or center registering
!  values = allocated space by the sub containing data: the caller is
!           in charge of the deallocate
!
      subroutine read_ncglob(cfile,cvar,iires,iores,lreg,values)
        use netcdf
        use mod_block
        implicit none
        character(len=*) , intent(in) :: cfile , cvar
        integer , intent(in) :: iires, iores
        logical , intent(in) :: lreg
        real(4) , dimension(:,:) , allocatable , intent(out) :: values
        integer :: ncid , ivar , istatus
        integer :: nlat , nlon , itl
        integer :: i , j , iosec , inpsec , iopsec , ifrac , ireg
        integer , dimension(2) :: istart, icount
        real(4) , dimension(:,:) , allocatable :: readbuf
        integer :: nlogb , nlagb , hnlogb , hnlagb
        integer , parameter :: secpd = 3600
        integer , parameter :: secpm = 60
        real(4) :: delta

        iosec = iores*secpm
        inpsec = secpd/iires
        nlogb = 360*inpsec
        nlagb = 180*inpsec
        if (mod(nlogb,iosec) /= 0) then
          print *, 'Sorry, subroutine read_ncglob does not support'// &
                   ' iores = ',iores
          stop
        end if
        hnlogb = nlogb/2
        hnlagb = nlagb/2
        iopsec = secpd/iosec
        ifrac = inpsec/iopsec

        ireg = 0
        if (lreg) ireg = 1
        delta = ((dble(iires)/2.0D0)/dble(secpd))*dble(ireg)

        grdltmn = floor(xminlat)  -delta
        grdltma = ceiling(xmaxlat)+delta
        if (grdltmn < -90+delta) grdltmn = -90.0+delta
        if (grdltma > +90-delta) grdltma =  90.0-delta
        if (lonwrap) then
          grdlnmn = -180.0+delta
          grdlnma =  180.0-delta
        else
          grdlnmn = floor(xminlon)  -delta
          grdlnma = ceiling(xmaxlon)+delta
        end if

        nlat = (grdltma-grdltmn)*inpsec
        if (lonwrap) then
          nlon = nlogb-ireg
        else if (lcrosstime) then
          nlon = ((180.0-delta-grdlnmn)+(180.0-delta+grdlnma))*inpsec
        else
          nlon = (grdlnma-grdlnmn)*inpsec
        end if

        nlatin = (nlat/ifrac)+1
        nlonin = (nlon/ifrac)
        nlat = nlat+1
        if (.not. lonwrap) then
          nlon = nlon+1
          nlonin = nlonin + 1
        end if

        allocate(values(nlonin,nlatin),stat=istatus)
        if (istatus /= 0) then
          print *, 'Memory error on allocating ', &
                    nlatin*nlonin*4,' bytes.'
          stop
        end if

        allocate(readbuf(nlon,nlat), stat=istatus)
        if (istatus /= 0) then
          print *, 'Memory error on allocating ',nlat*nlon*4,' bytes.'
          stop
        end if

        istatus = nf90_open(cfile, nf90_nowrite, ncid)
        call checkerr(istatus)
        istatus = nf90_inq_varid(ncid, cvar, ivar)
        call checkerr(istatus)

        istart(2) = hnlagb+(grdltmn*inpsec)+1
        if (lonwrap) then
          istart(1) = 1
        else
          istart(1) = hnlogb+(grdlnmn*inpsec)+1
        end if
        if (.not. lcrosstime) then
          ! Simple case: not crossing timeline
          icount(2) = nlat
          icount(1) = nlon
          istatus = nf90_get_var(ncid, ivar, readbuf, istart, icount)
          call checkerr(istatus)
        else
          ! Crossing timeline
          itl = nlogb - istart(1)
          icount(2) = nlat
          icount(1) = itl
          istatus = nf90_get_var(ncid, ivar, readbuf(1:itl,:), &
                                 istart, icount)
          call checkerr(istatus)
          istart(1) = 1
          icount(1) = nlon-itl
          istatus = nf90_get_var(ncid, ivar, readbuf(itl+1:nlon,:), &
                                 istart, icount)
          call checkerr(istatus)
        end if

        istatus = nf90_close(ncid)
        call checkerr(istatus)

        do i = 1 , nlatin
          do j = 1 , nlonin
            values(j,i) = readbuf((j-1)*ifrac+1,(i-1)*ifrac+1)
          end do
        end do
        deallocate(readbuf)

      end subroutine read_ncglob
!
      subroutine checkerr(ierr)
        use netcdf
        implicit none
        integer , intent(in) :: ierr
        if ( ierr /= nf90_noerr ) then
          write (6, *) nf90_strerror(ierr)
          stop
        end if
      end subroutine checkerr
!
      end module mod_rdldtr

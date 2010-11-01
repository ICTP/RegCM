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
!  imeth = flag for selecting resampling method
!     default = nearest
!          1  = mean
!          2  = median
!          3  = most present
!  values = allocated space by the sub containing data: the caller is
!           in charge of the deallocate
!
      subroutine read_ncglob(cfile,cvar,iires,iores,lreg,imeth,values)
        use netcdf
        use mod_block
        implicit none
        character(len=*) , intent(in) :: cfile , cvar
        integer , intent(in) :: iires , iores , imeth
        logical , intent(in) :: lreg
        real(4) , dimension(:,:) , allocatable , intent(out) :: values
        integer :: ncid , ivar , istatus
        integer :: nlat , nlon , itl
        integer :: i , j , iosec , inpsec , iopsec , ifrac , ireg
        integer , dimension(2) :: istart, icount
        real(4) , dimension(:,:) , allocatable :: readbuf
        real(4) , dimension(:) , allocatable :: copybuf
        integer :: nlogb , nlagb , hnlogb , hnlagb , nfrac
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

        print *, 'Opening '//trim(cfile)
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

        print '(a$)', ' Resampling'
        nfrac = ifrac*ifrac
        allocate(copybuf(nfrac), stat=istatus)
        if (istatus /= 0) then
          print *, 'Memory error on allocating ',nfrac*4,' bytes'
          stop
        end if
        select case (imeth)
          case (1)
            do i = 1 , nlatin
              if (mod(i,10) == 0) print '(a$)', '.'
              do j = 1 , nlonin
                call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                             (i-1)*ifrac+1,ifrac,lcrosstime)
                values(j,i) = sum(copybuf)/size(copybuf)
              end do
            end do
          case (2)
            do i = 1 , nlatin
              if (mod(i,10) == 0) print '(a$)', '.'
              do j = 1 , nlonin
                call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                             (i-1)*ifrac+1,ifrac,lcrosstime)
                call qsort(copybuf)
                values(j,i) = copybuf(nfrac/2)
              end do
            end do
          case (3)
            do i = 1 , nlatin
              if (mod(i,10) == 0) print '(a$)', '.'
              do j = 1 , nlonin
                call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                             (i-1)*ifrac+1,ifrac,lcrosstime)
                values(j,i) = mpindex(copybuf)
              end do
            end do
          case default
            do i = 1 , nlatin
              if (mod(i,10) == 0) print '(a$)', '.'
              do j = 1 , nlonin
                values(j,i) = readbuf((j-1)*ifrac+1,(i-1)*ifrac+1)
              end do
            end do
        end select
        print *, 'Done.'
        deallocate(readbuf)
        deallocate(copybuf)

      end subroutine read_ncglob
!
      subroutine checkerr(ierr)
        use netcdf
        implicit none
        integer , intent(in) :: ierr
        if ( ierr /= nf90_noerr ) then
          write (6, *) 'NetCDF library error.'
          write (6, *) nf90_strerror(ierr)
          stop
        end if
      end subroutine checkerr
!
      subroutine fillbuf(copybuf,readbuf,ni,nj,i,j,isize,lwrap)
        implicit none
        integer , intent(in) :: ni , nj , isize 
        real(4) , dimension(isize*isize) , intent(out) :: copybuf
        real(4) , dimension(ni,nj) , intent(in) :: readbuf
        integer , intent(in) :: i , j
        logical , intent(in) :: lwrap
        integer :: hsize , imin , imax , jmin , jmax , icnt , jcnt , ip
        integer :: ib , jb

        hsize = isize / 2
        imin = i - hsize
        imax = i + hsize
        jmin = j - hsize
        jmax = j + hsize

        ip = 1
        do icnt = 1 , isize
          do jcnt = 1 , isize
            ib = imin+icnt-1
            jb = jmin+jcnt-1
            if (lwrap) then
              if (ib < 1) then
                ib = ni - ib
              else if (ib > ni) then
                ib = ib - ni
              end if
            end if
            ib = min(max(ib,1),ni)
            jb = min(max(jb,1),nj)
            copybuf(ip) = readbuf(ib,jb)
            ip = ip + 1
          end do
        end do
      end subroutine fillbuf
!
      function mpindex(x)
        implicit none
        real(4) , dimension(:) , intent(in) :: x
        integer :: mpindex
        integer , dimension(32) :: cnt
        integer :: i

        do i = 1 , 32
          cnt(i) = sum(int(x)/i, int(x) == i)
        end do
        mpindex = 0
        do i = 1 , 32
          if (cnt(i) > mpindex) mpindex = i
        end do
      end function mpindex
!
      recursive subroutine qsort(a)
        implicit none 
        real(4) , dimension(:) , intent(in out) :: a
        integer :: np , isplit
 
        np = size(a)
        if (np > 1) then
         call partition(a, isplit)
         call qsort(a(:isplit-1))
         call qsort(a(isplit:))
        end if
      end subroutine qsort
 
      subroutine partition(a, marker)
        implicit none 
        real(4) , dimension(:) , intent(inout) :: a
        integer , intent(out) :: marker
        integer :: np , left , right
        real(4) :: temp , pivot
 
        np = size(a)
        pivot = (a(1) + a(np))/2.0E0
        left = 0
        right = np + 1
 
        do while (left < right)
          right = right - 1
          do while (a(right) > pivot)
            right = right-1
          end do
          left = left + 1
          do while (a(left) < pivot)
            left = left + 1
          end do
          if (left < right) then 
            temp = a(left)
            a(left) = a(right)
            a(right) = temp
          end if
        end do
        if (left == right) then
          marker = left + 1
        else
          marker = left
        end if
      end subroutine partition

      end module mod_rdldtr

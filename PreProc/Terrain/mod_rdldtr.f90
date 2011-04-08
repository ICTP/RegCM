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
  use m_stdio
  use m_die
  use m_mall
  use mod_constants

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
  subroutine read_ncglob(cfile,cvar,iires,iores,lreg,imeth)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: cfile , cvar
    integer , intent(in) :: iires , iores , imeth
    logical , intent(in) :: lreg
    integer :: ncid , ivar , istatus
    integer :: nlat , nlon , itl
    integer :: i , j , iosec , inpsec , iopsec , ifrac , ireg
    integer , dimension(2) :: istart, icount
    real(SP) , dimension(:,:) , allocatable :: readbuf
    real(SP) , dimension(:) , allocatable :: copybuf
    integer :: nlogb , nlagb , hnlogb , hnlagb , nfrac
    integer , parameter :: secpd = 3600
    integer , parameter :: secpm = 60
    real(DP) :: delta

    iosec = iores*secpm
    inpsec = secpd/iires
    nlogb = 360*inpsec
    nlagb = 180*inpsec
    if (mod(nlogb,iosec) /= 0) then
      write(stderr,*) 'Subroutine read_ncglob do not support'// &
               ' iores = ',iores
      call die('read_ncglob')
    end if
    hnlogb = nlogb/2
    hnlagb = nlagb/2
    iopsec = secpd/iosec
    ifrac = inpsec/iopsec

    ireg = 0
    if (lreg) ireg = 1
    delta = ((dble(iires)/d_two)/dble(secpd))*dble(ireg)

    grdltmn = dble(floor(xminlat))  -delta
    grdltma = dble(ceiling(xmaxlat))+delta
    if (grdltmn < -deg90+delta) grdltmn = -deg90+delta
    if (grdltma > +deg90-delta) grdltma =  deg90-delta
    if (lonwrap) then
      grdlnmn = -deg180+delta
      grdlnma =  deg180-delta
    else
      grdlnmn = dble(floor(xminlon))  -delta
      grdlnma = dble(ceiling(xmaxlon))+delta
    end if

    nlat = idnint((grdltma-grdltmn)*dble(inpsec))
    if (lonwrap) then
      nlon = nlogb + 1 - ireg
    else if (lcrosstime) then
      nlon = idnint(((deg180-delta-grdlnmn)+ &
                     (deg180-delta+grdlnma))*dble(inpsec))
    else
      nlon = idnint((grdlnma-grdlnmn)*dble(inpsec))
    end if

    nlatin = (nlat/ifrac)+1
    nlonin = (nlon/ifrac)
    nlat = nlat+1
    if (.not. lonwrap) then
      nlon = nlon+1
      nlonin = nlonin + 1
    end if

    call getspace

    allocate(readbuf(nlon,nlat), stat=istatus)
    if (istatus /= 0) then
      write(stderr,*) 'Memory error on allocating ', &
                      nlat*nlon*4,' bytes.'
      call die('read_ncglob')
    end if
    call mall_mci(readbuf,'read_ncglob')

    write(stdout,*) 'Opening '//trim(cfile)
    istatus = nf90_open(cfile, nf90_nowrite, ncid)
    call checkerr(istatus)
    istatus = nf90_inq_varid(ncid, cvar, ivar)
    call checkerr(istatus)

    istart(2) = hnlagb+idnint(grdltmn*dble(inpsec))+1
    if (lonwrap) then
      istart(1) = 1
    else
      istart(1) = hnlogb+idnint(grdlnmn*dble(inpsec))+1
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

    ! Fix South pole missing in dataset
    if (istart(2) == 1) then
      readbuf(:,1) = readbuf(:,2)
    end if

    write(stdout,'(a)',advance='no') ' Resampling'
    nfrac = ifrac*ifrac
    allocate(copybuf(nfrac), stat=istatus)
    if (istatus /= 0) then
      write(stderr,*) 'Memory error on allocating ',nfrac*4,' bytes'
      call die('read_ncglob')
    end if
    call mall_mci(copybuf,'read_ncglob')
    select case (imeth)
      case (1)
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                         (i-1)*ifrac+1,ifrac,lcrosstime)
            values(j,i) = sum(copybuf)/real(size(copybuf))
          end do
        end do
      case (2)
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                         (i-1)*ifrac+1,ifrac,lcrosstime)
            call qsort(copybuf)
            values(j,i) = 0.5*(copybuf(nfrac/2)+copybuf(nfrac/2+1))
          end do
        end do
      case (3)
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                         (i-1)*ifrac+1,ifrac,lcrosstime)
            values(j,i) = real(mpindex(copybuf))
          end do
        end do
      case default
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            values(j,i) = readbuf((j-1)*ifrac+1,(i-1)*ifrac+1)
          end do
        end do
    end select
    write(stdout,'(a)') 'Done.'
    deallocate(readbuf)
    call mall_mco(readbuf,'read_ncglob')
    deallocate(copybuf)
    call mall_mco(copybuf,'read_ncglob')

  end subroutine read_ncglob
!
  subroutine checkerr(ierr)
    use netcdf
    use m_die
    implicit none
    integer , intent(in) :: ierr
    if ( ierr /= nf90_noerr ) then
      call die('read_ncglob','NetCDF library error.', &
               ierr,nf90_strerror(ierr),0)
    end if
  end subroutine checkerr
!
  subroutine fillbuf(copybuf,readbuf,ni,nj,i,j,isize,lwrap)
    implicit none
    integer , intent(in) :: ni , nj , isize 
    real(SP) , dimension(isize*isize) , intent(out) :: copybuf
    real(SP) , dimension(ni,nj) , intent(in) :: readbuf
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
    real(SP) , dimension(:) , intent(in) :: x
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
    real(SP) , dimension(:) , intent(in out) :: a
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
    real(SP) , dimension(:) , intent(inout) :: a
    integer , intent(out) :: marker
    integer :: np , left , right
    real(SP) :: temp , pivot
 
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

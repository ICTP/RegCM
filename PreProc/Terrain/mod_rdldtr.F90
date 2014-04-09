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
  use mod_stdio
  use mod_constants
  use mod_memutil
  use mod_nchelper

  private

  real(rk8) , dimension(:,:) , pointer :: readbuf
  real(rk8) , dimension(:) , pointer :: copybuf
  real(rk8) , dimension(:,:) , pointer :: values

  public :: values , read_ncglob

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
    integer(ik4) , intent(in) :: iires , iores , imeth
    logical , intent(in) :: lreg
    integer(ik4) :: ncid , ivar , istatus
    integer(ik4) :: nlat , nlon , itl
    integer(ik4) :: i , j , iosec , inpsec , iopsec , ifrac , ireg
    integer(ik4) , dimension(2) :: istart, icount
    integer(ik4) :: nlogb , nlagb , hnlogb , hnlagb , nfrac
    integer(ik4) , parameter :: secpd = 3600
    integer(ik4) , parameter :: secpm = 60
    real(rk8) :: delta

    if ( iores <= 0 ) then
      iosec = iires
    else
      iosec = iores*secpm
    end if
    inpsec = secpd/iires
    nlogb = 360*inpsec
    nlagb = 180*inpsec
    if (mod(nlogb,iosec) /= 0) then
      write(stderr,*) 'Subroutine read_ncglob do not support iores = ',iores
      call die('read_ncglob')
    end if
    hnlogb = nlogb/2
    hnlagb = nlagb/2
    iopsec = secpd/iosec
    ifrac = inpsec/iopsec
#ifdef DEBUG
    write(stderr,*) 'INPSEC = ', inpsec
    write(stderr,*) 'IOPSEC = ', iopsec
#endif

    ireg = 0
    if (lreg) ireg = 1
    delta = ((dble(iires)/d_two)/dble(secpd))*dble(ireg)

    grdltmn = floor(xminlat)  -delta
    grdltma = ceiling(xmaxlat)+delta
    if (grdltmn < -deg90+delta) grdltmn = -deg90+delta
    if (grdltma > +deg90-delta) grdltma =  deg90-delta
    if (lonwrap) then
      grdlnmn = -deg180+delta
      grdlnma =  deg180-delta
    else
      grdlnmn = floor(xminlon)  -delta
      grdlnma = ceiling(xmaxlon)+delta
    end if

#ifdef DEBUG
    write(stderr,*) 'BOUNDS IN LAT: ', grdltmn , grdltma
    write(stderr,*) 'BOUNDS IN LON: ', grdlnmn , grdlnma
#endif

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

#ifdef DEBUG
    write(stderr,*) 'WILL READ ', nlon , 'x', nlat, ' points'
    write(stderr,*) 'WILL GIVE ', nlonin , 'x', nlatin, ' points'
#endif

    call getmem2d(values,1,nlonin,1,nlatin,'rdldtr:values')
    call getmem2d(readbuf,1,nlon,1,nlat,'rdldtr:readbuf')

    write(stdout,*) 'Opening '//trim(cfile)
    istatus = nf90_open(cfile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inq_varid(ncid, cvar, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    istart(2) = hnlagb+idnint(grdltmn*dble(inpsec))+1
    if (lonwrap) then
      istart(1) = 1
    else
      istart(1) = hnlogb+idnint(grdlnmn*dble(inpsec))+1
    end if
    if (.not. lcrosstime) then
      ! Simple case: not crossing timeline
      icount(2) = nlat
      if (icount(2)+istart(2) > nlagb) icount(2) = icount(2) - 1
      icount(1) = nlon
      istatus = nf90_get_var(ncid, ivar, readbuf, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    else
      ! Crossing timeline
      itl = nlogb - istart(1)
      icount(2) = nlat
      icount(1) = itl
      istatus = nf90_get_var(ncid, ivar, readbuf(1:itl,:), &
                             istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      istart(1) = 1
      icount(1) = nlon-itl
      istatus = nf90_get_var(ncid, ivar, readbuf(itl+1:nlon,:), &
                             istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    end if

    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    ! Fix Poles for interpolations
    if (istart(2) == 1) then
      write (stdout,*) 'Correcting South pole.'
      readbuf(:,1) = readbuf(:,2)
    else if (istart(2)+icount(2)-1 == nlagb) then
      write (stdout,*) 'Correcting North pole.'
      readbuf(:,nlat) = readbuf(:,nlat-1)
    end if

    write(stdout,'(a)',advance='no') ' Resampling'
    nfrac = ifrac*ifrac
    call getmem1d(copybuf,1,nfrac,'rdldtr:copybuf')
    select case (imeth)
      case (1)
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                         (i-1)*ifrac+1,ifrac,lcrosstime)
            values(j,i) = sum(copybuf)/dble(size(copybuf))
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
            values(j,i) = dble(mpindex(copybuf))
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
    call relmem2d(readbuf)
    call relmem1d(copybuf)
    write(stdout,'(a)') 'Done.'

  end subroutine read_ncglob
!
  subroutine fillbuf(copybuf,readbuf,ni,nj,i,j,isize,lwrap)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , isize 
    real(rk8) , dimension(isize*isize) , intent(out) :: copybuf
    real(rk8) , dimension(ni,nj) , intent(in) :: readbuf
    integer(ik4) , intent(in) :: i , j
    logical , intent(in) :: lwrap
    integer(ik4) :: hsize , imin , imax , jmin , jmax , icnt , jcnt , ip
    integer(ik4) :: ib , jb

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
    real(rk8) , dimension(:) , intent(in) :: x
    integer(ik4) :: mpindex
    integer(ik4) , dimension(32) :: cnt
    integer(ik4) :: i
    cnt = 0
    do i = 1 , 32
      cnt(i) = count(int(x) == i)
    end do
    mpindex = maxloc(cnt,1,cnt>0)
  end function mpindex
!
  recursive subroutine qsort(a)
    implicit none 
    real(rk8) , dimension(:) , intent(in out) :: a
    integer(ik4) :: np , isplit
 
    np = size(a)
    if (np > 1) then
     call partition(a, isplit)
     call qsort(a(:isplit-1))
     call qsort(a(isplit:))
    end if
  end subroutine qsort
 
  subroutine partition(a, marker)
    implicit none 
    real(rk8) , dimension(:) , intent(inout) :: a
    integer(ik4) , intent(out) :: marker
    integer(ik4) :: np , left , right
    real(rk8) :: temp , pivot
 
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

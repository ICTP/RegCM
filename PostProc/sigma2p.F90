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

program sigma2p
  use netcdf
  use mod_constants

  implicit none

  integer , parameter :: np = 11
  real(4) , dimension(np) :: plevs

  character(256) :: prgname , ncsfile , ncpfile
  character(128) :: attname , dimname , varname
  integer :: numarg , istatus , ncid , ncout

  integer , allocatable , dimension(:) :: dimids , dimlen
  real(4) , allocatable , dimension(:) :: sigma
  real(4) , allocatable , dimension(:,:,:) :: xvar
  real(4) , allocatable , dimension(:,:,:) :: pvar
  real(4) , allocatable , dimension(:,:) :: ps
  real(4) , allocatable , dimension(:) :: avar
  real(4) , allocatable , dimension(:) :: apvar
  real(8) , allocatable , dimension(:) :: times
  logical , allocatable , dimension(:) :: lkvarflag , ltvarflag
  integer , allocatable , dimension(:) :: varsize
  integer , allocatable , dimension(:) :: intscheme
  integer , allocatable , dimension(:) :: nvdims
  integer , allocatable , dimension(:,:) :: dimsize
  integer , allocatable , dimension(:) :: istart , icount
  integer :: ndims , nvars , natts , udimid , nvatts
  integer :: ivarid , idimid , xtype
  integer :: jxdimid , iydimid , kzdimid , itdimid , itvarid , ikvarid
  integer :: ipsvarid
  integer :: jx , iy , kz , nt
  real(4) :: ptop
  integer :: i , j , it , iv , iid1 , iid2 , ii , i3d , p3d
  integer :: n3d , ip3d
#ifdef __PGI
  integer , external :: iargc
#endif
#ifdef IBM 
  integer , external :: iargc
#endif

  data plevs /1000.,925.,850.,700.,500.,400.,300.,250.,200.,150.,100./

  call getarg(0, prgname)
  numarg = iargc( )
  if (numarg < 1) then
    write (6,*) 'Not enough arguments.'
    write (6,*) ' '
    write (6,*) 'Usage : ', trim(prgname), ' Rcmfile.nc'
    write (6,*) ' '
    stop
  end if

  call getarg(1, ncsfile)
  iid1 = scan(ncsfile, '/', .true.)
  iid2 = scan(ncsfile, '.', .true.)
  ncpfile = trim(ncsfile(iid1+1:iid2-1))//'_pressure.nc'

  istatus = nf90_open(ncsfile, nf90_nowrite, ncid)
  call check_ok('Error Opening Input file '//trim(ncsfile))

#ifdef NETCDF4_HDF5
  istatus = nf90_create(ncpfile, &
            ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model), &
            ncout)
#else
  istatus = nf90_create(ncpfile, nf90_clobber, ncout)
#endif
  call check_ok('Error Opening Output file '//trim(ncpfile))

  istatus = nf90_inquire(ncid,ndims,nvars,natts,udimid)
  call check_ok('Error Reading Netcdf file '//trim(ncsfile))

  do i = 1 , natts
    istatus = nf90_inq_attname(ncid, nf90_global, i, attname)
    call check_ok('Error read global attribute')
    istatus = nf90_copy_att(ncid, nf90_global, attname, ncout, nf90_global)
    call check_ok('Error copy global attribute '//trim(attname))
  end do

  allocate(dimlen(ndims), stat=istatus)
  call check_alloc('Memory error allocating dimlen')
 
  kz = 0
  nt = 0
  do i = 1 , ndims
    istatus = nf90_inquire_dimension(ncid, i, dimname, dimlen(i))
    call check_ok('Error reading dimension info')
    if (dimname == 'iy' .or. dimname == 'y') then
      iy = dimlen(i)
      iydimid = i
    else if (dimname == 'jx' .or. dimname == 'x') then
      jx = dimlen(i)
      jxdimid = i
    else if (dimname == 'kz' .or. dimname == 'z') then
      kz = dimlen(i)
      kzdimid = i
    end if
    if (dimname == 'time') then
      itdimid = i
      nt = dimlen(i)
      istatus = nf90_def_dim(ncout, dimname, nf90_unlimited, idimid)
    else if (dimname == 'kz') then
      istatus = nf90_def_dim(ncout, 'plev', np, idimid)
    else
      istatus = nf90_def_dim(ncout, dimname, dimlen(i), idimid)
    end if
    call check_ok('Error define dimension '//trim(dimname))
  end do

  i3d = jx*iy*kz
  p3d = jx*iy*np

  if (kz == 0 .or. nt == 0) then
    write (6,*) 'Nothing to do: no vertical dimension or no time'
    istatus = nf90_close(ncout)
    istatus = nf90_close(ncid)
    stop
  end if

  allocate(dimids(ndims), stat=istatus)
  call check_alloc('Memory error allocating dimids')
  allocate(lkvarflag(nvars), stat=istatus)
  call check_alloc('Memory error allocating lkvarflag')
  allocate(ltvarflag(nvars), stat=istatus)
  call check_alloc('Memory error allocating ltvarflag')
  allocate(varsize(nvars), stat=istatus)
  call check_alloc('Memory error allocating varsize')
  allocate(intscheme(nvars), stat=istatus)
  call check_alloc('Memory error allocating intscheme')
  allocate(nvdims(nvars), stat=istatus)
  call check_alloc('Memory error allocating nvdims')
  allocate(dimsize(ndims,nvars), stat=istatus)
  call check_alloc('Memory error allocating dimsize')
  allocate(istart(ndims), stat=istatus)
  call check_alloc('Memory error allocating istart')
  allocate(icount(ndims), stat=istatus)
  call check_alloc('Memory error allocating icount')
  allocate(ps(jx,iy), stat=istatus)
  call check_alloc('Memory error allocating ps')
  allocate(xvar(jx,iy,kz), stat=istatus)
  call check_alloc('Memory error allocating xvar')
  allocate(pvar(jx,iy,np), stat=istatus)
  call check_alloc('Memory error allocating pvar')
  allocate(sigma(kz), stat=istatus)
  call check_alloc('Memory error allocating sigma')
  allocate(times(nt), stat=istatus)
  call check_alloc('Memory error allocating times')

  itvarid = 0
  do i = 1 , nvars
    lkvarflag(i) = .false.
    ltvarflag(i) = .false.
    varsize(i) = 1
    intscheme(i) = 1
    dimsize(:,i) = 1
    istatus = nf90_inquire_variable(ncid,i,varname,xtype,nvdims(i), &
                                    dimids,nvatts)
    call check_ok('Error inquire variable')
    if (varname == 'time') then
      itvarid = i
    else if (varname == 'sigma') then
      varname = 'plev'
      ikvarid = i
    else if (varname == 'ptop') then
      istatus = nf90_get_var(ncid,i,ptop)
      call check_ok('Error read variable ptop')
    else if (varname == 'ps') then
      ipsvarid = i
    else if (varname == 't') then
      intscheme(i) = 2
    end if
    iv = 1
    do j = 1 , nvdims(i)
      if (dimids(j) == kzdimid) then
        lkvarflag(i) = .true.
      else if (dimids(j) == itdimid) then
        ltvarflag(i) = .true.
        iv = iv + 1
        cycle
      end if
      dimsize(iv,i) = dimlen(dimids(j))
      iv = iv + 1
      varsize(i) = varsize(i) * dimlen(dimids(j))
    end do
    istatus = nf90_def_var(ncout, varname, xtype, dimids(1:nvdims(i)), ivarid)
    call check_ok('Error define variable '//trim(varname))
#ifdef NETCDF4_HDF5
    if (nvdims(i) > 2) then
      istatus = nf90_def_var_deflate(ncout, ivarid, 1, 1, 9)
      call check_ok('Error setting compression for '//trim(varname))
    end if
#endif
    if (varname == 'plev') then
      istatus = nf90_put_att(ncout, ivarid, 'standard_name', 'pressure')
      call check_ok('Error adding pressure standard name')
      istatus = nf90_put_att(ncout, ivarid, 'long_name', &
                             'Interpolated pressure')
      call check_ok('Error adding pressure long name')
      istatus = nf90_put_att(ncout, ivarid, 'units', 'hPa')
      call check_ok('Error adding pressure units')
      istatus = nf90_put_att(ncout, ivarid, 'axis', 'Z')
      call check_ok('Error adding pressure axis')
      istatus = nf90_put_att(ncout, ivarid, 'positive', 'down')
      call check_ok('Error adding pressure positive')
      cycle
    end if
    do j = 1 , nvatts
      istatus = nf90_inq_attname(ncid, i, j, attname)
      call check_ok('Error read attribute for '//trim(varname))
      istatus = nf90_copy_att(ncid, i, attname, ncout, ivarid)
      call check_ok('Error copy attribute '//trim(attname)//' for '// &
                    trim(varname))
    end do
  end do

  istatus = nf90_inq_varid(ncid, "sigma", ivarid)
  call check_ok('Error reading variable sigma.')
  istatus = nf90_get_var(ncid, ivarid, sigma)
  call check_ok('Error reading variable sigma.')

  istatus = nf90_inq_varid(ncid, "time", ivarid)
  call check_ok('Error reading variable time.')
  istatus = nf90_get_var(ncid, ivarid, times)
  call check_ok('Error reading variable time.')

  istatus = nf90_enddef(ncout)
  call check_ok('Error preparing for write on output')

! Write time independent variables

  do i = 1 , nvars
    if (i == itvarid) cycle
    if (i == ikvarid) then
      istatus = nf90_put_var(ncout, i, plevs)
      call check_ok('Error writing variable plev')
      cycle
    end if
    if (.not. ltvarflag(i)) then
      allocate(avar(varsize(i)), stat=istatus)
      call check_alloc('Memory error allocating avar')
      iv = nvdims(i)
      istart(:) = 1
      icount(1:iv) = dimsize(1:iv,i)
      istatus = nf90_get_var(ncid, i, avar, istart(1:iv), icount(1:iv))
      call check_ok('Error reading variable.')
      istatus = nf90_put_var(ncout, i, avar, istart(1:iv), icount(1:iv))
      call check_ok('Error writing variable.')
      deallocate(avar)
    end if
  end do

! Write time dependent variables

  do it = 1 , nt
    istart(1) = it
    icount(1) = 1
    istatus = nf90_put_var(ncout, itvarid, times(it:it), &
                           istart(1:1), icount(1:1))
    call check_ok('Error writing time.')
    istart(1) = 1
    istart(2) = 1
    istart(3) = it
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1
    istatus = nf90_get_var(ncid, ipsvarid, ps, istart(1:3), icount(1:3))
    call check_ok('Error reading ps.')
    istatus = nf90_put_var(ncout, ipsvarid, ps, istart(1:3), icount(1:3))
    call check_ok('Error writing ps.')
    do i = 1 , nvars
      if (.not. ltvarflag(i)) cycle
      if (i == itvarid) cycle
      if (i == ipsvarid) cycle

      if (lkvarflag(i)) then

!       Do interpolation

        iv = nvdims(i)
        istart(iv) = it
        icount(iv) = 1
        istart(1:iv-1) = 1
        icount(1:iv-1) = dimsize(1:iv-1,i)
        allocate(avar(varsize(i)),stat=istatus)
        call check_alloc('Memory error allocating avar')
        istatus = nf90_get_var(ncid, i, avar, istart(1:iv), icount(1:iv))
        call check_ok('Error reading variable to interpolate.')

        n3d = varsize(i) / i3d
        ip3d = p3d*n3d

        allocate(apvar(ip3d),stat=istatus)
        call check_alloc('Memory error allocating apvar')
        do ii = 1 , n3d
          xvar = reshape(avar((ii-1)*i3d+1:ii*i3d),(/jx,iy,kz/))
          if (intscheme(i) == 1) then
            call intlin(pvar,xvar,ps,sigma,jx,iy,kz,plevs,np)
          else if (intscheme(i) == 2) then
            call intlog(pvar,xvar,ps,sigma,jx,iy,kz,plevs,np)
          end if
          apvar((ii-1)*ip3d+1:ii*ip3d) = reshape(pvar,(/ip3d/))
        end do

        icount(3) = np
        istatus = nf90_put_var(ncout, i, apvar, istart(1:iv), icount(1:iv))
        call check_ok('Error writing interpolated variable.')

        deallocate(avar)
        deallocate(apvar)

      else

!       No interpolation

        iv = nvdims(i)
        istart(iv) = it
        icount(iv) = 1
        istart(1:iv-1) = 1
        icount(1:iv-1) = dimsize(1:iv-1,i)
        allocate(avar(varsize(i)),stat=istatus)
        call check_alloc('Memory error allocating avar')
        istatus = nf90_get_var(ncid, i, avar, istart(1:iv), &
                               icount(1:iv))
        call check_ok('Error reading variable to passthrough')
        istatus = nf90_put_var(ncout, i, avar, istart(1:iv), &
                               icount(1:iv))
        call check_ok('Error writing variable to passthrough')
        deallocate(avar)

      end if

    end do
  end do

  deallocate(dimlen)
  deallocate(dimids)
  deallocate(lkvarflag)
  deallocate(ltvarflag)
  deallocate(varsize)
  deallocate(nvdims)
  deallocate(dimsize)
  deallocate(sigma)
  deallocate(times)
  deallocate(ps)
  deallocate(xvar)
  deallocate(pvar)

  istatus = nf90_close(ncid)
  call check_ok('Error close input file '//trim(ncsfile))

  istatus = nf90_close(ncout)
  call check_ok('Error close output file '//trim(ncpfile))

contains
!
  subroutine intlog(fp,f,pstar,sig,im,jm,km,p,kp)
    implicit none
    integer , intent(in) :: im , jm , km , kp
    real(4) , intent(out) , dimension(im,jm,kp) :: fp
    real(4) , intent(in) , dimension(im,jm) :: pstar
    real(4) , intent(in) , dimension(km) :: sig
    real(4) , intent(in) , dimension(kp) :: p
    real(4) , intent(in) , dimension(im,jm,km) :: f
    integer :: i , j , k , n
    integer :: k1 , k1p , kbc
    real(4) :: sigp , wp , w1
!
! intlog is for vertical interpolation of t.  the interpolation is
! linear in log p.  where extrapolation upward is necessary,
! the t field is considered to have 0 vertical derivative.
! where extrapolation downward is necessary, the t field is
! considered to have a lapse rate of lrate (k/m), and the
! thickness is determined hydrostatically from the mean of the
! two extreme temperatures in the layer.

!   find first sigma level above boundary layer (less than sig=bltop)
    do k = 1 , km
      if (sig(k).lt.bltop) kbc = k
    end do
    do j = 1 , jm
      do i = 1 , im
        do n = 1 , kp
          sigp = (p(n)-ptop) / (pstar(i,j)-ptop)
          k1 = 0
          do k = 1 , km
            if (sigp.gt.sig(k)) k1 = k
          end do
          if (sigp.le.sig(1)) then
            fp(i,j,n) = f(i,j,1)
          else if ((sigp.gt.sig(1)) .and. (sigp.lt.sig(km))) then
            k1p = k1 + 1
            wp  = log(sigp/sig(k1)) / log(sig(k1p)/sig(k1))
            w1  = 1. - wp
            fp(i,j,n)= w1*f(i,j,k1) + wp*f(i,j,k1p)
          else if ((sigp.ge.sig(km)) .and. (sigp.le.1.))then
            fp(i,j,n) = f(i,j,km)
          else if (sigp.gt.1.) then
            fp(i,j,n) = f(i,j,kbc) &
              * exp(-rgas*lrate*log(sigp/sig(kbc))*rgti)
!           ***** from r. errico, see routine height *****
          end if
        end do
      end do
    end do
  end subroutine intlog
!
  subroutine intlin(fp,f,pstar,sig,im,jm,km,p,kp)
    implicit none
    integer , intent(in) :: im , jm , km , kp
    real(4) , intent(out) , dimension(im,jm,kp) :: fp
    real(4) , intent(in) , dimension(im,jm) :: pstar
    real(4) , intent(in) , dimension(km) :: sig
    real(4) , intent(in) , dimension(kp) :: p
    real(4) , intent(in) , dimension(im,jm,km) :: f
    integer :: i , j , k , n
    integer :: k1 , k1p
    real(4) :: sigp , wp , w1
!
!  intlin is for vertical interpolation of u, v, and relative humidity.
!   the interpolation is linear in p.  where extrapolation is
!   necessary, fields are considered to have 0 vertical derivative.

    do j = 1 , jm
      do i = 1 , im
        do n = 1 , kp
          sigp = (p(n)-ptop) / (pstar(i,j)-ptop)
          k1 = 0
          do k = 1 , km
            if (sigp.gt.sig(k)) k1 = k
          end do
          if (sigp.le.sig(1)) then
            fp(i,j,n) = f(i,j,1)
          else if ((sigp.gt.sig(1)) .and. (sigp.lt.sig(km))) then
            k1p = k1 + 1
            wp  = (sigp-sig(k1))/(sig(k1p)-sig(k1))
            w1  = 1.-wp
            fp(i,j,n)  = w1*f(i,j,k1) + wp*f(i,j,k1p)
          else if (sigp.ge.sig(km)) then
            fp(i,j,n)  = f(i,j,km)
          end if
        end do
      end do
    end do
  end subroutine intlin
!
  subroutine check_ok(message)
    use netcdf
    implicit none
    character(*) :: message
    if (istatus /= nf90_noerr) then
      write (6,*) message
      write (6,*) nf90_strerror(istatus)
      stop
    end if
  end subroutine check_ok

  subroutine check_alloc(message)
    use netcdf
    implicit none
    character(*) :: message
    if (istatus /= 0) then
      write (6,*) message
      stop
    end if
  end subroutine check_alloc

end program sigma2p

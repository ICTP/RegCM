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

! Unit test of Share/mod_rdldtr.F90. A dataset on a regular 30 arc second
! grid is written in the current directory; read_ncglob reads the window
! around a domain and resamples it on blocks of 4x4 cells (2 arc minutes):
! each method keeps a constant, the mean, median and trimmed mean of a
! linear ramp are its value at the centre of the block, and the bounds of
! the resampled grid are the centres of the blocks. gfopen and gfread
! interpolate the dataset with mod_kdinterp, with and without a land mask.

program test_rdldtr
  use netcdf
  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_rdldtr
  use mod_kdinterp, only : h_missing_value
  use share_test_utils
  implicit none
  character(len=*), parameter :: fname = 'test_rdldtr.nc'
  integer(ik4), parameter :: ng = 600
  real(rk8), parameter :: dg = 1.0_rk8/120.0_rk8
  integer(ik4), parameter :: nj = 6, ni = 5
  real(rkx), dimension(nj,ni) :: xlat, xlon, var
  real(rkx), pointer, contiguous, dimension(:,:) :: values => null()
  real(rk8) :: lnma, lnmn, ltma, ltmn, err, x0, y0
  integer(ik4) :: nlatin, nlonin, i, j, m
  type(globalfile) :: gfile
  character(len=48) :: msg

  call memory_init
  do i = 1, ni
    do j = 1, nj
      xlat(j,i) = 1.0_rkx + 0.5_rkx*real(i,rkx)
      xlon(j,i) = 1.0_rkx + 0.5_rkx*real(j,rkx)
    end do
  end do

  ! A constant field: every resampling method keeps it.
  call write_dataset(.false.,.false.)
  do m = 0, 4
    call read_ncglob(fname,'z',2,m,0,xlat,xlon,lnma,lnmn,ltma,ltmn, &
                     nlatin,nlonin,values)
    write(msg,'(a,i0,a)') 'resampling method ', m, ' keeps a constant'
    if ( m == 3 ) then
      call check(all(values == 7.0_rkx),trim(msg))
    else
      call check(all(abs(values-7.0_rkx) < 1.0e-5_rkx),trim(msg))
    end if
  end do

  ! A ramp: the block statistics are the centre value, and
  ! the resampled grid starts at the centre of the first block.
  call write_dataset(.true.,.false.)
  do m = 1, 4
    if ( m == 3 ) cycle
    call read_ncglob(fname,'z',2,m,0,xlat,xlon,lnma,lnmn,ltma,ltmn, &
                     nlatin,nlonin,values)
    err = 0.0_rk8
    do i = 1, nlatin
      y0 = ltmn + real(i-1,rk8)*2.0_rk8/60.0_rk8
      do j = 1, nlonin
        x0 = lnmn + real(j-1,rk8)*2.0_rk8/60.0_rk8
        err = max(err,abs(real(values(j,i),rk8)-ramp(x0,y0)))
      end do
    end do
    write(msg,'(a,i0,a)') 'resampling method ', m, ' of a ramp'
    call check(err < 1.0e-3_rk8,trim(msg))
  end do
  call check(abs((lnma-lnmn)-real(nlonin-1,rk8)*2.0_rk8/60.0_rk8) &
             < 1.0e-9_rk8 .and. abs((ltma-ltmn) - &
             real(nlatin-1,rk8)*2.0_rk8/60.0_rk8) < 1.0e-9_rk8, &
    'bounds of the resampled grid')
  call read_ncglob(fname,'z',0,0,0,xlat,xlon,lnma,lnmn,ltma,ltmn, &
                   nlatin,nlonin,values)
  err = 0.0_rk8
  do j = 1, nlonin
    x0 = lnmn + real(j-1,rk8)*dg
    err = max(err,abs(real(values(j,1),rk8)-ramp(x0,ltmn)))
  end do
  call check(err < 1.0e-3_rk8,'window read without resampling')

  ! Interpolation with gfread: the ramp, then a masked dataset.
  call gfopen(gfile,fname,xlat,xlon,50.0_rkx,1.0_rkx,0)
  call gfread(gfile,'z',var,0.0_rkx)
  err = 0.0_rk8
  do i = 1, ni
    do j = 1, nj
      err = max(err,abs(real(var(j,i),rk8) - &
                ramp(real(xlon(j,i),rk8),real(xlat(j,i),rk8))))
    end do
  end do
  call check(err < 0.05_rk8,'gfread of a ramp')
  call gfclose(gfile)
  call write_dataset(.false.,.true.)
  call gfopen(gfile,fname,xlat,xlon,50.0_rkx,1.0_rkx,0)
  call gfread(gfile,'z',var,0.0_rkx)
  call check(all(abs(var-7.0_rkx) < 1.0e-4_rkx .or. &
                 var == h_missing_value),'gfread with a land mask')
  call check(any(var == h_missing_value) .and. &
             any(abs(var-7.0_rkx) < 1.0e-4_rkx), &
    'masked points are missing')
  call gfclose(gfile)

  call memory_destroy
  call test_end('test_rdldtr')

  contains

  ! Distinct values in a block, symmetric about its centre
  real(rk8) function ramp(x,y)
    implicit none
    real(rk8), intent(in) :: x, y
    ramp = 100.0_rk8 + 10.0_rk8*x + 1.3_rk8*y
  end function ramp

  ! Dataset from 0 to 5 degrees in latitude and longitude, 30 arc
  ! seconds: a constant 7, or a ramp; with a land mask
  ! the western half is sea.
  subroutine write_dataset(lramp,lmask)
    implicit none
    logical, intent(in) :: lramp, lmask
    integer(ik4) :: ncid, ilat, ilon, ivlat, ivlon, ivz, ivm, istat
    real(rk8), dimension(ng) :: lat, lon
    real(rk4), dimension(ng,ng) :: z
    integer(ik4), dimension(ng,ng) :: mask
    integer(ik4) :: k, l
    do k = 1, ng
      lat(k) = real(k-1,rk8)*dg
      lon(k) = real(k-1,rk8)*dg
    end do
    do l = 1, ng
      do k = 1, ng
        if ( lramp ) then
          z(k,l) = real(ramp(lon(k),lat(l)),rk4)
        else
          z(k,l) = 7.0_rk4
        end if
        mask(k,l) = merge(0,1,lon(k) < 2.5_rk8)
      end do
    end do
    istat = nf90_create(fname,nf90_clobber,ncid)
    istat = nf90_def_dim(ncid,'lat',ng,ilat)
    istat = nf90_def_dim(ncid,'lon',ng,ilon)
    istat = nf90_def_var(ncid,'lat',nf90_double,[ilat],ivlat)
    istat = nf90_def_var(ncid,'lon',nf90_double,[ilon],ivlon)
    istat = nf90_def_var(ncid,'z',nf90_float,[ilon,ilat],ivz)
    if ( lmask ) istat = nf90_def_var(ncid,'landmask',nf90_int, &
                                      [ilon,ilat],ivm)
    istat = nf90_enddef(ncid)
    istat = nf90_put_var(ncid,ivlat,lat)
    istat = nf90_put_var(ncid,ivlon,lon)
    istat = nf90_put_var(ncid,ivz,z)
    if ( lmask ) istat = nf90_put_var(ncid,ivm,mask)
    istat = nf90_close(ncid)
  end subroutine write_dataset

end program test_rdldtr

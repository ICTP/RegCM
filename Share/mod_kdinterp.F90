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

module mod_kdinterp

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_stdio
  use mod_message
  use mod_earth
  use mod_spbarcoord
  use mod_dynparam, only : idynamic
  use mod_kdtree2

  private

  public :: h_interpolator
  public :: h_interpolator_create
  public :: h_interpolate_cont, h_interpolate_class, h_interpolate_nn
  public :: h_interpolator_destroy
  public :: h_missing_value

  interface h_interpolator_create
    module procedure interp_create_ll_g
    module procedure interp_create_g_g
    module procedure interp_create_g_ll
    module procedure interp_create_ll_ll
    module procedure interp_create_ll_ll_1d
  end interface h_interpolator_create

  interface h_interpolate_cont
    module procedure interp_1d
    module procedure interp_2d
    module procedure interp_3d
    module procedure interp_4d
    module procedure interp_5d
  end interface h_interpolate_cont

  interface h_interpolate_nn
    module procedure interp_1d_nn
    module procedure interp_2d_nn
  end interface h_interpolate_nn

  interface h_interpolate_class
    module procedure interp_class_2dr
    module procedure interp_class_3dr
    module procedure interp_class_4dr
    module procedure interp_class_5dr
    module procedure interp_class_i
    module procedure interp_class_ld
  end interface h_interpolate_class

  interface compwgt_genlin
    module procedure compwgt_genlin_1d
    module procedure compwgt_genlin_2d
  end interface compwgt_genlin

  interface maxedis
    module procedure maxedis1
    module procedure maxedis2
  end interface maxedis

  ! Need at least three point to triangulate, 4 for a quasi-linear
  integer(ik4), parameter :: minp = 4

  ! Try not to chocke memory...
  integer(ik4), parameter :: maxp = 32

  real(rkx), parameter :: missl = -9999.0_rkx
  real(rkx), parameter :: h_missing_value = missl
  real(rkx), parameter :: missc = -9990.0_rkx
  real(rk8), parameter :: mindis = 1.0e-4_rkx

  type pwgt
    integer(ik4) :: i, j
    real(rk8) :: wgt
  end type pwgt

  type ijwgt
    integer(ik4) :: np
    type(pwgt), dimension(:), pointer, contiguous :: wgt => null( )
  end type ijwgt

  type ftarget
    integer(ik4), dimension(2) :: tshape
    type(ijwgt), dimension(:,:), pointer, contiguous :: ft => null( )
  end type ftarget

  type h_interpolator
    integer(ik4), dimension(2) :: sshape
    type(ftarget) :: tg
  end type h_interpolator

  contains

  real(rk8) function maxedis1(lat,lon) result(maxdis)
    implicit none
    real(rk8), dimension(:), intent(in) :: lat
    real(rk8), dimension(:), intent(in) :: lon
    integer(ik4) :: i, j
    real(rk8) :: r1, r2
    maxdis = 0.1_rk8
    do i = 2, size(lat)
      do j = 2, size(lon)
        r1 = gcdist_simple(lat(i-1),lon(j-1),lat(i),lon(j))
        r2 = gcdist_simple(lat(i-1),lon(j),lat(i),lon(j-1))
        if ( r1 > maxdis ) maxdis = r1
        if ( r2 > maxdis ) maxdis = r2
      end do
    end do
  end function maxedis1

  real(rk8) function maxedis2(lat,lon) result(maxdis)
    implicit none
    real(rk8), dimension(:,:), intent(in) :: lat
    real(rk8), dimension(:,:), intent(in) :: lon
    integer(ik4) :: i, j
    integer(ik4), dimension(2) :: nij
    real(rk8) :: r1, r2
    maxdis = 0.1_rk8
    nij = shape(lat)
    do i = 2, nij(2)
      do j = 2, nij(1)
        r1 = gcdist_simple(lat(j-1,i-1),lon(j-1,i-1),lat(j,i),lon(j,i))
        r2 = gcdist_simple(lat(j,i-1),lon(j,i-1),lat(j-1,i),lon(j-1,i))
        if ( r1 > maxdis ) maxdis = r1
        if ( r2 > maxdis ) maxdis = r2
      end do
    end do
  end function maxedis2

  subroutine interp_create_ll_g(h_i,slat,slon,tlat,tlon,ds,roi)
    implicit none
    type(h_interpolator), intent(out) :: h_i
    real(rkx), dimension(:), intent(in) :: slat
    real(rkx), dimension(:), intent(in) :: slon
    real(rkx), dimension(:,:), intent(in) :: tlat
    real(rkx), dimension(:,:), intent(in) :: tlon
    real(rkx), intent(in), optional :: ds
    real(rkx), intent(in), optional :: roi
    real(rk8), dimension(:,:), allocatable :: x
    real(kdkind), dimension(3) :: p
    type(kdtree2), pointer :: mr
    type(kdtree2_result), pointer, contiguous, dimension(:) :: results
    integer(ik4) :: n1, n2, np, ni, nj, nf, i, j, i10
    real(rk8) :: dx, r2, rin
    logical :: lbil

    if ( any(shape(tlat) /= shape(tlon)) ) then
      call die('interp_create_ll_g','Target shapes non conforming',1)
    end if
    n1 = size(slat)
    n2 = size(slon)
    h_i%sshape(1) = n2
    h_i%sshape(2) = n1
    np = n1 * n2
    allocate(x(3,np))
    call ll2xyz(slat,slon,x)
    mr => kdtree2_create(x,sort=.true.,rearrange=.true.)
    rin = maxedis(slat,slon)
    dx = rin/erkm
    lbil = .true.
    if ( present(ds) ) then
      if ( rin < ds ) then
        dx = ds/erkm
        lbil = .false.
      end if
    end if
    if ( present(roi) ) then
      r2 = (dx*dx*roi*roi)
    else
      r2 = (dx*dx)
    end if
    nj = size(tlat,1)
    ni = size(tlat,2)
    h_i%tg%tshape = shape(tlat)
    allocate(h_i%tg%ft(nj,ni))
    write(stdout,'(a,f16.8,a)',advance='no') &
      ' Computing weights, R = ',sqrt(r2)*erkm, ' '
    i10 = ni/10
    do i = 1, ni
      if (mod(i,i10) == 0) write(stdout,'(a)',advance='no') '.'
      do j = 1, nj
        call ll2xyz(tlat(j,i),tlon(j,i),p)
        if ( lbil ) then
          np = minp
          allocate(results(minp))
          call kdtree2_n_nearest(mr,p,np,results)
        else
          np = kdtree2_r_count(mr,p,r2)
          if ( np < minp ) then
            h_i%tg%ft(j,i)%np = 0
            cycle
          else
            if ( np > maxp ) then
              np = maxp
              allocate(results(maxp))
              call kdtree2_n_nearest(mr,p,np,results)
            else
              allocate(results(np))
              call kdtree2_r_nearest(mr,p,r2,nf,np,results)
              np = nf
            end if
          end if
        end if
        call compwgt_genlin(np,n2,p,x,results,h_i%tg%ft(j,i)%wgt)
        !call compwgt_distwei(np,n2,results,h_i%tg%ft(j,i)%wgt)
        h_i%tg%ft(j,i)%np = np
        deallocate(results)
      end do
    end do
    deallocate(x)
    call kdtree2_destroy(mr)
    write(stdout,'(a)') ' Done.'
  end subroutine interp_create_ll_g

  subroutine interp_create_ll_ll(h_i,slat,slon,tlat,tlon,roi)
    implicit none
    type(h_interpolator), intent(out) :: h_i
    real(rkx), dimension(:), intent(in) :: slat
    real(rkx), dimension(:), intent(in) :: slon
    real(rkx), dimension(:), intent(in) :: tlat
    real(rkx), dimension(:), intent(in) :: tlon
    real(rkx), intent(in), optional :: roi
    real(rk8), dimension(:,:), allocatable :: x
    real(kdkind), dimension(3) :: p
    type(kdtree2), pointer :: mr
    type(kdtree2_result), pointer, contiguous, dimension(:) :: results
    integer(ik4) :: n1, n2, np, ni, nj, nf, i, j, i10
    real(rk8) :: dx, dx1, r2, rin, rin1
    logical :: lbil

    n1 = size(slat)
    n2 = size(slon)
    h_i%sshape(1) = n2
    h_i%sshape(2) = n1
    np = n1 * n2
    allocate(x(3,np))
    call ll2xyz(slat,slon,x)
    mr => kdtree2_create(x,sort=.true.,rearrange=.true.)
    rin = maxedis(slat,slon)
    dx = rin/erkm
    rin1 = maxedis(tlat,tlon)
    dx1 = rin1/erkm
    lbil = .true.
    if ( rin < rin1 ) then
      dx = dx1
      lbil = .false.
    end if
    if ( present(roi) ) then
      r2 = (dx*dx*roi*roi)
    else
      r2 = (dx*dx)
    end if
    ni = size(tlat)
    nj = size(tlon)
    h_i%tg%tshape(1) = nj
    h_i%tg%tshape(2) = ni
    allocate(h_i%tg%ft(nj,ni))
    write(stdout,'(a,f16.8,a)',advance='no') &
      ' Computing weights, R = ',sqrt(r2)*erkm, ' '
    i10 = ni/10
    do i = 1, ni
      if (mod(i,i10) == 0) write(stdout,'(a)',advance='no') '.'
      do j = 1, nj
        call ll2xyz(tlat(i),tlon(j),p)
        if ( lbil ) then
          np = minp
          allocate(results(minp))
          call kdtree2_n_nearest(mr,p,np,results)
        else
          np = kdtree2_r_count(mr,p,r2)
          if ( np < minp ) then
            h_i%tg%ft(j,i)%np = 0
            cycle
          else
            if ( np > maxp ) then
              np = maxp
              allocate(results(maxp))
              call kdtree2_n_nearest(mr,p,np,results)
            else
              allocate(results(np))
              call kdtree2_r_nearest(mr,p,r2,nf,np,results)
              np = nf
            end if
          end if
        end if
        call compwgt_genlin(np,n2,p,x,results,h_i%tg%ft(j,i)%wgt)
        !call compwgt_distwei(np,n2,results,h_i%tg%ft(j,i)%wgt)
        h_i%tg%ft(j,i)%np = np
        deallocate(results)
      end do
    end do
    deallocate(x)
    call kdtree2_destroy(mr)
    write(stdout,'(a)') ' Done.'
  end subroutine interp_create_ll_ll

  subroutine interp_create_ll_ll_1d(h_i,ni,slat,slon,no,tlat,tlon)
    implicit none
    integer, intent(in) :: ni, no
    type(h_interpolator), intent(out) :: h_i
    real(rkx), dimension(ni), intent(in) :: slat
    real(rkx), dimension(ni), intent(in) :: slon
    real(rkx), dimension(no), intent(in) :: tlat
    real(rkx), dimension(no), intent(in) :: tlon
    real(rk8), dimension(:,:), allocatable :: x
    real(kdkind), dimension(3) :: p
    type(kdtree2), pointer :: mr
    type(kdtree2_result), pointer, contiguous, dimension(:) :: results
    integer(ik4) :: i, i10, np

    h_i%sshape(1) = ni
    h_i%sshape(2) = 1
    allocate(x(3,ni))
    call ll2xyz(ni,slat,slon,x)
    mr => kdtree2_create(x,sort=.true.,rearrange=.true.)
    h_i%tg%tshape(1) = no
    h_i%tg%tshape(2) = 1
    allocate(h_i%tg%ft(no,1))
    write(stdout,'(a)',advance='no') ' Computing weights'
    i10 = no/10
    do i = 1, no
      if (mod(i,i10) == 0) write(stdout,'(a)',advance='no') '.'
      call ll2xyz(tlat(i),tlon(i),p)
      np = minp
      allocate(results(np))
      call kdtree2_n_nearest(mr,p,np,results)
      call compwgt_genlin(np,p,x,results,h_i%tg%ft(i,1)%wgt)
      h_i%tg%ft(i,1)%np = np
      deallocate(results)
    end do
    deallocate(x)
    call kdtree2_destroy(mr)
    write(stdout,'(a)') ' Done.'
  end subroutine interp_create_ll_ll_1d

  subroutine interp_create_g_g(h_i,slat,slon,tlat,tlon,ds,roi)
    implicit none
    type(h_interpolator), intent(out) :: h_i
    real(rkx), dimension(:,:), intent(in) :: slat
    real(rkx), dimension(:,:), intent(in) :: slon
    real(rkx), dimension(:,:), intent(in) :: tlat
    real(rkx), dimension(:,:), intent(in) :: tlon
    real(rkx), intent(in), optional :: ds
    real(rkx), intent(in), optional :: roi
    real(rk8), dimension(:,:), allocatable :: x
    real(kdkind), dimension(3) :: p
    type(kdtree2), pointer :: mr
    type(kdtree2_result), pointer, contiguous, dimension(:) :: results
    integer(ik4) :: n1, n2, np, ni, nj, nf, i, j, i10
    real(rk8) :: dx, r2, rin
    logical :: lbil

    if ( any(shape(slat) /= shape(slon)) ) then
      call die('interp_create_g_g','Source shapes non conforming',1)
    end if
    if ( any(shape(tlat) /= shape(tlon)) ) then
      call die('interp_create_g_g','Target shapes non conforming',1)
    end if
    n1 = size(slat,2)
    n2 = size(slat,1)
    h_i%sshape(1) = n2
    h_i%sshape(2) = n1
    np = n1 * n2
    allocate(x(3,np))
    call ll2xyz(slat,slon,x)
    mr => kdtree2_create(x,sort=.true.,rearrange=.true.)
    rin = maxedis(slat,slon)
    dx = rin/erkm
    lbil = .true.
    if ( present(ds) ) then
      if ( rin < ds ) then
        dx = ds/erkm
        lbil = .false.
      end if
    end if
    if ( present(roi) ) then
      r2 = (dx*dx*roi*roi)
    else
      r2 = (dx*dx)
    end if
    ni = size(tlat,2)
    nj = size(tlat,1)
    h_i%tg%tshape = shape(tlat)
    allocate(h_i%tg%ft(nj,ni))
    write(stdout,'(a,f16.8,a)',advance='no') &
      ' Computing weights, R = ',sqrt(r2)*erkm, ' '
    i10 = ni/10
    do i = 1, ni
      if (mod(i,i10) == 0) write(stdout,'(a)',advance='no') '.'
      do j = 1, nj
        call ll2xyz(tlat(j,i),tlon(j,i),p)
        if ( lbil ) then
          np = minp
          allocate(results(minp))
          call kdtree2_n_nearest(mr,p,np,results)
        else
          np = kdtree2_r_count(mr,p,r2)
          if ( np < minp ) then
            h_i%tg%ft(j,i)%np = 0
            cycle
          else
            if ( np > maxp ) then
              np = maxp
              allocate(results(maxp))
              call kdtree2_n_nearest(mr,p,np,results)
            else
              allocate(results(np))
              call kdtree2_r_nearest(mr,p,r2,nf,np,results)
              np = nf
            end if
          end if
        end if
        call compwgt_genlin(np,n2,p,x,results,h_i%tg%ft(j,i)%wgt)
        !call compwgt_distwei(np,n2,results,h_i%tg%ft(j,i)%wgt)
        h_i%tg%ft(j,i)%np = np
        deallocate(results)
      end do
    end do
    deallocate(x)
    call kdtree2_destroy(mr)
    write(stdout,'(a)') ' Done.'
  end subroutine interp_create_g_g

  subroutine interp_create_g_ll(h_i,slat,slon,tlat,tlon,roi)
    implicit none
    type(h_interpolator), intent(out) :: h_i
    real(rkx), dimension(:,:), intent(in) :: slat
    real(rkx), dimension(:,:), intent(in) :: slon
    real(rkx), dimension(:), intent(in) :: tlat
    real(rkx), dimension(:), intent(in) :: tlon
    real(rkx), intent(in), optional :: roi
    real(rk8), dimension(:,:), allocatable :: x
    real(kdkind), dimension(3) :: p
    type(kdtree2), pointer :: mr
    type(kdtree2_result), pointer, contiguous, dimension(:) :: results
    integer(ik4) :: n1, n2, np, ni, nj, nf, i, j, i10
    real(rk8) :: dx, dx1, r2, rin, rin1
    logical :: lbil

    if ( any(shape(slat) /= shape(slon)) ) then
      call die('interp_create_g_g','Source shapes non conforming',1)
    end if
    n1 = size(slat,2)
    n2 = size(slat,1)
    h_i%sshape(1) = n2
    h_i%sshape(2) = n1
    np = n1 * n2
    allocate(x(3,np))
    call ll2xyz(slat,slon,x)
    mr => kdtree2_create(x,sort=.true.,rearrange=.true.)
    rin = maxedis(slat,slon)
    dx = rin/erkm
    rin1 = maxedis(tlat,tlon)
    dx1 = rin1/erkm
    lbil = .true.
    if ( rin < rin1 ) then
      dx = dx1
      lbil = .false.
    end if
    if ( present(roi) ) then
      r2 = (dx*dx*roi*roi)
    else
      r2 = (dx*dx)
    end if
    ni = size(tlat)
    nj = size(tlon)
    h_i%tg%tshape(1) = nj
    h_i%tg%tshape(2) = ni
    allocate(h_i%tg%ft(nj,ni))
    write(stdout,'(a,f16.8,a)',advance='no') &
      ' Computing weights, R = ',sqrt(r2)*erkm, ' '
    i10 = ni/10
    do i = 1, ni
      if (mod(i,i10) == 0) write(stdout,'(a)',advance='no') '.'
      do j = 1, nj
        call ll2xyz(tlat(i),tlon(j),p)
        if ( lbil ) then
          np = minp
          allocate(results(minp))
          call kdtree2_n_nearest(mr,p,np,results)
        else
          np = kdtree2_r_count(mr,p,r2)
          if ( np < minp ) then
            h_i%tg%ft(j,i)%np = 0
            cycle
          else
            if ( np > maxp ) then
              np = maxp
              allocate(results(maxp))
              call kdtree2_n_nearest(mr,p,np,results)
            else
              allocate(results(np))
              call kdtree2_r_nearest(mr,p,r2,nf,np,results)
              np = nf
            end if
          end if
        end if
        call compwgt_genlin(np,n2,p,x,results,h_i%tg%ft(j,i)%wgt)
        !call compwgt_distwei(np,n2,results,h_i%tg%ft(j,i)%wgt)
        h_i%tg%ft(j,i)%np = np
        deallocate(results)
      end do
    end do
    deallocate(x)
    call kdtree2_destroy(mr)
    write(stdout,'(a)') ' Done.'
  end subroutine interp_create_g_ll

  subroutine h_interpolator_destroy(h_i)
    implicit none
    type(h_interpolator), intent(inout) :: h_i
    integer :: ni, nj, j, i
    if ( associated(h_i%tg%ft) ) then
      nj = h_i%tg%tshape(1)
      ni = h_i%tg%tshape(2)
      do i = 1, ni
        do j = 1, nj
          if ( h_i%tg%ft(j,i)%np > 0 ) then
            deallocate(h_i%tg%ft(j,i)%wgt)
          end if
        end do
      end do
      deallocate(h_i%tg%ft)
    end if
  end subroutine h_interpolator_destroy

  subroutine interp_1d_nn(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:), intent(in) :: g
    real(rkx), dimension(:), intent(out) :: f
    integer(ik4) :: i, ni, n, si, gni
    real(rk8) :: gmax
    if ( size(g) /= h_i%sshape(1) ) then
      write(stderr,*) 'SOURCE SHAPE INTERP = ',h_i%sshape(1),' /= ',size(g)
      call die('interp_1d_nn','Non conforming shape for source',1)
    end if
    if ( size(f) /= h_i%tg%tshape(1) ) then
      write(stderr,*) 'TARGET SHAPE INTERP = ',h_i%tg%tshape(1),' /= ',size(f)
      call die('interp_1d_nn','Non conforming shape for target',1)
    end if
    ni = size(f)
    gni = size(g)
    do i = 1, ni
      gmax = -1.0_rkx
      f(i) = missl
      do n = 1, h_i%tg%ft(i,1)%np
        si = h_i%tg%ft(i,1)%wgt(n)%i
        if ( g(si) > missc .and. gmax < h_i%tg%ft(i,1)%wgt(n)%wgt ) then
          gmax = h_i%tg%ft(i,1)%wgt(n)%wgt
          f(i) = g(si)
        end if
      end do
    end do
  end subroutine interp_1d_nn

  subroutine interp_1d(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:), intent(in) :: g
    real(rkx), dimension(:), intent(out) :: f
    integer(ik4) :: i, ni, n, si, gni
    real(rkx) :: gmax, gmin
    real(rk8) :: gsum, gwgt
    if ( size(g) /= h_i%sshape(1) ) then
      write(stderr,*) 'SOURCE SHAPE INTERP = ',h_i%sshape(1),' /= ',size(g)
      call die('interp_1d','Non conforming shape for source',1)
    end if
    if ( size(f) /= h_i%tg%tshape(1) ) then
      write(stderr,*) 'TARGET SHAPE INTERP = ',h_i%tg%tshape(1),' /= ',size(f)
      call die('interp_1d','Non conforming shape for target',1)
    end if
    ni = size(f)
    gni = size(g)
    gmax = -1e20_rkx
    gmin = 1e20_rkx
    do i = 1, gni
      if ( g(i) > missc ) then
        if ( gmax < g(i) ) gmax = g(i)
        if ( gmin > g(i) ) gmin = g(i)
      end if
    end do
    do i = 1, ni
      gsum = d_zero
      gwgt = d_zero
      do n = 1, h_i%tg%ft(i,1)%np
        si = h_i%tg%ft(i,1)%wgt(n)%i
        if ( g(si) > missc ) then
          gsum = gsum + g(si) * h_i%tg%ft(i,1)%wgt(n)%wgt
          gwgt = gwgt + h_i%tg%ft(i,1)%wgt(n)%wgt
        end if
      end do
      if ( gwgt > d_zero ) then
        f(i) = real(gsum / gwgt,rkx)
      else
        f(i) = missl
      end if
    end do
    f = max(gmin,min(gmax,f))
  end subroutine interp_1d

  subroutine interp_2d(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:,:), intent(in) :: g
    real(rkx), dimension(:,:), intent(out) :: f
    integer(ik4) :: i, j, ni, nj, n, si, sj, gni, gnj
    real(rk8) :: gsum, gwgt
    real(rkx) :: gmax, gmin
    if ( any(shape(g) /= h_i%sshape) ) then
      write(stderr,*) 'SOURCE SHAPE INTERP = ',h_i%sshape,' /= ',shape(g)
      call die('interp_2d','Non conforming shape for source',1)
    end if
    if ( any(shape(f) /= h_i%tg%tshape) ) then
      write(stderr,*) 'TARGET SHAPE INTERP = ',h_i%tg%tshape,' /= ',shape(f)
      call die('interp_2d','Non conforming shape for target',1)
    end if
    nj = size(f,1)
    ni = size(f,2)
    gnj = size(g,1)
    gni = size(g,2)
    gmax = -1e20_rkx
    gmin = 1e20_rkx
    do i = 1, gni
      do j = 1, gnj
        if ( g(j,i) > missc ) then
          if ( gmax < g(j,i) ) gmax = g(j,i)
          if ( gmin > g(j,i) ) gmin = g(j,i)
        end if
      end do
    end do
    do i = 1, ni
      do j = 1, nj
        gsum = d_zero
        gwgt = d_zero
        do n = 1, h_i%tg%ft(j,i)%np
          si = h_i%tg%ft(j,i)%wgt(n)%i
          sj = h_i%tg%ft(j,i)%wgt(n)%j
          if ( g(sj,si) > missc ) then
            gsum = gsum + g(sj,si) * h_i%tg%ft(j,i)%wgt(n)%wgt
            gwgt = gwgt + h_i%tg%ft(j,i)%wgt(n)%wgt
          end if
        end do
        if ( gwgt > d_zero ) then
          f(j,i) = real(gsum / gwgt,rkx)
        else
          f(j,i) = missl
        end if
      end do
    end do
    call smtdsmt(f)
    do i = 1, ni
      do j = 1, nj
        if ( f(j,i) > missc ) then
          f(j,i) = max(gmin,min(gmax,f(j,i)))
        end if
      end do
    end do
  end subroutine interp_2d

  subroutine interp_2d_nn(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:,:), intent(in) :: g
    real(rkx), dimension(:,:), intent(out) :: f
    integer(ik4) :: i, j, ni, nj, n, si, sj, gni, gnj
    real(rk8) :: gmax
    if ( any(shape(g) /= h_i%sshape) ) then
      write(stderr,*) 'SOURCE SHAPE INTERP = ',h_i%sshape,' /= ',shape(g)
      call die('interp_2d_nn','Non conforming shape for source',1)
    end if
    if ( any(shape(f) /= h_i%tg%tshape) ) then
      write(stderr,*) 'TARGET SHAPE INTERP = ',h_i%tg%tshape,' /= ',shape(f)
      call die('interp_2d_nn','Non conforming shape for target',1)
    end if
    nj = size(f,1)
    ni = size(f,2)
    gnj = size(g,1)
    gni = size(g,2)
    do i = 1, ni
      do j = 1, nj
        gmax = -1.0_rkx
        f(j,i) = missl
        do n = 1, h_i%tg%ft(j,i)%np
          si = h_i%tg%ft(j,i)%wgt(n)%i
          sj = h_i%tg%ft(j,i)%wgt(n)%j
          if ( g(sj,si) > missc .and. gmax < h_i%tg%ft(j,i)%wgt(n)%wgt ) then
            gmax = h_i%tg%ft(j,i)%wgt(n)%wgt
            f(j,i) = g(sj,si)
          end if
        end do
      end do
    end do
  end subroutine interp_2d_nn

  subroutine interp_3d(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:,:,:), intent(in) :: g
    real(rkx), dimension(:,:,:), intent(out) :: f
    integer(ik4) :: n3, n
    n3 = size(g,3)
    if ( n3 /= size(f,3) ) then
      write(stderr,*) 'DIMENSION 3 g = ',size(g,3)
      write(stderr,*) 'DIMENSION 3 f = ',size(f,3)
      call die('interp_3d','Non conforming shapes',1)
    end if
!$OMP PARALLEL DO
    do n = 1, n3
      call interp_2d(h_i,g(:,:,n),f(:,:,n))
    end do
!$OMP END PARALLEL DO
  end subroutine interp_3d

  subroutine interp_4d(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:,:,:,:), intent(in) :: g
    real(rkx), dimension(:,:,:,:), intent(out) :: f
    integer(ik4) :: n4, n
    n4 = size(g,4)
    if ( n4 /= size(f,4) ) then
      write(stderr,*) 'DIMENSION 4 g = ',size(g,4)
      write(stderr,*) 'DIMENSION 4 f = ',size(f,4)
      call die('interp_4d','Non conforming shapes',1)
    end if
    do n = 1, n4
      call interp_3d(h_i,g(:,:,:,n),f(:,:,:,n))
    end do
  end subroutine interp_4d

  subroutine interp_5d(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:,:,:,:,:), intent(in) :: g
    real(rkx), dimension(:,:,:,:,:), intent(out) :: f
    integer(ik4) :: n5, n
    n5 = size(g,5)
    if ( n5 /= size(f,5) ) then
      write(stderr,*) 'DIMENSION 5 g = ',size(g,5)
      write(stderr,*) 'DIMENSION 5 f = ',size(f,5)
      call die('interp_5d','Non conforming shapes',1)
    end if
    do n = 1, n5
      call interp_4d(h_i,g(:,:,:,:,n),f(:,:,:,:,n))
    end do
  end subroutine interp_5d

  subroutine interp_class_2dr(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:,:), intent(in) :: g
    real(rkx), dimension(:,:), intent(out) :: f
    real(rk8) :: wgtm
    integer(ik4) :: i, j, ni, nj, n, si, sj
    if ( any(shape(g) /= h_i%sshape) ) then
      write(stderr,*) 'SOURCE SHAPE INTERP = ',h_i%sshape,' /= ',shape(g)
      call die('interp_class','Non conforming shape for source',1)
    end if
    if ( any(shape(f) /= h_i%tg%tshape) ) then
      write(stderr,*) 'TARGET SHAPE INTERP = ',h_i%tg%tshape,' /= ',shape(f)
      call die('interp_class','Non conforming shape for target',1)
    end if
    nj = size(f,1)
    ni = size(f,2)
    do i = 1, ni
      do j = 1, nj
        if ( h_i%tg%ft(j,i)%np > 0 ) then
          wgtm = h_i%tg%ft(j,i)%wgt(1)%wgt
          si = h_i%tg%ft(j,i)%wgt(1)%i
          sj = h_i%tg%ft(j,i)%wgt(1)%j
          f(j,i) = g(sj,si)
          do n = 2, h_i%tg%ft(j,i)%np
            if ( wgtm < h_i%tg%ft(j,i)%wgt(n)%wgt ) then
              wgtm = h_i%tg%ft(j,i)%wgt(n)%wgt
              si = h_i%tg%ft(j,i)%wgt(n)%i
              sj = h_i%tg%ft(j,i)%wgt(n)%j
              f(j,i) = g(sj,si)
            end if
          end do
        else
          f(j,i) = missl
        end if
      end do
    end do
  end subroutine interp_class_2dr

  subroutine interp_class_3dr(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:,:,:), intent(in) :: g
    real(rkx), dimension(:,:,:), intent(out) :: f
    integer(ik4) :: n3, n
    n3 = size(g,3)
    if ( n3 /= size(f,3) ) then
      write(stderr,*) 'DIMENSION 3 g = ',size(g,3)
      write(stderr,*) 'DIMENSION 3 f = ',size(f,3)
      call die('interp_class_3dr','Non conforming shapes',1)
    end if
!$OMP PARALLEL DO
    do n = 1, n3
      call interp_class_2dr(h_i,g(:,:,n),f(:,:,n))
    end do
!$OMP END PARALLEL DO
  end subroutine interp_class_3dr

  subroutine interp_class_4dr(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:,:,:,:), intent(in) :: g
    real(rkx), dimension(:,:,:,:), intent(out) :: f
    integer(ik4) :: n4, n
    n4 = size(g,4)
    if ( n4 /= size(f,4) ) then
      write(stderr,*) 'DIMENSION 4 g = ',size(g,4)
      write(stderr,*) 'DIMENSION 4 f = ',size(f,4)
      call die('interp_4d','Non conforming shapes',1)
    end if
    do n = 1, n4
      call interp_class_3dr(h_i,g(:,:,:,n),f(:,:,:,n))
    end do
  end subroutine interp_class_4dr

  subroutine interp_class_5dr(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    real(rkx), dimension(:,:,:,:,:), intent(in) :: g
    real(rkx), dimension(:,:,:,:,:), intent(out) :: f
    integer(ik4) :: n5, n
    n5 = size(g,5)
    if ( n5 /= size(f,5) ) then
      write(stderr,*) 'DIMENSION 5 g = ',size(g,5)
      write(stderr,*) 'DIMENSION 5 f = ',size(f,5)
      call die('interp_5d','Non conforming shapes',1)
    end if
    do n = 1, n5
      call interp_class_4dr(h_i,g(:,:,:,:,n),f(:,:,:,:,n))
    end do
  end subroutine interp_class_5dr

  subroutine interp_class_i(h_i,g,f)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    integer(ik4), dimension(:,:), intent(in) :: g
    integer(ik4), dimension(:,:), intent(out) :: f
    integer(ik4) :: i, j, ni, nj, n, si, sj, iv, nc, n1, n2
    integer(ik4), dimension(1) :: v
    real(rk8), dimension(:), allocatable :: gvals
    if ( any(shape(g) /= h_i%sshape) ) then
      write(stderr,*) 'SOURCE SHAPE INTERP = ',h_i%sshape,' /= ',shape(g)
      call die('interp_class','Non conforming shape for source',1)
    end if
    if ( any(shape(f) /= h_i%tg%tshape) ) then
      write(stderr,*) 'TARGET SHAPE INTERP = ',h_i%tg%tshape,' /= ',shape(f)
      call die('interp_class','Non conforming shape for target',1)
    end if
    n1 = minval(g)
    n2 = maxval(g)
    nc = n2 - n1 + 1
    if ( nc <= 0 ) then
      write(stderr,*) 'INCONSISTENCY IN CLASS NUMBER = ',nc
      call die('interp_class','CLASS NUMBER <= 0',1)
    end if
    nj = size(f,1)
    ni = size(f,2)
    allocate(gvals(n1:n2))
    do i = 1, ni
      do j = 1, nj
        gvals(:) = d_zero
        do n = 1, h_i%tg%ft(j,i)%np
          si = h_i%tg%ft(j,i)%wgt(n)%i
          sj = h_i%tg%ft(j,i)%wgt(n)%j
          iv = g(sj,si)
          gvals(iv) = gvals(iv) + 1
        end do
        v = maxloc(gvals) - 1 + n1
        f(j,i) = v(1)
      end do
    end do
    deallocate(gvals)
  end subroutine interp_class_i

  subroutine interp_class_ld(h_i,g,f,iw,pct)
    implicit none
    type(h_interpolator), intent(in) :: h_i
    integer(ik4), dimension(:,:), intent(in) :: g
    integer(ik4), intent(in) :: iw
    real(rkx), dimension(:,:), intent(out) :: f
    real(rkx), intent(in) :: pct
    integer(ik4) :: i, j, ni, nj, n, si, sj, iv, nc, n1, n2
    real(rk8) :: wgt
    integer(ik4), dimension(1) :: v
    integer(ik4), dimension(:), allocatable :: gvals
    if ( any(shape(g) /= h_i%sshape) ) then
      write(stderr,*) 'SOURCE SHAPE INTERP = ',h_i%sshape,' /= ',shape(g)
      call die('interp_class','Non conforming shape for source',1)
    end if
    if ( any(shape(f) /= h_i%tg%tshape) ) then
      write(stderr,*) 'TARGET SHAPE INTERP = ',h_i%tg%tshape,' /= ',shape(f)
      call die('interp_class','Non conforming shape for target',1)
    end if
    n1 = minval(g)
    n2 = maxval(g)
    nc = n2 - n1 + 1
    if ( nc <= 0 ) then
      write(stderr,*) 'INCONSISTENCY IN CLASS NUMBER = ',nc
      call die('interp_class','CLASS NUMBER <= 0',1)
    end if
    nj = size(f,1)
    ni = size(f,2)
    allocate(gvals(n1:n2))
    do i = 1, ni
      do j = 1, nj
        gvals(:) = 0
        do n = 1, h_i%tg%ft(j,i)%np
          si = h_i%tg%ft(j,i)%wgt(n)%i
          sj = h_i%tg%ft(j,i)%wgt(n)%j
          iv = g(sj,si)
          gvals(iv) = gvals(iv) + 1
        end do
        v = maxloc(gvals) - 1 + n1
        if ( v(1) == iw ) then
          wgt = real(gvals(v(1)),rkx)/real(sum(gvals),rkx)
          if ( wgt < pct/d_100 ) then
            gvals(iw) = 0
            v = maxloc(gvals) - 1 + n1
          end if
        end if
        f(j,i) = real(v(1),rkx)
      end do
    end do
    deallocate(gvals)
  end subroutine interp_class_ld

  subroutine smther(f)
    implicit none
    real(rkx), intent(inout), dimension(:,:) :: f
    integer(ik4) :: i1, i2, j1, j2
    integer(ik4) :: i, is, ie, j, js, je
    real(rkx), allocatable, dimension(:,:) :: tmp
    i1 = lbound(f,2)
    i2 = ubound(f,2)
    j1 = lbound(f,1)
    j2 = ubound(f,1)
    allocate(tmp(j1:j2,i1:i2))
    ie = i2-1
    je = j2-1
    is = i1+1
    js = j1+1
    tmp(:,:) = f(:,:)
    do i = is, ie
      do j = js, je
        if ( all(tmp(j-1:j+1,i-1:i+1) > missc) ) then
          f(j,i) = (tmp(j-1,i-1)+tmp(j-1,i)+tmp(j-1,i+1) + &
                    tmp(j+1,i-1)+tmp(j+1,i)+tmp(j+1,i+1) + &
                    tmp(j,i-1)+ tmp(j,i+1)+6.0_rkx*tmp(j,i))/14.0_rkx
        end if
      end do
    end do
    do i = is, ie
      if ( all(tmp(je:j2,i-1:i+1) > missc) ) then
        f(j2,i) = (tmp(je,i-1)+tmp(je,i)+tmp(je,i+1) + &
                   tmp(j2,i-1)+ tmp(j2,i+1)+3.0_rkx*tmp(j2,i))/8.0_rkx
      end if
      if ( all(tmp(j1:js,i-1:i+1) > missc) ) then
        f(j1,i) = (tmp(js,i-1)+tmp(js,i)+tmp(js,i+1) + &
                   tmp(j1,i-1)+ tmp(j1,i+1)+3.0_rkx*tmp(j1,i))/8.0_rkx
      end if
    end do
    do j = js, je
      if ( all(tmp(j-1:j+1,ie:i2) > missc) ) then
        f(j,i2) = (tmp(j-1,ie)+tmp(j,ie)+tmp(j+1,ie) + &
                   tmp(j+1,i2)+tmp(j-1,i2)+3.0_rkx* tmp(j,i2))/8.0_rkx
      end if
      if ( all(tmp(j-1:j+1,i1:is) > missc) ) then
        f(j,i1) = (tmp(j-1,is)+tmp(j,is)+tmp(j+1,is) + &
                   tmp(j+1,i1)+tmp(j-1,i1)+3.0_rkx*tmp(j,i1))/8.0_rkx
      end if
    end do
    if ( all(tmp(j1:js,i1:is) > missc) ) then
      f(j1,i1) = (2.0_rkx*tmp(j1,i1)+tmp(j1,is)+tmp(js,i1)+tmp(js,js))/5.0_rkx
    end if
    if ( all(tmp(je:j2,i1:is) > missc) ) then
      f(j2,i1) = (2.0_rkx*tmp(j2,i1)+tmp(j2,is)+tmp(je,i1)+tmp(je,is))/5.0_rkx
    end if
    if ( all(tmp(je:j2,ie:i2) > missc) ) then
      f(j2,i2) = (2.0_rkx*tmp(j2,i2)+tmp(j2,ie)+tmp(je,i2)+tmp(je,ie))/5.0_rkx
    end if
    if ( all(tmp(j1:js,ie:i2) > missc) ) then
      f(j1,i2) = (2.0_rkx*tmp(j1,i2)+tmp(j1,ie)+tmp(js,i2)+tmp(js,ie))/5.0_rkx
    end if
    deallocate(tmp)
  end subroutine smther

  subroutine smtdsmt(f)
    implicit none
    real(rkx), intent(inout), dimension(:,:) :: f
    real(rkx) :: aplus, asv, cell
    integer(ik4) :: i1, i2, j1, j2
    integer(ik4) :: i, is, ie, j, js, je, kp, np
    real(rkx), dimension(2) :: xnu
    integer(ik4), parameter :: npass = 4
    !
    ! purpose: spatially smooth data in f to dampen short
    ! wavelength components
    !
    i1 = lbound(f,2)
    i2 = ubound(f,2)
    j1 = lbound(f,1)
    j2 = ubound(f,1)
    ie = i2-1
    je = j2-1
    is = i1+1
    js = j1+1
    xnu(1) =  0.50_rkx
    xnu(2) = -0.52_rkx
    do np = 1, npass
      do kp = 1, 2
        ! first smooth in the ni direction
        do i = i1, i2
          asv = f(j1,i)
          do j = js, je
            cell = f(j,i)
            aplus = f(j+1,i)
            if ( asv > missc .and. aplus > missc .and. cell > missc ) then
              f(j,i) = cell + xnu(kp)*( (asv+aplus)/d_two - cell)
            end if
            asv = cell
          end do
        end do
        ! smooth in the nj direction
        do j = j1, j2
          asv = f(j,i1)
          do i = is, ie
            cell = f(j,i)
            aplus = f(j,i+1)
            if ( asv > missc .and. aplus > missc .and. cell > missc ) then
              f(j,i) = cell + xnu(kp)*((asv+aplus)/d_two - cell)
            end if
            asv = cell
          end do
        end do
      end do
    end do
  end subroutine smtdsmt

  subroutine compwgt_distwei(np,n2,r,w)
    implicit none
    integer(ik4), intent(inout) :: np
    integer(ik4), intent(in) :: n2
    type(kdtree2_result), pointer, contiguous, dimension(:), intent(in) :: r
    type(pwgt), dimension(:), pointer, contiguous, intent(inout) :: w
    integer(ik4) :: n
    real(rk8) :: rx, rmax

    rmax = r(1)%dis
    do n = 2, np
      if ( rmax < r(n)%dis ) rmax = r(n)%dis
    end do
    allocate(w(np))
    do n = 1, np
      if ( r(n)%dis < epsilon(rx) ) then
        np = 1
        deallocate(w)
        allocate(w(1))
        w(1)%i = (r(n)%idx-1)/n2 + 1
        w(1)%j = r(n)%idx - n2*(w(1)%i-1)
        w(1)%wgt = d_one
        return
      end if
      rx = sqrt(r(n)%dis/rmax)
      w(n)%i = (r(n)%idx-1)/n2 + 1
      w(n)%j = r(n)%idx - n2*(w(n)%i-1)
      w(n)%wgt = d_one/rx
    end do
  end subroutine compwgt_distwei

  subroutine compwgt_genlin_1d(np,p,xp,r,w)
    implicit none
    integer(ik4), intent(inout) :: np
    real(rk8), dimension(:,:), intent(in) :: xp
    real(rk8), dimension(3), intent(in) :: p
    type(kdtree2_result), pointer, contiguous, dimension(:), intent(in) :: r
    type(pwgt), dimension(:), pointer, contiguous, intent(inout) :: w
    real(rk8), dimension(3,np) :: v
    real(rk8), dimension(np) :: lambda
    real(rk8) :: rx
    integer(ik4) :: i, n

    ! Check perfect match
    do n = 1, np
      if ( abs(r(n)%dis) < epsilon(rx) ) then
        np = 1
        allocate(w(1))
        w(1)%i = r(n)%idx
        w(1)%j = 1
        w(1)%wgt = d_one
        return
      end if
    end do
    do n = 1, np
      do i = 1, 3
        v(i,n) = xp(i,r(n)%idx)
      end do
    end do
    allocate(w(np))
    call spherical_barycentric(np,p,v,lambda)
    do n = 1, np
      w(n)%i = r(n)%idx
      w(n)%j = 1
      w(n)%wgt = lambda(n)
    end do
  end subroutine compwgt_genlin_1d

  subroutine compwgt_genlin_2d(np,n2,p,xp,r,w)
    implicit none
    integer(ik4), intent(in) :: n2
    integer(ik4), intent(inout) :: np
    real(rk8), dimension(:,:), intent(in) :: xp
    real(rk8), dimension(3), intent(in) :: p
    type(kdtree2_result), pointer, contiguous, dimension(:), intent(in) :: r
    type(pwgt), dimension(:), pointer, contiguous, intent(inout) :: w
    real(rk8), dimension(3,np) :: v
    real(rk8), dimension(np) :: lambda
    real(rk8) :: rx
    integer(ik4) :: i, n

    ! Check perfect match
    do n = 1, np
      if ( abs(r(n)%dis) < epsilon(rx) ) then
        np = 1
        allocate(w(1))
        w(1)%i = (r(n)%idx-1)/n2 + 1
        w(1)%j = r(n)%idx - n2*(w(1)%i-1)
        w(1)%wgt = d_one
        return
      end if
    end do
    do n = 1, np
      do i = 1, 3
        v(i,n) = xp(i,r(n)%idx)
      end do
    end do
    allocate(w(np))
    call spherical_barycentric(np,p,v,lambda)
    do n = 1, np
      w(n)%i = (r(n)%idx-1)/n2 + 1
      w(n)%j = r(n)%idx - n2*(w(n)%i-1)
      w(n)%wgt = lambda(n)
    end do
  end subroutine compwgt_genlin_2d

end module mod_kdinterp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

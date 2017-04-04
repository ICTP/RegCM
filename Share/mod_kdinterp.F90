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

module mod_kdinterp

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_stdio
  use mod_message
  use mod_earth
  use mod_kdtree2

  private

  public :: h_interpolator
  public :: h_interpolator_create
  public :: h_interpolate_cont , h_interpolate_class
  public :: h_interpolator_destroy
  public :: h_missing_value

  interface h_interpolator_create
    module procedure interp_create_ll_g
    module procedure interp_create_g_g
    module procedure interp_create_g_ll
    module procedure interp_create_ll_ll
  end interface h_interpolator_create

  interface h_interpolate_cont
    module procedure interp_2d
    module procedure interp_3d
    module procedure interp_4d
    module procedure interp_5d
  end interface h_interpolate_cont

  interface h_interpolate_class
    module procedure interp_class_r
    module procedure interp_class_i
  end interface h_interpolate_class

  integer(ik4) , parameter :: minp = 5

  real(rkx) , parameter :: missl = -9999.0_rkx
  real(rkx) , parameter :: h_missing_value = missl
  real(rkx) , parameter :: missc = -9990.0_rkx
  real(rkx) , parameter :: mindis = 1.0e-6_rkx

  type pwgt
    integer(ik4) :: i , j
    real(rkx) :: wgt
  end type pwgt

  type ijwgt
    integer(ik4) :: np
    type(pwgt) , dimension(:) , allocatable :: wgt
  end type ijwgt

  type ftarget
    integer(ik4) , dimension(2) :: tshape
    type(ijwgt) , dimension(:,:) , allocatable :: ft
  end type ftarget

  type h_interpolator
    integer(ik4) , dimension(2) :: sshape
    type(ftarget) :: tg
    real(rkx) :: ds
  end type h_interpolator

  contains

  subroutine interp_create_ll_g(h_i,slat,slon,tlat,tlon,ds,roi)
    implicit none
    type(h_interpolator) , intent(out) :: h_i
    real(rkx) , dimension(:) , intent(in) :: slat
    real(rkx) , dimension(:) , intent(in) :: slon
    real(rkx) , dimension(:,:) , intent(in) :: tlat
    real(rkx) , dimension(:,:) , intent(in) :: tlon
    real(rkx) , intent(in) :: ds
    real(rkx) , intent(in) , optional :: roi
    real(rkx) , dimension(:,:) , allocatable :: x
    real(kdkind) , dimension(3) :: p
    type(kdtree2) , pointer :: mr
    type(kdtree2_result) , pointer , dimension(:) :: results
    integer(ik4) :: n1 , n2 , np , ni , nj , nf , i , j , n
    real(rkx) :: r2 , rx

    if ( any(shape(tlat) /= shape(tlon)) ) then
      call die('interp_create_ll_g','Target shapes non conforming',1)
    end if
    n1 = size(slat)
    n2 = size(slon)
    h_i%ds = ds
    h_i%sshape(1) = n2
    h_i%sshape(2) = n1
    np = n1 * n2
    allocate(x(3,np))
    call ll2xyz(slat,slon,x)
    mr => kdtree2_create(x,sort=.true.,rearrange=.true.)
    deallocate(x)
    if ( present(roi) ) then
      r2 = ds*ds*roi*roi
    else
      r2 = ds*ds
    end if
    nj = size(tlat,1)
    ni = size(tlat,2)
    h_i%tg%tshape = shape(tlat)
    allocate(h_i%tg%ft(nj,ni))
    do i = 1 , ni
      do j = 1 , nj
        call ll2xyz(tlat(j,i),tlon(j,i),p)
        np = kdtree2_r_count(mr,p,r2)
        if ( np < minp ) then
          np = minp
          allocate(results(minp))
          call kdtree2_n_nearest(mr,p,np,results)
          h_i%tg%ft(j,i)%np = np
        else
          allocate(results(np))
          call kdtree2_r_nearest(mr,p,r2,nf,np,results)
          h_i%tg%ft(j,i)%np = nf
        end if
        np = h_i%tg%ft(j,i)%np
        allocate(h_i%tg%ft(j,i)%wgt(h_i%tg%ft(j,i)%np))
        do n = 1 , np
          h_i%tg%ft(j,i)%wgt(n)%i = (results(n)%idx-1)/n2 + 1
          h_i%tg%ft(j,i)%wgt(n)%j = results(n)%idx - &
                                    n2*(h_i%tg%ft(j,i)%wgt(n)%i - 1)
          rx = max(sqrt(results(n)%dis)/ds,mindis)
          h_i%tg%ft(j,i)%wgt(n)%wgt = d_one/(rx*rx)
        end do
        deallocate(results)
      end do
    end do
    call kdtree2_destroy(mr)
  end subroutine interp_create_ll_g

  subroutine interp_create_ll_ll(h_i,slat,slon,tlat,tlon,ds,roi)
    implicit none
    type(h_interpolator) , intent(out) :: h_i
    real(rkx) , dimension(:) , intent(in) :: slat
    real(rkx) , dimension(:) , intent(in) :: slon
    real(rkx) , dimension(:) , intent(in) :: tlat
    real(rkx) , dimension(:) , intent(in) :: tlon
    real(rkx) , intent(in) :: ds
    real(rkx) , intent(in) , optional :: roi
    real(rkx) , dimension(:,:) , allocatable :: x
    real(kdkind) , dimension(3) :: p
    type(kdtree2) , pointer :: mr
    type(kdtree2_result) , pointer , dimension(:) :: results
    integer(ik4) :: n1 , n2 , np , ni , nj , nf , i , j , n
    real(rkx) :: r2 , rx

    n1 = size(slat)
    n2 = size(slon)
    h_i%ds = ds
    h_i%sshape(1) = n2
    h_i%sshape(2) = n1
    np = n1 * n2
    allocate(x(3,np))
    call ll2xyz(slat,slon,x)
    mr => kdtree2_create(x,sort=.true.,rearrange=.true.)
    deallocate(x)
    if ( present(roi) ) then
      r2 = ds*ds*roi*roi
    else
      r2 = ds*ds
    end if
    nj = size(tlat)
    ni = size(tlon)
    h_i%tg%tshape(1) = nj
    h_i%tg%tshape(2) = ni
    allocate(h_i%tg%ft(nj,ni))
    do i = 1 , ni
      do j = 1 , nj
        call ll2xyz(tlat(i),tlon(j),p)
        np = kdtree2_r_count(mr,p,r2)
        if ( np < minp ) then
          np = minp
          allocate(results(minp))
          call kdtree2_n_nearest(mr,p,np,results)
          h_i%tg%ft(j,i)%np = np
        else
          allocate(results(np))
          call kdtree2_r_nearest(mr,p,r2,nf,np,results)
          h_i%tg%ft(j,i)%np = nf
        end if
        np = h_i%tg%ft(j,i)%np
        allocate(h_i%tg%ft(j,i)%wgt(h_i%tg%ft(j,i)%np))
        do n = 1 , np
          h_i%tg%ft(j,i)%wgt(n)%i = (results(n)%idx-1)/n2 + 1
          h_i%tg%ft(j,i)%wgt(n)%j = results(n)%idx - &
                                    n2*(h_i%tg%ft(j,i)%wgt(n)%i-1)
          rx = max(sqrt(results(n)%dis)/ds,mindis)
          h_i%tg%ft(j,i)%wgt(n)%wgt = d_one/(rx*rx)
        end do
        deallocate(results)
      end do
    end do
    call kdtree2_destroy(mr)
  end subroutine interp_create_ll_ll

  subroutine interp_create_g_g(h_i,slat,slon,tlat,tlon,ds,roi)
    implicit none
    type(h_interpolator) , intent(out) :: h_i
    real(rkx) , dimension(:,:) , intent(in) :: slat
    real(rkx) , dimension(:,:) , intent(in) :: slon
    real(rkx) , dimension(:,:) , intent(in) :: tlat
    real(rkx) , dimension(:,:) , intent(in) :: tlon
    real(rkx) , intent(in) :: ds
    real(rkx) , intent(in) , optional :: roi
    real(rkx) , dimension(:,:) , allocatable :: x
    real(kdkind) , dimension(3) :: p
    type(kdtree2) , pointer :: mr
    type(kdtree2_result) , pointer , dimension(:) :: results
    integer(ik4) :: n1 , n2 , np , ni , nj , nf , i , j , n
    real(rkx) :: r2 , rx

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
    h_i%ds = ds
    np = n1 * n2
    allocate(x(3,np))
    call ll2xyz(slat,slon,x)
    mr => kdtree2_create(x,sort=.true.,rearrange=.true.)
    deallocate(x)
    if ( present(roi) ) then
      r2 = ds*ds*roi*roi
    else
      r2 = ds*ds
    end if
    ni = size(tlat,2)
    nj = size(tlat,1)
    h_i%tg%tshape = shape(tlat)
    allocate(h_i%tg%ft(nj,ni))
    do i = 1 , ni
      do j = 1 , nj
        call ll2xyz(tlat(j,i),tlon(j,i),p)
        np = kdtree2_r_count(mr,p,r2)
        if ( np < minp ) then
          np = minp
          allocate(results(minp))
          call kdtree2_n_nearest(mr,p,np,results)
          h_i%tg%ft(j,i)%np = np
        else
          allocate(results(np))
          call kdtree2_r_nearest(mr,p,r2,nf,np,results)
          h_i%tg%ft(j,i)%np = nf
        end if
        np = h_i%tg%ft(j,i)%np
        allocate(h_i%tg%ft(j,i)%wgt(h_i%tg%ft(j,i)%np))
        do n = 1 , np
          h_i%tg%ft(j,i)%wgt(n)%i = (results(n)%idx-1)/n2 + 1
          h_i%tg%ft(j,i)%wgt(n)%j = results(n)%idx - &
                                    n2*(h_i%tg%ft(j,i)%wgt(n)%i-1)
          rx = max(sqrt(results(n)%dis)/ds,mindis)
          h_i%tg%ft(j,i)%wgt(n)%wgt = d_one/(rx*rx)
        end do
        deallocate(results)
      end do
    end do
    call kdtree2_destroy(mr)
  end subroutine interp_create_g_g

  subroutine interp_create_g_ll(h_i,slat,slon,tlat,tlon,ds,roi)
    implicit none
    type(h_interpolator) , intent(out) :: h_i
    real(rkx) , dimension(:,:) , intent(in) :: slat
    real(rkx) , dimension(:,:) , intent(in) :: slon
    real(rkx) , dimension(:) , intent(in) :: tlat
    real(rkx) , dimension(:) , intent(in) :: tlon
    real(rkx) , intent(in) :: ds
    real(rkx) , intent(in) , optional :: roi
    real(rkx) , dimension(:,:) , allocatable :: x
    real(kdkind) , dimension(3) :: p
    type(kdtree2) , pointer :: mr
    type(kdtree2_result) , pointer , dimension(:) :: results
    integer(ik4) :: n1 , n2 , np , ni , nj , nf , i , j , n
    real(rkx) :: r2 , rx

    if ( any(shape(slat) /= shape(slon)) ) then
      call die('interp_create_g_g','Source shapes non conforming',1)
    end if
    n1 = size(slat,2)
    n2 = size(slat,1)
    h_i%sshape(1) = n2
    h_i%sshape(2) = n1
    h_i%ds = ds
    np = n1 * n2
    allocate(x(3,np))
    call ll2xyz(slat,slon,x)
    mr => kdtree2_create(x,sort=.true.,rearrange=.true.)
    deallocate(x)
    if ( present(roi) ) then
      r2 = ds*ds*roi*roi
    else
      r2 = ds*ds
    end if
    ni = size(tlat)
    nj = size(tlon)
    h_i%tg%tshape(1) = nj
    h_i%tg%tshape(2) = ni
    allocate(h_i%tg%ft(nj,ni))
    do i = 1 , ni
      do j = 1 , nj
        call ll2xyz(tlat(i),tlon(j),p)
        np = kdtree2_r_count(mr,p,r2)
        if ( np < minp ) then
          np = minp
          allocate(results(minp))
          call kdtree2_n_nearest(mr,p,np,results)
          h_i%tg%ft(j,i)%np = np
        else
          allocate(results(np))
          call kdtree2_r_nearest(mr,p,r2,nf,np,results)
          h_i%tg%ft(j,i)%np = nf
        end if
        np = h_i%tg%ft(j,i)%np
        allocate(h_i%tg%ft(j,i)%wgt(h_i%tg%ft(j,i)%np))
        do n = 1 , np
          h_i%tg%ft(j,i)%wgt(n)%i = (results(n)%idx-1)/n2 + 1
          h_i%tg%ft(j,i)%wgt(n)%j = results(n)%idx - &
                                    n2*(h_i%tg%ft(j,i)%wgt(n)%i-1)
          rx = max(sqrt(results(n)%dis)/ds,mindis)
          h_i%tg%ft(j,i)%wgt(n)%wgt = d_one/(rx*rx)
        end do
        deallocate(results)
      end do
    end do
    call kdtree2_destroy(mr)
  end subroutine interp_create_g_ll

  subroutine h_interpolator_destroy(h_i)
    implicit none
    type(h_interpolator) , intent(inout) :: h_i
    integer :: ni , nj , j , i
    if ( allocated(h_i%tg%ft) ) then
      nj = h_i%tg%tshape(1)
      ni = h_i%tg%tshape(2)
      do i = 1 , ni
        do j = 1 , nj
          deallocate(h_i%tg%ft(j,i)%wgt)
        end do
      end do
      deallocate(h_i%tg%ft)
    end if
  end subroutine h_interpolator_destroy

  subroutine interp_2d(h_i,g,f)
    implicit none
    type(h_interpolator) , intent(in) :: h_i
    real(rkx) , dimension(:,:) , intent(in) :: g
    real(rkx) , dimension(:,:) , intent(out) :: f
    integer(ik4) :: i , j , ni , nj , n , si , sj
    real(rkx) :: gsum , gwgt
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
    do i = 1 , ni
      do j = 1 , nj
        gsum = d_zero
        gwgt = d_zero
        do n = 1 , h_i%tg%ft(j,i)%np
          si = h_i%tg%ft(j,i)%wgt(n)%i
          sj = h_i%tg%ft(j,i)%wgt(n)%j
          if ( g(sj,si) > missc ) then
            gsum = gsum + g(sj,si) * h_i%tg%ft(j,i)%wgt(n)%wgt
            gwgt = gwgt + h_i%tg%ft(j,i)%wgt(n)%wgt
          end if
        end do
        if ( gwgt > d_zero ) then
          f(j,i) = gsum / gwgt
        else
          f(j,i) = missl
        end if
      end do
    end do
    call smtdsmt(f)
  end subroutine interp_2d

  subroutine interp_3d(h_i,g,f)
    implicit none
    type(h_interpolator) , intent(in) :: h_i
    real(rkx) , dimension(:,:,:) , intent(in) :: g
    real(rkx) , dimension(:,:,:) , intent(out) :: f
    integer(ik4) :: n3 , n
    n3 = size(g,3)
    if ( n3 /= size(f,3) ) then
      write(stderr,*) 'DIMENSION 3 g = ',size(g,3)
      write(stderr,*) 'DIMENSION 3 f = ',size(f,3)
      call die('interp_3d','Non conforming shapes',1)
    end if
    do n = 1 , n3
      call interp_2d(h_i,g(:,:,n),f(:,:,n))
    end do
  end subroutine interp_3d

  subroutine interp_4d(h_i,g,f)
    implicit none
    type(h_interpolator) , intent(in) :: h_i
    real(rkx) , dimension(:,:,:,:) , intent(in) :: g
    real(rkx) , dimension(:,:,:,:) , intent(out) :: f
    integer(ik4) :: n4 , n
    n4 = size(g,4)
    if ( n4 /= size(f,4) ) then
      write(stderr,*) 'DIMENSION 4 g = ',size(g,4)
      write(stderr,*) 'DIMENSION 4 f = ',size(f,4)
      call die('interp_4d','Non conforming shapes',1)
    end if
    do n = 1 , n4
      call interp_3d(h_i,g(:,:,:,n),f(:,:,:,n))
    end do
  end subroutine interp_4d

  subroutine interp_5d(h_i,g,f)
    implicit none
    type(h_interpolator) , intent(in) :: h_i
    real(rkx) , dimension(:,:,:,:,:) , intent(in) :: g
    real(rkx) , dimension(:,:,:,:,:) , intent(out) :: f
    integer(ik4) :: n5 , n
    n5 = size(g,5)
    if ( n5 /= size(f,5) ) then
      write(stderr,*) 'DIMENSION 5 g = ',size(g,5)
      write(stderr,*) 'DIMENSION 5 f = ',size(f,5)
      call die('interp_5d','Non conforming shapes',1)
    end if
    do n = 1 , n5
      call interp_4d(h_i,g(:,:,:,:,n),f(:,:,:,:,n))
    end do
  end subroutine interp_5d

  subroutine interp_class_r(h_i,g,f)
    implicit none
    type(h_interpolator) , intent(in) :: h_i
    real(rkx) , dimension(:,:) , intent(in) :: g
    real(rkx) , dimension(:,:) , intent(out) :: f
    integer(ik4) :: i , j , ni , nj , n , si , sj , iv , nc , n1 , n2
    integer(ik4) , dimension(1) :: v
    real(rkx) , dimension(:) , allocatable :: gvals
    if ( any(shape(g) /= h_i%sshape) ) then
      write(stderr,*) 'SOURCE SHAPE INTERP = ',h_i%sshape,' /= ',shape(g)
      call die('interp_class','Non conforming shape for source',1)
    end if
    if ( any(shape(f) /= h_i%tg%tshape) ) then
      write(stderr,*) 'TARGET SHAPE INTERP = ',h_i%tg%tshape,' /= ',shape(f)
      call die('interp_class','Non conforming shape for target',1)
    end if
    n1 = int(minval(g))
    n2 = int(maxval(g))
    nc = n2 - n1 + 1
    if ( nc <= 0 ) then
      write(stderr,*) 'INCONSISTENCY IN CLASS NUMBER = ',nc
      call die('interp_class','CLASS NUMBER <= 0',1)
    end if
    nj = size(f,1)
    ni = size(f,2)
    allocate(gvals(n1:n2))
    do i = 1 , ni
      do j = 1 , nj
        gvals(:) = d_zero
        do n = 1 , h_i%tg%ft(j,i)%np
          si = h_i%tg%ft(j,i)%wgt(n)%i
          sj = h_i%tg%ft(j,i)%wgt(n)%j
          iv = nint(g(sj,si))
          gvals(iv) = gvals(iv) + sqrt(h_i%tg%ft(j,i)%wgt(n)%wgt)
        end do
        v = maxloc(gvals)
        f(j,i) = real(v(1),rkx)
      end do
    end do
    deallocate(gvals)
  end subroutine interp_class_r

  subroutine interp_class_i(h_i,g,f)
    implicit none
    type(h_interpolator) , intent(in) :: h_i
    integer(ik4) , dimension(:,:) , intent(in) :: g
    integer(ik4) , dimension(:,:) , intent(out) :: f
    integer(ik4) :: i , j , ni , nj , n , si , sj , iv , nc , n1 , n2
    integer(ik4) , dimension(1) :: v
    real(rkx) , dimension(:) , allocatable :: gvals
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
    do i = 1 , ni
      do j = 1 , nj
        gvals(:) = d_zero
        do n = 1 , h_i%tg%ft(j,i)%np
          si = h_i%tg%ft(j,i)%wgt(n)%i
          sj = h_i%tg%ft(j,i)%wgt(n)%j
          iv = g(sj,si)
          gvals(iv) = gvals(iv) + sqrt(h_i%tg%ft(j,i)%wgt(n)%wgt)
        end do
        v = maxloc(gvals)
        f(j,i) = v(1)
      end do
    end do
    deallocate(gvals)
  end subroutine interp_class_i

  subroutine smtdsmt(f)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:) :: f
    real(rkx) :: aplus , asv , cell
    integer(ik4) :: i1 , i2 , j1 , j2
    integer(ik4) :: i , is , ie , j , js , je , k , kp
    real(rkx) , dimension(2) :: xnu
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
    do k = 1 , 4
      do kp = 1 , 2
        ! first smooth in the ni direction
        do i = is , ie
          asv = f(j1,i)
          do j = js , je
            cell = f(j,i)
            aplus = f(j+1,i)
            if ( asv > missc .and. aplus > missc .and. cell > missc ) then
              f(j,i) = cell + xnu(kp)*( (asv+aplus)/d_two - cell)
            end if
            asv = cell
          end do
        end do
        ! smooth in the nj direction
        do j = js , je
          asv = f(j,i1)
          do i = is , ie
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

end module mod_kdinterp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
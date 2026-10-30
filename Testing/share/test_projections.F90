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

! Unit test of Share/mod_projections.F90. For each projection: the domain
! centre (clat,clon) is at the middle of the grid, ijll and llij are
! inverse, the map factor matches the distance between grid points, the
! wind rotation gives the components along the grid axes, and the back
! rotation restores the input.

program test_projections
  use mod_intkinds
  use mod_realkinds
  use mod_projections
  use mod_memutil
  use share_test_utils
  implicit none
  real(rk8), parameter :: d2r = acos(-1.0_rk8)/180.0_rk8
  call memory_init
  call run('LAMCON NH','LAMCON',45.0_rk8,10.0_rk8, &
           0.0_rk8,0.0_rk8,30.0_rk8,60.0_rk8,61,51)
  call run('LAMCON SH','LAMCON',-35.0_rk8,-60.0_rk8, &
           0.0_rk8,0.0_rk8,-30.0_rk8,-60.0_rk8,61,51)
  call run('LAMCON tangent','LAMCON',45.0_rk8,10.0_rk8, &
           0.0_rk8,0.0_rk8,45.0_rk8,45.0_rk8,61,51)
  call run('POLSTR NH','POLSTR',75.0_rk8,20.0_rk8, &
           0.0_rk8,0.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('POLSTR SH','POLSTR',-75.0_rk8,120.0_rk8, &
           0.0_rk8,0.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('POLSTR SH2','POLSTR',-60.0_rk8,0.0_rk8, &
           0.0_rk8,0.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('NORMER equator','NORMER',0.0_rk8,20.0_rk8, &
           0.0_rk8,0.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('NORMER 30N','NORMER',30.0_rk8,20.0_rk8, &
           0.0_rk8,0.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('ROTMER pole=centre','ROTMER',45.0_rk8,10.0_rk8, &
           45.0_rk8,10.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('ROTMER pole/=centre','ROTMER',45.0_rk8,10.0_rk8, &
           40.0_rk8,5.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('ROTMER SH','ROTMER',-30.0_rk8,-60.0_rk8, &
           -30.0_rk8,-60.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('ROTLLR AUS','ROTLLR',-25.0_rk8,135.0_rk8, &
           60.31_rk8,141.38_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('ROTLLR EUR','ROTLLR',50.0_rk8,10.0_rk8, &
           39.25_rk8,-162.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('ROTLLR SA','ROTLLR',-20.0_rk8,-60.0_rk8, &
           70.0_rk8,120.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('ROTLLR dateline','ROTLLR',60.0_rk8,178.0_rk8, &
           30.0_rk8,-2.0_rk8,0.0_rk8,0.0_rk8,61,51)
  call run('ROTMER tall','ROTMER',30.0_rk8,10.0_rk8, &
           30.0_rk8,10.0_rk8,0.0_rk8,0.0_rk8,61,201)
  call run('NORMER tall','NORMER',30.0_rk8,10.0_rk8, &
           0.0_rk8,0.0_rk8,0.0_rk8,0.0_rk8,61,201)
  call run('LAMCON even grid','LAMCON',45.0_rk8,10.0_rk8, &
           0.0_rk8,0.0_rk8,30.0_rk8,60.0_rk8,60,50)
  call run('ROTLLR even grid','ROTLLR',50.0_rk8,10.0_rk8, &
           39.25_rk8,-162.0_rk8,0.0_rk8,0.0_rk8,60,50)
  call memory_destroy
  call test_end('test_projections')

  contains

  subroutine run(label,code,clat,clon,plat,plon,t1,t2,nx,ny)
    implicit none
    character(len=*), intent(in) :: label, code
    real(rk8), intent(in) :: clat, clon, plat, plon, t1, t2
    integer(ik4), intent(in) :: nx, ny
    real(rk8), parameter :: ds = 25000.0_rk8
    type(anyprojparams) :: pp
    type(regcm_projection) :: pj
    real(rkx) :: ri, rj, la1, lo1, la2, lo2
    real(rkx), pointer, contiguous, dimension(:,:) :: xl => null()
    real(rkx), pointer, contiguous, dimension(:,:) :: xo => null()
    real(rkx), pointer, contiguous, dimension(:,:) :: xm => null()
    real(rkx), pointer, contiguous, dimension(:,:) :: u => null()
    real(rkx), pointer, contiguous, dimension(:,:) :: v => null()
    real(rkx), pointer, contiguous, dimension(:,:) :: u0 => null()
    real(rkx), pointer, contiguous, dimension(:,:) :: v0 => null()
    real(rk8), dimension(2) :: ex, ey
    real(rk8) :: errij, errm, errw, errb, mnum, un, vn
    integer(ik4) :: i, j, k

    pp%pcode = code
    pp%ds = ds
    pp%clat = clat
    pp%clon = clon
    pp%plat = plat
    pp%plon = plon
    pp%trlat1 = t1
    pp%trlat2 = t2
    pp%nlat = ny
    pp%nlon = nx
    pp%staggerx = .false.
    pp%staggery = .false.
    pp%rotparam = .true.
    call pj%initialize(pp)

    call pj%llij(real(clat,rkx),real(clon,rkx),ri,rj)
    call check(abs(ri-real(nx+1,rkx)*0.5_rkx) < 1.0e-9_rkx .and. &
               abs(rj-real(ny+1,rkx)*0.5_rkx) < 1.0e-9_rkx, &
               label//': centre at the middle of the grid')

    call getmem(xl,1,nx,1,ny,'xl')
    call getmem(xo,1,nx,1,ny,'xo')
    call getmem(xm,1,nx,1,ny,'xm')
    call getmem(u,1,nx,1,ny,'u')
    call getmem(v,1,nx,1,ny,'v')
    call getmem(u0,1,nx,1,ny,'u0')
    call getmem(v0,1,nx,1,ny,'v0')
    errij = 0.0_rk8
    do j = 1, ny
      do i = 1, nx
        call pj%ijll(real(i,rkx),real(j,rkx),xl(i,j),xo(i,j))
        call pj%llij(xl(i,j),xo(i,j),ri,rj)
        errij = max(errij,abs(ri-real(i,rkx)),abs(rj-real(j,rkx)))
      end do
    end do
    call check(errij < 1.0e-9_rk8,label//': ijll and llij are inverse')

    ! The map factor is the grid distance over the true distance between
    ! the two points half a grid step on either side.
    call pj%mapfac(xl,xo,xm)
    errm = 0.0_rk8
    do j = 2, ny-1
      do i = 2, nx-1
        call pj%ijll(real(i,rkx)-0.5_rkx,real(j,rkx),la1,lo1)
        call pj%ijll(real(i,rkx)+0.5_rkx,real(j,rkx),la2,lo2)
        mnum = ds/gcd(la1,lo1,la2,lo2)
        errm = max(errm,abs(xm(i,j)/mnum-1.0_rk8))
      end do
    end do
    ! The reference distance uses a 6371 km radius, the model 6371.229 km.
    call check(errm < 1.0e-4_rk8,label//': map factor')

    ! An eastward (k = 1) and a northward (k = 2) unit wind: the rotated
    ! components are the projections on the grid axes.
    errw = 0.0_rk8
    errb = 0.0_rk8
    do k = 1, 2
      u = merge(1.0_rkx,0.0_rkx,k == 1)
      v = merge(0.0_rkx,1.0_rkx,k == 1)
      u0 = u
      v0 = v
      call pj%wind2_rotate(u,v)
      do j = 2, ny-1
        do i = 2, nx-1
          ex = axis(pj,real(i,rkx),real(j,rkx),1)
          ey = axis(pj,real(i,rkx),real(j,rkx),2)
          un = u0(i,j)*ex(1) + v0(i,j)*ex(2)
          vn = u0(i,j)*ey(1) + v0(i,j)*ey(2)
          errw = max(errw,abs(u(i,j)-un),abs(v(i,j)-vn))
        end do
      end do
      call pj%wind2_antirotate(u,v)
      errb = max(errb,maxval(abs(u-u0)),maxval(abs(v-v0)))
    end do
    call check(errw < 1.0e-6_rk8,label//': wind rotation')
    call check(errb < 1.0e-12_rk8,label//': wind back rotation')

    call relmem(xl)
    call relmem(xo)
    call relmem(xm)
    call relmem(u)
    call relmem(v)
    call relmem(u0)
    call relmem(v0)
    call pj%destruct( )
  end subroutine run

  ! Unit vector (east, north) along the grid i (n = 1) or j (n = 2) axis.
  function axis(pj,ri,rj,n) result(e)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: ri, rj
    integer(ik4), intent(in) :: n
    real(rk8), dimension(2) :: e
    real(rkx), parameter :: h = 0.01_rkx
    real(rkx) :: la1, lo1, la2, lo2
    if ( n == 1 ) then
      call pj%ijll(ri-h,rj,la1,lo1)
      call pj%ijll(ri+h,rj,la2,lo2)
    else
      call pj%ijll(ri,rj-h,la1,lo1)
      call pj%ijll(ri,rj+h,la2,lo2)
    end if
    e(1) = dwrap(real(lo2-lo1,rk8))*cos(0.5_rk8*real(la1+la2,rk8)*d2r)
    e(2) = real(la2-la1,rk8)
    e = e/sqrt(sum(e*e))
  end function axis

  real(rk8) function dwrap(x)
    implicit none
    real(rk8), intent(in) :: x
    dwrap = modulo(x+180.0_rk8,360.0_rk8) - 180.0_rk8
  end function dwrap

  ! Great circle distance (m) on a sphere of 6371 km radius.
  real(rk8) function gcd(la1,lo1,la2,lo2)
    implicit none
    real(rkx), intent(in) :: la1, lo1, la2, lo2
    real(rk8) :: a
    a = sin(real(la2-la1,rk8)*d2r*0.5_rk8)**2 + &
        cos(real(la1,rk8)*d2r)*cos(real(la2,rk8)*d2r) * &
        sin(dwrap(real(lo2-lo1,rk8))*d2r*0.5_rk8)**2
    gcd = 2.0_rk8*6.371e6_rk8*asin(sqrt(a))
  end function gcd

end program test_projections

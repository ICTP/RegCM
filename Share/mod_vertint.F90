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

module mod_vertint

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_constants
  use mod_message
  use mod_stdatm
  use mod_interp , only : interp1d

  implicit none

  private

  real(rkx) , parameter :: missl = -9999.0_rkx

  real(rkx) , parameter :: rgas2 = rgas/d_two
  ! lrate is defined as a positive constant.
  real(rkx) , parameter :: rglrog = rgas*lrate*regrav
  real(rkx) , parameter :: b1 = -egrav/lrate

  interface intlin
    module procedure intlin_double
    module procedure intlin_single
    module procedure intlin_o_double
    module procedure intlin_o_single
  end interface intlin

  interface intlinz
    module procedure intlin_z_o_single
    module procedure intlin_z_o_double
  end interface intlinz

  interface intlog
    module procedure intlog_double
    module procedure intlog_single
    module procedure intlog_o_double
    module procedure intlog_o_single
  end interface intlog

  interface intlinreg
    module procedure intlinreg_p
    module procedure intlinreg_z
  end interface intlinreg

  interface intzps
    module procedure intzps1
    module procedure intzps2
  end interface intzps

  public :: intlin , intgtb , intlog , intlinz
  public :: intpsn , intv0 , intv1 , intvp , intv2 , intv3
  public :: intlinreg , intlinprof
  public :: intz1 , intp1 , intp3 , intz3 , intzps

  contains

  ! The subroutine vertically interpolates from regular grid of pressure
  ! levels p to a non regular pressure levels grid p3d. Surface pressure
  ! is used to create a dimensionless vertical coordinate.
  ! Input:
  !   f(im,jm,kp)   - the field to be interpolated from
  !   p(kp)         - the pressure levels to interpolate from
  !   p3d(im,jm,km) - the pressure levels to interpolate to
  !   ps(im,jm)     - the surface pressure
  ! Output:
  !   fp(im,jm,km)  - the field interpolate on regular pressure grid p
  !
  subroutine intlinreg_p(fp,f,ps,p,im1,im2,jm1,jm2,kp,p3d,km)
    implicit none
    integer(ik4) , intent(in) :: im1 , im2 , jm1 , jm2 , km , kp
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: f
    real(rkx) , pointer , dimension(:,:) , intent(in) :: ps
    real(rkx) , pointer , dimension(:) , intent(in) :: p
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: p3d
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rkx) , dimension(kp) :: sig
    real(rkx) :: sigp , w1 , wp

    fp(:,:,:) = missl
    if ( p(1) > p(kp) ) then
      do j = jm1 , jm2
        do i = im1 , im2
          if ( ps(i,j) < 1.0_rkx ) cycle
          do k = 1 , kp
            sig(k) = (p(k)-p(kp))/(ps(i,j)-p(kp))
          end do
          do n = 1 , km
            sigp = (p3d(i,j,n)-p(kp))/(ps(i,j)-p(kp))
            if ( sigp <= sig(kp) ) then
              fp(i,j,n) = f(i,j,kp)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            end if
            kx = 2
            do k = 2 , kp
              kx = k
              if ( sig(k) < sigp ) exit
            end do
            knx = kx - 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    else
      do j = jm1 , jm2
        do i = im1 , im2
          if ( ps(i,j) < 1.0_rkx ) cycle
          do k = 1 , kp
            sig(k) = (p(k)-p(1))/(ps(i,j)-p(1))
          end do
          do n = 1 , km
            sigp = (p3d(i,j,n)-p(1))/(ps(i,j)-p(1))
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            else if ( sigp >= sig(kp) ) then
              fp(i,j,n) = f(i,j,kp)
              cycle
            end if
            kx = 2
            do k = 2 , kp
              kx = k
              if ( sig(k) > sigp ) exit
            end do
            knx = kx - 1
            wp = (sig(kx)-sigp)/(sig(kx)-sig(knx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlinreg_p

  ! The subroutine vertically interpolates from regular grid of height
  ! levels z to a non regular height levels grid hz.
  ! Input:
  !   f(km,jm,im)   - the field to be interpolated from
  !   z(km)         - the height levels the field is defined at
  !   hz(kz,jm,im)  - the height levels to interpolate to
  ! Output:
  !   fz(kz,jm,im)  - the field interpolate on regular pressure grid p
  !
  subroutine intlinreg_z(fz,f,z,i1,i2,j1,j2,km,hz,kz)
    implicit none
    integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , km , kz
    real(rkx) , dimension(i1:i2,j1:j2,km) , intent(in) ::  f
    real(rkx) , dimension(i1:i2,j1:j2,kz) , intent(in) ::  hz
    real(rkx) , dimension(km) , intent(in) :: z
    real(rkx) , dimension(i1:i2,j1:j2,kz) , intent(out) :: fz
    integer(ik4) :: i , j , k , kx , knx , n
    real(rkx) :: w1 , w2
    if ( z(1) < z(km) ) then
      do j = j1 , j2
        do i = i1 , i2
          do n = 1 , kz
            if ( z(km) >= hz(i,j,n) ) then
              fz(i,j,n) = f(i,j,km)
              cycle
            else if ( z(1) <= hz(i,j,n) ) then
              fz(i,j,n) = f(i,j,1)
              cycle
            end if
            do k = 1 , km
              kx = k
              if ( hz(i,j,n) > z(k) ) exit
            end do
            knx = kx - 1
            w1 = (hz(i,j,n)-z(knx))/(z(kx)-z(knx))
            w2 = 1.0 - w1
            fz(i,j,n) = w1*f(i,j,kx) + w2*f(i,j,knx)
          end do
        end do
      end do
    else
      do j = j1 , j2
        do i = i1 , i2
          do n = 1 , kz
            if ( z(1) >= hz(i,j,n) ) then
              fz(i,j,n) = f(i,j,1)
              cycle
            else if ( z(km) <= hz(i,j,km) ) then
              fz(i,j,n) = f(i,j,km)
              cycle
            end if
            do k = km , 1 , -1
              kx = k
              if ( hz(i,j,n) > z(k) ) exit
            end do
            knx = kx + 1
            w1 = (hz(i,j,n)-z(knx))/(z(kx)-z(knx))
            w2 = 1.0 - w1
            fz(i,j,n) = w1*f(i,j,kx) + w2*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlinreg_z

  ! The subroutine vertically interpolates a vertical profile f on p
  ! levels onto a grid with pressure levels p3d. Surface pressure
  ! is used to create a dimensionless vertical coordinate.
  ! Input:
  !   f(kp)         - the field to be interpolated from
  !   p(kp)         - the pressure levels the field is defined at
  !   p3d(im,jm,km) - the pressure levels to interpolate to
  !   ps(im,jm)     - the surface pressure
  ! Output:
  !   fp(im,jm,km)  - the field interpolate on regular pressure grid p3d
  !
  subroutine intlinprof(fp,f,ps,p,im1,im2,jm1,jm2,kp,p3d,km)
    implicit none
    integer(ik4) , intent(in) :: im1 , im2 , jm1 , jm2 , km , kp
    real(rkx) , pointer , dimension(:) , intent(in) :: f
    real(rkx) , pointer , dimension(:,:) , intent(in) :: ps
    real(rkx) , pointer , dimension(:) , intent(in) :: p
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: p3d
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rkx) , dimension(kp) :: sig
    real(rkx) :: sigp , w1 , wp

    fp = missl

    if ( p(1) > p(kp) ) then
      do j = jm1 , jm2
        do i = im1 , im2
          if ( ps(i,j) < 1.0_rkx ) cycle
          do k = 1 , kp
            sig(k) = (p(k)-p(kp))/(ps(i,j)-p(kp))
          end do
          do n = 1 , km
            sigp = (p3d(i,j,n)-p(kp))/(ps(i,j)-p(kp))
            if ( sigp <= sig(kp) ) then
              fp(i,j,n) = f(kp)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(1)
              cycle
            end if
            kx = 2
            do k = 2 , kp
              kx = k
              if ( sig(k) < sigp ) exit
            end do
            knx = kx - 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(kx) + wp*f(knx)
          end do
        end do
      end do
    else
      do j = jm1 , jm2
        do i = im1 , im2
          if ( ps(i,j) < 1.0_rkx ) cycle
          do k = 1 , kp
            sig(k) = (p(k)-p(1))/(ps(i,j)-p(1))
          end do
          do n = 1 , km
            sigp = (p3d(i,j,n)-p(1))/(ps(i,j)-p(1))
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(1)
              cycle
            else if ( sigp >= sig(kp) ) then
              fp(i,j,n) = f(kp)
              cycle
            end if
            kx = 2
            do k = 2 , kp
              kx = k
              if ( sig(k) > sigp ) exit
            end do
            knx = kx - 1
            wp = (sig(kx)-sigp)/(sig(kx)-sig(knx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(kx) + wp*f(knx)
          end do
        end do
      end do
    end if
  end subroutine intlinprof

  ! The subroutine vertically interpolates on regular grid of pressure
  ! levels p from a non regular pressure levels grid p3d.
  ! The single and double precision version are one after the other.
  ! Input:
  !   f(im,jm,km)   - the field to be interpolated from
  !   p3d(im,jm,km) - the pressure levels the field is defined at
  !   p(kp)         - the pressure levels to interpolate to
  ! Output:
  !   fp(im,jm,kp)  - the field interpolate on regular pressure grid p
  !
  subroutine intlin_double(fp,f,p3d,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , dimension(im,jm,km) , intent(in) :: f , p3d
    real(rk8) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk8) , dimension(kp) , intent(in) :: p
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) , dimension(km) :: sig
    real(rk8) :: sigp , w1 , wp , tp , bp
    if ( p3d(1,1,1) > p3d(1,1,km) ) then
      do j = 1 , jm
        do i = 1 , im
          tp = p3d(i,j,km)
          bp = p3d(i,j,1)-tp
          do k = 1 , km
            sig(k) = (p3d(i,j,k)-tp)/bp
          end do
          do n = 1 , kp
            sigp = (p(n)-tp)/bp
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
            else
              kx = 2
              do k = 2 , km
                kx = k
                if ( sig(k) < sigp ) exit
              end do
              knx = kx - 1
              wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          tp = p3d(i,j,1)
          bp = p3d(i,j,km)-tp
          do k = 1 , km
            sig(k) = (p3d(i,j,k)-tp)/bp
          end do
          do n = 1 , kp
            sigp = (p(n)-tp)/bp
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
            else if ( sigp >= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
            else
              kx = 2
              do k = 2 , km
                kx = k
                if ( sig(k) > sigp ) exit
              end do
              knx = kx - 1
              wp = (sig(kx)-sigp)/(sig(kx)-sig(knx))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    end if
  end subroutine intlin_double

  subroutine intlin_single(fp,f,p3d,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk4) , dimension(im,jm,km) , intent(in) :: f , p3d
    real(rk4) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk4) , dimension(kp) , intent(in) :: p
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk4) , dimension(km) :: sig
    real(rk4) :: sigp , w1 , wp , tp , bp
    if ( p3d(1,1,1) > p3d(1,1,km) ) then
      do j = 1 , jm
        do i = 1 , im
          tp = p3d(i,j,km)
          bp = p3d(i,j,1) - tp
          do k = 1 , km
            sig(k) = (p3d(i,j,k)-tp)/bp
          end do
          do n = 1 , kp
            sigp = (p(n)-tp)/bp
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
            else
              kx = 2
              do k = 2 , km
                kx = k
                if ( sig(k) < sigp ) exit
              end do
              knx = kx - 1
              wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          tp = p3d(i,j,1)
          bp = p3d(i,j,km) - tp
          do k = 1 , km
            sig(k) = (p3d(i,j,k)-tp)/bp
          end do
          do n = 1 , kp
            sigp = (p(n)-tp)/bp
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
            else if ( sigp >= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
            else
              kx = 2
              do k = 2 , km
                kx = k
                if ( sig(k) > sigp ) exit
              end do
              knx = kx - 1
              wp = (sig(kx)-sigp)/(sig(kx)-sig(knx))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    end if
  end subroutine intlin_single

  ! The subroutine vertically interpolates on regular grid of pressure
  ! levels p from hydrostatic vertical levels at sigma oordinates.
  ! The single and double precision version are one after the other.
  ! Input:
  !   f(im,jm,km)   - the field to be interpolated from
  !   ps(im,jm)     - the surface pressure
  !   sig(km)       - the sigma coordinate of input levels
  !   ptop          - the model top pressure
  !   p(kp)         - the pressure levels to interpolate to
  ! Output:
  !   fp(im,jm,kp)  - the field interpolate on regular pressure grid p
  !
  subroutine intlin_o_double(fp,f,ps,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , intent(in) :: ptop
    real(rk8) , dimension(im,jm,km) , intent(in) :: f
    real(rk8) , dimension(kp) , intent(in) :: p
    real(rk8) , dimension(im,jm) , intent(in) :: ps
    real(rk8) , dimension(km) , intent(in) :: sig
    real(rk8) , dimension(im,jm,kp) , intent(out) :: fp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) :: sigp , w1 , wp

    if ( sig(1) > sig(2) ) then
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kp
            sigp = (p(n)-ptop)/(ps(i,j)-ptop)
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            end if
            do k = 2 , km
              kx = k
              if ( sig(k) < sigp ) exit
            end do
            knx = kx - 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kp
            sigp = (p(n)-ptop)/(ps(i,j)-ptop)
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            else if ( sigp >= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            end if
            do k = 2 , km
              kx = k
              if ( sig(k) > sigp ) exit
            end do
            knx = kx - 1
            wp = (sig(kx)-sigp)/(sig(kx)-sig(knx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlin_o_double

  subroutine intlin_o_single(fp,f,ps,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rkx) , intent(in) :: ptop
    real(rk4) , dimension(im,jm,km) , intent(in) :: f
    real(rk4) , dimension(kp) , intent(in) :: p
    real(rk4) , dimension(im,jm) , intent(in) :: ps
    real(rk4) , dimension(km) , intent(in) :: sig
    real(rk4) , dimension(im,jm,kp) , intent(out) :: fp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk4) :: sigp , w1 , wp , pt
    pt = real(ptop)
    if ( sig(1) > sig(2) ) then
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kp
            sigp = (p(n)-pt)/(ps(i,j)-pt)
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            end if
            do k = 2 , km
              kx = k
              if ( sig(k) < sigp ) exit
            end do
            knx = kx - 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = 1.0 - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kp
            sigp = (p(n)-pt)/(ps(i,j)-pt)
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            else if ( sigp >= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            end if
            do k = 2 , km
              kx = k
              if ( sig(k) > sigp ) exit
            end do
            knx = kx - 1
            wp = (sig(kx)-sigp)/(sig(kx)-sig(knx))
            w1 = 1.0 - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlin_o_single

  ! The subroutine vertically interpolates from non regular grid of height
  ! levels hz to a regular height levels grid z.
  ! Input:
  !   f(km,jm,im)   - the field to be interpolated from
  !   hz(km,jm,im)  - the height levels the field is defined at
  !   z(kz)         - the height levels to interpolate to
  ! Output:
  !   fz(kz,jm,im)  - the field interpolate on regular pressure grid p
  !
  subroutine intlin_z_o_single(fz,f,hz,im,jm,km,z,kz)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kz
    real(rk4) , dimension(im,jm,km) , intent(in) :: f , hz
    real(rk4) , dimension(kz) , intent(in) :: z
    real(rk4) , dimension(im,jm,kz) , intent(out) :: fz
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk4) :: w1 , wz
    if ( hz(1,1,1) < hz(1,1,km) ) then
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kz
            if ( z(n) >= hz(i,j,km) ) then
              fz(i,j,n) = f(i,j,km)
              cycle
            else if ( z(n) <= hz(i,j,1) ) then
              fz(i,j,n) = f(i,j,1)
              cycle
            end if
            kx = 2
            do k = 2 , km
              kx = k
              if ( z(n) < hz(i,j,k) ) exit
            end do
            knx = kx - 1
            wz = (hz(i,j,kx)-z(n))/(hz(i,j,kx)-hz(i,j,knx))
            w1 = 1.0 - wz
            fz(i,j,n) = w1*f(i,j,kx) + wz*f(i,j,knx)
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kz
            if ( z(n) >= hz(i,j,1) ) then
              fz(i,j,n) = f(i,j,1)
              cycle
            else if ( z(n) <= hz(i,j,km) ) then
              fz(i,j,n) = f(i,j,km)
              cycle
            end if
            kx = 2
            do k = 2 , km
              kx = k
              if ( hz(i,j,k) < z(n) ) exit
            end do
            knx = kx - 1
            wz = (z(n)-hz(i,j,kx))/(hz(i,j,knx)-hz(i,j,kx))
            w1 = 1.0 - wz
            fz(i,j,n) = w1*f(i,j,kx) + wz*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlin_z_o_single

  subroutine intlin_z_o_double(fz,f,hz,im,jm,km,z,kz)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kz
    real(rk8) , dimension(im,jm,km) , intent(in) :: f , hz
    real(rk8) , dimension(kz) , intent(in) :: z
    real(rk8) , dimension(im,jm,kz) , intent(out) :: fz
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) :: w1 , wz
    if ( hz(1,1,1) < hz(1,1,km) ) then
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kz
            if ( z(n) >= hz(i,j,km) ) then
              fz(i,j,n) = f(i,j,km)
              cycle
            else if ( z(n) <= hz(i,j,1) ) then
              fz(i,j,n) = f(i,j,1)
              cycle
            end if
            kx = 2
            do k = 2 , km
              kx = k
              if ( z(n) < hz(i,j,k) ) exit
            end do
            knx = kx - 1
            wz = (hz(i,j,kx)-z(n))/(hz(i,j,kx)-hz(i,j,knx))
            w1 = 1.0 - wz
            fz(i,j,n) = w1*f(i,j,kx) + wz*f(i,j,knx)
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kz
            if ( z(n) >= hz(i,j,1) ) then
              fz(i,j,n) = f(i,j,1)
              cycle
            else if ( z(n) <= hz(i,j,km) ) then
              fz(i,j,n) = f(i,j,km)
              cycle
            end if
            kx = 2
            do k = 2 , km
              kx = k
              if ( hz(i,j,k) < z(n) ) exit
            end do
            knx = kx - 1
            wz = (z(n)-hz(i,j,kx))/(hz(i,j,knx)-hz(i,j,kx))
            w1 = 1.0 - wz
            fz(i,j,n) = w1*f(i,j,kx) + wz*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlin_z_o_double

  ! The subroutine vertically interpolates or extrapolates if needed
  ! lowermost pressure, temperature and elevation from vertical profiles
  ! at fixed pressure levels.
  !
  ! THE ORDERING OF LEVELS MUST BE FROM BOTTOM TO TOP.
  !
  ! input:    tp      temps on global model pressure levels
  !           zp      heights of global model pressure levels
  !           zrcm    rcm topography
  !           sccm    global pressure levels (divided by surface)
  !           pss     global surface pressure
  ! output:   tlayer  mean layer temp above rcm surface
  !           pa      pressure at top of layer
  !           za      height at pressure pa
  !
  ! NOTES
  !        Ordering is BOTTOM -> TOP for the P coordinate
  !        Sigma is thus 1 -> 0
  !        Pressure coordinates for pss , psrccm and ptop MUST match.
  !        Lowermost level is 1
  !
  subroutine intgtb(pa,za,tlayer,zrcm,tp,zp,pss,sccm,pst,ni,nj,nz)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nz
    real(rkx) :: pss , pst
    real(rkx) , dimension(ni,nj) , intent(in) :: zrcm
    real(rkx) , dimension(nz) , intent(in) :: sccm
    real(rkx) , dimension(ni,nj,nz) , intent(in) :: tp , zp
    real(rkx) , dimension(ni,nj) , intent(out) :: tlayer , pa , za
    integer(ik4) :: i , j , k , kb , kt
    real(rkx) :: wu , wl

    do j = 1 , nj
      do i = 1 , ni
        if ( zrcm(i,j) < zp(i,j,1) ) then
          tlayer(i,j) = tp(i,j,1)
          za(i,j) = zp(i,j,1)
          pa(i,j) = pss + pst
        else if ( zrcm(i,j) > zp(i,j,nz) ) then
          write(stderr,*) 'REGIONAL MODEL ELEVATION HIGHER THAN GCM TOP'
          tlayer(i,j) = missl
          za(i,j) = missl
          pa(i,j) = missl
        else
          kb = 0
          do k = 1 , nz - 1
            if ( zrcm(i,j) >   zp(i,j,k)    .and. &
                 zrcm(i,j) <=  zp(i,j,k+1) ) then
              kb = k
              exit
            end if
          end do
          kt = kb + 1
          wu = (zrcm(i,j)-zp(i,j,kb))/(zp(i,j,kt)-zp(i,j,kb))
          wl = d_one - wu
          tlayer(i,j) = tp(i,j,kt) * wl + tp(i,j,kb) * wu
          tlayer(i,j) = (tp(i,j,kt) + tlayer(i,j))/d_two
          za(i,j) = zp(i,j,kt)
          pa(i,j) = pss*sccm(kt) + pst
        end if
      end do
    end do
  end subroutine intgtb

  ! The subroutine vertically interpolates on regular grid of pressure
  ! levels p from a non regular pressure levels grid p3d.
  ! The single and double precision version are one after the other.
  !
  ! The interpolation is logarithmic.
  !
  ! Input:
  !   f(im,jm,km)   - the field to be interpolated from
  !   p3d(im,jm,km) - the pressure levels the field is defined at
  !   p(kp)         - the pressure levels to interpolate to
  ! Output:
  !   fp(im,jm,kp)  - the field interpolate on regular pressure grid p
  !
  subroutine intlog_double(fp,f,p3d,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , dimension(im,jm,km) , intent(in) :: f , p3d
    real(rk8) , dimension(kp) , intent(in) :: p
    real(rk8) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk8) :: sigp , w1 , wp , tp , bp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) , dimension(km) :: sig

    if ( p3d(1,1,1) > p3d(1,1,km) ) then
      do j = 1 , jm
        do i = 1 , im
          tp = p3d(i,j,km)
          bp = p3d(i,j,1)-tp
          do k = 1 , km
            sig(k) = (p3d(i,j,k)-tp)/bp
          end do
          do n = 1 , kp
            sigp = (p(n)-tp)/bp
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
            else if ( sigp > sig(1) ) then
              fp(i,j,n) = d_half*(f(i,j,1)+f(i,j,2)) * &
                        dexp(rglrog*dlog(sigp/sig(1)))
            else
              kx = 2
              do k = 2 , km
                kx = k
                if ( sig(k) < sigp ) exit
              end do
              knx = kx - 1
              wp = (dlog(sigp)-dlog(sig(kx)))/(dlog(sig(knx))-dlog(sig(kx)))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          tp = p3d(i,j,1)
          bp = p3d(i,j,km)-tp
          do k = 1 , km
            sig(k) = (p3d(i,j,k)-tp)/bp
          end do
          do n = 1 , kp
            sigp = (p(n)-tp)/bp
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
            else if ( sigp > sig(km) ) then
              fp(i,j,n) = d_half*(f(i,j,km)+f(i,j,km-1)) * &
                        dexp(rglrog*dlog(sigp/sig(km)))
            else
              kx = 2
              do k = 2 , km
                kx = k
                if ( sig(k) > sigp ) exit
              end do
              knx = kx - 1
              wp = dlog(sig(kx)/sigp)/dlog(sig(kx)/sig(knx))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    end if
  end subroutine intlog_double

  subroutine intlog_single(fp,f,p3d,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk4) , dimension(im,jm,km) , intent(in) :: f , p3d
    real(rk4) , dimension(kp) , intent(in) :: p
    real(rk4) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk4) :: sigp , w1 , wp , tp , bp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk4) , dimension(km) :: sig

    if ( p3d(1,1,1) > p3d(1,1,km) ) then
      do j = 1 , jm
        do i = 1 , im
          tp = p3d(i,j,km)
          bp = p3d(i,j,1)-tp
          do k = 1 , km
            sig(k) = (p3d(i,j,k)-tp)/bp
          end do
          do n = 1 , kp
            sigp = (p(n)-tp)/bp
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
            else if ( sigp > sig(1) ) then
              fp(i,j,n) = 0.5*(f(i,j,1)+f(i,j,2)) * &
                        exp(real(rglrog,rk4)*log(sigp/sig(1)))
            else
              kx = 2
              do k = 2 , km
                kx = k
                if ( sig(k) < sigp ) exit
              end do
              knx = kx - 1
              wp = log(sigp/sig(kx))/log(sig(knx)/sig(kx))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          tp = p3d(i,j,1)
          bp = p3d(i,j,km)-tp
          do k = 1 , km
            sig(k) = (p3d(i,j,k)-tp)/bp
          end do
          do n = 1 , kp
            sigp = (p(n)-tp)/bp
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
            else if ( sigp > sig(km) ) then
              fp(i,j,n) = 0.5*(f(i,j,km)+f(i,j,km-1)) * &
                        exp(real(rglrog,rk4)*log(sigp/sig(km)))
            else
              kx = 2
              do k = 2 , km
                kx = k
                if ( sig(k) > sigp ) exit
              end do
              knx = kx - 1
              wp = (log(sigp)-log(sig(kx)))/(log(sig(knx))-log(sig(kx)))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    end if
  end subroutine intlog_single

  ! The subroutine vertically interpolates on regular grid of pressure
  ! levels p from a non regular hydrostatic pressure levels.
  ! The single and double precision version are one after the other.
  !
  ! The interpolation is logarithmic.
  !
  ! Input:
  !   f(im,jm,km)   - the field to be interpolated from
  !   ps(im,jm)     - the surface pressure values
  !   sig(km)       - the sigma coordinate values
  !   ptop          - the model top pressue
  !   p(kp)         - the pressure levels to interpolate to
  ! Output:
  !   fp(im,jm,kp)  - the field interpolate on regular pressure grid p
  !
  subroutine intlog_o_double(fp,f,ps,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , intent(in) :: ptop
    real(rk8) , dimension(im,jm,km) , intent(in) :: f
    real(rk8) , dimension(kp) , intent(in) :: p
    real(rk8) , dimension(im,jm) , intent(in) :: ps
    real(rk8) , dimension(km) , intent(in) :: sig
    real(rk8) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk8) :: sigp , w1 , wp
    integer(ik4) :: i , j , k , kx , knx , n
    if ( sig(1) > sig(2) ) then
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kp
            sigp = (p(n)-ptop)/(ps(i,j)-ptop)
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
            else if ( sigp > sig(1) ) then
              fp(i,j,n) = d_half*(f(i,j,1)+f(i,j,2)) * &
                             dexp(rglrog*dlog(sigp/sig(1)))
            else
              do k = 2 , km
                kx = k
                if ( sig(k) < sigp ) exit
              end do
              knx = kx - 1
              wp = dlog(sigp/sig(kx))/dlog(sig(knx)/sig(kx))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kp
            sigp = (p(n)-ptop)/(ps(i,j)-ptop)
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
            else if ( sigp > sig(km) ) then
              fp(i,j,n) = d_half*(f(i,j,km)+f(i,j,km-1)) * &
                        dexp(rglrog*dlog(sigp/sig(km)))
            else
              do k = 2 , km
                kx = k
                if ( sig(k) > sigp ) exit
              end do
              knx = kx - 1
              wp = dlog(sig(kx)/sigp)/dlog(sig(kx)/sig(knx))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    end if
  end subroutine intlog_o_double

  subroutine intlog_o_single(fp,f,ps,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rkx) , intent(in) :: ptop
    real(rk4) , dimension(im,jm,km) , intent(in) :: f
    real(rk4) , dimension(kp) , intent(in) :: p
    real(rk4) , dimension(im,jm) , intent(in) :: ps
    real(rk4) , dimension(km) , intent(in) :: sig
    real(rk4) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk4) :: sigp , w1 , wp , pt
    integer(ik4) :: i , j , k , kx , knx , n

    pt = real(ptop)
    if ( sig(1) > sig(2) ) then
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kp
            sigp = (p(n)-pt)/(ps(i,j)-pt)
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
            else if ( sigp > sig(1) ) then
              fp(i,j,n) = 0.5*(f(i,j,1)+f(i,j,2)) * &
                        exp(real(rglrog,rk4)*log(sigp/sig(1)))
            else
              do k = 2 , km
                kx = k
                if ( sig(k) < sigp ) exit
              end do
              knx = kx - 1
              wp = log(sigp/sig(kx))/log(sig(knx)/sig(kx))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    else
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kp
            sigp = (p(n)-pt)/(ps(i,j)-pt)
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
            else if ( sigp > sig(km) ) then
              fp(i,j,n) = 0.5*(f(i,j,km)+f(i,j,km-1)) * &
                        exp(real(rglrog,rk4)*log(sigp/sig(km)))
            else
              do k = 2 , km
                kx = k
                if ( sig(k) > sigp ) exit
              end do
              knx = kx - 1
              wp = log(sig(kx)/sigp)/log(sig(kx)/sig(knx))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end if
          end do
        end do
      end do
    end if
  end subroutine intlog_o_single

  ! This subroutine extrapolates surface pressure from closest level above.
  ! Output is pstar = surface pressure - ptop
  !
  subroutine intpsn(psrcm,zrcm,pa,za,tlayer,pt,ni,nj)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    real(rkx) , intent(in) :: pt
    real(rkx) , dimension(ni,nj) , intent(in) :: pa , tlayer , za , zrcm
    real(rkx) , dimension(ni,nj) , intent(out) :: psrcm
    integer(ik4) :: i , j

    do i = 1 , ni
      do j = 1 , nj
        psrcm(i,j) = pa(i,j)*exp(-govr*(zrcm(i,j)-za(i,j))/tlayer(i,j)) - pt
      end do
    end do
  end subroutine intpsn

  !
  ! INTV0 is for vertical interpolation with the same vertical ordering
  ! of regcm. The interpolation is linear in p.  Where extrapolation
  ! is necessary, fields are considered to have 0 vertical derivative.
  !
  subroutine intv0(frcm,fccm,psrcm,srcm,pss,sccm,pt,pst,ni,nj,krcm,kccm)
    implicit none
    integer(ik4) , intent(in) :: kccm , krcm , ni , nj
    real(rkx) , intent(in) :: pt , pss , pst
    real(rkx) , dimension(ni,nj,kccm) , intent(in) :: fccm
    real(rkx) , dimension(ni,nj) , intent(in) :: psrcm
    real(rkx) , dimension(kccm) , intent(in) :: sccm
    real(rkx) , dimension(krcm) , intent(in) :: srcm
    real(rkx) , dimension(ni,nj,krcm) , intent(out) :: frcm
    real(rkx) :: rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1 , n
    do i = 1 , ni
      do j = 1 , nj
        do n = 1 , krcm
          sc = (srcm(n)*psrcm(i,j) + pt - pst)/pss
          k1 = 0
          do k = 1 , kccm
            if ( sc > sccm(k) ) k1 = k
          end do
          if ( k1 == 0 ) then
            frcm(i,j,n) = fccm(i,j,1)
          else if ( k1 /= kccm ) then
            kp1 = k1+1
            rc = (sc-sccm(k1))/(sccm(kp1)-sccm(k1))
            rc1 = d_one-rc
            frcm(i,j,n) = rc1*fccm(i,j,k1)+rc*fccm(i,j,kp1)
          else
            frcm(i,j,n) = fccm(i,j,kccm)
          end if
        end do
      end do
    end do
  end subroutine intv0

  !
  ! INTVP is for vertical interpolation with the same vertical ordering
  ! of regcm. The interpolation is linear in p.  Where extrapolation
  ! is necessary, fields are considered to have 0 vertical derivative.
  ! For a particular sigma value only.
  !
  subroutine intvp(frcm,fccm,psrcm,srcm,pss,sccm,pt,pst,ni,nj,kccm)
    implicit none
    integer(ik4) , intent(in) :: kccm , ni , nj
    real(rkx) , intent(in) :: pt , pss , pst
    real(rkx) , dimension(ni,nj,kccm) , intent(in) :: fccm
    real(rkx) , dimension(ni,nj) , intent(in) :: psrcm
    real(rkx) , dimension(kccm) , intent(in) :: sccm
    real(rkx) , intent(in) :: srcm
    real(rkx) , dimension(ni,nj) , intent(out) :: frcm
    real(rkx) :: rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1

    do i = 1 , ni
      do j = 1 , nj
        sc = (srcm * psrcm(i,j) + pt - pst)/pss
        k1 = 0
        do k = 1 , kccm
          if ( sc > sccm(k) ) k1 = k
        end do
        if ( k1 == 0 ) then
          frcm(i,j) = fccm(i,j,kccm)
        else if ( k1 /= kccm ) then
          kp1 = k1 + 1
          rc = (sccm(k1)-sc)/(sccm(k1)-sccm(kp1))
          rc1 = d_one - rc
          frcm(i,j) = rc1*fccm(i,j,kccm-k1+1)+rc*fccm(i,j,kccm-kp1+1)
        else
          frcm(i,j) = fccm(i,j,1)
        end if
      end do
    end do
  end subroutine intvp

  ! Vertical interpolation to Z levels.
  !
  subroutine intz1(frcm,fccm,zrcm,zccm,trcm,ni,nj,krcm,kccm,a,e1,e2)
    implicit none
    integer(ik4) , intent(in) :: kccm , krcm , ni , nj
    real(rkx) , intent(in) :: a , e1 , e2
    real(rkx) , dimension(ni,nj,kccm) , intent(in) :: fccm , zccm
    real(rkx) , dimension(ni,nj,krcm) , intent(in) :: zrcm
    real(rkx) , dimension(ni,nj) , intent(in) :: trcm
    real(rkx) , dimension(ni,nj,krcm) , intent(out) :: frcm
    real(rkx) , dimension(kccm) :: xc , fc
    real(rkx) , dimension(krcm) :: xr , fr
    integer(ik4) :: i , j
    do j = 1 , nj
      do i = 1 , ni
        xc(:) = zccm(i,j,:)
        fc(:) = fccm(i,j,:)
        xr(:) = zrcm(i,j,:) + trcm(i,j)
        call interp1d(xc,fc,xr,fr,a,e1,e2)
        frcm(i,j,:) = fr(:)
      end do
    end do
  end subroutine intz1

  ! Vertical interpolation to P levels.
  !
  subroutine intp1(frcm,fccm,prcm,pccm,ni,nj,krcm,kccm,a,e1,e2)
    implicit none
    integer(ik4) , intent(in) :: kccm , krcm , ni , nj
    real(rkx) , intent(in) :: a , e1 , e2
    real(rkx) , dimension(ni,nj,kccm) , intent(in) :: fccm , pccm
    real(rkx) , dimension(ni,nj,krcm) , intent(in) :: prcm
    real(rkx) , dimension(ni,nj,krcm) , intent(out) :: frcm
    real(rkx) , dimension(kccm) :: xc , fc
    real(rkx) , dimension(krcm) :: xr , fr
    integer(ik4) :: i , j , kt , kb
    if ( pccm(1,1,1) > pccm(1,1,kccm) ) then
      kt = kccm
      kb = 1
    else
      kt = 1
      kb = kccm
    end if
    do j = 1 , nj
      do i = 1 , ni
        xc(:) = (pccm(i,j,:)-pccm(i,j,kt))/(pccm(i,j,kb)-pccm(i,j,kt))
        fc(:) = fccm(i,j,:)
        xr(:) = (prcm(i,j,:)-pccm(i,j,kt))/(pccm(i,j,kb)-pccm(i,j,kt))
        call interp1d(xc,fc,xr,fr,a,e1,e2)
        frcm(i,j,:) = fr(:)
      end do
    end do
  end subroutine intp1

  !
  ! INTV1 is for vertical interpolation of U, V, and RH
  ! The interpolation is linear in P.  Where extrapolation
  ! is necessary, fields are considered to have 0 vertical derivative.
  !
  ! Use method 2 for Q to set 10e-8 as minimum value.
  !
  ! NOTES
  !        INPUT Ordering is BOTTOM -> TOP for the P coordinate
  !        SCCM is thus 1 -> 0
  !        Pressure coordinates for pss , psrccm and ptop MUST match.
  !        Lowermost input level is 1
  !
  subroutine intv1(frcm,fccm,psrcm,srcm,pss,sccm,pt,pst,ni,nj,krcm,kccm,imeth)
    implicit none
    integer(ik4) , intent(in) :: kccm , krcm , ni , nj , imeth
    real(rkx) , intent(in) :: pt , pss , pst
    real(rkx) , dimension(ni,nj,kccm) , intent(in) :: fccm
    real(rkx) , dimension(ni,nj) , intent(in) :: psrcm
    real(rkx) , dimension(kccm) , intent(in) :: sccm
    real(rkx) , dimension(krcm) , intent(in) :: srcm
    real(rkx) , dimension(ni,nj,krcm) , intent(out) :: frcm
    real(rkx) :: rc1 , rc2 , sc
    integer(ik4) :: i , j , k , k1 , km1 , n

    if ( imeth == 1 ) then
      do j = 1 , nj
        do i = 1 , ni
          do n = 1 , krcm
            sc = ((srcm(n)*psrcm(i,j) + pt) - pst)/pss
            if ( sc > sccm(1) ) then
              frcm(i,j,n) = fccm(i,j,1)
            else if ( sc < sccm(kccm) ) then
              frcm(i,j,n) = fccm(i,j,kccm)
            else
              k1 = 0
              do k = 1 , kccm
                if ( sc > sccm(k) ) then
                  k1 = k
                  exit
                end if
              end do
              km1 = k1 - 1
              rc1 = (sc-sccm(km1))/(sccm(k1)-sccm(km1))
              rc2 = d_one - rc1
              frcm(i,j,n) = rc1*fccm(i,j,k1)+rc2*fccm(i,j,km1)
            end if
          end do
        end do
      end do
    else
      do j = 1 , nj
        do i = 1 , ni
          do n = 1 , krcm
            sc = ((srcm(n)*psrcm(i,j) + pt) - pst)/pss
            if ( sc > sccm(1) ) then
              frcm(i,j,n) = max(fccm(i,j,1),1.0e-8_rkx)
            else if ( sc < sccm(kccm) ) then
              frcm(i,j,n) = max(fccm(i,j,kccm),1.0e-8_rkx)
            else
              k1 = 0
              do k = 1 , kccm
                if ( sc > sccm(k) ) then
                  k1 = k
                  exit
                end if
              end do
              km1 = k1 - 1
              rc1 = (sc-sccm(km1))/(sccm(k1)-sccm(km1))
              rc2 = d_one - rc1
              frcm(i,j,n) = rc1 * max(fccm(i,j,k1),1.0e-8_rkx) + &
                            rc2 * max(fccm(i,j,km1),1.0e-8_rkx)
            end if
          end do
        end do
      end do
    end if
  end subroutine intv1

  subroutine intv2(frcm,fccm,psrcm,srcm,pss,sccm,pt,pst,ni,nj,krcm,kccm)
    implicit none
    integer(ik4) , intent(in) :: kccm , krcm , ni , nj
    real(rkx) , intent(in) :: pt , pss , pst
    real(rkx) , dimension(ni,nj,kccm) , intent(in) :: fccm
    real(rkx) , dimension(ni,nj) , intent(in) :: psrcm
    real(rkx) , dimension(kccm) , intent(in) :: sccm
    real(rkx) , dimension(krcm) , intent(in) :: srcm
    real(rkx) , dimension(ni,nj,krcm) , intent(out) :: frcm
    real(rkx) :: a1 , rc1 , rc2 , sc
    integer(ik4) :: i , j , k , k1 , km1 , n

    do j = 1 , nj
      do i = 1 , ni
        do n = 1 , krcm
          sc = ((srcm(n)*psrcm(i,j) + pt) - pst)/pss
          if ( sc > sccm(1) ) then
            a1 = rgas2*log(sc/sccm(1))
            frcm(i,j,n) = fccm(i,j,1)*(b1-a1)/(b1+a1)
          else if ( sc < sccm(kccm) ) then
            frcm(i,j,n) = fccm(i,j,kccm)
          else
            k1 = 0
            do k = 1 , kccm
              if ( sc > sccm(k) ) then
                k1 = k
                exit
              end if
            end do
            km1 = k1 - 1
            rc1 = log(sc/sccm(km1))/log(sccm(k1)/sccm(km1))
            rc2 = d_one - rc1
            frcm(i,j,n) = rc1*fccm(i,j,k1)+rc2*fccm(i,j,km1)
          end if
        end do
      end do
    end do
  end subroutine intv2

  subroutine intz3(fsrcm,fccm,zccm,zrcm,ni,nj,kccm,a,e1,e2)
    implicit none
    integer(ik4) , intent(in) :: kccm , ni , nj
    real(rkx) , intent(in) :: a , e1 , e2
    real(rkx) , dimension(ni,nj,kccm) , intent(in) :: fccm , zccm
    real(rkx) , dimension(ni,nj) , intent(in) :: zrcm
    real(rkx) , dimension(ni,nj) , intent(out) :: fsrcm
    real(rkx) , dimension(kccm) :: zc , fc
    real(rkx) , dimension(1) :: rx , rc
    integer(ik4) :: i , j
    do j = 1 , nj
      do i = 1 , ni
        zc = zccm(i,j,:)
        fc = fccm(i,j,:)
        rx(1) = zrcm(i,j)
        call interp1d(zc,fc,rx,rc,a,e1,e2)
        fsrcm(i,j) = rc(1)
      end do
    end do
  end subroutine intz3

  subroutine intp3(fsrcm,fccm,pccm,psrcm,ni,nj,kccm,a,e1,e2)
    implicit none
    integer(ik4) , intent(in) :: kccm , ni , nj
    real(rkx) , intent(in) :: a , e1 , e2
    real(rkx) , dimension(ni,nj,kccm) , intent(in) :: fccm , pccm
    real(rkx) , dimension(ni,nj) , intent(in) :: psrcm
    real(rkx) , dimension(ni,nj) , intent(out) :: fsrcm
    real(rkx) , dimension(kccm) :: zc , fc
    real(rkx) , dimension(1) :: rx , rc
    integer(ik4) :: i , j , kt , kb
    if ( pccm(1,1,1) > pccm(1,1,kccm) ) then
      kt = kccm
      kb = 1
    else
      kt = 1
      kb = kccm
    end if
    do j = 1 , nj
      do i = 1 , ni
        zc(:) = (pccm(i,j,:)-pccm(i,j,kt))/(pccm(i,j,kb)-pccm(i,j,kt))
        fc(:) = fccm(i,j,:)
        rx(1) = (psrcm(i,j)-pccm(i,j,kt))/(pccm(i,j,kb)-pccm(i,j,kt))
        call interp1d(zc,fc,rx,rc,a,e1,e2)
        fsrcm(i,j) = rc(1)
      end do
    end do
  end subroutine intp3

  !
  ! INTV3 is for vertical interpolation.  The interpolation
  ! is linear in log P.  Where extrapolation upward is necessary,
  ! the T field is considered to have 0 vertical derivative.
  ! Where extrapolation downward is necessary, the T field is
  ! considered to have a lapse rate of rlapse (k/m), and the
  ! thickness is determined hydrostatically from the mean of the
  ! two extreme temperatues in the layer.
  !
  ! NOTES
  !        Ordering is BOTTOM -> TOP for the P coordinate
  !        Sigma is thus 1 -> 0
  !        Pressure coordinates for pss , psrccm and ptop MUST match.
  !        Lowermost level is 1
  !
  subroutine intv3(fsccm,fccm,psrccm,pss,sccm,ptop,pst,ni,nj,kccm)
    implicit none
    integer(ik4) , intent(in) :: kccm , ni , nj
    real(rkx) , intent(in) :: ptop , pss , pst
    real(rkx) , dimension(ni,nj,kccm) , intent(in) :: fccm
    real(rkx) , dimension(ni,nj) , intent(in) :: psrccm
    real(rkx) , dimension(ni,nj) , intent(out) :: fsccm
    real(rkx) , dimension(kccm) , intent(in) :: sccm
    real(rkx) :: a1 , rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , km1

    do j = 1 , nj
      do i = 1 , ni
        sc = ((psrccm(i,j)+ptop) - pst)/pss
        if ( sc < sccm(kccm) ) then
          fsccm(i,j) = fccm(i,j,kccm)
        else if ( sc > sccm(1) ) then
          ! If the surface is below the GCM's lowest level, extrapolate
          a1 = rgas2*log(sc/sccm(1))
          fsccm(i,j) = fccm(i,j,1)*(b1-a1)/(b1+a1)
        else
          k1 = 0
          do k = 1 , kccm
            if ( sc > sccm(k) ) then
              k1 = k
              exit
            end if
          end do
          ! Interpolate the temperature between the two adjacent GCM levels
          km1 = k1 - 1
          rc = log(sc/sccm(km1))/log(sccm(k1)/sccm(km1))
          rc1 = d_one - rc
          fsccm(i,j) = rc*fccm(i,j,k1)+rc1*fccm(i,j,km1)
        end if
      end do
    end do
  end subroutine intv3

  subroutine intzps1(psrcm,zrcm,tp,zp,pss,sccm,pst,lat,jday,ni,nj,nz)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nz
    real(rkx) , intent(in) :: pss , jday , pst
    real(rkx) , dimension(ni,nj) , intent(in) :: zrcm , lat
    real(rkx) , dimension(ni,nj) , intent(out) :: psrcm
    real(rkx) , dimension(nz) , intent(in) :: sccm
    real(rkx) , dimension(ni,nj,nz) , intent(in) :: tp , zp
    integer(ik4) :: i , j , k , kb , kt
    real(rkx) :: wu , wl , tlayer , tva , tvb , pa , pb , za , zb , dz , lrt

    if ( zp(1,1,1) < zp(1,1,nz) ) then
      do j = 1 , nj
        do i = 1 , ni
          if ( zp(i,j,1) > zrcm(i,j) ) then
            pb = pss + pst
            za = zp(i,j,2)
            zb = zp(i,j,1)
            tva = tp(i,j,2)
            tvb = tp(i,j,1)
            dz = zb - zrcm(i,j)
            lrt = (tva-tvb)/(za-zb)
            !lrt = 0.65_rkx*lrt + 0.35_rkx*stdlrate(jday,lat(i,j))
            lrt = 0.65_rkx*lrt - 0.35_rkx*lrate
            tlayer = tvb - dz*lrt
            tlayer = (tvb + tlayer) * 0.5_rkx
            psrcm(i,j) = pb * exp(govr*dz/tlayer)
          else if ( zp(1,1,nz) < zrcm(i,j) ) then
            write(stderr,*) 'REGIONAL MODEL ELEVATION HIGHER THAN GCM TOP'
            psrcm(i,j) = missl
          else
            kb = 1
            do k = 1 , nz - 1
              if ( zrcm(i,j) <= zp(i,j,k+1) .and. &
                   zrcm(i,j) > zp(i,j,k) ) then
                kb = k
                exit
              end if
            end do
            kt = kb + 1
            pa = pss*sccm(kt) + pst
            tva = tp(i,j,kt)
            tvb = tp(i,j,kb)
            za = zp(i,j,kt)
            zb = zp(i,j,kb)
            dz = za - zrcm(i,j)
            wu = (zrcm(i,j)-zb)/(za-zb)
            wl = 1.0_rkx - wu
            tlayer = tva * wu + tvb * wl
            tlayer = (tva + tlayer) * 0.5_rkx
            psrcm(i,j) = pa * exp(govr*dz/tlayer)
          end if
        end do
      end do
    else
      do j = 1 , nj
        do i = 1 , ni
          if ( zp(i,j,nz) > zrcm(i,j) ) then
            pb = pss + pst
            za = zp(i,j,nz-1)
            zb = zp(i,j,nz)
            tva = tp(i,j,nz-1)
            tvb = tp(i,j,nz)
            dz = za - zrcm(i,j)
            lrt = (tva-tvb)/(za-zb)
            !lrt = 0.65_rkx*lrt + 0.35_rkx*stdlrate(jday,lat(i,j))
            lrt = 0.65_rkx*lrt - 0.35_rkx*lrate
            tlayer = tvb - dz*lrt
            tlayer = (tvb + tlayer) * 0.5_rkx
            psrcm(i,j) = pb * exp(govr*dz/tlayer)
          else if ( zp(1,1,1) < zrcm(i,j) ) then
            write(stderr,*) 'REGIONAL MODEL ELEVATION HIGHER THAN GCM TOP'
            psrcm(i,j) = missl
          else
            kb = nz
            do k = nz , 2 , -1
              if ( zrcm(i,j) <= zp(i,j,k-1) .and. &
                   zrcm(i,j) > zp(i,j,k) ) then
                kb = k
                exit
              end if
            end do
            kt = kb - 1
            pa = pss*sccm(kt) + pst
            za = zp(i,j,kt)
            zb = zp(i,j,kb)
            tva = tp(i,j,kt)
            tvb = tp(i,j,kb)
            dz = za - zrcm(i,j)
            wu = (zrcm(i,j)-zb)/(za-zb)
            wl = 1.0_rkx - wu
            tlayer = tva * wu + tvb * wl
            tlayer = (tva + tlayer) * 0.5_rkx
            psrcm(i,j) = pa * exp(govr*dz/tlayer)
          end if
        end do
      end do
    end if
  end subroutine intzps1

  subroutine intzps2(psrcm,zrcm,tp,zp,pp,lat,jday,ni,nj,nz)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nz
    real(rkx) , intent(in) :: jday
    real(rkx) , dimension(ni,nj) , intent(in) :: zrcm , lat
    real(rkx) , dimension(ni,nj) , intent(out) :: psrcm
    real(rkx) , dimension(ni,nj,nz) , intent(in) :: tp , zp , pp
    integer(ik4) :: i , j , k , kb , kt
    real(rkx) :: wu , wl , tlayer , pa , za , dz , lrt

    if ( zp(1,1,1) < zp(1,1,nz) ) then
      do j = 1 , nj
        do i = 1 , ni
          kb = 0
          do k = 1 , nz - 1
            if ( zrcm(i,j) <= zp(i,j,k+1) .and. &
                 zrcm(i,j) > zp(i,j,k) ) then
              kb = k
              exit
            end if
          end do
          if ( kb /= 0 ) then
            kt = kb + 1
            pa = pp(i,j,kt)
            za = zp(i,j,kt)
            dz = zrcm(i,j)-za
            wu = (zrcm(i,j)-zp(i,j,kb))/(zp(i,j,kt)-zp(i,j,kb))
            wl = d_one - wu
            tlayer = tp(i,j,kt) * wu + tp(i,j,kb) * wl
            tlayer = (tp(i,j,kt) + tlayer)/d_two
          else
            pa = pp(i,j,1)
            za = zp(i,j,1)
            dz = zrcm(i,j)-za
            lrt = (tp(i,j,2)-tp(i,j,1))/(zp(i,j,2)-zp(i,j,1))
            !lrt = 0.65_rkx*lrt + 0.35_rkx*stdlrate(jday,lat(i,j))
            lrt = 0.65_rkx*lrt - 0.35_rkx*lrate
            tlayer = tp(i,j,1) - 0.5_rkx*dz*lrt
          end if
          psrcm(i,j) = pa * exp(-govr*dz/tlayer)
        end do
      end do
    else
      do j = 1 , nj
        do i = 1 , ni
          kb = 0
          do k = 2 , nz - 1
            if ( zrcm(i,j) <= zp(i,j,k-1) .and. &
                 zrcm(i,j) > zp(i,j,k) ) then
              kb = k
              exit
            end if
          end do
          if ( kb /= 0 ) then
            kt = kb - 1
            pa = pp(i,j,kt)
            za = zp(i,j,kt)
            dz = zrcm(i,j)-za
            wu = (zrcm(i,j)-zp(i,j,kb))/(zp(i,j,kt)-zp(i,j,kb))
            wl = d_one - wu
            tlayer = tp(i,j,kt) * wu + tp(i,j,kb) * wl
            tlayer = (tp(i,j,kt) + tlayer)/d_two
          else
            pa = pp(i,j,nz)
            za = zp(i,j,nz)
            dz = zrcm(i,j)-za
            lrt = (tp(i,j,nz-1)-tp(i,j,nz))/(zp(i,j,nz-1)-zp(i,j,nz))
            !lrt = 0.65_rkx*lrt + 0.35_rkx*stdlrate(jday,lat(i,j))
            lrt = 0.65_rkx*lrt - 0.35_rkx*lrate
            tlayer = tp(i,j,nz) - 0.5_rkx*dz*lrt
          end if
          psrcm(i,j) = pa * exp(-govr*dz/tlayer)
        end do
      end do
    end if
  end subroutine intzps2

end module mod_vertint

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

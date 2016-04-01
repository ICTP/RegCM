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

  implicit none

  private

  real(rk8) , parameter :: rgas2 = rgas/d_two
  ! lrate is defined as a positive constant.
  real(rk8) , parameter :: rglrog = rgas*lrate*regrav
  real(rk8) , parameter :: b1 = -egrav/lrate

  interface intlin
    module procedure intlin_double
    module procedure intlin_single
    module procedure intlin_o_double
    module procedure intlin_o_single
    module procedure intlin_z_o_single
  end interface intlin

  interface intlog
    module procedure intlog_double
    module procedure intlog_single
    module procedure intlog_o_double
    module procedure intlog_o_single
  end interface intlog

  public :: intlin , intgtb , intlog
  public :: intpsn , intv0 , intv1 , intv2 , intv3
  public :: intlinreg , intlinprof

  contains

  subroutine intlinreg(fp,f,ps,p3d,im1,im2,jm1,jm2,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im1 , im2 , jm1 , jm2 , km , kp
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: f
    real(rk8) , pointer , dimension(:,:) , intent(in) :: ps
    real(rk8) , pointer , dimension(:) , intent(in) :: p
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: p3d
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: fp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) , dimension(im1:im2,jm1:jm2,kp) :: ff
    real(rk8) , dimension(kp) :: pp
    real(rk8) , dimension(kp) :: sig
    real(rk8) :: sigp , w1 , wp
    logical :: same_order
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    !
    same_order = .true.
    if ( (p(1) > p(kp) .and. p3d(im1,jm1,1) < p3d(im1,jm1,km)) .or. &
         (p(1) < p(kp) .and. p3d(im1,jm1,1) > p3d(im1,jm1,km)) ) then
      same_order = .false.
    end if
    if ( same_order ) then
      ff(:,:,:) = f(:,:,:)
      pp(:) = p(:)
    else
      do k = 1 , kp
        ff(:,:,kp-k+1) = f(:,:,k)
        pp(kp-k+1) = p(k)
      end do
    end if
    !
    ! HERE BOTTOM TO TOP
    !
    if ( pp(1) > pp(kp) ) then
      !
      ! Loop over points
      !
      do j = jm1 , jm2
        do i = im1 , im2
          if ( ps(i,j) < 1.0D0 ) cycle
          !
          ! Sigma values in this point
          !
          do k = 1 , kp
            sig(k) = pp(k)/ps(i,j)
          end do
          !
          ! For each of the requested levels
          !
          do n = 1 , km
            !
            ! The searched sigma value
            !
            sigp = p3d(i,j,n)/ps(i,j)
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(kp) ) then
              fp(i,j,n) = ff(i,j,kp)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = ff(i,j,1)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = 0
            do k = 1 , kp-1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx + 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*ff(i,j,kx) + wp*ff(i,j,knx)
          end do
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Loop over points
      !
      do j = jm1 , jm2
        do i = im1 , im2
          if ( ps(i,j) < 1.0D0 ) cycle
          !
          ! Sigma values in this point
          !
          do k = 1 , kp
            sig(k) = pp(k)/ps(i,j)
          end do
          !
          ! For each of the requested levels
          !
          do n = 1 , km
            !
            ! The searched sigma value
            !
            sigp = p3d(i,j,n)/ps(i,j)
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = ff(i,j,1)
              cycle
            else if ( sigp >= sig(kp) ) then
              fp(i,j,n) = ff(i,j,kp)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = kp + 1
            do k = kp , 2 , -1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx - 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*ff(i,j,kx) + wp*ff(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlinreg

  subroutine intlinprof(fp,f,ps,p3d,im1,im2,jm1,jm2,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im1 , im2 , jm1 , jm2 , km , kp
    real(rk8) , pointer , dimension(:) , intent(in) :: f
    real(rk8) , pointer , dimension(:,:) , intent(in) :: ps
    real(rk8) , pointer , dimension(:) , intent(in) :: p
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: p3d
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: fp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) , dimension(kp) :: ff
    real(rk8) , dimension(kp) :: pp
    real(rk8) , dimension(kp) :: sig
    real(rk8) :: sigp , w1 , wp
    logical :: same_order
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    !
    same_order = .true.
    if ( (p(1) > p(kp) .and. p3d(im1,jm1,1) < p3d(im1,jm1,km)) .or. &
         (p(1) < p(kp) .and. p3d(im1,jm1,1) > p3d(im1,jm1,km)) ) then
      same_order = .false.
    end if
    if ( same_order ) then
      ff(:) = f(:)
      pp(:) = p(:)
    else
      do k = 1 , kp
        ff(kp-k+1) = f(k)
        pp(kp-k+1) = p(k)
      end do
    end if
    !
    ! HERE BOTTOM TO TOP
    !
    if ( pp(1) > pp(kp) ) then
      !
      ! Loop over points
      !
      do j = jm1 , jm2
        do i = im1 , im2
          if ( ps(i,j) < 1.0D0 ) cycle
          !
          ! Sigma values in this point
          !
          do k = 1 , kp
            sig(k) = pp(k)/ps(i,j)
          end do
          !
          ! For each of the requested levels
          !
          do n = 1 , km
            !
            ! The searched sigma value
            !
            sigp = p3d(i,j,n)/ps(i,j)
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(kp) ) then
              fp(i,j,n) = ff(kp)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = ff(1)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = 0
            do k = 1 , kp-1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx + 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*ff(kx) + wp*ff(knx)
          end do
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Loop over points
      !
      do j = jm1 , jm2
        do i = im1 , im2
          if ( ps(i,j) < 1.0D0 ) cycle
          !
          ! Sigma values in this point
          !
          do k = 1 , kp
            sig(k) = pp(k)/ps(i,j)
          end do
          !
          ! For each of the requested levels
          !
          do n = 1 , km
            !
            ! The searched sigma value
            !
            sigp = p3d(i,j,n)/ps(i,j)
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = ff(1)
              cycle
            else if ( sigp >= sig(kp) ) then
              fp(i,j,n) = ff(kp)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = kp + 1
            do k = kp , 2 , -1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx - 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*ff(kx) + wp*ff(knx)
          end do
        end do
      end do
    end if
  end subroutine intlinprof

  subroutine intlin_double(fp,f,ps,p3d,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , dimension(im,jm,km) , intent(in) :: f , p3d
    real(rk8) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk8) , dimension(kp) , intent(in) :: p
    real(rk8) , dimension(im,jm) , intent(in) :: ps
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) , dimension(km) :: sig
    real(rk8) , dimension(kp) :: pp1
    real(rk8) :: sigp , w1 , wp
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    !
    ! HERE BOTTOM TO TOP
    !
    if ( p3d(1,1,1) > p3d(1,1,km) ) then
      !
      ! Assure same order requested in output
      !
      if ( p(1) > p(kp) ) then
        pp1 = p
      else
        do k = 1 , kp
          pp1(k) = p(kp-k+1)
        end do
      end if
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Discard missing values
          !
          if ( ps(i,j) > -9995.0D0 ) then
            !
            ! Sigma values in this point
            !
            do k = 1 , km
              sig(k) = p3d(i,j,k)/ps(i,j)
            end do
            !
            ! For each of the requested levels
            !
            do n = 1 , kp
              !
              ! The searched sigma value
              !
              sigp = pp1(n)/ps(i,j)
              !
              ! Over the top or below bottom level
              !
              if ( sigp <= sig(km) ) then
                fp(i,j,n) = f(i,j,km)
                cycle
              else if ( sigp >= sig(1) ) then
                fp(i,j,n) = f(i,j,1)
                cycle
              end if
              !
              ! Search k level below the requested one
              !
              kx = 0
              do k = 1 , km-1
                if ( sigp > sig(k) ) exit
                kx = k
              end do
              !
              ! This is the above level
              !
              knx = kx + 1
              wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end do
          else
            call die('intlin','Missing value in surface pressure',1)
          end if
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Assure same order requested in output
      !
      if ( p(1) > p(kp) ) then
        do k = 1 , kp
          pp1(k) = p(kp-k+1)
        end do
      else
        pp1 = p
      end if
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Discard missing values
          !
          if ( ps(i,j) > -9995.0D0 ) then
            !
            ! Sigma values in this point
            !
            do k = 1 , km
              sig(k) = p3d(i,j,k)/ps(i,j)
            end do
            !
            ! For each of the requested levels
            !
            do n = 1 , kp
              !
              ! The searched sigma value
              !
              sigp = pp1(n)/ps(i,j)
              !
              ! Over the top or below bottom level
              !
              if ( sigp <= sig(1) ) then
                fp(i,j,n) = f(i,j,1)
                cycle
              else if ( sigp >= sig(km) ) then
                fp(i,j,n) = f(i,j,km)
                cycle
              end if
              !
              ! Search k level below the requested one
              !
              kx = km + 1
              do k = km , 2 , -1
                if ( sigp > sig(k) ) exit
                kx = k
              end do
              !
              ! This is the above level
              !
              knx = kx - 1
              wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end do
          else
            call die('intlin','Missing value in surface pressure',1)
          end if
        end do
      end do
    end if
  end subroutine intlin_double

  subroutine intlin_single(fp,f,ps,p3d,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk4) , dimension(im,jm,km) , intent(in) :: f , p3d
    real(rk4) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk4) , dimension(kp) , intent(in) :: p
    real(rk4) , dimension(im,jm) , intent(in) :: ps
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk4) , dimension(km) :: sig
    real(rk4) , dimension(kp) :: pp1
    real(rk4) :: sigp , w1 , wp
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    !
    ! HERE BOTTOM TO TOP
    !
    if ( p3d(1,1,1) > p3d(1,1,km) ) then
      !
      ! Assure same order requested in output
      !
      if ( p(1) > p(kp) ) then
        pp1 = p
      else
        do k = 1 , kp
          pp1(k) = p(kp-k+1)
        end do
      end if
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Discard missing values
          !
          if ( ps(i,j) > -9995.0 ) then
            !
            ! Sigma values in this point
            !
            do k = 1 , km
              sig(k) = p3d(i,j,k)/ps(i,j)
            end do
            !
            ! For each of the requested levels
            !
            do n = 1 , kp
              !
              ! The searched sigma value
              !
              sigp = pp1(n)/ps(i,j)
              !
              ! Over the top or below bottom level
              !
              if ( sigp <= sig(km) ) then
                fp(i,j,n) = f(i,j,km)
                cycle
              else if ( sigp >= sig(1) ) then
                fp(i,j,n) = f(i,j,1)
                cycle
              end if
              !
              ! Search k level below the requested one
              !
              kx = 0
              do k = 1 , km-1
                if ( sigp > sig(k) ) exit
                kx = k
              end do
              !
              ! This is the above level
              !
              knx = kx + 1
              wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end do
          else
            call die('intlin','Missing value in surface pressure',1)
          end if
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Assure same order requested in output
      !
      if ( p(1) > p(kp) ) then
        do k = 1 , kp
          pp1(k) = p(kp-k+1)
        end do
      else
        pp1 = p
      end if
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Discard missing values
          !
          if ( ps(i,j) > -9995.0 ) then
            !
            ! Sigma values in this point
            !
            do k = 1 , km
              sig(k) = p3d(i,j,k)/ps(i,j)
            end do
            !
            ! For each of the requested levels
            !
            do n = 1 , kp
              !
              ! The searched sigma value
              !
              sigp = pp1(n)/ps(i,j)
              !
              ! Over the top or below bottom level
              !
              if ( sigp <= sig(1) ) then
                fp(i,j,n) = f(i,j,1)
                cycle
              else if ( sigp >= sig(km) ) then
                fp(i,j,n) = f(i,j,km)
                cycle
              end if
              !
              ! Search k level below the requested one
              !
              kx = km + 1
              do k = km , 2 , -1
                if ( sigp > sig(k) ) exit
                kx = k
              end do
              !
              ! This is the above level
              !
              knx = kx - 1
              wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end do
          else
            call die('intlin','Missing value in surface pressure',1)
          end if
        end do
      end do
    end if
  end subroutine intlin_single
  !
  !-----------------------------------------------------------------------
  !
  subroutine intlin_o_double(fp,f,pstar,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , intent(in) :: ptop
    real(rk8) , dimension(im,jm,km) , intent(in) :: f
    real(rk8) , dimension(kp) , intent(in) :: p
    real(rk8) , dimension(im,jm) , intent(in) :: pstar
    real(rk8) , dimension(km) , intent(in) :: sig
    real(rk8) , dimension(im,jm,kp) , intent(out) :: fp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) :: sigp , w1 , wp
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    !
    ! HERE BOTTOM TO TOP
    !
    if ( sig(1) > sig(2) ) then
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! For each of the requested levels
          !
          do n = 1 , kp
            !
            ! The searched sigma value
            !
            sigp = (p(n)-ptop)/(pstar(i,j)-ptop)
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = 0
            do k = 1 , km-1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx + 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! For each of the requested levels
          !
          do n = 1 , kp
            !
            ! The searched sigma value
            !
            sigp = (p(n)-ptop)/(pstar(i,j)-ptop)
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            else if ( sigp >= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = km + 1
            do k = km , 2 , -1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx - 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlin_o_double

  subroutine intlin_o_single(fp,f,pstar,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , intent(in) :: ptop
    real(rk4) , dimension(im,jm,km) , intent(in) :: f
    real(rk4) , dimension(kp) , intent(in) :: p
    real(rk4) , dimension(im,jm) , intent(in) :: pstar
    real(rk4) , dimension(km) , intent(in) :: sig
    real(rk4) , dimension(im,jm,kp) , intent(out) :: fp
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk4) :: sigp , w1 , wp , pt
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    !
    pt = real(ptop)
    !
    ! HERE BOTTOM TO TOP
    !
    if ( sig(1) > sig(2) ) then
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! For each of the requested levels
          !
          do n = 1 , kp
            !
            ! The searched sigma value
            !
            sigp = (p(n)-pt)/(pstar(i,j)-pt)
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = 0
            do k = 1 , km-1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx + 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = 1.0 - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! For each of the requested levels
          !
          do n = 1 , kp
            !
            ! The searched sigma value
            !
            sigp = (p(n)-pt)/(pstar(i,j)-pt)
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            else if ( sigp >= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = km + 1
            do k = km , 2 , -1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx - 1
            wp = (sigp-sig(kx))/(sig(knx)-sig(kx))
            w1 = 1.0 - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlin_o_single

  subroutine intlin_z_o_single(fz,f,hz,sig,im,jm,km,z,kz)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kz
    real(rk4) , dimension(im,jm,km) , intent(in) :: f , hz
    real(rk4) , dimension(kz) , intent(in) :: z
    real(rk4) , dimension(km) , intent(in) :: sig
    real(rk4) , dimension(im,jm,kz) , intent(out) :: fz
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk4) :: w1 , wz
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN Z.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    !
    ! HERE BOTTOM TO TOP
    !
    if ( sig(1) < sig(2) ) then
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! For each of the requested levels
          !
          do n = 1 , kz
            !
            ! Over the top or below bottom level
            !
            if ( z(n) <= hz(i,j,km) ) then
              fz(i,j,n) = f(i,j,km)
              cycle
            else if ( z(n) >= hz(i,j,1) ) then
              fz(i,j,n) = f(i,j,1)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = 0
            do k = 1 , km-1
              if ( z(n) > hz(i,j,k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx + 1
            wz = (z(n)-hz(i,j,kx))/(hz(i,j,knx)-hz(i,j,kx))
            w1 = 1.0 - wz
            fz(i,j,n) = w1*f(i,j,kx) + wz*f(i,j,knx)
          end do
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! For each of the requested levels
          !
          do n = 1 , kz
            !
            ! Over the top or below bottom level
            !
            if ( z(n) <= hz(i,j,1) ) then
              fz(i,j,n) = f(i,j,1)
              cycle
            else if ( z(n) >= hz(i,j,km) ) then
              fz(i,j,n) = f(i,j,km)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = km + 1
            do k = km , 2 , -1
              if ( z(n) > hz(i,j,k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx - 1
            wz = (z(n)-hz(i,j,kx))/(hz(i,j,knx)-hz(i,j,kx))
            w1 = 1.0 - wz
            fz(i,j,n) = w1*f(i,j,kx) + wz*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlin_z_o_single
  !
  !-----------------------------------------------------------------------
  !
  subroutine intgtb(pa,za,tlayer,zrcm,tp,zp,pss,sccm,ni,nj,nlev1)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nlev1
    real(rk8) :: pss
    real(rk8) , dimension(ni,nj) , intent(in) :: zrcm
    real(rk8) , dimension(nlev1) , intent(in) :: sccm
    real(rk8) , dimension(ni,nj,nlev1) , intent(in) :: tp , zp
    real(rk8) , dimension(ni,nj) , intent(out) :: tlayer , pa , za
    integer(ik4) :: i , j , k , kb , kt
    real(rk8) :: wu , wl
    !
    ! INTGTB CALCULATES ALL VARIABLES NEEDED TO COMPUTE P* ON THE RCM
    ! TOPOGRAPHY.  THE MEAN TEMPERATURE IN THE LAYER BETWEEN
    ! THE TOPOGRAPHY AND THE PRESSURE LEVEL ABOVE IS CALULATED
    ! BY LINEARLY INTERPOLATING WITH HEIGHT THE TEMPS ON
    ! PRESSURE LEVELS.
    !
    ! INPUT:    TP      TEMPS ON GLOBAL MODEL PRESSURE LEVELS
    !           ZP      HEIGHTS OF GLOBAL MODEL PRESSURE LEVELS
    !           ZRCM    RCM TOPOGRAPHY
    !           SCCM    GLOBAL PRESSURE LEVELS (DIVIDED BY SURFACE)
    ! OUTPUT:   TLAYER  MEAN LAYER TEMP ABOVE RCM SURFACE
    !           PA      PRESSURE AT TOP OF LAYER
    !           ZA      HEIGHT AT PRESSURE PA
    !
    do i = 1 , ni
      do j = 1 , nj
        kt = 0
        do k = 1 , nlev1 - 1
          if ( zrcm(i,j) <= zp(i,j,nlev1+1-k) .and. &
               zrcm(i,j) > zp(i,j,nlev1-k) ) kt = k
        end do
        kb = kt + 1
        if ( kt /= 0 ) then
          wu = ( zrcm(i,j) - zp(i,j,nlev1+1-kb) ) / &
               ( zp(i,j,nlev1+1-kt) - zp(i,j,nlev1+1-kb) )
          wl = d_one - wu
          tlayer(i,j) = tp(i,j,nlev1+1-kt) * wu + tp(i,j,nlev1+1-kb) * wl
          tlayer(i,j) = (tp(i,j,nlev1+1-kt) + tlayer(i,j))/d_two
          za(i,j) = zp(i,j,nlev1+1-kt)
          pa(i,j) = pss*sccm(kt)
        else
          tlayer(i,j) = tp(i,j,1)
          za(i,j) = zp(i,j,1)
          pa(i,j) = pss
        end if
      end do
    end do
  end subroutine intgtb
  !
  !-----------------------------------------------------------------------
  !
  subroutine intlog_double(fp,f,ps,p3d,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , dimension(im,jm,km) , intent(in) :: f , p3d
    real(rk8) , dimension(kp) , intent(in) :: p
    real(rk8) , dimension(im,jm) , intent(in) :: ps
    real(rk8) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk8) :: sigp , w1 , wp
    integer(ik4) :: i , j , k , kx , knx , n , kbc
    real(rk8) , dimension(km) :: sig
    real(rk8) , dimension(kp) :: pp1
    !
    ! INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
    ! LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
    ! THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    ! WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
    ! CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
    ! THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
    ! TWO EXTREME TEMPERATURES IN THE LAYER.
    !
    ! HERE BOTTOM TO TOP
    !
    if ( p3d(1,1,1) > p3d(1,1,km) ) then
      !
      ! Assure same order requested in output
      !
      if ( p(1) > p(kp) ) then
        pp1 = p
      else
        do k = 1 , kp
          pp1(k) = p(kp-k+1)
        end do
      end if
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Discard missing values
          !
          if ( ps(i,j) > -9995.0D0 ) then
            !
            ! Sigma values in this point , and find boundary layer
            !
            kbc = 1
            do k = 1 , km
              sig(k) = p3d(i,j,k)/ps(i,j)
              if ( sig(k) >= bltop ) kbc = k
            end do
            !
            ! For each of the requested levels
            !
            do n = 1 , kp
              !
              ! The searched sigma value
              !
              sigp = pp1(n)/ps(i,j)
              !
              ! Extrapolation
              !
              if ( sigp > d_one ) then
                fp(i,j,n) = f(i,j,kbc)*dexp(rglrog*dlog(sigp/sig(kbc)))
                cycle
              end if
              !
              ! Over the top or below bottom level
              !
              if ( sigp <= sig(km) ) then
                fp(i,j,n) = f(i,j,km)
                cycle
              else if ( sigp >= sig(1) ) then
                fp(i,j,n) = f(i,j,1)
                cycle
              end if
              !
              ! Search k level below the requested one
              !
              kx = 0
              do k = 1 , km
                if ( sigp > sig(k) ) exit
                kx = k
              end do
              !
              ! This is the above level
              !
              knx = kx + 1
              wp = dlog(sigp/sig(kx))/dlog(sig(knx)/sig(kx))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end do
          else
            call die('intlog','Missing value in surface pressure',1)
          end if
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Assure same order requested in output
      !
      if ( p(1) > p(kp) ) then
        do k = 1 , kp
          pp1(k) = p(kp-k+1)
        end do
      else
        pp1 = p
      end if
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Discard missing values
          !
          if ( ps(i,j) > -9995.0D0 ) then
            !
            ! Sigma values in this point , and find boundary layer
            !
            kbc = km
            do k = km , 1 , -1
              sig(k) = p3d(i,j,k)/ps(i,j)
              if ( sig(k) >= bltop ) kbc = k
            end do
            !
            ! For each of the requested levels
            !
            do n = 1 , kp
              !
              ! The searched sigma value
              !
              sigp = pp1(n)/ps(i,j)
              !
              ! Extrapolation
              !
              if ( sigp > d_one ) then
                fp(i,j,n) = f(i,j,kbc)*dexp(rglrog*dlog(sigp/sig(kbc)))
                cycle
              end if
              !
              ! Over the top or below bottom level
              !
              if ( sigp <= sig(1) ) then
                fp(i,j,n) = f(i,j,1)
                cycle
              else if ( sigp >= sig(km) ) then
                fp(i,j,n) = f(i,j,km)
                cycle
              end if
              !
              ! Search k level below the requested one
              !
              kx = km + 1
              do k = km , 2 , -1
                if ( sigp > sig(k) ) exit
                kx = k
              end do
              !
              ! This is the above level
              !
              knx = kx - 1
              wp = dlog(sigp/sig(kx))/dlog(sig(knx)/sig(kx))
              w1 = d_one - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end do
          else
            call die('intlog','Missing value in surface pressure',1)
          end if
        end do
      end do
    end if
  end subroutine intlog_double

  subroutine intlog_single(fp,f,ps,p3d,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk4) , dimension(im,jm,km) , intent(in) :: f , p3d
    real(rk4) , dimension(kp) , intent(in) :: p
    real(rk4) , dimension(im,jm) , intent(in) :: ps
    real(rk4) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk4) :: sigp , w1 , wp
    integer(ik4) :: i , j , k , kx , knx , n , kbc
    real(rk4) , dimension(km) :: sig
    real(rk4) , dimension(kp) :: pp1
    !
    ! INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
    ! LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
    ! THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    ! WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
    ! CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
    ! THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
    ! TWO EXTREME TEMPERATURES IN THE LAYER.
    !
    ! HERE BOTTOM TO TOP
    !
    if ( p3d(1,1,1) > p3d(1,1,km) ) then
      !
      ! Assure same order requested in output
      !
      if ( p(1) > p(kp) ) then
        pp1 = p
      else
        do k = 1 , kp
          pp1(k) = p(kp-k+1)
        end do
      end if
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Discard missing values
          !
          if ( ps(i,j) > -9995.0 ) then
            !
            ! Sigma values in this point , and find boundary layer
            !
            kbc = 1
            do k = 1 , km
              sig(k) = p3d(i,j,k)/ps(i,j)
              if ( sig(k) >= bltop ) kbc = k
            end do
            !
            ! For each of the requested levels
            !
            do n = 1 , kp
              !
              ! The searched sigma value
              !
              sigp = pp1(n)/ps(i,j)
              !
              ! Extrapolation
              !
              if ( sigp > 1.0 ) then
                fp(i,j,n) = f(i,j,kbc)*exp(real(rglrog)*log(sigp/sig(kbc)))
                cycle
              end if
              !
              ! Over the top or below bottom level
              !
              if ( sigp <= sig(km) ) then
                fp(i,j,n) = f(i,j,km)
                cycle
              else if ( sigp >= sig(1) ) then
                fp(i,j,n) = f(i,j,1)
                cycle
              end if
              !
              ! Search k level below the requested one
              !
              kx = 0
              do k = 1 , km
                if ( sigp > sig(k) ) exit
                kx = k
              end do
              !
              ! This is the above level
              !
              knx = kx + 1
              wp = log(sigp/sig(kx))/log(sig(knx)/sig(kx))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end do
          else
            call die('intlog','Missing value in surface pressure',1)
          end if
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Assure same order requested in output
      !
      if ( p(1) > p(kp) ) then
        do k = 1 , kp
          pp1(k) = p(kp-k+1)
        end do
      else
        pp1 = p
      end if
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Discard missing values
          !
          if ( ps(i,j) > -9995.0 ) then
            !
            ! Sigma values in this point , and find boundary layer
            !
            kbc = km
            do k = km , 1 , -1
              sig(k) = p3d(i,j,k)/ps(i,j)
              if ( sig(k) >= bltop ) kbc = k
            end do
            !
            ! For each of the requested levels
            !
            do n = 1 , kp
              !
              ! The searched sigma value
              !
              sigp = pp1(n)/ps(i,j)
              !
              ! Extrapolation
              !
              if ( sigp > 1.0 ) then
                fp(i,j,n) = f(i,j,kbc)*exp(real(rglrog)*log(sigp/sig(kbc)))
                cycle
              end if
              !
              ! Over the top or below bottom level
              !
              if ( sigp <= sig(1) ) then
                fp(i,j,n) = f(i,j,1)
                cycle
              else if ( sigp >= sig(km) ) then
                fp(i,j,n) = f(i,j,km)
                cycle
              end if
              !
              ! Search k level below the requested one
              !
              kx = km + 1
              do k = km , 2 , -1
                if ( sigp > sig(k) ) exit
                kx = k
              end do
              !
              ! This is the above level
              !
              knx = kx - 1
              wp = log(sigp/sig(kx))/log(sig(knx)/sig(kx))
              w1 = 1.0 - wp
              fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
            end do
          else
            call die('intlog','Missing value in surface pressure',1)
          end if
        end do
      end do
    end if
  end subroutine intlog_single
  !
  !-----------------------------------------------------------------------
  !
  subroutine intlog_o_double(fp,f,pstar,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , intent(in) :: ptop
    real(rk8) , dimension(im,jm,km) , intent(in) :: f
    real(rk8) , dimension(kp) , intent(in) :: p
    real(rk8) , dimension(im,jm) , intent(in) :: pstar
    real(rk8) , dimension(km) , intent(in) :: sig
    real(rk8) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk8) :: sigp , w1 , wp
    integer(ik4) :: i , j , k , kx , knx , kbc , n
    !
    ! INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
    ! LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
    ! THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    ! WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
    ! CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
    ! THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
    ! TWO EXTREME TEMPERATURES IN THE LAYER.
    !

    !
    ! HERE BOTTOM TO TOP
    !
    if ( sig(1) > sig(2) ) then
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! find boundary layer top
          !
          kbc = 1
          do k = 1 , km
            if ( sig(k) >= bltop ) kbc = k
          end do
          !
          ! For each of the requested levels
          !
          do n = 1 , kp
            !
            ! The searched sigma value
            !
            sigp = (p(n)-ptop)/(pstar(i,j)-ptop)
            !
            ! Extrapolation
            !
            if ( sigp > d_one ) then
              fp(i,j,n) = f(i,j,kbc)*dexp(rglrog*dlog(sigp/sig(kbc)))
              cycle
            end if
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = 0
            do k = 1 , km
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx + 1
            wp = dlog(sigp/sig(kx))/dlog(sig(knx)/sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Find boundary layer Top
          !
          kbc = km
          do k = km , 1 , -1
            if ( sig(k) >= bltop ) kbc = k
          end do
          !
          ! For each of the requested levels
          !
          do n = 1 , kp
            !
            ! The searched sigma value
            !
            sigp = (p(n)-ptop)/(pstar(i,j)-ptop)
            !
            ! Extrapolation
            !
            if ( sigp > d_one ) then
              fp(i,j,n) = f(i,j,kbc)*dexp(rglrog*dlog(sigp/sig(kbc)))
              cycle
            end if
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            else if ( sigp >= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = km + 1
            do k = km , 2 , -1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx - 1
            wp = dlog(sigp/sig(kx))/dlog(sig(knx)/sig(kx))
            w1 = d_one - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlog_o_double

  subroutine intlog_o_single(fp,f,pstar,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , intent(in) :: ptop
    real(rk4) , dimension(im,jm,km) , intent(in) :: f
    real(rk4) , dimension(kp) , intent(in) :: p
    real(rk4) , dimension(im,jm) , intent(in) :: pstar
    real(rk4) , dimension(km) , intent(in) :: sig
    real(rk4) , dimension(im,jm,kp) , intent(out) :: fp
    real(rk4) :: sigp , w1 , wp , pt
    integer(ik4) :: i , j , k , kx , knx , kbc , n
    !
    ! INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
    ! LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
    ! THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    ! WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
    ! CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
    ! THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
    ! TWO EXTREME TEMPERATURES IN THE LAYER.
    !
    pt = real(ptop)
    !
    ! HERE BOTTOM TO TOP
    !
    if ( sig(1) > sig(2) ) then
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! find boundary layer top
          !
          kbc = 1
          do k = 1 , km
            if ( sig(k) >= bltop ) kbc = k
          end do
          !
          ! For each of the requested levels
          !
          do n = 1 , kp
            !
            ! The searched sigma value
            !
            sigp = (p(n)-pt)/(pstar(i,j)-pt)
            !
            ! Extrapolation
            !
            if ( sigp > 1.0 ) then
              fp(i,j,n) = f(i,j,kbc)*exp(real(rglrog)*log(sigp/sig(kbc)))
              cycle
            end if
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            else if ( sigp >= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = 0
            do k = 1 , km
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx + 1
            wp = log(sigp/sig(kx))/log(sig(knx)/sig(kx))
            w1 = 1.0 - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    !
    ! HERE TOP TO BOTTOM
    !
    else
      !
      ! Loop over points
      !
      do j = 1 , jm
        do i = 1 , im
          !
          ! Find boundary layer Top
          !
          kbc = km
          do k = km , 1 , -1
            if ( sig(k) >= bltop ) kbc = k
          end do
          !
          ! For each of the requested levels
          !
          do n = 1 , kp
            !
            ! The searched sigma value
            !
            sigp = (p(n)-pt)/(pstar(i,j)-pt)
            !
            ! Extrapolation
            !
            if ( sigp > 1.0 ) then
              fp(i,j,n) = f(i,j,kbc)*exp(real(rglrog)*log(sigp/sig(kbc)))
              cycle
            end if
            !
            ! Over the top or below bottom level
            !
            if ( sigp <= sig(1) ) then
              fp(i,j,n) = f(i,j,1)
              cycle
            else if ( sigp >= sig(km) ) then
              fp(i,j,n) = f(i,j,km)
              cycle
            end if
            !
            ! Search k level below the requested one
            !
            kx = km + 1
            do k = km , 2 , -1
              if ( sigp > sig(k) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx - 1
            wp = log(sigp/sig(kx))/log(sig(knx)/sig(kx))
            w1 = 1.0 - wp
            fp(i,j,n) = w1*f(i,j,kx) + wp*f(i,j,knx)
          end do
        end do
      end do
    end if
  end subroutine intlog_o_single
  !
  !-----------------------------------------------------------------------
  !
  subroutine intpsn(psrcm,zrcm,pa,za,tlayer,pt,ni,nj)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    real(rk8) , intent(in) :: pt
    real(rk8) , dimension(ni,nj) , intent(in) :: pa , tlayer , za , zrcm
    real(rk8) , dimension(ni,nj) , intent(out) :: psrcm
    integer(ik4) :: i , j
    !
    ! EXTRAPOLATE SURFACE PRESSURE FROM CLOSEST PRESSURE LEVEL ABOVE.
    ! USE TLAYER CALCULATED IN INTGTB.
    ! PSRCM = SURFACE PRESSURE - PTOP
    !
    do i = 1 , ni
      do j = 1 , nj
        psrcm(i,j) = pa(i,j)*exp(-govr*(zrcm(i,j)-za(i,j))/tlayer(i,j)) - pt
      end do
    end do
  end subroutine intpsn
  !
  !-----------------------------------------------------------------------
  !
  subroutine intv0(frcm,fccm,psrcm,srcm,pss,sccm,pt,ni,nj,krcm,kccm)
    implicit none
    integer(ik4) , intent(in) :: kccm , krcm , ni , nj
    real(rk8) , intent(in) :: pt , pss
    real(rk8) , dimension(ni,nj,kccm) , intent(in) :: fccm
    real(rk8) , dimension(ni,nj) , intent(in) :: psrcm
    real(rk8) , dimension(kccm) , intent(in) :: sccm
    real(rk8) , dimension(krcm) , intent(in) :: srcm
    real(rk8) , dimension(ni,nj,krcm) , intent(out) :: frcm
    real(rk8) :: dp1 , pt1 , rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1 , n
    !
    ! INTV0 is for vertical interpolation of tracer where the tracer has
    ! the same vertical ordering of regcm.
    ! The interpolation is linear in p.  Where extrapolation
    ! is necessary, fields are considered to have 0 vertical derivative.
    !
    pt1 = pt/pss
    do i = 1 , ni
      do j = 1 , nj
        dp1 = psrcm(i,j)/pss
        do n = 1 , krcm
          sc = srcm(n)*dp1 + pt1
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

  subroutine intv1(frcm,fccm,psrcm,srcm,pss,sccm,pt,ni,nj,krcm,kccm)
    implicit none
    integer(ik4) , intent(in) :: kccm , krcm , ni , nj
    real(rk8) , intent(in) :: pt , pss
    real(rk8) , dimension(ni,nj,kccm) , intent(in) :: fccm
    real(rk8) , dimension(ni,nj) , intent(in) :: psrcm
    real(rk8) , dimension(kccm) , intent(in) :: sccm
    real(rk8) , dimension(krcm) , intent(in) :: srcm
    real(rk8) , dimension(ni,nj,krcm) , intent(out) :: frcm
    real(rk8) :: dp1 , pt1 , rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1 , n
    !
    ! INTV1 is for vertical interpolation of U, V, and RH
    ! The interpolation is linear in P.  Where extrapolation
    ! is necessary, fields are considered to have 0 vertical derivative.
    !
    pt1 = pt/pss
    do i = 1 , ni
      do j = 1 , nj
        dp1 = psrcm(i,j)/pss
        do n = 1 , krcm
          sc = srcm(n)*dp1 + pt1
          k1 = 0
          do k = 1 , kccm
            if ( sc > sccm(k) ) k1 = k
          end do
          if ( k1 == 0 ) then
            frcm(i,j,n) = fccm(i,j,kccm)
          else if ( k1 /= kccm ) then
            kp1 = k1 + 1
            rc = (sccm(k1)-sc)/(sccm(k1)-sccm(kp1))
            rc1 = d_one - rc
            frcm(i,j,n) = rc1*fccm(i,j,kccm-k1+1)+rc*fccm(i,j,kccm-kp1+1)
          else
            frcm(i,j,n) = fccm(i,j,1)
          end if
        end do
      end do
    end do
  end subroutine intv1
  !
  !-----------------------------------------------------------------------
  !
  subroutine intv2(frcm,fccm,psrcm,srcm,pss,sccm,pt,ni,nj,krcm,kccm)
    implicit none
    integer(ik4) , intent(in) :: kccm , krcm , ni , nj
    real(rk8) , intent(in) :: pt , pss
    real(rk8) , dimension(ni,nj,kccm) , intent(in) :: fccm
    real(rk8) , dimension(ni,nj) , intent(in) :: psrcm
    real(rk8) , dimension(kccm) , intent(in) :: sccm
    real(rk8) , dimension(krcm) , intent(in) :: srcm
    real(rk8) , dimension(ni,nj,krcm) , intent(out) :: frcm
    real(rk8) :: a1 , dp1 , pt1 , rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1 , n
    !
    ! INTV2 is for vertical interpolation.  The interpolation is
    ! linear in log P.  Where extrapolation upward is necessary,
    ! the T field is considered to have 0 vertical derivative.
    ! Where extrapolation downward is necessary, the T field is
    ! considered to have a lapse rate of rlapse (k/m), and the
    ! thickness is determined hydrostatically from the mean of the
    ! two extreme temperatues in the layer.
    !
    pt1 = pt/pss
    do i = 1 , ni
      do j = 1 , nj
        dp1 = psrcm(i,j)/pss
        do n = 1 , krcm
          sc = srcm(n)*dp1 + pt1
          k1 = 0
          do k = 1 , kccm
            if ( sc > sccm(k) ) k1 = k
          end do
          if ( k1 == 0 ) then
            frcm(i,j,n) = fccm(i,j,kccm)
          else if ( k1 /= kccm ) then
            kp1 = k1 + 1
            rc = log(sccm(k1)/sc)/log(sccm(k1)/sccm(kp1))
            rc1 = d_one - rc
            frcm(i,j,n) = rc1*fccm(i,j,kccm-k1+1)+rc*fccm(i,j,kccm-kp1+1)
          else
            a1 = rgas2*log(sc/sccm(kccm))
            frcm(i,j,n) = fccm(i,j,1)*(b1-a1)/(b1+a1)
          end if
        end do
      end do
    end do
  end subroutine intv2
  !
  !-----------------------------------------------------------------------
  !
  subroutine intv3(fsccm,fccm,psrccm,pss,sccm,ptop,ni,nj,kccm)
    implicit none
    integer(ik4) , intent(in) :: kccm , ni , nj
    real(rk8) , intent(in) :: ptop , pss
    real(rk8) , dimension(ni,nj,kccm) , intent(in) :: fccm
    real(rk8) , dimension(ni,nj) , intent(in) :: psrccm
    real(rk8) , dimension(ni,nj) , intent(out) :: fsccm
    real(rk8) , dimension(kccm) , intent(in) :: sccm
    real(rk8) :: a1 , rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1
    !
    ! INTV3 is for vertical interpolation.  The interpolation
    ! is linear in log P.  Where extrapolation upward is necessary,
    ! the T field is considered to have 0 vertical derivative.
    ! Where extrapolation downward is necessary, the T field is
    ! considered to have a lapse rate of rlapse (k/m), and the
    ! thickness is determined hydrostatically from the mean of the
    ! two extreme temperatues in the layer.
    !
    do i = 1 , ni
      do j = 1 , nj
        sc = (psrccm(i,j)+ptop)/pss
        k1 = 0
        do k = 1 , kccm
          if ( sc > sccm(k) ) k1 = k
        end do
        if ( k1 == 0 ) then
          fsccm(i,j) = fccm(i,j,1)
        else if ( k1 /= kccm ) then
          ! Interpolate the temperature between
          ! the two adjacent GCM levels
          kp1 = k1 + 1
          rc = log(sccm(k1)/sc)/log(sccm(k1)/sccm(kp1))
          rc1 = d_one - rc
          fsccm(i,j) = rc1*fccm(i,j,kccm+1-k1)+rc*fccm(i,j,kccm+1-kp1)
        else
          ! If the surface is below the GCM's lowest level,
          ! extrapolate
          a1 = rgas2*log(sc/sccm(kccm))
          fsccm(i,j) = fccm(i,j,kccm+1-kccm)*(b1-a1)/(b1+a1)
        end if
      end do
    end do
  end subroutine intv3

end module mod_vertint

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

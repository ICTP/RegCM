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
!
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_constants
  use mod_message
!
  private

  real(rk8) , parameter :: rgas2 = rgas/d_two
  ! lrate is defined as a positive constant.
  real(rk8) , parameter :: rglrog = rgas*lrate*regrav
  real(rk8) , parameter :: b1 = -egrav/lrate
  real(rk8) , parameter :: psccm = 100.0D0
!
  interface intlin_o
    module procedure intlin_o_double
    module procedure intlin_o_single
  end interface intlin_o

  interface intlog_o
    module procedure intlog_o_double
    module procedure intlog_o_single
  end interface intlog_o

  public :: intlin , intlin_o , intgtb , intlog , intlog_o
  public :: intpsn , intv0 , intv1 , intv2 , intv3
!
  contains

  subroutine intlin(fp,f,ps,p3d,im,jm,km,p,kp)
    implicit none
    integer(ik4) :: im , jm , km , kp
    real(rk8) , dimension(im,jm,km) :: f , p3d
    real(rk8) , dimension(im,jm,kp) :: fp
    real(rk8) , dimension(kp) :: p
    real(rk8) , dimension(im,jm) :: ps
    intent (in) f , im , jm , km , kp , p , ps , p3d
    intent (out) fp
!
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) , dimension(km) :: sig
    real(rk8) , dimension(kp) :: pp1
    real(rk8) :: sigp , w1 , wp

    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    !

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
  end subroutine intlin
!
!-----------------------------------------------------------------------
!
  subroutine intlin_o_double(fp,f,pstar,sig,ptop,im,jm,km,p,kp)
    implicit none
!
    integer(ik4) :: im , jm , km , kp
    real(rk8) :: ptop
    real(rk8) , dimension(im,jm,km) :: f
    real(rk8) , dimension(im,jm,kp) :: fp
    real(rk8) , dimension(kp) :: p
    real(rk8) , dimension(im,jm) :: pstar
    real(rk8) , dimension(km) :: sig
    intent (in) f , im , jm , km , kp , p , pstar , ptop , sig
    intent (out) fp
!
    integer(ik4) :: i , j , k , kx , knx , n
    real(rk8) :: sigp , w1 , wp
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
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
!
    integer(ik4) :: im , jm , km , kp
    real(rk8) :: ptop
    real(rk4) , dimension(im,jm,km) :: f
    real(rk4) , dimension(im,jm,kp) :: fp
    real(rk4) , dimension(kp) :: p
    real(rk4) , dimension(im,jm) :: pstar
    real(rk4) , dimension(km) :: sig
    intent (in) f , im , jm , km , kp , p , pstar , ptop , sig
    intent (out) fp
!
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
!
!-----------------------------------------------------------------------
!
  subroutine intgtb(pa,za,tlayer,zrcm,tp,zp,sccm,ni,nj,nlev1)
    implicit none
!
    integer(ik4) :: ni , nj , nlev1
    real(rk8) , dimension(ni,nj) :: pa , tlayer , za , zrcm
    real(rk8) , dimension(nlev1) :: sccm
    real(rk8) , dimension(ni,nj,nlev1) :: tp , zp
    intent (in) ni , nj , nlev1 , sccm , tp , zp , zrcm
    intent (out) pa , za
    intent (inout) tlayer
!
    integer(ik4) :: i , j , k , kb , kt
    !
    !     INTGTB CALCULATES ALL VARIABLES NEEDED TO COMPUTE P* ON THE RCM
    !     TOPOGRAPHY.  THE MEAN TEMPERATURE IN THE LAYER BETWEEN
    !     THE TOPOGRAPHY AND THE PRESSURE LEVEL ABOVE IS CALULATED
    !     BY LINEARLY INTERPOLATING WITH HEIGHT THE TEMPS ON
    !     PRESSURE LEVELS.
    !     INPUT:    TP        TEMPS ON ECMWF PRESSURE LEVELS
    !     ZP        HEIGHTS OF ECMWF PRESSURE LEVELS
    !     ZRCM      RCM TOPOGRAPHY
    !     SCCM      ECMWF PRESSURE LEVELS (DIVIDED BY 1000.)
    !     OUTPUT:   TLAYER    MEAN LAYER TEMP ABOVE RCM SURFACE
    !     PA        PRESSURE AT TOP OF LAYER
    !     ZA        HEIGHT AT PRESSURE PA
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
          tlayer(i,j) = (tp(i,j,nlev1+1-kt)*(zrcm(i,j)-          &
                         zp(i,j,nlev1+1-kb))+tp(i,j,nlev1+1-kb)  &
                        *(zp(i,j,nlev1+1-kt)-zrcm(i,j)))         &
                        /(zp(i,j,nlev1+1-kt)-zp(i,j,nlev1+1-kb))
          tlayer(i,j) = (tp(i,j,nlev1+1-kt)+tlayer(i,j))/d_two
          za(i,j) = zp(i,j,nlev1+1-kt)
          pa(i,j) = d_100*sccm(kt)
        else
          tlayer(i,j) = tp(i,j,1)
          za(i,j) = zp(i,j,1)
          pa(i,j) = d_100
        end if
      end do
    end do
  end subroutine intgtb
!
!-----------------------------------------------------------------------
!
  subroutine intlog(fp,f,ps,p3d,im,jm,km,p,kp)
    implicit none
!
    integer(ik4) :: im , jm , km , kp
    real(rk8) , dimension(im,jm,km) :: f , p3d
    real(rk8) , dimension(im,jm,kp) :: fp
    real(rk8) , dimension(kp) :: p
    real(rk8) , dimension(im,jm) :: ps
    intent (in) f , im , jm , km , kp , p , p3d , ps
    intent (out) fp
!
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
  end subroutine intlog
!
!-----------------------------------------------------------------------
!
  subroutine intlog_o_double(fp,f,pstar,sig,ptop,im,jm,km,p,kp)
    implicit none
!
    integer(ik4) :: im , jm , km , kp
    real(rk8) :: ptop
    real(rk8) , dimension(im,jm,km) :: f
    real(rk8) , dimension(im,jm,kp) :: fp
    real(rk8) , dimension(kp) :: p
    real(rk8) , dimension(im,jm) :: pstar
    real(rk8) , dimension(km) :: sig
    intent (in) f , im , jm , km , kp , p , pstar , ptop , sig
    intent (out) fp
!
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
!
    integer(ik4) :: im , jm , km , kp
    real(rk8) :: ptop
    real(rk4) , dimension(im,jm,km) :: f
    real(rk4) , dimension(im,jm,kp) :: fp
    real(rk4) , dimension(kp) :: p
    real(rk4) , dimension(im,jm) :: pstar
    real(rk4) , dimension(km) :: sig
    intent (in) f , im , jm , km , kp , p , pstar , ptop , sig
    intent (out) fp
!
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
    integer(ik4) :: ni , nj
    real(rk8) :: pt
    real(rk8) , dimension(ni,nj) :: pa , psrcm , tlayer , za , zrcm
    intent (in) ni , nj , pa , pt , tlayer , za , zrcm
    intent (out) psrcm
!
    real(rk8) :: tb
    integer(ik4) :: i , j
    !
    ! EXTRAPOLATE SURFACE PRESSURE FROM CLOSEST PRESSURE LEVEL ABOVE.
    ! USE TLAYER CALCULATED IN INTGTB.
    ! PSRCM = SURFACE PRESSURE - PTOP
    !
    do i = 1 , ni
      do j = 1 , nj
        tb = tlayer(i,j)
        psrcm(i,j) = pa(i,j)*dexp(-govr*(zrcm(i,j)-za(i,j))/tb)-pt
      end do
    end do
  end subroutine intpsn
!
!-----------------------------------------------------------------------
!
  subroutine intv0(frcm,fccm,psrcm,srcm,sccm,pt,ni,nj,krcm,kccm)
    implicit none
!
    integer(ik4) :: kccm , krcm , ni , nj
    real(rk8) :: pt
    real(rk8) , dimension(ni,nj,kccm) :: fccm
    real(rk8) , dimension(ni,nj,krcm) :: frcm
    real(rk8) , dimension(ni,nj) :: psrcm
    real(rk8) , dimension(kccm) :: sccm
    real(rk8) , dimension(krcm) :: srcm
    intent (in) fccm , kccm , krcm , ni , nj , psrcm , pt , sccm , srcm
    intent (out) frcm
!
    real(rk8) :: dp1 , pt1 , rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1 , n

!
!   INTV0 IS FOR VERTICAL INTERPOLATION OF TRACER WHERE THE TRACER HAS
!   THE SAME VERTICAL ORDERING OF REGCM.
!   THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
!   IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!
    pt1 = pt/psccm
    do i = 1 , ni
      do j = 1 , nj
        dp1 = psrcm(i,j)/psccm
        do n = 1 , krcm
          sc = srcm(n)*dp1 + pt1
          k1 = 0
          do k = 1 , kccm
            if ( sc > sccm(k) ) then
              k1 = k
            end if
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

  subroutine intv1(frcm,fccm,psrcm,srcm,sccm,pt,ni,nj,krcm,kccm)
    implicit none
!
    integer(ik4) :: kccm , krcm , ni , nj
    real(rk8) :: pt
    real(rk8) , dimension(ni,nj,kccm) :: fccm
    real(rk8) , dimension(ni,nj,krcm) :: frcm
    real(rk8) , dimension(ni,nj) :: psrcm
    real(rk8) , dimension(kccm) :: sccm
    real(rk8) , dimension(krcm) :: srcm
    intent (in) fccm , kccm , krcm , ni , nj , psrcm , pt , sccm , srcm
    intent (out) frcm
!
    real(rk8) :: dp1 , pt1 , rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1 , n

!
!   INTV1 IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
!   HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
!   IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
! 
    pt1 = pt/psccm
    do i = 1 , ni
      do j = 1 , nj
        dp1 = psrcm(i,j)/psccm
        do n = 1 , krcm
          sc = srcm(n)*dp1 + pt1
          k1 = 0
          do k = 1 , kccm
            if ( sc > sccm(k) ) then
              k1 = k
            end if
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
  subroutine intv2(frcm,fccm,psrcm,srcm,sccm,pt,ni,nj,krcm,kccm)
    implicit none
!
    integer(ik4) :: kccm , krcm , ni , nj
    real(rk8) :: pt
    real(rk8) , dimension(ni,nj,kccm) :: fccm
    real(rk8) , dimension(ni,nj,krcm) :: frcm
    real(rk8) , dimension(ni,nj) :: psrcm
    real(rk8) , dimension(kccm) :: sccm
    real(rk8) , dimension(krcm) :: srcm
    intent (in) fccm , kccm , krcm , ni , nj , psrcm , pt , sccm , srcm
    intent (out) frcm
!
    real(rk8) :: a1 , dp1 , pt1 , rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1 , n
!
!   INTV2 IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!   LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!   THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!   WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!   CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
!   THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!   TWO EXTREME TEMPERATUES IN THE LAYER.
!
    pt1 = pt/psccm
    do i = 1 , ni
      do j = 1 , nj
        dp1 = psrcm(i,j)/psccm
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
  subroutine intv3(fsccm,fccm,psrccm,sccm,ptop,ni,nj,kccm)
    implicit none
!
    integer(ik4) :: kccm , ni , nj
    real(rk8) :: ptop
    real(rk8) , dimension(ni,nj,kccm) :: fccm
    real(rk8) , dimension(ni,nj) :: fsccm , psrccm
    real(rk8) , dimension(kccm) :: sccm
    intent (in) fccm , kccm , ni , nj , psrccm , ptop , sccm
    intent (out) fsccm
!
    real(rk8) :: a1 , rc , rc1 , sc
    integer(ik4) :: i , j , k , k1 , kp1
!
!   INTV3 IS FOR VERTICAL INTERPOLATION OF TSCCM.  THE INTERPOLATION
!   IS LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!   THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!   WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!   CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
!   THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!   TWO EXTREME TEMPERATUES IN THE LAYER.
!
    do i = 1 , ni
      do j = 1 , nj
        sc = (psrccm(i,j)+ptop)*d_r100
        k1 = 0
        do k = 1 , kccm - 1
          if ( sc <= sccm(k+1) .and. sc >= sccm(k) ) k1 = k
        end do
        !If the surface is below the GCM's lowest level,
        !then extrapolate temperature
        if (sc > sccm(kccm) ) then
          a1 = rgas2*log(sc/sccm(kccm))
          fsccm(i,j) = fccm(i,j,kccm+1-kccm)*(b1-a1)/(b1+a1)
        !Otherwise, interpolate the surface temperature between
        !the two adjacent GCM levels
        else if ( k1 == 0 ) then
          fsccm(i,j) = fccm(i,j,1)
        else
          kp1 = k1 + 1
          rc = log(sccm(k1)/sc)/log(sccm(k1)/sccm(kp1))
          rc1 = d_one - rc
          fsccm(i,j) = rc1*fccm(i,j,kccm+1-k1)+rc*fccm(i,j,kccm+1-kp1)
        end if
      end do
    end do
  end subroutine intv3
!
end module mod_vertint

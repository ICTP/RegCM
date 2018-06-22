module mod_vertint

  implicit none

  private

  ! Standard Gravity (m/sec**2) 3rd CGPM
  real(4) , parameter :: egrav = 9.80665
  ! Gas constant for dry air in Joules/kg/K
  real(4) , parameter :: rgas = 287.05823
  ! Standard atmosphere ICAO 1993
  real(4) , parameter :: lrate = 0.00649

  real(4) , parameter :: bltop = 0.960

  real(4) , parameter :: rglrog = rgas*lrate/egrav

  public :: intlin , intlog , intlin_z

  contains

  subroutine intlin(im,jm,km,f,pstar,sig,ptop,p,fp)
    implicit none
    integer , intent(in) :: im , jm , km
    real(8) , intent(in) :: ptop
    real(4) , intent(in) :: p
    real(4) , intent(in) , dimension(km,jm,im) :: f
    real(4) , intent(in) , dimension(jm,im) :: pstar
    real(4) , intent(in) , dimension(km) :: sig
    real(4) , intent(out) , dimension(jm,im) :: fp

    integer :: i , j , k , kx , knx
    real(4) :: sigp , w1 , wp , pt
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
          ! The searched sigma value
          !
          sigp = (p-pt)/(pstar(j,i)-pt)
          !
          ! Over the top or below bottom level
          !
          if ( sigp <= sig(km) ) then
            fp(j,i) = f(km,j,i)
          else if ( sigp >= sig(1) ) then
            fp(j,i) = f(1,j,i)
          else
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
            fp(j,i) = w1*f(kx,j,i) + wp*f(knx,j,i)
          end if
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
          ! The searched sigma value
          !
          sigp = (p-pt)/(pstar(j,i)-pt)
          !
          ! Over the top or below bottom level
          !
          if ( sigp <= sig(1) ) then
            fp(j,i) = f(1,j,i)
          else if ( sigp >= sig(km) ) then
            fp(j,i) = f(km,j,i)
          else
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
            fp(j,i) = w1*f(kx,j,i) + wp*f(knx,j,i)
          end if
        end do
      end do
    end if
  end subroutine intlin

  subroutine intlin_z(im,jm,km,f,hz,sig,z,fz)
    implicit none

    integer , intent(in) :: im , jm , km
    real(4) , intent(in) , dimension(km,jm,im) :: f , hz
    real(4) , intent(out) , dimension(jm,im) :: fz
    real(4) , intent(in) , dimension(km) :: sig
    real(4) , intent(in) :: z
!
    integer :: i , j , k , kx , knx
    real(4) :: w1 , wz
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN Z.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    !
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
          ! Over the top or below bottom level
          !
          if ( z <= hz(km,j,i) ) then
            fz(j,i) = f(km,j,i)
          else if ( z >= hz(1,j,i) ) then
            fz(j,i) = f(1,j,i)
          else
            !
            ! Search k level below the requested one
            !
            kx = 0
            do k = 1 , km-1
              if ( z > hz(k,j,i) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx + 1
            wz = (z-hz(kx,j,i))/(hz(knx,j,i)-hz(kx,j,i))
            w1 = 1.0 - wz
            fz(j,i) = w1*f(kx,j,i) + wz*f(knx,j,i)
          end if
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
          ! Over the top or below bottom level
          !
          if ( z <= hz(1,j,i) ) then
            fz(j,i) = f(1,j,i)
          else if ( z >= hz(km,j,i) ) then
            fz(j,i) = f(km,j,i)
          else
            !
            ! Search k level below the requested one
            !
            kx = km + 1
            do k = km , 2 , -1
              if ( z > hz(k,j,i) ) exit
              kx = k
            end do
            !
            ! This is the above level
            !
            knx = kx - 1
            wz = (z-hz(kx,j,i))/(hz(knx,j,i)-hz(kx,j,i))
            w1 = 1.0 - wz
            fz(j,i) = w1*f(kx,j,i) + wz*f(knx,j,i)
          end if
        end do
      end do
    end if
  end subroutine intlin_z

  subroutine intlog(im,jm,km,f,pstar,sig,ptop,p,fp)
    implicit none
    integer , intent(in) :: im , jm , km
    real(8) , intent(in) :: ptop
    real(4) , intent(in) , dimension(km,jm,im) :: f
    real(4) , intent(out) , dimension(jm,im) :: fp
    real(4) , intent(in) :: p
    real(4) , intent(in) , dimension(jm,im) :: pstar
    real(4) , intent(in) , dimension(km) :: sig

    real(4) :: sigp , w1 , wp , pt
    integer :: i , j , k , kx , knx , kbc
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
          ! The searched sigma value
          !
          sigp = (p-pt)/(pstar(j,i)-pt)
          !
          ! Extrapolation
          !
          if ( sigp > 1.0 ) then
            fp(j,i) = f(kbc,j,i)*exp(rglrog*log(sigp/sig(kbc)))
          !
          ! Over the top or below bottom level
          !
          else if ( sigp <= sig(km) ) then
            fp(j,i) = f(km,j,i)
          else if ( sigp >= sig(1) ) then
            fp(j,i) = f(1,j,i)
          else
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
            fp(j,i) = w1*f(kx,j,i) + wp*f(knx,j,i)
          end if
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
          ! The searched sigma value
          !
          sigp = (p-pt)/(pstar(j,i)-pt)
          !
          ! Extrapolation
          !
          if ( sigp > 1.0 ) then
            fp(j,i) = f(kbc,j,i)*exp(rglrog*log(sigp/sig(kbc)))
          !
          ! Over the top or below bottom level
          !
          else if ( sigp <= sig(1) ) then
            fp(j,i) = f(1,j,i)
          else if ( sigp >= sig(km) ) then
            fp(j,i) = f(km,j,i)
          else
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
            fp(j,i) = w1*f(kx,j,i) + wp*f(knx,j,i)
          end if
        end do
      end do
    end if
  end subroutine intlog

end module mod_vertint
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

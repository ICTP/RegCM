module mod_vertint

  implicit none

  private

  ! Standard Gravity (m/sec**2) 3rd CGPM
  real(4) , parameter :: egrav = 9.80665
  ! Gas constant for dry air in Joules/kg/K
  real(4) , parameter :: rgas = 287.05823
  ! Standard atmosphere ICAO 1993
  real(4) , parameter :: lrate = 0.00649

  real(4) , parameter :: rglrog = rgas*lrate/egrav

  real(4) , parameter :: bltop = 0.960

  public :: intlin_hy , intlog_hy
  public :: intlin_nonhy , intlog_nonhy

  contains

  subroutine intlin_nonhy(km,jm,im,f,mp,p,fp)
    implicit none
    integer , intent(in) :: im , jm , km
    real(4) , intent(in) :: p
    real(4) , intent(in) , dimension(km,jm,im) :: f
    real(4) , intent(in) , dimension(km,jm,im) :: mp
    real(4) , intent(out) , dimension(jm,im) :: fp

    integer :: i , j , k , knx
    integer :: kx = 0
    real(4) :: w1 , wp , sp , tp , bp
    real(4) , dimension(km) :: spp
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    ! HERE WE ASSUME PRESSURE ARE TOP TO BOTTOM ON KM
    !
    ! Loop over points
    !
    do j = 1 , jm
      do i = 1 , im
        !
        ! Over the top or below bottom level
        !
        tp = mp(1,j,i)
        bp = mp(km,j,i)-tp
        do k = 1 , km
          spp(k) = (mp(k,j,i)-tp) / bp
        end do
        sp = (p-tp) / bp
        if ( sp <= spp(1) ) then
          fp(j,i) = f(1,j,i)
        else if ( sp >= spp(km) ) then
          fp(j,i) = f(km,j,i)
        else
          !
          ! Search p level below
          !
          do k = 2 , km
            kx = k
            if ( spp(k) > sp ) exit
          end do
          knx = kx - 1
          wp = (spp(kx)-sp)/(spp(kx)-spp(knx))
          w1 = 1.0 - wp
          fp(j,i) = w1*f(kx,j,i) + wp*f(knx,j,i)
        end if
      end do
    end do
  end subroutine intlin_nonhy

  subroutine intlog_nonhy(km,jm,im,f,mp,p,fp)
    implicit none
    integer , intent(in) :: im , jm , km
    real(4) , intent(in) , dimension(km,jm,im) :: f
    real(4) , intent(in) , dimension(km,jm,im) :: mp
    real(4) , intent(out) , dimension(jm,im) :: fp
    real(4) , intent(in) :: p

    real(4) :: w1 , wp , sp , bp , tp
    real(4) , dimension(km) :: spp
    integer :: i , j , k , knx
    integer :: kx = 0
    !
    ! INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
    ! LINEAR IN LOG P.
    !
    do j = 1 , jm
      do i = 1 , im
        !
        ! Over the top or below bottom level
        !
        tp = mp(1,j,i)
        bp = mp(km,j,i)-tp
        do k = 1 , km
          spp(k) = (mp(k,j,i) - tp) / bp
        end do
        sp = (p - tp) / bp
        if ( sp <= spp(1) ) then
          fp(j,i) = f(1,j,i)
        else if ( sp > spp(km) ) then
          !
          ! Extrapolate to surface
          !
          fp(j,i) = 0.5*(f(km,j,i)+f(km-1,j,i))*exp(rglrog*log(sp/spp(km)))
        else
          !
          ! Search p level below the requested one
          !
          do k = 2 , km
            kx = k
            if ( spp(k) > sp ) exit
          end do
          knx = kx - 1
          wp = log(spp(kx)/sp)/log(spp(kx)/spp(knx))
          w1 = 1.0 - wp
          fp(j,i) = w1*f(kx,j,i) + wp*f(knx,j,i)
        end if
      end do
    end do
  end subroutine intlog_nonhy

  subroutine intlin_hy(km,jm,im,f,pstar,sig,ptop,p,fp)
    implicit none
    integer , intent(in) :: im , jm , km
    real(4) , intent(in) :: p
    real(8) , intent(in) :: ptop
    real(4) , intent(in) , dimension(km,jm,im) :: f
    real(4) , intent(in) , dimension(jm,im) :: pstar
    real(4) , intent(in) , dimension(km) :: sig
    real(4) , intent(out) , dimension(jm,im) :: fp

    integer :: i , j , k , knx
    integer :: kx = 0
    real(4) :: sigp , w1 , wp , ptp
    !
    ! INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
    ! HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
    ! IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    ! HERE WE ASSUME PRESSURE ARE TOP TO BOTTOM ON KM
    !
    ! Loop over points
    !
    ptp = real(ptop) * 100.0
    do j = 1 , jm
      do i = 1 , im
        !
        ! The searched sigma value
        !
        sigp = (p-ptp)/(pstar(j,i)-ptp)
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
          do k = 2 , km
            kx = k
            if ( sig(k) > sigp ) exit
          end do
          knx = kx - 1
          wp = (sig(kx)-sigp)/(sig(kx)-sig(knx))
          w1 = 1.0 - wp
          fp(j,i) = w1*f(kx,j,i) + wp*f(knx,j,i)
        end if
      end do
    end do
  end subroutine intlin_hy

  subroutine intlog_hy(km,jm,im,f,pstar,sig,ptop,p,fp)
    implicit none
    integer , intent(in) :: im , jm , km
    real(4) , intent(in) , dimension(km,jm,im) :: f
    real(4) , intent(out) , dimension(jm,im) :: fp
    real(8) , intent(in) :: ptop
    real(4) , intent(in) :: p
    real(4) , intent(in) , dimension(jm,im) :: pstar
    real(4) , intent(in) , dimension(km) :: sig

    real(4) :: sigp , w1 , wp , ptp
    integer :: i , j , k , knx
    integer :: kx = 0
    !
    ! INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
    ! LINEAR IN LOG P.
    !
    ! Loop over points
    !
    ptp = real(ptop) * 100.0
    do j = 1 , jm
      do i = 1 , im
        !
        ! Find boundary layer Top
        !
        ! The searched sigma value
        !
        sigp = (p-ptp)/(pstar(j,i)-ptp)
        !
        ! Over the top or below bottom level
        !
        if ( sigp <= sig(1) ) then
          fp(j,i) = f(1,j,i)
        else if ( sigp > sig(km) ) then
          fp(j,i) = 0.5*(f(km,j,i)+f(km-1,j,i))*exp(rglrog*log(sigp/sig(km)))
        else
          !
          ! Search k level above
          !
          do k = 2 , km
            kx = k
            if ( sig(k) > sigp ) exit
          end do
          knx = kx - 1
          wp = log(sig(kx)/sigp)/log(sig(kx)/sig(knx))
          w1 = 1.0 - wp
          fp(j,i) = w1*f(kx,j,i) + wp*f(knx,j,i)
        end if
      end do
    end do
  end subroutine intlog_hy

end module mod_vertint
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

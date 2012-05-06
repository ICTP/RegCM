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
  use mod_realkinds
  use mod_stdio
  use mod_constants
  use mod_message
!
  private

  real(sp) , parameter :: rgas2 = real(rgas/2.0D0)
  real(sp) , parameter :: rglr = real(rgas*lrate)
  real(sp) , parameter :: b1 = -real(egrav/lrate)
  real(sp) , parameter :: rbltop = real(bltop)
  real(sp) , parameter :: psccm = 100.0
!
  integer , parameter :: maxnlev = 100

  public :: intlin , intlin_o , intgtb , intlog , intlog_o
  public :: intpsn , intv0 , intv1 , intv2 , intv3
!
  contains

  subroutine intlin(fp,f,ps,p3d,im,jm,km,p,kp)
  implicit none
!
  integer :: im , jm , km , kp
  real(sp) , dimension(im,jm,km) :: f , p3d
  real(sp) , dimension(im,jm,kp) :: fp
  real(sp) , dimension(kp) :: p
  real(sp) , dimension(im,jm) :: ps
  intent (in) f , im , jm , km , kp , p , p3d , ps
  intent (out) fp
!
  integer :: i , j , k , k1 , kp1 , n
  real(sp) , dimension(maxnlev) :: sig
  real(sp) :: sigp , w1 , wp
!
!     INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
!     HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
!     IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!
  if ( km > maxnlev ) then
    call fatal(__FILE__,__LINE__,'Unrecoverable error, increase maxnlev')
  end if
  do j = 1 , jm
    do i = 1 , im
      if ( ps(i,j) > -9995.0 ) then
        do k = 1 , km
          sig(k) = p3d(i,j,k)/ps(i,j)
        end do
        do n = 1 , kp
          sigp = p(n)/ps(i,j)
          k1 = 0
          do k = 1 , km
            if ( sigp > sig(k) ) k1 = k
          end do
          if ( sigp <= sig(1) ) then
            fp(i,j,n) = f(i,j,1)
          else if ( (sigp > sig(1)) .and. (sigp < sig(km)) ) then
            kp1 = k1 + 1
            wp = (sigp-sig(k1))/(sig(kp1)-sig(k1))
            w1 = 1. - wp
            fp(i,j,n) = w1*f(i,j,k1) + wp*f(i,j,kp1)
          else if ( sigp >= sig(km) ) then
            fp(i,j,n) = f(i,j,km)
          else
          end if
        end do
      else
        do n = 1 , kp
          fp(i,j,n) = -9999.0
        end do
      end if
    end do
  end do
  end subroutine intlin
!
!-----------------------------------------------------------------------
!
  subroutine intlin_o(fp,f,pstar,sig,ptop,im,jm,km,p,kp)
  implicit none
!
  integer :: im , jm , km , kp
  real(dp) :: ptop
  real(sp) , dimension(im,jm,km) :: f
  real(sp) , dimension(im,jm,kp) :: fp
  real(sp) , dimension(kp) :: p
  real(sp) , dimension(im,jm) :: pstar
  real(sp) , dimension(km) :: sig
  intent (in) f , im , jm , km , kp , p , pstar , ptop , sig
  intent (out) fp
!
  integer :: i , j , k , k1 , kp1 , n
  real(sp) :: sigp , w1 , wp
!
!     INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
!     HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
!     IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!
  do j = 1 , jm
    do i = 1 , im
      do n = 1 , kp
        sigp = (p(n)-real(ptop))/(pstar(i,j)-real(ptop))
        k1 = 0
        do k = 1 , km
          if ( sigp > sig(k) ) k1 = k
        end do
        if ( sigp <= sig(1) ) then
          fp(i,j,n) = f(i,j,1)
        else if ( (sigp > sig(1)) .and. (sigp < sig(km)) ) then
          kp1 = k1 + 1
          wp = (sigp-sig(k1))/(sig(kp1)-sig(k1))
          w1 = 1. - wp
          fp(i,j,n) = w1*f(i,j,k1) + wp*f(i,j,kp1)
        else if ( sigp >= sig(km) ) then
          fp(i,j,n) = f(i,j,km)
        else
        end if
      end do
    end do
  end do
  end subroutine intlin_o
!
!-----------------------------------------------------------------------
!
  subroutine intgtb(pa,za,tlayer,zrcm,tp,zp,sccm,ni,nj,nlev1)
  implicit none
!
  integer :: ni , nj , nlev1
  real(sp) , dimension(ni,nj) :: pa , tlayer , za , zrcm
  real(sp) , dimension(nlev1) :: sccm
  real(sp) , dimension(ni,nj,nlev1) :: tp , zp
  intent (in) ni , nj , nlev1 , sccm , tp , zp , zrcm
  intent (out) pa , za
  intent (inout) tlayer
!
  integer :: i , j , k , kb , kt
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
        tlayer(i,j) = (tp(i,j,nlev1+1-kt)+tlayer(i,j))/2.
        za(i,j) = zp(i,j,nlev1+1-kt)
        pa(i,j) = 100.*sccm(kt)
      else
        tlayer(i,j) = tp(i,j,1)
        za(i,j) = zp(i,j,1)
        pa(i,j) = 100.
      end if
 
    end do
  end do
 
!     PRINT *, 'ZRCM, ZP(6)   = ', ZRCM(5,5), ZP(5,5,NLEV1+1-6)
!     PRINT *, '      TP(6)   = ',            TP(5,5,NLEV1+1-6)
!     PRINT *, 'TLAYER, ZA, PA = ', TLAYER(5,5), ZA(5,5), PA(5,5)
 
  end subroutine intgtb
!
!-----------------------------------------------------------------------
!
  subroutine intlog(fp,f,ps,p3d,im,jm,km,p,kp)
  implicit none
!
  integer :: im , jm , km , kp
  real(sp) , dimension(im,jm,km) :: f , p3d
  real(sp) , dimension(im,jm,kp) :: fp
  real(sp) , dimension(kp) :: p
  real(sp) , dimension(im,jm) :: ps
  intent (in) f , im , jm , km , kp , p , p3d , ps
  intent (out) fp
!
  real(sp) :: sigp , w1 , wp
  integer :: i , j , k , k1 , kp1 , kbc , n
  real(sp) , dimension(maxnlev) :: sig
!
!     INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!     LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!     THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!     WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!     CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
!     THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!     TWO EXTREME TEMPERATURES IN THE LAYER.
 
!
!**   FIND FIRST SIGMA LEVEL ABOVE BOUNDARY LAYER (LESS THAN SIG = BLTOP)
  if ( km > maxnlev ) then
    call fatal(__FILE__,__LINE__,'Unrecoverable error, increase maxnlev')
  end if
  do j = 1 , jm
    do i = 1 , im
      if ( ps(i,j) > -9995.0 ) then
        kbc = 1
        do k = 1 , km
          sig(k) = p3d(i,j,k)/ps(i,j)
          if ( sig(k) < rbltop ) kbc = k
        end do
        do n = 1 , kp
          sigp = p(n)/ps(i,j)
          k1 = 0
          do k = 1 , km
            if ( sigp > sig(k) ) k1 = k
          end do
          if ( sigp <= sig(1) ) then
            fp(i,j,n) = f(i,j,1)
          else if ( (sigp > sig(1)) .and. (sigp < sig(km)) ) then
            kp1 = k1 + 1
            wp = log(sigp/sig(k1))/log(sig(kp1)/sig(k1))
            w1 = 1. - wp
            fp(i,j,n) = w1*f(i,j,k1) + wp*f(i,j,kp1)
          else if ( (sigp >= sig(km)) .and. (sigp <= 1.) ) then
            fp(i,j,n) = f(i,j,km)
          else if ( sigp > 1. ) then
            fp(i,j,n) = f(i,j,kbc)*exp(+rglr*log(sigp/sig(kbc))*real(regrav))
!               ***** FROM R. ERRICO, SEE ROUTINE HEIGHT *****
          else
          end if
        end do
      else
        do n = 1 , kp
          fp(i,j,n) = -9999.0
        end do
      end if
    end do
  end do
 
  end subroutine intlog
!
!-----------------------------------------------------------------------
!
  subroutine intlog_o(fp,f,pstar,sig,ptop,im,jm,km,p,kp)
  implicit none
!
  integer :: im , jm , km , kp
  real(dp) :: ptop
  real(sp) , dimension(im,jm,km) :: f
  real(sp) , dimension(im,jm,kp) :: fp
  real(sp) , dimension(kp) :: p
  real(sp) , dimension(im,jm) :: pstar
  real(sp) , dimension(km) :: sig
  intent (in) f , im , jm , km , kp , p , pstar , ptop , sig
  intent (out) fp
!
  real(sp) :: sigp , w1 , wp
  integer :: i , j , k , k1 , kp1 , kbc , n
!
!     INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!     LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!     THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!     WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!     CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
!     THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!     TWO EXTREME TEMPERATURES IN THE LAYER.
 
!
!**   FIND FIRST SIGMA LEVEL ABOVE BOUNDARY LAYER (LESS THAN SIG = BLTOP)
  kbc = 1
  do k = 1 , km
    if ( sig(k) < rbltop ) kbc = k
  end do
  do j = 1 , jm
    do i = 1 , im
      do n = 1 , kp
        sigp = (p(n)-real(ptop))/(pstar(i,j)-real(ptop))
        k1 = 0
        do k = 1 , km
          if ( sigp > sig(k) ) k1 = k
        end do
        if ( sigp <= sig(1) ) then
          fp(i,j,n) = f(i,j,1)
        else if ( (sigp > sig(1)) .and. (sigp < sig(km)) ) then
          kp1 = k1 + 1
          wp = log(sigp/sig(k1))/log(sig(kp1)/sig(k1))
          w1 = 1. - wp
          fp(i,j,n) = w1*f(i,j,k1) + wp*f(i,j,kp1)
        else if ( (sigp >= sig(km)) .and. (sigp <= 1.) ) then
          fp(i,j,n) = f(i,j,km)
        else if ( sigp > 1. ) then
          fp(i,j,n) = f(i,j,kbc)*exp(rglr*log(sigp/sig(kbc))*real(regrav))
!             ***** FROM R. ERRICO, SEE ROUTINE HEIGHT *****
        else
        end if
      end do
    end do
  end do
 
  end subroutine intlog_o
!
!-----------------------------------------------------------------------
!
  subroutine intpsn(psrcm,zrcm,pa,za,tlayer,pt,ni,nj)
  implicit none
!
  integer :: ni , nj
  real(dp) :: pt
  real(sp) , dimension(ni,nj) :: pa , psrcm , tlayer , za , zrcm
  intent (in) ni , nj , pa , pt , tlayer , za , zrcm
  intent (out) psrcm
!
  real(sp) :: tb
  integer :: i , j
!
!     EXTRAPOLATE SURFACE PRESSURE FROM CLOSEST PRESSURE LEVEL ABOVE.
!     USE TLAYER CALCULATED IN INTGTB.
!     PSRCM = SURFACE PRESSURE - PTOP
!
  do i = 1 , ni
    do j = 1 , nj
      tb = tlayer(i,j)
      psrcm(i,j) = pa(i,j)*exp(-real(govr)*(zrcm(i,j)-za(i,j))/tb)-real(pt)
    end do
  end do
 
!     PRINT *, 'ZRCM, ZA, PA, PT = ', ZRCM(5,5), ZA(5,5), PA(5,5), PT
!     PRINT *, 'TLAYER(5,5), PSRCM(5,5) = ', TLAYER(5,5), PSRCM(5,5)
 
  end subroutine intpsn
!
!-----------------------------------------------------------------------
!
  subroutine intv0(frcm,fccm,psrcm,srcm,sccm,pt,ni,nj,krcm,kccm)
    implicit none
!
    integer :: kccm , krcm , ni , nj
    real(dp) :: pt
    real(sp) , dimension(ni,nj,kccm) :: fccm
    real(sp) , dimension(ni,nj,krcm) :: frcm
    real(sp) , dimension(ni,nj) :: psrcm
    real(sp) , dimension(kccm) :: sccm
    real(sp) , dimension(krcm) :: srcm
    intent (in) fccm , kccm , krcm , ni , nj , psrcm , pt , sccm , srcm
    intent (out) frcm
!
    real(sp) :: dp1 , pt1 , rc , rc1 , sc
    integer :: i , j , k , k1 , kp1 , n

!
!   INTV0 IS FOR VERTICAL INTERPOLATION OF TRACER WHERE THE TRACER HAS
!   THE SAME VERTICAL ORDERING OF REGCM.
!   THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
!   IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!
    pt1 = real(pt)/psccm
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
            rc1 = 1.0-rc
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
    integer :: kccm , krcm , ni , nj
    real(dp) :: pt
    real(sp) , dimension(ni,nj,kccm) :: fccm
    real(sp) , dimension(ni,nj,krcm) :: frcm
    real(sp) , dimension(ni,nj) :: psrcm
    real(sp) , dimension(kccm) :: sccm
    real(sp) , dimension(krcm) :: srcm
    intent (in) fccm , kccm , krcm , ni , nj , psrcm , pt , sccm , srcm
    intent (out) frcm
!
    real(sp) :: dp1 , pt1 , rc , rc1 , sc
    integer :: i , j , k , k1 , kp1 , n

!
!   INTV1 IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
!   HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
!   IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
! 
    pt1 = real(pt)/psccm
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
            rc1 = 1.0 - rc
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
    integer :: kccm , krcm , ni , nj
    real(dp) :: pt
    real(sp) , dimension(ni,nj,kccm) :: fccm
    real(sp) , dimension(ni,nj,krcm) :: frcm
    real(sp) , dimension(ni,nj) :: psrcm
    real(sp) , dimension(kccm) :: sccm
    real(sp) , dimension(krcm) :: srcm
    intent (in) fccm , kccm , krcm , ni , nj , psrcm , pt , sccm , srcm
    intent (out) frcm
!
    real(sp) :: a1 , dp1 , pt1 , rc , rc1 , sc
    integer :: i , j , k , k1 , kp1 , n
!
!   INTV2 IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!   LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!   THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!   WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!   CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
!   THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!   TWO EXTREME TEMPERATUES IN THE LAYER.
!
    pt1 = real(pt)/psccm
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
            rc1 = 1.0 - rc
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
    integer :: kccm , ni , nj
    real(dp) :: ptop
    real(sp) , dimension(ni,nj,kccm) :: fccm
    real(sp) , dimension(ni,nj) :: fsccm , psrccm
    real(sp) , dimension(kccm) :: sccm
    intent (in) fccm , kccm , ni , nj , psrccm , ptop , sccm
    intent (out) fsccm
!
    real(sp) :: a1 , rc , rc1 , sc
    integer :: i , j , k , k1 , kp1
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
        sc = (psrccm(i,j)+real(ptop))/100.
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
          rc1 = 1.0 - rc
          fsccm(i,j) = rc1*fccm(i,j,kccm+1-k1)+rc*fccm(i,j,kccm+1-kp1)
        end if
      end do
    end do
  end subroutine intv3
!
end module mod_vertint

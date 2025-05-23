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

module mod_hgt

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_message

  implicit none

  private

  integer(ik4), parameter :: maxnlev = 100

  real(rk4), parameter :: srovg = real(rovg,rk4)
  real(rk4), parameter :: slrate = real(lrate,rk4)
  real(rk4), parameter :: segrav = real(egrav,rk4)
  real(rk4), parameter :: srgas = real(rgas,rk4)

  interface htsig_o
    module procedure htsig_o_double
    module procedure htsig_o_single
  end interface htsig_o

  interface height_o
    module procedure height_o_double
    module procedure height_o_double_nonhy
    module procedure height_o_single
    module procedure height_o_single_nonhy
  end interface height_o

  interface nonhydrost
    module procedure nonhydrost_double
    module procedure nonhydrost_single
    module procedure nonhydrost_single_2
  end interface nonhydrost

  interface htsig
    module procedure htsig_1
    module procedure htsig_2
    module procedure htsig_3
  end interface htsig

  public :: hydrost, nonhydrost, mslp2ps
  public :: height, height_o
  public :: htsig, htsig_o
  public :: psig, psig1, mslp, gs_filter
  public :: crc_zeta_from_p, crc_p_from_zeta
  public :: p2pai, pai2p

  contains

  subroutine hydrost(h,t,topo,ps,ptop,sigmah,ni,nj,nk)
    implicit none
    integer(ik4), intent(in) :: ni, nj, nk
    real(rkx), intent(in) :: ptop
    real(rkx), intent(in), dimension(nk) :: sigmah
    real(rkx), intent(in), dimension(ni,nj,nk) :: t
    real(rkx), intent(in), dimension(ni,nj) :: topo, ps
    real(rkx), intent(out), dimension(ni,nj,nk) :: h

    integer(ik4) :: i, j, k
    real(rkx), dimension(nk+1) :: sigmaf
    real(rkx), dimension(nk) :: dsigma
    real(rkx) :: pf, tbar
    !
    ! ROUTINE TO COMPUTE HEIGHT USING THE HYDROSTATIC RELATION.
    ! THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
    ! HEIGHT IS COMPUTED IN THE RCM MODEL.
    !
    sigmaf(1) = d_zero
    sigmaf(nk+1) = d_one
    do k = 2, nk
      sigmaf(k) = d_half*(sigmah(k-1)+sigmah(k))
    end do
    do k = 1, nk
      dsigma(k) = sigmaf(k+1) - sigmaf(k)
    end do
    !
    ! SET BOUNDARY VALUES TO ZERO AT ALL LEVELS SINCE THE HEIGHT IS
    ! DEFINED AT CROSS POINTS AND AT HALF LEVELS.
    !
    do concurrent ( j = 1:nj, i = 1:ni )
      pf = ptop/ps(i,j)
      h(i,j,nk) = topo(i,j) + rovg*t(i,j,nk)*log((d_one+pf)/(sigmah(nk)+pf))
      do k = nk - 1, 1, -1
        tbar = (t(i,j,k)*dsigma(k) + &
          t(i,j,k+1)*dsigma(k+1))/(dsigma(k)+dsigma(k+1))
        h(i,j,k) = h(i,j,k+1)+rovg*tbar*log((sigmah(k+1)+pf)/(sigmah(k)+pf))
      end do
    end do
  end subroutine hydrost

  subroutine nonhydrost_double(h,t0,p0,ps,topo,ni,nj,nk)
    implicit none
    integer(ik4), intent(in) :: ni, nj, nk
    real(rk8), intent(in), dimension(ni,nj,nk) :: t0, p0
    real(rk8), intent(in), dimension(ni,nj) :: ps, topo
    real(rk8), intent(out), dimension(ni,nj,nk) :: h
    real(rk8) :: tbar
    integer(ik4) :: i, j, k
    !
    ! ROUTINE TO COMPUTE HEIGHT FOR THE NON-HYDROSTATIC CORE
    ! THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
    ! HEIGHT IS COMPUTED IN THE RCM MODEL.
    !
    do concurrent ( i = 1:ni, j = 1:nj )
      h(i,j,nk) = topo(i,j) + rovg * t0(i,j,nk) * log(ps(i,j)/p0(i,j,nk))
    end do
    do k = nk-1, 1, -1
      do concurrent ( i = 1:ni, j = 1:nj )
        tbar = d_half*(t0(i,j,k)+t0(i,j,k+1))
        h(i,j,k) = h(i,j,k+1) + rovg * tbar * log(p0(i,j,k+1)/p0(i,j,k))
      end do
    end do
  end subroutine nonhydrost_double

  subroutine nonhydrost_single(h,t0,p0,ps,topo,ni,nj,nk)
    implicit none
    integer(ik4), intent(in) :: ni, nj, nk
    real(rk4), intent(in), dimension(ni,nj,nk) :: t0, p0
    real(rk4), intent(in), dimension(ni,nj) :: ps, topo
    real(rk4), intent(out), dimension(ni,nj,nk) :: h
    real(rk4) :: tbar
    integer(ik4) :: i, j, k
    !
    ! ROUTINE TO COMPUTE HEIGHT FOR THE NON-HYDROSTATIC CORE
    ! THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
    ! HEIGHT IS COMPUTED IN THE RCM MODEL.
    !
    do concurrent ( i = 1:ni, j = 1:nj )
      h(i,j,nk) = topo(i,j) + srovg * t0(i,j,nk) * log(ps(i,j)/p0(i,j,nk))
    end do
    do k = nk-1, 1, -1
      do concurrent ( i = 1:ni, j = 1:nj )
        tbar = 0.5*(t0(i,j,k)+t0(i,j,k+1))
        h(i,j,k) = h(i,j,k+1) + srovg * tbar * log(p0(i,j,k+1)/p0(i,j,k))
      end do
    end do
  end subroutine nonhydrost_single

  subroutine nonhydrost_single_2(h,t0,p0,ps,topo,i1,i2,j1,j2,nk)
    implicit none
    integer(ik4), intent(in) :: i1, i2, j1, j2, nk
    real(rk8), intent(in), dimension(i1:i2,j1:j2,nk) :: t0, p0
    real(rk8), intent(in), dimension(i1:i2,j1:j2) :: ps, topo
    real(rk8), intent(out), dimension(i1:i2,j1:j2,nk) :: h
    real(rk8) :: tbar
    integer(ik4) :: i, j, k
    !
    ! ROUTINE TO COMPUTE HEIGHT FOR THE NON-HYDROSTATIC CORE
    ! THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
    ! HEIGHT IS COMPUTED IN THE RCM MODEL.
    !
    do concurrent ( i = i1:i2, j = j1:j2 )
      h(i,j,nk) = topo(i,j) + srovg * t0(i,j,nk) * log(ps(i,j)/p0(i,j,nk))
    end do
    do k = nk-1, 1, -1
      do concurrent ( i = i1:i2, j = j1:j2 )
        tbar = d_half*(t0(i,j,k)+t0(i,j,k+1))
        h(i,j,k) = h(i,j,k+1) + srovg * tbar * log(p0(i,j,k+1)/p0(i,j,k))
      end do
    end do
  end subroutine nonhydrost_single_2

  subroutine height(hp,h,t,ps,p3d,ht,im,jm,km,p,kp)
    implicit none
    integer(ik4), intent(in) :: im, jm, km, kp
    real(rkx), intent(in), dimension(im,jm,km) :: h, p3d, t
    real(rkx), intent(in), dimension(im,jm) :: ht, ps
    real(rkx), intent(in), dimension(kp) :: p
    real(rkx), intent(out), dimension(im,jm,kp) :: hp

    real(rkx) :: psfc, temp, wb, wt, pt, pb
    integer(ik4) :: i, j, k, kb, kt, n, ipb, ipt, ipi
    real(rkx), dimension(km) :: psig
    !
    !  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
    !     ON INPUT:
    !        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
    !        PS = SURFACE PRESSURE
    !        PTOP = MODEL TOP PRESSURE.
    !        SIG = SIGMA LEVELS.
    !        P = PRESSURE LEVELS DESIRED.
    !     ON OUTPUT:
    !        ALL FIELDS EXCEPT H ARE UNCHANGED.
    !        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
    !
    !  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
    !  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
    !     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
    !     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
    !     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
    !     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
    !      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
    !
    if ( p3d(1,1,1) < p3d(1,1,km) ) then
      ipt = 1
      ipb = km
      ipi = +1
    else
      ipt = km
      ipb = 1
      ipi = -1
    end if
#ifdef STDPAR
    do concurrent ( i = 1:im, j = 1:jm ) local(psig)
#else
    do j = 1, jm
      do i = 1, im
#endif
        psfc = ps(i,j)
        if ( psfc > -9995.0_rkx ) then
          do k = 1, km
            psig(k) = p3d(i,j,k)
          end do
          pt = psig(ipt)
          pb = psig(ipb)
          do n = 1, kp
            if ( p(n) < pt ) then
              temp = t(i,j,ipt)
              hp(i,j,n) = h(i,j,ipt) + rovg*temp*log(psig(ipt)/p(n))
            else if ( p(n) > psfc ) then
              temp = 0.5 * ( t(i,j,ipb) + t(i,j,ipb-ipi))
              hp(i,j,n) = ht(i,j) + &
                    (temp/lrate)*(d_one-exp(+rovg*lrate*log(p(n)/psfc)))
            else if ( p(n) >= pb ) then
              temp = t(i,j,ipb)
              hp(i,j,n) = ht(i,j) + rovg*temp*log(psfc/p(n))
            else
              kt = ipt
              do k = ipt, ipb, ipi
                if ( psig(k) > p(n) ) exit
                kt = k
              end do
              kb = kt + ipi
              wt = log(psig(kb)/p(n))/log(psig(kb)/psig(kt))
              wb = 1.0_rkx - wt
              temp = wt*t(i,j,kt) + wb*t(i,j,kb)
              temp = (temp+t(i,j,kb))/d_two
              hp(i,j,n) = h(i,j,kb) + rovg*temp*log(psig(kb)/p(n))
            end if
          end do
        else
          do n = 1, kp
            hp(i,j,n) = -9999.0_rkx
          end do
        end if
#ifndef STDPAR
      end do
#endif
    end do
  end subroutine height
!
!-----------------------------------------------------------------------
!
  subroutine height_o_double(hp,h,t,ps,ht,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4), intent(in) :: im, jm, km, kp
    real(rkx), intent(in) :: ptop
    real(rk8), intent(in), dimension(im,jm,km) :: h, t
    real(rk8), intent(out), dimension(im,jm,kp) :: hp
    real(rk8), intent(in), dimension(im,jm) :: ht, ps
    real(rk8), intent(in), dimension(kp) :: p
    real(rk8), intent(in), dimension(km) :: sig

    real(rk8) :: psfc, temp, wb, wt
    integer(ik4) :: i, j, k, kb, kt, n
    real(rk8), dimension(km) :: psig
    !
    !  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
    !     ON INPUT:
    !        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
    !        PSTAR = SURFACE PRESSURE - MODEL TOP PRESSURE.
    !        SIG = SIGMA LEVELS.
    !        P = PRESSURE LEVELS DESIRED.
    !     ON OUTPUT:
    !        ALL FIELDS EXCEPT H ARE UNCHANGED.
    !        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
    !
    !  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
    !  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
    !     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
    !     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
    !     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
    !     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
    !      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
    !
#ifdef STDPAR
    do concurrent ( i = 1:im, j = 1:jm ) local(psig)
#else
    do j = 1, jm
      do i = 1, im
#endif
        do k = 1, km
          psig(k) = sig(k)*(ps(i,j)-ptop) + ptop
        end do
        psfc = ps(i,j)
        do n = 1, kp
          if ( p(n) <= psig(1) ) then
            temp = t(i,j,1)
            hp(i,j,n) = h(i,j,1) + rovg*temp*log(psig(1)/p(n))
          else if ( p(n) > psfc ) then
            temp = 0.5 * (t(i,j,km) + t(i,j,km-1))
            hp(i,j,n) = ht(i,j) + &
                  (temp/lrate)*(d_one-exp(+rovg*lrate*log(p(n)/psfc)))
          else if ( p(n) >= psig(km) ) then
            temp = t(i,j,km)
            hp(i,j,n) = ht(i,j) + rovg*temp*log(psfc/p(n))
          else
            kt = 1
            do k = 1, km
              if ( psig(k) > p(n) ) exit
              kt = k
            end do
            kb = kt + 1
            wt = log(psig(kb)/p(n))/log(psig(kb)/psig(kt))
            wb = 1.0_rk8 - wt
            temp = wt*t(i,j,kt) + wb*t(i,j,kb)
            temp = (temp+t(i,j,kb))/d_two
            hp(i,j,n) = h(i,j,kb) + rovg*temp*log(psig(kb)/p(n))
          end if
        end do
#ifndef STDPAR
      end do
#endif
    end do
  end subroutine height_o_double

  subroutine height_o_double_nonhy(hp,h,t,ps,ht,p3,im,jm,km,p,kp)
    implicit none
    integer(ik4), intent(in) :: im, jm, km, kp
    real(rk8), intent(in), dimension(im,jm,km) :: h, t, p3
    real(rk8), intent(out), dimension(im,jm,kp) :: hp
    real(rk8), intent(in), dimension(im,jm) :: ht, ps
    real(rk8), intent(in), dimension(kp) :: p

    real(rk8) :: psfc, temp, wb, wt
    integer(ik4) :: i, j, k, kb, kt, n
    !
    !  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
    !     ON INPUT:
    !        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
    !        PSTAR = SURFACE PRESSURE - MODEL TOP PRESSURE.
    !        SIG = SIGMA LEVELS.
    !        P = PRESSURE LEVELS DESIRED.
    !     ON OUTPUT:
    !        ALL FIELDS EXCEPT H ARE UNCHANGED.
    !        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
    !
    !  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
    !  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
    !     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
    !     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
    !     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
    !     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
    !      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
    !
    do concurrent ( i = 1:im, j = 1:jm )
      psfc = ps(i,j)
      do n = 1, kp
        if ( p(n) <= p3(i,j,1) ) then
          temp = t(i,j,1)
          hp(i,j,n) = h(i,j,1) + rovg*temp*log(p3(i,j,1)/p(n))
        else if ( p(n) > psfc ) then
          temp = 0.5 * (t(i,j,km) + t(i,j,km-1))
          hp(i,j,n) = ht(i,j) + &
                (temp/lrate)*(d_one-exp(+rovg*lrate*log(p(n)/psfc)))
        else if ( p(n) >= p3(i,j,km) ) then
          temp = t(i,j,km)
          hp(i,j,n) = ht(i,j) + rovg*temp*log(psfc/p(n))
        else
          kt = 1
          do k = 1, km
            if ( p3(i,j,k) > p(n) ) exit
            kt = k
          end do
          kb = kt + 1
          wt = log(p3(i,j,kb)/p(n))/log(p3(i,j,kb)/p3(i,j,kt))
          wb = 1.0_rk8 - wt
          temp = wt*t(i,j,kt) + wb*t(i,j,kb)
          temp = (temp+t(i,j,kb))/d_two
          hp(i,j,n) = h(i,j,kb) + rovg*temp*log(p3(i,j,kb)/p(n))
        end if
      end do
    end do
  end subroutine height_o_double_nonhy

  subroutine height_o_single(hp,h,t,ps,ht,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4), intent(in) :: im, jm, km, kp
    real(rkx), intent(in) :: ptop
    real(rk4), intent(in), dimension(im,jm,km) :: h, t
    real(rk4), intent(out), dimension(im,jm,kp) :: hp
    real(rk4), intent(in), dimension(im,jm) :: ht, ps
    real(rk4), intent(in), dimension(kp) :: p
    real(rk4), intent(in), dimension(km) :: sig

    real(rk4) :: psfc, temp, wb, wt, ptp
    integer(ik4) :: i, j, k, kb, kt, n
    real(rk4), dimension(km) :: psig
    !
    !  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
    !     ON INPUT:
    !        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
    !        PSTAR = SURFACE PRESSURE - MODEL TOP PRESSURE.
    !        SIG = SIGMA LEVELS.
    !        P = PRESSURE LEVELS DESIRED.
    !     ON OUTPUT:
    !        ALL FIELDS EXCEPT H ARE UNCHANGED.
    !        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
    !
    !  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
    !  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
    !     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
    !     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
    !     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
    !     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
    !      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
    !
    ptp = real(ptop)
#ifdef STDPAR
    do concurrent ( i = 1:im, j = 1:jm ) local(psig)
#else
    do j = 1, jm
      do i = 1, im
#endif
        do k = 1, km
          psig(k) = sig(k)*(ps(i,j)-ptp) + ptp
        end do
        psfc = ps(i,j)
        do n = 1, kp
          if ( p(n) <= psig(1) ) then
            temp = t(i,j,1)
            hp(i,j,n) = h(i,j,1) + srovg*temp*log(psig(1)/p(n))
          else if ( p(n) > psfc ) then
            temp = 0.5 * (t(i,j,km) + t(i,j,km-1))
            hp(i,j,n) = ht(i,j) + &
                  (temp/slrate)*(1.0-exp(+srovg*slrate*log(p(n)/psfc)))
          else if ( p(n) >= psig(km) ) then
            temp = t(i,j,km)
            hp(i,j,n) = ht(i,j) + srovg*temp*log(psfc/p(n))
          else
            kt = 1
            do k = 1, km
              if ( psig(k) > p(n) ) exit
              kt = k
            end do
            kb = kt + 1
            wt = log(psig(kb)/p(n))/log(psig(kb)/psig(kt))
            wb = 1.0 - wt
            temp = wt*t(i,j,kt) + wb*t(i,j,kb)
            temp = (temp+t(i,j,kb))/2.0
            hp(i,j,n) = h(i,j,kb) + srovg*temp*log(psig(kb)/p(n))
          end if
        end do
#ifndef STDPAR
      end do
#endif
    end do
  end subroutine height_o_single

  subroutine height_o_single_nonhy(hp,h,t,ps,ht,p3,im,jm,km,p,kp)
    implicit none
    integer(ik4), intent(in) :: im, jm, km, kp
    real(rk4), intent(in), dimension(im,jm,km) :: h, t, p3
    real(rk4), intent(out), dimension(im,jm,kp) :: hp
    real(rk4), intent(in), dimension(im,jm) :: ht, ps
    real(rk4), intent(in), dimension(kp) :: p

    real(rk4) :: psfc, temp, wb, wt
    integer(ik4) :: i, j, k, kb, kt, n
    !
    !  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
    !     ON INPUT:
    !        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
    !        PSTAR = SURFACE PRESSURE - MODEL TOP PRESSURE.
    !        SIG = SIGMA LEVELS.
    !        P = PRESSURE LEVELS DESIRED.
    !     ON OUTPUT:
    !        ALL FIELDS EXCEPT H ARE UNCHANGED.
    !        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
    !
    !  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
    !  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
    !     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
    !     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
    !     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
    !     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
    !      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
    !
    do concurrent ( i = 1:im, j = 1:jm )
      psfc = ps(i,j)
      do n = 1, kp
        if ( p(n) > psfc ) then
          temp = 0.5 * (t(i,j,km) + t(i,j,km-1))
          hp(i,j,n) = ht(i,j) + &
                (temp/slrate)*(1.0-exp(+srovg*slrate*log(p(n)/psfc)))
        else if ( p(n) <= p3(i,j,1) ) then
          temp = t(i,j,1)
          hp(i,j,n) = h(i,j,1) + srovg*temp*log(p3(i,j,1)/p(n))
        else if ( p(n) >= p3(i,j,km) ) then
          temp = t(i,j,km)
          hp(i,j,n) = ht(i,j) + srovg*temp*log(psfc/p(n))
        else
          kt = 1
          do k = 1, km
            if ( p3(i,j,k) > p(n) ) exit
            kt = k
          end do
          kb = kt + 1
          wt = log(p3(i,j,kb)/p(n))/log(p3(i,j,kb)/p3(i,j,kt))
          wb = 1.0 - wt
          temp = wt*t(i,j,kt) + wb*t(i,j,kb)
          temp = (temp+t(i,j,kb))/2.0
          hp(i,j,n) = h(i,j,kb) + srovg*temp*log(p3(i,j,kb)/p(n))
        end if
      end do
    end do
  end subroutine height_o_single_nonhy

  subroutine htsig_1(t,h,p3d,ps,ht,im,jm,km)
    implicit none
    integer(ik4), intent(in) :: im, jm, km
    real(rkx), intent(in), dimension(im,jm,km) :: p3d, t
    real(rkx), intent(out), dimension(im,jm,km) :: h
    real(rkx), intent(in), dimension(im,jm) :: ht, ps

    real(rkx) :: tbar
    integer(ik4) :: i, j, k

    if ( p3d(1,1,km) > p3d(1,1,1) ) then
      do concurrent ( i = 1:im, j = 1:jm )
        if ( ps(i,j) > -9995.0_rkx ) then
          h(i,j,km) = ht(i,j) + rovg*t(i,j,km)*log(ps(i,j)/p3d(i,j,km))
        else
          h(i,j,km) = -9999.0_rkx
        end if
      end do
      do k = km - 1, 1, -1
        do concurrent ( i = 1:im, j = 1:jm )
          if ( h(i,j,k+1) > -9995.0_rkx ) then
            tbar = d_half*(t(i,j,k)+t(i,j,k+1))
            h(i,j,k) = h(i,j,k+1)+rovg*tbar*log(p3d(i,j,k+1)/p3d(i,j,k))
          else
            h(i,j,k) = -9999.0_rkx
          end if
        end do
      end do
    else
      do concurrent ( i = 1:im, j = 1:jm )
        if ( ps(i,j) > -9995.0_rkx ) then
          h(i,j,1) = ht(i,j) + rovg*t(i,j,km)*log(ps(i,j)/p3d(i,j,1))
        else
          h(i,j,1) = -9999.0_rkx
        end if
      end do
      do k = 2, km
        do concurrent ( i = 1:im, j = 1:jm )
          if ( h(i,j,k-1) > -9995.0_rkx ) then
            tbar = d_half*(t(i,j,k)+t(i,j,k-1))
            h(i,j,k) = h(i,j,k-1)+rovg*tbar*log(p3d(i,j,k-1)/p3d(i,j,k))
          else
            h(i,j,k) = -9999.0_rkx
          end if
        end do
      end do
    end if
  end subroutine htsig_1

  subroutine htsig_2(t,h,pstar,ht,sig,ptop,i1,i2,j1,j2,km)
    implicit none
    integer(ik4), intent(in) :: i1, i2, j1, j2, km
    real(rk8), intent(in) :: ptop
    real(rk8), intent(in), dimension(i1:i2,j1:j2,km) :: t
    real(rk8), intent(in), dimension(i1:i2,j1:j2) :: ht, pstar
    real(rk8), intent(in), dimension(km) :: sig
    real(rk8), intent(out), dimension(i1:i2,j1:j2,km) :: h
    real(rk8) :: tbar
    integer(ik4) :: i, j, k
    do concurrent ( i = i1:i2, j = j1:j2 )
      h(i,j,km) = ht(i,j) + &
          rovg*t(i,j,km)*log(pstar(i,j)/((pstar(i,j)-ptop)*sig(km)+ptop))
    end do
    do k = km - 1, 1, -1
      do concurrent ( i = i1:i2, j = j1:j2 )
        tbar = d_half*(t(i,j,k)+t(i,j,k+1))
        h(i,j,k) = h(i,j,k+1) + &
          rovg*tbar*log(((pstar(i,j)-ptop)*sig(k+1)+ptop)/ &
                        ((pstar(i,j)-ptop)*sig(k)+ptop))
      end do
    end do
  end subroutine htsig_2

  subroutine htsig_3(z,t,p,q,ps,ht)
    implicit none
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: t, p, q
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: ps, ht
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(inout) :: z
    integer(ik4) :: ni, nj, nk
    integer(ik4) :: i, j, k
    real(rkx) :: tv0, h0, p0, tv1
    ni = size(z,1)
    nj = size(z,2)
    nk = size(z,3)
    if ( p(1,1,1) > p(1,1,nk) ) then
      do concurrent( i = 1:ni, j = 1:nj )
        p0 = ps(i,j)
        h0 = ht(i,j)
        tv0 = t(i,j,1) * (1.0_rkx + ep1 * q(i,j,1))
        z(i,j,1) = h0 + rovg * tv0 * log(p0/p(i,j,1))
      end do
      do k = 2, nk
        do concurrent( i = 1:ni, j = 1:nj )
          tv0 = t(i,j,k-1) * (1.0_rkx + ep1 * q(i,j,k-1))
          tv1 = t(i,j,k) * (1.0_rkx + ep1 * q(i,j,k))
          z(i,j,k) = z(i,j,k-1) + &
            rovg * 0.5_rkx * (tv0+tv1) * log(p(i,j,k-1)/p(i,j,k))
        end do
      end do
    else
      do concurrent( i = 1:ni, j = 1:nj )
        tv0 = t(i,j,nk) * (1.0_rkx + ep1 * q(i,j,nk))
        z(i,j,nk) = ht(i,j) + rovg * tv0 * log(ps(i,j)/p(i,j,nk))
      end do
      do k = nk-1, 1, -1
        do concurrent( i = 1:ni, j = 1:nj )
          tv0 = t(i,j,k+1) * (1.0_rkx + ep1 * q(i,j,k+1))
          tv1 = t(i,j,k) * (1.0_rkx + ep1 * q(i,j,k))
          z(i,j,k) = z(i,j,k+1) + &
            rovg * 0.5_rkx * (tv0+tv1) * log(p(i,j,k+1)/p(i,j,k))
        end do
      end do
    end if
  end subroutine htsig_3

  subroutine htsig_o_double(t,h,pstar,ht,sig,ptop,im,jm,km)
    implicit none
    integer(ik4), intent(in) :: im, jm, km
    real(rkx), intent(in) :: ptop
    real(rk8), intent(in), dimension(im,jm,km) :: t
    real(rk8), intent(in), dimension(im,jm) :: ht, pstar
    real(rk8), intent(in), dimension(km) :: sig
    real(rk8), intent(out), dimension(im,jm,km) :: h
    real(rk8) :: tbar
    integer(ik4) :: i, j, k
    do concurrent( i = 1:im, j = 1:jm )
      h(i,j,km) = ht(i,j) + &
          rovg*t(i,j,km)*log(pstar(i,j)/((pstar(i,j)-ptop)*sig(km)+ptop))
    end do
    do k = km - 1, 1, -1
      do concurrent( i = 1:im, j = 1:jm )
        tbar = d_half*(t(i,j,k)+t(i,j,k+1))
        h(i,j,k) = h(i,j,k+1) + &
          rovg*tbar*log(((pstar(i,j)-ptop)*sig(k+1)+ptop)/ &
                         ((pstar(i,j)-ptop)*sig(k)+ptop))
      end do
    end do
  end subroutine htsig_o_double

  subroutine htsig_o_single(t,h,pstar,ht,sig,ptop,im,jm,km)
    implicit none
    integer(ik4), intent(in) :: im, jm, km
    real(rkx), intent(in) :: ptop
    real(rk4), intent(in), dimension(im,jm,km) :: t
    real(rk4), intent(in), dimension(im,jm) :: ht, pstar
    real(rk4), intent(in), dimension(km) :: sig
    real(rk4), intent(out), dimension(im,jm,km) :: h
    real(rk4) :: tbar, rpt
    integer(ik4) :: i, j, k
    rpt = real(ptop)
    do concurrent( i = 1:im, j = 1:jm )
      h(i,j,km) = ht(i,j) + &
          srovg*t(i,j,km)*log(pstar(i,j)/((pstar(i,j)-rpt)*sig(km)+rpt))
    end do
    do k = km - 1, 1, -1
      do concurrent( i = 1:im, j = 1:jm )
        tbar = 0.5*(t(i,j,k)+t(i,j,k+1))
        h(i,j,k) = h(i,j,k+1) + &
          srovg*tbar*log(((pstar(i,j)-rpt)*sig(k+1)+rpt)/ &
                         ((pstar(i,j)-rpt)*sig(k)+rpt))
      end do
    end do
  end subroutine htsig_o_single

  pure elemental real(rkx) function p2pai(p) result(pai)
    implicit none
    real(rkx), intent(in) :: p
    pai = (p/p00)**rovcp
  end function p2pai

  pure elemental real(rkx) function pai2p(pai) result(p)
    implicit none
    real(rkx), intent(in) :: pai
    p = (pai**cpovr) * p00
  end function pai2p

  subroutine mslp2ps(h,t,slp,ht,ps,im,jm,km)
    implicit none
    integer(ik4), intent(in) :: im, jm, km
    real(rkx), dimension(im,jm,km), intent(in) :: h, t
    real(rkx), dimension(im,jm), intent(in) :: ht, slp
    real(rkx), dimension(im,jm), intent(out) :: ps
    integer(ik4) :: kbc, i, j , k
    real(rkx) :: tsfc
    real(rkx), parameter :: blhgt = 1000.0_rkx

    do j = 1, jm
      kbc = km
      do i = 1, im
        do k = 1, km
          if ( h(i,j,k) > blhgt ) then
            kbc = k
            exit
          end if
        end do
        tsfc = t(i,j,kbc)-lrate*(h(i,j,kbc)-ht(i,j))
        ps(i,j) = slp(i,j) / exp(-egrav/(rgas*lrate)* &
                             log(d_one-ht(i,j)*lrate/tsfc))
      end do
    end do
  end subroutine mslp2ps

  subroutine mslp(t,ps,ht,slp,im,jm,kz)
    implicit none
    integer(ik4), intent(in) :: im, jm, kz
    real(rk4), dimension(im,jm,kz), intent(in) :: t
    real(rk4), dimension(im,jm), intent(in) :: ht, ps
    real(rk4), dimension(im,jm), intent(out) :: slp
    integer(ik4) :: i, j
    real(rk4) :: tstar, hstar, alpha, sraval

    ! Follow Kallen 1996
    alpha = real(lrate*rgas/egrav,rk4)
    do concurrent ( i = 1:im, j = 1:jm )
      tstar = t(i,j,kz)
      if ( tstar < 255.0 ) then
        tstar = (tstar+255.0)*0.5
      else if ( tstar > 290.5 ) then
        tstar = 290.5 + (0.005*(tstar-290.5))**2
      end if
      hstar = ht(i,j)*segrav/(srgas*tstar)
      sraval = 0.5*alpha*hstar
      slp(i,j) = ps(i,j) * exp(hstar*(1.0 - sraval + (sraval*sraval)/3.0))
    end do
  end subroutine mslp

  subroutine psig(t,h,p3d,ps,ht,im,jm,km)
    implicit none
    integer(ik4), intent(in) :: im, jm, km
    real(rkx), intent(in), dimension(im,jm,km) :: h, t
    real(rkx), intent(in), dimension(im,jm) :: ht, ps
    real(rkx), intent(out), dimension(im,jm,km) :: p3d
    real(rkx) :: tbar
    integer(ik4) :: i, j, k

    do concurrent ( i = 1:im, j = 1:jm )
      p3d(i,j,1) = ps(i,j)*exp(-(h(i,j,1)-ht(i,j))/rovg/t(i,j,km))
    end do
    do k = 2, km
      do concurrent ( i = 1:im, j = 1:jm )
        tbar = d_half*(t(i,j,k)+t(i,j,k-1))
        p3d(i,j,k) = p3d(i,j,k-1)*exp(-(h(i,j,k)-h(i,j,k-1))/rovg/tbar)
      end do
    end do
  end subroutine psig

  subroutine psig1(t,h,p3d,ps,ht,im,jm,km)
    implicit none
    integer(ik4), intent(in) :: im, jm, km
    real(rkx), intent(in), dimension(im,jm,km) :: h, t
    real(rkx), intent(in), dimension(im,jm) :: ht, ps
    real(rkx), intent(out), dimension(im,jm,km) :: p3d
    real(rkx) :: tbar
    integer(ik4) :: i, j, k

    do concurrent ( i = 1:im, j = 1:jm )
      p3d(i,j,1) = ps(i,j)*exp(-(h(i,j,1)-ht(i,j))/rovg/t(i,j,km))
    end do
    do k = 2, km
      do concurrent ( i = 1:im, j = 1:jm )
        tbar = t(i,j,k)
        p3d(i,j,k) = p3d(i,j,k-1)*exp(-(h(i,j,k)-h(i,j,k-1))/rovg/tbar)
      end do
    end do
  end subroutine psig1

  ! Gauss Siedel Filtering
  subroutine gs_filter(v,vm,im,jm)
    implicit none
    integer(ik4), intent(in) :: im, jm
    real(rk4), dimension(im,jm), intent(in) :: vm
    real(rk4), dimension(im,jm), intent(inout) :: v
    integer(ik4) :: i, j, n
    integer(ik4), parameter :: niter = 20
    real(rk4), dimension(im,jm) :: v1, mask
    real(rk4) :: mval
    v1(:,1) = v(:,1)
    v1(1,:) = v(1,:)
    v1(im,:) = v(im,:)
    v1(:,jm) = v(:,jm)
    mask(:,:) = 0.0
    mval = 0.5*(maxval(vm)-minval(vm))
    do concurrent ( i = 2:im-1, j = 2:jm-1 )
      mask(i,j) = (vm(i-1,j)+vm(i+1,j)+vm(i,j-1) + &
                   vm(i,j+1)-4.0*vm(i,j))/mval
    end do
    do n = 1, niter
      do j = 2, jm-1
        do i = 2, im-1
          v1(i,j) = 0.25*(v1(i-1,j)+v(i+1,j)+v1(i,j-1)+v(i,j+1)-mask(i,j))
        end do
      end do
      v(:,:) = v1(:,:)
    end do
  end subroutine gs_filter

  pure real(rkx) function crc_zeta_from_p(p) result(z)
    implicit none
    real(rkx), intent(in) :: p
    real(rk8), parameter :: a1 = 44330.8_rk8
    real(rk8), parameter :: a2 = 4946.54_rk8
    real(rk8), parameter :: ex = 0.1902632_rk8
    real(rk8) :: p8, z8
    p8 = real(p, rk8)
    z8 = a1 - a2 * p8**ex
    z = real(z8,rkx)
  end function crc_zeta_from_p

  pure real(rkx) function crc_p_from_zeta(z) result(p)
    implicit none
    real(rkx), intent(in) :: z
    real(rk8) :: p8, z8
    z8 = z
    p8 = 100.0_rk8*((44331.514_rk8-z)/11880.516_rk8)**(1.0_rk8/0.1902632_rk8)
    p = real(p8,rkx)
  end function crc_p_from_zeta

end module mod_hgt
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

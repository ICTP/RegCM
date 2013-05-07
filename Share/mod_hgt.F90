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

module mod_hgt

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_message

  private

  integer(ik4) , parameter :: maxnlev = 100

  real(rk4) , parameter :: srovg = real(rovg)
  real(rk4) , parameter :: slrate = real(lrate)
  real(rk4) , parameter :: segrav = real(egrav)
  real(rk4) , parameter :: srgas = real(rgas)

  interface htsig_o
    module procedure htsig_o_double
    module procedure htsig_o_single
  end interface htsig_o

  interface height_o
    module procedure height_o_double
    module procedure height_o_single
  end interface height_o

  public :: hydrost , mslp2ps
  public :: height , height_o
  public :: htsig , htsig_o
  public :: psig , mslp , gs_filter

  contains

  subroutine hydrost(h,t,topo,ps,ptop,sigmah,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk8) , intent(in) :: ptop
    real(rk8) , intent(in) , dimension(nk) :: sigmah
    real(rk8) , intent(in) , dimension(ni,nj,nk) :: t
    real(rk8) , intent(out) , dimension(ni,nj,nk) :: h
    real(rk8) , intent(in) , dimension(ni,nj) :: topo , ps
!
    integer(ik4) :: i , j , k
    real(rk8) , dimension(nk+1) :: sigmaf
    real(rk8) , dimension(nk) :: dsigma
    real(rk8) :: pf , tbar
    !
    ! ROUTINE TO COMPUTE HEIGHT USING THE HYDROSTATIC RELATION.
    ! THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
    ! HEIGHT IS COMPUTED IN THE RCM MODEL.
    !
    sigmaf(1) = d_zero
    sigmaf(nk+1) = d_one
    do k = 2 , nk
      sigmaf(k) = d_half*(sigmah(k-1)+sigmah(k))
    end do
    do k = 1 , nk
      dsigma(k) = sigmaf(k+1) - sigmaf(k)
    end do
    !
    ! SET BOUNDARY VALUES TO ZERO AT ALL LEVELS SINCE THE HEIGHT IS
    ! DEFINED AT CROSS POINTS AND AT HALF LEVELS.
    !
    do i = 1 , ni
      do j = 1 , nj
        pf = ptop/ps(i,j)
        h(i,j,nk) = topo(i,j) + rovg*t(i,j,nk)*dlog((d_one+pf)/(sigmah(nk)+pf))
        do k = nk - 1 , 1 , -1
          tbar = (t(i,j,k)*dsigma(k) + &
            t(i,j,k+1)*dsigma(k+1))/(dsigma(k)+dsigma(k+1))
          h(i,j,k) = h(i,j,k+1)+rovg*tbar*dlog((sigmah(k+1)+pf)/(sigmah(k)+pf))
        end do
      end do
    end do
  end subroutine hydrost

  subroutine height(hp,h,t,ps,p3d,ht,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , intent(in) , dimension(im,jm,km) :: h , p3d , t
    real(rk8) , intent(in) , dimension(im,jm) :: ht , ps
    real(rk8) , intent(in) , dimension(kp) :: p
    real(rk8) , intent(out) , dimension(im,jm,kp) :: hp
!
    real(rk8) :: psfc , temp , wb , wt
    integer(ik4) :: i , j , k , kb , kbc , kt , n
    real(rk8) , dimension(km+1) :: psig
    real(rk8) , dimension(km) :: sig
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
    do j = 1 , jm
      do i = 1 , im
        psfc = ps(i,j)
        kbc = 1
        if ( psfc > -9995.0D0 ) then
          do k = 1 , km
            sig(k) = p3d(i,j,k)/ps(i,j)
            if ( sig(k) < bltop ) kbc = k
            psig(k) = p3d(i,j,k)
          end do
          do n = 1 , kp
            kt = 1
            do k = 1 , km
              if ( psig(k) < p(n) ) kt = k
            end do
            kb = kt + 1
            if ( p(n) <= psig(1) ) then
              temp = t(i,j,1)
              hp(i,j,n) = h(i,j,1) + rovg*temp*dlog(psig(1)/p(n))
            else if ( (p(n) > psig(1)) .and. (p(n) < psig(km)) ) then
              wt = dlog(psig(kb)/p(n))/dlog(psig(kb)/psig(kt))
              wb = dlog(p(n)/psig(kt))/dlog(psig(kb)/psig(kt))
              temp = wt*t(i,j,kt) + wb*t(i,j,kb)
              temp = (temp+t(i,j,kb))/d_two
              hp(i,j,n) = h(i,j,kb) + rovg*temp*dlog(psig(kb)/p(n))
            else if ( (p(n) >= psig(km)) .and. (p(n) <= psfc) ) then
              temp = t(i,j,km)
              hp(i,j,n) = ht(i,j) + rovg*temp*dlog(psfc/p(n))
            else if ( p(n) > psfc ) then
              temp = t(i,j,kbc) + lrate*(h(i,j,kbc)-ht(i,j))
              hp(i,j,n) = ht(i,j)+ &
                      (temp/lrate)*(d_one-dexp(+rovg*lrate*dlog(p(n)/psfc)))
            end if
          end do
        else
          do n = 1 , kp
            hp(i,j,n) = -9999.0D0
          end do
        end if
      end do
    end do
  end subroutine height
!
!-----------------------------------------------------------------------
!
  subroutine height_o_double(hp,h,t,pstar,ht,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , intent(in) :: ptop
    real(rk8) , intent(in) , dimension(im,jm,km) :: h , t
    real(rk8) , intent(out) , dimension(im,jm,kp) :: hp
    real(rk8) , intent(in) , dimension(im,jm) :: ht , pstar
    real(rk8) , intent(in) , dimension(kp) :: p
    real(rk8) , intent(in) , dimension(km) :: sig
!
    real(rk8) :: psfc , temp , wb , wt
    integer(ik4) :: i , j , k , kb , kbc , kt , n
    real(rk8) , dimension(km) :: psig
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
    kbc = 1
    do k = 1 , km
      if ( sig(k) < bltop ) kbc = k
    end do
    do j = 1 , jm
      do i = 1 , im
        do k = 1 , km
          psig(k) = sig(k)*(pstar(i,j)-ptop) + ptop
        end do
        psfc = pstar(i,j)
        do n = 1 , kp
          kt = 1
          do k = 1 , km
            if ( psig(k) < p(n) ) kt = k
          end do
          kb = kt + 1
          if ( p(n) <= psig(1) ) then
            temp = t(i,j,1)
            hp(i,j,n) = h(i,j,1) + rovg*temp*dlog(psig(1)/p(n))
          else if ( (p(n) > psig(1)) .and. (p(n) < psig(km)) ) then
            wt = dlog(psig(kb)/p(n))/dlog(psig(kb)/psig(kt))
            wb = dlog(p(n)/psig(kt))/dlog(psig(kb)/psig(kt))
            temp = wt*t(i,j,kt) + wb*t(i,j,kb)
            temp = (temp+t(i,j,kb))/d_two
            hp(i,j,n) = h(i,j,kb) + rovg*temp*dlog(psig(kb)/p(n))
          else if ( (p(n) >= psig(km)) .and. (p(n) <= psfc) ) then
            temp = t(i,j,km)
            hp(i,j,n) = ht(i,j) + rovg*temp*dlog(psfc/p(n))
          else if ( p(n) > psfc ) then
            temp = t(i,j,kbc) + lrate*(h(i,j,kbc)-ht(i,j))
            hp(i,j,n) = ht(i,j) + &
                  (temp/lrate)*(d_one-dexp(+rovg*lrate*dlog(p(n)/psfc)))
          end if
        end do
      end do
    end do
  end subroutine height_o_double

  subroutine height_o_single(hp,h,t,pstar,ht,sig,ptop,im,jm,km,p,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km , kp
    real(rk8) , intent(in) :: ptop
    real(rk4) , intent(in) , dimension(im,jm,km) :: h , t
    real(rk4) , intent(out) , dimension(im,jm,kp) :: hp
    real(rk4) , intent(in) , dimension(im,jm) :: ht , pstar
    real(rk4) , intent(in) , dimension(kp) :: p
    real(rk4) , intent(in) , dimension(km) :: sig
!
    real(rk4) :: psfc , temp , wb , wt , ptp
    integer(ik4) :: i , j , k , kb , kbc , kt , n
    real(rk4) , dimension(km) :: psig
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
    kbc = 1
    do k = 1 , km
      if ( sig(k) < bltop ) kbc = k
    end do
    do j = 1 , jm
      do i = 1 , im
        do k = 1 , km
          psig(k) = sig(k)*(pstar(i,j)-ptp) + ptp
        end do
        psfc = pstar(i,j)
        do n = 1 , kp
          kt = 1
          do k = 1 , km
            if ( psig(k) < p(n) ) kt = k
          end do
          kb = kt + 1
          if ( p(n) <= psig(1) ) then
            temp = t(i,j,1)
            hp(i,j,n) = h(i,j,1) + srovg*temp*log(psig(1)/p(n))
          else if ( (p(n) > psig(1)) .and. (p(n) < psig(km)) ) then
            wt = log(psig(kb)/p(n))/log(psig(kb)/psig(kt))
            wb = log(p(n)/psig(kt))/log(psig(kb)/psig(kt))
            temp = wt*t(i,j,kt) + wb*t(i,j,kb)
            temp = (temp+t(i,j,kb))/2.0
            hp(i,j,n) = h(i,j,kb) + srovg*temp*log(psig(kb)/p(n))
          else if ( (p(n) >= psig(km)) .and. (p(n) <= psfc) ) then
            temp = t(i,j,km)
            hp(i,j,n) = ht(i,j) + srovg*temp*log(psfc/p(n))
          else if ( p(n) > psfc ) then
            temp = t(i,j,kbc) + slrate*(h(i,j,kbc)-ht(i,j))
            hp(i,j,n) = ht(i,j) + &
                  (temp/slrate)*(1.0-exp(+srovg*slrate*log(p(n)/psfc)))
          end if
        end do
      end do
    end do
  end subroutine height_o_single
!
!-----------------------------------------------------------------------
!
  subroutine htsig(t,h,p3d,ps,ht,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rk8) , intent(in) , dimension(im,jm,km) :: p3d , t
    real(rk8) , intent(out) , dimension(im,jm,km) :: h
    real(rk8) , intent(in) , dimension(im,jm) :: ht , ps
!
    real(rk8) :: tbar
    integer(ik4) :: i , j , k
!
    do j = 1 , jm
      do i = 1 , im
        if ( ps(i,j) > -9995.0D0 ) then
          h(i,j,km) = ht(i,j) + rovg*t(i,j,km)*dlog(ps(i,j)/p3d(i,j,km))
        else
          h(i,j,km) = -9999.0D0
        end if
      end do
    end do
    do k = km - 1 , 1 , -1
      do j = 1 , jm
        do i = 1 , im
          if ( h(i,j,k+1) > -9995.0D0 ) then
            tbar = d_half*(t(i,j,k)+t(i,j,k+1))
            h(i,j,k) = h(i,j,k+1)+rovg*tbar*dlog(p3d(i,j,k+1)/p3d(i,j,k))
          else
            h(i,j,k) = -9999.0D0
          end if
        end do
      end do
    end do
  end subroutine htsig
!
!-----------------------------------------------------------------------
!
  subroutine htsig_o_double(t,h,pstar,ht,sig,ptop,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rk8) , intent(in) :: ptop
    real(rk8) , intent(in) , dimension(im,jm,km) :: t
    real(rk8) , intent(in) , dimension(im,jm) :: ht , pstar
    real(rk8) , intent(in) , dimension(km) :: sig
    real(rk8) , intent(out) , dimension(im,jm,km) :: h
    real(rk8) :: tbar
    integer(ik4) :: i , j , k
    do j = 1 , jm
      do i = 1 , im
        h(i,j,km) = ht(i,j) + &
          rovg*t(i,j,km)*dlog(pstar(i,j)/((pstar(i,j)-ptop)*sig(km)+ptop))
      end do
    end do
    do k = km - 1 , 1 , -1
      do j = 1 , jm
        do i = 1 , im
          tbar = d_half*(t(i,j,k)+t(i,j,k+1))
          h(i,j,k) = h(i,j,k+1) + &
            rovg*tbar*dlog(((pstar(i,j)-ptop)*sig(k+1)+ptop)/ &
                           ((pstar(i,j)-ptop)*sig(k)+ptop))
        end do
      end do
    end do
  end subroutine htsig_o_double
!
  subroutine htsig_o_single(t,h,pstar,ht,sig,ptop,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rk8) , intent(in) :: ptop
    real(rk4) , intent(in) , dimension(im,jm,km) :: t
    real(rk4) , intent(in) , dimension(im,jm) :: ht , pstar
    real(rk4) , intent(in) , dimension(km) :: sig
    real(rk4) , intent(out) , dimension(im,jm,km) :: h
    real(rk4) :: tbar , rpt
    integer(ik4) :: i , j , k
    rpt = real(ptop)
    do j = 1 , jm
      do i = 1 , im
        h(i,j,km) = ht(i,j) + &
          srovg*t(i,j,km)*log(pstar(i,j)/((pstar(i,j)-rpt)*sig(km)+rpt))
      end do
    end do
    do k = km - 1 , 1 , -1
      do j = 1 , jm
        do i = 1 , im
          tbar = 0.5*(t(i,j,k)+t(i,j,k+1))
          h(i,j,k) = h(i,j,k+1) + &
            srovg*tbar*log(((pstar(i,j)-rpt)*sig(k+1)+rpt)/ &
                           ((pstar(i,j)-rpt)*sig(k)+rpt))
        end do
      end do
    end do
  end subroutine htsig_o_single

  subroutine mslp2ps(h,t,slp,ht,ps,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rk8) , dimension(im,jm,km) , intent(in) :: h , t
    real(rk8) , dimension(im,jm) , intent(in) :: ht , slp
    real(rk8) , dimension(im,jm) , intent(out) :: ps
    integer(ik4) :: kbc , i , j  , k
    real(rk8) :: tsfc
    real(rk8) , parameter :: blhgt = 1000.0D0

    do j = 1 , jm
      kbc = km
      do i = 1 , im
        do k = 1 , km
          if ( h(i,j,k) > blhgt ) then
            kbc = k
            exit
          end if
        end do
        tsfc = t(i,j,kbc)-lrate*(h(i,j,kbc)-ht(i,j))
        ps(i,j) = slp(i,j) / dexp(-egrav/(rgas*lrate)* &
                             dlog(d_one-ht(i,j)*lrate/tsfc))
      end do
    end do
  end subroutine mslp2ps
!
  subroutine mslp(t,ps,ht,slp,im,jm,kz)
    implicit none
    integer(ik4) , intent(in) :: im , jm , kz
    real(rk4) , dimension(im,jm,kz) , intent(in) :: t
    real(rk4) , dimension(im,jm) , intent(in) :: ht , ps
    real(rk4) , dimension(im,jm) , intent(out) :: slp
    integer(ik4) :: i , j
    real(rk4) :: tstar , hstar , alpha , sraval
!
    ! Follow Kallen 1996
    alpha = real(lrate*rgas/egrav)
    do j = 1 , jm
      do i = 1 , im
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
    end do
  end subroutine mslp

  subroutine psig(t,h,p3d,ps,ht,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rk8) , intent(in) , dimension(im,jm,km) :: h , t
    real(rk8) , intent(in) , dimension(im,jm) :: ht , ps
    real(rk8) , intent(out) , dimension(im,jm,km) :: p3d
!
    real(rk8) :: tbar
    integer(ik4) :: i , j , k
!
    do j = 1 , jm
      do i = 1 , im
        p3d(i,j,1) = ps(i,j)*dexp(-(h(i,j,1)-ht(i,j))/rovg/t(i,j,km))
      end do
    end do
    do k = 2 , km
      do j = 1 , jm
        do i = 1 , im
          tbar = d_half*(t(i,j,k)+t(i,j,k-1))
          p3d(i,j,k) = p3d(i,j,k-1)*dexp(-(h(i,j,k)-h(i,j,k-1))/rovg/tbar)
        end do
      end do
    end do
  end subroutine psig
!
  ! Gauss Siedel Filtering
  subroutine gs_filter(v,vm,im,jm)
    implicit none
    integer(ik4) , intent(in) :: im , jm
    real(rk4) , dimension(im,jm) , intent(in) :: vm
    real(rk4) , dimension(im,jm) , intent(inout) :: v
    integer(ik4) :: i , j , n
    integer(ik4) , parameter :: niter = 20
    real(rk4) , dimension(im,jm) :: v1 , mask
    real(rk4) :: mval
    v1(:,1) = v(:,1)
    v1(1,:) = v(1,:)
    v1(im,:) = v(im,:)
    v1(:,jm) = v(:,jm)
    mask(:,:) = 0.0
    mval = 0.5*(maxval(vm)-minval(vm))
    do j = 2 , jm-1
      do i = 2 , im-1
        mask(i,j) = (vm(i-1,j)+vm(i+1,j)+vm(i,j-1) + &
                     vm(i,j+1)-4.0*vm(i,j))/mval
      end do
    end do
    do n = 1 , niter
      do j = 2 , jm-1
        do i = 2 , im-1
          v1(i,j) = 0.25*(v1(i-1,j)+v(i+1,j)+v1(i,j-1)+v(i,j+1)-mask(i,j))
        end do
      end do
      v(:,:) = v1(:,:)
    end do
  end subroutine gs_filter

end module mod_hgt

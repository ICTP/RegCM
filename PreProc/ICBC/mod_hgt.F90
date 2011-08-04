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

  use mod_constants
  use mod_message
  use m_realkinds

  real(sp) , private , parameter :: srovg = real(rovg)
  real(sp) , private , parameter :: slrate = real(lrate)

  integer , private , parameter :: maxnlev = 100
  contains

  subroutine hydrost(h,t,phis,ps,ptop,sigmaf,sigmah,dsigma,ni,nj,nk)
  implicit none
!
  integer :: ni , nj , nk
  real(dp) :: ptop
  real(sp) , dimension(nk) :: dsigma , sigmah
  real(sp) , dimension(ni,nj,nk) :: h , t
  real(sp) , dimension(ni,nj) :: phis , ps
  real(sp) , dimension(nk+1) :: sigmaf
  intent (in) ni , nj , nk , phis , ps , ptop , sigmaf , sigmah , t
  intent (inout) dsigma , h
!
  integer :: i , j , k , k1 , k2
  real(sp) :: pf , tbar , pt
!
!     ROUTINE TO COMPUTE HEIGHT USING THE HYDROSTATIC RELATION.
!     THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
!     HEIGHT IS COMPUTED IN THE RCM MODEL.
!
  pt = real(ptop)
  do k = 1 , nk
    dsigma(k) = sigmaf(k+1) - sigmaf(k)
  end do
!
!     SET BOUNDARY VALUES TO ZERO AT ALL LEVELS SINCE THE HEIGHT IS
!     DEFINED AT CROSS POINTS AND AT HALF LEVELS.
!
  do i = 1 , ni
    do j = 1 , nj
      pf = pt/ps(i,j)
      h(i,j,nk) = phis(i,j) + srovg*t(i,j,nk)*log((1.+pf)/(sigmah(nk)+pf))
      do k2 = 1 , nk - 1
        k = nk - k2
        k1 = k + 1
        tbar = (t(i,j,k)*dsigma(k)+t(i,j,k1)*dsigma(k1))/(dsigma(k)+dsigma(k1))
        h(i,j,k) = h(i,j,k1)+srovg*tbar*log((sigmah(k1)+pf)/(sigmah(k)+pf))
      end do
    end do
  end do
  end subroutine hydrost
!
!-----------------------------------------------------------------------
!
  subroutine height(hp,h,t,ps,p3d,ht,im,jm,km,p,kp)
 
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
  implicit none
!
  integer :: im , jm , km , kp
  real(sp) , dimension(im,jm,km) :: h , p3d , t
  real(sp) , dimension(im,jm,kp) :: hp
  real(sp) , dimension(im,jm) :: ht , ps
  real(sp) , dimension(kp) :: p
  intent (in) h , ht , im , jm , km , kp , p , p3d , ps , t
  intent (out) hp
!
  real(sp) :: psfc , temp , wb , wt
  integer :: i , j , k , kb , kbc , kt , n
  real(sp) , dimension(maxnlev+1) :: psig
  real(sp) , dimension(maxnlev) :: sig
!
  if ( km > maxnlev ) then
    call fatal(__FILE__,__LINE__,'Unrecoverable error, increase maxnlev')
  end if
  do j = 1 , jm
    do i = 1 , im
      psfc = ps(i,j)
      kbc = 1
      if ( psfc > -9995.0 ) then
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
            hp(i,j,n) = h(i,j,1) + srovg*temp*log(psig(1)/p(n))
          else if ( (p(n) > psig(1)) .and. (p(n) < psig(km)) ) then
            wt = log(psig(kb)/p(n))/log(psig(kb)/psig(kt))
            wb = log(p(n)/psig(kt))/log(psig(kb)/psig(kt))
            temp = wt*t(i,j,kt) + wb*t(i,j,kb)
            temp = (temp+t(i,j,kb))/2.
            hp(i,j,n) = h(i,j,kb) + srovg*temp*log(psig(kb)/p(n))
          else if ( (p(n) >= psig(km)) .and. (p(n) <= psfc) ) then
            temp = t(i,j,km)
            hp(i,j,n) = ht(i,j) + srovg*temp*log(psfc/p(n))
          else if ( p(n) > psfc ) then
            temp = t(i,j,kbc) + slrate*(h(i,j,kbc)-ht(i,j))
            hp(i,j,n) = ht(i,j)+ &
                    (temp/slrate)*(1.-exp(+srovg*slrate*log(p(n)/psfc)))
!
          else
          end if
        end do
      else
        do n = 1 , kp
          hp(i,j,n) = -9999.0
        end do
      end if
    end do
  end do
  end subroutine height
!
!-----------------------------------------------------------------------
!
  subroutine height_o(hp,h,t,pstar,ht,sig,ptop,im,jm,km,p,kp)
 
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
  implicit none
!
  integer :: im , jm , km , kp
  real(dp) :: ptop
  real(sp) , dimension(im,jm,km) :: h , t
  real(sp) , dimension(im,jm,kp) :: hp
  real(sp) , dimension(im,jm) :: ht , pstar
  real(sp) , dimension(kp) :: p
  real(sp) , dimension(km) :: sig
  intent (in) h , ht , im , jm , km , kp , p , pstar , ptop , sig , t
  intent (out) hp
!
  real(sp) :: psfc , temp , wb , wt , pt
  integer :: i , j , k , kb , kbc , kt , n
  real(sp) , dimension(100) :: psig
!
  pt = real(ptop)
  kbc = 1
  do k = 1 , km
    if ( sig(k) < bltop ) kbc = k
  end do
!     PRINT *,'FIRST SIGMA LEVEL ABOVE BNDY LAYER:', SIG(KBC)
!
  do j = 1 , jm
    do i = 1 , im
      do k = 1 , km
        psig(k) = sig(k)*(pstar(i,j)-pt) + pt
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
          temp = (temp+t(i,j,kb))/2.
          hp(i,j,n) = h(i,j,kb) + srovg*temp*log(psig(kb)/p(n))
        else if ( (p(n) >= psig(km)) .and. (p(n) <= psfc) ) then
          temp = t(i,j,km)
          hp(i,j,n) = ht(i,j) + srovg*temp*log(psfc/p(n))
        else if ( p(n) > psfc ) then
          temp = t(i,j,kbc) + slrate*(h(i,j,kbc)-ht(i,j))
          hp(i,j,n) = ht(i,j) + &
                (temp/slrate)*(1.-exp(+srovg*slrate*log(p(n)/psfc)))
!
        else
        end if
      end do
    end do
  end do
  end subroutine height_o
!
!-----------------------------------------------------------------------
!
  subroutine htsig(t,h,p3d,ps,ht,im,jm,km)
  implicit none
!
  integer :: im , jm , km
  real(sp) , dimension(im,jm,km) :: h , p3d , t
  real(sp) , dimension(im,jm) :: ht , ps
  intent (in) ht , im , jm , km , p3d , ps , t
  intent (inout) h
!
  real(sp) :: tbar
  integer :: i , j , k
!
  do j = 1 , jm
    do i = 1 , im
      if ( ps(i,j) > -9995.0 ) then
        h(i,j,km) = ht(i,j) + srovg*t(i,j,km)*log(ps(i,j)/p3d(i,j,km))
      else
        h(i,j,km) = -9999.0
      end if
    end do
  end do
  do k = km - 1 , 1 , -1
    do j = 1 , jm
      do i = 1 , im
        if ( h(i,j,k+1) > -9995.0 ) then
          tbar = 0.5*(t(i,j,k)+t(i,j,k+1))
          h(i,j,k) = h(i,j,k+1)+srovg*tbar*log(p3d(i,j,k+1)/p3d(i,j,k))
        else
          h(i,j,k) = -9999.0
        end if
      end do
    end do
  end do
  end subroutine htsig
!
!-----------------------------------------------------------------------
!
  subroutine htsig_o(t,h,pstar,ht,sig,ptop,im,jm,km)
  implicit none
!
  integer :: im , jm , km
  real(dp) :: ptop
  real(sp) , dimension(im,jm,km) :: h , t
  real(sp) , dimension(im,jm) :: ht , pstar
  real(sp) , dimension(km) :: sig
  intent (in) ht , im , jm , km , pstar , ptop , sig , t
  intent (inout) h
!
  real(sp) :: tbar , pt
  integer :: i , j , k
!
  pt = real(ptop)
  do j = 1 , jm
    do i = 1 , im
      h(i,j,km) = ht(i,j) + &
              srovg*t(i,j,km)*log(pstar(i,j)/((pstar(i,j)-pt)*sig(km)+pt))
    end do
  end do
  do k = km - 1 , 1 , -1
    do j = 1 , jm
      do i = 1 , im
        tbar = 0.5*(t(i,j,k)+t(i,j,k+1))
        h(i,j,k) = h(i,j,k+1) + &
                srovg*tbar*log(((pstar(i,j)-pt)*sig(k+1)+pt)/ &
                               ((pstar(i,j)-pt)*sig(k)+pt))
      end do
    end do
  end do
  end subroutine htsig_o
!
end module mod_hgt

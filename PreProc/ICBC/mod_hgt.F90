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

  integer(ik4) , private , parameter :: maxnlev = 100

  contains

  subroutine hydrost(h,t,phis,ps,ptop,sigmaf,sigmah,dsigma,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk8) , intent(in) :: ptop
    real(rk8) , intent(in) , dimension(nk) :: sigmah
    real(rk8) , intent(inout) , dimension(nk) :: dsigma
    real(rk8) , intent(in) , dimension(ni,nj,nk) :: t
    real(rk8) , intent(inout) , dimension(ni,nj,nk) :: h
    real(rk8) , intent(in) , dimension(ni,nj) :: phis , ps
    real(rk8) , intent(in) , dimension(nk+1) :: sigmaf
!
    integer(ik4) :: i , j , k , k1 , k2
    real(rk8) :: pf , tbar
    !
    ! ROUTINE TO COMPUTE HEIGHT USING THE HYDROSTATIC RELATION.
    ! THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
    ! HEIGHT IS COMPUTED IN THE RCM MODEL.
    !
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
        h(i,j,nk) = phis(i,j) + rovg*t(i,j,nk)*dlog((d_one+pf)/(sigmah(nk)+pf))
        do k2 = 1 , nk - 1
          k = nk - k2
          k1 = k + 1
          tbar = (t(i,j,k)*dsigma(k) + &
            t(i,j,k1)*dsigma(k1))/(dsigma(k)+dsigma(k1))
          h(i,j,k) = h(i,j,k1)+rovg*tbar*dlog((sigmah(k1)+pf)/(sigmah(k)+pf))
        end do
      end do
    end do
  end subroutine hydrost
!
!-----------------------------------------------------------------------
!
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
  subroutine height_o(hp,h,t,pstar,ht,sig,ptop,im,jm,km,p,kp)
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
  end subroutine height_o
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
  subroutine htsig_o(t,h,pstar,ht,sig,ptop,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rk8) , intent(in) :: ptop
    real(rk8) , intent(in) , dimension(im,jm,km) :: t
    real(rk8) , intent(in) , dimension(im,jm) :: ht , pstar
    real(rk8) , intent(in) , dimension(km) :: sig
    real(rk8) , intent(out) , dimension(im,jm,km) :: h
!
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
  end subroutine htsig_o
!
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
end module mod_hgt

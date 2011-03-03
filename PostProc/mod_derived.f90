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

module mod_derived

  use mod_constants

  contains

  subroutine humid2(t,q,preslv,im,jm,kp)
    implicit none
    integer im,jm,kp
    real    t(im,jm,kp),q(im,jm,kp)
    real    preslv(kp)
    integer i,j,k
    real    hl,satvp,qs
    real , parameter :: qmin = 0.0 ! minimum value of specific humidity
!
!  this routine replaces specific humidity by relative humidity
!  data on sigma levels
!
    do k = 1 , kp
      do j = 1 , jm
        do i = 1 , im
          hl = 597.3-.566*(t(i,j,k)-273.16)
          satvp = 6.11*exp(9.045*hl*(rtzero-1./t(i,j,k)))
          qs = .622*satvp/(preslv(k)-satvp)           ! preslv (hpa)
          if (q(i,j,k).lt.qmin) q(i,j,k) = qmin       ! specified minimum
          q(i,j,k) = q(i,j,k)*qs
        end do
      end do
    end do
  end subroutine humid2

  subroutine humid1(t,q,ps,sigma,im,jm,km,ptop)
    implicit none
    integer im,jm,km
    real    ptop
    real    t(im,jm,km),q(im,jm,km)
    real    ps(im,jm)
    real    sigma(km)
    integer i,j,k
    real    hl,satvp,qs,p
    real , parameter :: qmin = 0.0 ! minimum value of specific humidity
!
!  this routine replaces specific humidity by relative humidity
!  data on sigma levels
!
    do k = 1 , km
      do j = 1 , jm
        do i = 1 , im
          p = sigma(k)*(ps(i,j)-ptop)+ptop
          hl = 597.3-.566*(t(i,j,k)-273.16)           ! latent heat of evap.
          satvp = 6.11*exp(9.045*hl*(rtzero-1./t(i,j,k))) ! saturation vap press.
          qs = .622*satvp/(p-satvp)                   ! sat. mixing ratio
          q(i,j,k) = q(i,j,k)/qs
        end do
      end do
    end do
  end subroutine humid1
!
  subroutine height(hp,h,t,pstar,ht,sig,im,jm,km,p,kp,ptop)

!  height determines the height of pressure levels.
!  on input:
!    h and t are height and temperature on sigma, respectively.
!    pstar = surface pressure - model top pressure.
!    sig = sigma levels.
!    p = pressure levels desired.
!  on output:
!    all fields except h are unchanged.
!  h has height fields at kp pressure levels.
!
!  for upward extrapolation, t is considered to have 0 veritcal deriv.
!  for downward extrapolation, t has lapse rate of lrate (k/km)
!  and extrapolation is done from the lowest sigma level above
!  the boundary layer (top arbitrarily taken at sigma = bltop).
!  equation used is exact solution to hydrostatic relation,
!  gotten from r. errico (also used in slpres routine):
!  z = z0 - (t0/lrate) * (1.-exp(-r*lrate*ln(p/p0)/g))
!
    implicit none
    integer im,jm,km,kp
    real    t(im,jm,km),h(im,jm,km),hp(im,jm,kp)
    real    pstar(im,jm),ht(im,jm)
    real    sig(km),p(kp)
    real    ptop
    real    psig(100)
    integer i,j,k,kbc,n,kt,kb
    real    psfc,temp,wt,wb
!
    do k = 1 , km
      if (sig(k).lt.bltop) then
        kbc = k
      end if
    end do
!
    do j = 1 , jm
      do i = 1 , im
        do k = 1 , km
          psig(k) = sig(k) * (pstar(i,j)-ptop) + ptop
        end do
        psfc = pstar(i,j)
        do n = 1 , kp
          kt = 1
          do k = 1 , km
            if (psig(k).lt.p(n)) kt = k
          end do
          kb = kt + 1
          if (p(n).le.psig(1)) then
            temp = t(i,j,1)
            hp(i,j,n) = h(i,j,1)+rgas*temp*log(psig(1)/p(n))*regrav
          else if((p(n).gt.psig(1)) .and. (p(n).lt.psig(km))) then
            wt = log(psig(kb)/p(n)) / log(psig(kb)/psig(kt))
            wb = log(p(n)/psig(kt)) / log(psig(kb)/psig(kt))
            temp = wt * t(i,j,kt) + wb * t(i,j,kb)
            temp = ( temp + t(i,j,kb) ) / 2.
            hp(i,j,n) = h(i,j,kb)+rgas*temp*log(psig(kb)/p(n))*regrav
          else if ((p(n).ge.psig(km)) .and. (p(n).le.psfc)) then
            temp = t(i,j,km)
            hp(i,j,n) = ht(i,j)+rgas*temp*log(psfc/p(n))*regrav
          else if (p(n).gt.psfc) then
            temp = t(i,j,kbc) - lrate * (h(i,j,kbc)-ht(i,j))
            hp(i,j,n) = ht(i,j)-(temp/lrate)  &
                    * ( 1.-exp(-rgas*lrate*log(p(n)/psfc)*regrav))
          end if
        end do
      end do
    end do
  end subroutine height
!
  subroutine slpres(h,t,pstar,ht,slp,sig,im,jm,km)
    implicit none
    integer im,jm,km
    real    t(im,jm,km),h(im,jm,km),pstar(im,jm),ht(im,jm)
!   real    tg(im,jm)
    real    slp(im,jm)
    real    sig(km)
    integer kbc,i,j,k
    real    tsfc
!
    do k=1,km
      if (sig(k).lt.bltop) then
        kbc=k
      end if
    end do
    do j = 1 , jm
      do i = 1 , im
        tsfc = t(i,j,kbc)-lrate*(h(i,j,kbc)-ht(i,j))
        slp(i,j) = pstar(i,j)  &
            * exp( -egrav/(rgas*lrate)*log(1.-ht(i,j)*lrate/tsfc))
      end do
    end do

!   do j=1,jm
!     do i=1,im
!       slp2(i,j) = pstar(i,j) * &
!             exp( egrav*ht(i,j)/(rgas*0.5*(tg(i,j)+288.15)))
!     end do
!   end do
  end subroutine slpres
!
  subroutine htsig(t,h,pstar,ht,sig,im,jm,km,ptop)
    implicit none
    integer im,jm,km
    real    t(im,jm,km),h(im,jm,km)
    real    pstar(im,jm),ht(im,jm)
    real    sig(km)
    real    ptop
    integer i,j,k
    real    tbar
!
    do j = 1 , jm
      do i = 1 , im
         h(i,j,km) = ht(i,j) + rgas*regrav*t(i,j,km) &
                   * log(pstar(i,j)/((pstar(i,j)-ptop)*sig(km)+ptop))
      end do
    end do
    do k = km - 1, 1, -1
      do j = 1 , jm
        do i = 1 , im
          tbar = 0.5*( t(i,j,k)+t(i,j,k+1) )
          h(i,j,k) = h(i,j,k+1) +rgas*regrav*tbar  &
                   * log(((pstar(i,j)-ptop)*sig(k+1)+ptop)  &
                        /((pstar(i,j)-ptop)*sig(k)+ptop))
        end do
      end do
    end do
  end subroutine htsig

end module mod_derived

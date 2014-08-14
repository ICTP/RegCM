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

module mod_humid

  implicit none

  private

  ! numbers
  real(8) , parameter :: d_zero = 0.0D+00
  real(8) , parameter :: d_one = 1.0D+00
  real(8) , parameter :: d_two = 2.0D+00
  real(8) , parameter :: d_half = 0.50D+00
  real(8) , parameter :: d_60 = 60.0D+00
  real(8) , parameter :: d_10 = 1.0D+01
  real(8) , parameter :: d_r10 = 1.0D-01
  real(8) , parameter :: d_100 = 1.0D+02
  real(8) , parameter :: d_r100 = 1.0D-02
  real(8) , parameter :: d_1000 = 1.0D+03
  real(8) , parameter :: d_r1000 = 1.0D-03

  ! Low/Hi values
  real(8) , parameter :: minqx   = 1.0D-14
  real(8) , parameter :: dlowval = 1.0D-30
  real(8) , parameter :: dhival  = 1.0D+30
  real(4) , parameter :: slowval = 1.0E-30
  real(4) , parameter :: shival  = 1.0E+30
  real(8) , parameter :: dmissval = 1.0D+20
  real(4) , parameter :: smissval = 1.0E+20

  ! Standard Gravity (m/sec**2) 3rd CGPM
  real(8) , parameter :: egrav = 9.80665D+00

  ! Boltzman Constant k CODATA 2007
  real(8) , parameter :: boltzk = 1.3806504D-23
  ! Avogadro Constant
  real(8) , parameter :: navgdr = 6.02214129D23
  ! Effective molecular weight of dry air (g/mol)
  real(8) , parameter :: amd = 28.9644D+00
  ! Effective molecular weight of water (g/mol)
  real(8) , parameter :: amw = 18.0153D+00
  ! Effective molecular weight of ozone (g/mol)
  real(8) , parameter :: amo = 47.9942D+00
  ! Effective molecular weight of carbon dioxide (g/mol)
  real(8) , parameter :: amco2 = 44.01D+00

  real(8) , parameter :: rgasmol = navgdr*boltzk
  ! Gas constant for dry air in Joules/kg/K
  real(8) , parameter :: rgas = (rgasmol/amd)*1000.0D+00
  real(8) , parameter :: rdry = rgas

  ! Specific heat at constant pressure for dry air J/kg/K
  real(8) , parameter :: cpd = 3.5D+00*rgas

  ! Various utility terms used in calculations
  real(8) , parameter :: tzero = 273.15D+00

  ! Ratio of mean molecular weight of water to that of dry air
  real(8) , parameter :: ep2 = amw/amd

  public :: humid1 , humid2
  public :: pfesat , pfesdt , pfqsat , pfqsdt , sig2p , msig2p , wlh

  contains

  subroutine humid1(im,jm,km,t,q,ps,sigma,ptop,p,rh)
    implicit none
    integer , intent(in) :: im , jm , km
    real(8) , intent(in) :: ptop
    real(8) , intent(in) :: p
    real(4) , intent(in) , dimension(jm,im) :: ps
    real(4) , intent(in) , dimension(km,jm,im) :: t , q
    real(4) , intent(in) , dimension(km) :: sigma
    real(4) , intent(out) , dimension(jm,im) :: rh

    real(4) :: qs , qa , sigp , rpt , wp , w1
    real(8) :: pt
    integer :: i , j , k , kx , knx
    real(4) , dimension(km) :: rhs
    !
    ! THIS ROUTINE COMPUTES SPECIFIC HUMIDITY 
    ! DATA ON PRESSURE LEVEL P
    !
    rpt = real(ptop)
    do j = 1 , jm
      do i = 1 , im
        do k = 1 , km
          pt = dble(sigma(k))*(dble(ps(j,i))-ptop) + ptop
          qs = real(pfqsat(dble(t(k,j,i)),pt*d_100))
          qa = q(k,j,i)/(1.0-q(k,j,i))
          rhs(k) = max(qa/qs,0.0)*100.0
        end do
        ! Now vertical interpolation on pressure level!
        sigp = (real(p)-rpt)/(ps(j,i)-rpt)
        !
        ! Over the top or below bottom level
        !
        if ( sigp <= sigma(km) ) then
          rh(j,i) = rhs(km)
        else if ( sigp >= sigma(1) ) then
          rh(j,i) = rhs(1)
        else
          !
          ! Search k level below the requested one
          !
          kx = 0
          do k = 1 , km-1
            if ( sigp > sigma(k) ) exit
            kx = k
          end do
          !
          ! This is the above level
          !
          knx = kx + 1
          wp = (sigp-sigma(kx))/(sigma(knx)-sigma(kx))
          w1 = 1.0 - wp
          rh(j,i) = w1*rhs(kx) + wp*rhs(knx)
        end if
      end do
    end do
  end subroutine humid1

  subroutine humid2(im,jm,tas,qas,ps,rh)
    implicit none
    integer , intent(in) :: im , jm
    real(4) , dimension(jm,im) , intent(in) :: tas , qas , ps
    real(4) , dimension(jm,im) , intent(out) :: rh
    real(4) :: qs , qa
    integer :: i , j
    do j = 1 , jm
      do i = 1 , im
        qs = real(pfqsat(dble(tas(j,i)),dble(ps(j,i)*100.0)))
        qa = qas(j,i)/(1.0-qas(j,i))
        rh(j,i) = (qa/qs)*100.0
      end do
    end do
  end subroutine humid2

  ! Computes saturation pressurre
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
  real(8) function pfesat(t)
    implicit none
    real(8) , intent(in) :: t     ! Temperature (K)

    real(8) :: td , t_limit
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(8) , parameter :: a0 =  6.11213476D0
    real(8) , parameter :: a1 =  0.444007856D0
    real(8) , parameter :: a2 =  0.143064234D-01
    real(8) , parameter :: a3 =  0.264461437D-03
    real(8) , parameter :: a4 =  0.305903558D-05
    real(8) , parameter :: a5 =  0.196237241D-07
    real(8) , parameter :: a6 =  0.892344772D-10
    real(8) , parameter :: a7 = -0.373208410D-12
    real(8) , parameter :: a8 =  0.209339997D-15
    !
    ! For ice (temperature range -75C-0C)
    !
    real(8) , parameter :: c0 =  6.11123516D0
    real(8) , parameter :: c1 =  0.503109514D0
    real(8) , parameter :: c2 =  0.188369801D-01
    real(8) , parameter :: c3 =  0.420547422D-03
    real(8) , parameter :: c4 =  0.614396778D-05
    real(8) , parameter :: c5 =  0.602780717D-07
    real(8) , parameter :: c6 =  0.387940929D-09
    real(8) , parameter :: c7 =  0.149436277D-11
    real(8) , parameter :: c8 =  0.262655803D-14

    t_limit = t - tzero
    if ( t_limit > 100.0D0 ) t_limit = 100.0D0
    if ( t_limit < -75.0D0 ) t_limit = -75.0D0
    td = t_limit
    if ( td >= 0.0D0 ) then
      pfesat = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
         + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    else
      pfesat = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
         + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    end if
    pfesat = pfesat * d_100 ! pa
  end function pfesat
  !
  ! Computes derivative in temperature of saturation pressure
  !
  real(8) function pfesdt(t)
    implicit none
    real(8), intent(in)  :: t     ! Temperature (K)

    real(8) :: td , t_limit
    !
    ! For derivative:water vapor
    !
    real(8), parameter :: b0 =  0.444017302D0
    real(8), parameter :: b1 =  0.286064092D-01
    real(8), parameter :: b2 =  0.794683137D-03
    real(8), parameter :: b3 =  0.121211669D-04
    real(8), parameter :: b4 =  0.103354611D-06
    real(8), parameter :: b5 =  0.404125005D-09
    real(8), parameter :: b6 = -0.788037859D-12
    real(8), parameter :: b7 = -0.114596802D-13
    real(8), parameter :: b8 =  0.381294516D-16
    !
    ! For derivative:ice
    !
    real(8), parameter :: d0 =  0.503277922D0
    real(8), parameter :: d1 =  0.377289173D-01
    real(8), parameter :: d2 =  0.126801703D-02
    real(8), parameter :: d3 =  0.249468427D-04
    real(8), parameter :: d4 =  0.313703411D-06
    real(8), parameter :: d5 =  0.257180651D-08
    real(8), parameter :: d6 =  0.133268878D-10
    real(8), parameter :: d7 =  0.394116744D-13
    real(8), parameter :: d8 =  0.498070196D-16

    t_limit = t - tzero
    if ( t_limit > 100.0D0 ) t_limit = 100.0D0
    if ( t_limit < -75.0D0 ) t_limit = -75.0D0
    td = t_limit
    if ( td >= 0.0D0 ) then
      pfesdt = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
           + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
      pfesdt = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
           + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    end if
    pfesdt = pfesdt * 100.D0 ! pa/K
  end function pfesdt

  real(8) function pfqsat(t,p,e)
    implicit none
    real(8) , intent(in) :: t             ! Temperature (K)
    real(8) , intent(in) :: p             ! Pressure (Pa)
    real(8) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(8) :: es , vp , vp1 , vp2
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    ! Bolton 1980
    vp  = 1.0D0 / (p - 0.378D0*es)
    vp1 = ep2 * vp
    vp2 = vp1 * vp
    pfqsat = max(es * vp1, minqx)  ! kg/kg
  end function pfqsat

  real(8) function pfqsdt(t,p,e,dedt)
    implicit none
    real(8) , intent(in) :: t             ! Temperature (K)
    real(8) , intent(in) :: p             ! Pressure (Pa)
    real(8) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(8) , intent(in) , optional :: dedt ! derivative of e in dt (Pa/K)
    real(8) :: es , esdt , vp , vp1 , vp2
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    if ( present(dedt) ) then
      esdt = dedt
    else
      esdt = pfesdt(t)
    end if
    vp  = 1.0D0 / (p - 0.378D0*es)
    vp1 = ep2 * vp
    vp2 = vp1 * vp
    pfqsdt = esdt * vp2 * p ! 1 / K
  end function pfqsdt

  real(8) function sig2p(ps,sigma,ptop)
    implicit none
    real(8) , intent(in) :: ps , sigma , ptop
    !!!!!!!!!! Assume input is in cbar !!!!!!!!!!
    sig2p = (sigma * ( ps - ptop ) + ptop) * d_1000 ! Pressure in Pa
    !!!!!!!!!! Assume input is in cbar !!!!!!!!!!
  end function sig2p

  real(8) function msig2p(ps,sigma,ptop)
    implicit none
    real(8) , intent(in) :: ps , sigma , ptop
    !!!!!!!!!!!!! Assume input is in cbar !!!!!!!!!!!!!
    msig2p = (sigma * ps + ptop) * d_1000 ! Pressure in Pa !
    !!!!!!!!! In the model, ps is ps-ptop !!!!!!!!!!!!!
    !!!!!!!!!!!!! Assume input is in cbar !!!!!!!!!!!!!
  end function msig2p

  real(8) function wlh(t)
    implicit none
    real(8) , intent(in) :: t
    if ( t > tzero ) then
      wlh = 2500.8D0 - 2.36D0*t + 0.0016D0*t*t - 0.00006D0*t*t*t
    else
      wlh = 2834.1D0 - 0.29D0*t - 0.004D0*t*t
    end if
  end function wlh

end module mod_humid

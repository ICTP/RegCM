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

module mod_capecin

  use mod_realkinds
  use mod_intkinds
  use mod_constants

  implicit none

  real(rkx) , parameter :: pinc = 100.0_rkx ! Pressure increment (Pa)
                                            ! smaller number yields more
                                            ! accurate results, larger
                                            ! number makes code go faster

  integer(ik4) , parameter :: source = 2 ! Source parcel:
                                         ! 1 = surface
                                         ! 2 = most unstable (max theta-e)
                                         ! 3 = mixed-layer (specify ml_depth)

  real(rkx) , parameter :: ml_depth =  200.0_rkx ! depth (m) of mixed layer
                                                 ! for source=3

  integer(ik4) , parameter :: adiabat = 3 ! Formulation of moist adiabat:
                                          ! 1 = pseudoadiabatic, liquid only
                                          ! 2 = reversible, liquid only
                                          ! 3 = pseudoadiabatic, with ice
                                          ! 4 = reversible, with ice

  public :: getcape

  contains
  !
  !
  !  getcape - a fortran90 subroutine to calculate Convective Available
  !            Potential Energy (CAPE) from a sounding.
  !
  !  Version 1.02                           Last modified:  10 October 2008
  !
  !  Author:  George H. Bryan
  !           Mesoscale and Microscale Meteorology Division
  !           National Center for Atmospheric Research
  !           Boulder, Colorado, USA
  !           gbryan@ucar.edu
  !
  !  Disclaimer:  This code is made available WITHOUT WARRANTY.
  !
  !  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
  !               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
  !
  !  Input:     nk - number of levels in the sounding (integer)
  !
  !              p - one-dimensional array of pressure (Pa) (real)
  !
  !              t - one-dimensional array of temperature (K) (real)
  !
  !             rh - one-dimensional array of relative humidity (0-1)
  !
  !  Output:  cape - Convective Available Potential Energy (J/kg) (real)
  !
  !            cin - Convective Inhibition (J/kg) (real)
  !
  subroutine getcape(nk,p,t,rh,cape,cin)
    implicit none

    integer(ik4) , intent(in) :: nk
    real(rkx) , dimension(nk) , intent(in) :: p , t , rh
    real(rkx) , intent(out) :: cape , cin

    logical :: doit , ice , cloud , not_converged
    integer(ik4) :: k , kmax , n , nloop , i
    real(rkx) , dimension(nk) :: td , pi , q , th , thv , z

    real(rkx) :: the , maxthe , parea , narea , lfc
    real(rkx) :: th1 , p1 , t1 , qv1 , ql1 , qi1 , b1 , pi1
    real(rkx) :: thv1 , qt , dp , dz , ps , frac
    real(rkx) :: th2 , p2 , t2 , qv2 , ql2 , qi2 , b2 , pi2 , thv2
    real(rkx) :: thlast , fliq , fice , tbar , qvbar , qlbar , qibar
    real(rkx) :: lhv , lhs , lhf , rm , cpm
    real(rkx) :: avgth , avgqv

    real(rkx), parameter :: lv1   = wlhv+(cpw-cpv)*tzero
    real(rkx), parameter :: lv2   = cpw-cpv
    real(rkx), parameter :: ls1   = wlhs+(cpi-cpv)*tzero
    real(rkx), parameter :: ls2   = cpi-cpv
    real(rkx), parameter :: rp00  = d_one/p00
    real(rkx), parameter :: reps  = rwat/rgas
    real(rkx), parameter :: rddcp = rgas/cpd
    real(rkx), parameter :: cpdg  = cpd*regrav

    real(rkx), parameter :: converge = 0.002_rkx

    ! Get td,pi,q,th,thv

    do k = 1 , nk
      pi(k) = (p(k)*rp00)**rddcp
      td(k) = getdewp(t(k)-tzero,rh(k))
      q(k) = getqvs(p(k),td(k))
      th(k) = t(k)/pi(k)
      thv(k) = th(k)*(d_one+reps*q(k))/(d_one+q(k))
    end do

    ! get height using the hydrostatic equation

    z(1) = d_zero
    do k = 2 , nk
      dz = -cpdg*0.5_rkx*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))
      z(k) = z(k-1) + dz
    end do

    ! Find source parcel

    kmax = 1

    select case (source)
      case (1)
        ! use surface parcel
        k    = kmax
        th2  = th(kmax)
        pi2  = pi(kmax)
        p2   = p(kmax)
        t2   = t(kmax)
        thv2 = thv(kmax)
        qv2  = q(kmax)
        b2   = d_zero
      case (2)
        ! use most unstable parcel (max theta-e)
        if ( p(1) < 50000.0_rkx ) then
          ! first report is above 500 mb ... just use the first level reported
          maxthe = getthe(p(1),t(1),td(1),q(1))
        else
          ! find max thetae below 500 mb
          maxthe = d_zero
          do k = 1 , nk
            if ( p(k) >= 50000.0_rkx ) then
              the = getthe(p(k),t(k),td(k),q(k))
              if ( the > maxthe ) then
                maxthe = the
                kmax = k
              end if
            end if
          end do
        end if
        k    = kmax
        th2  = th(kmax)
        pi2  = pi(kmax)
        p2   = p(kmax)
        t2   = t(kmax)
        thv2 = thv(kmax)
        qv2  = q(kmax)
        b2   = d_zero
      case default
        ! use mixed layer
        if ( (z(2)-z(1)) > ml_depth ) then
          ! the second level is above the mixed-layer depth:  just use the
          ! lowest level
          avgth = th(1)
          avgqv = q(1)
        else if ( z(nk) < ml_depth ) then
          ! the top-most level is within the mixed layer:  just use the
          ! upper-most level
          avgth = th(nk)
          avgqv = q(nk)
          kmax = nk
        else
          ! calculate the mixed-layer properties:
          avgth = d_zero
          avgqv = d_zero
          k = 2
          do while ( (z(k) <= ml_depth) .and. (k <= nk) )
            avgth = avgth + 0.5_rkx*(z(k)-z(k-1))*(th(k)+th(k-1))
            avgqv = avgqv + 0.5_rkx*(z(k)-z(k-1))*(q(k)+q(k-1))
            k = k + 1
          end do
          th2 = th(k-1)+(th(k)-th(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
          qv2 =  q(k-1)+( q(k)- q(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
          avgth = avgth + 0.5_rkx*(ml_depth-z(k-1))*(th2+th(k-1))
          avgqv = avgqv + 0.5_rkx*(ml_depth-z(k-1))*(qv2+q(k-1))
          avgth = avgth/ml_depth
          avgqv = avgqv/ml_depth
        end if
        k    = kmax
        th2  = avgth
        qv2  = avgqv
        thv2 = th2*(d_one+reps*qv2)/(d_one+qv2)
        pi2  = pi(kmax)
        p2   = p(kmax)
        t2   = th2*pi2
        b2   = egrav*( thv2-thv(kmax) )/thv(kmax)
    end select

    ! Define parcel properties at initial location
    narea = d_zero
    ql2 = d_zero
    qi2 = d_zero
    qt  = qv2

    cape = d_zero
    cin  = d_zero
    lfc  = d_zero

    doit = .true.
    cloud = .false.

    if ( adiabat == 1 .or. adiabat == 2 ) then
      ice = .false.
    else
      ice = .true.
    end if

    the = getthe(p2,t2,t2,qv2)

    ! Begin ascent of parcel

    do while ( doit .and. (k < nk) )
      k = k+1
      b1 =  b2
      dp = p(k-1)-p(k)
      if ( dp < pinc ) then
        nloop = 1
      else
        nloop = 1 + int( dp/pinc )
        dp = dp/real(nloop,rkx)
      end if
      do n = 1 , nloop
        p1 =  p2
        t1 =  t2
        pi1 = pi2
        th1 = th2
        qv1 = qv2
        ql1 = ql2
        qi1 = qi2
        thv1 = thv2
        p2 = p2 - dp
        pi2 = (p2*rp00)**rddcp
        thlast = th1
        i = 0
        not_converged = .true.
        do while ( not_converged .and. i < 100 )
          i = i + 1
          t2 = thlast*pi2
          if ( ice ) then
            fliq = max(min((t2-233.15_rkx)/(tzero-233.15),d_one),d_zero)
            fice = d_one-fliq
          else
            fliq = d_one
            fice = d_zero
          endif
          qv2 = min(qt, fliq*getqvs(p2,t2) + fice*getqvi(p2,t2) )
          qi2 = max(fice*(qt-qv2), d_zero)
          ql2 = max(qt-qv2-qi2, d_zero)
          tbar  = d_half*(t1+t2)
          qvbar = d_half*(qv1+qv2)
          qlbar = d_half*(ql1+ql2)
          qibar = d_half*(qi1+qi2)

          lhv = lv1-lv2*tbar
          lhs = ls1-ls2*tbar
          lhf = lhs-lhv

          rm = rgas+rwat*qvbar
          cpm = cpd+cpv*qvbar+cpw*qlbar+cpi*qibar
          th2 = th1*exp( lhv*(ql2-ql1)/(cpm*tbar) +   &
                         lhs*(qi2-qi1)/(cpm*tbar) +   &
                         (rm/cpm-rgas/cpd)*log(p2/p1) )
          if ( abs(th2-thlast) < converge ) then
            thlast = thlast+0.3_rkx*(th2-thlast)
          else
            not_converged = .false.
          end if
        end do

        ! Latest pressure increment is complete.  Calculate some
        ! important stuff:

        if ( ql2 > 1.0e-10_rkx ) cloud = .true.

        if ( adiabat == 1 .or. adiabat == 3 ) then
          ! pseudoadiabat
          qt  = qv2
          ql2 = d_zero
          qi2 = d_zero
        end if
      end do

      thv2 = th2*(d_one+reps*qv2)/(d_one+qv2+ql2+qi2)
      b2 = egrav*(thv2-thv(k))/thv(k)
      dz = -cpdg*d_half*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))

      the = getthe(p2,t2,t2,qv2)

      ! Get contributions to CAPE and CIN:

      if ( (b2 >= d_zero) .and. (b1 <= d_zero) ) then
        ! first trip into positive area
        ps = p(k-1)+(p(k)-p(k-1))*(d_zero-b1)/(b2-b1)
        frac = b2/(b2-b1)
        parea = d_half*b2*dz*frac
        narea = narea-d_half*b1*dz*(d_one-frac)
        cin = cin + narea
        narea = d_zero
      else if ( (b2 < d_zero) .and. (b1 > d_zero) ) then
        ! first trip into neg area
        ps = p(k-1)+(p(k)-p(k-1))*(d_zero-b1)/(b2-b1)
        frac = b1/(b1-b2)
        parea = d_half*b1*dz*frac
        narea = -d_half*b2*dz*(d_one-frac)
      else if ( b2 < d_zero ) then
        ! still collecting negative buoyancy
        parea = d_zero
        narea = narea-d_half*dz*(b1+b2)
      else
        ! still collecting positive buoyancy
        parea = d_half*dz*(b1+b2)
        narea = d_zero
      endif

      cape = cape + max(d_zero,parea)

      if ( (p(k) <= 10000.0_rkx) .and. (b2 < d_zero) ) then
        ! stop if b < 0 and p < 100 mb
        doit = .false.
      end if

    end do

    contains

    pure real(rkx) function getdewp(tc,rh)
      implicit none
      real(rkx) , intent(in) :: tc , rh
      real(rkx) , parameter :: b = 18.678_rkx
      real(rkx) , parameter :: c = 257.14_rkx ! [C]
      real(rkx) , parameter :: d = 234.50_rkx ! [C]
      real(rkx) :: gm
      gm = log(rh * exp((b - tc/d)*(tc/(c+tc))))
      getdewp = tzero + c * gm/(b-gm)
    end function getdewp

    pure real(rkx) function getqvs(p,t)
      implicit none
      real(rkx) , intent(in) :: p , t
      real(rkx) :: es
      es = 611.2_rkx*exp(17.67_rkx*(t-273.15_rkx)/(t-29.65_rkx))
      getqvs = ep2*es/(p-es)
    end function getqvs

    pure real(rkx) function getqvi(p,t)
      implicit none
      real(rkx) , intent(in) :: p , t
      real(rkx) :: es
      es = 611.2_rkx*exp(21.8745584_rkx*(t-273.15_rkx)/(t-7.66_rkx))
      getqvi = ep2*es/(p-es)
    end function getqvi

    pure real(rkx) function getthe(p,t,td,q)
      implicit none
      real(rkx) , intent(in) :: p , t , td , q
      real(rkx) :: tlcl
      if ( (td-t) >= -0.1_rkx ) then
        tlcl = t
      else
        tlcl = 56.0_rkx + ( (td-56.0_rkx)**(-1.0_rkx) + &
               0.00125_rkx * log(t/td) )**(-1.0_rkx)
      end if
      getthe = t * ( (100000.0_rkx/p)**(0.2854_rkx*(1.0_rkx-0.28_rkx*q)) ) * &
               exp( ((3376.0_rkx/tlcl)-2.54_rkx)*q*(1.0_rkx+0.81_rkx*q) )
    end function getthe

    end subroutine getcape

end module mod_capecin
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

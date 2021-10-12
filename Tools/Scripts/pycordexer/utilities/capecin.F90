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

  implicit none

  real , parameter :: rovg = 29.2716599

  real , parameter :: pinc = 100.0 ! Pressure increment (Pa)
                                   ! smaller number yields more
                                   ! accurate results, larger
                                   ! number makes code go faster

  integer , parameter :: source = 2 ! Source parcel:
                                    ! 1 = surface
                                    ! 2 = most unstable (max theta-e)
                                    ! 3 = mixed-layer (specify ml_depth)

  real , parameter :: ml_depth =  200.0 ! depth (m) of mixed layer
                                        ! for source=3

  integer , parameter :: adiabat = 1 ! Formulation of moist adiabat:
                                     ! 1 = pseudoadiabatic, liquid only
                                     ! 2 = reversible, liquid only
                                     ! 3 = pseudoadiabatic, with ice
                                     ! 4 = reversible, with ice

  public :: getcape
  public :: getcape_hy
  public :: getcape_nhy
  public :: getcape_moloch

  contains

  subroutine getcape(im,jm,km,p,t,rh,cape,cin)
    implicit none
    integer , intent(in) :: im , jm , km
    real(4) , intent(in) , dimension(km,jm,im) :: p , t , rh
    real(4) , intent(out) , dimension(jm,im) :: cape , cin
    real(4) , dimension(km) :: pa , ta , rha
    integer :: i , j , k , kk , iloop , icount
    real(4) :: z1
    icount = im * jm
!$OMP PARALLEL DO PRIVATE (i,j,k,kk,pa,ta,rha,z1)
    do iloop = 1 , icount
      i = iloop/jm + 1
      j = iloop - (i-1)*jm
      do k = 1 , km
        kk = km - k + 1
        pa(kk) = p(k,j,i)
        ta(kk) = t(k,j,i)
        rha(kk) = rh(k,j,i) * 0.01
      end do
      z1 = 40.0 ! Arbitrary here: Assume lower model level to be at z1 in meters
      call capecin(km,z1,pa,ta,rha,cape(j,i),cin(j,i))
    end do
!$OMP END PARALLEL DO
  end subroutine getcape

  subroutine getcape_moloch(im,jm,km,ps,t,rh,pai,cape,cin)
    implicit none
    integer , intent(in) :: im , jm , km
    real(4) , intent(in) , dimension(jm,im) :: ps
    real(4) , intent(in) , dimension(km,jm,im) :: pai , t , rh
    real(4) , intent(out) , dimension(jm,im) :: cape , cin
    real(4) , dimension(km) :: pa , ta , rha
    integer :: i , j , k , kk , iloop , icount
    real(4) :: z1
    icount = im * jm
!$OMP PARALLEL DO PRIVATE (i,j,k,kk,pa,ta,rha,z1)
    do iloop = 1 , icount
      i = iloop/jm + 1
      j = iloop - (i-1)*jm
      do k = 1 , km
        kk = km - k + 1
        pa(kk) = 100000.0*pai(k,j,i)**3.5
        ta(kk) = t(k,j,i)
        rha(kk) = rh(k,j,i) * 0.01
      end do
      z1 = rovg * ta(1) * log(ps(j,i)/pa(1))
      call capecin(km,z1,pa,ta,rha,cape(j,i),cin(j,i))
    end do
!$OMP END PARALLEL DO
  end subroutine getcape_moloch

  subroutine getcape_hy(im,jm,km,ps,t,rh,sigma,ptop,cape,cin)
    implicit none
    integer , intent(in) :: im , jm , km
    real(4) , intent(in) , dimension(jm,im) :: ps
    real(4) , intent(in) , dimension(km,jm,im) :: t , rh
    real(8) , intent(in) :: ptop
    real(4) , intent(in) , dimension(km) :: sigma
    real(4) , intent(out) , dimension(jm,im) :: cape , cin
    real(4) , dimension(km) :: pa , ta , rha
    integer :: i , j , k , kk , iloop , icount
    real(4) :: ptp , z1
    ptp = real(ptop) * 100.0
    icount = im * jm
!$OMP PARALLEL DO PRIVATE (i,j,k,kk,pa,ta,rha,z1)
    do iloop = 1 , icount
      i = iloop/jm + 1
      j = iloop - (i-1)*jm
      do k = 1 , km
        kk = km - k + 1
        pa(kk) = sigma(k)*(ps(j,i)-ptp) + ptp
        ta(kk) = t(k,j,i)
        rha(kk) = rh(k,j,i) * 0.01
      end do
      z1 = rovg * ta(1) * log(ps(j,i)/pa(1))
      call capecin(km,z1,pa,ta,rha,cape(j,i),cin(j,i))
    end do
!$OMP END PARALLEL DO
  end subroutine getcape_hy

  subroutine getcape_nhy(im,jm,km,ps,t,p0,rh,sigma,ptop,pp,cape,cin)
    implicit none
    integer , intent(in) :: im , jm , km
    real(4) , intent(in) , dimension(jm,im) :: ps , p0
    real(4) , intent(in) , dimension(km,jm,im) :: t , rh , pp
    real(8) , intent(in) :: ptop
    real(4) , intent(in) , dimension(km) :: sigma
    real(4) , intent(out) , dimension(jm,im) :: cape , cin
    real(4) , dimension(km) :: pa , ta , rha
    integer :: i , j , k , kk , iloop , icount
    real(4) :: ptp , z1
    real(4) , dimension(jm,im) :: pstar
    ptp = real(ptop) * 100.0
    pstar = p0 - ptp
    icount = im * jm
!$OMP PARALLEL DO PRIVATE (i,j,k,kk,pa,ta,rha,z1)
    do iloop = 1 , icount
      i = iloop/jm + 1
      j = iloop - (i-1)*jm
      do k = 1 , km
        kk = km - k + 1
        pa(kk) = sigma(k)*pstar(j,i) + ptp + pp(k,j,i)
        ta(kk) = t(k,j,i)
        rha(kk) = rh(k,j,i) * 0.01
      end do
      z1 = rovg * ta(1) * log(ps(j,i)/pa(1))
      call capecin(km,z1,pa,ta,rha,cape(j,i),cin(j,i))
    end do
!$OMP END PARALLEL DO
  end subroutine getcape_nhy
  !
  !
  !  capecin - a fortran90 subroutine to calculate Convective Available
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
  subroutine capecin(nk,z1,p,t,rh,cape,cin)
    implicit none

    integer , intent(in) :: nk
    real , dimension(nk) , intent(in) :: p , t , rh
    real , intent(in) :: z1
    real , intent(out) :: cape , cin

    logical :: doit , ice , cloud , not_converged
    integer :: k , kmax , n , nloop , i
    real , dimension(nk) :: td , pi , q , th , thv , z

    real :: the , maxthe , parea , narea , lfc
    real :: th1 , p1 , t1 , qv1 , ql1 , qi1 , b1 , pi1
    real :: thv1 , qt , dp , dz , ps , frac
    real :: th2 , p2 , t2 , qv2 , ql2 , qi2 , b2 , pi2 , thv2
    real :: thlast , fliq , fice , tbar , qvbar , qlbar , qibar
    real :: lhv , lhs , lhf , rm , cpm
    real :: avgth , avgqv

    real , parameter :: wlhv = 2.50080e6
    real , parameter :: wlhf = 0.33355e6
    real , parameter :: wlhs = 2.83435e6
    real , parameter :: cpv = 1846.0932676
    real , parameter :: cpw = 4186.95
    real , parameter :: cpi = 2117.27
    real , parameter :: rgas = 287.0569248
    real , parameter :: egrav = 9.80665
    real , parameter :: regrav = 1.0/egrav
    real , parameter :: tzero = 273.15
    real , parameter :: cpd = 1004.6992368
    real , parameter :: rwat = 461.5233169
    real , parameter :: p00 = 1.000000e5
    real , parameter :: ep2 = 0.6219770795
    real , parameter :: lv1   = wlhv+(cpw-cpv)*tzero
    real , parameter :: lv2   = cpw-cpv
    real , parameter :: ls1   = wlhs+(cpi-cpv)*tzero
    real , parameter :: ls2   = cpi-cpv
    real , parameter :: rp00  = 1.0/p00
    real , parameter :: reps  = 1.0/ep2
    real , parameter :: rddcp = rgas/cpd
    real , parameter :: cpdg  = cpd*regrav

    real, parameter :: converge = 0.0002

    ! Get td,pi,q,th,thv

    do k = 1 , nk
      pi(k) = (p(k)*rp00)**rddcp
      td(k) = getdewp(t(k),rh(k))
      q(k) = getqvs(p(k),td(k))
      th(k) = t(k)/pi(k)
      thv(k) = th(k)*(1.0+reps*q(k))/(1.0+q(k))
    end do

    ! get height using the hydrostatic equation

    z(1) = z1
    do k = 2 , nk
      dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))
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
        b2   = 0.0
      case (2)
        ! use most unstable parcel (max theta-e)
        if ( p(1) < 50000.0 ) then
          ! first report is above 500 mb ... just use the first level reported
          maxthe = getthe(p(1),t(1),td(1),q(1))
        else
          ! find max thetae below 500 mb
          maxthe = 0.0
          do k = 1 , nk
            if ( p(k) >= 50000.0 ) then
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
        b2   = 0.0
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
          avgth = 0.0
          avgqv = 0.0
          k = 2
          do while ( (z(k) <= ml_depth) .and. (k <= nk) )
            avgth = avgth + 0.5*(z(k)-z(k-1))*(th(k)+th(k-1))
            avgqv = avgqv + 0.5*(z(k)-z(k-1))*(q(k)+q(k-1))
            k = k + 1
          end do
          th2 = th(k-1)+(th(k)-th(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
          qv2 =  q(k-1)+( q(k)- q(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
          avgth = avgth + 0.5*(ml_depth-z(k-1))*(th2+th(k-1))
          avgqv = avgqv + 0.5*(ml_depth-z(k-1))*(qv2+q(k-1))
          avgth = avgth/ml_depth
          avgqv = avgqv/ml_depth
        end if
        k    = kmax
        th2  = avgth
        qv2  = avgqv
        thv2 = th2*(1.0+reps*qv2)/(1.0+qv2)
        pi2  = pi(kmax)
        p2   = p(kmax)
        t2   = th2*pi2
        b2   = egrav*( thv2-thv(kmax) )/thv(kmax)
    end select

    ! Define parcel properties at initial location
    narea = 0.0
    ql2 = 0.0
    qi2 = 0.0
    qt  = qv2

    cape = 0.0
    cin  = 0.0
    lfc  = 0.0

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
        dp = dp/real(nloop)
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
            fliq = max(min((t2-233.15)/(tzero-233.15),1.0),0.0)
            fice = 1.0-fliq
          else
            fliq = 1.0
            fice = 0.0
          endif
          qv2 = min(qt, fliq*getqvs(p2,t2) + fice*getqvi(p2,t2) )
          qi2 = max(fice*(qt-qv2), 0.0)
          ql2 = max(qt-qv2-qi2, 0.0)
          tbar  = 0.5*(t1+t2)
          qvbar = 0.5*(qv1+qv2)
          qlbar = 0.5*(ql1+ql2)
          qibar = 0.5*(qi1+qi2)

          lhv = lv1-lv2*tbar
          lhs = ls1-ls2*tbar
          lhf = lhs-lhv

          rm = rgas+rwat*qvbar
          cpm = cpd+cpv*qvbar+cpw*qlbar+cpi*qibar
          th2 = th1*exp( lhv*(ql2-ql1)/(cpm*tbar) +   &
                         lhs*(qi2-qi1)/(cpm*tbar) +   &
                         (rm/cpm-rgas/cpd)*log(p2/p1) )
          if ( abs(th2-thlast) < converge ) then
            thlast = thlast+0.3*(th2-thlast)
          else
            not_converged = .false.
          end if
        end do

        ! Latest pressure increment is complete.  Calculate some
        ! important stuff:

        if ( ql2 > 1.0e-10 ) cloud = .true.

        if ( adiabat == 1 .or. adiabat == 3 ) then
          ! pseudoadiabat
          qt  = qv2
          ql2 = 0.0
          qi2 = 0.0
        end if
      end do

      thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2+qi2)
      b2 = egrav*(thv2-thv(k))/thv(k)
      dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))

      the = getthe(p2,t2,t2,qv2)

      ! Get contributions to CAPE and CIN:

      if ( (b2 >= 0.0) .and. (b1 <= 0.0) ) then
        ! first trip into positive area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b2/(b2-b1)
        parea = 0.5*b2*dz*frac
        narea = narea-0.5*b1*dz*(1.0-frac)
        cin = cin + narea
        narea = 0.0
      else if ( (b2 < 0.0) .and. (b1 > 0.0) ) then
        ! first trip into neg area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b1/(b1-b2)
        parea = 0.5*b1*dz*frac
        narea = -0.5*b2*dz*(1.0-frac)
      else if ( b2 < 0.0 ) then
        ! still collecting negative buoyancy
        parea = 0.0
        narea = narea-0.5*dz*(b1+b2)
      else
        ! still collecting positive buoyancy
        parea = 0.5*dz*(b1+b2)
        narea = 0.0
      endif

      cape = cape + max(0.0,parea)

      if ( (p(k) <= 10000.0) .and. (b2 < 0.0) ) then
        ! stop if b < 0 and p < 100 mb
        doit = .false.
      end if

    end do

    contains

    pure real function getdewp(t,rh)
      implicit none
      real , intent(in) :: t , rh
      real , parameter :: b = 18.678
      real , parameter :: c = 257.14 ! [C]
      real , parameter :: d = 234.50 ! [C]
      real :: tc , gm
      tc = t - tzero
      gm = log(rh * exp((b - tc/d)*(tc/(c+tc))))
      getdewp = tzero + c * gm/(b-gm)
    end function getdewp

    pure real function getqvs(p,t)
      implicit none
      real , intent(in) :: p , t
      real :: es
      es = 611.2*exp(17.67*(t-273.15)/(t-29.65))
      getqvs = ep2*es/(p-es)
    end function getqvs

    pure real function getqvi(p,t)
      implicit none
      real , intent(in) :: p , t
      real :: es
      es = 611.2*exp(21.8745584*(t-273.15)/(t-7.66))
      getqvi = ep2*es/(p-es)
    end function getqvi

    pure real function getthe(p,t,td,q)
      implicit none
      real , intent(in) :: p , t , td , q
      real :: tlcl
      if ( (td-t) >= -0.1 ) then
        tlcl = t
      else
        tlcl = 56.0 + ( (td-56.0)**(-1.0) + &
               0.00125 * log(t/td) )**(-1.0)
      end if
      getthe = t * ( (100000.0/p)**(0.2854*(1.0-0.28*q)) ) * &
               exp( ((3376.0/tlcl)-2.54)*q*(1.0+0.81*q) )
    end function getthe

  end subroutine capecin

end module mod_capecin
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

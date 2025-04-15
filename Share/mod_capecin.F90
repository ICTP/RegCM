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

module mod_capecin

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_spline

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

  logical :: table_empty = .true.

  integer(ik4) , parameter :: itb = 076
  integer(ik4) , parameter :: jtb = 134

  real(rkx) :: pl , thl , rdq , rdth , rdp , rdthe , plq , rdpq , rdtheq
  real(rkx) , dimension(jtb) :: qs0 , sqs
  real(rkx) , dimension(itb) :: the0 , sthe
  real(rkx) , dimension(itb,jtb) :: ptbl
  real(rkx) , dimension(jtb,itb) :: ttbl

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

    ! This routine computes a surface to 500mb lifted index.
    ! The lifted parcel is from the first atmpspheric ETA
    ! layer (ie, the ETA layer closest to the model ground).
    ! The lifted index is the difference between this parcel's
    ! temperature at 500mb and the ambient 500mb temperature.
    ! Russ Treadon W/NP2 @date 1993-03-10

    subroutine otlift(slindx,t,q,p,t500,ista,iend,jsta,jend,kk)
      implicit none
      integer(ik4) , intent(in) :: ista , iend , jsta , jend , kk
      real(rkx) , dimension(:,:) , pointer , intent(inout) :: slindx
      real(rkx) , dimension(:,:) , pointer , intent(in) :: t500
      real(rkx) , dimension(:,:,:) , pointer , intent(in) :: t , q , p

      real(rkx) , parameter :: d8202 = 0.820231e+00_rkx
      real(rkx) , parameter :: h5e4 = 5.e4_rkx
      real(rkx) , parameter :: p500 = 50000.0_rkx
      real(rkx) , parameter :: elivw = 2.72e6_rkx
      real(rkx) , parameter :: elocp = elivw/cpd
      real(rkx) , parameter :: oneps = 1.0_rkx-ep2
      real(rkx) , parameter :: pt = 1.0_rkx
      real(rkx) , parameter :: thl = 210.0_rkx

      real(rkx) :: tvp , esatp , qsatp
      real(rkx) :: tth , tp , apesp , partmp , thesp , tpsp
      real(rkx) :: bqs00 , sqs00 , bqs10 , sqs10 , bq , sq , tq
      real(rkx) :: pp00 , pp10 , pp01 , pp11 , t00 , t10 , t01 , t11
      real(rkx) :: bthe00 , sthe00 , bthe10 , sthe10 , bth , sth
      real(rkx) :: tqq , qq , qbt , tthbt , tbt , apebt , ppq , pp
      integer(ik4) :: i , j , ittbk , iq , it , iptbk
      integer(ik4) :: ith , ip , iqtb
      integer(ik4) :: ittb , iptb , ithtb
      !
      if ( table_empty ) then
        call table_fill( )
        table_empty = .false.
      end if
      !
      !**********************************************************
      ! Start otlift here
      !
      ! Initialize lifted index array to zero.
      do concurrent ( i = ista:iend, j = jsta:jend )
        slindx(i,j) = d_zero
      end do
      ! Find Exner at lowest level-------------------------------
#ifdef STDPAR
      do concurrent ( i = ista:iend, j = jsta:jend ) &
        local(tbt,qbt,apebt,tthbt,tth,tqq,ittb,ittbk,bqs00,sqs00,&
        bqs10,sqs10,bq,sq,tq,ppq,iqtb,iq,it,pp00,pp10,pp01,pp11, &
        tpsp,apesp,thesp,tp,qq,iptb,iptbk,bthe00,sthe00,bthe10,  &
        sthe10,bth,sth,pp,ithtb,ith,ip,t00,t10,t01,t11,partmp,   &
        esatp,qsatp,tvp)
#else
      do j = jsta , jend
        do i = ista , iend
#endif
          tbt = t(i,j,kk)
          ! Specific Humidity expected.
          qbt = q(i,j,kk)/(1.0_rkx+q(i,j,kk))
          apebt = (p00/p(i,j,kk))**rovcp
          ! Scaling potential temperature & table index----------
          tthbt = tbt*apebt
          tth = (tthbt-thl)*rdth
          tqq = tth-aint(tth)
          ittb = int(tth) + 1
          ! Keeping indices within the table---------------------
          if ( ittb < 1 ) then
            ittb = 1
            tqq = d_zero
          end if
          if ( ittb >= jtb ) then
            ittb = jtb - 1
            tqq = d_zero
          end if
          ! Base and scaling factor for spec. humidity-----------
          ittbk = ittb
          bqs00 = qs0(ittbk)
          sqs00 = sqs(ittbk)
          bqs10 = qs0(ittbk+1)
          sqs10 = sqs(ittbk+1)
          ! Scaling spec. humidity & table index-----------------
          bq = (bqs10-bqs00)*tqq + bqs00
          sq = (sqs10-sqs00)*tqq + sqs00
          tq = (qbt-bq)/sq*rdq
          ppq = tq - aint(tq)
          iqtb = int(tq) + 1
          ! Keeping indices within the table---------------------
          if ( iqtb < 1 ) then
            iqtb = 1
            ppq = d_zero
          end if
          if ( iqtb >= itb ) then
            iqtb = itb-1
            ppq = d_zero
          end if
          ! Saturation pressure at four surrounding table pts.---
          iq = iqtb
          it = ittb
          pp00 = ptbl(iq,it)
          pp10 = ptbl(iq+1,it)
          pp01 = ptbl(iq,it+1)
          pp11 = ptbl(iq+1,it+1)
          ! Saturation point variables at the bottom------------
          tpsp = pp00+(pp10-pp00)*ppq+(pp01-pp00)*tqq + &
                (pp00-pp10-pp01+pp11)*ppq*tqq
          if ( tpsp <= d_zero ) tpsp = p00
          apesp = (p00/tpsp)**rovcp
          thesp = tthbt*exp(elocp*qbt*apesp/tthbt)
          ! Scaling pressure & tt table index------------------
          tp = (h5e4-pl)*rdp
          qq = tp - aint(tp)
          iptb = int(tp)+1
          ! Keeping indices within the table-------------------
          if ( iptb < 1 ) then
            iptb = 1
            qq = d_zero
          end if
          if ( iptb >= itb ) then
            iptb = itb-1
            qq = d_zero
          end if
          ! Base and scaling factor for the-------------------
          iptbk = iptb
          bthe00 = the0(iptbk)
          sthe00 = sthe(iptbk)
          bthe10 = the0(iptbk+1)
          sthe10 = sthe(iptbk+1)
          ! Scaling the & tt table index----------------------
          bth = (bthe10-bthe00)*qq + bthe00
          sth = (sthe10-sthe00)*qq + sthe00
          tth = (thesp-bth)/sth*rdthe
          pp = tth-aint(tth)
          ithtb = int(tth) + 1
          ! Keeping indices within the table------------------
          if ( ithtb < 1 ) then
            ithtb = 1
            pp = d_zero
          end if
          if ( ithtb >= jtb ) then
            ithtb = jtb-1
            pp = d_zero
          end if
          ! Temperature at four surrounding tt table pts.----
          ith = ithtb
          ip = iptb
          t00 = ttbl(ith,ip)
          t10 = ttbl(ith+1,ip)
          t01 = ttbl(ith,ip+1)
          t11 = ttbl(ith+1,ip+1)
          ! Parcel temperature at 500mb----------------------
          if ( tpsp >= h5e4 ) then
            partmp=(t00+(t10-t00)*pp + (t01-t00)*qq + &
                   (t00-t10-t01+t11)*pp*qq)
          else
            partmp = tbt*apebt*d8202
          end if
          ! Lifted Index-------------------------------------
          !
          ! The parcel temperature at 500 mb has been
          ! computed, and we find the mixing ratio at that
          ! level which will be the saturation value since
          ! we're following a moist adiabat. Note that the
          ! ambient 500 mb should probably be virtualized,
          ! but the impact of moisture at that level is
          ! quite small
          !
          esatp = pfesat(partmp,p500)
          qsatp = ep2*esatp/(p500-esatp*oneps)
          tvp = partmp*(1.0_rkx+ep1*qsatp)
          slindx(i,j) = t500(i,j)-tvp
#ifndef STDPAR
        end do
#endif
      end do

      contains

#include <pfesat.inc>

      subroutine table_fill( )
        implicit none
        ! ****************************************************************
        ! *                                                              *
        ! *             GENERATE VALUES FOR LOOK-UP TABLES               *
        ! *                                                              *
        ! ****************************************************************
        real(rkx) , parameter :: thh = 365.0_rkx
        real(rkx) , parameter :: ph = 105000.0_rkx
        real(rkx) , parameter :: pq0 = 379.90516_rkx
        real(rkx) , parameter :: a2 = 17.2693882_rkx
        real(rkx) , parameter :: a3 = 273.16_rkx
        real(rkx) , parameter :: a4 = 35.86_rkx
        real(rkx) , parameter :: eliwv = 2.683e+6_rkx
        real(rkx) , parameter :: eps = 1.E-9_rkx
        real(rkx) , dimension(jtb) :: qsold , pold, qsnew , pnew , tnew
        real(rkx) , dimension(jtb) :: told , theold , thenew
        real(rkx) , dimension(jtb) :: app , apt , aqp , aqt , y2p , y2t
        real(rkx) :: dth , dp , th , p , ape , denom , qs0k , sqsk , dqs
        real(rkx) :: qs , sthek , the0k , dthe
        integer(ik4) :: kpm , kthm1 , kpm1 , kp , kthm , kth

        ! Coarse look-up table for saturation point----------------
        kthm  = jtb
        kpm   = itb
        kthm1 = kthm-1
        kpm1  = kpm-1
        pl = pt
        dth = (thh-thl) / real(kthm-1, rkx)
        dp  = (ph -pl ) / real(kpm -1, rkx)
        rdth = 1.0_rkx/dth
        rdp  = 1.0_rkx/dp
        rdq  = kpm-1
        th = thl - dth
        !-----------------------------------------------------------
        do kth = 1 , kthm
          th = th + dth
          p  = pl - dp
          do kp = 1 , kpm
            p = p + dp
            if ( p <= 0.0_rkx ) then
              pold(1)  = 0.0_rkx
              qsold(1) = 0.0_rkx
            else
              ape = (100000.0_rkx/p)**rovcp
              denom = th - a4*ape
              if ( denom > eps ) then
                qsold(kp) = pq0 / p*exp(a2*(th-a3*ape)/denom)
              else
                qsold(kp) = 0.0_rkx
              end if
              ! qsold(kp) = pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
              pold(kp) = p
            end if
          end do
          qs0k       = qsold(1)
          sqsk       = qsold(kpm) - qsold(1)
          qsold(1  ) = 0.0_rkx
          qsold(kpm) = 1.0_rkx
          do kp = 2 , kpm1
            qsold(kp) = (qsold(kp)-qs0k)/sqsk
            if ( (qsold(kp)-qsold(kp-1)) < eps ) then
              qsold(kp) = qsold(kp-1)+eps
            end if
          end do
          qs0(kth) = qs0k
          sqs(kth) = sqsk
          !-----------------------------------------------------------
          qsnew(1  ) = 0.0_rkx
          qsnew(kpm) = 1.0_rkx
          dqs = 1.0_rkx/real(kpm-1,rkx)
          do kp = 2 , kpm1
            qsnew(kp) = qsnew(kp-1) + dqs
          end do
          y2p(1   ) = 0.0_rkx
          y2p(kpm ) = 0.0_rkx
          call spline(jtb,kpm,qsold,pold,y2p,kpm,qsnew,pnew,app,aqp)
          do kp = 1 , kpm
            ptbl(kp,kth) = pnew(kp)
          end do
          !-----------------------------------------------------------
        end do
        ! Coarse look-up table for t(p) from constant the----------
        p = pl - dp
        do kp = 1 , kpm
          p  = p + dp
          th = thl - dth
          do kth = 1 , kthm
            th    = th + dth
            if ( p <= 0.0_rkx ) then
              told(kth)   = th
              theold(kth) = th
            else
              ape   = (100000.0_rkx/p)**rovcp
              denom = th - a4*ape
              if ( denom > eps ) then
                qs = pq0/p*exp(a2*(th-a3*ape)/denom)
              else
                qs = 0.0_rkx
              end if
              ! qs = pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
              told(kth) = th / ape
              theold(kth) = th*exp(eliwv*qs/(cpd*told(kth)))
            end if
          end do
          the0k = theold(1)
          sthek = theold(kthm) - theold(1)
          theold(1   ) = 0.0_rkx
          theold(kthm) = 1.0_rkx
          do kth = 2 , kthm1
            theold(kth) = (theold(kth)-the0k)/sthek
            if ( (theold(kth)-theold(kth-1)) < eps ) then
              theold(kth) = theold(kth-1) +  eps
            end if
          end do
          the0(kp) = the0k
          sthe(kp) = sthek
          !-----------------------------------------------------
          thenew(1  )  = 0.0_rkx
          thenew(kthm) = 1.0_rkx
          dthe         = 1.0_rkx/real(kthm-1,rkx)
          rdthe        = 1.0_rkx/dthe
          do kth = 2 , kthm1
            thenew(kth) = thenew(kth-1) + dthe
          end do
          y2t(1   ) = 0.0_rkx
          y2t(kthm) = 0.0_rkx
          call spline(jtb,kthm,theold,told,y2t,kthm,thenew,tnew,apt,aqt)
          do kth = 1 , kthm
            ttbl(kth,kp) = tnew(kth)
          end do
          !------------------------------------------------------
        end do
      end subroutine table_fill

    end subroutine otlift

end module mod_capecin
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

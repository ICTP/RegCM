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

module mod_cu_shallow
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_cu_common
  use mod_dynparam
  use mod_runparams , only : iqv , dtcum , dt
  use mod_regcm_types

  implicit none

  private

  real(rkx) , parameter :: rads = 50.0_rkx
  real(rkx) , parameter :: pcut = 400.0_rkx
  real(rkx) , parameter :: c0 = 0.0_rkx

  public :: shallcu

  contains

  subroutine shallcu(m2c)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c

    real(rkx) , dimension(kz) :: t , q , p , outts , tns , qns , outqs
    real(rkx) :: xmb
    real(rkx) :: dtime , ter11 , pret , psur

    integer(ik4) :: i , j , kbmax , k , kk , ier

    kbmax = ((kz*3)/4-3)
    dtime = dtcum

    do i = ici1 , ici2
      do j = jci1 , jci2
        !
        ! Prepare input, erase output
        !
        do k = 1 , kz
          kk = kzp1-k
          t(k) = m2c%tas(j,i,kk)
          q(k) = m2c%qxas(j,i,kk,iqv)
          if ( idynamic == 3 ) then
            tns(k) = t(k) + m2c%tten(j,i,kk)*dt
            qns(k) = q(k) + m2c%qxten(j,i,kk,iqv)*dt
          else
            tns(k) = t(k) + m2c%tten(j,i,kk)/m2c%psb(j,i)*dt
            qns(k) = q(k) + m2c%qxten(j,i,kk,iqv)/m2c%psb(j,i)*dt
          end if
          p(k) = m2c%pas(j,i,kk)*d_r100
          psur = m2c%psf(j,i)
          outts(k) = d_zero
          outqs(k) = d_zero
        end do
        pret = d_zero
        ter11 = m2c%ht(j,i)*regrav
        if ( ter11 <= d_zero ) ter11 = 1.e-05_rkx
        !
        ! Call cumulus parameterization
        !
        ier = 0
        call shallow(t,q,ter11,tns,qns,p,kz,pret,p,outts,outqs,dtcum, &
                     kbmax,pcut,c0,psur,ier,rads,xmb)
        if ( xmb > d_zero ) then
          total_precip_points = total_precip_points+1
          cu_prate(j,i) = pret/dt
        end if
        do k = 1 , kz
          kk = kzp1 - k
          cu_tten(j,i,kk) = outts(k)
          cu_qten(j,i,kk,iqv) = outqs(k)
        end do
      end do
    end do
  end subroutine shallcu

  subroutine shallow(ti,qi,z1,tio,qio,pio,klev,pre,pi,outtem,outq, &
                  dtime,kbmax,pcut,c0,psur,ier,rads,xmb)
    implicit none
      integer(ik4) , intent(in) :: klev , kbmax
      real(rkx) , intent(in) :: z1 , dtime , pcut , c0 , psur , rads
      real(rkx) , intent(out) :: xmb
      integer , intent(out) :: ier
      real(rkx) , intent(out) :: pre
      !
      ! Environmental properties before large scale forcing,
      ! only t,qe and p are needed for input!
      !
      real(rkx) , dimension(kz) , intent(in) :: pi , ti , qi
      !
      ! Environmental properties after large scale forcing,
      ! to,qo,po are needed as input!
      !
      real(rkx) , dimension(klev) , intent(in) :: tio , qio , pio
      real(rkx) , dimension(klev) , intent(out) :: outtem , outq
      !
      ! Environmental properties before large scale forcing
      !
      real(rkx) , dimension(kz) :: he , hes , qes , z , tv
      real(rkx) , dimension(kz) :: p , t , qe
      !
      ! Cloud updraft properties , which are calculated for
      ! every cloud (i) and stored. kz clouds are possible.
      !
      real(rkx) , dimension(kz) :: hc , zu , qrc , pw , zuo , pwo , qsc
      !
      ! Changed properties after simulated cloud influence,
      ! necessary for kernel calculations!
      !
      real(rkx) , dimension(kz) :: xt , xhe , xqe , xhes , xqes , xz
      real(rkx) , dimension(kz) :: xtv , xhh , xqsc
      real(rkx) , dimension(kz) :: xzu , xpw , xhc , xqrc
      real(rkx) :: xa
      !
      ! Various work dimensions, cloudworkfunctions (a,aold),
      ! kernels (xk), detrainment of updraft (cd)
      !
      real(rkx) :: aa , xx , f , aold
      real(rkx) :: xk
      real(rkx) , dimension(kz) :: hh , della
      real(rkx) , dimension(kz) :: cdd
      !
      ! Downdraft dimensions
      !
      real(rkx) :: edt
      real(rkx) , dimension(kz) :: zd , qrcd , hcd , pwd
      real(rkx) , dimension(kz) :: zdo , pwdo
      integer(ik4) :: kmin , itest
      integer(ik4) :: ktop , ktopx
      !
      ! DIMENSIONS FOR FEEDBACK ROUTINE, INCLUDING FEEDBACK
      ! OUTPUT TO DRIVER ROUTINE (OUTTEM,OUTQ,XMB,XMFLUX)!
      !
      real(rkx) :: rad
      real(rkx) , dimension(kz) :: dellt , dellq
      real(rkx) , dimension(kz) :: xmflux , het
      real(rkx) :: ax , xax , dh , dq , e , entnet , hkb
      real(rkx) :: mbdt , xhkb , xqkb
      real(rkx) :: toshall
      real(rkx) :: psurf , pbcdif , qkb
      integer(ik4) :: k , lloop , iph
      integer(ik4) :: kbhe , k22 , kb , kbcon , kbxx

      real(rkx) , parameter :: ht1 = wlhv/cpd
      real(rkx) , parameter :: ht2 = wlhs/cpd
      real(rkx) , parameter , dimension(2) :: be = [ ep2*ht1/c287 , &
                                                      ep2*ht2/c287 ]
      real(rkx) , parameter , dimension(2) :: ae = [ be(1)/tzero+log(c1es) , &
                                                      be(2)/tzero+log(c1es) ]

      rad = rads
      xx = 0.2_rkx/rad

      toshall = 50.0_rkx
      lloop = 0
      ier = 0
      mbdt = dtime*5.0e-3_rkx
      pre = d_zero

      do lloop = 1 , 2
        if ( lloop == 2 ) then
          psurf = psur
          aold = aa
          aa = d_zero
          xmb = d_zero
          f = d_zero
          edt = d_zero
          kmin = -1
          itest = -1
          do k = 1 , kz
            xmflux(k) = d_zero
            z(k) = d_zero
            t(k) = tio(k)
            qe(k) = qio(k)
            p(k) = pio(k)
          end do
          do k = 1 , kz
            zu(k) = d_zero
            zd(k) = d_zero
            hcd(k) = d_zero
            qrcd(k) = d_zero
            hc(k) = d_zero
            qrc(k) = d_zero
          end do
        else
          psurf = psur
          aa = d_zero
          f = d_zero
          xmb = d_zero
          edt = d_zero
          xk = d_zero

          do k = 1 , kz
            t(k) = ti(k)
            qe(k) = qi(K)
            p(k) = pi(k)
            della(k) = d_zero
            outtem(k) = d_zero
            outq(k) = d_zero
            hh(k) = d_zero
            cdd(k) = 1.0_rkx*xx
          end do

          do k = 1 , kz
            pwd(k) = d_zero
            pw(k) = d_zero
            dellt(k) = d_zero
            dellq(k) = d_zero
            zuo(k) = d_zero
            pwo(k) = d_zero
            pwdo(k) = d_zero
            zdo(k) = d_zero
            zu(k) = d_zero
            zd(k) = d_zero
            hcd(k) = d_zero
            qrcd(k) = d_zero
            hc(k) = d_zero
            qrc(k) = d_zero
          end do
        end if
        !
        ! Calculate environmental conditions
        !
        do k = 1 , kz
          iph = 1
          if ( t(k) <= toshall ) iph = 2
          e = exp(ae(iph)-be(iph)/t(k))
          qes(k) = ep2 * e/(d_100*p(k)-(d_one-ep2)*e)
          if ( qes(k) <= minqq ) qes(k) = minqq
          if ( qe(k) > qes(k) )qe(k) = qes(k)
          tv(k) = t(k)+ep1*qe(k)*t(k)
        end do
        !
        ! Heights
        !
        call heipre(p,z,tv,z1,psurf)
        !
        ! Determine moist static energy
        !
        call moiene(t,qe,z,he)
        call moiene(t,qes,z,hes)
        do k = 1 , kz
          if ( he(k) >= hes(k) ) he(k) = hes(k)
          if ( k < kz ) het(k) = d_half*(he(k)+he(k+1))
        end do
        het(kz) = he(kz)
        !
        ! Determine level with highest moist static energy content.
        ! First estimate a level in lower troposphere with high
        ! moist static energy.
        !
        kbhe = 1
        k22 = 1
        call minim(hes,k22,kbhe)
        call maxim(het,kbhe,k22)
        !
        ! Decide for convective cloud base
        !
        kbloop: do
          hkb = d_half*(he(k22)+he(k22+1))
          qkb = d_half*(qe(k22)+qe(k22+1))
          kb = k22
          kbcon = kb
          findbase: do
            dh = d_half * (hes(kbcon)+hes(kbcon+1))
            if ( hkb < dh ) then
              kbcon = kbcon+1
              if ( kbcon >= kz-1 ) then
                pbcdif = p(kb)-p(kbcon)
                pbcdif = d_half*(p(kb)+p(kb+1))-p(kbcon)
                if ( pbcdif > 50.0_rkx ) then
                  k22 = k22+1
                  if ( k22 > kbmax ) return
                  cycle kbloop
                end if
              end if
              cycle findbase
            end if
            exit findbase
          end do findbase
          exit kbloop
        end do kbloop
        !
        ! Static control
        !
        ! Moist static energy profile inside cloud
        ! and entrainment rate, negative entrainment rates
        ! are not possible.
        !
        call entrs(kbcon,he,hh,hes,xx,kz,z,hkb,ktop)
        if ( ktop < 2 ) return
        !
        ! Final moist static energy inside cloud
        !
        do k = 1 , kz
          hc(k) = hh(k)
        end do
        !
        ! Normalized mass flux profile
        !
        !  entnet is net entrainment minus detrainment which is zero for
        !  shallow convection
        !
        entnet = d_zero
        call zunc(k22,zu,kbcon,entnet,z,kz,ktop)
        !
        ! Moisture properties
        !
        call precip(kb,kbcon,kz,xx,hh,hes,t,qe,qes,pw,qsc,qrc, &
                    z,p,qkb,pcut,c0,zu,ktop)
        do k = 1 , kz
          hh(k) = d_zero
        end do
        !
        ! Store values, if lloop equals one!!!
        !
        if ( lloop == 1 ) then
          do k = 1 , kz
            zuo(k) = zu(k)
            pwo(k) = pw(k)
          end do
        end if

        call cloudws(hc,qes,hes,zu,z,t,kz,ax,kbcon,ktop)
        if ( ax <= d_zero ) then
          aa = d_zero
          ktop = 1
        else
          aa = ax
        end if
        if ( lloop > 1 ) then
          if ( (aold < dlowval ) .or. (aa < dlowval) ) then
            f = d_zero
          else
            f = (aa-aold)/dtime
          end if
        end if
        !
        ! Loop to calculate kernels. for every cloud calculate:
        !
        !    1.  Changes that cloud would do per unit mass flux
        !        if in existence to t,qe,he fields,
        !    2.  New environmental conditions,
        !    3.  New cloud work functions with modified environment,
        !    4.  Kernels
        !
        if ( lloop == 1 ) then
          if ( ktop < 2 ) return
          !
          ! Initialize old values first
          !
          xa = d_zero
          do k = 1 , kz
            xt(k) = t(k)
            xqe(k) = qe(k)
            xhe(k) = he(k)
            xhh(k) = d_zero
          end do
          xqrc = d_zero
          xpw = d_zero
          xzu = d_zero
          xhc = d_zero
          !
          ! Changes in moist static energy and moisture
          !
          call kerhels(he,xx,zu,hkb,hc,hh,p,z,kbcon,kz,xhe, &
                         k22,xhkb,1,cdd,psurf,ktop)
          call kerhels(qe,xx,zu,qkb,qrc,della,p,z,kbcon,kz,xqe, &
                         k22,xqkb,2,cdd,psurf,ktop)
          !
          ! Temperature
          !
          do k = 1 , kz
            dh = hh(k)
            dq = della(k)
            xt(k) = (mbdt/cpd)*(dh-wlhv*dq)+t(k)
            dellq(k) = dq
            dellt(k) = (d_one/cpd)*(dh-wlhv*dq)
            hh(k) = d_zero
            della(k) = d_zero
          end do
          !
          ! Virtual temperature for heights
          !
          do k = 1 , kz
            if ( xqe(k) < minqq ) xqe(k) = minqq
            iph = 1
            if ( xt(k) <= toshall ) iph = 2
            e = exp(ae(iph)-be(iph)/xt(k))
            xqes(k) = ep2*e/(d_100*p(k)-(d_one-ep2)*e)
            if ( xqes(k) < minqq ) xqes(k) = minqq
            if ( xqe(k) > xqes(k) ) xqe(k) = xqes(k)
            xtv(k) = xt(k) + ep1*xqe(k)*xt(k)
          end do
          !
          ! Heights
          !
          call heipre(p,xz,xtv,z1,psurf)
          !
          ! Determine saturation moist static energy
          !
          call moiene(xt,xqes,xz,xhes)
          do k = 1 , kz
            if ( xhe(k) >= xhes(k) ) xhe(k) = xhes(k)
          end do
          !
          ! Moist static energy profile inside cloud
          ! and entrainment rate, negative entrainment rates
          ! are not possible.
          !
          call entrs(kbcon,xhe,xhh,xhes,xx,kz,xz,xhkb,ktopx)
          if ( ktopx < 2 ) return
          !
          ! Final moist static energy inside cloud
          !
          do k = 1 , kz
            xhc(k) = xhh(k)
          end do
          !
          ! Normalized mass flux profile
          !
          call zunc(k22,xzu,kbcon,entnet,xz,kz,ktopx)
          !
          ! Moisture calculations necessary for downdrafts
          !
          call precip(kb,kbcon,kz,xx,xhh,xhes,xt,xqe,xqes,xpw, &
                      xqsc,xqrc,xz,p,xqkb,pcut,c0,xzu,ktopx)
          do k = 1 , kz
            xhh(k) = d_zero
          end do
          !
          ! Cloud work function
          !
          if ( ktopx > 2 ) then
            call cloudws(xhc,xqes,xhes,xzu,xz,xt,kz,xax,kbcon,ktopx)
            !
            ! No negative generation of buoyant energy, PLEEEAAAASSSEEE!
            !
            if ( xax <= d_zero ) then
              ktopx = 1
              xa = d_zero
            else
              xa = xax
              xax = d_zero
            end if
          end if
          if ( xa < dlowval .and. aa < dlowval ) then
            xk = d_zero
            return
          else
            xk = (xa-aa)/mbdt
            if ( xk > d_zero ) then
              xk = -1.0_rkx
            end if
          end if
          kbxx = kb
          kb = 1
       else
         !
         ! Calculate massfluxes
         !
         ier = 0
         if ( abs(xk) > dlowval ) xmb = -f/xk
         if ( xmb < d_zero ) xmb = d_zero
         if  ( kbxx == 0 .or. kbxx > kbmax ) then
           ier = 999
           return
         end if
       end if
     end do
     !
     ! Calculate output
     !
     call araouts(xmflux,xmb,zuo,dellt,dellq,kz,outtem,outq,ier,pre)
     if ( ier > 0 ) then
       do k = 1 , kz
         outtem(k) = d_zero
         outq(k) = d_zero
       end do
       pre = d_zero
       ier = 0
       return
     end if

  contains

    subroutine araouts(xmc,xmb,zu,delt,delq,kz,outtem,outq,ier,pre)
      implicit none
      integer(ik4) , intent(in) :: kz
      real(rkx) , intent(inout) :: xmb
      real(rkx) , dimension(kz) , intent(in) :: delt , delq , zu
      real(rkx) , dimension(kz) , intent(out) :: xmc , outtem , outq
      integer(ik4) , intent(out) :: ier
      real(rkx) , intent(out) :: pre
      integer(ik4) :: k , ibcout
      real(rkx) :: outtes

      ier = 0
      ibcout = 0
      pre = d_zero
      do k = 1 , kz
        outtem(k) = d_zero
        outq(k) = d_zero
        xmc(k) = d_zero
      end do
      if ( xmb > d_zero ) ibcout = 1
      if ( ibcout == 0 ) then
        ier = 888
        return
      end if
      !
      ! First calculate temperature and moisture changes
      !
      do k = 1 , kz
        outtes = delt(k)*xmb*86400.0_rkx
        if ( (outtes > 500.0_rkx) .or. (outtes < -200.0_rkx) ) then
          xmb = d_zero
          ier = 777
          return
        endif
        outtem(k) = outtem(k) + delt(k)*xmb
        outq(k) = outq(k) + delq(k)*xmb
      end do
      !
      ! How about them mass fluxes??
      !
      do k = 1 , kz
        xmc(k) = xmc(k) + zu(k)*xmb
      end do
    end subroutine araouts

    subroutine cloudws(hc,qes,hes,zu,z,tempp,kz,ax,kbcon,ktop)
      implicit none
      integer(ik4) , intent(in) :: kz , kbcon
      real(rkx) , dimension(kz) , intent(in) :: z , qes , hes , tempp
      real(rkx) , dimension(kz) , intent(in) :: hc , zu
      integer(ik4) , intent(in) :: ktop
      real(rkx) , intent(out) :: ax
      integer(ik4) :: k
      real(rkx) :: gamma1 , gamma2 , dhh , zzu , dt , dg , dh , dz , aa

      ax = d_zero
      !
      ! Calculate cloud work function for updraft
      !
      do  k = kbcon , ktop
        gamma1 = (wlhv/cpd)*(wlhv/(rwat*(tempp(k)**2)))*qes(k)
        gamma2 = (wlhv/cpd)*(wlhv/(rwat*(tempp(k+1)**2)))*qes(k+1)
        dhh = hc(k)
        zzu = zu(k)
        dt = d_half*(tempp(k)+tempp(k+1))
        dg = d_half*(gamma1+gamma2)
        dh = d_half*(hes(k)+hes(k+1))
        dz = z(k+1)-z(k)
        aa = dz*(egrav/(cpd*dt)) * zzu * ((dhh-(dh))/(d_one+dg))
        ax = ax+aa
      end do
    end subroutine cloudws
    !
    ! This subroutine calculates incloud moist static energy
    !
    subroutine entrs(kbc,h,hc,hsat,ent,kz,p,hkb,ktop)
      implicit none
      integer(ik4) , intent(in) :: kbc , kz
      real(rkx) , dimension(kz) , intent(in) :: h , hsat , p
      real(rkx) , dimension(kz) , intent(out) :: hc
      real(rkx) , intent(in) :: hkb , ent
      integer(ik4) , intent(inout) :: ktop
      real(rkx) :: dz , depth , dby
      integer(ik4) :: k , kend

      ktop = 1
      do k = 1 , kbc
        hc(k) = hkb
      end do
      do k = kbc+1 , kz-1
        hc(k) = d_half*(hsat(k)+hsat(k+1))
      end do
      hc(kz) = hsat(kz)

      kend = kz-1
      if ( kbc+1 > kend ) return
      do k = kbc+1 , kend
        dz = d_half*(p(k-1)+p(k+1))
        hc(k) = (hc(k-1)*(d_one - d_half*dz*ent) + ent*dz*h(k)) / &
                (d_one + d_half*ent*dz)
        dby = hc(k) - d_half*(hsat(k)+hsat(k+1))
        depth = p(k) - p(kbc+1)
        if ( dby < d_zero .or. depth > 3000.0_rkx ) then
          ktop = k-1
          exit
        end if
      end do
    end subroutine entrs
    !
    ! Subroutine kerhel calculates the changes a cloud with
    ! top height lp, base kb and entrainment rate r would do to
    ! the environment per unit mass (del*dp/g).
    !
    subroutine kerhels(var,r,zu,hkb,hc,della,p,z,kb,kz, &
                       xvar,kbeg,xhkb,ich,cd,psu,ktop)
      implicit none
      integer(ik4) , intent(in) :: kz , kbeg , kb , ich
      real(rkx) , intent(in) :: hkb , r , psu
      real(rkx) , intent(out) :: xhkb
      real(rkx) , dimension(kz) , intent(in) :: p , z , var
      real(rkx) , dimension(kz) , intent(out) :: xvar , della
      real(rkx) , dimension(kz) , intent(in) :: hc , zu , cd
      integer(ik4) , intent(in) :: ktop
      integer(ik4) :: k , lpt
      real(rkx) :: dp , dz , dv1 , dv2 , dv3 , zu1 , zu2 , detup

      do k = 1 , kz
        xvar(k) = var(k)
        della(k) = d_zero
      end do
      !
      ! Calculate changes, first surface, last top
      !
      if ( kbeg == 1 ) then
        dp = 50.0_rkx*(psu-p(2))
        dv1 = d_half*(var(1)+var(2))
        dv2 = var(1)
        della(kbeg) = (dv1-dv2)*egrav/dp
      else
        k = kbeg
        dp = 50.0_rkx*(p(k-1)-p(k+1))
        dv1 = d_half*(var(k)+var(k+1))
        dv2 = var(k)
        della(kbeg) = (dv1-dv2)*egrav/dp
      endif
      della(kbeg) = d_zero

      if ( kbeg+1 < ktop ) then
        do k = kbeg+1 , ktop-1
          zu1 = zu(k)
          zu2 = zu(k-1)
          dv1 = d_half*(var(k)+var(k+1))
          dv2 = var(k)
          dv3 = d_half*(var(k)+var(k-1))
          !
          ! Specifiy detrainment of downdraft, has to be consistent
          ! with zd calculations in soundd.
          !
          dz = d_half*(z(k+1)-z(k-1))
          !
          ! Changed due to subsidence and entrainment
          !
          detup = (d_half*(hc(k-1)+hc(k))-dv2)*r*cd(k)*dz*zu2
          if ( k <= kb ) detup = d_zero
          dp = 50.0_rkx*(p(k-1)-p(k+1))
          della(k) = ((zu1)*(dv1-dv2)+(zu2)*(dv2-dv3)+detup)*egrav/dp
        end do
      end if
      !
      ! Cloud top
      !
      lpt = ktop
      zu1 = zu(lpt-1)
      dz = d_half*(z(lpt)-z(lpt-1))
      dp = d_100*(p(lpt-1)-p(lpt))
      dv1 = d_half*(var(lpt)+var(lpt-1))
      della(lpt) = zu1*(hc(lpt)-dv1)*egrav/dp
      !
      ! Final changed variable per unit mass flux
      !
      do k = 1 , ktop
        xvar(k) = della(k)*mbdt + var(k)
        if ( xvar(k) <= d_zero ) then
          xvar(k) = minqq
          della(k) = d_zero
        end if
      end do
      xhkb = della(kbeg)*mbdt + hkb
      if ( ich <= 0 ) then
        if ( della(kbeg) >= d_zero ) then
          della(kbeg) = d_zero
          xhkb = della(kbeg)*mbdt + hkb
        end if
      end if
    end subroutine kerhels

    subroutine heipre(p,z,t,z1,psurf)
      implicit none
      real(rkx) , dimension(kz) , intent(in) :: p , t
      real(rkx) , dimension(kz) , intent(out) :: z
      real(rkx) , intent(in) :: psurf , z1
      real(rkx) :: zstart , tvbar
      integer(ik4) :: k
      zstart = z1-(log(p(1))-log(psurf))*rgas*t(1)/egrav
      z(1) = zstart
      do k = 2 , kz
        tvbar = d_half * (t(k)+t(k-1))
        z(k) = z(k-1) - (log(p(k))-log(p(k-1)))*rgas*tvbar/egrav
      end do
    end subroutine heipre

    subroutine moiene(t,q,z,h)
      implicit none
      real(rkx) , dimension(kz) , intent(in) :: t , q , z
      real(rkx) , dimension(kz) , intent(out) :: h
      integer(ik4) :: k
      do k = 1 , kz
        h(k) = egrav*z(k) + cpd*t(k) + wlhv*q(k)
      end do
    end subroutine moiene

    subroutine minim(array,ks,kmin)
      implicit none
      integer(ik4) , intent(in) :: ks
      integer(ik4) , intent(out) :: kmin
      real(rkx) , dimension(kz) :: array
      integer(ik4) :: k
      real(rkx) :: x
      kmin = ks
      x = array(ks)
      do k = ks+1 , kz
        if ( array(k) < x ) then
          x = array(k)
          kmin = k
        end if
      end do
    end subroutine minim

    subroutine maxim(array,ke,kmax)
      implicit none
      integer(ik4) , intent(in) :: ke
      integer(ik4) , intent(out) :: kmax
      real(rkx) , dimension(kz) :: array
      integer(ik4) :: k
      real(rkx) :: x
      kmax = 1
      x = array(1)
      do k = 2 , ke
        if ( array(k) >= x ) then
          x = array(k)
          kmax = k
        end if
      end do
    end subroutine maxim
    !
    ! Normalized updraft mass flux function
    !
    subroutine zunc(kbeg,zu,kb,r,z,kz,ktop)
      implicit none
      integer(ik4) , intent(in) :: kz , kbeg , kb
      real(rkx) , dimension(kz) , intent(out) :: zu
      real(rkx) , dimension(kz) , intent(in) :: z
      real(rkx) , intent(in) :: r
      integer(ik4) :: ktop
      integer(ik4) :: k , lpt
      real(rkx) :: dz
      !
      ! Below updraft air originating level (zu)
      !
      do k = 1 , kbeg
        zu(k) = d_zero
      end do
      !
      ! Between zu and level of free convectioN
      !
      do k = kbeg , kb
        zu(k) = d_one
      end do
      !
      ! Between level of free convection and cloud top
      !
      do k = kb+1 , ktop-1
        dz = d_half * (z(k+1)-z(k-1))
        zu(k) = zu(k-1) * (d_one + r*dz)
      end do
      !
      ! Cloud top
      !
      lpt = ktop
      dz = d_half * (z(lpt)-z(lpt-1))
      zu(lpt) = zu(lpt-1) * (d_one + r*dz)
    end subroutine zunc
    !
    ! Calculate precipitation stepwise:
    !
    !     1. Use steady state plume equation to calculate qc
    !     2. With saturation condition, calculate how much moisture
    !        could be in updraft (qrc).
    !     3. Calculate liquid water and precipitation (pw) at that
    !        level (qc-qrc).
    !     4. Determine new qc.
    !     5. Go back to 1..
    !
    subroutine precip(kb,kbcon,kz,r,hc,hes,t,qe, &
                      qes,pw,qc,qrc,z,p,qkb,pcut,c0,zu,ktop)
      implicit none
      integer(ik4) , intent(in) :: kbcon , kz , kb
      real(rkx) , dimension(kz) , intent(in) :: hc , hes , t , z , p
      real(rkx) , dimension(kz) , intent(in) :: qe , qes
      real(rkx) , dimension(kz) , intent(out) :: qc
      real(rkx) , dimension(kz) , intent(in) :: zu
      real(rkx) , dimension(kz) , intent(out) :: pw , qrc
      real(rkx) , intent(in) :: c0 , pcut , qkb , r
      integer(ik4) :: ktop
      real(rkx) :: dh , dz , agamma , qrch
      integer :: k
      !
      ! Erase
      !
      do k = 1 , kz
        pw(k) = d_zero
        qc(k) = d_zero
      end do
      !
      ! Initialize
      !
      do k = 1 , kb
        qrc(k) = d_half*(qe(k)+qe(k+1))
      end do
      do k = kb , kbcon
        qrc(k) = qkb
      end do
      do k = ktop , kz
        qrc(k) = qes(k)
      end do
      !
      ! Boundary condition
      !
      qc(kbcon) = qkb
      pw(kbcon) = d_zero
      !
      ! Start loop
      !
      do k = kbcon+1 , ktop
        if ( k == ktop ) then
          dh = hc(k)-hes(k)
          dz = d_half*(z(k)-z(k-1))
        else
          dh = hc(k)-d_half*(hes(k)+hes(k+1))
          dz = -d_half*z(k-1)+d_half*z(k+1)
        end if
        agamma = (wlhv/cpd)*(wlhv/(rwat*(t(k)**2)))*qe(k)
        !
        ! 1.
        !
        qc(k) = (qc(k-1)*(d_one-r*dz*d_half) + r*dz*qe(k))/(d_one+d_half*dz*r)
        !
        ! 2.
        !
        qrch = qes(k)+(d_one/wlhv)*(agamma/(d_one+agamma))*dh
        !
        ! Liquid water content+water vapor (qrc)
        !
        qrc(k) = (qc(k)-qrch)/(d_one+c0*dz)+qrch
        !
        ! 3. Precip
        !
        pw(k) = c0*dz*(qrc(k)-qrch)*zu(k)
        !
        ! If p<pcut, all water is rained out!
        !
        if ( k == ktop ) then
          qrc(k) = (qc(k)-qes(k))/(d_one+c0*dz)+qes(k)
          pw(k) = (qrc(k)-qes(k))*c0*dz*zu(k)
          if ( p(ktop) < pcut ) then
            qrc(k) = qes(k)
            pw(k) = (qc(k)-qrc(k))*zu(k)
          end if
        end if
        !
        ! 4.+5.
        !
        qc(k) = qrc(k)
      end do
    end subroutine precip

  end subroutine shallow

end module mod_cu_shallow

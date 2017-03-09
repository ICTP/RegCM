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
  use mod_constants
  use mod_cu_common
  use mod_dynparam
  use mod_runparams , only : iqv , dtcum , dt

  implicit none

  private

  integer(ik4) , parameter :: knum = 1
  real(rkx) , parameter :: rads = 50.0_rkx
  real(rkx) , parameter :: pcut = 400.0_rkx
  real(rkx) , parameter :: c0 = 0.0_rkx

  public :: shallcu

  contains

  subroutine shallcu(m2c)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c

    real(rkx) , dimension(kz) :: t , q , p , outts , tns , qns , outqs
    real(rkx) , dimension(knum) :: xmb
    real(rkx) :: mbdt , dtime , ter11 , pret , psur

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
          tns(k) = t(k) + avg_tten(j,i,kk)*dt
          qns(k) = q(k) + avg_qten(j,i,kk,iqv)*dt
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
        if ( xmb(1) > d_zero ) then
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
      real(rkx) , dimension(knum) , intent(out) :: xmb
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
      real(rkx) , dimension(kz,knum) :: hc , zu , qrc , pw , zuo , pwo , qsc
      !
      ! Changed properties after simulated cloud influence,
      ! necessary for kernel calculations!
      !
      real(rkx) , dimension(kz) :: xt , xhe , xqe , xhes , xqes , xz
      real(rkx) , dimension(kz) :: xtv , xhh , xqsc
      real(rkx) , dimension(kz,knum) :: xzu , xpw , xhc , xqrc
      real(rkx) , dimension(knum) :: xa
      !
      ! Various work dimensions, cloudworkfunctions (a,aold),
      ! kernels (xk), detrainment of updraft (cd)
      !
      real(rkx) , dimension(knum) :: aa , xx , xxo , f , aold
      real(rkx) , dimension(knum,knum) :: xk
      real(rkx) , dimension(kz) :: hh , della
      real(rkx) , dimension(kz,knum) :: cdd
      !
      ! Downdraft dimensions
      !
      real(rkx) , dimension(knum) :: edt , edto , edtoo
      real(rkx) , dimension(kz,knum) :: zd , qrcd , hcd , pwd , xzd , xqrcd
      real(rkx) , dimension(kz,knum) :: xhcd , xpwd , zdo , pwdo
      integer(ik4) , dimension(knum) :: kmin , itest , itestx , itest2
      integer(ik4) , dimension(knum) :: kminx , ktop , ktopx
      !
      ! DIMENSIONS FOR FEEDBACK ROUTINE, INCLUDING FEEDBACK
      ! OUTPUT TO DRIVER ROUTINE (OUTTEM,OUTQ,XMB,XMFLUX)!
      !
      real(rkx) , dimension(knum) :: rad
      real(rkx) , dimension(kz,knum) :: dellt , dellq
      real(rkx) , dimension(kz) :: xmflux , het
      real(rkx) :: ax , xax , dh , dq , e , entnet , hkb
      real(rkx) :: mbdt , xhkb , xqkb
      real(rkx) :: toshall
      real(rkx) :: psurf , pbcdif , qkb
      integer(ik4) :: k , kk , ll , lloop , iph , lp , lpx
      integer(ik4) :: kbhe , k22 , kb , kbcon , kbxx , knumen

      real(rkx) , parameter :: ht1 = wlhv/cpd
      real(rkx) , parameter :: ht2 = wlhs/cpd
      real(rkx) , parameter , dimension(2) :: be = (/ ep2*ht1/c287 , &
                                                      ep2*ht2/c287 /)
      real(rkx) , parameter , dimension(2) :: ae = (/ be(1)/tzero+log(c1es) , &
                                                      be(2)/tzero+log(c1es) /)

      rad(1) = rads
      xx(1) = 0.2_rkx/rad(1)

      toshall = 50.0_rkx
      lloop = 0
      ier = 0
      mbdt = dtime*5.0e-3_rkx
      pre = d_zero

      do lloop = 1 , 2
        if ( lloop == 2 ) then
          psurf = psur
          do k = 1 , knum
            aold(k) = aa(k)
            aa(k) = d_zero
            xmb(k) = d_zero
            f(k) = d_zero
            edt(k) = d_zero
            kmin(k) = -1
            itest(k) = -1
          end do
          do k = 1 , kz
            xmflux(k) = d_zero
            z(k) = d_zero
            t(k) = tio(k)
            qe(k) = qio(k)
            p(k) = pio(k)
          end do
          do ll = 1 , knum
            do k = 1 , kz
              zu(k,ll) = d_zero
              zd(k,ll) = d_zero
              hcd(k,ll) = d_zero
              qrcd(k,ll) = d_zero
              hc(k,ll) = d_zero
              qrc(k,ll) = d_zero
            end do
          end do
        else
          psurf = psur
          do k = 1 , knum
            aa(k) = d_zero
            f(k) = d_zero
            xmb(k) = d_zero
            edt(k) = d_zero
            do ll = 1 , knum
              xk(k,ll) = d_zero
            end do
          end do

          do k = 1 , kz
            t(k) = ti(k)
            qe(k) = qi(K)
            p(k) = pi(k)
            della(k) = d_zero
            outtem(k) = d_zero
            outq(k) = d_zero
            hh(k) = d_zero
            do ll = 1 , knum
              cdd(k,ll) = 1.0_rkx*xx(ll)
            end do
          end do

          do ll = 1 , knum
            do k = 1 , kz
              pwd(k,ll) = d_zero
              pw(k,ll) = d_zero
              dellt(k,ll) = d_zero
              dellq(k,ll) = d_zero
              zuo(k,ll) = d_zero
              pwo(k,ll) = d_zero
              pwdo(k,ll) = d_zero
              zdo(k,ll) = d_zero
              zu(k,ll) = d_zero
              zd(k,ll) = d_zero
              hcd(k,ll) = d_zero
              qrcd(k,ll) = d_zero
              hc(k,ll) = d_zero
              qrc(k,ll) = d_zero
            end do
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
        call heipre(p,z,tv,z1,kz,1,psurf)
        !
        ! Determine moist static energy
        !
        call moiene(t,qe,z,he,kz,1)
        call moiene(t,qes,z,hes,kz,1)
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
        call minim(hes,kz,k22,kz,kbhe)
        call maxim(het,kz,1,kbhe,k22)

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
        knumen = knum
        !
        ! Static control
        !
        do lp = 1 , knumen
          !
          ! Moist static energy profile inside cloud
          ! and entrainment rate, negative entrainment rates
          ! are not possible.
          !
          call entrs(kbcon,he,hh,hes,xx(lp),kz,lp,kb,z,hkb,knum,ktop)
          !
          ! Final moist static energy inside cloud
          !
          do k = 1 , kz
            hc(k,lp) = hh(k)
          end do
          !
          ! Normalized mass flux profile
          !
          !  entnet is net entrainment minus detrainment which is zero for
          !  shallow convection
          !
          entnet = d_zero
          call zunc(k22,zu,kbcon,entnet,z,kz,lp,knum,ktop)
          !
          ! Moisture properties
          !
          call precip(cdd,kb,kbcon,kz,xx(lp),hh,hes,t,qe,qes,pw,qsc,qrc, &
                      z,lp,p,qkb,pcut,c0,zu,knum,ktop)
          do k = 1 , kz
            hh(k) = d_zero
          end do
        end do
        !
        ! Store values, if lloop equals one!!!
        !
        if ( lloop == 1 ) then
          do kk = 1 , knum
            do k = 1 , kz
              zuo(k,kk) = zu(k,kk)
              pwo(k,kk) = pw(k,kk)
            end do
          end do
        end if

        do lp = 1 , knumen
          if ( ktop(lp) < 2 ) exit
          call cloudws(lp,hc,qes,hes,zu,z,t,kz,ax,kbcon,knum,ktop)
          if ( ax <= d_zero ) then
            aa(lp) = d_zero
            ktop(lp) = 1
          else
            aa(lp) = ax
          end if
          if ( lloop > 1 ) then
            if ( (aold(lp) < dlowval ) .or. (aa(lp) < dlowval) ) then
              f(lp) = d_zero
            else
              f(lp) = (aa(lp)-aold(lp))/dtime
            end if
          end if
        end do
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
          do lp = 1 , knumen
            if ( ktop(lp) < 2 ) continue
            !
            ! Initialize old values first
            !
            do k = 1 , knum
              xa(k) = d_zero
            end do
            do k = 1 , kz
              xt(k) = t(k)
              xqe(k) = qe(k)
              xhe(k) = he(k)
              xhh(k) = d_zero
            end do
            do ll = 1 , kz
              do k = 1 , knum
                xqrc(ll,k) = d_zero
                xpw(ll,k) = d_zero
                xzu(ll,k) = d_zero
                xhc(ll,k) = d_zero
              end do
            end do
            !
            ! Changes in moist static energy and moisture
            !
            call kerhels(he,xx(lp),lp,zu,hkb,hc,hh,p,z,kbcon,kz,xhe,mbdt, &
                         k22,xhkb,1,cdd,psurf,knum,ktop)
            call kerhels(qe,xx(lp),lp,zu,qkb,qrc,della,p,z,kbcon,kz,xqe, &
                         mbdt,k22,xqkb,2,cdd,psurf,knum,ktop)
            !
            ! Temperature
            !
            do k = 1 , kz
              dh = hh(k)
              dq = della(k)
              xt(k) = (mbdt/cpd)*(dh-wlhv*dq)+t(k)
              dellq(k,lp) = dq
              dellt(k,lp) = (d_one/cpd)*(dh-wlhv*dq)
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
            call heipre(p,xz,xtv,z1,kz,1,psurf)
            !
            ! Determine saturation moist static energy
            !
            call moiene(xt,xqes,xz,xhes,kz,1)
            do k = 1 , kz
              if ( xhe(k) >= xhes(k) ) xhe(k) = xhes(k)
            end do
            do lpx = 1 , knumen
              !
              ! Moist static energy profile inside cloud
              ! and entrainment rate, negative entrainment rates
              ! are not possible.
              !
              call entrs(kbcon,xhe,xhh,xhes,xx(lpx),kz,lpx,kb,xz,xhkb, &
                         knum,ktopx)
              if ( ktopx(lpx) < 2 ) continue
              !
              ! Final moist static energy inside cloud
              !
              do k = 1 , kz
                xhc(k,lpx) = xhh(k)
              end do
              !
              ! Normalized mass flux profile
              !
              call zunc(k22,xzu,kbcon,entnet,xz,kz,lpx,knum,ktopx)
              !
              ! Moisture calculations necessary for downdrafts
              !
              call precip(cdd,kb,kbcon,kz,xx(lpx),xhh,xhes,xt,xqe,xqes,xpw, &
                          xqsc,xqrc,xz,lpx,p,xqkb,pcut,c0,xzu,knum,ktopx)
              do k = 1 , kz
                xhh(k) = d_zero
              end do
            end do
            !
            ! Cloud work function
            !
            do lpx = 1 , knumen
              if ( ktopx(lpx) > 2 ) then
                call cloudws(lpx,xhc,xqes,xhes,xzu,xz,xt,kz,xax, &
                             kbcon,knum,ktopx)
                !
                ! No negative generation of buoyant energy, PLEEEAAAASSSEEE!
                !
                if ( xax <= d_zero ) then
                  ktopx(lpx) = 1
                  xa(lpx) = d_zero
                else
                  xa(lpx) = xax
                  xax = d_zero
                end if
              end if
              if ( xa(lpx) < dlowval .and. aa(lpx) < dlowval ) then
                xk(lpx,lp) = d_zero
                exit
              else
                xk(lpx,lp) = (xa(lpx)-aa(lpx))/mbdt
                if ( lpx <= lp ) then
                  if ( xk(lpx,lp) > d_zero ) then
                    xk(lpx,lp) = -1.0_rkx
                  end if
                end if
              end if
            end do
          end do

          kbxx = kb
          kb = 1
       else
         !
         ! Calculate massfluxes
         !
         ier = 0
         if ( abs(xk(1,1)) > dlowval ) xmb(1) = -f(1)/xk(1,1)
         if ( xmb(1) < d_zero ) xmb(1) = d_zero
         if  ( kbxx == 0 .or. kbxx > kbmax ) then
           ier = 999
           return
         end if
       end if
     end do
     !
     ! Calculate output
     !
     call araouts(xmflux,xmb,zuo,dellt,dellq,kz,outtem,outq,ier,pre,knum)
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

    subroutine araouts(xmc,xmb,zu,delt,delq,kz,outtem,outq,ier,pre,knum)
      implicit none
      integer(ik4) , intent(in) :: kz , knum
      real(rkx) , dimension(knum) , intent(inout) :: xmb
      real(rkx) , dimension(kz,knum) , intent(in) :: delt , delq , zu
      real(rkx) , dimension(kz) , intent(out) :: xmc , outtem , outq
      integer(ik4) , intent(out) :: ier
      real(rkx) , intent(out) :: pre
      integer(ik4) :: k , lp , ibcout
      real(rkx) :: outtes

      ier = 0
      ibcout = 0
      pre = d_zero
      do k = 1 , kz
        outtem(k) = d_zero
        outq(k) = d_zero
        xmc(k) = d_zero
      end do
      do k = 1 , knum
        if ( xmb(k) > d_zero ) ibcout = 1
      end do
      if ( ibcout == 0 ) then
        ier = 888
        return
      end if
      !
      ! First calculate temperature and moisture changes
      !
      do k = 1 , kz
        do lp = 1 , knum
          outtes = delt(k,lp)*xmb(lp)*86400.0_rkx
          if ( (outtes > 500.0_rkx) .or. (outtes < -200.0_rkx) ) then
            xmb(lp) = d_zero
            ier = 777
            return
          endif
          outtem(k) = outtem(k) + delt(k,lp)*xmb(lp)
          outq(k) = outq(k) + delq(k,lp)*xmb(lp)
        end do
      end do
      !
      ! How about them mass fluxes??
      !
      do k = 1 , kz
        xmc(k) = d_zero
        do lp = 1 , knum
          xmc(k) = xmc(k) + zu(k,lp)*xmb(lp)
        end do
      end do
    end subroutine araouts

    subroutine cloudws(l1,hc,qes,hes,zu,z,tempp,kz,ax,kbcon,knum,ktop)
      implicit none
      integer(ik4) , intent(in) :: kz , kbcon , knum , l1
      real(rkx) , dimension(kz) , intent(in) :: z , qes , hes , tempp
      real(rkx) , dimension(kz,knum) , intent(in) :: hc , zu
      integer(ik4) , dimension(knum) , intent(in) :: ktop
      real(rkx) , intent(out) :: ax
      integer(ik4) :: k
      real(rkx) :: gamma1 , gamma2 , dhh , zzu , dt , dg , dh , dz , aa

      ax = d_zero
      !
      ! Calculate cloud work function for updraft
      !
      do  k = kbcon , ktop(l1)
        gamma1 = (wlhv/cpd)*(wlhv/(rwat*(tempp(k)**2)))*qes(k)
        gamma2 = (wlhv/cpd)*(wlhv/(rwat*(tempp(k+1)**2)))*qes(k+1)
        dhh = hc(k,l1)
        zzu = zu(k,l1)
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
    subroutine entrs(kbc,h,hc,hsat,ent,kz,lp,kb,p,hkb,knum,ktop)
      implicit none
      integer(ik4) , intent(in) :: kbc , kz , lp , kb , knum
      real(rkx) , dimension(kz) , intent(in) :: h , hsat , p
      real(rkx) , dimension(kz) , intent(out) :: hc
      real(rkx) , intent(in) :: hkb , ent
      integer(ik4) , dimension(knum) , intent(inout) :: ktop
      real(rkx) :: dz , depth , dby
      integer(ik4) :: k , kend

      ktop(lp) = 1

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
          ktop(lp) = k-1
          exit
        end if
      end do
    end subroutine entrs
    !
    ! Subroutine kerhel calculates the changes a cloud with
    ! top height lp, base kb and entrainment rate r would do to
    ! the environment per unit mass (del*dp/g).
    !
    subroutine kerhels(var,r,lp,zu,hkb,hc,della,p,z,kb,kz, &
                       xvar,mbdt,kbeg,xhkb,ich,cd,psu,knum,ktop)
      implicit none
      integer(ik4) , intent(in) :: kz , knum , kbeg , lp , kb , ich
      real(rkx) , intent(in) :: hkb , mbdt , r , psu
      real(rkx) , intent(out) :: xhkb
      real(rkx) , dimension(kz) , intent(in) :: p , z , var
      real(rkx) , dimension(kz) , intent(out) :: xvar , della
      real(rkx) , dimension(kz,knum) , intent(in) :: hc , zu , cd
      integer(ik4) , dimension(knum) , intent(in) :: ktop
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

      if ( kbeg+1 < ktop(lp) ) then
        do k = kbeg+1 , ktop(lp)-1
          zu1 = zu(k,lp)
          zu2 = zu(k-1,lp)
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
          detup = (d_half*(hc(k-1,lp)+hc(k,lp))-dv2)*r*cd(k,lp)*dz*zu2
          if ( k <= kb ) detup = d_zero
          dp = 50.0_rkx*(p(k-1)-p(k+1))
          della(k) = ((zu1)*(dv1-dv2)+(zu2)*(dv2-dv3)+detup)*egrav/dp
        end do
      end if
      !
      ! Cloud top
      !
      lpt = ktop(lp)
      zu1 = zu(lpt-1,lp)
      dz = d_half*(z(lpt)-z(lpt-1))
      dp = d_100*(p(lpt-1)-p(lpt))
      dv1 = (var(lpt))
      dv1 = d_half*(var(lpt)+var(lpt-1))
      della(lpt) = zu1*(hc(lpt,lp)-dv1)*egrav/dp
      !
      ! Final changed variable per unit mass flux
      !
      do k = 1 , ktop(lp)
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

    subroutine heipre(p,z,t,z1,kz,ks,psurf)
      implicit none
      integer(ik4) , intent(in) :: kz , ks
      real(rkx) , dimension(kz) , intent(in) :: p , t
      real(rkx) , dimension(kz) , intent(out) :: z
      real(rkx) , intent(in) :: psurf , z1
      real(rkx) :: zstart , tvbar
      integer(ik4) :: k
      zstart = z1-(log(p(1))-log(psurf))*rgas*t(1)/egrav
      z(ks) = zstart
      do k = ks+1 , kz
        tvbar = d_half * (t(k)+t(k-1))
        z(k) = z(k-1) - (log(p(k))-log(p(k-1)))*rgas*tvbar/egrav
      end do
    end subroutine heipre

    subroutine moiene(t,q,z,h,kz,ks)
      implicit none
      integer(ik4) , intent(in) :: kz , ks
      real(rkx) , dimension(kz) , intent(in) :: t , q , z
      real(rkx) , dimension(kz) , intent(out) :: h
      integer(ik4) :: k
      do k = ks , kz
        h(k) = egrav*z(k) + cpd*t(k) + wlhv*q(k)
      end do
    end subroutine moiene

    subroutine minim(array,kz,ks,ke,kmin)
      implicit none
      integer(ik4) , intent(in) :: kz , ks , ke
      integer(ik4) , intent(out) :: kmin
      real(rkx) , dimension(kz) :: array
      integer(ik4) :: k
      real(rkx) :: x
      kmin = ks
      x = array(ks)
      do k = ks+1 , ke
        if ( array(k) < x ) then
          x = array(k)
          kmin = k
        end if
      end do
    end subroutine minim

    subroutine maxim(array,kz,ks,ke,kmax)
      implicit none
      integer(ik4) , intent(in) :: kz , ks , ke
      integer(ik4) , intent(out) :: kmax
      real(rkx) , dimension(kz) :: array
      integer(ik4) :: k
      real(rkx) :: x
      kmax = ks
      x = array(ks)
      do k = ks+1 , ke
        if ( array(k) >= x ) then
          x = array(k)
          kmax = k
        end if
      end do
    end subroutine maxim
    !
    ! Normalized updraft mass flux function
    !
    subroutine zunc(kbeg,zu,kb,r,z,kz,lp,knum,ktop)
      implicit none
      integer(ik4) , intent(in) :: kz , kbeg , kb , lp , knum
      real(rkx) , dimension(kz,knum) , intent(out) :: zu
      real(rkx) , dimension(kz) , intent(in) :: z
      real(rkx) , intent(in) :: r
      integer(ik4) , dimension(knum) :: ktop
      integer(ik4) :: k , lpt
      real(rkx) :: dz
      !
      ! Below updraft air originating level (zu)
      !
      do k = 1 , kbeg
        zu(k,lp) = d_zero
      end do
      !
      ! Between zu and level of free convectioN
      !
      do k = kbeg , kb
        zu(k,lp) = d_one
      end do
      !
      ! Between level of free convection and cloud top
      !
      do k = kb+1 , ktop(lp)-1
        dz = d_half * (z(k+1)-z(k-1))
        zu(k,lp) = zu(k-1,lp) * (d_one + r*dz)
      end do
      !
      ! Cloud top
      !
      lpt = ktop(lp)
      dz = d_half * (z(lpt)-z(lpt-1))
      zu(lpt,lp) = zu(lpt-1,lp) * (d_one + r*dz)
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
    subroutine precip(cd,kb,kbcon,kz,r,hc,hes,t,qe, &
                      qes,pw,qc,qrc,z,lp,p,qkb,pcut,c0,zu,knum,ktop)
      implicit none
      integer(ik4) , intent(in) :: kbcon , kz , kb , lp , knum
      real(rkx) , dimension(kz) , intent(in) :: hc , hes , t , z , p
      real(rkx) , dimension(kz) , intent(in) :: qe , qes
      real(rkx) , dimension(kz) , intent(out) :: qc
      real(rkx) , dimension(kz,knum) , intent(in) :: cd , zu
      real(rkx) , dimension(kz,knum) , intent(out) :: pw , qrc
      real(rkx) , intent(in) :: c0 , pcut , qkb , r
      integer(ik4) , dimension(knum) :: ktop
      real(rkx) :: dh , dz , agamma , qrch
      integer :: k
      !
      ! Erase
      !
      do k = 1 , kz
        pw(k,lp) = d_zero
        qc(k) = d_zero
      end do
      !
      ! Initialize
      !
      do k = 1 , kb
        qrc(k,lp) = d_half*(qe(k)+qe(k+1))
      end do
      do k = kb , kbcon
        qrc(k,lp) = qkb
      end do
      do k = ktop(lp) , kz
        qrc(k,lp) = qes(k)
      end do
      !
      ! Boundary condition
      !
      qc(kbcon) = qkb
      pw(kbcon,lp) = d_zero
      !
      ! Start loop
      !
      do k = kbcon+1 , ktop(lp)
        if ( k == ktop(lp) ) then
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
        qrc(k,lp) = (qc(k)-qrch)/(d_one+c0*dz)+qrch
        !
        ! 3. Precip
        !
        pw(k,lp) = c0*dz*(qrc(k,lp)-qrch)*zu(k,lp)
        !
        ! If p<pcut, all water is rained out!
        !
        if ( k == ktop(lp) ) then
          qrc(k,lp) = (qc(k)-qes(k))/(d_one+c0*dz)+qes(k)
          pw(k,lp) = (qrc(k,lp)-qes(k))*c0*dz*zu(k,lp)
          if ( p(ktop(lp)) < pcut ) then
            qrc(k,lp) = qes(k)
            pw(k,lp) = (qc(k)-qrc(k,lp))*zu(k,lp)
          end if
        end if
        !
        ! 4.+5.
        !
        qc(k) = qrc(k,lp)
      end do
    end subroutine precip

  end subroutine shallow

end module mod_cu_shallow

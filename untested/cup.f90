!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine cup(qcrit,t,q,z1,tn,qo,po,pre,p,outtem,outq,dtime,psur,&
                   & vsp,istart,iend,kdet,jslc)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      use mod_regcm_param
      use mod_param2
      use mod_pmoist
      use mod_rad
      use mod_trachem
      use mod_constants , only : gti , rgti , cpd , tmelt , wlhv ,      &
                               & rwat , rcpd , wlhvocp
      implicit none
!
! Dummy arguments
!
      real(8) :: dtime
      integer :: iend , istart , jslc
      integer , dimension(ix) :: kdet
      real(8) , dimension(ix,kx) :: outq , outtem , p , po , q , qo ,   &
                                  & t , tn , vsp
      real(8) , dimension(ix) :: pre , psur , qcrit , z1
      intent (in) dtime , jslc , kdet , p , po , psur , qcrit , t , tn ,&
                & z1
      intent (inout) outq , outtem , pre , q , qo
!
! Local variables
!
      real(8) , dimension(ix) :: aa0 , aa1 , bu , buo , edt , edto ,    &
                               & edtx , hcd , hcdo , hkb , hkbo , pwav ,&
                               & pwavo , pwev , pwevo , qcd , qcdo ,    &
                               & qck , qcko , qkb , qkbo , vshear ,     &
                               & xaa0 , xhcd , xhkb , xmb , xpwav ,     &
                               & xpwev , xqcd , xqck , xqkb
      real(8) :: adw , akclth , alsixt , aup , c0 , detdo ,             &
               & detdoq , dg , dh , dhh , dp , dq , dt , dv1 , dv1q ,   &
               & dv2 , dv2q , dv3 , dv3q , dz , dz1 , dz2 , dzo , e ,   &
               & eo , f , agamma , gamma0 , gamma1 , gamma2 , gammo ,   &
               & gammo0 , mbdt , outtes , pbcdif , qrch , qrcho ,       &
               & tcrit , tfinv , tvbar , tvbaro , xk
      real(8) , dimension(2) :: ae , be , ht
      real(8) , dimension(ix,kx) :: dby , dbyo , dellah , dellaq ,      &
                                  & dellat , dkk , he , heo , hes ,     &
                                  & heso , pw , pwd , pwdo , pwo , qc , &
                                  & qco , qes , qeso , qrcd , qrcdo ,   &
                                  & tv , tvo , xdby , xhe , xhes , xpw ,&
                                  & xpwd , xq , xqc , xqes , xqrcd ,    &
                                  & xt , xtv , xz , z , zo
      integer :: i , iph , ipho , k , kbcono , kclth , kk , lpt
      integer , dimension(ix) :: jmin , k22 , kb , kbcon , kds , ktop
!
      tcrit = 50.
 
      tfinv = 1./tmelt
      alsixt = dlog(610.71D0)
      ht(1) = wlhvocp
      ht(2) = 2.834E6*rcpd
      be(1) = .622*ht(1)*3.50
      ae(1) = be(1)*tfinv + alsixt
      be(2) = .622*ht(2)*3.50
      ae(2) = be(2)*tfinv + alsixt
      mbdt = dtime*5.E-03
      c0 = .002
      f = -1.
      xk = -1.
!
!---  environmental conditions, first heights
!
      do k = 1 , kx
        do i = istart , iend
          dkk(i,k) = 1.
          iph = 1
          ipho = 1
          if ( t(i,k).le.tcrit ) iph = 2
          if ( tn(i,k).le.tcrit ) ipho = 2
          e = dexp(ae(iph)-be(iph)/t(i,k))
          eo = dexp(ae(ipho)-be(ipho)/tn(i,k))
          qes(i,k) = .622*e/(100.*p(i,k)-(1.-.622)*e)
          qeso(i,k) = .622*eo/(100.*po(i,k)-(1.-.622)*eo)
          if ( qes(i,k).le.1.E-08 ) qes(i,k) = 1.E-08
          if ( q(i,k).gt.qes(i,k) ) q(i,k) = qes(i,k)
          if ( qeso(i,k).le.1.E-08 ) qeso(i,k) = 1.E-08
          if ( qo(i,k).gt.qeso(i,k) ) qo(i,k) = qeso(i,k)
          tv(i,k) = t(i,k) + .608*q(i,k)*t(i,k)
          tvo(i,k) = tn(i,k) + .608*qo(i,k)*tn(i,k)
        end do
      end do
      do i = 1 , ix
        hkb(i) = 0. ! EES
        qkb(i) = 0.
        hkbo(i) = 0.
        qkbo(i) = 0.
        xhkb(i) = 0.
        xqkb(i) = 0.
        edt(i) = 0.
        edto(i) = 0.
        edtx(i) = 0.
      end do
      do i = istart , iend
!       hkb(i)=0.
!       qkb(i)=0.
!       hkbo(i)=0.
!       qkbo(i)=0.
!       xhkb(i)=0.
!       xqkb(i)=0.
        aa1(i) = 0.
        aa0(i) = 0.
        if ( qcrit(i).le.0. ) aa0(i) = -1.
        xaa0(i) = 0.
        xpwav(i) = 0.
        xpwev(i) = 0.
        pwav(i) = 0.
        pwev(i) = 0.
        pwavo(i) = 0.
        pwevo(i) = 0.
        k22(i) = 1
        ktop(i) = 1
        kbcon(i) = 1
        kb(i) = 1
        kds(i) = 1
        jmin(i) = 1
!       edt(i)=0.
!       edto(i)=0.
!       edtx(i)=0.
        xmb(i) = 0.
        vshear(i) = 0.
        z(i,1) = z1(i) - (dlog(p(i,1))-dlog(psur(i)))*287.*tv(i,1)*rgti
        zo(i,1) = z1(i) - (dlog(po(i,1))-dlog(psur(i)))*287.*tvo(i,1)   &
                & *rgti
      end do
      do k = 2 , kx
        do i = istart , iend
          tvbar = .5*tv(i,k) + .5*tv(i,k-1)
          z(i,k) = z(i,k-1) - (dlog(p(i,k))-dlog(p(i,k-1)))             &
                 & *287.*tvbar*rgti
          tvbaro = .5*tvo(i,k) + .5*tvo(i,k-1)
          zo(i,k) = zo(i,k-1) - (dlog(po(i,k))-dlog(po(i,k-1)))         &
                  & *287.*tvbaro*rgti
        end do
      end do
!
!---  moist static energy
!
      do k = 1 , kx
        do i = istart , iend
          cldlwc(i,k) = 0.
          cldfra(i,k) = 0.
          pw(i,k) = 0.
          xpw(i,k) = 0.
          pwo(i,k) = 0.
          qc(i,k) = 0.
          xqc(i,k) = 0.
          qco(i,k) = 0.
          pwd(i,k) = 0.
          pwdo(i,k) = 0.
          xpwd(i,k) = 0.
          dellah(i,k) = 0.
          dellaq(i,k) = 0.
          dellat(i,k) = 0.
          he(i,k) = gti*z(i,k) + cpd*t(i,k) + 2.5E06*q(i,k)
          hes(i,k) = gti*z(i,k) + cpd*t(i,k) + 2.5E06*qes(i,k)
          if ( he(i,k).ge.hes(i,k) ) he(i,k) = hes(i,k)
          heo(i,k) = gti*zo(i,k) + cpd*tn(i,k) + 2.5E06*qo(i,k)
          heso(i,k) = gti*zo(i,k) + cpd*tn(i,k) + 2.5E06*qeso(i,k)
          if ( heo(i,k).ge.heso(i,k) ) heo(i,k) = heso(i,k)
          xt(i,k) = t(i,k)
          xq(i,k) = q(i,k)
          xhe(i,k) = he(i,k)
          if ( k.ne.kx ) qrcd(i,k) = .5*(qes(i,k)+qes(i,k+1))
          if ( k.ne.kx ) qrcdo(i,k) = .5*(qeso(i,k)+qeso(i,k+1))
        end do
      end do
!
!------- determine level with highest moist static energy content.
!
      call maximi(he,ix,kx,1,kbmax2d(i,jslc),k22,istart,iend)
      do i = istart , iend
        if ( aa0(i).ge.0. ) then
          if ( k22(i).ge.kbmax2d(i,jslc) ) then
            aa0(i) = -1.
            go to 100
          end if
          hkb(i) = he(i,k22(i))
          qkb(i) = q(i,k22(i))
          hkbo(i) = heo(i,k22(i))
          qkbo(i) = qo(i,k22(i))
          qck(i) = qkb(i)
          qcko(i) = qkbo(i)
        end if
        100 continue
      end do
!
!---  decide for convective cloud base
!
      do i = istart , iend
        if ( aa0(i).ge.0. ) then
          do k = 1 , kdet(i)
            kk = kdet(i) - k + 1
!           dkk(i,kk)=.75*dkk(i,kk+1)
            dkk(i,k) = 1. - dble(kk)/dble(kdet(i))
          end do
 120      continue
          kb(i) = k22(i)
!----------------------------------
          kbcon(i) = kb(i)
 140      continue
          dh = .5*hes(i,kbcon(i)) + .5*hes(i,kbcon(i)+1)
          if ( hkb(i).lt.dh ) then
            kbcon(i) = kbcon(i) + 1
            if ( kbcon(i).gt.kbmax2d(i,jslc) ) then
              aa0(i) = -1.
              go to 200
            end if
            go to 140
          else
!
!---        after large-scale forcing is applied, possible lid should be
!---        removed!!!
!
            kbcono = kb(i)
!ictp
 150        continue
            if ( kbcono.gt.kbmax2d(i,jslc) ) then
              aa0(i) = -1.
              go to 200
            end if
!ictp_
            dh = .5*heso(i,kbcono) + .5*heso(i,kbcono+1)
            if ( hkbo(i).lt.dh ) then
              kbcono = kbcono + 1
              go to 150
            else
              pbcdif = -p(i,kbcono) + p(i,kb(i))
!----------------------------- below was commented out
!as           uncommenting the following lines for experiment 2/5/95
              if ( pbcdif.gt.pbcmax2d(i,jslc) ) then
                                            !this is where typo was (pbdcdif)
                k22(i) = k22(i) + 1
                if ( k22(i).ge.kbmax2d(i,jslc) ) then
                  aa0(i) = -1.
                  go to 200
                end if
                hkb(i) = he(i,k22(i))
                qkb(i) = q(i,k22(i))
                hkbo(i) = heo(i,k22(i))
                qkbo(i) = qo(i,k22(i))
                qck(i) = qkb(i)
                qcko(i) = qkbo(i)
                go to 120
              end if
            end if
          end if
!as
        end if
 200    continue
      end do
 
!
!---  downdraft originating level
!
      call minimi(he,ix,kx,kb,kx,jmin,istart,iend)
      call maximi(vsp,ix,kx,1,kx,kds,istart,iend)
!
!**************************** static control
!
!
!---  determine cloud top
!
      do i = istart , iend
        if ( aa0(i).ge.0 ) then
          if ( jmin(i).le.3 ) then
            aa0(i) = -1.
            go to 300
          end if
          if ( kds(i).ge.kx ) kds(i) = kx - 1
          if ( kds(i).le.kbcon(i) ) kds(i) = kbcon(i)
          dby(i,kx) = hkb(i) - hes(i,kx)
          dbyo(i,kx) = hkbo(i) - heso(i,kx)
        end if
 300    continue
      end do
      do k = 1 , kx - 1
        do i = istart , iend
          if ( aa0(i).ne.-1. ) then
            dby(i,k) = hkb(i) - .5*(hes(i,k)+hes(i,k+1))
            dbyo(i,k) = hkbo(i) - .5*(heso(i,k)+heso(i,k+1))
          end if
        end do
      end do
      do i = istart , iend
        if ( aa0(i).ne.-1. ) then
          do k = 2 , kx - kbcon(i) - 1
            kk = kx - k + 1
            if ( dby(i,kk).ge.0. ) then
              ktop(i) = kk + 1
              go to 320
            end if
          end do
          aa0(i) = -1.
          go to 400
 320      continue
          if ( ktop(i).gt.kx ) ktop(i) = kx
          if ( p(i,kbcon(i))-p(i,ktop(i)).lt.mincld2d(i,jslc) ) aa0(i)  &
             & = -1.
        end if
 400    continue
      end do
!
 
!------- moisture and cloud work functions
!
      do k = 2 , kx - 1
        do i = istart , iend
          if ( aa0(i).ne.-1. ) then
            if ( k.gt.kbcon(i) ) then
              if ( k.lt.ktop(i) ) then
                dz = -.5*z(i,k-1) + .5*z(i,k+1)
                dz1 = z(i,k) - z(i,k-1)
                agamma = (wlhvocp)*(wlhv/(rwat*(t(i,k)**2)))*qes(i,k)
                gamma0 = (wlhvocp)*(wlhv/(rwat*(t(i,k-1)**2)))*         &
                        & qes(i,k-1)
                qrch = qes(i,k) + (1./wlhv)*(agamma/(1.+agamma))*       &
                        & dby(i,k)
                qc(i,k) = (qck(i)-qrch)/(1.+c0*dz) + qrch
                pw(i,k) = c0*dz*(qc(i,k)-qrch)
                qck(i) = qc(i,k)
                pwav(i) = pwav(i) + pw(i,k)
                dz1 = z(i,k) - z(i,k-1)
                aa0(i) = aa0(i)                                         &
                       & + dz1*(gti/(cpd*(.5*(t(i,k)+t(i,k-1)))))       &
                       & *dby(i,k-1)/(1.+.5*agamma+.5*gamma0)
                dzo = -.5*zo(i,k-1) + .5*zo(i,k+1)
                dz2 = zo(i,k) - zo(i,k-1)
                gammo = (wlhvocp)*(wlhv/(rwat*(tn(i,k)**2)))*qeso(i,k)
                gammo0 = (wlhvocp)*(wlhv/(rwat*(tn(i,k-1)**2)))*        &
                       & qeso(i,k-1)
                qrcho = qeso(i,k) + (1./wlhv)*(gammo/(1.+gammo))*       &
                       & dbyo(i,k)
                qco(i,k) = (qcko(i)-qrcho)/(1.+c0*dzo) + qrcho
                pwo(i,k) = c0*dzo*(qco(i,k)-qrcho)
                qcko(i) = qco(i,k)
                pwavo(i) = pwavo(i) + pwo(i,k)
                aa1(i) = aa1(i)                                         &
                       & + dz2*(gti/(cpd*(.5*(tn(i,k)+tn(i,k-1)))))     &
                       & *dbyo(i,k-1)/(1.+.5*gammo+.5*gammo0)
              end if
            end if
          end if
        end do
      end do
!
!
      do i = istart , iend
        if ( aa0(i).ne.-1. ) then
          k = ktop(i)
          dz = -.5*z(i,k-1) + .5*z(i,k)
          agamma = (wlhvocp)*(wlhv/(rwat*(t(i,k)**2)))*qes(i,k)
          qrch = qes(i,k) + (1./wlhv)*(agamma/(1.+agamma))*dby(i,k)
          qc(i,k) = qes(i,k)
          pw(i,k) = (qrch-qes(i,k))
          pwav(i) = pwav(i) + pw(i,k)
!
          dz = -.5*zo(i,k-1) + .5*zo(i,k)
          agamma = (wlhvocp)*(wlhv/(rwat*(tn(i,k)**2)))*qeso(i,k)
          qrcho = qeso(i,k) + (1./wlhv)*(agamma/(1.+agamma))*dbyo(i,k)
          qco(i,k) = qeso(i,k)
          pwo(i,k) = (qrcho-qeso(i,k))
          pwavo(i) = pwavo(i) + pwo(i,k)
        end if
      end do
!
!------- downdraft calculations
!
!
!---  determine downdraft strength in terms of windshear
!
      do kk = 1 , kx/2
        do i = istart , iend
          if ( aa0(i).ne.-1. ) vshear(i) = vshear(i)                    &
             & + dabs((vsp(i,kk+1)-vsp(i,kk))/(z(i,kk+1)-z(i,kk)))
        end do
      end do
      do i = istart , iend
        if ( aa0(i).ne.-1. ) then
          vshear(i) = 1.E3*vshear(i)/dble(kx/2)
          edt(i) = 1. - (1.591-.639*vshear(i)+.0953*(vshear(i)**2)      &
                 & -.00496*(vshear(i)**3))
 
          if ( edt(i).gt.shrmax2d(i,jslc) ) edt(i) = shrmax2d(i,jslc)
          if ( edt(i).lt.shrmin2d(i,jslc) ) edt(i) = shrmin2d(i,jslc)
 
          edto(i) = edt(i)
          edtx(i) = edt(i)
          qrcd(i,kx) = qes(i,kx)
          hcd(i) = .5*(he(i,jmin(i))+he(i,jmin(i)+1))
          qcd(i) = .5*(q(i,jmin(i))+q(i,jmin(i)+1))
          qrcdo(i,kx) = qeso(i,kx)
          hcdo(i) = heso(i,kx)
          hcdo(i) = .5*(heo(i,jmin(i))+heo(i,jmin(i)+1))
          qcdo(i) = .5*(qo(i,jmin(i))+qo(i,jmin(i)+1))
          bu(i) = 0.
          buo(i) = 0.
        end if
      end do
      do k = 1 , kx - 1
        do i = istart , iend
          if ( aa0(i).ne.-1. ) then
            if ( k.lt.jmin(i) ) then
              kk = jmin(i) - k
              dz = -(z(i,kk)-z(i,kk+2))*.5
              bu(i) = bu(i) + dz*(hcd(i)-.5*(hes(i,kk)+hes(i,kk+1)))
              dq = (qes(i,kk)+qes(i,kk+1))*.5
              dt = (t(i,kk)+t(i,kk+1))*.5
              agamma = (wlhvocp)*(wlhv/(rwat*(dt**2)))*dq
              dh = hcd(i) - .5*(hes(i,kk)+hes(i,kk+1))
              qrcd(i,kk) = (dq+(1./wlhv)*(agamma/(1.+agamma))*dh)
              pwd(i,kk) = dkk(i,kk)*(qcd(i)-qrcd(i,kk))
              qcd(i) = qrcd(i,kk)
              pwev(i) = pwev(i) + pwd(i,kk)
!
              dz = -(zo(i,kk)-zo(i,kk+2))*.5
              buo(i) = buo(i) + dz*(hcdo(i)-.5*(heso(i,kk)+heso(i,kk+1))&
                     & )
              dq = (qeso(i,kk)+qeso(i,kk+1))*.5
              dt = (tn(i,kk)+tn(i,kk+1))*.5
              agamma = (wlhvocp)*(wlhv/(rwat*(dt**2)))*dq
              dh = hcdo(i) - .5*(heso(i,kk)+heso(i,kk+1))
              qrcdo(i,kk) = (dq+(1./wlhv)*(agamma/(1.+agamma))*dh)
              pwdo(i,kk) = dkk(i,kk)*(qcdo(i)-qrcdo(i,kk))
              qcdo(i) = qrcdo(i,kk)
              pwevo(i) = pwevo(i) + pwdo(i,kk)
            end if
          end if
        end do
      end do
!
      do i = istart , iend
        if ( aa0(i).ne.-1. ) then
          if ( bu(i).ge.0 .or. buo(i).ge.0 .or. pwev(i).ge.0 .or.       &
             & pwevo(i).ge.0. ) aa0(i) = -1.
          edt(i) = -edt(i)*pwav(i)/pwev(i)
          if ( edt(i).gt.edtmax2d(i,jslc) ) edt(i) = edtmax2d(i,jslc)
          if ( edt(i).lt.edtmin2d(i,jslc) ) edt(i) = edtmin2d(i,jslc)
          edto(i) = -edto(i)*pwavo(i)/pwevo(i)
          if ( edto(i).gt.edtmaxo2d(i,jslc) ) edto(i)                   &
             & = edtmaxo2d(i,jslc)
          if ( edto(i).lt.edtmino2d(i,jslc) ) edto(i)                   &
             & = edtmino2d(i,jslc)
        end if
      end do
!
!---  what would the change be?
!
      do i = istart , iend
        if ( aa0(i).ne.-1. ) then
          k = 1
          dz = .5*(z(i,2)-z(i,1))
          dp = 50.*(psur(i)-p(i,2))
          dellah(i,1) = edt(i)                                          &
                      & *(dkk(i,1)*hcd(i)-dkk(i,1)*.5*(he(i,1)+he(i,2)))&
                      & *gti/dp
          dellaq(i,1) = edt(i)                                          &
                      & *(dkk(i,1)*qrcd(i,1)-dkk(i,1)*.5*(q(i,1)+q(i,2))&
                      & )*gti/dp
          xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
          xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
          dellat(i,k) = rcpd*(dellah(i,k)-2.5E06*dellaq(i,k))
          xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-2.5E06*dellaq(i,k))+t(i,k)
          if ( xq(i,k).le.0. ) xq(i,k) = 1.E-08
        end if
      end do
!
      do k = 1 , kx - 1
        do i = istart , iend
          if ( aa0(i).ne.-1. ) then
            if ( k.ne.1 .and. k.lt.ktop(i) ) then
              dv1 = .5*(he(i,k)+he(i,k+1))
              dv2 = he(i,k)
              dv3 = .5*(he(i,k)+he(i,k-1))
              dv1q = .5*(q(i,k)+q(i,k+1))
              dv2q = q(i,k)
              dv3q = .5*(q(i,k)+q(i,k-1))
!
!---          specifiy detrainment of downdraft, has to be consistent
!---          with zd calculations in soundd.
!
              detdo = (1.-dkk(i,k))*(hcd(i)-dv2)
              detdoq = (1.-dkk(i,k))*(qrcd(i,k)-dv2q)
              dz = .5*(z(i,k+1)-z(i,k-1))
!
!---          changed due to subsidence and entrainment
!
              aup = 1.
              if ( k.le.k22(i) ) aup = 0.
              adw = 1.
              if ( k.gt.jmin(i) ) adw = 0.
              dp = +50.*(p(i,k-1)-p(i,k+1))
              dellah(i,k) = ((aup-adw*edt(i))*(dv1-dv2)+(aup-adw*edt(i))&
                          & *(dv2-dv3))*gti/dp + adw*edt(i)*detdo*gti/dp
              dellaq(i,k) = ((aup-adw*edt(i))*(dv1q-dv2q)+(aup-adw*edt(i&
                          & ))*(dv2q-dv3q))*gti/dp + adw*edt(i)         &
                          & *detdoq*gti/dp
              xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
              xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
              dellat(i,k) = rcpd*(dellah(i,k)-2.5E06*dellaq(i,k))
              xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-2.5E06*dellaq(i,k))    &
                      & + t(i,k)
              if ( xq(i,k).le.0. ) xq(i,k) = 1.E-08
            end if
          end if
        end do
      end do
!
!------- cloud top
!
      do i = istart , iend
        if ( aa0(i).ne.-1. ) then
          lpt = ktop(i)
          dp = 100.*(p(i,lpt-1)-p(i,lpt))
          dv1 = .5*(he(i,lpt)+he(i,lpt-1))
          dellah(i,lpt) = (hkb(i)-dv1)*gti/dp
          dv1 = .5*(q(i,lpt)+q(i,lpt-1))
          dellaq(i,lpt) = (qes(i,lpt)-dv1)*gti/dp
          k = lpt
          xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
          xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
          dellat(i,k) = rcpd*(dellah(i,k)-2.5E06*dellaq(i,k))
          xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-2.5E06*dellaq(i,k))        &
                  & + t(i,k)
          if ( xq(i,k).le.0. ) xq(i,k) = 1.E-08
          xhkb(i) = dellah(i,kbcon(i))*mbdt + hkb(i)
          xqkb(i) = dellaq(i,kbcon(i))*mbdt + qkb(i)
          if ( xqkb(i).le.0. ) xqkb(i) = 1.E-08
        end if
      end do
!
!---  environmental conditions, first heights
!
      do k = 1 , kx
        do i = istart , iend
          if ( aa0(i).ne.-1. ) then
            iph = 1
            if ( xt(i,k).le.tcrit ) iph = 2
            e = dexp(ae(iph)-be(iph)/xt(i,k))
            xqes(i,k) = .622*e/(100.*p(i,k)-(1.-.622)*e)
            if ( xqes(i,k).le.1.E-08 ) xqes(i,k) = 1.E-08
            if ( xq(i,k).gt.xqes(i,k) ) xq(i,k) = xqes(i,k)
            xtv(i,k) = xt(i,k) + .608*xq(i,k)*xt(i,k)
          end if
        end do
      end do
!     bug fix
      do k = 1 , kx - 1
        do i = istart , iend
          if ( aa0(i).ne.-1 ) xqrcd(i,k) = .5*(xqes(i,k)+xqes(i,k+1))
        end do
      end do
!
      do i = istart , iend
        if ( aa0(i).ne.-1. ) xz(i,1) = z1(i)                            &
                                     & - (dlog(p(i,1))-dlog(psur(i)))   &
                                     & *287.*xtv(i,1)*rgti
      end do
      do k = 2 , kx
        do i = istart , iend
          if ( aa0(i).ne.-1. ) then
            tvbar = .5*xtv(i,k) + .5*xtv(i,k-1)
            xz(i,k) = xz(i,k-1) - (dlog(p(i,k))-dlog(p(i,k-1)))         &
                    & *287.*tvbar*rgti
          end if
        end do
      end do
!
!---  moist static energy
!
      do k = 1 , kx
        do i = istart , iend
          if ( aa0(i).ne.-1. ) then
            xhes(i,k) = gti*xz(i,k) + cpd*xt(i,k) + 2.5E06*xqes(i,k)
            if ( xhe(i,k).ge.xhes(i,k) ) xhe(i,k) = xhes(i,k)
          end if
        end do
      end do
!
!
!**************************** static control
!
      do i = istart , iend
        if ( aa0(i).ne.-1. ) then
          xqck(i) = xqkb(i)
          xdby(i,kx) = xhkb(i) - xhes(i,kx)
        end if
      end do
!
!------- moisture and cloud work functions
!
      do k = 1 , kx - 1
        do i = istart , iend
          if ( aa0(i).ge.0. ) then
            xdby(i,k) = xhkb(i) - .5*(xhes(i,k)+xhes(i,k+1))
            if ( k.gt.kbcon(i) .and. k.lt.ktop(i) ) then
              dz = -.5*xz(i,k-1) + .5*xz(i,k+1)
              dz1 = xz(i,k) - xz(i,k-1)
              agamma = (wlhvocp)*(wlhv/(rwat*(xt(i,k)**2)))*xqes(i,k)
              gamma0 = (wlhvocp)*(wlhv/(rwat*(xt(i,k-1)**2)))*          &
                      & xqes(i,k-1)
              qrch = xqes(i,k) + (1./wlhv)*(agamma/(1.+agamma))*        &
                      & xdby(i,k)
              xqc(i,k) = (xqck(i)-qrch)/(1.+c0*dz) + qrch
              xpw(i,k) = c0*dz*(xqc(i,k)-qrch)
              xqck(i) = xqc(i,k)
              xpwav(i) = xpwav(i) + xpw(i,k)
              xaa0(i) = xaa0(i)                                         &
                      & + dz1*(gti/(cpd*(.5*(xt(i,k)+xt(i,k-1)))))      &
                      & *xdby(i,k-1)/(1.+.5*agamma+.5*gamma0)
            end if
          end if
        end do
      end do
      do i = istart , iend
        if ( aa0(i).ge.0. ) then
          k = ktop(i)
          dz = -.5*xz(i,k-1) + .5*xz(i,k)
          agamma = (wlhvocp)*(wlhv/(rwat*(xt(i,k)**2)))*xqes(i,k)
          qrch = xqes(i,k) + (1./wlhv)*(agamma/(1.+agamma))*xdby(i,k)
          xqc(i,k) = xqes(i,k)
          xpw(i,k) = (qrch-xqes(i,k))
          xpwav(i) = xpwav(i) + xpw(i,k)
          xqrcd(i,kx) = xqes(i,kx)
          xhcd(i) = .5*(xhe(i,jmin(i))+xhe(i,jmin(i)+1))
          xqcd(i) = .5*(xq(i,jmin(i))+xq(i,jmin(i)+1))
          xpwev(i) = 0.
          bu(i) = 0.
        end if
      end do
!
!------- downdraft calculations
!
!
!---  downdraft moisture properties
!
      do k = 1 , kx - 1
        do i = istart , iend
          if ( aa0(i).ge.0. ) then
            if ( k.lt.jmin(i) ) then
              kk = jmin(i) - k
              dz = -(xz(i,kk)-xz(i,kk+2))*.5
              bu(i) = bu(i) + dz*(xhcd(i)-.5*(xhes(i,kk)+xhes(i,kk+1)))
              dq = (xqes(i,kk)+xqes(i,kk+1))*.5
              dt = (xt(i,kk)+xt(i,kk+1))*.5
              agamma = (wlhvocp)*(wlhv/(rwat*(dt**2)))*dq
              dh = xhcd(i) - .5*(xhes(i,kk)+xhes(i,kk+1))
              xqrcd(i,kk) = (dq+(1./wlhv)*(agamma/(1.+agamma))*dh)
              xpwd(i,kk) = dkk(i,kk)*(xqcd(i)-xqrcd(i,kk))
              xqcd(i) = xqrcd(i,kk)
              xpwev(i) = xpwev(i) + xpwd(i,kk)
            end if
          end if
        end do
      end do
      do i = istart , iend
        if ( aa0(i).ge.0. ) then
          if ( bu(i).ge.0. ) then
            aa0(i) = -1.
            go to 500
          end if
          if ( xpwev(i).ne.0. ) edtx(i) = -edtx(i)*xpwav(i)/xpwev(i)
          if ( edtx(i).gt.edtmaxx2d(i,jslc) ) edtx(i)                   &
             & = edtmaxx2d(i,jslc)
          if ( edtx(i).lt.edtminx2d(i,jslc) ) edtx(i)                   &
             & = edtminx2d(i,jslc)
        end if
 500    continue
      end do
!
!
!---  downdraft cloudwork functions
!
!
      do k = 1 , kx - 1
        do i = istart , iend
          if ( aa0(i).ge.0. ) then
            if ( k.lt.jmin(i) ) then
              kk = jmin(i) - k
!
!---          original
!
              gamma1 = (wlhvocp)*(wlhv/(rwat*(t(i,kk)**2)))*qes(i,kk)
              gamma2 = (wlhvocp)*(wlhv/(rwat*(t(i,kk+1)**2)))*          &
                     & qes(i,kk+1)
              dhh = hcd(i)
              dt = .5*(t(i,kk)+t(i,kk+1))
              dg = .5*(gamma1+gamma2)
              dh = .5*(hes(i,kk)+hes(i,kk+1))
              dz = (z(i,kk)-z(i,kk+1))*dkk(i,kk)
              aa0(i) = aa0(i) + edt(i)*dz*(gti/(cpd*dt))*((dhh-dh)/     &
                     & (1.+dg))
!
!---          modified by larger scale
!
              gamma1 = (wlhvocp)*(wlhv/(rwat*(tn(i,kk)**2)))*qeso(i,kk)
              gamma2 = (wlhvocp)*(wlhv/(rwat*(tn(i,kk+1)**2)))*         &
                      & qeso(i,kk+1)
              dhh = hcdo(i)
              dt = .5*(tn(i,kk)+tn(i,kk+1))
              dg = .5*(gamma1+gamma2)
              dh = .5*(heso(i,kk)+heso(i,kk+1))
              dz = (zo(i,kk)-zo(i,kk+1))*dkk(i,kk)
              aa1(i) = aa1(i) + edto(i)*dz*(gti/(cpd*dt))               &
                     & *((dhh-dh)/(1.+dg))
!
!---          modified by cloud
!
              gamma1 = (wlhvocp)*(wlhv/(rwat*(xt(i,kk)**2)))*xqes(i,kk)
              gamma2 = (wlhvocp)*(wlhv/(rwat*(xt(i,kk+1)**2)))*         &
                     & xqes(i,kk+1)
              dhh = xhcd(i)
              dt = .5*(xt(i,kk)+xt(i,kk+1))
              dg = .5*(gamma1+gamma2)
              dh = .5*(xhes(i,kk)+xhes(i,kk+1))
              dz = (xz(i,kk)-xz(i,kk+1))*dkk(i,kk)
              xaa0(i) = xaa0(i) + edtx(i)*dz*(gti/(cpd*dt))             &
                      & *((dhh-dh)/(1.+dg))
            end if
          end if
        end do
      end do
!
!---  large scale forcing
!
      do i = istart , iend
        if ( aa0(i).ge.0. ) then
          if ( igcc.eq.1 ) then
            f = (aa1(i)-aa0(i))/dtime  ! Arakawa-Schubert closure
          else if ( igcc.eq.2 ) then
            f = aa0(i)/dtauc2d(i,jslc)   ! Fritsch-Chappell closure
          else
          end if
          xk = (xaa0(i)-aa0(i))/mbdt
          xmb(i) = -f/xk
          if ( f.le.0 .or. xk.ge.0. ) xmb(i) = 0.
        end if
      end do
!chem2
      mflx(i,1) = xmb(i)
      mflx(i,2) = xmb(i)*edt(i)
!chem2_
!
!---  feedback
!
      do k = 1 , kx
        do i = istart , iend
          if ( aa0(i).ge.0. ) then
            if ( k.le.ktop(i) ) then
              outtes = dellat(i,k)*xmb(i)*86400.
              if ( (outtes.gt.htmax2d(i,jslc)) .or.                     &
                 & (outtes.lt.htmin2d(i,jslc)) ) then
                xmb(i) = 0.
                aa0(i) = -1.
              else
                outtem(i,k) = outtem(i,k) + dellat(i,k)*xmb(i)
                outq(i,k) = outq(i,k) + dellaq(i,k)*xmb(i)
                pre(i) = pre(i) + (pw(i,k)+edt(i)*pwd(i,k))*xmb(i)
              end if
            end if
          end if
        end do
      end do
!
!     calculate cloud fraction and water content
!
      do i = istart , iend
!chem2
        icumtop(i,jslc) = 0
        icumbot(i,jslc) = 0
        icumdwd(i,jslc) = 0
!chem2_
        if ( aa0(i).ge.0. ) then
 
          if ( ktop(i).gt.1 .and. kbcon(i).gt.1 ) then
            kclth = ktop(i) - kbcon(i) + 1
            akclth = 1./dble(kclth)
            do k = kbcon(i) , ktop(i)
              kk = kx - k + 1
              cldlwc(i,kk) = cllwcv
              cldfra(i,kk) = 1. - (1.-clfrcv)**akclth
            end do
!chem2
!chem2      define convection  base and top for tracers
            if ( ichem.eq.1 ) then
              if ( ktop(i).gt.1 .and. k22(i).ge.1 ) then
                icumtop(i,jslc) = kxp1 - ktop(i)
                icumbot(i,jslc) = kxp1 - k22(i)
                icumdwd(i,jslc) = kxp1 - jmin(i)
              end if
            end if
!chem2_
          end if
        end if
 
      end do
      end subroutine cup

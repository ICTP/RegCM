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

module mod_cu_grell

  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_main
  use mod_bats
  use mod_pmoist
  use mod_rad
  use mod_trachem
  use mod_date
  use mod_service 
 
contains

  subroutine cuparan(tten,qten,j)

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    implicit none
    !
    ! Dummy arguments
    !
    integer :: j
    real(8) , dimension(iy,kz) :: qten , tten
    intent (inout) qten , tten
    !
    ! Local variables
    !
    real(8) :: aprdiv , calc , dtime , pkdcut , pkk , prainx , us , vs
    integer :: i , iconj , icut , iend , istart , k , kk
    integer , dimension(iy) :: kdet
    real(8) , dimension(iy,kz) :: outq , outt , p , po , q , qo , t , &
         & tn , vsp
    real(8) , dimension(iy) :: pret , psur , qcrit , ter11
    !
    character (len=50) :: subroutine_name='cuparan'
    integer :: idindx=0
    !
    call time_begin(subroutine_name,idindx)
    !     zero out radiative clouds
    !
    cldlwc = 0.0D0
    cldfra = 0.0D0

    icut = 0
    dtime = dt
    pkdcut = 75.0D0
    istart = 2 + icut
    iend = iym2 - icut
    !
    !---  prepare input, erase output
    !
    do i = istart , iend
       kdet(i) = 2
       qcrit(i) = 0.0D0
       pret(i) = 0.0D0
    end do

    do k = 1 , kz
       do i = 2 + icut , iym2 - icut
          kk = kz - k + 1

          us = gwnd%usk(i,kk)
          vs = gwnd%vsk(i,kk)

          t(i,k) = atm2%t(i,kk,j)/sps2%ps(i,j)
          q(i,k) = atm2%qv(i,kk,j)/sps2%ps(i,j)
          if ( q(i,k) < 1.D-08 ) q(i,k) = 1.D-08
          tn(i,k) = t(i,k) + (tten(i,kk))/sps2%ps(i,j)*dtime
          qo(i,k) = q(i,k) + (qten(i,kk))/sps2%ps(i,j)*dtime
          p(i,k) = 10.0D0*sps2%ps(i,j)*a(kk) + 10.0D0*r8pt
          vsp(i,k) = dsqrt(us**2.0D0+vs**2.0D0)
          if ( qo(i,k) < 1.D-08 ) qo(i,k) = 1.D-08
          !
          po(i,k) = p(i,k)
          psur(i) = 10.0D0*sps2%ps(i,j) + 10.0D0*r8pt
          outt(i,k) = 0.0D0
          pkk = psur(i) - po(i,k)
          if ( pkk <= pkdcut ) kdet(i) = kdet(i) + 1
          outq(i,k) = 0.0D0
          ter11(i) = mddom%ht(i,j)*rgti
          if ( ter11(i) <= 0.0D0 ) ter11(i) = 1.D-05
          qcrit(i) = qcrit(i) + qten(i,kk)
       end do
    end do
    !
    !---  call cumulus parameterization
    !
    call cup(qcrit,t,q,ter11,tn,qo,po,pret,p,outt,outq,dtime,psur,vsp,&
         & istart,iend,kdet,j)
    do k = 1 , kz
       do i = istart , iend
          if ( pret(i) > 0.0D0 ) then
             kk = kz - k + 1
             tten(i,kk) = sps2%ps(i,j)*outt(i,k) + tten(i,kk)
             qten(i,kk) = sps2%ps(i,j)*outq(i,k) + qten(i,kk)
          end if
       end do
    end do
    !
    !---  rain in cm.
    !
    calc = 0.5D0
    iconj = 0
    do i = istart , iend
       if ( pret(i) > 0.0D0 ) then
          sfsta%rainc(i,j) = sfsta%rainc(i,j) + pret(i)*calc*dt
          !         print *,'sfsta%rainc(',i,j,')=',sfsta%rainc(i,j)
          iconj = iconj + 1
          !.....................precipitation rate for bats (mm/s)
          aprdiv = dble(nbatst)
          if ( jyear == jyear0 .and. ktau == 0 ) aprdiv = 1.0D0
          prainx = pret(i)*calc*dt
          pptc(i,j) = pptc(i,j) + prainx/(dtmin*60.0D0)/aprdiv
          !.......................................................
       end if
    end do
    icon(j) = iconj
    !
    call time_end(subroutine_name,idindx)
  end subroutine cuparan
  !
  !    GRELL CUMULUS SCHEME
  !
  subroutine cup(qcrit,t,q,z1,tn,qo,po,pre,p,outtem,outq,dtime,psur,&
       & vsp,istart,iend,kdet,jslc)

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    implicit none
    !
    ! Dummy arguments
    !
    real(8) :: dtime
    integer :: iend , istart , jslc
    integer , dimension(iy) :: kdet
    real(8) , dimension(iy,kz) :: outq , outtem , p , po , q , qo ,   &
         & t , tn , vsp
    real(8) , dimension(iy) :: pre , psur , qcrit , z1
    intent (in) dtime , jslc , kdet , p , po , psur , qcrit , t , tn ,&
         & z1
    intent (inout) outq , outtem , pre , q , qo
    !
    ! Local variables
    !
    real(8) , dimension(iy) :: aa0 , aa1 , bu , buo , edt , edto ,    &
         & edtx , hcd , hcdo , hkb , hkbo , pwav ,&
         & pwavo , pwev , pwevo , qcd , qcdo ,    &
         & qck , qcko , qkb , qkbo , vshear ,     &
         & xaa0 , xhcd , xhkb , xmb , xpwav ,     &
         & xpwev , xqcd , xqck , xqkb
    real(8) :: adw , akclth , alsixt , aup , c0 , detdo ,             &
         & detdoq , dg , dh , dhh , dp_s , dq , xdt , dv1 , dv1q ,  &
         & dv2 , dv2q , dv3 , dv3q , dz , dz1 , dz2 , dzo , e ,   &
         & eo , f , agamma , gamma0 , gamma1 , gamma2 , gammo ,   &
         & gammo0 , mbdt , outtes , pbcdif , qrch , qrcho ,       &
         & tcrit , tvbar , tvbaro , xk
    real(8) , dimension(2) :: ae , be , ht
    real(8) , dimension(iy,kz) :: dby , dbyo , dellah , dellaq ,      &
         & dellat , dkk , he , heo , hes ,     &
         & heso , pw , pwd , pwdo , pwo , qc , &
         & qco , qes , qeso , qrcd , qrcdo ,   &
         & tv , tvo , xdby , xhe , xhes , xpw ,&
         & xpwd , xq , xqc , xqes , xqrcd ,    &
         & xt , xtv , xz , z , zo
    integer :: i , iph , ipho , k , kbcono , kclth , kk , lpt
    integer , dimension(iy) :: jmin , k22 , kb , kbcon , kds , ktop
    !
    character (len=50) :: subroutine_name='cup'
    integer :: idindx=0
    !
    call time_begin(subroutine_name,idindx)
    tcrit = 50.0D0

    alsixt = dlog(610.71D0)
    ht(1) = wlhvocp
    ht(2) = 2.834D6*rcpd
    be(1) = ep2*ht(1)*3.50D0
    ae(1) = be(1)*rtzero + alsixt
    be(2) = ep2*ht(2)*3.50D0
    ae(2) = be(2)*rtzero + alsixt
    mbdt = dtime*5.D-03
    c0 = 0.002D0
    f = -1.0D0
    xk = -1.0D0
    !
    !---  environmental conditions, first heights
    !
    do k = 1 , kz
       do i = istart , iend
          dkk(i,k) = 1.0D0
          iph = 1
          ipho = 1
          if ( t(i,k) <= tcrit ) iph = 2
          if ( tn(i,k) <= tcrit ) ipho = 2
          e = dexp(ae(iph)-be(iph)/t(i,k))
          eo = dexp(ae(ipho)-be(ipho)/tn(i,k))
          qes(i,k) = ep2*e/(100.0D0*p(i,k)-(1.0D0-ep2)*e)
          qeso(i,k) = ep2*eo/(100.*po(i,k)-(1.-ep2)*eo)
          if ( qes(i,k) <= 1.D-08 ) qes(i,k) = 1.D-08
          if ( q(i,k) > qes(i,k) ) q(i,k) = qes(i,k)
          if ( qeso(i,k) <= 1.D-08 ) qeso(i,k) = 1.D-08
          if ( qo(i,k) > qeso(i,k) ) qo(i,k) = qeso(i,k)
          tv(i,k) = t(i,k) + 0.608D0*q(i,k)*t(i,k)
          tvo(i,k) = tn(i,k) + 0.608D0*qo(i,k)*tn(i,k)
       end do
    end do
    do i = 1 , iy
       hkb(i) = 0.0D0 ! EES
       qkb(i) = 0.0D0
       hkbo(i) = 0.0D0
       qkbo(i) = 0.0D0
       xhkb(i) = 0.0D0
       xqkb(i) = 0.0D0
       edt(i) = 0.0D0
       edto(i) = 0.0D0
       edtx(i) = 0.0D0
    end do
    do i = istart , iend
       !       hkb(i)=0.0D0
       !       qkb(i)=0.0D0
       !       hkbo(i)=0.0D0
       !       qkbo(i)=0.0D0
       !       xhkb(i)=0.0D0
       !       xqkb(i)=0.0D0
       aa1(i) = 0.0D0
       aa0(i) = 0.0D0
       if ( icup==99 .or. icup==98 ) then
         if (cumcon%cuscheme(i,jslc)/=2 ) then
           aa0(i) = -1.0D0
         end if
       end if
       if ( qcrit(i) <= 0.0D0 ) then
         aa0(i) = -1.0D0
       end if
       xaa0(i) = 0.0D0
       xpwav(i) = 0.0D0
       xpwev(i) = 0.0D0
       pwav(i) = 0.0D0
       pwev(i) = 0.0D0
       pwavo(i) = 0.0D0
       pwevo(i) = 0.0D0
       k22(i) = 1
       ktop(i) = 1
       kbcon(i) = 1
       kb(i) = 1
       kds(i) = 1
       jmin(i) = 1
       !       edt(i)=0.0D0
       !       edto(i)=0.0D0
       !       edtx(i)=0.0D0
       xmb(i) = 0.0D0
       vshear(i) = 0.0D0
       z(i,1) = z1(i) - (dlog(p(i,1))-dlog(psur(i)))*rgas*tv(i,1)*rgti
       zo(i,1) = z1(i) - (dlog(po(i,1))-dlog(psur(i)))*rgas*tvo(i,1)   &
            & *rgti
    end do
    do k = 2 , kz
       do i = istart , iend
          tvbar = 0.5D0*tv(i,k) + 0.5D0*tv(i,k-1)
          z(i,k) = z(i,k-1) - (dlog(p(i,k))-dlog(p(i,k-1)))             &
               & *rgas*tvbar*rgti
          tvbaro = 0.5D0*tvo(i,k) + 0.5D0*tvo(i,k-1)
          zo(i,k) = zo(i,k-1) - (dlog(po(i,k))-dlog(po(i,k-1)))         &
               & *rgas*tvbaro*rgti
       end do
    end do
    !
    !---  moist static energy
    !
    do k = 1 , kz
       do i = istart , iend
          cldlwc(i,k) = 0.D0
          cldfra(i,k) = 0.D0
          pw(i,k) = 0.D0
          xpw(i,k) = 0.D0
          pwo(i,k) = 0.D0
          qc(i,k) = 0.D0
          xqc(i,k) = 0.D0
          qco(i,k) = 0.D0
          pwd(i,k) = 0.D0
          pwdo(i,k) = 0.D0
          xpwd(i,k) = 0.D0
          dellah(i,k) = 0.D0
          dellaq(i,k) = 0.D0
          dellat(i,k) = 0.D0
          he(i,k) = gti*z(i,k) + cpd*t(i,k) + wlhv*q(i,k)
          hes(i,k) = gti*z(i,k) + cpd*t(i,k) + wlhv*qes(i,k)
          if ( he(i,k) >= hes(i,k) ) he(i,k) = hes(i,k)
          heo(i,k) = gti*zo(i,k) + cpd*tn(i,k) + wlhv*qo(i,k)
          heso(i,k) = gti*zo(i,k) + cpd*tn(i,k) + wlhv*qeso(i,k)
          if ( heo(i,k) >= heso(i,k) ) heo(i,k) = heso(i,k)
          xt(i,k) = t(i,k)
          xq(i,k) = q(i,k)
          xhe(i,k) = he(i,k)
          if ( k /= kz ) qrcd(i,k) = 0.5D0*(qes(i,k)+qes(i,k+1))
          if ( k /= kz ) qrcdo(i,k) = 0.5D0*(qeso(i,k)+qeso(i,k+1))
       end do
    end do
    !
    !------- determine level with highest moist static energy content.
    !
    call maximi(he,iy,kz,1,kbmax2d(i,jslc),k22,istart,iend)
    do i = istart , iend
       if ( aa0(i) >= 0.0D0 ) then
          if ( k22(i) >= kbmax2d(i,jslc) ) then
             aa0(i) = -1.0D0
             go to 100
          end if
          hkb(i) = he(i,k22(i))
          qkb(i) = q(i,k22(i))
          hkbo(i) = heo(i,k22(i))
          qkbo(i) = qo(i,k22(i))
          qck(i) = qkb(i)
          qcko(i) = qkbo(i)
       end if
100    continue
    end do
    !
    !---  decide for convective cloud base
    !
    do i = istart , iend
       if ( aa0(i) >= 0.0D0 ) then
          do k = 1 , kdet(i)
             kk = kdet(i) - k + 1
             !           dkk(i,kk)=.75*dkk(i,kk+1)
             dkk(i,k) = 1.0D0 - dble(kk)/dble(kdet(i))
          end do
120       continue
          kb(i) = k22(i)
          !----------------------------------
          kbcon(i) = kb(i)
140       continue
          dh = 0.5D0*hes(i,kbcon(i)) + 0.5D0*hes(i,kbcon(i)+1)
          if ( hkb(i) < dh ) then
             kbcon(i) = kbcon(i) + 1
             if ( kbcon(i) > kbmax2d(i,jslc) ) then
                aa0(i) = -1.0D0
                go to 200
             end if
             go to 140
          else
             !
             !-after large-scale forcing is applied, possible lid should be
             !-removed!!!
             !
             kbcono = kb(i)
             !ictp
150          continue
             if ( kbcono > kbmax2d(i,jslc) ) then
                aa0(i) = -1.0D0
                go to 200
             end if
             !ictp_
             dh = 0.5D0*heso(i,kbcono) + 0.5D0*heso(i,kbcono+1)
             if ( hkbo(i) < dh ) then
                kbcono = kbcono + 1
                go to 150
             else
                pbcdif = -p(i,kbcono) + p(i,kb(i))
                !-below was commented out
                !as uncommenting the following lines for experiment 2/5/95
                if ( pbcdif > pbcmax2d(i,jslc) ) then
                   !this is where typo was (pbdcdif)
                   k22(i) = k22(i) + 1
                   if ( k22(i) >= kbmax2d(i,jslc) ) then
                      aa0(i) = -1.0D0
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
    call minimi(he,iy,kz,kb,kz,jmin,istart,iend)
    call maximi(vsp,iy,kz,1,kz,kds,istart,iend)
    !
    !**************************** static control
    !
    !
    !---  determine cloud top
    !
    do i = istart , iend
       if ( aa0(i) >= 0 ) then
          if ( jmin(i) <= 3 ) then
             aa0(i) = -1.0D0
             go to 300
          end if
          if ( kds(i) >= kz ) kds(i) = kz - 1
          if ( kds(i) <= kbcon(i) ) kds(i) = kbcon(i)
          dby(i,kz) = hkb(i) - hes(i,kz)
          dbyo(i,kz) = hkbo(i) - heso(i,kz)
       end if
300    continue
    end do
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(i) > -1.0D0 ) then
             dby(i,k) = hkb(i) - 0.5D0*(hes(i,k)+hes(i,k+1))
             dbyo(i,k) = hkbo(i) - 0.5D0*(heso(i,k)+heso(i,k+1))
          end if
       end do
    end do
    do i = istart , iend
       if ( aa0(i) > -1.0D0 ) then
          do k = 2 , kz - kbcon(i) - 1
             kk = kz - k + 1
             if ( dby(i,kk) >= 0.0D0 ) then
                ktop(i) = kk + 1
                go to 320
             end if
          end do
          aa0(i) = -1.0D0
          go to 400
320       continue
          if ( ktop(i) > kz ) ktop(i) = kz
          if ( p(i,kbcon(i))-p(i,ktop(i)) < mincld2d(i,jslc) ) &
            aa0(i) = -1.0D0
       end if
400    continue
    end do
    !

    !------- moisture and cloud work functions
    !
    do k = 2 , kz - 1
       do i = istart , iend
          if ( aa0(i) > -1.0D0 ) then
             if ( k > kbcon(i) ) then
                if ( k < ktop(i) ) then
                   dz = -0.5D0*z(i,k-1) + 0.5D0*z(i,k+1)
                   dz1 = z(i,k) - z(i,k-1)
                   agamma = (wlhvocp)* &
                      (wlhv/(rwat*(t(i,k)**2.0D0)))*qes(i,k)
                   gamma0 = (wlhvocp)*(wlhv/(rwat*(t(i,k-1)**2.0D0)))* &
                        & qes(i,k-1)
                   qrch = qes(i,k) + (1.0D0/wlhv)* &
                       (agamma/(1.0D0+agamma))*dby(i,k)
                   qc(i,k) = (qck(i)-qrch)/(1.0D0+c0*dz) + qrch
                   pw(i,k) = c0*dz*(qc(i,k)-qrch)
                   qck(i) = qc(i,k)
                   pwav(i) = pwav(i) + pw(i,k)
                   dz1 = z(i,k) - z(i,k-1)
                   aa0(i) = aa0(i)                                     &
                        & + dz1*(gti/(cpd*(0.5D0*(t(i,k)+t(i,k-1)))))  &
                        & *dby(i,k-1)/(1.0D0+0.5D0*agamma+0.5D0*gamma0)
                   dzo = -0.5D0*zo(i,k-1) + 0.5D0*zo(i,k+1)
                   dz2 = zo(i,k) - zo(i,k-1)
                   gammo = (wlhvocp)*(wlhv/ &
                             (rwat*(tn(i,k)**2.0D0)))*qeso(i,k)
                   gammo0 = (wlhvocp)*(wlhv/(rwat*(tn(i,k-1)**2.0D0)))*  &
                        & qeso(i,k-1)
                   qrcho = qeso(i,k) + & 
                     (1.0D0/wlhv)*(gammo/(1.0D0+gammo))*dbyo(i,k)
                   qco(i,k) = (qcko(i)-qrcho)/(1.0D0+c0*dzo) + qrcho
                   pwo(i,k) = c0*dzo*(qco(i,k)-qrcho)
                   qcko(i) = qco(i,k)
                   pwavo(i) = pwavo(i) + pwo(i,k)
                   aa1(i) = aa1(i)                                       &
                        & + dz2*(gti/(cpd*(0.5D0*(tn(i,k)+tn(i,k-1)))))  &
                        & *dbyo(i,k-1)/(1.0D0+0.5D0*gammo+0.5D0*gammo0)
                end if
             end if
          end if
       end do
    end do
    !
    !
    do i = istart , iend
       if ( aa0(i) > -1.0D0 ) then
          k = ktop(i)
          dz = -0.5D0*z(i,k-1) + 0.5D0*z(i,k)
          agamma = (wlhvocp)*(wlhv/(rwat*(t(i,k)**2.0D0)))*qes(i,k)
          qrch = qes(i,k) + (1.0D0/wlhv)*(agamma/(1.0D0+agamma))*dby(i,k)
          qc(i,k) = qes(i,k)
          pw(i,k) = (qrch-qes(i,k))
          pwav(i) = pwav(i) + pw(i,k)
          !
          dz = -0.5D0*zo(i,k-1) + 0.5D0*zo(i,k)
          agamma = (wlhvocp)*(wlhv/(rwat*(tn(i,k)**2.0D0)))*qeso(i,k)
          qrcho = qeso(i,k) + (1.0D0/wlhv)* &
               (agamma/(1.0D0+agamma))*dbyo(i,k)
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
    do kk = 1 , kz/2
       do i = istart , iend
          if ( aa0(i) > -1.0D0 ) vshear(i) = vshear(i)                    &
               & + dabs((vsp(i,kk+1)-vsp(i,kk))/(z(i,kk+1)-z(i,kk)))
       end do
    end do
    do i = istart , iend
       if ( aa0(i) > -1.0D0 ) then
          vshear(i) = 1.D3*vshear(i)/dble(kz/2)
          edt(i) = 1.0D0 - (1.591D0-0.639D0*vshear(i)+  &
                            0.0953D0*(vshear(i)**2.0D0) &
                           -0.00496D0*(vshear(i)**3.0D0))

          if ( edt(i) > shrmax2d(i,jslc) ) edt(i) = shrmax2d(i,jslc)
          if ( edt(i) < shrmin2d(i,jslc) ) edt(i) = shrmin2d(i,jslc)

          edto(i) = edt(i)
          edtx(i) = edt(i)
          qrcd(i,kz) = qes(i,kz)
          hcd(i) = 0.5D0*(he(i,jmin(i))+he(i,jmin(i)+1))
          qcd(i) = 0.5D0*(q(i,jmin(i))+q(i,jmin(i)+1))
          qrcdo(i,kz) = qeso(i,kz)
          hcdo(i) = heso(i,kz)
          hcdo(i) = 0.5D0*(heo(i,jmin(i))+heo(i,jmin(i)+1))
          qcdo(i) = 0.5D0*(qo(i,jmin(i))+qo(i,jmin(i)+1))
          bu(i) = 0.0D0
          buo(i) = 0.0D0
       end if
    end do
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(i) > -1.0D0 ) then
             if ( k < jmin(i) ) then
                kk = jmin(i) - k
                dz = -(z(i,kk)-z(i,kk+2))*0.5D0
                bu(i) = bu(i) + dz*(hcd(i)-0.5D0*(hes(i,kk)+hes(i,kk+1)))
                dq = (qes(i,kk)+qes(i,kk+1))*0.5D0
                xdt = (t(i,kk)+t(i,kk+1))*0.5D0
                agamma = (wlhvocp)*(wlhv/(rwat*(xdt**2.0D0)))*dq
                dh = hcd(i) - 0.5D0*(hes(i,kk)+hes(i,kk+1))
                qrcd(i,kk) = (dq+(1.0D0/wlhv)*(agamma/(1.0D0+agamma))*dh)
                pwd(i,kk) = dkk(i,kk)*(qcd(i)-qrcd(i,kk))
                qcd(i) = qrcd(i,kk)
                pwev(i) = pwev(i) + pwd(i,kk)
                !
                dz = -(zo(i,kk)-zo(i,kk+2))*0.5D0
                buo(i) = buo(i) + dz*(hcdo(i)- &
                         0.5D0*(heso(i,kk)+heso(i,kk+1)))
                dq = (qeso(i,kk)+qeso(i,kk+1))*0.5D0
                xdt = (tn(i,kk)+tn(i,kk+1))*0.5D0
                agamma = (wlhvocp)*(wlhv/(rwat*(xdt**2.0D0)))*dq
                dh = hcdo(i) - 0.5D0*(heso(i,kk)+heso(i,kk+1))
                qrcdo(i,kk) = (dq+(1.0D0/wlhv)*(agamma/(1.0D0+agamma))*dh)
                pwdo(i,kk) = dkk(i,kk)*(qcdo(i)-qrcdo(i,kk))
                qcdo(i) = qrcdo(i,kk)
                pwevo(i) = pwevo(i) + pwdo(i,kk)
             end if
          end if
       end do
    end do
    !
    do i = istart , iend
       if ( aa0(i) > -1.0D0 ) then
          if ( bu(i) >= 0.0D0 .or. buo(i) >= 0.0D0 .or. &
             pwev(i) >= 0.0D0 .or. pwevo(i) >= 0.0D0 ) then
            aa0(i) = -1.0D0
          end if
          edt(i) = -edt(i)*pwav(i)/pwev(i)
          if ( edt(i) > edtmax2d(i,jslc) ) edt(i) = edtmax2d(i,jslc)
          if ( edt(i) < edtmin2d(i,jslc) ) edt(i) = edtmin2d(i,jslc)
          edto(i) = -edto(i)*pwavo(i)/pwevo(i)
          if ( edto(i) > edtmaxo2d(i,jslc) ) &
            edto(i) = edtmaxo2d(i,jslc)
          if ( edto(i) < edtmino2d(i,jslc) ) &
            edto(i) = edtmino2d(i,jslc)
       end if
    end do
    !
    !---  what would the change be?
    !
    do i = istart , iend
       if ( aa0(1) > -1.0D0 ) then
          k = 1
          dz = 0.5D0*(z(i,2)-z(i,1))
          dp_s = 50.0D0*(psur(i)-p(i,2))
          dellah(i,1) = edt(i)                                      &
               & *(dkk(i,1)*hcd(i)-dkk(i,1)* &
                 0.5D0*(he(i,1)+he(i,2)))*gti/dp_s
          dellaq(i,1) = edt(i)                                      &
               & *(dkk(i,1)*qrcd(i,1)-dkk(i,1)* &
                 0.5D0*(q(i,1)+q(i,2)))*gti/dp_s
          xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
          xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
          dellat(i,k) = rcpd*(dellah(i,k)-wlhv*dellaq(i,k))
          xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-wlhv*dellaq(i,k))+t(i,k)
          if ( xq(i,k) <= 0.0D0 ) xq(i,k) = 1.D-08
       end if
    end do
    !
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(1) > -1.0D0 ) then
             if ( k /= 1 .and. k < ktop(i) ) then
                dv1 = 0.5D0*(he(i,k)+he(i,k+1))
                dv2 = he(i,k)
                dv3 = 0.5D0*(he(i,k)+he(i,k-1))
                dv1q = 0.5D0*(q(i,k)+q(i,k+1))
                dv2q = q(i,k)
                dv3q = 0.5D0*(q(i,k)+q(i,k-1))
                !
                !---   specifiy detrainment of downdraft, has to be consistent
                !---   with zd calculations in soundd.
                !
                detdo = (1.0D0-dkk(i,k))*(hcd(i)-dv2)
                detdoq = (1.0D0-dkk(i,k))*(qrcd(i,k)-dv2q)
                dz = 0.5D0*(z(i,k+1)-z(i,k-1))
                !
                !---    changed due to subsidence and entrainment
                !
                aup = 1.0D0
                if ( k <= k22(i) ) aup = 0.0D0
                adw = 1.0D0
                if ( k > jmin(i) ) adw = 0.0D0
                dp_s = +50.0D0*(p(i,k-1)-p(i,k+1))
                dellah(i,k) = ((aup-adw*edt(i))*(dv1-dv2)+(aup-adw*edt(i))&
                     & *(dv2-dv3))*gti/dp_s + adw*edt(i)*detdo*gti/dp_s
                dellaq(i,k) = ((aup-adw*edt(i))* &
                    (dv1q-dv2q)+(aup-adw*edt(i))* &
                    (dv2q-dv3q))*gti/dp_s + adw*edt(i)*detdoq*gti/dp_s
                xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
                xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
                dellat(i,k) = rcpd*(dellah(i,k)-wlhv*dellaq(i,k))
                xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-wlhv*dellaq(i,k))    &
                     & + t(i,k)
                if ( xq(i,k) <= 0.0D0 ) xq(i,k) = 1.D-08
             end if
          end if
       end do
    end do
    !
    !------- cloud top
    !
    do i = istart , iend
       if ( aa0(1) > -1.0D0 ) then
          lpt = ktop(i)
          dp_s = 100.0D0*(p(i,lpt-1)-p(i,lpt))
          dv1 = 0.5D0*(he(i,lpt)+he(i,lpt-1))
          dellah(i,lpt) = (hkb(i)-dv1)*gti/dp_s
          dv1 = 0.5D0*(q(i,lpt)+q(i,lpt-1))
          dellaq(i,lpt) = (qes(i,lpt)-dv1)*gti/dp_s
          k = lpt
          xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
          xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
          dellat(i,k) = rcpd*(dellah(i,k)-wlhv*dellaq(i,k))
          xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-wlhv*dellaq(i,k))        &
               & + t(i,k)
          if ( xq(i,k) <= 0.0D0 ) xq(i,k) = 1.D-08
          xhkb(i) = dellah(i,kbcon(i))*mbdt + hkb(i)
          xqkb(i) = dellaq(i,kbcon(i))*mbdt + qkb(i)
          if ( xqkb(i) <= 0.0D0 ) xqkb(i) = 1.D-08
       end if
    end do
    !
    !---  environmental conditions, first heights
    !
    do k = 1 , kz
       do i = istart , iend
          if ( aa0(1) > -1.0D0 ) then
             iph = 1
             if ( xt(i,k) <= tcrit ) iph = 2
             e = dexp(ae(iph)-be(iph)/xt(i,k))
             xqes(i,k) = ep2*e/(100.0D0*p(i,k)-(1.0D0-ep2)*e)
             if ( xqes(i,k) <= 1.D-08 ) xqes(i,k) = 1.D-08
             if ( xq(i,k) > xqes(i,k) ) xq(i,k) = xqes(i,k)
             xtv(i,k) = xt(i,k) + 0.608D0*xq(i,k)*xt(i,k)
          end if
       end do
    end do
    !     bug fix
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(1) > -1.0D0 ) &
            xqrcd(i,k) = 0.5D0*(xqes(i,k)+xqes(i,k+1))
       end do
    end do
    !
    do i = istart , iend
       if ( aa0(1) > -1.0D0 ) xz(i,1) = z1(i)                            &
            & - (dlog(p(i,1))-dlog(psur(i)))   &
            & *rgas*xtv(i,1)*rgti
    end do
    do k = 2 , kz
       do i = istart , iend
          if ( aa0(1) > -1.0D0 ) then
             tvbar = 0.5D0*xtv(i,k) + 0.5D0*xtv(i,k-1)
             xz(i,k) = xz(i,k-1) - (dlog(p(i,k))-dlog(p(i,k-1)))         &
                  & *rgas*tvbar*rgti
          end if
       end do
    end do
    !
    !---  moist static energy
    !
    do k = 1 , kz
       do i = istart , iend
          if ( aa0(1) > -1.0D0 ) then
             xhes(i,k) = gti*xz(i,k) + cpd*xt(i,k) + wlhv*xqes(i,k)
             if ( xhe(i,k) >= xhes(i,k) ) xhe(i,k) = xhes(i,k)
          end if
       end do
    end do
    !
    !
    !**************************** static control
    !
    do i = istart , iend
       if ( aa0(1) > -1.0D0 ) then
          xqck(i) = xqkb(i)
          xdby(i,kz) = xhkb(i) - xhes(i,kz)
       end if
    end do
    !
    !------- moisture and cloud work functions
    !
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(i) >= 0.0D0 ) then
             xdby(i,k) = xhkb(i) - 0.5D0*(xhes(i,k)+xhes(i,k+1))
             if ( k > kbcon(i) .and. k < ktop(i) ) then
                dz = -0.5D0*xz(i,k-1) + 0.5D0*xz(i,k+1)
                dz1 = xz(i,k) - xz(i,k-1)
                agamma = (wlhvocp)*(wlhv/(rwat*(xt(i,k)**2.0D0)))*xqes(i,k)
                gamma0 = (wlhvocp)*(wlhv/(rwat*(xt(i,k-1)**2.0D0)))* &
                     & xqes(i,k-1)
                qrch = xqes(i,k) + (1.0D0/wlhv)*(agamma/(1.0D0+agamma))*   &
                     & xdby(i,k)
                xqc(i,k) = (xqck(i)-qrch)/(1.0D0+c0*dz) + qrch
                xpw(i,k) = c0*dz*(xqc(i,k)-qrch)
                xqck(i) = xqc(i,k)
                xpwav(i) = xpwav(i) + xpw(i,k)
                xaa0(i) = xaa0(i)                                         &
                     & + dz1*(gti/(cpd*(0.5D0*(xt(i,k)+xt(i,k-1)))))      &
                     & *xdby(i,k-1)/(1.0D0+0.5D0*agamma+0.5D0*gamma0)
             end if
          end if
       end do
    end do
    do i = istart , iend
       if ( aa0(i) >= 0.0D0 ) then
          k = ktop(i)
          dz = -0.5D0*xz(i,k-1) + 0.5D0*xz(i,k)
          agamma = (wlhvocp)*(wlhv/(rwat*(xt(i,k)**2.0D0)))*xqes(i,k)
          qrch = xqes(i,k) + (1.0D0/wlhv)*(agamma/(1.0D0+agamma))*xdby(i,k)
          xqc(i,k) = xqes(i,k)
          xpw(i,k) = (qrch-xqes(i,k))
          xpwav(i) = xpwav(i) + xpw(i,k)
          xqrcd(i,kz) = xqes(i,kz)
          xhcd(i) = 0.5D0*(xhe(i,jmin(i))+xhe(i,jmin(i)+1))
          xqcd(i) = 0.5D0*(xq(i,jmin(i))+xq(i,jmin(i)+1))
          xpwev(i) = 0.0D0
          bu(i) = 0.0D0
       end if
    end do
    !
    !------- downdraft calculations
    !
    !
    !---  downdraft moisture properties
    !
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(i) >= 0.0D0 ) then
             if ( k < jmin(i) ) then
                kk = jmin(i) - k
                dz = -(xz(i,kk)-xz(i,kk+2))*0.5D0
                bu(i) = bu(i) + dz*(xhcd(i)-0.5D0*(xhes(i,kk)+xhes(i,kk+1)))
                dq = (xqes(i,kk)+xqes(i,kk+1))*0.5D0
                xdt = (xt(i,kk)+xt(i,kk+1))*0.5D0
                agamma = (wlhvocp)*(wlhv/(rwat*(xdt**2.0D0)))*dq
                dh = xhcd(i) - 0.5D0*(xhes(i,kk)+xhes(i,kk+1))
                xqrcd(i,kk) = (dq+(1.0D0/wlhv)*(agamma/(1.0D0+agamma))*dh)
                xpwd(i,kk) = dkk(i,kk)*(xqcd(i)-xqrcd(i,kk))
                xqcd(i) = xqrcd(i,kk)
                xpwev(i) = xpwev(i) + xpwd(i,kk)
             end if
          end if
       end do
    end do
    do i = istart , iend
       if ( aa0(i) >= 0.0D0 ) then
          if ( bu(i) >= 0.0D0 ) then
             aa0(i) = -1.0D0
             go to 500
          end if
          if ( xpwev(i) /= 0.0D0 ) edtx(i) = -edtx(i)*xpwav(i)/xpwev(i)
          if ( edtx(i) > edtmaxx2d(i,jslc) ) edtx(i)                   &
               & = edtmaxx2d(i,jslc)
          if ( edtx(i) < edtminx2d(i,jslc) ) edtx(i)                   &
               & = edtminx2d(i,jslc)
       end if
500    continue
    end do
    !
    !
    !---  downdraft cloudwork functions
    !
    !
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(i) >= 0.0D0 ) then
             if ( k < jmin(i) ) then
                kk = jmin(i) - k
                !
                !---          original
                !
                gamma1 = (wlhvocp)*(wlhv/(rwat*(t(i,kk)**2.0D0)))*qes(i,kk)
                gamma2 = (wlhvocp)*(wlhv/(rwat*(t(i,kk+1)**2.0D0)))*          &
                     & qes(i,kk+1)
                dhh = hcd(i)
                xdt = 0.5D0*(t(i,kk)+t(i,kk+1))
                dg = 0.5D0*(gamma1+gamma2)
                dh = 0.5D0*(hes(i,kk)+hes(i,kk+1))
                dz = (z(i,kk)-z(i,kk+1))*dkk(i,kk)
                aa0(i) = aa0(i) + edt(i)*dz*(gti/(cpd*xdt))*((dhh-dh)/    &
                     & (1.0D0+dg))
                !
                !---          modified by larger scale
                !
                gamma1 = (wlhvocp)*(wlhv/(rwat*(tn(i,kk)**2.0D0)))*qeso(i,kk)
                gamma2 = (wlhvocp)*(wlhv/(rwat*(tn(i,kk+1)**2.0D0)))*         &
                     & qeso(i,kk+1)
                dhh = hcdo(i)
                xdt = 0.5D0*(tn(i,kk)+tn(i,kk+1))
                dg = 0.5D0*(gamma1+gamma2)
                dh = 0.5D0*(heso(i,kk)+heso(i,kk+1))
                dz = (zo(i,kk)-zo(i,kk+1))*dkk(i,kk)
                aa1(i) = aa1(i) + edto(i)*dz*(gti/(cpd*xdt))              &
                     & *((dhh-dh)/(1.0D0+dg))
                !
                !---          modified by cloud
                !
                gamma1 = (wlhvocp)*(wlhv/(rwat*(xt(i,kk)**2.0D0)))*xqes(i,kk)
                gamma2 = (wlhvocp)*(wlhv/(rwat*(xt(i,kk+1)**2.0D0)))*         &
                     & xqes(i,kk+1)
                dhh = xhcd(i)
                xdt = 0.5D0*(xt(i,kk)+xt(i,kk+1))
                dg = 0.5D0*(gamma1+gamma2)
                dh = 0.5D0*(xhes(i,kk)+xhes(i,kk+1))
                dz = (xz(i,kk)-xz(i,kk+1))*dkk(i,kk)
                xaa0(i) = xaa0(i) + edtx(i)*dz*(gti/(cpd*xdt))            &
                     & *((dhh-dh)/(1.0D0+dg))
             end if
          end if
       end do
    end do
    !
    !---  large scale forcing
    !
    do i = istart , iend
       if ( aa0(i) >= 0.0D0 ) then
          if ( igcc == 1 ) then
             f = (aa1(i)-aa0(i))/dtime  ! Arakawa-Schubert closure
          else if ( igcc == 2 ) then
             f = aa0(i)/dtauc2d(i,jslc)   ! Fritsch-Chappell closure
          end if
          xk = (xaa0(i)-aa0(i))/mbdt
          xmb(i) = -f/xk
          if ( f <= 0.0D0 .or. xk >= 0.0D0 ) xmb(i) = 0.0D0
       end if
    end do
    !chem2
    mflx(i,1) = xmb(i)
    mflx(i,2) = xmb(i)*edt(i)
    !chem2_
    !
    !---  feedback
    !
    do k = 1 , kz
       do i = istart , iend
          if ( aa0(i) >= 0.0D0 ) then
             if ( k <= ktop(i) ) then
                outtes = dellat(i,k)*xmb(i)*86400.0D0
                if ( (outtes > htmax2d(i,jslc)) .or.  &
                     & (outtes < htmin2d(i,jslc)) ) then
                   xmb(i) = 0.0D0
                   aa0(i) = -1.0D0
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
       if ( aa0(i) >= 0.0D0 ) then

          if ( ktop(i) > 1 .and. kbcon(i) > 1 ) then
             kclth = ktop(i) - kbcon(i) + 1
             akclth = 1.0D0/dble(kclth)
             do k = kbcon(i) , ktop(i)
                kk = kz - k + 1
                cldlwc(i,kk) = cllwcv
                cldfra(i,kk) = 1.0D0 - (1.0D0-clfrcv)**akclth
             end do
             !chem2
             !chem2      define convection  base and top for tracers
             if ( ichem == 1 ) then
                if ( ktop(i) > 1 .and. k22(i) >= 1 ) then
                   icumtop(i,jslc) = kzp1 - ktop(i)
                   icumbot(i,jslc) = kzp1 - k22(i)
                   icumdwd(i,jslc) = kzp1 - jmin(i)
                end if
             end if
             !chem2_
          end if
       end if

    end do
    call time_end(subroutine_name,idindx)
  end subroutine cup
  !
  !
  !
  subroutine minimi(array,iy,kz,ks,kend,kt,istart,iend)
    !
    implicit none
    !
    integer :: iend , istart , iy , kend , kz
    real(8) , dimension(iy,kz) :: array
    integer , dimension(iy) :: ks , kt
    intent (in) array , iend , istart , iy , kend , ks , kz
    intent (out) kt
    !
    integer :: i , k
    real(8) :: x
    character (len=50) :: subroutine_name='minimi'
    integer :: idindx=0
    !
    call time_begin(subroutine_name,idindx)
    !
    do i = istart , iend
       kt(i) = ks(i)
       x = array(i,ks(i))
       do k = ks(i) + 1 , kend
          if ( array(i,k) < x ) then
             x = array(i,k)
             kt(i) = k
          end if
       end do
    end do
    !
    call time_end(subroutine_name,idindx) 
  end subroutine minimi
  !
  !
  !
  subroutine maximi(array,iy,kz,ks,ke,imax,istart,iend)

    implicit none
    !
    integer :: iend , istart , iy , ke , ks , kz
    real(8) , dimension(iy,kz) :: array
    integer , dimension(iy) :: imax
    intent (in) array , iend , istart , iy , ke , ks , kz
    intent (out) imax
    !
    integer :: i , k
    real(8) :: x , xar
    !
    character (len=50) :: subroutine_name='maximi'
    integer :: idindx=0
    !
    call time_begin(subroutine_name,idindx)
    do i = istart , iend
       imax(i) = ks
       x = array(i,ks)
       do k = ks , ke
          xar = array(i,k)
          if ( xar >= x ) then
             x = xar
             imax(i) = k
          end if
       end do
    end do
    !
    call time_end(subroutine_name,idindx) 
  end subroutine maximi

end module mod_cu_grell

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
    cldlwc = d_zero
    cldfra = d_zero

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
       qcrit(i) = d_zero
       pret(i) = d_zero
    end do

    do k = 1 , kz
       do i = 2 + icut , iym2 - icut
          kk = kz - k + 1

          us = gwnd%usk(i,kk)
          vs = gwnd%vsk(i,kk)

          t(i,k) = atm2%t(i,kk,j)/sps2%ps(i,j)
          q(i,k) = atm2%qv(i,kk,j)/sps2%ps(i,j)
          if ( q(i,k) < 1.0D-08 ) q(i,k) = 1.0D-08
          tn(i,k) = t(i,k) + (tten(i,kk))/sps2%ps(i,j)*dtime
          qo(i,k) = q(i,k) + (qten(i,kk))/sps2%ps(i,j)*dtime
          p(i,k) = d_10*sps2%ps(i,j)*a(kk) + d_10*r8pt
          vsp(i,k) = dsqrt(us**d_two+vs**d_two)
          if ( qo(i,k) < 1.0D-08 ) qo(i,k) = 1.0D-08
          !
          po(i,k) = p(i,k)
          psur(i) = d_10*sps2%ps(i,j) + d_10*r8pt
          outt(i,k) = d_zero
          pkk = psur(i) - po(i,k)
          if ( pkk <= pkdcut ) kdet(i) = kdet(i) + 1
          outq(i,k) = d_zero
          ter11(i) = mddom%ht(i,j)*rgti
          if ( ter11(i) <= d_zero ) ter11(i) = 1.0D-05
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
          if ( pret(i) > d_zero ) then
             kk = kz - k + 1
             tten(i,kk) = sps2%ps(i,j)*outt(i,k) + tten(i,kk)
             qten(i,kk) = sps2%ps(i,j)*outq(i,k) + qten(i,kk)
          end if
       end do
    end do
    !
    !---  rain in cm.
    !
    calc = d_half
    iconj = 0
    do i = istart , iend
       if ( pret(i) > d_zero ) then
          sfsta%rainc(i,j) = sfsta%rainc(i,j) + pret(i)*calc*dt
          !         print *,'sfsta%rainc(',i,j,')=',sfsta%rainc(i,j)
          iconj = iconj + 1
          !.....................precipitation rate for bats (mm/s)
          aprdiv = dble(nbatst)
          if ( jyear == jyear0 .and. ktau == 0 ) aprdiv = d_one
          prainx = pret(i)*calc*dt
          pptc(i,j) = pptc(i,j) + prainx/(dtmin*minph)/aprdiv
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
    mbdt = dtime*5.0D-03
    c0 = 0.002D0
    f = -d_one
    xk = -d_one
    !
    !---  environmental conditions, first heights
    !
    do k = 1 , kz
       do i = istart , iend
          dkk(i,k) = d_one
          iph = 1
          ipho = 1
          if ( t(i,k) <= tcrit ) iph = 2
          if ( tn(i,k) <= tcrit ) ipho = 2
          e = dexp(ae(iph)-be(iph)/t(i,k))
          eo = dexp(ae(ipho)-be(ipho)/tn(i,k))
          qes(i,k) = ep2*e/(d_100*p(i,k)-(d_one-ep2)*e)
          qeso(i,k) = ep2*eo/(d_100*po(i,k)-(d_one-ep2)*eo)
          if ( qes(i,k) <= 1.0D-08 ) qes(i,k) = 1.0D-08
          if ( q(i,k) > qes(i,k) ) q(i,k) = qes(i,k)
          if ( qeso(i,k) <= 1.0D-08 ) qeso(i,k) = 1.0D-08
          if ( qo(i,k) > qeso(i,k) ) qo(i,k) = qeso(i,k)
          tv(i,k) = t(i,k) + 0.608D0*q(i,k)*t(i,k)
          tvo(i,k) = tn(i,k) + 0.608D0*qo(i,k)*tn(i,k)
       end do
    end do
    do i = 1 , iy
       hkb(i) = d_zero ! EES
       qkb(i) = d_zero
       hkbo(i) = d_zero
       qkbo(i) = d_zero
       xhkb(i) = d_zero
       xqkb(i) = d_zero
       edt(i) = d_zero
       edto(i) = d_zero
       edtx(i) = d_zero
    end do
    do i = istart , iend
       !       hkb(i)=d_zero
       !       qkb(i)=d_zero
       !       hkbo(i)=d_zero
       !       qkbo(i)=d_zero
       !       xhkb(i)=d_zero
       !       xqkb(i)=d_zero
       aa1(i) = d_zero
       aa0(i) = d_zero
       if ( icup==99 .or. icup==98 ) then
         if (cumcon%cuscheme(i,jslc)/=2 ) then
           aa0(i) = -d_one
         end if
       end if
       if ( qcrit(i) <= d_zero ) then
         aa0(i) = -d_one
       end if
       xaa0(i) = d_zero
       xpwav(i) = d_zero
       xpwev(i) = d_zero
       pwav(i) = d_zero
       pwev(i) = d_zero
       pwavo(i) = d_zero
       pwevo(i) = d_zero
       k22(i) = 1
       ktop(i) = 1
       kbcon(i) = 1
       kb(i) = 1
       kds(i) = 1
       jmin(i) = 1
       !       edt(i)=d_zero
       !       edto(i)=d_zero
       !       edtx(i)=d_zero
       xmb(i) = d_zero
       vshear(i) = d_zero
       z(i,1) = z1(i) - (dlog(p(i,1))-dlog(psur(i)))*rgas*tv(i,1)*rgti
       zo(i,1) = z1(i) - (dlog(po(i,1))-dlog(psur(i)))*rgas*tvo(i,1)   &
            & *rgti
    end do
    do k = 2 , kz
       do i = istart , iend
          tvbar = d_half*tv(i,k) + d_half*tv(i,k-1)
          z(i,k) = z(i,k-1) - (dlog(p(i,k))-dlog(p(i,k-1)))             &
               & *rgas*tvbar*rgti
          tvbaro = d_half*tvo(i,k) + d_half*tvo(i,k-1)
          zo(i,k) = zo(i,k-1) - (dlog(po(i,k))-dlog(po(i,k-1)))         &
               & *rgas*tvbaro*rgti
       end do
    end do
    !
    !---  moist static energy
    !
    do k = 1 , kz
       do i = istart , iend
          cldlwc(i,k) = d_zero
          cldfra(i,k) = d_zero
          pw(i,k) = d_zero
          xpw(i,k) = d_zero
          pwo(i,k) = d_zero
          qc(i,k) = d_zero
          xqc(i,k) = d_zero
          qco(i,k) = d_zero
          pwd(i,k) = d_zero
          pwdo(i,k) = d_zero
          xpwd(i,k) = d_zero
          dellah(i,k) = d_zero
          dellaq(i,k) = d_zero
          dellat(i,k) = d_zero
          he(i,k) = gti*z(i,k) + cpd*t(i,k) + wlhv*q(i,k)
          hes(i,k) = gti*z(i,k) + cpd*t(i,k) + wlhv*qes(i,k)
          if ( he(i,k) >= hes(i,k) ) he(i,k) = hes(i,k)
          heo(i,k) = gti*zo(i,k) + cpd*tn(i,k) + wlhv*qo(i,k)
          heso(i,k) = gti*zo(i,k) + cpd*tn(i,k) + wlhv*qeso(i,k)
          if ( heo(i,k) >= heso(i,k) ) heo(i,k) = heso(i,k)
          xt(i,k) = t(i,k)
          xq(i,k) = q(i,k)
          xhe(i,k) = he(i,k)
          if ( k /= kz ) qrcd(i,k) = d_half*(qes(i,k)+qes(i,k+1))
          if ( k /= kz ) qrcdo(i,k) = d_half*(qeso(i,k)+qeso(i,k+1))
       end do
    end do
    !
    !------- determine level with highest moist static energy content.
    !
    call maximi(he,iy,kz,1,kbmax2d(i,jslc),k22,istart,iend)
    do i = istart , iend
       if ( aa0(i) >= d_zero ) then
          if ( k22(i) >= kbmax2d(i,jslc) ) then
             aa0(i) = -d_one
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
       if ( aa0(i) >= d_zero ) then
          do k = 1 , kdet(i)
             kk = kdet(i) - k + 1
             !           dkk(i,kk)=.75*dkk(i,kk+1)
             dkk(i,k) = d_one - dble(kk)/dble(kdet(i))
          end do
120       continue
          kb(i) = k22(i)
          !----------------------------------
          kbcon(i) = kb(i)
140       continue
          dh = d_half*hes(i,kbcon(i)) + d_half*hes(i,kbcon(i)+1)
          if ( hkb(i) < dh ) then
             kbcon(i) = kbcon(i) + 1
             if ( kbcon(i) > kbmax2d(i,jslc) ) then
                aa0(i) = -d_one
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
                aa0(i) = -d_one
                go to 200
             end if
             !ictp_
             dh = d_half*heso(i,kbcono) + d_half*heso(i,kbcono+1)
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
                      aa0(i) = -d_one
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
             aa0(i) = -d_one
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
          if ( aa0(i) > -d_one ) then
             dby(i,k) = hkb(i) - d_half*(hes(i,k)+hes(i,k+1))
             dbyo(i,k) = hkbo(i) - d_half*(heso(i,k)+heso(i,k+1))
          end if
       end do
    end do
    do i = istart , iend
       if ( aa0(i) > -d_one ) then
          do k = 2 , kz - kbcon(i) - 1
             kk = kz - k + 1
             if ( dby(i,kk) >= d_zero ) then
                ktop(i) = kk + 1
                go to 320
             end if
          end do
          aa0(i) = -d_one
          go to 400
320       continue
          if ( ktop(i) > kz ) ktop(i) = kz
          if ( p(i,kbcon(i))-p(i,ktop(i)) < mincld2d(i,jslc) ) &
            aa0(i) = -d_one
       end if
400    continue
    end do
    !

    !------- moisture and cloud work functions
    !
    do k = 2 , kz - 1
       do i = istart , iend
          if ( aa0(i) > -d_one ) then
             if ( k > kbcon(i) ) then
                if ( k < ktop(i) ) then
                   dz = -d_half*z(i,k-1) + d_half*z(i,k+1)
                   dz1 = z(i,k) - z(i,k-1)
                   agamma = (wlhvocp)* &
                      (wlhv/(rwat*(t(i,k)**d_two)))*qes(i,k)
                   gamma0 = (wlhvocp)*(wlhv/(rwat*(t(i,k-1)**d_two)))* &
                        & qes(i,k-1)
                   qrch = qes(i,k) + (d_one/wlhv)* &
                       (agamma/(d_one+agamma))*dby(i,k)
                   qc(i,k) = (qck(i)-qrch)/(d_one+c0*dz) + qrch
                   pw(i,k) = c0*dz*(qc(i,k)-qrch)
                   qck(i) = qc(i,k)
                   pwav(i) = pwav(i) + pw(i,k)
                   dz1 = z(i,k) - z(i,k-1)
                   aa0(i) = aa0(i)                                     &
                        & + dz1*(gti/(cpd*(d_half*(t(i,k)+t(i,k-1)))))  &
                        & *dby(i,k-1)/(d_one+d_half*agamma+d_half*gamma0)
                   dzo = -d_half*zo(i,k-1) + d_half*zo(i,k+1)
                   dz2 = zo(i,k) - zo(i,k-1)
                   gammo = (wlhvocp)*(wlhv/ &
                             (rwat*(tn(i,k)**d_two)))*qeso(i,k)
                   gammo0 = (wlhvocp)*(wlhv/(rwat*(tn(i,k-1)**d_two)))*  &
                        & qeso(i,k-1)
                   qrcho = qeso(i,k) + & 
                     (d_one/wlhv)*(gammo/(d_one+gammo))*dbyo(i,k)
                   qco(i,k) = (qcko(i)-qrcho)/(d_one+c0*dzo) + qrcho
                   pwo(i,k) = c0*dzo*(qco(i,k)-qrcho)
                   qcko(i) = qco(i,k)
                   pwavo(i) = pwavo(i) + pwo(i,k)
                   aa1(i) = aa1(i)                                       &
                        & + dz2*(gti/(cpd*(d_half*(tn(i,k)+tn(i,k-1)))))  &
                        & *dbyo(i,k-1)/(d_one+d_half*gammo+d_half*gammo0)
                end if
             end if
          end if
       end do
    end do
    !
    !
    do i = istart , iend
       if ( aa0(i) > -d_one ) then
          k = ktop(i)
          dz = -d_half*z(i,k-1) + d_half*z(i,k)
          agamma = (wlhvocp)*(wlhv/(rwat*(t(i,k)**d_two)))*qes(i,k)
          qrch = qes(i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*dby(i,k)
          qc(i,k) = qes(i,k)
          pw(i,k) = (qrch-qes(i,k))
          pwav(i) = pwav(i) + pw(i,k)
          !
          dz = -d_half*zo(i,k-1) + d_half*zo(i,k)
          agamma = (wlhvocp)*(wlhv/(rwat*(tn(i,k)**d_two)))*qeso(i,k)
          qrcho = qeso(i,k) + (d_one/wlhv)* &
               (agamma/(d_one+agamma))*dbyo(i,k)
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
          if ( aa0(i) > -d_one ) vshear(i) = vshear(i)                    &
               & + dabs((vsp(i,kk+1)-vsp(i,kk))/(z(i,kk+1)-z(i,kk)))
       end do
    end do
    do i = istart , iend
       if ( aa0(i) > -d_one ) then
          vshear(i) = d_1000*vshear(i)/dble(kz/2)
          edt(i) = d_one - (1.591D0-0.639D0*vshear(i)+  &
                            0.0953D0*(vshear(i)**d_two) &
                           -0.00496D0*(vshear(i)**d_three))

          if ( edt(i) > shrmax2d(i,jslc) ) edt(i) = shrmax2d(i,jslc)
          if ( edt(i) < shrmin2d(i,jslc) ) edt(i) = shrmin2d(i,jslc)

          edto(i) = edt(i)
          edtx(i) = edt(i)
          qrcd(i,kz) = qes(i,kz)
          hcd(i) = d_half*(he(i,jmin(i))+he(i,jmin(i)+1))
          qcd(i) = d_half*(q(i,jmin(i))+q(i,jmin(i)+1))
          qrcdo(i,kz) = qeso(i,kz)
          hcdo(i) = heso(i,kz)
          hcdo(i) = d_half*(heo(i,jmin(i))+heo(i,jmin(i)+1))
          qcdo(i) = d_half*(qo(i,jmin(i))+qo(i,jmin(i)+1))
          bu(i) = d_zero
          buo(i) = d_zero
       end if
    end do
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(i) > -d_one ) then
             if ( k < jmin(i) ) then
                kk = jmin(i) - k
                dz = -(z(i,kk)-z(i,kk+2))*d_half
                bu(i) = bu(i) + dz*(hcd(i)-d_half*(hes(i,kk)+hes(i,kk+1)))
                dq = (qes(i,kk)+qes(i,kk+1))*d_half
                xdt = (t(i,kk)+t(i,kk+1))*d_half
                agamma = (wlhvocp)*(wlhv/(rwat*(xdt**d_two)))*dq
                dh = hcd(i) - d_half*(hes(i,kk)+hes(i,kk+1))
                qrcd(i,kk) = (dq+(d_one/wlhv)*(agamma/(d_one+agamma))*dh)
                pwd(i,kk) = dkk(i,kk)*(qcd(i)-qrcd(i,kk))
                qcd(i) = qrcd(i,kk)
                pwev(i) = pwev(i) + pwd(i,kk)
                !
                dz = -(zo(i,kk)-zo(i,kk+2))*d_half
                buo(i) = buo(i) + dz*(hcdo(i)- &
                         d_half*(heso(i,kk)+heso(i,kk+1)))
                dq = (qeso(i,kk)+qeso(i,kk+1))*d_half
                xdt = (tn(i,kk)+tn(i,kk+1))*d_half
                agamma = (wlhvocp)*(wlhv/(rwat*(xdt**d_two)))*dq
                dh = hcdo(i) - d_half*(heso(i,kk)+heso(i,kk+1))
                qrcdo(i,kk) = (dq+(d_one/wlhv)*(agamma/(d_one+agamma))*dh)
                pwdo(i,kk) = dkk(i,kk)*(qcdo(i)-qrcdo(i,kk))
                qcdo(i) = qrcdo(i,kk)
                pwevo(i) = pwevo(i) + pwdo(i,kk)
             end if
          end if
       end do
    end do
    !
    do i = istart , iend
       if ( aa0(i) > -d_one ) then
          if ( bu(i) >= d_zero .or. buo(i) >= d_zero .or. &
             pwev(i) >= d_zero .or. pwevo(i) >= d_zero ) then
            aa0(i) = -d_one
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
       if ( aa0(1) > -d_one ) then
          k = 1
          dz = d_half*(z(i,2)-z(i,1))
          dp_s = 50.0D0*(psur(i)-p(i,2))
          dellah(i,1) = edt(i)                                      &
               & *(dkk(i,1)*hcd(i)-dkk(i,1)* &
                 d_half*(he(i,1)+he(i,2)))*gti/dp_s
          dellaq(i,1) = edt(i)                                      &
               & *(dkk(i,1)*qrcd(i,1)-dkk(i,1)* &
                 d_half*(q(i,1)+q(i,2)))*gti/dp_s
          xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
          xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
          dellat(i,k) = rcpd*(dellah(i,k)-wlhv*dellaq(i,k))
          xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-wlhv*dellaq(i,k))+t(i,k)
          if ( xq(i,k) <= d_zero ) xq(i,k) = 1.0D-08
       end if
    end do
    !
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(1) > -d_one ) then
             if ( k /= 1 .and. k < ktop(i) ) then
                dv1 = d_half*(he(i,k)+he(i,k+1))
                dv2 = he(i,k)
                dv3 = d_half*(he(i,k)+he(i,k-1))
                dv1q = d_half*(q(i,k)+q(i,k+1))
                dv2q = q(i,k)
                dv3q = d_half*(q(i,k)+q(i,k-1))
                !
                !---   specifiy detrainment of downdraft, has to be consistent
                !---   with zd calculations in soundd.
                !
                detdo = (d_one-dkk(i,k))*(hcd(i)-dv2)
                detdoq = (d_one-dkk(i,k))*(qrcd(i,k)-dv2q)
                dz = d_half*(z(i,k+1)-z(i,k-1))
                !
                !---    changed due to subsidence and entrainment
                !
                aup = d_one
                if ( k <= k22(i) ) aup = d_zero
                adw = d_one
                if ( k > jmin(i) ) adw = d_zero
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
                if ( xq(i,k) <= d_zero ) xq(i,k) = 1.0D-08
             end if
          end if
       end do
    end do
    !
    !------- cloud top
    !
    do i = istart , iend
       if ( aa0(1) > -d_one ) then
          lpt = ktop(i)
          dp_s = d_100*(p(i,lpt-1)-p(i,lpt))
          dv1 = d_half*(he(i,lpt)+he(i,lpt-1))
          dellah(i,lpt) = (hkb(i)-dv1)*gti/dp_s
          dv1 = d_half*(q(i,lpt)+q(i,lpt-1))
          dellaq(i,lpt) = (qes(i,lpt)-dv1)*gti/dp_s
          k = lpt
          xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
          xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
          dellat(i,k) = rcpd*(dellah(i,k)-wlhv*dellaq(i,k))
          xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-wlhv*dellaq(i,k))        &
               & + t(i,k)
          if ( xq(i,k) <= d_zero ) xq(i,k) = 1.0D-08
          xhkb(i) = dellah(i,kbcon(i))*mbdt + hkb(i)
          xqkb(i) = dellaq(i,kbcon(i))*mbdt + qkb(i)
          if ( xqkb(i) <= d_zero ) xqkb(i) = 1.0D-08
       end if
    end do
    !
    !---  environmental conditions, first heights
    !
    do k = 1 , kz
       do i = istart , iend
          if ( aa0(1) > -d_one ) then
             iph = 1
             if ( xt(i,k) <= tcrit ) iph = 2
             e = dexp(ae(iph)-be(iph)/xt(i,k))
             xqes(i,k) = ep2*e/(d_100*p(i,k)-(d_one-ep2)*e)
             if ( xqes(i,k) <= 1.0D-08 ) xqes(i,k) = 1.0D-08
             if ( xq(i,k) > xqes(i,k) ) xq(i,k) = xqes(i,k)
             xtv(i,k) = xt(i,k) + 0.608D0*xq(i,k)*xt(i,k)
          end if
       end do
    end do
    !     bug fix
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(1) > -d_one ) &
            xqrcd(i,k) = d_half*(xqes(i,k)+xqes(i,k+1))
       end do
    end do
    !
    do i = istart , iend
       if ( aa0(1) > -d_one ) xz(i,1) = z1(i)                            &
            & - (dlog(p(i,1))-dlog(psur(i)))   &
            & *rgas*xtv(i,1)*rgti
    end do
    do k = 2 , kz
       do i = istart , iend
          if ( aa0(1) > -d_one ) then
             tvbar = d_half*xtv(i,k) + d_half*xtv(i,k-1)
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
          if ( aa0(1) > -d_one ) then
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
       if ( aa0(1) > -d_one ) then
          xqck(i) = xqkb(i)
          xdby(i,kz) = xhkb(i) - xhes(i,kz)
       end if
    end do
    !
    !------- moisture and cloud work functions
    !
    do k = 1 , kz - 1
       do i = istart , iend
          if ( aa0(i) >= d_zero ) then
             xdby(i,k) = xhkb(i) - d_half*(xhes(i,k)+xhes(i,k+1))
             if ( k > kbcon(i) .and. k < ktop(i) ) then
                dz = -d_half*xz(i,k-1) + d_half*xz(i,k+1)
                dz1 = xz(i,k) - xz(i,k-1)
                agamma = (wlhvocp)*(wlhv/(rwat*(xt(i,k)**d_two)))*xqes(i,k)
                gamma0 = (wlhvocp)*(wlhv/(rwat*(xt(i,k-1)**d_two)))* &
                     & xqes(i,k-1)
                qrch = xqes(i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*   &
                     & xdby(i,k)
                xqc(i,k) = (xqck(i)-qrch)/(d_one+c0*dz) + qrch
                xpw(i,k) = c0*dz*(xqc(i,k)-qrch)
                xqck(i) = xqc(i,k)
                xpwav(i) = xpwav(i) + xpw(i,k)
                xaa0(i) = xaa0(i)                                         &
                     & + dz1*(gti/(cpd*(d_half*(xt(i,k)+xt(i,k-1)))))      &
                     & *xdby(i,k-1)/(d_one+d_half*agamma+d_half*gamma0)
             end if
          end if
       end do
    end do
    do i = istart , iend
       if ( aa0(i) >= d_zero ) then
          k = ktop(i)
          dz = -d_half*xz(i,k-1) + d_half*xz(i,k)
          agamma = (wlhvocp)*(wlhv/(rwat*(xt(i,k)**d_two)))*xqes(i,k)
          qrch = xqes(i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*xdby(i,k)
          xqc(i,k) = xqes(i,k)
          xpw(i,k) = (qrch-xqes(i,k))
          xpwav(i) = xpwav(i) + xpw(i,k)
          xqrcd(i,kz) = xqes(i,kz)
          xhcd(i) = d_half*(xhe(i,jmin(i))+xhe(i,jmin(i)+1))
          xqcd(i) = d_half*(xq(i,jmin(i))+xq(i,jmin(i)+1))
          xpwev(i) = d_zero
          bu(i) = d_zero
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
          if ( aa0(i) >= d_zero ) then
             if ( k < jmin(i) ) then
                kk = jmin(i) - k
                dz = -(xz(i,kk)-xz(i,kk+2))*d_half
                bu(i) = bu(i) + dz*(xhcd(i)-d_half*(xhes(i,kk)+xhes(i,kk+1)))
                dq = (xqes(i,kk)+xqes(i,kk+1))*d_half
                xdt = (xt(i,kk)+xt(i,kk+1))*d_half
                agamma = (wlhvocp)*(wlhv/(rwat*(xdt**d_two)))*dq
                dh = xhcd(i) - d_half*(xhes(i,kk)+xhes(i,kk+1))
                xqrcd(i,kk) = (dq+(d_one/wlhv)*(agamma/(d_one+agamma))*dh)
                xpwd(i,kk) = dkk(i,kk)*(xqcd(i)-xqrcd(i,kk))
                xqcd(i) = xqrcd(i,kk)
                xpwev(i) = xpwev(i) + xpwd(i,kk)
             end if
          end if
       end do
    end do
    do i = istart , iend
       if ( aa0(i) >= d_zero ) then
          if ( bu(i) >= d_zero ) then
             aa0(i) = -d_one
             go to 500
          end if
          if ( xpwev(i) /= d_zero ) edtx(i) = -edtx(i)*xpwav(i)/xpwev(i)
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
          if ( aa0(i) >= d_zero ) then
             if ( k < jmin(i) ) then
                kk = jmin(i) - k
                !
                !---          original
                !
                gamma1 = (wlhvocp)*(wlhv/(rwat*(t(i,kk)**d_two)))*qes(i,kk)
                gamma2 = (wlhvocp)*(wlhv/(rwat*(t(i,kk+1)**d_two)))*          &
                     & qes(i,kk+1)
                dhh = hcd(i)
                xdt = d_half*(t(i,kk)+t(i,kk+1))
                dg = d_half*(gamma1+gamma2)
                dh = d_half*(hes(i,kk)+hes(i,kk+1))
                dz = (z(i,kk)-z(i,kk+1))*dkk(i,kk)
                aa0(i) = aa0(i) + edt(i)*dz*(gti/(cpd*xdt))*((dhh-dh)/    &
                     & (d_one+dg))
                !
                !---          modified by larger scale
                !
                gamma1 = (wlhvocp)*(wlhv/(rwat*(tn(i,kk)**d_two)))*qeso(i,kk)
                gamma2 = (wlhvocp)*(wlhv/(rwat*(tn(i,kk+1)**d_two)))*         &
                     & qeso(i,kk+1)
                dhh = hcdo(i)
                xdt = d_half*(tn(i,kk)+tn(i,kk+1))
                dg = d_half*(gamma1+gamma2)
                dh = d_half*(heso(i,kk)+heso(i,kk+1))
                dz = (zo(i,kk)-zo(i,kk+1))*dkk(i,kk)
                aa1(i) = aa1(i) + edto(i)*dz*(gti/(cpd*xdt))              &
                     & *((dhh-dh)/(d_one+dg))
                !
                !---          modified by cloud
                !
                gamma1 = (wlhvocp)*(wlhv/(rwat*(xt(i,kk)**d_two)))*xqes(i,kk)
                gamma2 = (wlhvocp)*(wlhv/(rwat*(xt(i,kk+1)**d_two)))*         &
                     & xqes(i,kk+1)
                dhh = xhcd(i)
                xdt = d_half*(xt(i,kk)+xt(i,kk+1))
                dg = d_half*(gamma1+gamma2)
                dh = d_half*(xhes(i,kk)+xhes(i,kk+1))
                dz = (xz(i,kk)-xz(i,kk+1))*dkk(i,kk)
                xaa0(i) = xaa0(i) + edtx(i)*dz*(gti/(cpd*xdt))            &
                     & *((dhh-dh)/(d_one+dg))
             end if
          end if
       end do
    end do
    !
    !---  large scale forcing
    !
    do i = istart , iend
       if ( aa0(i) >= d_zero ) then
          if ( igcc == 1 ) then
             f = (aa1(i)-aa0(i))/dtime  ! Arakawa-Schubert closure
          else if ( igcc == 2 ) then
             f = aa0(i)/dtauc2d(i,jslc)   ! Fritsch-Chappell closure
          end if
          xk = (xaa0(i)-aa0(i))/mbdt
          xmb(i) = -f/xk
          if ( f <= d_zero .or. xk >= d_zero ) xmb(i) = d_zero
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
          if ( aa0(i) >= d_zero ) then
             if ( k <= ktop(i) ) then
                outtes = dellat(i,k)*xmb(i)*secpd
                if ( (outtes > htmax2d(i,jslc)) .or.  &
                     & (outtes < htmin2d(i,jslc)) ) then
                   xmb(i) = d_zero
                   aa0(i) = -d_one
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
       if ( aa0(i) >= d_zero ) then

          if ( ktop(i) > 1 .and. kbcon(i) > 1 ) then
             kclth = ktop(i) - kbcon(i) + 1
             akclth = d_one/dble(kclth)
             do k = kbcon(i) , ktop(i)
                kk = kz - k + 1
                cldlwc(i,kk) = cllwcv
                cldfra(i,kk) = d_one - (d_one-clfrcv)**akclth
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

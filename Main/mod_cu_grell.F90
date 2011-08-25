!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   This file is part of ICTP RegCM.
!
!   ICTP RegCM is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   ICTP RegCM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_cu_grell

  use mod_runparams
  use mod_atm_interface
  use mod_cu_common
  use mod_bats
  use mod_rad
  use mod_che_trac
  use mod_service 
  use mod_memutil
 
  private

  integer , parameter :: icut = 0
  real(8) , parameter :: xacact = -0.99999D0
  real(8) , parameter :: tcrit = 50.0D0
  real(8) , parameter :: c0 = 0.002D0
  real(8) :: alsixt
  real(8) , dimension(2) :: ae , be
  real(8) , pointer , dimension(:,:) :: outq , outt , p , po ,  &
                              q , qo , t , tn , vsp
  real(8) , pointer , dimension(:,:) :: dby , dbyo , dellah ,       &
              dellaq , dellat , dkk , he , heo , hes , heso , pwc , &
              pwcd , pwcdo , pwco , qc , qco , qes , qeso , qrcd ,  &
              qrcdo , tv , tvo , xdby , xhe , xhes , xpwc , xpwcd , &
              xq , xqc , xqes , xqrcd , xt , xtv , xz , z , zo
  real(8) , pointer , dimension(:) :: pret , psur , qcrit , ter11
  integer , pointer , dimension(:) :: kdet
  integer , pointer , dimension(:) :: kmin , k22 , kb , kbcon , kds , ktop
  real(8) , pointer , dimension(:) :: xac , xao , bu , buo , edt ,  &
              edto , edtx , hcd , hcdo , hkb , hkbo , pwcav ,       &
              pwcavo , pwcev , pwcevo , qcd , qcdo , qck , qcko ,   &
              qkb , qkbo , vshear , xxac , xhcd , xhkb , xmb ,      &
              xpwcav , xpwcev , xqcd , xqck , xqkb
!
  real(8) , public , pointer , dimension(:,:) :: dtauc2d , pbcmax2d ,   &
              mincld2d , shrmax2d , shrmin2d , edtmax2d , edtmin2d ,    &
              edtmaxo2d , edtmaxx2d , edtmino2d , edtminx2d , htmax2d , &
              htmin2d
  integer , public , pointer , dimension(:,:) :: kbmax2d
!
  public :: allocate_mod_cu_grell , cuparan
!
  contains

  subroutine allocate_mod_cu_grell
    implicit none

    call getmem2d(dtauc2d,1,iy,1,jxp,'cu_grell:dtauc2d')
    call getmem2d(pbcmax2d,1,iy,1,jxp,'cu_grell:pbcmax2d')
    call getmem2d(mincld2d,1,iy,1,jxp,'cu_grell:mincld2d')
    call getmem2d(shrmax2d,1,iy,1,jxp,'cu_grell:shrmax2d')
    call getmem2d(shrmin2d,1,iy,1,jxp,'cu_grell:shrmin2d')
    call getmem2d(edtmax2d,1,iy,1,jxp,'cu_grell:edtmax2d')
    call getmem2d(edtmin2d,1,iy,1,jxp,'cu_grell:edtmin2d')
    call getmem2d(edtmaxo2d,1,iy,1,jxp,'cu_grell:edtmaxo2d')
    call getmem2d(edtmino2d,1,iy,1,jxp,'cu_grell:edtmino2d')
    call getmem2d(edtmaxx2d,1,iy,1,jxp,'cu_grell:edtmaxx2d')
    call getmem2d(edtminx2d,1,iy,1,jxp,'cu_grell:edtminx2d')
    call getmem2d(htmax2d,1,iy,1,jxp,'cu_grell:htmax2d')
    call getmem2d(htmin2d,1,iy,1,jxp,'cu_grell:htmin2d')
    call getmem2d(kbmax2d,1,iy,1,jxp,'cu_grell:kbmax2d')

    call getmem2d(outq,1,iy,1,kz,'cu_grell:outq')
    call getmem2d(outt,1,iy,1,kz,'cu_grell:outt')
    call getmem2d(p,1,iy,1,kz,'cu_grell:p')
    call getmem2d(po,1,iy,1,kz,'cu_grell:po')
    call getmem2d(q,1,iy,1,kz,'cu_grell:q')
    call getmem2d(qo,1,iy,1,kz,'cu_grell:qo')
    call getmem2d(t,1,iy,1,kz,'cu_grell:t')
    call getmem2d(tn,1,iy,1,kz,'cu_grell:tn')
    call getmem2d(vsp,1,iy,1,kz,'cu_grell:vsp')
!
    call getmem1d(pret,1,iy,'cu_grell:pret')
    call getmem1d(psur,1,iy,'cu_grell:psur')
    call getmem1d(qcrit,1,iy,'cu_grell:qcrit')
    call getmem1d(ter11,1,iy,'cu_grell:ter11')
!
    call getmem1d(kdet,1,iy,'cu_grell:kdet')
    call getmem1d(kmin,1,iy,'cu_grell:kmin')
    call getmem1d(k22,1,iy,'cu_grell:k22')
    call getmem1d(kb,1,iy,'cu_grell:kb')
    call getmem1d(kbcon,1,iy,'cu_grell:kbcon')
    call getmem1d(kds,1,iy,'cu_grell:kds')
    call getmem1d(ktop,1,iy,'cu_grell:ktop')
!
    call getmem2d(dby,1,iy,1,kz,'cu_grell:dby')
    call getmem2d(dbyo,1,iy,1,kz,'cu_grell:dbyo')
    call getmem2d(dellah,1,iy,1,kz,'cu_grell:dellah')
    call getmem2d(dellaq,1,iy,1,kz,'cu_grell:dellaq')
    call getmem2d(dellat,1,iy,1,kz,'cu_grell:dellat')
    call getmem2d(dkk,1,iy,1,kz,'cu_grell:dkk')
    call getmem2d(he,1,iy,1,kz,'cu_grell:he')
    call getmem2d(heo,1,iy,1,kz,'cu_grell:heo')
    call getmem2d(hes,1,iy,1,kz,'cu_grell:hes')
    call getmem2d(heso,1,iy,1,kz,'cu_grell:heso')
    call getmem2d(pwc,1,iy,1,kz,'cu_grell:pwc')
    call getmem2d(pwco,1,iy,1,kz,'cu_grell:pwco')
    call getmem2d(pwcd,1,iy,1,kz,'cu_grell:pwcd')
    call getmem2d(pwcdo,1,iy,1,kz,'cu_grell:pwcdo')
    call getmem2d(qc,1,iy,1,kz,'cu_grell:qc')
    call getmem2d(qco,1,iy,1,kz,'cu_grell:qco')
    call getmem2d(qes,1,iy,1,kz,'cu_grell:qes')
    call getmem2d(qeso,1,iy,1,kz,'cu_grell:qeso')
    call getmem2d(qrcd,1,iy,1,kz,'cu_grell:qrcd')
    call getmem2d(qrcdo,1,iy,1,kz,'cu_grell:qrcdo')
    call getmem2d(tv,1,iy,1,kz,'cu_grell:tv')
    call getmem2d(tvo,1,iy,1,kz,'cu_grell:tvo')
    call getmem2d(xdby,1,iy,1,kz,'cu_grell:xdby')
    call getmem2d(xhe,1,iy,1,kz,'cu_grell:xhe')
    call getmem2d(xhes,1,iy,1,kz,'cu_grell:xhes')
    call getmem2d(xpwc,1,iy,1,kz,'cu_grell:xpwc')
    call getmem2d(xpwcd,1,iy,1,kz,'cu_grell:xpwcd')
    call getmem2d(xq,1,iy,1,kz,'cu_grell:xq')
    call getmem2d(xqc,1,iy,1,kz,'cu_grell:xqc')
    call getmem2d(xqes,1,iy,1,kz,'cu_grell:xqes')
    call getmem2d(xqrcd,1,iy,1,kz,'cu_grell:xqrcd')
    call getmem2d(xt,1,iy,1,kz,'cu_grell:xt')
    call getmem2d(xtv,1,iy,1,kz,'cu_grell:xtv')
    call getmem2d(xz,1,iy,1,kz,'cu_grell:xz')
    call getmem2d(z,1,iy,1,kz,'cu_grell:z')
    call getmem2d(zo,1,iy,1,kz,'cu_grell:zo')
!
    call getmem1d(xac,1,iy,'cu_grell:xac')
    call getmem1d(xao,1,iy,'cu_grell:xao')
    call getmem1d(bu,1,iy,'cu_grell:bu')
    call getmem1d(buo,1,iy,'cu_grell:buo')
    call getmem1d(edt,1,iy,'cu_grell:edt')
    call getmem1d(edto,1,iy,'cu_grell:edto')
    call getmem1d(edtx,1,iy,'cu_grell:edtx')
    call getmem1d(hcd,1,iy,'cu_grell:hcd')
    call getmem1d(hcdo,1,iy,'cu_grell:hcdo')
    call getmem1d(hkb,1,iy,'cu_grell:hkb')
    call getmem1d(hkbo,1,iy,'cu_grell:hkbo')
    call getmem1d(pwcav,1,iy,'cu_grell:pwcav')
    call getmem1d(pwcavo,1,iy,'cu_grell:pwcavo')
    call getmem1d(pwcev,1,iy,'cu_grell:pwcev')
    call getmem1d(pwcevo,1,iy,'cu_grell:pwcevo')
    call getmem1d(qcd,1,iy,'cu_grell:qcd')
    call getmem1d(qcdo,1,iy,'cu_grell:qcdo')
    call getmem1d(qck,1,iy,'cu_grell:qck')
    call getmem1d(qcko,1,iy,'cu_grell:qcko')
    call getmem1d(qkb,1,iy,'cu_grell:qkb')
    call getmem1d(qkbo,1,iy,'cu_grell:qkbo')
    call getmem1d(vshear,1,iy,'cu_grell:vshear')
    call getmem1d(xxac,1,iy,'cu_grell:xxac')
    call getmem1d(xhcd,1,iy,'cu_grell:xhcd')
    call getmem1d(xhkb,1,iy,'cu_grell:xhkb')
    call getmem1d(xmb,1,iy,'cu_grell:xmb')
    call getmem1d(xpwcav,1,iy,'cu_grell:xpwcav')
    call getmem1d(xpwcev,1,iy,'cu_grell:xpwcev')
    call getmem1d(xqcd,1,iy,'cu_grell:xqcd')
    call getmem1d(xqck,1,iy,'cu_grell:xqck')
    call getmem1d(xqkb,1,iy,'cu_grell:xqkb')
!
    alsixt = dlog(610.71D0)
    be(1) = ep2*wlhvocp*3.50D0
    be(2) = ep2*2.834D6*rcpd*3.50D0
    ae(1) = be(1)*rtzero + alsixt
    ae(2) = be(2)*rtzero + alsixt
  end subroutine allocate_mod_cu_grell

  subroutine cuparan(j)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    implicit none
!
    integer , intent(in) :: j
!
    real(8) :: aprdiv , calc , pkdcut , pkk , prainx , us , vs
    integer :: i , jp1 , iconj , iend , istart , k , kk
!
    character (len=64) :: subroutine_name='cuparan'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!   zero out radiative clouds
!
    cldlwc = d_zero
    cldfra = d_zero

    pkdcut = 75.0D0
    istart = 2 + icut
    iend = iym2 - icut
!
!   prepare input, erase output
!
    outq(:,:) = d_zero
    outt(:,:) = d_zero
    p(:,:)    = d_zero
    po(:,:)   = d_zero
    q(:,:)    = d_zero
    qo(:,:)   = d_zero
    t(:,:)    = d_zero
    tn(:,:)   = d_zero
    vsp(:,:)  = d_zero
!
    pret(:)  = d_zero
    psur(:)  = d_zero
    qcrit(:) = d_zero
    ter11(:) = d_zero
!
    kdet(:)  = 2
    k22(:)   = 1
    ktop(:)  = 1
    kbcon(:) = 1
    kb(:)    = 1
    kds(:)   = 1
    kmin(:)  = 1
!
    dby(:,:) = d_zero
    dbyo(:,:) = d_zero
    dellah(:,:) = d_zero
    dellaq(:,:) = d_zero
    dellat(:,:) = d_zero
    dkk(:,:) = d_one
    he(:,:) = d_zero
    heo(:,:) = d_zero
    hes(:,:) = d_zero
    heso(:,:) = d_zero
    pwc(:,:) = d_zero
    pwco(:,:) = d_zero
    pwcd(:,:) = d_zero
    pwcdo(:,:) = d_zero
    qc(:,:) = d_zero
    qco(:,:) = d_zero
    qes(:,:) = d_zero
    qeso(:,:) = d_zero
    qrcd(:,:) = d_zero
    qrcdo(:,:) = d_zero
    tv(:,:) = d_zero
    tvo(:,:) = d_zero
    xdby(:,:) = d_zero
    xhe(:,:) = d_zero
    xhes(:,:) = d_zero
    xpwc(:,:) = d_zero
    xpwcd(:,:) = d_zero
    xq(:,:) = d_zero
    xqc(:,:) = d_zero
    xqes(:,:) = d_zero
    xqrcd(:,:) = d_zero
    xt(:,:) = d_zero
    xtv(:,:) = d_zero
    xz(:,:) = d_zero
    z(:,:) = d_zero
    zo(:,:) = d_zero
!
    xac(:) = d_zero
    xao(:) = d_zero
    bu(:) = d_zero
    buo(:) = d_zero
    edt(:) = d_zero
    edto(:) = d_zero
    edtx(:) = d_zero
    hcd(:) = d_zero
    hcdo(:) = d_zero
    hkb(:) = d_zero
    hkbo(:) = d_zero
    pwcav(:) = d_zero
    pwcavo(:) = d_zero
    pwcev(:) = d_zero
    pwcevo(:) = d_zero
    qcd(:) = d_zero
    qcdo(:) = d_zero
    qck(:) = d_zero
    qcko(:) = d_zero
    qkb(:) = d_zero
    qkbo(:) = d_zero
    vshear(:) = d_zero
    xxac(:) = d_zero
    xhcd(:) = d_zero
    xhkb(:) = d_zero
    xmb(:) = d_zero
    xpwcav(:) = d_zero
    xpwcev(:) = d_zero
    xqcd(:) = d_zero
    xqck(:) = d_zero
    xqkb(:) = d_zero
!
    do k = 1 , kz
      do i = istart , iend
        kk = kz - k + 1
        jp1 = j + 1
        us = (atm1%u(i,kk,j)/sps2%ps(i,j)+       &
              atm1%u(i+1,kk,j)/sps2%ps(i+1,j)+   &
              atm1%u(i,kk,jp1)/sps2%ps(i,jp1)+   &
              atm1%u(i+1,kk,jp1)/sps2%ps(i+1,jp1))*d_rfour
        vs = (atm1%v(i,kk,j)/sps2%ps(i,j)+       &
              atm1%v(i+1,kk,j)/sps2%ps(i+1,j)+   &
              atm1%v(i,kk,jp1)/sps2%ps(i,jp1)+   &
              atm1%v(i+1,kk,jp1)/sps2%ps(i+1,jp1))*d_rfour
        t(i,k) = atm2%t(i,kk,j)/sps2%ps(i,j)
        q(i,k) = atm2%qv(i,kk,j)/sps2%ps(i,j)
        if ( q(i,k) < 1.0D-08 ) q(i,k) = 1.0D-08
        tn(i,k) = t(i,k) + (aten%t(i,kk,j))/sps2%ps(i,j)*dt
        qo(i,k) = q(i,k) + (aten%qv(i,kk,j))/sps2%ps(i,j)*dt
        p(i,k) = d_10*sps2%ps(i,j)*a(kk) + d_10*ptop
        vsp(i,k) = dsqrt(us**d_two+vs**d_two)
        if ( qo(i,k) < 1.0D-08 ) qo(i,k) = 1.0D-08
!
        po(i,k) = p(i,k)
        psur(i) = d_10*sps2%ps(i,j) + d_10*ptop
        outt(i,k) = d_zero
        pkk = psur(i) - po(i,k)
        if ( pkk <= pkdcut ) kdet(i) = kdet(i) + 1
        outq(i,k) = d_zero
        ter11(i) = mddom%ht(i,j)*regrav
        if ( ter11(i) <= d_zero ) ter11(i) = 1.0D-05
        qcrit(i) = qcrit(i) + aten%qv(i,kk,j)
      end do
    end do
!
!   call cumulus parameterization
!
    call cup(istart,iend,j)
!
!   return cumulus parameterization
!
    do k = 1 , kz
      do i = istart , iend
        if ( pret(i) > d_zero ) then
          kk = kz - k + 1
          aten%t(i,kk,j) = sps2%ps(i,j)*outt(i,k) + aten%t(i,kk,j)
          aten%qv(i,kk,j) = sps2%ps(i,j)*outq(i,k) + aten%qv(i,kk,j)
        end if
      end do
    end do
!
!   rain in cm.
!
    calc = d_half
    iconj = 0
    do i = istart , iend
      prainx = pret(i)*calc*dt
      if ( prainx > dlowval ) then
        sfsta%rainc(i,j) = sfsta%rainc(i,j) + prainx
!       precipitation rate for bats (mm/s)
        aprdiv = dble(ntsrf)
        if ( ktau == 0 ) aprdiv = d_one
        pptc(i,j) = pptc(i,j) + prainx/dtsec/aprdiv
        iconj = iconj + 1
      end if
    end do
    icon(j) = iconj
!
    call time_end(subroutine_name,idindx)
  end subroutine cuparan
!
!   GRELL CUMULUS SCHEME
!
  subroutine cup(istart,iend,j)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    implicit none
!
    integer :: iend , istart , j
    intent (in) j , istart , iend
!
    real(8) :: adw , akclth , aup , detdo , detdoq , dg , dh ,   &
               dhh , dp_s , dq , xdt , dv1 , dv1q , dv2 , dv2q , &
               dv3 , dv3q , dz , dz1 , dz2 , dzo , e , eo , f ,  &
               agamma , agamma0 , agamma1 , agamma2 , agammo ,   &
               agammo0 , mbdt , outtes , pbcdif , qrch , qrcho , &
               tvbar , tvbaro , xk
    integer :: i , iph , ipho , k , kbcono , kclth , kk , lpt
!
    character (len=64) :: subroutine_name='cup'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

    mbdt = dt*5.0D-03
    f  = -d_one
    xk = -d_one
!
!   environmental conditions, first heights
!
    do k = 1 , kz
      do i = istart , iend
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

    do i = istart , iend
      if ( qcrit(i) <= d_zero ) then
        xac(i) = -d_one
      end if
      if ( icup /= 2 ) then
        if (cucontrol(i,j) /= 2 ) then
          xac(i) = -d_one
        end if
      end if
      z(i,1)  = ter11(i) - (dlog(p(i,1))-dlog(psur(i)))*rgas*tv(i,1)*regrav
      zo(i,1) = ter11(i) - (dlog(po(i,1))-dlog(psur(i)))*rgas*tvo(i,1)*regrav
    end do

    do k = 2 , kz
      do i = istart , iend
        tvbar = d_half*(tv(i,k)+tv(i,k-1))
        z(i,k) = z(i,k-1)-(dlog(p(i,k))-dlog(p(i,k-1)))*rgas*tvbar*regrav
        tvbaro = d_half*(tvo(i,k)+tvo(i,k-1))
        zo(i,k) = zo(i,k-1)-(dlog(po(i,k))-dlog(po(i,k-1)))*rgas*tvbaro*regrav
      end do
    end do
!
!   moist static energy
!
    do k = 1 , kz
      do i = istart , iend
        he(i,k) = egrav*z(i,k) + cpd*t(i,k) + wlhv*q(i,k)
        hes(i,k) = egrav*z(i,k) + cpd*t(i,k) + wlhv*qes(i,k)
        if ( he(i,k) >= hes(i,k) ) he(i,k) = hes(i,k)
        heo(i,k) = egrav*zo(i,k) + cpd*tn(i,k) + wlhv*qo(i,k)
        heso(i,k) = egrav*zo(i,k) + cpd*tn(i,k) + wlhv*qeso(i,k)
        if ( heo(i,k) >= heso(i,k) ) heo(i,k) = heso(i,k)
        xt(i,k) = t(i,k)
        xq(i,k) = q(i,k)
        xhe(i,k) = he(i,k)
        if ( k /= kz ) qrcd(i,k) = d_half*(qes(i,k)+qes(i,k+1))
        if ( k /= kz ) qrcdo(i,k) = d_half*(qeso(i,k)+qeso(i,k+1))
      end do
    end do
!
!   determine level with highest moist static energy content.
!
    call maximi(he,iy,kz,1,kbmax2d(i,j),k22,istart,iend)
    do i = istart , iend
      if ( xac(i) >= d_zero ) then
        if ( k22(i) >= kbmax2d(i,j) ) then
          xac(i) = -d_one
          cycle
        end if
        hkb(i) = he(i,k22(i))
        qkb(i) = q(i,k22(i))
        hkbo(i) = heo(i,k22(i))
        qkbo(i) = qo(i,k22(i))
        qck(i) = qkb(i)
        qcko(i) = qkbo(i)
      end if
    end do
!
!   decide for convective cloud base
!
    do i = istart , iend
      if ( xac(i) >= d_zero ) then
        do k = 1 , kdet(i)
          kk = kdet(i) - k + 1
!         dkk(i,kk) = 0.75D0*dkk(i,kk+1)
          dkk(i,k) = d_one - dble(kk)/dble(kdet(i))
        end do

120     continue

        kb(i) = k22(i)
        kbcon(i) = kb(i)

140     continue

        dh = d_half*(hes(i,kbcon(i))+hes(i,kbcon(i)+1))
        if ( hkb(i) < dh ) then
          kbcon(i) = kbcon(i) + 1
          if ( kbcon(i) > kbmax2d(i,j) ) then
            xac(i) = -d_one
            cycle
          end if
          go to 140
        else
!
!         after large-scale forcing is applied, possible lid should be removed!!!
!
          kbcono = kb(i)
!ictp
150       continue
          if ( kbcono > kbmax2d(i,j) ) then
            xac(i) = -d_one
            cycle
          end if
!ictp_
          dh = d_half*(heso(i,kbcono)+heso(i,kbcono+1))
          if ( hkbo(i) < dh ) then
            kbcono = kbcono + 1
            go to 150
          else
            pbcdif = -p(i,kbcono) + p(i,kb(i))
!           below was commented out
!           as uncommenting the following lines for experiment 2/5/95
            if ( pbcdif > pbcmax2d(i,j) ) then
!             this is where typo was (pbdcdif)
              k22(i) = k22(i) + 1
              if ( k22(i) >= kbmax2d(i,j) ) then
                xac(i) = -d_one
                cycle
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
      end if
    end do
!
!   downdraft originating level
!
    call minimi(he,iy,kz,kb,kz,kmin,istart,iend)
    call maximi(vsp,iy,kz,1,kz,kds,istart,iend)
!
!   static control
!
!   determine cloud top
!
    do i = istart , iend
      if ( xac(i) >= 0 ) then
        if ( kmin(i) <= 3 ) then
          xac(i) = -d_one
          cycle
        end if
        if ( kds(i) >= kz ) kds(i) = kz - 1
        if ( kds(i) <= kbcon(i) ) kds(i) = kbcon(i)
        dby(i,kz) = hkb(i) - hes(i,kz)
        dbyo(i,kz) = hkbo(i) - heso(i,kz)
      end if
    end do
    do k = 1 , kz - 1
      do i = istart , iend
        if ( xac(i) > xacact ) then
          dby(i,k) = hkb(i) - d_half*(hes(i,k)+hes(i,k+1))
          dbyo(i,k) = hkbo(i) - d_half*(heso(i,k)+heso(i,k+1))
        end if
      end do
    end do
    do i = istart , iend
      if ( xac(i) > xacact ) then
        do k = 2 , kz - kbcon(i) - 1
          kk = kz - k + 1
          if ( dby(i,kk) >= d_zero ) then
            ktop(i) = kk + 1
            go to 320
          end if
        end do
        xac(i) = -d_one
        go to 400
320     continue
        if ( ktop(i) > kz ) ktop(i) = kz
        if ( p(i,kbcon(i))-p(i,ktop(i)) < mincld2d(i,j) ) xac(i) = -d_one
      end if
400   continue
    end do
!
!   moisture and cloud work functions
!
    do k = 2 , kz - 1
      do i = istart , iend
        if ( xac(i) > xacact ) then
          if ( k > kbcon(i) ) then
            if ( k < ktop(i) ) then
              dz = d_half*(z(i,k+1)-z(i,k-1))
              dz1 = z(i,k) - z(i,k-1)
              agamma = (wlhvocp)*(wlhv/(rwat*(t(i,k)**d_two)))*qes(i,k)
              agamma0 = (wlhvocp)*(wlhv/(rwat*(t(i,k-1)**d_two)))*qes(i,k-1)
              qrch = qes(i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*dby(i,k)
              qc(i,k) = (qck(i)-qrch)/(d_one+c0*dz) + qrch
              pwc(i,k) = c0*dz*(qc(i,k)-qrch)
              qck(i) = qc(i,k)
              pwcav(i) = pwcav(i) + pwc(i,k)
              dz1 = z(i,k) - z(i,k-1)
              xac(i) = xac(i)+dz1*(egrav/(cpd*(d_half*(t(i,k)+t(i,k-1))))) * &
                       dby(i,k-1)/(d_one+d_half*(agamma+agamma0))
              dzo = d_half*(zo(i,k+1)-zo(i,k-1))
              dz2 = zo(i,k) - zo(i,k-1)
              agammo = (wlhvocp)*(wlhv/(rwat*(tn(i,k)**d_two)))*qeso(i,k)
              agammo0 = (wlhvocp)*(wlhv/(rwat*(tn(i,k-1)**d_two)))*qeso(i,k-1)
              qrcho = qeso(i,k)+(d_one/wlhv)*(agammo/(d_one+agammo))*dbyo(i,k)
              qco(i,k) = (qcko(i)-qrcho)/(d_one+c0*dzo) + qrcho
              pwco(i,k) = c0*dzo*(qco(i,k)-qrcho)
              qcko(i) = qco(i,k)
              pwcavo(i) = pwcavo(i) + pwco(i,k)
              xao(i) = xao(i) + dz2*(egrav/(cpd*((tn(i,k)+tn(i,k-1))*d_half))) * &
                       dbyo(i,k-1)/(d_one+d_half*(agammo+agammo0))
            end if
          end if
        end if
      end do
    end do
!
    do i = istart , iend
      if ( xac(i) > xacact ) then
        k = ktop(i)
        dz = d_half*(z(i,k)-z(i,k-1))
        agamma = (wlhvocp)*(wlhv/(rwat*(t(i,k)**d_two)))*qes(i,k)
        qrch = qes(i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*dby(i,k)
        qc(i,k) = qes(i,k)
        pwc(i,k) = (qrch-qes(i,k))
        pwcav(i) = pwcav(i) + pwc(i,k)
!
        dz = d_half*(zo(i,k)-zo(i,k-1))
        agamma = (wlhvocp)*(wlhv/(rwat*(tn(i,k)**d_two)))*qeso(i,k)
        qrcho = qeso(i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*dbyo(i,k)
        qco(i,k) = qeso(i,k)
        pwco(i,k) = (qrcho-qeso(i,k))
        pwcavo(i) = pwcavo(i) + pwco(i,k)
      end if
    end do
!
!   downdraft calculations
!
!   determine downdraft strength in terms of windshear
!
    do kk = 1 , kz/2
      do i = istart , iend
        if ( xac(i) > xacact ) then
          vshear(i) = vshear(i) + dabs((vsp(i,kk+1)-vsp(i,kk))/(z(i,kk+1)-z(i,kk)))
        end if
      end do
    end do
    do i = istart , iend
      if ( xac(i) > xacact ) then
        vshear(i) = d_1000*vshear(i)/dble(kz/2)
        edt(i) = d_one - (1.591D0-0.639D0*vshear(i)+0.0953D0*(vshear(i)**d_two) - &
                          0.00496D0*(vshear(i)**d_three))

        if ( edt(i) > shrmax2d(i,j) ) edt(i) = shrmax2d(i,j)
        if ( edt(i) < shrmin2d(i,j) ) edt(i) = shrmin2d(i,j)

        edto(i) = edt(i)
        edtx(i) = edt(i)
        qrcd(i,kz) = qes(i,kz)
        hcd(i) = d_half*(he(i,kmin(i))+he(i,kmin(i)+1))
        qcd(i) = d_half*(q(i,kmin(i))+q(i,kmin(i)+1))
        qrcdo(i,kz) = qeso(i,kz)
        hcdo(i) = heso(i,kz)
        hcdo(i) = d_half*(heo(i,kmin(i))+heo(i,kmin(i)+1))
        qcdo(i) = d_half*(qo(i,kmin(i))+qo(i,kmin(i)+1))
        bu(i) = d_zero
        buo(i) = d_zero
      end if
    end do
    do k = 1 , kz - 1
      do i = istart , iend
        if ( xac(i) > xacact ) then
          if ( k < kmin(i) ) then
            kk = kmin(i) - k
            dz = -(z(i,kk)-z(i,kk+2))*d_half
            bu(i) = bu(i) + dz*(hcd(i)-d_half*(hes(i,kk)+hes(i,kk+1)))
            dq = (qes(i,kk)+qes(i,kk+1))*d_half
            xdt = (t(i,kk)+t(i,kk+1))*d_half
            agamma = (wlhvocp)*(wlhv/(rwat*(xdt**d_two)))*dq
            dh = hcd(i) - d_half*(hes(i,kk)+hes(i,kk+1))
            qrcd(i,kk) = (dq+(d_one/wlhv)*(agamma/(d_one+agamma))*dh)
            pwcd(i,kk) = dkk(i,kk)*(qcd(i)-qrcd(i,kk))
            qcd(i) = qrcd(i,kk)
            pwcev(i) = pwcev(i) + pwcd(i,kk)
!
            dz = d_half*(zo(i,kk+2)-zo(i,kk))
            buo(i) = buo(i) + dz*(hcdo(i)-(heso(i,kk)+heso(i,kk+1))*d_half)
            dq = (qeso(i,kk)+qeso(i,kk+1))*d_half
            xdt = (tn(i,kk)+tn(i,kk+1))*d_half
            agamma = (wlhvocp)*(wlhv/(rwat*(xdt**d_two)))*dq
            dh = hcdo(i) - d_half*(heso(i,kk)+heso(i,kk+1))
            qrcdo(i,kk) = (dq+(d_one/wlhv)*(agamma/(d_one+agamma))*dh)
            pwcdo(i,kk) = dkk(i,kk)*(qcdo(i)-qrcdo(i,kk))
            qcdo(i) = qrcdo(i,kk)
            pwcevo(i) = pwcevo(i) + pwcdo(i,kk)
          end if
        end if
      end do
    end do
!
    do i = istart , iend
      if ( xac(i) > xacact ) then
        if ( bu(i) >= d_zero .or. buo(i) >= d_zero .or. &
          pwcev(i) >= d_zero .or. pwcevo(i) >= d_zero ) then
          xac(i) = -d_one
        end if
        edt(i) = -edt(i)*pwcav(i)/pwcev(i)
        if ( edt(i) > edtmax2d(i,j) ) edt(i) = edtmax2d(i,j)
        if ( edt(i) < edtmin2d(i,j) ) edt(i) = edtmin2d(i,j)
        edto(i) = -edto(i)*pwcavo(i)/pwcevo(i)
        if ( edto(i) > edtmaxo2d(i,j) ) edto(i) = edtmaxo2d(i,j)
        if ( edto(i) < edtmino2d(i,j) ) edto(i) = edtmino2d(i,j)
      end if
    end do
!
!   what would the change be?
!
    do i = istart , iend
      if ( xac(i) > xacact ) then
        k = 1
        dz = d_half*(z(i,2)-z(i,1))
        dp_s = 50.0D0*(psur(i)-p(i,2))
        dellah(i,1) = edt(i)*(dkk(i,1)*hcd(i)-dkk(i,1) * &
                      d_half*(he(i,1)+he(i,2)))*egrav/dp_s
        dellaq(i,1) = edt(i)*(dkk(i,1)*qrcd(i,1)-dkk(i,1) * &
                      d_half*(q(i,1)+q(i,2)))*egrav/dp_s
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
        if ( xac(i) > xacact ) then
          if ( k /= 1 .and. k < ktop(i) ) then
            dv1 = d_half*(he(i,k)+he(i,k+1))
            dv2 = he(i,k)
            dv3 = d_half*(he(i,k)+he(i,k-1))
            dv1q = d_half*(q(i,k)+q(i,k+1))
            dv2q = q(i,k)
            dv3q = d_half*(q(i,k)+q(i,k-1))
!
!           specifiy detrainment of downdraft, has to be consistent
!           with zd calculations in soundd.
!
            detdo = (d_one-dkk(i,k))*(hcd(i)-dv2)
            detdoq = (d_one-dkk(i,k))*(qrcd(i,k)-dv2q)
            dz = d_half*(z(i,k+1)-z(i,k-1))
!
!            changed due to subsidence and entrainment
!
            aup = d_one
            if ( k <= k22(i) ) aup = d_zero
            adw = d_one
            if ( k > kmin(i) ) adw = d_zero
            dp_s = +50.0D0*(p(i,k-1)-p(i,k+1))
            dellah(i,k) = ((aup-adw*edt(i))*(dv1-dv2)+(aup-adw*edt(i)) * &
                  (dv2-dv3))*egrav/dp_s + adw*edt(i)*detdo*egrav/dp_s
            dellaq(i,k) = ((aup-adw*edt(i))*(dv1q-dv2q)+(aup-adw*edt(i)) * &
               (dv2q-dv3q))*egrav/dp_s + adw*edt(i)*detdoq*egrav/dp_s
            xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
            xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
            dellat(i,k) = rcpd*(dellah(i,k)-wlhv*dellaq(i,k))
            xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-wlhv*dellaq(i,k)) + t(i,k)
            if ( xq(i,k) <= d_zero ) xq(i,k) = 1.0D-08
          end if
        end if
      end do
    end do
!
!   cloud top
!
    do i = istart , iend
      if ( xac(i) > xacact ) then
        lpt = ktop(i)
        dp_s = d_100*(p(i,lpt-1)-p(i,lpt))
        dv1 = d_half*(he(i,lpt)+he(i,lpt-1))
        dellah(i,lpt) = (hkb(i)-dv1)*egrav/dp_s
        dv1 = d_half*(q(i,lpt)+q(i,lpt-1))
        dellaq(i,lpt) = (qes(i,lpt)-dv1)*egrav/dp_s
        k = lpt
        xhe(i,k) = dellah(i,k)*mbdt + he(i,k)
        xq(i,k) = dellaq(i,k)*mbdt + q(i,k)
        dellat(i,k) = rcpd*(dellah(i,k)-wlhv*dellaq(i,k))
        xt(i,k) = (mbdt*rcpd)*(dellah(i,k)-wlhv*dellaq(i,k)) + t(i,k)
        if ( xq(i,k) <= d_zero ) xq(i,k) = 1.0D-08
        xhkb(i) = dellah(i,kbcon(i))*mbdt + hkb(i)
        xqkb(i) = dellaq(i,kbcon(i))*mbdt + qkb(i)
        if ( xqkb(i) <= d_zero ) xqkb(i) = 1.0D-08
      end if
    end do
!
!   environmental conditions, first heights
!
    do k = 1 , kz
      do i = istart , iend
        if ( xac(i) > xacact ) then
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
!   bug fix
    do k = 1 , kz - 1
      do i = istart , iend
        if ( xac(i) > xacact ) xqrcd(i,k) = d_half*(xqes(i,k)+xqes(i,k+1))
      end do
    end do
!
    do i = istart , iend
      if ( xac(i) > xacact ) then
        xz(i,1) = ter11(i) - (dlog(p(i,1))-dlog(psur(i)))*rgas*xtv(i,1)*regrav
       end if
    end do
    do k = 2 , kz
      do i = istart , iend
        if ( xac(i) > xacact ) then
          tvbar = d_half*(xtv(i,k)+xtv(i,k-1))
          xz(i,k) = xz(i,k-1) - (dlog(p(i,k))-dlog(p(i,k-1)))*rgas*tvbar*regrav
        end if
      end do
    end do
!
!   moist static energy
!
    do k = 1 , kz
      do i = istart , iend
        if ( xac(i) > xacact ) then
          xhes(i,k) = egrav*xz(i,k) + cpd*xt(i,k) + wlhv*xqes(i,k)
          if ( xhe(i,k) >= xhes(i,k) ) xhe(i,k) = xhes(i,k)
        end if
      end do
    end do
!
!   static control
!
    do i = istart , iend
      if ( xac(i) > xacact ) then
        xqck(i) = xqkb(i)
        xdby(i,kz) = xhkb(i) - xhes(i,kz)
      end if
    end do
!
!   moisture and cloud work functions
!
    do k = 1 , kz - 1
      do i = istart , iend
        if ( xac(i) >= d_zero ) then
          xdby(i,k) = xhkb(i) - d_half*(xhes(i,k)+xhes(i,k+1))
          if ( k > kbcon(i) .and. k < ktop(i) ) then
            dz = d_half*(xz(i,k+1)-xz(i,k-1))
            dz1 = xz(i,k) - xz(i,k-1)
            agamma = (wlhvocp)*(wlhv/(rwat*(xt(i,k)**d_two)))*xqes(i,k)
            agamma0 = (wlhvocp)*(wlhv/(rwat*(xt(i,k-1)**d_two)))*xqes(i,k-1)
            qrch = xqes(i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*xdby(i,k)
            xqc(i,k) = (xqck(i)-qrch)/(d_one+c0*dz) + qrch
            xpwc(i,k) = c0*dz*(xqc(i,k)-qrch)
            xqck(i) = xqc(i,k)
            xpwcav(i) = xpwcav(i) + xpwc(i,k)
            xxac(i) = xxac(i) + dz1*(egrav/(cpd*(d_half*(xt(i,k)+xt(i,k-1))))) * &
                      xdby(i,k-1)/(d_one+d_half*(agamma+agamma0))
          end if
        end if
      end do
    end do
    do i = istart , iend
      if ( xac(i) >= d_zero ) then
        k = ktop(i)
        dz = d_half*(xz(i,k)-xz(i,k-1))
        agamma = (wlhvocp)*(wlhv/(rwat*(xt(i,k)**d_two)))*xqes(i,k)
        qrch = xqes(i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*xdby(i,k)
        xqc(i,k) = xqes(i,k)
        xpwc(i,k) = (qrch-xqes(i,k))
        xpwcav(i) = xpwcav(i) + xpwc(i,k)
        xqrcd(i,kz) = xqes(i,kz)
        xhcd(i) = d_half*(xhe(i,kmin(i))+xhe(i,kmin(i)+1))
        xqcd(i) = d_half*(xq(i,kmin(i))+xq(i,kmin(i)+1))
        xpwcev(i) = d_zero
        bu(i) = d_zero
      end if
    end do
!
!   downdraft calculations
!
!   downdraft moisture properties
!
    do k = 1 , kz - 1
      do i = istart , iend
        if ( xac(i) >= d_zero ) then
          if ( k < kmin(i) ) then
            kk = kmin(i) - k
            dz = -(xz(i,kk)-xz(i,kk+2))*d_half
            bu(i) = bu(i)+dz*(xhcd(i)-d_half*(xhes(i,kk)+xhes(i,kk+1)))
            dq = d_half*(xqes(i,kk)+xqes(i,kk+1))
            xdt = d_half*(xt(i,kk)+xt(i,kk+1))
            agamma = (wlhvocp)*(wlhv/(rwat*(xdt**d_two)))*dq
            dh = xhcd(i) - d_half*(xhes(i,kk)+xhes(i,kk+1))
            xqrcd(i,kk) = (dq+(d_one/wlhv)*(agamma/(d_one+agamma))*dh)
            xpwcd(i,kk) = dkk(i,kk)*(xqcd(i)-xqrcd(i,kk))
            xqcd(i) = xqrcd(i,kk)
            xpwcev(i) = xpwcev(i) + xpwcd(i,kk)
          end if
        end if
      end do
    end do
    do i = istart , iend
      if ( xac(i) >= d_zero ) then
        if ( bu(i) >= d_zero ) then
          xac(i) = -d_one
          cycle
        end if
        if ( dabs(xpwcev(i)) > dlowval ) then
          edtx(i) = -edtx(i)*xpwcav(i)/xpwcev(i)
        end if
        if ( edtx(i) > edtmaxx2d(i,j) ) then
          edtx(i) = edtmaxx2d(i,j)
        end if
        if ( edtx(i) < edtminx2d(i,j) ) then
          edtx(i) = edtminx2d(i,j)
        end if
      end if
    end do
!
!   downdraft cloudwork functions
!
    do k = 1 , kz - 1
      do i = istart , iend
        if ( xac(i) >= d_zero ) then
          if ( k < kmin(i) ) then
            kk = kmin(i) - k
!
!           original
!
            agamma1 = (wlhvocp)*(wlhv/(rwat*(t(i,kk)**d_two)))*qes(i,kk)
            agamma2 = (wlhvocp)*(wlhv/(rwat*(t(i,kk+1)**d_two)))*qes(i,kk+1)
            dhh = hcd(i)
            xdt = d_half*(t(i,kk)+t(i,kk+1))
            dg = d_half*(agamma1+agamma2)
            dh = d_half*(hes(i,kk)+hes(i,kk+1))
            dz = (z(i,kk)-z(i,kk+1))*dkk(i,kk)
            xac(i) = xac(i) + edt(i)*dz*(egrav/(cpd*xdt))*((dhh-dh)/(d_one+dg))
!
!           modified by larger scale
!
            agamma1 = (wlhvocp)*(wlhv/(rwat*(tn(i,kk)**d_two)))*qeso(i,kk)
            agamma2 = (wlhvocp)*(wlhv/(rwat*(tn(i,kk+1)**d_two)))*qeso(i,kk+1)
            dhh = hcdo(i)
            xdt = d_half*(tn(i,kk)+tn(i,kk+1))
            dg = d_half*(agamma1+agamma2)
            dh = d_half*(heso(i,kk)+heso(i,kk+1))
            dz = (zo(i,kk)-zo(i,kk+1))*dkk(i,kk)
            xao(i) = xao(i) + edto(i)*dz*(egrav/(cpd*xdt))*((dhh-dh)/(d_one+dg))
!
!           modified by cloud
!
            agamma1 = (wlhvocp)*(wlhv/(rwat*(xt(i,kk)**d_two)))*xqes(i,kk)
            agamma2 = (wlhvocp)*(wlhv/(rwat*(xt(i,kk+1)**d_two)))*xqes(i,kk+1)
            dhh = xhcd(i)
            xdt = d_half*(xt(i,kk)+xt(i,kk+1))
            dg = d_half*(agamma1+agamma2)
            dh = d_half*(xhes(i,kk)+xhes(i,kk+1))
            dz = (xz(i,kk)-xz(i,kk+1))*dkk(i,kk)
            xxac(i) = xxac(i) + edtx(i)*dz*(egrav/(cpd*xdt))*((dhh-dh)/(d_one+dg))
          end if
        end if
      end do
    end do
!
!   large scale forcing
!
    do i = istart , iend
      if ( xac(i) >= d_zero ) then
        if ( igcc == 1 ) then
          f = (xao(i)-xac(i))/dt   ! Arakawa-Schubert closure
        else if ( igcc == 2 ) then
          f = xac(i)/dtauc2d(i,j)  ! Fritsch-Chappell closure
        end if
        xk = (xxac(i)-xac(i))/mbdt
        xmb(i) = -f/xk
        if ( f <= d_zero .or. xk >= d_zero ) xmb(i) = d_zero
      end if
    end do
!
    mflx(i,1) = xmb(i)
    mflx(i,2) = xmb(i)*edt(i)
!
!   feedback
!
    do k = 1 , kz
      do i = istart , iend
        if ( xac(i) >= d_zero ) then
          if ( k <= ktop(i) ) then
            outtes = dellat(i,k)*xmb(i)*secpd
            if ( (outtes > htmax2d(i,j)) .or. (outtes < htmin2d(i,j)) ) then
              xmb(i) = d_zero
              xac(i) = -d_one
            else
              outt(i,k) = outt(i,k) + dellat(i,k)*xmb(i)
              outq(i,k) = outq(i,k) + dellaq(i,k)*xmb(i)
              pret(i) = pret(i) + (pwc(i,k)+edt(i)*pwcd(i,k))*xmb(i)
            end if
          end if
        end if
      end do
    end do
!
!   calculate cloud fraction and water content
!
    do i = istart , iend
!
      icumtop(i,j) = 0
      icumbot(i,j) = 0
      icumdwd(i,j) = 0
!
      if ( xac(i) >= d_zero ) then

        if ( ktop(i) > 1 .and. kbcon(i) > 1 ) then
          kclth = ktop(i) - kbcon(i) + 1
          akclth = d_one/dble(kclth)
          do k = kbcon(i) , ktop(i)
            kk = kz - k + 1
            cldlwc(i,kk) = cllwcv
            cldfra(i,kk) = d_one - (d_one-clfrcv)**akclth
          end do
!
!         define convection  base and top for tracers
!
          if ( ichem == 1 ) then
            if ( ktop(i) > 1 .and. k22(i) >= 1 ) then
              icumtop(i,j) = kzp1 - ktop(i)
              icumbot(i,j) = kzp1 - k22(i)
              icumdwd(i,j) = kzp1 - kmin(i)
            end if
          end if
        end if
      end if

    end do
    call time_end(subroutine_name,idindx)

    contains
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
       character (len=64) :: subroutine_name='minimi'
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
     subroutine maximi(array,iy,kz,ks,ke,imax,istart,iend)
!
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
      character (len=64) :: subroutine_name='maximi'
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

  end subroutine cup

end module mod_cu_grell

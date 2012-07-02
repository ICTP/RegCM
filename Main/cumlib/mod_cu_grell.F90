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

  use mod_dynparam
  use mod_memutil
  use mod_service 
  use mod_cu_common
 
  private

  real(dp) , parameter :: xacact = -0.99999D0
  real(dp) , parameter :: tcrit = 50.0D0
  real(dp) , parameter :: c0 = 0.002D0
  real(dp) :: alsixt
  real(dp) , dimension(2) :: ae , be
  real(dp) , pointer , dimension(:,:,:) :: outq , outt , p , po ,  &
                              q , qo , t , tn , vsp
  real(dp) , pointer , dimension(:,:,:) :: dby , dbyo , dellah ,       &
              dellaq , dellat , dkk , he , heo , hes , heso , pwc , &
              pwcd , pwcdo , pwco , qc , qco , qes , qeso , qrcd ,  &
              qrcdo , tv , tvo , xdby , xhe , xhes , xpwc , xpwcd , &
              xq , xqc , xqes , xqrcd , xt , xtv , xz , z , zo
  real(dp) , pointer , dimension(:,:) :: pret , psur , qcrit , ter11
  integer , pointer , dimension(:,:) :: kdet
  integer , pointer , dimension(:,:) :: kmin , k22 , kb , kbcon , kds , ktop
  real(dp) , pointer , dimension(:,:) :: xac , xao , bu , buo , edt ,  &
              edto , edtx , hcd , hcdo , hkb , hkbo , pwcav ,       &
              pwcavo , pwcev , pwcevo , qcd , qcdo , qck , qcko ,   &
              qkb , qkbo , vshear , xxac , xhcd , xhkb , xmb ,      &
              xpwcav , xpwcev , xqcd , xqck , xqkb
!
  real(dp) , public , pointer , dimension(:,:,:) :: mflx
!
  real(dp) , public , pointer , dimension(:,:) :: dtauc2d , pbcmax2d ,   &
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

    call getmem3d(mflx,jci1,jci2,ici1,ici2,1,2,'cu_grell:mflx')

    call getmem2d(dtauc2d,jci1,jci2,ici1,ici2,'cu_grell:dtauc2d')
    call getmem2d(pbcmax2d,jci1,jci2,ici1,ici2,'cu_grell:pbcmax2d')
    call getmem2d(mincld2d,jci1,jci2,ici1,ici2,'cu_grell:mincld2d')
    call getmem2d(shrmax2d,jci1,jci2,ici1,ici2,'cu_grell:shrmax2d')
    call getmem2d(shrmin2d,jci1,jci2,ici1,ici2,'cu_grell:shrmin2d')
    call getmem2d(edtmax2d,jci1,jci2,ici1,ici2,'cu_grell:edtmax2d')
    call getmem2d(edtmin2d,jci1,jci2,ici1,ici2,'cu_grell:edtmin2d')
    call getmem2d(edtmaxo2d,jci1,jci2,ici1,ici2,'cu_grell:edtmaxo2d')
    call getmem2d(edtmino2d,jci1,jci2,ici1,ici2,'cu_grell:edtmino2d')
    call getmem2d(edtmaxx2d,jci1,jci2,ici1,ici2,'cu_grell:edtmaxx2d')
    call getmem2d(edtminx2d,jci1,jci2,ici1,ici2,'cu_grell:edtminx2d')
    call getmem2d(htmax2d,jci1,jci2,ici1,ici2,'cu_grell:htmax2d')
    call getmem2d(htmin2d,jci1,jci2,ici1,ici2,'cu_grell:htmin2d')
    call getmem2d(kbmax2d,jci1,jci2,ici1,ici2,'cu_grell:kbmax2d')

    call getmem3d(outq,jci1,jci2,ici1,ici2,1,kz,'cu_grell:outq')
    call getmem3d(outt,jci1,jci2,ici1,ici2,1,kz,'cu_grell:outt')
    call getmem3d(p,jci1,jci2,ici1,ici2,1,kz,'cu_grell:p')
    call getmem3d(po,jci1,jci2,ici1,ici2,1,kz,'cu_grell:po')
    call getmem3d(q,jci1,jci2,ici1,ici2,1,kz,'cu_grell:q')
    call getmem3d(qo,jci1,jci2,ici1,ici2,1,kz,'cu_grell:qo')
    call getmem3d(t,jci1,jci2,ici1,ici2,1,kz,'cu_grell:t')
    call getmem3d(tn,jci1,jci2,ici1,ici2,1,kz,'cu_grell:tn')
    call getmem3d(vsp,jci1,jci2,ici1,ici2,1,kz,'cu_grell:vsp')
!
    call getmem2d(pret,jci1,jci2,ici1,ici2,'cu_grell:pret')
    call getmem2d(psur,jci1,jci2,ici1,ici2,'cu_grell:psur')
    call getmem2d(qcrit,jci1,jci2,ici1,ici2,'cu_grell:qcrit')
    call getmem2d(ter11,jci1,jci2,ici1,ici2,'cu_grell:ter11')
!
    call getmem2d(kdet,jci1,jci2,ici1,ici2,'cu_grell:kdet')
    call getmem2d(kmin,jci1,jci2,ici1,ici2,'cu_grell:kmin')
    call getmem2d(k22,jci1,jci2,ici1,ici2,'cu_grell:k22')
    call getmem2d(kb,jci1,jci2,ici1,ici2,'cu_grell:kb')
    call getmem2d(kbcon,jci1,jci2,ici1,ici2,'cu_grell:kbcon')
    call getmem2d(kds,jci1,jci2,ici1,ici2,'cu_grell:kds')
    call getmem2d(ktop,jci1,jci2,ici1,ici2,'cu_grell:ktop')
!
    call getmem3d(dby,jci1,jci2,ici1,ici2,1,kz,'cu_grell:dby')
    call getmem3d(dbyo,jci1,jci2,ici1,ici2,1,kz,'cu_grell:dbyo')
    call getmem3d(dellah,jci1,jci2,ici1,ici2,1,kz,'cu_grell:dellah')
    call getmem3d(dellaq,jci1,jci2,ici1,ici2,1,kz,'cu_grell:dellaq')
    call getmem3d(dellat,jci1,jci2,ici1,ici2,1,kz,'cu_grell:dellat')
    call getmem3d(dkk,jci1,jci2,ici1,ici2,1,kz,'cu_grell:dkk')
    call getmem3d(he,jci1,jci2,ici1,ici2,1,kz,'cu_grell:he')
    call getmem3d(heo,jci1,jci2,ici1,ici2,1,kz,'cu_grell:heo')
    call getmem3d(hes,jci1,jci2,ici1,ici2,1,kz,'cu_grell:hes')
    call getmem3d(heso,jci1,jci2,ici1,ici2,1,kz,'cu_grell:heso')
    call getmem3d(pwc,jci1,jci2,ici1,ici2,1,kz,'cu_grell:pwc')
    call getmem3d(pwco,jci1,jci2,ici1,ici2,1,kz,'cu_grell:pwco')
    call getmem3d(pwcd,jci1,jci2,ici1,ici2,1,kz,'cu_grell:pwcd')
    call getmem3d(pwcdo,jci1,jci2,ici1,ici2,1,kz,'cu_grell:pwcdo')
    call getmem3d(qc,jci1,jci2,ici1,ici2,1,kz,'cu_grell:qc')
    call getmem3d(qco,jci1,jci2,ici1,ici2,1,kz,'cu_grell:qco')
    call getmem3d(qes,jci1,jci2,ici1,ici2,1,kz,'cu_grell:qes')
    call getmem3d(qeso,jci1,jci2,ici1,ici2,1,kz,'cu_grell:qeso')
    call getmem3d(qrcd,jci1,jci2,ici1,ici2,1,kz,'cu_grell:qrcd')
    call getmem3d(qrcdo,jci1,jci2,ici1,ici2,1,kz,'cu_grell:qrcdo')
    call getmem3d(tv,jci1,jci2,ici1,ici2,1,kz,'cu_grell:tv')
    call getmem3d(tvo,jci1,jci2,ici1,ici2,1,kz,'cu_grell:tvo')
    call getmem3d(xdby,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xdby')
    call getmem3d(xhe,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xhe')
    call getmem3d(xhes,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xhes')
    call getmem3d(xpwc,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xpwc')
    call getmem3d(xpwcd,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xpwcd')
    call getmem3d(xq,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xq')
    call getmem3d(xqc,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xqc')
    call getmem3d(xqes,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xqes')
    call getmem3d(xqrcd,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xqrcd')
    call getmem3d(xt,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xt')
    call getmem3d(xtv,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xtv')
    call getmem3d(xz,jci1,jci2,ici1,ici2,1,kz,'cu_grell:xz')
    call getmem3d(z,jci1,jci2,ici1,ici2,1,kz,'cu_grell:z')
    call getmem3d(zo,jci1,jci2,ici1,ici2,1,kz,'cu_grell:zo')
!
    call getmem2d(xac,jci1,jci2,ici1,ici2,'cu_grell:xac')
    call getmem2d(xao,jci1,jci2,ici1,ici2,'cu_grell:xao')
    call getmem2d(bu,jci1,jci2,ici1,ici2,'cu_grell:bu')
    call getmem2d(buo,jci1,jci2,ici1,ici2,'cu_grell:buo')
    call getmem2d(edt,jci1,jci2,ici1,ici2,'cu_grell:edt')
    call getmem2d(edto,jci1,jci2,ici1,ici2,'cu_grell:edto')
    call getmem2d(edtx,jci1,jci2,ici1,ici2,'cu_grell:edtx')
    call getmem2d(hcd,jci1,jci2,ici1,ici2,'cu_grell:hcd')
    call getmem2d(hcdo,jci1,jci2,ici1,ici2,'cu_grell:hcdo')
    call getmem2d(hkb,jci1,jci2,ici1,ici2,'cu_grell:hkb')
    call getmem2d(hkbo,jci1,jci2,ici1,ici2,'cu_grell:hkbo')
    call getmem2d(pwcav,jci1,jci2,ici1,ici2,'cu_grell:pwcav')
    call getmem2d(pwcavo,jci1,jci2,ici1,ici2,'cu_grell:pwcavo')
    call getmem2d(pwcev,jci1,jci2,ici1,ici2,'cu_grell:pwcev')
    call getmem2d(pwcevo,jci1,jci2,ici1,ici2,'cu_grell:pwcevo')
    call getmem2d(qcd,jci1,jci2,ici1,ici2,'cu_grell:qcd')
    call getmem2d(qcdo,jci1,jci2,ici1,ici2,'cu_grell:qcdo')
    call getmem2d(qck,jci1,jci2,ici1,ici2,'cu_grell:qck')
    call getmem2d(qcko,jci1,jci2,ici1,ici2,'cu_grell:qcko')
    call getmem2d(qkb,jci1,jci2,ici1,ici2,'cu_grell:qkb')
    call getmem2d(qkbo,jci1,jci2,ici1,ici2,'cu_grell:qkbo')
    call getmem2d(vshear,jci1,jci2,ici1,ici2,'cu_grell:vshear')
    call getmem2d(xxac,jci1,jci2,ici1,ici2,'cu_grell:xxac')
    call getmem2d(xhcd,jci1,jci2,ici1,ici2,'cu_grell:xhcd')
    call getmem2d(xhkb,jci1,jci2,ici1,ici2,'cu_grell:xhkb')
    call getmem2d(xmb,jci1,jci2,ici1,ici2,'cu_grell:xmb')
    call getmem2d(xpwcav,jci1,jci2,ici1,ici2,'cu_grell:xpwcav')
    call getmem2d(xpwcev,jci1,jci2,ici1,ici2,'cu_grell:xpwcev')
    call getmem2d(xqcd,jci1,jci2,ici1,ici2,'cu_grell:xqcd')
    call getmem2d(xqck,jci1,jci2,ici1,ici2,'cu_grell:xqck')
    call getmem2d(xqkb,jci1,jci2,ici1,ici2,'cu_grell:xqkb')
!
    alsixt = dlog(610.71D0)
    be(1) = ep2*wlhvocp*3.50D0
    be(2) = ep2*2.834D6*rcpd*3.50D0
    ae(1) = be(1)*rtzero + alsixt
    ae(2) = be(2)*rtzero + alsixt
  end subroutine allocate_mod_cu_grell

  subroutine cuparan(ktau)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    implicit none
!
    integer(8) , intent(in) :: ktau
!
    real(dp) :: pkdcut , pkk , prainx , us , vs
    integer :: i , j , k , jp1 , kk
!
    character (len=64) :: subroutine_name='cuparan'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    pkdcut = 75.0D0
!
!   prepare input, erase output
!
    outq(:,:,:) = d_zero
    outt(:,:,:) = d_zero
    p(:,:,:)    = d_zero
    po(:,:,:)   = d_zero
    q(:,:,:)    = d_zero
    qo(:,:,:)   = d_zero
    t(:,:,:)    = d_zero
    tn(:,:,:)   = d_zero
    vsp(:,:,:)  = d_zero
!
    pret(:,:)  = d_zero
    psur(:,:)  = d_zero
    qcrit(:,:) = d_zero
    ter11(:,:) = d_zero
    if ( lchem ) convpr(:,:,:) = d_zero 
!
    kdet(:,:)  = 2
    k22(:,:)   = 1
    ktop(:,:)  = 1
    kbcon(:,:) = 1
    kb(:,:)    = 1
    kds(:,:)   = 1
    kmin(:,:)  = 1
!
    dby(:,:,:) = d_zero
    dbyo(:,:,:) = d_zero
    dellah(:,:,:) = d_zero
    dellaq(:,:,:) = d_zero
    dellat(:,:,:) = d_zero
    dkk(:,:,:) = d_one
    he(:,:,:) = d_zero
    heo(:,:,:) = d_zero
    hes(:,:,:) = d_zero
    heso(:,:,:) = d_zero
    pwc(:,:,:) = d_zero
    pwco(:,:,:) = d_zero
    pwcd(:,:,:) = d_zero
    pwcdo(:,:,:) = d_zero
    qc(:,:,:) = d_zero
    qco(:,:,:) = d_zero
    qes(:,:,:) = d_zero
    qeso(:,:,:) = d_zero
    qrcd(:,:,:) = d_zero
    qrcdo(:,:,:) = d_zero
    tv(:,:,:) = d_zero
    tvo(:,:,:) = d_zero
    xdby(:,:,:) = d_zero
    xhe(:,:,:) = d_zero
    xhes(:,:,:) = d_zero
    xpwc(:,:,:) = d_zero
    xpwcd(:,:,:) = d_zero
    xq(:,:,:) = d_zero
    xqc(:,:,:) = d_zero
    xqes(:,:,:) = d_zero
    xqrcd(:,:,:) = d_zero
    xt(:,:,:) = d_zero
    xtv(:,:,:) = d_zero
    xz(:,:,:) = d_zero
    z(:,:,:) = d_zero
    zo(:,:,:) = d_zero
!
    xac(:,:) = d_zero
    xao(:,:) = d_zero
    bu(:,:) = d_zero
    buo(:,:) = d_zero
    edt(:,:) = d_zero
    edto(:,:) = d_zero
    edtx(:,:) = d_zero
    hcd(:,:) = d_zero
    hcdo(:,:) = d_zero
    hkb(:,:) = d_zero
    hkbo(:,:) = d_zero
    pwcav(:,:) = d_zero
    pwcavo(:,:) = d_zero
    pwcev(:,:) = d_zero
    pwcevo(:,:) = d_zero
    qcd(:,:) = d_zero
    qcdo(:,:) = d_zero
    qck(:,:) = d_zero
    qcko(:,:) = d_zero
    qkb(:,:) = d_zero
    qkbo(:,:) = d_zero
    vshear(:,:) = d_zero
    xxac(:,:) = d_zero
    xhcd(:,:) = d_zero
    xhkb(:,:) = d_zero
    xmb(:,:) = d_zero
    xpwcav(:,:) = d_zero
    xpwcev(:,:) = d_zero
    xqcd(:,:) = d_zero
    xqck(:,:) = d_zero
    xqkb(:,:) = d_zero
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          kk = kz - k + 1
          jp1 = j + 1
          us = (puatm(j,i,kk)/sfcps(j,i)+       &
                puatm(j,i+1,kk)/sfcps(j,i+1)+   &
                puatm(jp1,i,kk)/sfcps(jp1,i)+   &
                puatm(jp1,i+1,kk)/sfcps(jp1,i+1))*d_rfour
          vs = (pvatm(j,i,kk)/sfcps(j,i)+       &
                pvatm(j,i+1,kk)/sfcps(j,i+1)+   &
                pvatm(jp1,i,kk)/sfcps(jp1,i)+   &
                pvatm(jp1,i+1,kk)/sfcps(jp1,i+1))*d_rfour
          t(j,i,k) = tas(j,i,kk)
          q(j,i,k) = qvas(j,i,kk)
          if ( q(j,i,k) < 1.0D-08 ) q(j,i,k) = 1.0D-08
          tn(j,i,k) = t(j,i,k) + (tten(j,i,kk))/sfcps(j,i)*dtcum
          qo(j,i,k) = q(j,i,k) + (qvten(j,i,kk))/sfcps(j,i)*dtcum
          p(j,i,k) = d_10*sfcps(j,i)*hlev(kk) + d_10*ptop
          vsp(j,i,k) = dsqrt(us**d_two+vs**d_two)
          if ( qo(j,i,k) < 1.0D-08 ) qo(j,i,k) = 1.0D-08
!
          po(j,i,k) = p(j,i,k)
          psur(j,i) = d_10*sfcps(j,i) + d_10*ptop
          outt(j,i,k) = d_zero
          pkk = psur(j,i) - po(j,i,k)
          if ( pkk <= pkdcut ) kdet(j,i) = kdet(j,i) + 1
          outq(j,i,k) = d_zero
          ter11(j,i) = sfhgt(j,i)*regrav
          if ( ter11(j,i) <= d_zero ) ter11(j,i) = 1.0D-05
          qcrit(j,i) = qcrit(j,i) + qvten(j,i,kk)
        end do
      end do
    end do
!
!   call cumulus parameterization
!
    call cup
!
!   return cumulus parameterization
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( pret(j,i) > d_zero ) then
            kk = kz - k + 1
            tten(j,i,kk) = sfcps(j,i)*outt(j,i,k) + tten(j,i,kk)
            qvten(j,i,kk) = sfcps(j,i)*outq(j,i,k) + qvten(j,i,kk)
          end if
        end do
      end do
    end do
!
!   rain in cm.
!
    total_precip_points = 0
    do i = ici1 , ici2
      do j = jci1 , jci2
        prainx = pret(j,i)*dtmdl
        if ( prainx > dlowval ) then
          rainc(j,i) = rainc(j,i) + prainx
!         precipitation rate for bats (mm/s)
          if ( ktau == 0 .and. debug_level > 2 ) then
            lmpcpc(j,i) = lmpcpc(j,i) + pret(j,i)
          else
            lmpcpc(j,i) = lmpcpc(j,i) + pret(j,i)*aprdiv
          end if
          total_precip_points = total_precip_points + 1
        end if
      end do
    end do
    call time_end(subroutine_name,idindx)
  end subroutine cuparan
!
!   GRELL CUMULUS SCHEME
!
  subroutine cup

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    implicit none
!
    real(dp) :: adw , akclth , aup , detdo , detdoq , dg , dh ,   &
               dhh , dp_s , dq , xdt , dv1 , dv1q , dv2 , dv2q , &
               dv3 , dv3q , dz , dz1 , dz2 , dzo , e , eo , f ,  &
               agamma , agamma0 , agamma1 , agamma2 , agammo ,   &
               agammo0 , mbdt , outtes , pbcdif , qrch , qrcho , &
               tvbar , tvbaro , xk
    integer :: i , j , k , iph , ipho , kbcono , kclth , kk , lpt
!
    character (len=64) :: subroutine_name='cup'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

    mbdt = dtcum*5.0D-03
    f  = -d_one
    xk = -d_one
!
!   environmental conditions, first heights
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          iph = 1
          ipho = 1
          if ( t(j,i,k) <= tcrit ) iph = 2
          if ( tn(j,i,k) <= tcrit ) ipho = 2
          e = dexp(ae(iph)-be(iph)/t(j,i,k))
          eo = dexp(ae(ipho)-be(ipho)/tn(j,i,k))
          qes(j,i,k) = ep2*e/(d_100*p(j,i,k)-(d_one-ep2)*e)
          qeso(j,i,k) = ep2*eo/(d_100*po(j,i,k)-(d_one-ep2)*eo)
          if ( qes(j,i,k) <= 1.0D-08 ) qes(j,i,k) = 1.0D-08
          if ( q(j,i,k) > qes(j,i,k) ) q(j,i,k) = qes(j,i,k)
          if ( qeso(j,i,k) <= 1.0D-08 ) qeso(j,i,k) = 1.0D-08
          if ( qo(j,i,k) > qeso(j,i,k) ) qo(j,i,k) = qeso(j,i,k)
          tv(j,i,k) = t(j,i,k) + 0.608D0*q(j,i,k)*t(j,i,k)
          tvo(j,i,k) = tn(j,i,k) + 0.608D0*qo(j,i,k)*tn(j,i,k)
        end do
      end do
    end do

    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( qcrit(j,i) <= d_zero ) then
          xac(j,i) = -d_one
        end if
        if ( icup /= 2 ) then
          if (cucontrol(j,i) /= 2 ) then
            xac(j,i) = -d_one
          end if
        end if
        z(j,i,1)  = ter11(j,i) - &
              (dlog(p(j,i,1))-dlog(psur(j,i)))*rgas*tv(j,i,1)*regrav
        zo(j,i,1) = ter11(j,i) - &
              (dlog(po(j,i,1))-dlog(psur(j,i)))*rgas*tvo(j,i,1)*regrav
      end do
    end do

    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tvbar = d_half*(tv(j,i,k)+tv(j,i,k-1))
          z(j,i,k) = z(j,i,k-1) - &
               (dlog(p(j,i,k))-dlog(p(j,i,k-1)))*rgas*tvbar*regrav
          tvbaro = d_half*(tvo(j,i,k)+tvo(j,i,k-1))
          zo(j,i,k) = zo(j,i,k-1) - &
               (dlog(po(j,i,k))-dlog(po(j,i,k-1)))*rgas*tvbaro*regrav
        end do
      end do
    end do
!
!   moist static energy
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          he(j,i,k) = egrav*z(j,i,k) + cpd*t(j,i,k) + wlhv*q(j,i,k)
          hes(j,i,k) = egrav*z(j,i,k) + cpd*t(j,i,k) + wlhv*qes(j,i,k)
          if ( he(j,i,k) >= hes(j,i,k) ) he(j,i,k) = hes(j,i,k)
          heo(j,i,k) = egrav*zo(j,i,k) + cpd*tn(j,i,k) + wlhv*qo(j,i,k)
          heso(j,i,k) = egrav*zo(j,i,k) + cpd*tn(j,i,k) + wlhv*qeso(j,i,k)
          if ( heo(j,i,k) >= heso(j,i,k) ) heo(j,i,k) = heso(j,i,k)
          xt(j,i,k) = t(j,i,k)
          xq(j,i,k) = q(j,i,k)
          xhe(j,i,k) = he(j,i,k)
          if ( k /= kz ) qrcd(j,i,k) = d_half*(qes(j,i,k)+qes(j,i,k+1))
          if ( k /= kz ) qrcdo(j,i,k) = d_half*(qeso(j,i,k)+qeso(j,i,k+1))
        end do
      end do
    end do
!
!   determine level with highest moist static energy content.
!
    call maximi2(he,1,kbmax2d,k22)

    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) >= d_zero ) then
          if ( k22(j,i) >= kbmax2d(j,i) ) then
            xac(j,i) = -d_one
            cycle
          end if
          hkb(j,i) = he(j,i,k22(j,i))
          qkb(j,i) = q(j,i,k22(j,i))
          hkbo(j,i) = heo(j,i,k22(j,i))
          qkbo(j,i) = qo(j,i,k22(j,i))
          qck(j,i) = qkb(j,i)
          qcko(j,i) = qkbo(j,i)
        end if
      end do
    end do
!
!   decide for convective cloud base
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) >= d_zero ) then
          do k = 1 , kdet(j,i)
            kk = kdet(j,i) - k + 1
!           dkk(j,i,kk) = 0.75D0*dkk(j,i,kk+1)
            dkk(j,i,k) = d_one - dble(kk)/dble(kdet(j,i))
          end do

120       continue

          kb(j,i) = k22(j,i)
          kbcon(j,i) = kb(j,i)

140       continue

          dh = d_half*(hes(j,i,kbcon(j,i))+hes(j,i,kbcon(j,i)+1))
          if ( hkb(j,i) < dh ) then
            kbcon(j,i) = kbcon(j,i) + 1
            if ( kbcon(j,i) > kbmax2d(j,i) ) then
              xac(j,i) = -d_one
              cycle
            end if
            go to 140
          else
!
!    after large-scale forcing is applied, possible lid should be removed!!!
!
            kbcono = kb(j,i)
!ictp
150         continue
            if ( kbcono > kbmax2d(j,i) ) then
              xac(j,i) = -d_one
              cycle
            end if
!ictp_
            dh = d_half*(heso(j,i,kbcono)+heso(j,i,kbcono+1))
            if ( hkbo(j,i) < dh ) then
              kbcono = kbcono + 1
              go to 150
            else
              pbcdif = -p(j,i,kbcono) + p(j,i,kb(j,i))
!             below was commented out
!             as uncommenting the following lines for experiment 2/5/95
              if ( pbcdif > pbcmax2d(j,i) ) then
!               this is where typo was (pbdcdif)
                k22(j,i) = k22(j,i) + 1
                if ( k22(j,i) >= kbmax2d(j,i) ) then
                  xac(j,i) = -d_one
                  cycle
                end if
                hkb(j,i) = he(j,i,k22(j,i))
                qkb(j,i) = q(j,i,k22(j,i))
                hkbo(j,i) = heo(j,i,k22(j,i))
                qkbo(j,i) = qo(j,i,k22(j,i))
                qck(j,i) = qkb(j,i)
                qcko(j,i) = qkbo(j,i)
                go to 120
              end if
            end if
          end if
        end if
      end do
    end do
!
!   downdraft originating level
!
    call minimi(he,kb,kz,kmin)
    call maximi1(vsp,1,kz,kds)
!
!   static control
!
!   determine cloud top
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) >= 0 ) then
          if ( kmin(j,i) <= 3 ) then
            xac(j,i) = -d_one
            cycle
          end if
          if ( kds(j,i) >= kz ) kds(j,i) = kz - 1
          if ( kds(j,i) <= kbcon(j,i) ) kds(j,i) = kbcon(j,i)
          dby(j,i,kz) = hkb(j,i) - hes(j,i,kz)
          dbyo(j,i,kz) = hkbo(j,i) - heso(j,i,kz)
        end if
      end do
    end do
    do k = 1 , kz - 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) > xacact ) then
            dby(j,i,k) = hkb(j,i) - d_half*(hes(j,i,k)+hes(j,i,k+1))
            dbyo(j,i,k) = hkbo(j,i) - d_half*(heso(j,i,k)+heso(j,i,k+1))
          end if
        end do
      end do
    end do
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) > xacact ) then
          do k = 2 , kz - kbcon(j,i) - 1
            kk = kz - k + 1
            if ( dby(j,i,kk) >= d_zero ) then
              ktop(j,i) = kk + 1
              go to 320
            end if
          end do
          xac(j,i) = -d_one
          go to 400
320       continue
          if ( ktop(j,i) > kz ) ktop(j,i) = kz
          if ( p(j,i,kbcon(j,i))-p(j,i,ktop(j,i)) < mincld2d(j,i) ) then
            xac(j,i) = -d_one
          end if
        end if
400     continue
      end do
    end do
!
!   moisture and cloud work functions
!
    do k = 2 , kz - 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) > xacact ) then
            if ( k > kbcon(j,i) ) then
              if ( k < ktop(j,i) ) then
                dz = d_half*(z(j,i,k+1)-z(j,i,k-1))
                dz1 = z(j,i,k) - z(j,i,k-1)
                agamma = (wlhvocp)*(wlhv/(rwat*(t(j,i,k)**d_two)))*qes(j,i,k)
                agamma0 = (wlhvocp) * &
                  (wlhv/(rwat*(t(j,i,k-1)**d_two)))*qes(j,i,k-1)
                qrch = qes(j,i,k) + &
                      (d_one/wlhv)*(agamma/(d_one+agamma))*dby(j,i,k)
                qc(j,i,k) = (qck(j,i)-qrch)/(d_one+c0*dz) + qrch
                pwc(j,i,k) = c0*dz*(qc(j,i,k)-qrch)
                qck(j,i) = qc(j,i,k)
                pwcav(j,i) = pwcav(j,i) + pwc(j,i,k)
                dz1 = z(j,i,k) - z(j,i,k-1)
                xac(j,i) = xac(j,i) + &
                       dz1*(egrav/(cpd*(d_half*(t(j,i,k)+t(j,i,k-1))))) * &
                       dby(j,i,k-1)/(d_one+d_half*(agamma+agamma0))
                dzo = d_half*(zo(j,i,k+1)-zo(j,i,k-1))
                dz2 = zo(j,i,k) - zo(j,i,k-1)
                agammo = (wlhvocp)*(wlhv/(rwat*(tn(j,i,k)**d_two)))*qeso(j,i,k)
                agammo0 = (wlhvocp) * &
                  (wlhv/(rwat*(tn(j,i,k-1)**d_two)))*qeso(j,i,k-1)
                qrcho = qeso(j,i,k) + &
                     (d_one/wlhv)*(agammo/(d_one+agammo))*dbyo(j,i,k)
                qco(j,i,k) = (qcko(j,i)-qrcho)/(d_one+c0*dzo) + qrcho
                pwco(j,i,k) = c0*dzo*(qco(j,i,k)-qrcho)
                qcko(j,i) = qco(j,i,k)
                pwcavo(j,i) = pwcavo(j,i) + pwco(j,i,k)
                xao(j,i) = xao(j,i) + &
                    dz2*(egrav/(cpd*((tn(j,i,k)+tn(j,i,k-1))*d_half))) * &
                    dbyo(j,i,k-1)/(d_one+d_half*(agammo+agammo0))
              end if
            end if
          end if
        end do
      end do
    end do
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) > xacact ) then
          k = ktop(j,i)
          dz = d_half*(z(j,i,k)-z(j,i,k-1))
          agamma = (wlhvocp)*(wlhv/(rwat*(t(j,i,k)**d_two)))*qes(j,i,k)
          qrch = qes(j,i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*dby(j,i,k)
          qc(j,i,k) = qes(j,i,k)
          pwc(j,i,k) = (qrch-qes(j,i,k))
          pwcav(j,i) = pwcav(j,i) + pwc(j,i,k)
!
          dz = d_half*(zo(j,i,k)-zo(j,i,k-1))
          agamma = (wlhvocp)*(wlhv/(rwat*(tn(j,i,k)**d_two)))*qeso(j,i,k)
          qrcho = qeso(j,i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*dbyo(j,i,k)
          qco(j,i,k) = qeso(j,i,k)
          pwco(j,i,k) = (qrcho-qeso(j,i,k))
          pwcavo(j,i) = pwcavo(j,i) + pwco(j,i,k)
        end if
      end do
    end do
!
!   downdraft calculations
!
!   determine downdraft strength in terms of windshear
!
    do kk = 1 , kz/2
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) > xacact ) then
            vshear(j,i) = vshear(j,i) + &
                 dabs((vsp(j,i,kk+1)-vsp(j,i,kk))/(z(j,i,kk+1)-z(j,i,kk)))
          end if
        end do
      end do
    end do
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) > xacact ) then
          vshear(j,i) = d_1000*vshear(j,i)/dble(kz/2)
          edt(j,i) = d_one - &
             (1.591D0-0.639D0*vshear(j,i)+0.0953D0*(vshear(j,i)**d_two) - &
              0.00496D0*(vshear(j,i)**d_three))
          if ( edt(j,i) > shrmax2d(j,i) ) edt(j,i) = shrmax2d(j,i)
          if ( edt(j,i) < shrmin2d(j,i) ) edt(j,i) = shrmin2d(j,i)

          edto(j,i) = edt(j,i)
          edtx(j,i) = edt(j,i)
          qrcd(j,i,kz) = qes(j,i,kz)
          hcd(j,i) = d_half*(he(j,i,kmin(j,i))+he(j,i,kmin(j,i)+1))
          qcd(j,i) = d_half*(q(j,i,kmin(j,i))+q(j,i,kmin(j,i)+1))
          qrcdo(j,i,kz) = qeso(j,i,kz)
          hcdo(j,i) = heso(j,i,kz)
          hcdo(j,i) = d_half*(heo(j,i,kmin(j,i))+heo(j,i,kmin(j,i)+1))
          qcdo(j,i) = d_half*(qo(j,i,kmin(j,i))+qo(j,i,kmin(j,i)+1))
          bu(j,i) = d_zero
          buo(j,i) = d_zero
        end if
      end do
    end do
    do k = 1 , kz - 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) > xacact ) then
            if ( k < kmin(j,i) ) then
              kk = kmin(j,i) - k
              dz = -(z(j,i,kk)-z(j,i,kk+2))*d_half
              bu(j,i) = bu(j,i) + &
                dz*(hcd(j,i)-d_half*(hes(j,i,kk)+hes(j,i,kk+1)))
              dq = (qes(j,i,kk)+qes(j,i,kk+1))*d_half
              xdt = (t(j,i,kk)+t(j,i,kk+1))*d_half
              agamma = (wlhvocp)*(wlhv/(rwat*(xdt**d_two)))*dq
              dh = hcd(j,i) - d_half*(hes(j,i,kk)+hes(j,i,kk+1))
              qrcd(j,i,kk) = (dq+(d_one/wlhv)*(agamma/(d_one+agamma))*dh)
              pwcd(j,i,kk) = dkk(j,i,kk)*(qcd(j,i)-qrcd(j,i,kk))
              qcd(j,i) = qrcd(j,i,kk)
              pwcev(j,i) = pwcev(j,i) + pwcd(j,i,kk)
!
              dz = d_half*(zo(j,i,kk+2)-zo(j,i,kk))
              buo(j,i) = buo(j,i) + &
                dz*(hcdo(j,i)-(heso(j,i,kk)+heso(j,i,kk+1))*d_half)
              dq = (qeso(j,i,kk)+qeso(j,i,kk+1))*d_half
              xdt = (tn(j,i,kk)+tn(j,i,kk+1))*d_half
              agamma = (wlhvocp)*(wlhv/(rwat*(xdt**d_two)))*dq
              dh = hcdo(j,i) - d_half*(heso(j,i,kk)+heso(j,i,kk+1))
              qrcdo(j,i,kk) = (dq+(d_one/wlhv)*(agamma/(d_one+agamma))*dh)
              pwcdo(j,i,kk) = dkk(j,i,kk)*(qcdo(j,i)-qrcdo(j,i,kk))
              qcdo(j,i) = qrcdo(j,i,kk)
              pwcevo(j,i) = pwcevo(j,i) + pwcdo(j,i,kk)
            end if
          end if
        end do
      end do
    end do
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) > xacact ) then
          if ( bu(j,i) >= d_zero .or. buo(j,i) >= d_zero .or. &
            pwcev(j,i) >= d_zero .or. pwcevo(j,i) >= d_zero ) then
            xac(j,i) = -d_one
          end if
          edt(j,i) = -edt(j,i)*pwcav(j,i)/pwcev(j,i)
          if ( edt(j,i) > edtmax2d(j,i) ) edt(j,i) = edtmax2d(j,i)
          if ( edt(j,i) < edtmin2d(j,i) ) edt(j,i) = edtmin2d(j,i)
          edto(j,i) = -edto(j,i)*pwcavo(j,i)/pwcevo(j,i)
          if ( edto(j,i) > edtmaxo2d(j,i) ) edto(j,i) = edtmaxo2d(j,i)
          if ( edto(j,i) < edtmino2d(j,i) ) edto(j,i) = edtmino2d(j,i)
        end if
      end do
    end do
!
!   what would the change be?
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) > xacact ) then
          k = 1
          dz = d_half*(z(j,i,2)-z(j,i,1))
          dp_s = 50.0D0*(psur(j,i)-p(j,i,2))
          dellah(j,i,1) = edt(j,i)*(dkk(j,i,1)*hcd(j,i)-dkk(j,i,1) * &
                        d_half*(he(j,i,1)+he(j,i,2)))*egrav/dp_s
          dellaq(j,i,1) = edt(j,i)*(dkk(j,i,1)*qrcd(j,i,1)-dkk(j,i,1) * &
                        d_half*(q(j,i,1)+q(j,i,2)))*egrav/dp_s
          xhe(j,i,k) = dellah(j,i,k)*mbdt + he(j,i,k)
          xq(j,i,k) = dellaq(j,i,k)*mbdt + q(j,i,k)
          dellat(j,i,k) = rcpd*(dellah(j,i,k)-wlhv*dellaq(j,i,k))
          xt(j,i,k) = (mbdt*rcpd)*(dellah(j,i,k)-wlhv*dellaq(j,i,k))+t(j,i,k)
          if ( xq(j,i,k) <= d_zero ) xq(j,i,k) = 1.0D-08
        end if
      end do
    end do
!
    do k = 1 , kz - 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) > xacact ) then
            if ( k /= 1 .and. k < ktop(j,i) ) then
              dv1 = d_half*(he(j,i,k)+he(j,i,k+1))
              dv2 = he(j,i,k)
              dv3 = d_half*(he(j,i,k)+he(j,i,k-1))
              dv1q = d_half*(q(j,i,k)+q(j,i,k+1))
              dv2q = q(j,i,k)
              dv3q = d_half*(q(j,i,k)+q(j,i,k-1))
!
!             specifiy detrainment of downdraft, has to be consistent
!             with zd calculations in soundd.
!
              detdo = (d_one-dkk(j,i,k))*(hcd(j,i)-dv2)
              detdoq = (d_one-dkk(j,i,k))*(qrcd(j,i,k)-dv2q)
              dz = d_half*(z(j,i,k+1)-z(j,i,k-1))
!
!              changed due to subsidence and entrainment
!
              aup = d_one
              if ( k <= k22(j,i) ) aup = d_zero
              adw = d_one
              if ( k > kmin(j,i) ) adw = d_zero
              dp_s = +50.0D0*(p(j,i,k-1)-p(j,i,k+1))
              dellah(j,i,k) = ((aup-adw*edt(j,i)) * &
                    (dv1-dv2)+(aup-adw*edt(j,i)) * &
                    (dv2-dv3))*egrav/dp_s + adw*edt(j,i)*detdo*egrav/dp_s
              dellaq(j,i,k) = ((aup-adw*edt(j,i)) * &
                 (dv1q-dv2q)+(aup-adw*edt(j,i)) * &
                 (dv2q-dv3q))*egrav/dp_s + adw*edt(j,i)*detdoq*egrav/dp_s
              xhe(j,i,k) = dellah(j,i,k)*mbdt + he(j,i,k)
              xq(j,i,k) = dellaq(j,i,k)*mbdt + q(j,i,k)
              dellat(j,i,k) = rcpd*(dellah(j,i,k)-wlhv*dellaq(j,i,k))
              xt(j,i,k) = (mbdt*rcpd)*(dellah(j,i,k) - &
                           wlhv*dellaq(j,i,k)) + t(j,i,k)
              if ( xq(j,i,k) <= d_zero ) xq(j,i,k) = 1.0D-08
            end if
          end if
        end do
      end do
    end do
!
!   cloud top
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) > xacact ) then
          lpt = ktop(j,i)
          dp_s = d_100*(p(j,i,lpt-1)-p(j,i,lpt))
          dv1 = d_half*(he(j,i,lpt)+he(j,i,lpt-1))
          dellah(j,i,lpt) = (hkb(j,i)-dv1)*egrav/dp_s
          dv1 = d_half*(q(j,i,lpt)+q(j,i,lpt-1))
          dellaq(j,i,lpt) = (qes(j,i,lpt)-dv1)*egrav/dp_s
          k = lpt
          xhe(j,i,k) = dellah(j,i,k)*mbdt + he(j,i,k)
          xq(j,i,k) = dellaq(j,i,k)*mbdt + q(j,i,k)
          dellat(j,i,k) = rcpd*(dellah(j,i,k)-wlhv*dellaq(j,i,k))
          xt(j,i,k) = (mbdt*rcpd)*(dellah(j,i,k)-wlhv*dellaq(j,i,k)) + t(j,i,k)
          if ( xq(j,i,k) <= d_zero ) xq(j,i,k) = 1.0D-08
          xhkb(j,i) = dellah(j,i,kbcon(j,i))*mbdt + hkb(j,i)
          xqkb(j,i) = dellaq(j,i,kbcon(j,i))*mbdt + qkb(j,i)
          if ( xqkb(j,i) <= d_zero ) xqkb(j,i) = 1.0D-08
        end if
      end do
    end do
!
!   environmental conditions, first heights
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) > xacact ) then
            iph = 1
            if ( xt(j,i,k) <= tcrit ) iph = 2
            e = dexp(ae(iph)-be(iph)/xt(j,i,k))
            xqes(j,i,k) = ep2*e/(d_100*p(j,i,k)-(d_one-ep2)*e)
            if ( xqes(j,i,k) <= 1.0D-08 ) xqes(j,i,k) = 1.0D-08
            if ( xq(j,i,k) > xqes(j,i,k) ) xq(j,i,k) = xqes(j,i,k)
            xtv(j,i,k) = xt(j,i,k) + 0.608D0*xq(j,i,k)*xt(j,i,k)
          end if
        end do
      end do
    end do
!   bug fix
    do k = 1 , kz - 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) > xacact ) then
            xqrcd(j,i,k) = d_half*(xqes(j,i,k)+xqes(j,i,k+1))
          end if
        end do
      end do
    end do
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) > xacact ) then
          xz(j,i,1) = ter11(j,i) - &
              (dlog(p(j,i,1))-dlog(psur(j,i)))*rgas*xtv(j,i,1)*regrav
         end if
       end do
    end do
    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) > xacact ) then
            tvbar = d_half*(xtv(j,i,k)+xtv(j,i,k-1))
            xz(j,i,k) = xz(j,i,k-1) - &
                (dlog(p(j,i,k))-dlog(p(j,i,k-1)))*rgas*tvbar*regrav
          end if
        end do
      end do
    end do
!
!   moist static energy
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) > xacact ) then
            xhes(j,i,k) = egrav*xz(j,i,k) + cpd*xt(j,i,k) + wlhv*xqes(j,i,k)
            if ( xhe(j,i,k) >= xhes(j,i,k) ) xhe(j,i,k) = xhes(j,i,k)
          end if
        end do
      end do
    end do
!
!   static control
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) > xacact ) then
          xqck(j,i) = xqkb(j,i)
          xdby(j,i,kz) = xhkb(j,i) - xhes(j,i,kz)
        end if
      end do
    end do
!
!   moisture and cloud work functions
!
    do k = 1 , kz - 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) >= d_zero ) then
            xdby(j,i,k) = xhkb(j,i) - d_half*(xhes(j,i,k)+xhes(j,i,k+1))
            if ( k > kbcon(j,i) .and. k < ktop(j,i) ) then
              dz = d_half*(xz(j,i,k+1)-xz(j,i,k-1))
              dz1 = xz(j,i,k) - xz(j,i,k-1)
              agamma = (wlhvocp)*(wlhv/(rwat*(xt(j,i,k)**d_two)))*xqes(j,i,k)
              agamma0 = (wlhvocp) * &
                (wlhv/(rwat*(xt(j,i,k-1)**d_two)))*xqes(j,i,k-1)
              qrch = xqes(j,i,k) + &
                 (d_one/wlhv)*(agamma/(d_one+agamma))*xdby(j,i,k)
              xqc(j,i,k) = (xqck(j,i)-qrch)/(d_one+c0*dz) + qrch
              xpwc(j,i,k) = c0*dz*(xqc(j,i,k)-qrch)
              xqck(j,i) = xqc(j,i,k)
              xpwcav(j,i) = xpwcav(j,i) + xpwc(j,i,k)
              xxac(j,i) = xxac(j,i) + &
                dz1*(egrav/(cpd*(d_half*(xt(j,i,k)+xt(j,i,k-1))))) * &
                xdby(j,i,k-1)/(d_one+d_half*(agamma+agamma0))
            end if
          end if
        end do
      end do
    end do
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) >= d_zero ) then
          k = ktop(j,i)
          dz = d_half*(xz(j,i,k)-xz(j,i,k-1))
          agamma = (wlhvocp)*(wlhv/(rwat*(xt(j,i,k)**d_two)))*xqes(j,i,k)
          qrch = xqes(j,i,k) + (d_one/wlhv)*(agamma/(d_one+agamma))*xdby(j,i,k)
          xqc(j,i,k) = xqes(j,i,k)
          xpwc(j,i,k) = (qrch-xqes(j,i,k))
          xpwcav(j,i) = xpwcav(j,i) + xpwc(j,i,k)
          xqrcd(j,i,kz) = xqes(j,i,kz)
          xhcd(j,i) = d_half*(xhe(j,i,kmin(j,i))+xhe(j,i,kmin(j,i)+1))
          xqcd(j,i) = d_half*(xq(j,i,kmin(j,i))+xq(j,i,kmin(j,i)+1))
          xpwcev(j,i) = d_zero
          bu(j,i) = d_zero
        end if
      end do
    end do
!
!   downdraft calculations
!
!   downdraft moisture properties
!
    do k = 1 , kz - 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) >= d_zero ) then
            if ( k < kmin(j,i) ) then
              kk = kmin(j,i) - k
              dz = -(xz(j,i,kk)-xz(j,i,kk+2))*d_half
              bu(j,i) = bu(j,i) + &
                dz*(xhcd(j,i)-d_half*(xhes(j,i,kk)+xhes(j,i,kk+1)))
              dq = d_half*(xqes(j,i,kk)+xqes(j,i,kk+1))
              xdt = d_half*(xt(j,i,kk)+xt(j,i,kk+1))
              agamma = (wlhvocp)*(wlhv/(rwat*(xdt**d_two)))*dq
              dh = xhcd(j,i) - d_half*(xhes(j,i,kk)+xhes(j,i,kk+1))
              xqrcd(j,i,kk) = (dq+(d_one/wlhv)*(agamma/(d_one+agamma))*dh)
              xpwcd(j,i,kk) = dkk(j,i,kk)*(xqcd(j,i)-xqrcd(j,i,kk))
              xqcd(j,i) = xqrcd(j,i,kk)
              xpwcev(j,i) = xpwcev(j,i) + xpwcd(j,i,kk)
            end if
          end if
        end do
      end do
    end do
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) >= d_zero ) then
          if ( bu(j,i) >= d_zero ) then
            xac(j,i) = -d_one
            cycle
          end if
          if ( dabs(xpwcev(j,i)) > dlowval ) then
            edtx(j,i) = -edtx(j,i)*xpwcav(j,i)/xpwcev(j,i)
          end if
          if ( edtx(j,i) > edtmaxx2d(j,i) ) then
            edtx(j,i) = edtmaxx2d(j,i)
          end if
          if ( edtx(j,i) < edtminx2d(j,i) ) then
            edtx(j,i) = edtminx2d(j,i)
          end if
        end if
      end do
    end do
!
!   downdraft cloudwork functions
!
    do k = 1 , kz - 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) >= d_zero ) then
            if ( k < kmin(j,i) ) then
              kk = kmin(j,i) - k
!
!             original
!
              agamma1 = (wlhvocp)*(wlhv/(rwat*(t(j,i,kk)**d_two)))*qes(j,i,kk)
              agamma2 = (wlhvocp) * &
                (wlhv/(rwat*(t(j,i,kk+1)**d_two)))*qes(j,i,kk+1)
              dhh = hcd(j,i)
              xdt = d_half*(t(j,i,kk)+t(j,i,kk+1))
              dg = d_half*(agamma1+agamma2)
              dh = d_half*(hes(j,i,kk)+hes(j,i,kk+1))
              dz = (z(j,i,kk)-z(j,i,kk+1))*dkk(j,i,kk)
              xac(j,i) = xac(j,i) + &
                 edt(j,i)*dz*(egrav/(cpd*xdt))*((dhh-dh)/(d_one+dg))
!
!             modified by larger scale
!
              agamma1 = (wlhvocp)*(wlhv/(rwat*(tn(j,i,kk)**d_two)))*qeso(j,i,kk)
              agamma2 = (wlhvocp) * &
                (wlhv/(rwat*(tn(j,i,kk+1)**d_two)))*qeso(j,i,kk+1)
              dhh = hcdo(j,i)
              xdt = d_half*(tn(j,i,kk)+tn(j,i,kk+1))
              dg = d_half*(agamma1+agamma2)
              dh = d_half*(heso(j,i,kk)+heso(j,i,kk+1))
              dz = (zo(j,i,kk)-zo(j,i,kk+1))*dkk(j,i,kk)
              xao(j,i) = xao(j,i) + &
                edto(j,i)*dz*(egrav/(cpd*xdt))*((dhh-dh)/(d_one+dg))
!
!             modified by cloud
!
              agamma1 = (wlhvocp)*(wlhv/(rwat*(xt(j,i,kk)**d_two)))*xqes(j,i,kk)
              agamma2 = (wlhvocp) * &
                (wlhv/(rwat*(xt(j,i,kk+1)**d_two)))*xqes(j,i,kk+1)
              dhh = xhcd(j,i)
              xdt = d_half*(xt(j,i,kk)+xt(j,i,kk+1))
              dg = d_half*(agamma1+agamma2)
              dh = d_half*(xhes(j,i,kk)+xhes(j,i,kk+1))
              dz = (xz(j,i,kk)-xz(j,i,kk+1))*dkk(j,i,kk)
              xxac(j,i) = xxac(j,i) + edtx(j,i)*dz * &
                        (egrav/(cpd*xdt))*((dhh-dh)/(d_one+dg))
            end if
          end if
        end do
      end do
    end do
!
!   large scale forcing
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( xac(j,i) >= d_zero ) then
          if ( igcc == 1 ) then
            f = (xao(j,i)-xac(j,i))/dtcum ! Arakawa-Schubert closure
          else if ( igcc == 2 ) then
            f = xac(j,i)/dtauc2d(j,i)    ! Fritsch-Chappell closure
          end if
          xk = (xxac(j,i)-xac(j,i))/mbdt
          xmb(j,i) = -f/xk
          if ( f <= d_zero .or. xk >= d_zero ) xmb(j,i) = d_zero
        end if
        mflx(j,i,1) = xmb(j,i)
        mflx(j,i,2) = xmb(j,i)*edt(j,i)
      end do
    end do
!
!   feedback
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( xac(j,i) >= d_zero ) then
            if ( k <= ktop(j,i) ) then
              outtes = dellat(j,i,k)*xmb(j,i)*secpd
              if ( (outtes > htmax2d(j,i)) .or. (outtes < htmin2d(j,i)) ) then
                xmb(j,i) = d_zero
                xac(j,i) = -d_one
              else
                outt(j,i,k) = outt(j,i,k) + dellat(j,i,k)*xmb(j,i)
                outq(j,i,k) = outq(j,i,k) + dellaq(j,i,k)*xmb(j,i)
                pret(j,i) = pret(j,i) + &
                  (pwc(j,i,k)+edt(j,i)*pwcd(j,i,k))*xmb(j,i)
              end if
            end if
          end if
        end do
      end do
    end do
    if ( lchem ) then
      do k = 2 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( xac(j,i) >= d_zero ) then
              if ( k <= ktop(j,i) ) then
                outtes = dellat(j,i,k)*xmb(j,i)*secpd
                if ( (outtes < htmax2d(j,i)) .and. &
                     (outtes > htmin2d(j,i)) ) then
                  !FAB save the layer rain rate for chem removal
                  convpr(j,i,k) = convpr(j,i,k-1) + &
                     (pwc(j,i,k)+edt(j,i)*pwcd(j,i,k))*xmb(j,i)
                else if ( outtes < htmin2d(j,i) )
                  convpr(j,i,k) = pret(j,i)
                end if
              end if
            end if
          end do
        end do
      end do
    end if
!
!   calculate cloud fraction and water content
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        icumtop(j,i) = 0
        icumbot(j,i) = 0
        icumdwd(j,i) = 0
!
        if ( xac(j,i) >= d_zero ) then

          if ( ktop(j,i) > 1 .and. kbcon(j,i) > 1 ) then
            kclth = ktop(j,i) - kbcon(j,i) + 1
            akclth = d_one/dble(kclth)
            do k = kbcon(j,i) , ktop(j,i)
              kk = kz - k + 1
              rcldlwc(j,i,kk) = cllwcv
              rcldfra(j,i,kk) = d_one - (d_one-clfrcv)**akclth
            end do
!
!           define convection  base and top for tracers
!
            if ( lchem ) then
              if ( ktop(j,i) > 1 .and. k22(j,i) >= 1 ) then
                icumtop(j,i) = kzp1 - ktop(j,i)
                icumbot(j,i) = kzp1 - k22(j,i)
                icumdwd(j,i) = kzp1 - kmin(j,i)
              end if
            end if
          end if
        end if
      end do
    end do

    call time_end(subroutine_name,idindx)

    contains
!
     subroutine minimi(array,ks,ke,kt)
!
       implicit none
!
       integer , intent (in) :: ke
       real(dp) , intent(in) , pointer , dimension(:,:,:) :: array
       integer , intent(in) , pointer , dimension(:,:) :: ks
       integer , intent(out) , pointer , dimension(:,:) :: kt
!
       integer :: i , j , k
       real(dp) :: x
       character (len=64) :: subroutine_name='minimi'
       integer :: idindx=0
!
       call time_begin(subroutine_name,idindx)
!
       do i = ici1 , ici2
         do j = jci1 , jci2
           kt(j,i) = ks(j,i)
           x = array(j,i,ks(j,i))
           do k = ks(j,i) + 1 , ke
             if ( array(j,i,k) < x ) then
               x = array(j,i,k)
               kt(j,i) = k
             end if
           end do
         end do
       end do
       call time_end(subroutine_name,idindx) 
     end subroutine minimi
!
     subroutine maximi1(array,ks,ke,imax)
!
      implicit none
!
      integer , intent (in) :: ks , ke
      real(dp) , intent(in) , pointer , dimension(:,:,:) :: array
      integer , intent(out) , pointer , dimension(:,:) :: imax
!
      integer :: i , j , k
      real(dp) :: x , xar
!
      character (len=64) :: subroutine_name='maximi1'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
      do i = ici1 , ici2
        do j = jci1 , jci2
          imax(j,i) = ks
          x = array(j,i,ks)
          do k = ks , ke
            xar = array(j,i,k)
            if ( xar >= x ) then
              x = xar
              imax(j,i) = k
            end if
          end do
        end do
      end do
!
      call time_end(subroutine_name,idindx) 
    end subroutine maximi1

     subroutine maximi2(array,ks,ke,imax)
!
      implicit none
!
      integer , intent (in) :: ks
      real(dp) , intent(in) , pointer , dimension(:,:,:) :: array
      integer , intent(in) , pointer , dimension(:,:) :: ke
      integer , intent(out) , pointer , dimension(:,:) :: imax
!
      integer :: i , j , k
      real(dp) :: x , xar
!
      character (len=64) :: subroutine_name='maximi2'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
      do i = ici1 , ici2
        do j = jci1 , jci2
          imax(j,i) = ks
          x = array(j,i,ks)
          do k = ks , ke(j,i)
            xar = array(j,i,k)
            if ( xar >= x ) then
              x = xar
              imax(j,i) = k
            end if
          end do
        end do
      end do
!
      call time_end(subroutine_name,idindx) 
    end subroutine maximi2

  end subroutine cup

end module mod_cu_grell

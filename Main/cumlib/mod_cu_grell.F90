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

  use mod_intkinds
  use mod_realkinds
  use mod_service
  use mod_memutil
  use mod_cu_common
  use mod_mpmessage
  use mod_runparams , only : iqv , dt , dtsec , igcc , ichem
  use mod_regcm_types

  implicit none

  private

  real(rk8) , parameter :: xacact = -0.99999D0
  real(rk8) , parameter :: tcrit = 50.0D0
  real(rk8) , parameter :: c0 = 0.002D0
  real(rk8) :: alsixt
  real(rk8) , dimension(2) :: ae , be
  real(rk8) , pointer , dimension(:,:) :: outq , outt , p , po ,  &
                              q , qo , t , tn , vsp
  real(rk8) , pointer , dimension(:,:) :: dby , dbyo , dellah ,       &
              dellaq , dellat , dkk , he , heo , hes , heso , pwc , &
              pwcd , pwcdo , pwco , qc , qco , qes , qeso , qrcd ,  &
              qrcdo , tv , tvo , xdby , xhe , xhes , xpwc , xpwcd , &
              xq , xqc , xqes , xqrcd , xt , xtv , xz , z , zo
  real(rk8) , pointer , dimension(:) :: pratec , psur , qcrit , ter11
  integer(ik4) , pointer , dimension(:) :: kdet
  integer(ik4) , pointer , dimension(:) :: kmin , k22 , kb , kbcon , &
              kds , ktop , iac , jac
  real(rk8) , pointer , dimension(:) :: xac , xao , bu , buo , edt ,  &
              edto , edtx , hcd , hcdo , hkb , hkbo , pwcav ,       &
              pwcavo , pwcev , pwcevo , qcd , qcdo , qck , qcko ,   &
              qkb , qkbo , vshear , xxac , xhcd , xhkb , xmb ,      &
              xpwcav , xpwcev , xqcd , xqck , xqkb

  real(rk8) , pointer , dimension(:) :: dtauc , pbcmax ,   &
              mincld , shrmax , shrmin , edtmax , edtmin ,    &
              edtmaxo , edtmaxx , edtmino , edtminx , htmax , &
              htmin
  integer(ik4) , public , pointer , dimension(:) :: kbmax
  real(rk8) , public , pointer , dimension(:,:) :: dtauc2d , pbcmax2d ,   &
              mincld2d , shrmax2d , shrmin2d , edtmax2d , edtmin2d ,    &
              edtmaxo2d , edtmaxx2d , edtmino2d , edtminx2d , htmax2d , &
              htmin2d
  integer(ik4) , public , pointer , dimension(:,:) :: kbmax2d

  integer :: ncp , nap

  public :: allocate_mod_cu_grell , cuparan

  contains

  subroutine allocate_mod_cu_grell
    implicit none

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

    ncp = (jci2-jci1+1)*(ici2-ici1+1)

    call getmem2d(outq,1,ncp,1,kz,'cu_grell:outq')
    call getmem2d(outt,1,ncp,1,kz,'cu_grell:outt')
    call getmem2d(p,1,ncp,1,kz,'cu_grell:p')
    call getmem2d(po,1,ncp,1,kz,'cu_grell:po')
    call getmem2d(q,1,ncp,1,kz,'cu_grell:q')
    call getmem2d(qo,1,ncp,1,kz,'cu_grell:qo')
    call getmem2d(t,1,ncp,1,kz,'cu_grell:t')
    call getmem2d(tn,1,ncp,1,kz,'cu_grell:tn')
    call getmem2d(vsp,1,ncp,1,kz,'cu_grell:vsp')

    call getmem1d(pratec,1,ncp,'cu_grell:pratec')
    call getmem1d(psur,1,ncp,'cu_grell:psur')
    call getmem1d(qcrit,1,ncp,'cu_grell:qcrit')
    call getmem1d(ter11,1,ncp,'cu_grell:ter11')

    call getmem1d(dtauc,1,ncp,'cu_grell:dtauc')
    call getmem1d(pbcmax,1,ncp,'cu_grell:pbcmax')
    call getmem1d(mincld,1,ncp,'cu_grell:mincld')
    call getmem1d(shrmax,1,ncp,'cu_grell:shrmax')
    call getmem1d(shrmin,1,ncp,'cu_grell:shrmin')
    call getmem1d(edtmax,1,ncp,'cu_grell:edtmax')
    call getmem1d(edtmin,1,ncp,'cu_grell:edtmin')
    call getmem1d(edtmaxo,1,ncp,'cu_grell:edtmaxo')
    call getmem1d(edtmino,1,ncp,'cu_grell:edtmino')
    call getmem1d(edtmaxx,1,ncp,'cu_grell:edtmaxx')
    call getmem1d(edtminx,1,ncp,'cu_grell:edtminx')
    call getmem1d(htmax,1,ncp,'cu_grell:htmax')
    call getmem1d(htmin,1,ncp,'cu_grell:htmin')
    call getmem1d(kbmax,1,ncp,'cu_grell:kbmax')

    call getmem1d(kdet,1,ncp,'cu_grell:kdet')
    call getmem1d(kmin,1,ncp,'cu_grell:kmin')
    call getmem1d(k22,1,ncp,'cu_grell:k22')
    call getmem1d(kb,1,ncp,'cu_grell:kb')
    call getmem1d(kbcon,1,ncp,'cu_grell:kbcon')
    call getmem1d(kds,1,ncp,'cu_grell:kds')
    call getmem1d(ktop,1,ncp,'cu_grell:ktop')
    call getmem1d(iac,1,ncp,'cu_grell:iac')
    call getmem1d(jac,1,ncp,'cu_grell:jac')

    call getmem2d(dby,1,ncp,1,kz,'cu_grell:dby')
    call getmem2d(dbyo,1,ncp,1,kz,'cu_grell:dbyo')
    call getmem2d(dellah,1,ncp,1,kz,'cu_grell:dellah')
    call getmem2d(dellaq,1,ncp,1,kz,'cu_grell:dellaq')
    call getmem2d(dellat,1,ncp,1,kz,'cu_grell:dellat')
    call getmem2d(dkk,1,ncp,1,kz,'cu_grell:dkk')
    call getmem2d(he,1,ncp,1,kz,'cu_grell:he')
    call getmem2d(heo,1,ncp,1,kz,'cu_grell:heo')
    call getmem2d(hes,1,ncp,1,kz,'cu_grell:hes')
    call getmem2d(heso,1,ncp,1,kz,'cu_grell:heso')
    call getmem2d(pwc,1,ncp,1,kz,'cu_grell:pwc')
    call getmem2d(pwco,1,ncp,1,kz,'cu_grell:pwco')
    call getmem2d(pwcd,1,ncp,1,kz,'cu_grell:pwcd')
    call getmem2d(pwcdo,1,ncp,1,kz,'cu_grell:pwcdo')
    call getmem2d(qc,1,ncp,1,kz,'cu_grell:qc')
    call getmem2d(qco,1,ncp,1,kz,'cu_grell:qco')
    call getmem2d(qes,1,ncp,1,kz,'cu_grell:qes')
    call getmem2d(qeso,1,ncp,1,kz,'cu_grell:qeso')
    call getmem2d(qrcd,1,ncp,1,kz,'cu_grell:qrcd')
    call getmem2d(qrcdo,1,ncp,1,kz,'cu_grell:qrcdo')
    call getmem2d(tv,1,ncp,1,kz,'cu_grell:tv')
    call getmem2d(tvo,1,ncp,1,kz,'cu_grell:tvo')
    call getmem2d(xdby,1,ncp,1,kz,'cu_grell:xdby')
    call getmem2d(xhe,1,ncp,1,kz,'cu_grell:xhe')
    call getmem2d(xhes,1,ncp,1,kz,'cu_grell:xhes')
    call getmem2d(xpwc,1,ncp,1,kz,'cu_grell:xpwc')
    call getmem2d(xpwcd,1,ncp,1,kz,'cu_grell:xpwcd')
    call getmem2d(xq,1,ncp,1,kz,'cu_grell:xq')
    call getmem2d(xqc,1,ncp,1,kz,'cu_grell:xqc')
    call getmem2d(xqes,1,ncp,1,kz,'cu_grell:xqes')
    call getmem2d(xqrcd,1,ncp,1,kz,'cu_grell:xqrcd')
    call getmem2d(xt,1,ncp,1,kz,'cu_grell:xt')
    call getmem2d(xtv,1,ncp,1,kz,'cu_grell:xtv')
    call getmem2d(xz,1,ncp,1,kz,'cu_grell:xz')
    call getmem2d(z,1,ncp,1,kz,'cu_grell:z')
    call getmem2d(zo,1,ncp,1,kz,'cu_grell:zo')

    call getmem1d(xac,1,ncp,'cu_grell:xac')
    call getmem1d(xao,1,ncp,'cu_grell:xao')
    call getmem1d(bu,1,ncp,'cu_grell:bu')
    call getmem1d(buo,1,ncp,'cu_grell:buo')
    call getmem1d(edt,1,ncp,'cu_grell:edt')
    call getmem1d(edto,1,ncp,'cu_grell:edto')
    call getmem1d(edtx,1,ncp,'cu_grell:edtx')
    call getmem1d(hcd,1,ncp,'cu_grell:hcd')
    call getmem1d(hcdo,1,ncp,'cu_grell:hcdo')
    call getmem1d(hkb,1,ncp,'cu_grell:hkb')
    call getmem1d(hkbo,1,ncp,'cu_grell:hkbo')
    call getmem1d(pwcav,1,ncp,'cu_grell:pwcav')
    call getmem1d(pwcavo,1,ncp,'cu_grell:pwcavo')
    call getmem1d(pwcev,1,ncp,'cu_grell:pwcev')
    call getmem1d(pwcevo,1,ncp,'cu_grell:pwcevo')
    call getmem1d(qcd,1,ncp,'cu_grell:qcd')
    call getmem1d(qcdo,1,ncp,'cu_grell:qcdo')
    call getmem1d(qck,1,ncp,'cu_grell:qck')
    call getmem1d(qcko,1,ncp,'cu_grell:qcko')
    call getmem1d(qkb,1,ncp,'cu_grell:qkb')
    call getmem1d(qkbo,1,ncp,'cu_grell:qkbo')
    call getmem1d(vshear,1,ncp,'cu_grell:vshear')
    call getmem1d(xxac,1,ncp,'cu_grell:xxac')
    call getmem1d(xhcd,1,ncp,'cu_grell:xhcd')
    call getmem1d(xhkb,1,ncp,'cu_grell:xhkb')
    call getmem1d(xmb,1,ncp,'cu_grell:xmb')
    call getmem1d(xpwcav,1,ncp,'cu_grell:xpwcav')
    call getmem1d(xpwcev,1,ncp,'cu_grell:xpwcev')
    call getmem1d(xqcd,1,ncp,'cu_grell:xqcd')
    call getmem1d(xqck,1,ncp,'cu_grell:xqck')
    call getmem1d(xqkb,1,ncp,'cu_grell:xqkb')

    alsixt = dlog(610.71D0)
    be(1) = ep2*wlhvocp*3.50D0
    be(2) = ep2*2.834D6*rcpd*3.50D0
    ae(1) = be(1)*rtzero + alsixt
    ae(2) = be(2)*rtzero + alsixt
  end subroutine allocate_mod_cu_grell

  subroutine cuparan(m2c,c2m)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c
    type(cum_2_mod) , intent(inout) :: c2m
    real(rk8) :: pkdcut , pkk , prainx
    integer(ik4) :: i , j , k , kk , n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cuparan'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    pkdcut = 75.0D0
    !
    ! prepare input, erase output
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

    pratec(:)  = d_zero
    psur(:)  = d_zero
    qcrit(:) = d_zero
    ter11(:) = d_zero
    if ( ichem == 1 ) c2m%convpr(:,:,:) = d_zero

    kdet(:)  = 2
    k22(:)   = 1
    ktop(:)  = 1
    kbcon(:) = 1
    kb(:)    = 1
    kds(:)   = 1
    kmin(:)  = 1

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

    nap = 0
    do i = ici1 , ici2
      do j = jci1 , jci2
        if (cuscheme(j,i) == 2 ) then
          nap = nap + 1
          iac(nap) = i
          jac(nap) = j
        end if
      end do
    end do

    if ( nap == 0 ) return

    ! Pressures in millibar here.

    do n = 1 , nap
      i = iac(n)
      j = jac(n)
      psur(n) = m2c%psf(j,i) * d_r100
      ter11(n) = m2c%ht(j,i) * regrav
      dtauc(n) = dtauc2d(j,i)
      pbcmax(n) = pbcmax2d(j,i)
      mincld(n) = mincld2d(j,i)
      shrmax(n) = shrmax2d(j,i)
      shrmin(n) = shrmin2d(j,i)
      edtmax(n) = edtmax2d(j,i)
      edtmin(n) = edtmin2d(j,i)
      edtmaxo(n) = edtmaxo2d(j,i)
      edtmino(n) = edtmino2d(j,i)
      edtmaxx(n) = edtmaxx2d(j,i)
      edtminx(n) = edtminx2d(j,i)
      htmax(n) = htmax2d(j,i)
      htmin(n) = htmin2d(j,i)
      kbmax(n) = kbmax2d(j,i)
    end do

    ! Grell scheme requires values bottom -> top

    do k = 1 , kz
      kk = kz - k + 1
      do n = 1 , nap
        i = iac(n)
        j = jac(n)
        vsp(n,k) = sqrt(m2c%uas(j,i,kk)**2 + m2c%vas(j,i,kk)**2)
        p(n,k) = m2c%pas(j,i,kk) * d_r100
        t(n,k) = m2c%tas(j,i,kk)
        q(n,k) = max(m2c%qxas(j,i,kk,iqv),1.0D-08)
        po(n,k) = p(n,k)
        tn(n,k) = t(n,k) + (c2m%tten(j,i,kk))/m2c%psb(j,i)*dt
        qo(n,k) = q(n,k) + (c2m%qxten(j,i,kk,iqv))/m2c%psb(j,i)*dt
        if ( qo(n,k) < 1.0D-08 ) qo(n,k) = 1.0D-08
        pkk = psur(n) - po(n,k)
        qcrit(n) = qcrit(n) + c2m%qxten(j,i,kk,iqv)/m2c%psb(j,i)
        if ( pkk <= pkdcut ) kdet(n) = kdet(n) + 1
      end do
    end do

    if ( maxval(kdet) > kz ) then
      write(stderr,*) 'At point ', iac(maxloc(kdet)), jac(maxloc(kdet))
      write(stderr,*) 'Convection has reached column top!'
      write(stderr,*) 'Convection scheme cannot solve.'
      write(stderr,*) 'Try decreasing dt.'
      call fatal(__FILE__,__LINE__,'GRELL INSTABILITY!')
    end if
    !
    ! call cumulus parameterization
    !
    call cup
    !
    ! return cumulus parameterization
    !
    do k = 1 , kz
      kk = kzp1 - k
      do n = 1 , nap
        i = iac(n)
        j = jac(n)
        c2m%tten(j,i,kk) = m2c%psb(j,i)*outt(n,k) + c2m%tten(j,i,kk)
        c2m%qxten(j,i,kk,iqv) = m2c%psb(j,i)*outq(n,k) + c2m%qxten(j,i,kk,iqv)
      end do
    end do

    ! build for chemistry 3d table of constant precipitation rate
    ! from the surface to the top of the convection
    if ( ichem == 1 ) then
      do n = 1 , nap
        i = iac(n)
        j = jac(n)
        do k = 1 , ktop(n)-1
          c2m%convpr(j,i,kz-k+1) = pratec(n)
        end do
      end do
    end if
    !
    ! calculate cloud fraction and water content
    !
    c2m%kcumtop(:,:) = 0
    c2m%kcumbot(:,:) = 0
    do n = 1 , nap
      i = iac(n)
      j = jac(n)
      if ( xac(n) >= d_zero ) then
        !
        ! define convection base and top
        !
        if ( ktop(n) > 1 .and. kbcon(n) > 1 ) then
          if ( ktop(n) > 1 .and. kbcon(n) >= 1 ) then
            c2m%kcumtop(j,i) = kzp1 - ktop(n)
            c2m%kcumbot(j,i) = kzp1 - kbcon(n)
          end if
        end if
      end if
    end do

    do n = 1 , nap
      i = iac(n)
      j = jac(n)
      prainx = pratec(n)*dtsec
      if ( prainx > dlowval ) then
        c2m%rainc(j,i) = c2m%rainc(j,i) + prainx
        ! precipitation rate for surface (kg m-2 s-1)
        c2m%pcratec(j,i) = c2m%pcratec(j,i) + pratec(n)
        total_precip_points = total_precip_points + 1
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cuparan
  !
  ! GRELL CUMULUS SCHEME
  !
  subroutine cup
    implicit none
    real(rk8) :: adw , aup , detdo , detdoq , dg , dh ,   &
               dhh , dp_s , dq , xdt , dv1 , dv1q , dv2 , dv2q , &
               dv3 , dv3q , dz , dz1 , dz2 , dzo , e , eo , f ,  &
               agamma , agamma0 , agamma1 , agamma2 , agammo ,   &
               agammo0 , mbdt , outtes , pbcdif , qrch , qrcho , &
               tvbar , tvbaro , xk
    integer(ik4) :: n , k , iph , ipho , kbcono , kk , lpt
    logical :: foundtop
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cup'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    mbdt = dt*5.0D-03
    f  = -d_one
    xk = -d_one
    !
    ! environmental conditions, first heights
    !
    do k = 1 , kz
      do n = 1 , nap
        iph = 1
        ipho = 1
        if ( t(n,k) <= tcrit ) iph = 2
        if ( tn(n,k) <= tcrit ) ipho = 2
        e = dexp(ae(iph)-be(iph)/t(n,k))
        eo = dexp(ae(ipho)-be(ipho)/tn(n,k))
        qes(n,k) = ep2*e/(d_100*p(n,k)-(d_one-ep2)*e)
        qeso(n,k) = ep2*eo/(d_100*po(n,k)-(d_one-ep2)*eo)
        if ( qes(n,k) <= 1.0D-08 ) qes(n,k) = 1.0D-08
        if ( q(n,k) > qes(n,k) ) q(n,k) = qes(n,k)
        if ( qeso(n,k) <= 1.0D-08 ) qeso(n,k) = 1.0D-08
        if ( qo(n,k) > qeso(n,k) ) qo(n,k) = qeso(n,k)
        tv(n,k) = t(n,k) + 0.608D0*q(n,k)*t(n,k)
        tvo(n,k) = tn(n,k) + 0.608D0*qo(n,k)*tn(n,k)
      end do
    end do

    do n = 1 , nap
      if ( qcrit(n) <= d_zero ) then
        xac(n) = -d_one
      end if
      z(n,1)  = ter11(n) - &
            (dlog(p(n,1))-dlog(psur(n)))*rgas*tv(n,1)*regrav
      zo(n,1) = ter11(n) - &
            (dlog(po(n,1))-dlog(psur(n)))*rgas*tvo(n,1)*regrav
    end do

    do k = 2 , kz
      do n = 1 , nap
        tvbar = d_half*(tv(n,k)+tv(n,k-1))
        z(n,k) = z(n,k-1) - &
             (dlog(p(n,k))-dlog(p(n,k-1)))*rgas*tvbar*regrav
        tvbaro = d_half*(tvo(n,k)+tvo(n,k-1))
        zo(n,k) = zo(n,k-1) - &
             (dlog(po(n,k))-dlog(po(n,k-1)))*rgas*tvbaro*regrav
      end do
    end do
    !
    ! moist static energy
    !
    do k = 1 , kz
      do n = 1 , nap
        he(n,k) = egrav*z(n,k) + cpd*t(n,k) + wlhv*q(n,k)
        hes(n,k) = egrav*z(n,k) + cpd*t(n,k) + wlhv*qes(n,k)
        if ( he(n,k) >= hes(n,k) ) he(n,k) = hes(n,k)
        heo(n,k) = egrav*zo(n,k) + cpd*tn(n,k) + wlhv*qo(n,k)
        heso(n,k) = egrav*zo(n,k) + cpd*tn(n,k) + wlhv*qeso(n,k)
        if ( heo(n,k) >= heso(n,k) ) heo(n,k) = heso(n,k)
        xt(n,k) = t(n,k)
        xq(n,k) = q(n,k)
        xhe(n,k) = he(n,k)
        if ( k /= kz ) qrcd(n,k) = d_half*(qes(n,k)+qes(n,k+1))
        if ( k /= kz ) qrcdo(n,k) = d_half*(qeso(n,k)+qeso(n,k+1))
      end do
    end do
    !
    ! determine level with highest moist static energy content.
    !
    call maximi2(he,1,kbmax,k22)

    do n = 1 , nap
      if ( xac(n) >= d_zero ) then
        if ( k22(n) >= kbmax(n) ) then
          xac(n) = -d_one
          cycle
        end if
        hkb(n) = he(n,k22(n))
        qkb(n) = q(n,k22(n))
        hkbo(n) = heo(n,k22(n))
        qkbo(n) = qo(n,k22(n))
        qck(n) = qkb(n)
        qcko(n) = qkbo(n)
      end if
    end do
    !
    ! decide for convective cloud base
    !
    pointloop: &
    do n = 1 , nap
      if ( xac(n) >= d_zero ) then
        do k = 1 , kdet(n)
          kk = kdet(n) - k + 1
          ! dkk(n,kk) = 0.75D0 * dkk(n,kk+1)
          dkk(n,k) = d_one - dble(kk)/dble(kdet(n))
        end do
        cloudbase: do
          kb(n) = k22(n)
          kbcon(n) = kb(n)
          gotonext: do
            dh = d_half*(hes(n,kbcon(n))+hes(n,kbcon(n)+1))
            if ( hkb(n) < dh ) then
              kbcon(n) = kbcon(n) + 1
              if ( kbcon(n) > kbmax(n) ) then
                xac(n) = -d_one
                cycle pointloop
              end if
              cycle gotonext
            else
              !
              ! after large-scale forcing is applied, possible lid should
              ! be removed!!!
              !
              kbcono = kb(n)
              removelid: do
                if ( kbcono > kbmax(n) ) then
                  xac(n) = -d_one
                  cycle pointloop
                end if
                dh = d_half*(heso(n,kbcono)+heso(n,kbcono+1))
                if ( hkbo(n) < dh ) then
                  kbcono = kbcono + 1
                  cycle removelid
                else
                  pbcdif = -p(n,kbcono) + p(n,kb(n))
                  ! below was commented out
                  ! as uncommenting the following lines for experiment 2/5/95
                  if ( pbcdif > pbcmax(n) ) then
                    ! this is where typo was (pbdcdif)
                    k22(n) = k22(n) + 1
                    if ( k22(n) >= kbmax(n) ) then
                      xac(n) = -d_one
                      cycle pointloop
                    end if
                    hkb(n) = he(n,k22(n))
                    qkb(n) = q(n,k22(n))
                    hkbo(n) = heo(n,k22(n))
                    qkbo(n) = qo(n,k22(n))
                    qck(n) = qkb(n)
                    qcko(n) = qkbo(n)
                    cycle cloudbase
                  end if
                end if
                exit removelid
              end do removelid
            end if
            exit gotonext
          end do gotonext
          exit cloudbase
        end do cloudbase
      end if
    end do pointloop
    !
    ! downdraft originating level
    !
    call minimi(he,kb,kz,kmin)
    call maximi1(vsp,1,kz,kds)
    !
    ! static control
    !
    ! determine cloud top
    !
    do n = 1 , nap
      if ( xac(n) >= 0 ) then
        if ( kmin(n) <= 3 ) then
          xac(n) = -d_one
          cycle
        end if
        if ( kds(n) >= kz ) kds(n) = kz - 1
        if ( kds(n) <= kbcon(n) ) kds(n) = kbcon(n)
        dby(n,kz) = hkb(n) - hes(n,kz)
        dbyo(n,kz) = hkbo(n) - heso(n,kz)
      end if
    end do
    do k = 1 , kz - 1
      do n = 1 , nap
        if ( xac(n) > xacact ) then
          dby(n,k) = hkb(n) - d_half*(hes(n,k)+hes(n,k+1))
          dbyo(n,k) = hkbo(n) - d_half*(heso(n,k)+heso(n,k+1))
        end if
      end do
    end do
    do n = 1 , nap
      if ( xac(n) > xacact ) then
        foundtop = .false.
        findtop: &
        do k = 2 , kz - kbcon(n) - 1
          kk = kz - k + 1
          if ( dby(n,kk) >= d_zero ) then
            ktop(n) = kk + 1
            foundtop = .true.
            exit findtop
          end if
        end do findtop
        if ( foundtop ) then
          if ( ktop(n) > kz ) ktop(n) = kz
          if ( p(n,kbcon(n))-p(n,ktop(n)) < mincld(n) ) then
            xac(n) = -d_one
          end if
        else
          xac(n) = -d_one
        end if
      end if
    end do
    !
    ! moisture and cloud work functions
    !
    do k = 2 , kz - 1
      do n = 1 , nap
        if ( xac(n) > xacact ) then
          if ( k > kbcon(n) ) then
            if ( k < ktop(n) ) then
              dz = d_half*(z(n,k+1)-z(n,k-1))
              dz1 = z(n,k) - z(n,k-1)
              agamma = (wlhvocp)*(wlhv/(rwat*(t(n,k)**2)))*qes(n,k)
              agamma0 = (wlhvocp) * &
                (wlhv/(rwat*(t(n,k-1)**2)))*qes(n,k-1)
              qrch = qes(n,k) + rwlhv*(agamma/(d_one+agamma))*dby(n,k)
              qc(n,k) = (qck(n)-qrch)/(d_one+c0*dz) + qrch
              pwc(n,k) = c0*dz*(qc(n,k)-qrch)
              qck(n) = qc(n,k)
              pwcav(n) = pwcav(n) + pwc(n,k)
              dz1 = z(n,k) - z(n,k-1)
              xac(n) = xac(n) + &
                     dz1*(egrav/(cpd*(d_half*(t(n,k)+t(n,k-1))))) * &
                     dby(n,k-1)/(d_one+d_half*(agamma+agamma0))
              dzo = d_half*(zo(n,k+1)-zo(n,k-1))
              dz2 = zo(n,k) - zo(n,k-1)
              agammo = (wlhvocp)*(wlhv/(rwat*(tn(n,k)**2)))*qeso(n,k)
              agammo0 = (wlhvocp) * &
                (wlhv/(rwat*(tn(n,k-1)**2)))*qeso(n,k-1)
              qrcho = qeso(n,k) + rwlhv*(agammo/(d_one+agammo))*dbyo(n,k)
              qco(n,k) = (qcko(n)-qrcho)/(d_one+c0*dzo) + qrcho
              pwco(n,k) = c0*dzo*(qco(n,k)-qrcho)
              qcko(n) = qco(n,k)
              pwcavo(n) = pwcavo(n) + pwco(n,k)
              xao(n) = xao(n) + &
                  dz2*(egrav/(cpd*((tn(n,k)+tn(n,k-1))*d_half))) * &
                  dbyo(n,k-1)/(d_one+d_half*(agammo+agammo0))
            end if
          end if
        end if
      end do
    end do

    do n = 1 , nap
      if ( xac(n) > xacact ) then
        k = ktop(n)
        dz = d_half*(z(n,k)-z(n,k-1))
        agamma = (wlhvocp)*(wlhv/(rwat*(t(n,k)**2)))*qes(n,k)
        qrch = qes(n,k) + rwlhv*(agamma/(d_one+agamma))*dby(n,k)
        qc(n,k) = qes(n,k)
        pwc(n,k) = (qrch-qes(n,k))
        pwcav(n) = pwcav(n) + pwc(n,k)
        dz = d_half*(zo(n,k)-zo(n,k-1))
        agamma = (wlhvocp)*(wlhv/(rwat*(tn(n,k)**2)))*qeso(n,k)
        qrcho = qeso(n,k) + rwlhv*(agamma/(d_one+agamma))*dbyo(n,k)
        qco(n,k) = qeso(n,k)
        pwco(n,k) = (qrcho-qeso(n,k))
        pwcavo(n) = pwcavo(n) + pwco(n,k)
      end if
    end do
    !
    ! downdraft calculations
    !
    ! determine downdraft strength in terms of windshear
    !
    do kk = 1 , kz/2
      do n = 1 , nap
        if ( xac(n) > xacact ) then
          vshear(n) = vshear(n) + &
               dabs((vsp(n,kk+1)-vsp(n,kk))/(z(n,kk+1)-z(n,kk)))
        end if
      end do
    end do
    do n = 1 , nap
      if ( xac(n) > xacact ) then
        vshear(n) = d_1000*vshear(n)/dble(kz/2)
        edt(n) = d_one - &
           (1.591D0-0.639D0*vshear(n)+0.0953D0*(vshear(n)**2) - &
            0.00496D0*(vshear(n)**3))
        if ( edt(n) > shrmax(n) ) edt(n) = shrmax(n)
        if ( edt(n) < shrmin(n) ) edt(n) = shrmin(n)

        edto(n) = edt(n)
        edtx(n) = edt(n)
        qrcd(n,kz) = qes(n,kz)
        hcd(n) = d_half*(he(n,kmin(n))+he(n,kmin(n)+1))
        qcd(n) = d_half*(q(n,kmin(n))+q(n,kmin(n)+1))
        qrcdo(n,kz) = qeso(n,kz)
        hcdo(n) = heso(n,kz)
        hcdo(n) = d_half*(heo(n,kmin(n))+heo(n,kmin(n)+1))
        qcdo(n) = d_half*(qo(n,kmin(n))+qo(n,kmin(n)+1))
        bu(n) = d_zero
        buo(n) = d_zero
      end if
    end do
    do k = 1 , kz - 1
      do n = 1 , nap
        if ( xac(n) > xacact ) then
          if ( k < kmin(n) ) then
            kk = kmin(n) - k
            dz = -(z(n,kk)-z(n,kk+2))*d_half
            bu(n) = bu(n) + &
              dz*(hcd(n)-d_half*(hes(n,kk)+hes(n,kk+1)))
            dq = (qes(n,kk)+qes(n,kk+1))*d_half
            xdt = (t(n,kk)+t(n,kk+1))*d_half
            agamma = (wlhvocp)*(wlhv/(rwat*(xdt**2)))*dq
            dh = hcd(n) - d_half*(hes(n,kk)+hes(n,kk+1))
            qrcd(n,kk) = dq + rwlhv*(agamma/(d_one+agamma))*dh
            pwcd(n,kk) = dkk(n,kk)*(qcd(n)-qrcd(n,kk))
            qcd(n) = qrcd(n,kk)
            pwcev(n) = pwcev(n) + pwcd(n,kk)
            dz = d_half*(zo(n,kk+2)-zo(n,kk))
            buo(n) = buo(n) + &
              dz*(hcdo(n)-(heso(n,kk)+heso(n,kk+1))*d_half)
            dq = (qeso(n,kk)+qeso(n,kk+1))*d_half
            xdt = (tn(n,kk)+tn(n,kk+1))*d_half
            agamma = (wlhvocp)*(wlhv/(rwat*(xdt**2)))*dq
            dh = hcdo(n) - d_half*(heso(n,kk)+heso(n,kk+1))
            qrcdo(n,kk) = dq + rwlhv*(agamma/(d_one+agamma))*dh
            pwcdo(n,kk) = dkk(n,kk)*(qcdo(n)-qrcdo(n,kk))
            qcdo(n) = qrcdo(n,kk)
            pwcevo(n) = pwcevo(n) + pwcdo(n,kk)
          end if
        end if
      end do
    end do

    do n = 1 , nap
      if ( xac(n) > xacact ) then
        if ( bu(n) >= d_zero .or. buo(n) >= d_zero .or. &
          pwcev(n) >= d_zero .or. pwcevo(n) >= d_zero ) then
          xac(n) = -d_one
        end if
        edt(n) = -edt(n)*pwcav(n)/pwcev(n)
        if ( edt(n) > edtmax(n) ) edt(n) = edtmax(n)
        if ( edt(n) < edtmin(n) ) edt(n) = edtmin(n)
        edto(n) = -edto(n)*pwcavo(n)/pwcevo(n)
        if ( edto(n) > edtmaxo(n) ) edto(n) = edtmaxo(n)
        if ( edto(n) < edtmino(n) ) edto(n) = edtmino(n)
      end if
    end do
    !
    ! what would the change be?
    !
    do n = 1 , nap
      if ( xac(n) > xacact ) then
        k = 1
        dz = d_half*(z(n,2)-z(n,1))
        dp_s = 50.0D0*(psur(n)-p(n,2))
        dellah(n,1) = edt(n)*(dkk(n,1)*hcd(n)-dkk(n,1) * &
                      d_half*(he(n,1)+he(n,2)))*egrav/dp_s
        dellaq(n,1) = edt(n)*(dkk(n,1)*qrcd(n,1)-dkk(n,1) * &
                      d_half*(q(n,1)+q(n,2)))*egrav/dp_s
        xhe(n,k) = dellah(n,k)*mbdt + he(n,k)
        xq(n,k) = dellaq(n,k)*mbdt + q(n,k)
        dellat(n,k) = rcpd*(dellah(n,k)-wlhv*dellaq(n,k))
        xt(n,k) = (mbdt*rcpd)*(dellah(n,k)-wlhv*dellaq(n,k))+t(n,k)
        if ( xq(n,k) <= d_zero ) xq(n,k) = 1.0D-08
      end if
    end do

    do k = 1 , kz - 1
      do n = 1 , nap
        if ( xac(n) > xacact ) then
          if ( k /= 1 .and. k < ktop(n) ) then
            dv1 = d_half*(he(n,k)+he(n,k+1))
            dv2 = he(n,k)
            dv3 = d_half*(he(n,k)+he(n,k-1))
            dv1q = d_half*(q(n,k)+q(n,k+1))
            dv2q = q(n,k)
            dv3q = d_half*(q(n,k)+q(n,k-1))
            !
            ! specifiy detrainment of downdraft, has to be consistent
            ! with zd calculations in soundd.
            !
            detdo = (d_one-dkk(n,k))*(hcd(n)-dv2)
            detdoq = (d_one-dkk(n,k))*(qrcd(n,k)-dv2q)
            dz = d_half*(z(n,k+1)-z(n,k-1))
            !
            ! changed due to subsidence and entrainment
            !
            aup = d_one
            if ( k <= k22(n) ) aup = d_zero
            adw = d_one
            if ( k > kmin(n) ) adw = d_zero
            dp_s = +50.0D0*(p(n,k-1)-p(n,k+1))
            dellah(n,k) = ((aup-adw*edt(n)) * &
                  (dv1-dv2)+(aup-adw*edt(n)) * &
                  (dv2-dv3))*egrav/dp_s + adw*edt(n)*detdo*egrav/dp_s
            dellaq(n,k) = ((aup-adw*edt(n)) * &
               (dv1q-dv2q)+(aup-adw*edt(n)) * &
               (dv2q-dv3q))*egrav/dp_s + adw*edt(n)*detdoq*egrav/dp_s
            xhe(n,k) = dellah(n,k)*mbdt + he(n,k)
            xq(n,k) = dellaq(n,k)*mbdt + q(n,k)
            dellat(n,k) = rcpd*(dellah(n,k)-wlhv*dellaq(n,k))
            xt(n,k) = (mbdt*rcpd)*(dellah(n,k) - &
                         wlhv*dellaq(n,k)) + t(n,k)
            if ( xq(n,k) <= d_zero ) xq(n,k) = 1.0D-08
          end if
        end if
      end do
    end do
    !
    ! cloud top
    !
    do n = 1 , nap
      if ( xac(n) > xacact ) then
        lpt = ktop(n)
        dp_s = d_100*(p(n,lpt-1)-p(n,lpt))
        dv1 = d_half*(he(n,lpt)+he(n,lpt-1))
        dellah(n,lpt) = (hkb(n)-dv1)*egrav/dp_s
        dv1 = d_half*(q(n,lpt)+q(n,lpt-1))
        dellaq(n,lpt) = (qes(n,lpt)-dv1)*egrav/dp_s
        k = lpt
        xhe(n,k) = dellah(n,k)*mbdt + he(n,k)
        xq(n,k) = dellaq(n,k)*mbdt + q(n,k)
        dellat(n,k) = rcpd*(dellah(n,k)-wlhv*dellaq(n,k))
        xt(n,k) = (mbdt*rcpd)*(dellah(n,k)-wlhv*dellaq(n,k)) + t(n,k)
        if ( xq(n,k) <= d_zero ) xq(n,k) = 1.0D-08
        xhkb(n) = dellah(n,kbcon(n))*mbdt + hkb(n)
        xqkb(n) = dellaq(n,kbcon(n))*mbdt + qkb(n)
        if ( xqkb(n) <= d_zero ) xqkb(n) = 1.0D-08
      end if
    end do
    !
    ! environmental conditions, first heights
    !
    do k = 1 , kz
      do n = 1 , nap
        if ( xac(n) > xacact ) then
          iph = 1
          if ( xt(n,k) <= tcrit ) iph = 2
          e = dexp(ae(iph)-be(iph)/xt(n,k))
          xqes(n,k) = ep2*e/(d_100*p(n,k)-(d_one-ep2)*e)
          if ( xqes(n,k) <= 1.0D-08 ) xqes(n,k) = 1.0D-08
          if ( xq(n,k) > xqes(n,k) ) xq(n,k) = xqes(n,k)
          xtv(n,k) = xt(n,k) + 0.608D0*xq(n,k)*xt(n,k)
        end if
      end do
    end do
    ! bug fix
    do k = 1 , kz - 1
      do n = 1 , nap
        if ( xac(n) > xacact ) then
          xqrcd(n,k) = d_half*(xqes(n,k)+xqes(n,k+1))
        end if
      end do
    end do

    do n = 1 , nap
      if ( xac(n) > xacact ) then
        xz(n,1) = ter11(n)-(dlog(p(n,1))-dlog(psur(n)))*rgas*xtv(n,1)*regrav
       end if
    end do
    do k = 2 , kz
      do n = 1 , nap
        if ( xac(n) > xacact ) then
          tvbar = d_half*(xtv(n,k)+xtv(n,k-1))
          xz(n,k) = xz(n,k-1)-(dlog(p(n,k))-dlog(p(n,k-1)))*rgas*tvbar*regrav
        end if
      end do
    end do
    !
    ! moist static energy
    !
    do k = 1 , kz
      do n = 1 , nap
        if ( xac(n) > xacact ) then
          xhes(n,k) = egrav*xz(n,k) + cpd*xt(n,k) + wlhv*xqes(n,k)
          if ( xhe(n,k) >= xhes(n,k) ) xhe(n,k) = xhes(n,k)
        end if
      end do
    end do
    !
    ! static control
    !
    do n = 1 , nap
      if ( xac(n) > xacact ) then
        xqck(n) = xqkb(n)
        xdby(n,kz) = xhkb(n) - xhes(n,kz)
      end if
    end do
    !
    ! moisture and cloud work functions
    !
    do k = 1 , kz - 1
      do n = 1 , nap
        if ( xac(n) >= d_zero ) then
          xdby(n,k) = xhkb(n) - d_half*(xhes(n,k)+xhes(n,k+1))
          if ( k > kbcon(n) .and. k < ktop(n) ) then
            dz = d_half*(xz(n,k+1)-xz(n,k-1))
            dz1 = xz(n,k) - xz(n,k-1)
            agamma = (wlhvocp)*(wlhv/(rwat*(xt(n,k)**2)))*xqes(n,k)
            agamma0 = (wlhvocp) * &
              (wlhv/(rwat*(xt(n,k-1)**2)))*xqes(n,k-1)
            qrch = xqes(n,k) + rwlhv*(agamma/(d_one+agamma))*xdby(n,k)
            xqc(n,k) = (xqck(n)-qrch)/(d_one+c0*dz) + qrch
            xpwc(n,k) = c0*dz*(xqc(n,k)-qrch)
            xqck(n) = xqc(n,k)
            xpwcav(n) = xpwcav(n) + xpwc(n,k)
            xxac(n) = xxac(n) + &
              dz1*(egrav/(cpd*(d_half*(xt(n,k)+xt(n,k-1))))) * &
              xdby(n,k-1)/(d_one+d_half*(agamma+agamma0))
          end if
        end if
      end do
    end do
    do n = 1 , nap
      if ( xac(n) >= d_zero ) then
        k = ktop(n)
        dz = d_half*(xz(n,k)-xz(n,k-1))
        agamma = (wlhvocp)*(wlhv/(rwat*(xt(n,k)**2)))*xqes(n,k)
        qrch = xqes(n,k) + rwlhv*(agamma/(d_one+agamma))*xdby(n,k)
        xqc(n,k) = xqes(n,k)
        xpwc(n,k) = (qrch-xqes(n,k))
        xpwcav(n) = xpwcav(n) + xpwc(n,k)
        xqrcd(n,kz) = xqes(n,kz)
        xhcd(n) = d_half*(xhe(n,kmin(n))+xhe(n,kmin(n)+1))
        xqcd(n) = d_half*(xq(n,kmin(n))+xq(n,kmin(n)+1))
        xpwcev(n) = d_zero
        bu(n) = d_zero
      end if
    end do
    !
    ! downdraft calculations
    !
    ! downdraft moisture properties
    !
    do k = 1 , kz - 1
      do n = 1 , nap
        if ( xac(n) >= d_zero ) then
          if ( k < kmin(n) ) then
            kk = kmin(n) - k
            dz = -(xz(n,kk)-xz(n,kk+2))*d_half
            bu(n) = bu(n) + dz*(xhcd(n)-d_half*(xhes(n,kk)+xhes(n,kk+1)))
            dq = d_half*(xqes(n,kk)+xqes(n,kk+1))
            xdt = d_half*(xt(n,kk)+xt(n,kk+1))
            agamma = (wlhvocp)*(wlhv/(rwat*(xdt**2)))*dq
            dh = xhcd(n) - d_half*(xhes(n,kk)+xhes(n,kk+1))
            xqrcd(n,kk) = dq + rwlhv*(agamma/(d_one+agamma))*dh
            xpwcd(n,kk) = dkk(n,kk)*(xqcd(n)-xqrcd(n,kk))
            xqcd(n) = xqrcd(n,kk)
            xpwcev(n) = xpwcev(n) + xpwcd(n,kk)
          end if
        end if
      end do
    end do
    do n = 1 , nap
      if ( xac(n) >= d_zero ) then
        if ( bu(n) >= d_zero ) then
          xac(n) = -d_one
          cycle
        end if
        if ( dabs(xpwcev(n)) > dlowval ) then
          edtx(n) = -edtx(n)*xpwcav(n)/xpwcev(n)
        end if
        if ( edtx(n) > edtmaxx(n) ) then
          edtx(n) = edtmaxx(n)
        end if
        if ( edtx(n) < edtminx(n) ) then
          edtx(n) = edtminx(n)
        end if
      end if
    end do
    !
    ! downdraft cloudwork functions
    !
    do k = 1 , kz - 1
      do n = 1 , nap
        if ( xac(n) >= d_zero ) then
          if ( k < kmin(n) ) then
            kk = kmin(n) - k
            !
            ! original
            !
            agamma1 = (wlhvocp)*(wlhv/(rwat*(t(n,kk)**2)))*qes(n,kk)
            agamma2 = (wlhvocp) * &
              (wlhv/(rwat*(t(n,kk+1)**2)))*qes(n,kk+1)
            dhh = hcd(n)
            xdt = d_half*(t(n,kk)+t(n,kk+1))
            dg = d_half*(agamma1+agamma2)
            dh = d_half*(hes(n,kk)+hes(n,kk+1))
            dz = (z(n,kk)-z(n,kk+1))*dkk(n,kk)
            xac(n) = xac(n) + &
               edt(n)*dz*(egrav/(cpd*xdt))*((dhh-dh)/(d_one+dg))
            !
            ! modified by larger scale
            !
            agamma1 = (wlhvocp)*(wlhv/(rwat*(tn(n,kk)**2)))*qeso(n,kk)
            agamma2 = (wlhvocp) * &
              (wlhv/(rwat*(tn(n,kk+1)**2)))*qeso(n,kk+1)
            dhh = hcdo(n)
            xdt = d_half*(tn(n,kk)+tn(n,kk+1))
            dg = d_half*(agamma1+agamma2)
            dh = d_half*(heso(n,kk)+heso(n,kk+1))
            dz = (zo(n,kk)-zo(n,kk+1))*dkk(n,kk)
            xao(n) = xao(n) + &
              edto(n)*dz*(egrav/(cpd*xdt))*((dhh-dh)/(d_one+dg))
            !
            ! modified by cloud
            !
            agamma1 = (wlhvocp)*(wlhv/(rwat*(xt(n,kk)**2)))*xqes(n,kk)
            agamma2 = (wlhvocp) * &
              (wlhv/(rwat*(xt(n,kk+1)**2)))*xqes(n,kk+1)
            dhh = xhcd(n)
            xdt = d_half*(xt(n,kk)+xt(n,kk+1))
            dg = d_half*(agamma1+agamma2)
            dh = d_half*(xhes(n,kk)+xhes(n,kk+1))
            dz = (xz(n,kk)-xz(n,kk+1))*dkk(n,kk)
            xxac(n) = xxac(n) + edtx(n)*dz * &
                      (egrav/(cpd*xdt))*((dhh-dh)/(d_one+dg))
          end if
        end if
      end do
    end do
    !
    ! large scale forcing
    !
    do n = 1 , nap
      if ( xac(n) >= d_zero ) then
        if ( igcc == 1 ) then
          f = (xao(n)-xac(n))/dt ! Arakawa-Schubert closure
        else if ( igcc == 2 ) then
          f = xac(n)/dtauc(n)  ! Fritsch-Chappell closure
        end if
        xk = (xxac(n)-xac(n))/mbdt
        xmb(n) = -f/xk
        if ( f <= d_zero .or. xk >= d_zero ) xmb(n) = d_zero
      end if
    end do
    !
    ! feedback
    !
    do k = 1 , kz
      do n = 1 , nap
        if ( xac(n) >= d_zero ) then
          if ( k <= ktop(n) ) then
            outtes = dellat(n,k)*xmb(n)*secpd
            if ( (outtes > htmax(n)) .or. (outtes < htmin(n)) ) then
              xmb(n) = d_zero
              xac(n) = -d_one
            else
              outt(n,k) = outt(n,k) + dellat(n,k)*xmb(n)
              outq(n,k) = outq(n,k) + dellaq(n,k)*xmb(n)
              pratec(n) = pratec(n) + &
                (pwc(n,k)+edt(n)*pwcd(n,k))*xmb(n)
            end if
          end if
        end if
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
    contains

     subroutine minimi(array,ks,ke,kt)
       implicit none
       integer(ik4) , intent (in) :: ke
       real(rk8) , intent(in) , pointer , dimension(:,:) :: array
       integer(ik4) , intent(in) , pointer , dimension(:) :: ks
       integer(ik4) , intent(inout) , pointer , dimension(:) :: kt
       integer(ik4) :: n , k
       real(rk8) :: x
#ifdef DEBUG
       character(len=dbgslen) :: subroutine_name = 'minimi'
       integer(ik4) , save :: idindx = 0
       call time_begin(subroutine_name,idindx)
#endif
       do n = 1 , nap
         kt(n) = ks(n)
         x = array(n,ks(n))
         do k = ks(n) + 1 , ke
           if ( array(n,k) < x ) then
             x = array(n,k)
             kt(n) = k
           end if
         end do
       end do
#ifdef DEBUG
       call time_end(subroutine_name,idindx)
#endif
     end subroutine minimi

     subroutine maximi1(array,ks,ke,imax)
      implicit none
      integer(ik4) , intent (in) :: ks , ke
      real(rk8) , intent(in) , pointer , dimension(:,:) :: array
      integer(ik4) , intent(inout) , pointer , dimension(:) :: imax
      integer(ik4) :: n , k
      real(rk8) :: x , xar
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'maximi1'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      do n = 1 , nap
        imax(n) = ks
        x = array(n,ks)
        do k = ks , ke
          xar = array(n,k)
          if ( xar >= x ) then
            x = xar
            imax(n) = k
          end if
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine maximi1

    subroutine maximi2(array,ks,ke,imax)
      implicit none
      integer(ik4) , intent (in) :: ks
      real(rk8) , intent(in) , pointer , dimension(:,:) :: array
      integer(ik4) , intent(in) , pointer , dimension(:) :: ke
      integer(ik4) , intent(inout) , pointer , dimension(:) :: imax
      integer(ik4) :: n , k
      real(rk8) :: x , xar
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'maximi2'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      do n = 1 , nap
        imax(n) = ks
        x = array(n,ks)
        do k = ks , ke(n)
          xar = array(n,k)
          if ( xar >= x ) then
            x = xar
            imax(n) = k
          end if
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine maximi2

  end subroutine cup

end module mod_cu_grell

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

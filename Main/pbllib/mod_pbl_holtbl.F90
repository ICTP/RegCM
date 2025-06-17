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

module mod_pbl_holtbl
  !
  ! Holtslag planetary boundary layer scheme
  ! Reference : Holtslag, De Bruijn and Pan - MWR - 8/90
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams, only : iqv, dt, rdt, ichem, ichdrdepo, &
        zhnew_fac, ifaholtth10, ifaholt, holtth10iter, ipptls, &
        iqc, iqi, dsigma
  use mod_mppparam
  use mod_memutil
  use mod_service
  use mod_mpmessage
  use mod_pbl_common
  use mod_regcm_types

  implicit none

  private

  public :: allocate_mod_pbl_holtbl, holtbl

  real(rkx), pointer, contiguous, dimension(:,:,:) :: cgh, cgs, kvc, kvh, &
                                          kvm, kvq
  real(rkx), pointer, contiguous, dimension(:,:) :: xhfx, xqfx, exns, pfcor
  real(rkx), pointer, contiguous, dimension(:,:) :: hfxv, obklen, thv10, ustr
  logical, pointer, contiguous, dimension(:,:) :: lunstb

  real(rkx), pointer, contiguous, dimension(:,:,:) :: alphak, betak, &
                        coefe, coeff1, coeff2, tpred1, tpred2, cfac, vv
  real(rkx), pointer, contiguous, dimension(:,:,:) :: kzm, ttnp
  real(rkx), pointer, contiguous, dimension(:,:,:) :: qtenv
  real(rkx), pointer, contiguous, dimension(:,:,:) :: qtenc
  real(rkx), pointer, contiguous, dimension(:,:,:) :: qteni
  real(rkx), pointer, contiguous, dimension(:,:) :: uvdrage
  real(rkx), pointer, contiguous, dimension(:,:,:) :: hydf

  real(rkx), pointer, contiguous, dimension(:,:,:) :: dza, thvx
  real(rkx), pointer, contiguous, dimension(:,:,:) :: akzz1, akzz2
  real(rkx), pointer, contiguous, dimension(:,:,:) :: rhohf

  real(rkx), pointer, contiguous, dimension(:,:,:) :: ri

  ! minimum eddy diffusivity ( background value )
  real(rkx), parameter :: kzo = 1.0_rkx ! m^2s-1
  real(rkx), parameter :: turbulent_lenght_scale = 100.0_rkx ! m
  real(rkx), parameter :: szkm = (turbulent_lenght_scale*vonkar)**2
  ! coef. of proportionality and lower % of bl in sfc layer
  real(rkx), parameter :: fak = 8.5_rkx
  real(rkx), parameter :: fakn = 7.2_rkx
  real(rkx), parameter :: sffrac = 0.1_rkx
  ! beta coefs. for momentum, stable conditions and heat
  real(rkx), parameter :: betam = 15.0_rkx
  real(rkx), parameter :: betas = 5.0_rkx
  real(rkx), parameter :: betah = 15.0_rkx
  real(rkx), parameter :: ccon = fak*sffrac*vonkar
  real(rkx), parameter :: gvk = egrav*vonkar
  real(rkx), parameter :: binm = betam*sffrac
  real(rkx), parameter :: binh = betah*sffrac
  real(rkx), parameter :: kzfrac = 0.8_rkx
  ! power in formula for k in critical ri for judging stability
  real(rkx), parameter :: pink = 2.0_rkx

  contains

  subroutine allocate_mod_pbl_holtbl
    implicit none
    call getmem3d(cfac,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:cfac')
    call getmem3d(vv,jci1,jci2,ici1,ici2,2,kz,'mod_holtbl:vv')
    call getmem3d(ri,1,kz,jci1,jci2,ici1,ici2,'mod_holtbl:ri')
    call getmem3d(kzm,jci1,jci2,ici1,ici2,2,kz,'mod_holtbl:kzm')
    call getmem3d(qtenv,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:qtenv')
    call getmem3d(qtenc,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:qtenc')
    if ( ipptls > 1 ) then
      call getmem3d(qteni,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:qteni')
    end if
    call getmem3d(ttnp,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:ttnp')
    call getmem3d(hydf,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:hydf')
    call getmem3d(cgh,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:cgh')
    call getmem3d(cgs,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:cgs')
    call getmem3d(kvh,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:kvh')
    call getmem3d(kvm,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:kvm')
    call getmem3d(kvq,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:kvq')
    call getmem2d(hfxv,jci1,jci2,ici1,ici2,'mod_holtbl:hfxv')
    call getmem2d(lunstb,jci1,jci2,ici1,ici2,'mod_holtbl:lunstb')
    call getmem2d(xhfx,jci1,jci2,ici1,ici2,'mod_holtbl:xhfx')
    call getmem2d(xqfx,jci1,jci2,ici1,ici2,'mod_holtbl:xqfx')
    call getmem2d(obklen,jci1,jci2,ici1,ici2,'mod_holtbl:obklen')
    call getmem2d(exns,jci1,jci2,ici1,ici2,'mod_holtbl:enxns')
    call getmem2d(pfcor,jci1,jci2,ici1,ici2,'mod_holtbl:pfcor')
    call getmem2d(thv10,jci1,jci2,ici1,ici2,'mod_holtbl:thv10')
    call getmem2d(ustr,jci1,jci2,ici1,ici2,'mod_holtbl:ustr')
    call getmem3d(thvx,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:thvx')
    call getmem3d(dza,jci1,jci2,ici1,ici2,1,kzm1,'mod_holtbl:dza')
    call getmem3d(rhohf,jci1,jci2,ici1,ici2,1,kzm1,'mod_holtbl:rhohf')
    call getmem3d(alphak,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:alphak')
    call getmem3d(betak,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:betak')
    call getmem3d(coefe,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:coefe')
    call getmem3d(coeff1,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:coeff1')
    call getmem3d(coeff2,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:coeff2')
    call getmem3d(tpred1,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:tpred1')
    call getmem3d(tpred2,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:tpred2')
    call getmem3d(akzz1,jci1ga,jci2,ici1ga,ici2,1,kz,'mod_holtbl:akzz1')
    call getmem3d(akzz2,jci1ga,jci2,ici1ga,ici2,1,kz,'mod_holtbl:akzz2')
    call getmem2d(uvdrage,jci1ga,jci2,ici1ga,ici2,'mod_holtbl:uvdrage')
    if ( ichem == 1 ) then
      call getmem3d(kvc,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:kvc')
    end if
  end subroutine allocate_mod_pbl_holtbl

  subroutine holtbl(m2p,p2m)
    !@acc use nvtx
    implicit none
    type(mod_2_pbl), intent(in) :: m2p
    type(pbl_2_mod), intent(inout) :: p2m
    integer(ik4) :: i, j, k, n
    real(rkx) :: dudz, dvdz, ss, n2, rin, fofri, kzmax
    real(rkx) :: rrho, uflxsfx, vflxsfx, uu
    real(rkx) :: sh10, oblen
    real(rkx) :: zlv, tlv, ulv, vlv, zkv, tkv, vvk
    real(rkx) :: xfmt, wsc, therm, phpblm, zpbl, xfht
    real(rkx) :: z, zm, zp, zh, zl, wstr
    real(rkx) :: zzh, zzhnew, zzhnew2
    real(rkx) :: pblk, pblk1, pblk2, pr
    real(rkx) :: fak1, fak2, fak3, term
    real(rkx) :: drgdot, uflxsf, vflxsf
    real(rkx) :: coef1, coef2, coef3
    integer(ik4) :: iter
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'holtbl'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !@acc call nvtxStartRange("holtbl")
    !
    ! density at surface is stored in rhox2d
    ! the full level density is stored in rhohf.
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzm1 )
      dza(j,i,k) = m2p%za(j,i,k) - m2p%za(j,i,k+1)
      rhohf(j,i,k) = (m2p%patm(j,i,k+1)-m2p%patm(j,i,k)) / &
                     (egrav*dza(j,i,k))
    end do
    !
    ! Compute the g/dp term, virtual temperature and the
    !    theta->temperature function
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      hydf(j,i,k) = egrav/(m2p%patmf(j,i,k+1)-m2p%patmf(j,i,k))
      thvx(j,i,k) = m2p%thatm(j,i,k) * (d_one+ep1*m2p%qxatm(j,i,k,iqv))
      cfac(j,i,k) = m2p%tatm(j,i,k)/m2p%thatm(j,i,k)
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      exns(j,i) = (m2p%patmf(j,i,kzp1)/p00)**rovcp
      pfcor(j,i) = max(abs(m2p%coriol(j,i)),2.546e-5_rkx)
    end do
    !
    ! Compute the diffusion coefficient using Blackadar scheme above boundary
    ! layer top.
    ! The free atmosphere diffusivities are based on standard mixing length
    ! forms for the neutral diffusivity multiplied by functns of Richardson
    ! number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for
    ! momentum, potential temperature, and constitutents.
    ! The stable Richardson num function (Ri>0) is taken from Holtslag and
    ! Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))
    ! The unstable Richardson number function (Ri<0) is taken from  CCM1.
    ! f = sqrt(1 - 18*Ri)
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
      vv(j,i,k) = max(m2p%uxatm(j,i,k)*m2p%uxatm(j,i,k) + &
                      m2p%vxatm(j,i,k)*m2p%vxatm(j,i,k),0.025_rkx)
      dudz = (m2p%uxatm(j,i,k-1)-m2p%uxatm(j,i,k))/dza(j,i,k-1)
      dvdz = (m2p%vxatm(j,i,k-1)-m2p%vxatm(j,i,k))/dza(j,i,k-1)
      ! Vertical wind shear (square)
      ss = dudz**2 + dvdz**2 + 1.0e-10_rkx
      ! Brunt-Vaissala frequency
      n2 = egrav * (m2p%thatm(j,i,k-1)-m2p%thatm(j,i,k)) / &
          (dza(j,i,k-1)*0.5_rkx*(m2p%thatm(j,i,k-1)+m2p%thatm(j,i,k)))
      ! Compute the gradient Richardson number
      rin = max(-5.0_rkx,min(10.0_rkx,n2/ss))
      if ( rin < 0.0_rkx ) then
        fofri = sqrt(max(1.0_rkx-18.0_rkx*rin,0.0_rkx))
      else
        fofri = 1.0_rkx/(1.0_rkx+10.0_rkx*rin*(1.0_rkx+8.0_rkx*rin))
      end if
      kzm(j,i,k) = szkm*sqrt(ss)*fofri
      kzmax = kzfrac*dza(j,i,k-1)*m2p%dzq(j,i,k)*rdt
      kzm(j,i,k) = max(min(kzm(j,i,k),kzmax),kzo)
    end do
    !
    ! Holtslag pbl
    !
    ! Initialize bl diffusion coefficients and counter-gradient terms
    ! with free atmosphere values
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
      ! counter gradient terms for heat and moisture
      cgs(j,i,k) = 0.0_rkx
      cgh(j,i,k) = 0.0_rkx
      ! eddy diffusivities for momentum, heat and moisture
      kvm(j,i,k) = kzm(j,i,k)
      kvh(j,i,k) = kzm(j,i,k)
      kvq(j,i,k) = kzm(j,i,k)
      if ( ichem == 1 ) then
        kvc(j,i,k) = kzm(j,i,k)
      end if
    end do

    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      ! compute friction velocity
      rrho = 1.0_rkx/m2p%rhox2d(j,i)
      uflxsfx = -m2p%uvdrag(j,i)*m2p%uxatm(j,i,kz)*rrho
      vflxsfx = -m2p%uvdrag(j,i)*m2p%vxatm(j,i,kz)*rrho
      ! Minimum allowed ustr = 0.01
      uu = max(uflxsfx*uflxsfx+vflxsfx*vflxsfx,0.00000001_rkx)
      ustr(j,i) = sqrt(sqrt(uu))
      ! convert surface fluxes to kinematic units
      xhfx(j,i) = m2p%hfx(j,i)*rrho*rcpd
      xqfx(j,i) = m2p%qfx(j,i)*rrho
      ! Compute virtual heat flux at surface (surface kinematic buoyancy flux)
      hfxv(j,i) = xhfx(j,i) + ep1 * m2p%thatm(j,i,kz) * xqfx(j,i)
      lunstb(j,i) = (hfxv(j,i) > 0.0_rkx)
    end do
    !
    ! estimate potential temperature at 10m via log temperature
    ! profile in the surface layer (brutsaert, p. 63).
    ! calculate mixing ratio at 10m by assuming a constant
    ! value from the surface to the lowest model level.
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      ! Compute specific humidity
      sh10 = m2p%qxatm(j,i,kz,iqv)/(m2p%qxatm(j,i,kz,iqv)+d_one)
      ! "virtual" potential temperature
      if ( lunstb(j,i) ) then
        if ( m2p%ldmsk(j,i) == 1 ) then
          thv10(j,i) = m2p%tg(j,i) * (d_one + ep1*sh10)
        else
          thv10(j,i) = m2p%tg(j,i) * &
            (d_one + ep1*pfqsat(m2p%tg(j,i),m2p%patmf(j,i,kzp1)))
        end if
      else
        ! first approximation for obhukov length
        if ( ifaholtth10 == 1 ) then
          thv10(j,i) = (0.25_rkx*m2p%thatm(j,i,kz) + &
                       0.75_rkx*m2p%tg(j,i))*(d_one+ep1*sh10)
        else if ( ifaholtth10 == 2 ) then
          thv10(j,i) = thvx(j,i,kz) + hfxv(j,i)/(vonkar*ustr(j,i)* &
                     log(m2p%za(j,i,kz)*d_r10))
        else
          thv10(j,i) = (d_half*(m2p%thatm(j,i,kz)+m2p%tg(j,i))) * &
                       (d_one + ep1*sh10)
        end if
      end if
    end do
    do iter = 1, holtth10iter
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        oblen = -(thv10(j,i)*ustr(j,i)**3) / &
          (gvk*(hfxv(j,i)+sign(1.e-10_rkx,hfxv(j,i))))
        if ( oblen >= m2p%za(j,i,kz) ) then
          thv10(j,i) = thvx(j,i,kz) + hfxv(j,i)/(vonkar*ustr(j,i))*  &
             (log(m2p%za(j,i,kz)*d_r10)+d_five/oblen*(m2p%za(j,i,kz)-d_10))
        else
          if ( oblen < m2p%za(j,i,kz) .and. oblen > d_10 ) then
            thv10(j,i) = thvx(j,i,kz) + hfxv(j,i)/(vonkar*ustr(j,i))*  &
                (log(oblen*d_r10)+d_five/oblen*(oblen-d_10)+         &
                6.0_rkx*log(m2p%za(j,i,kz)/oblen))
          else
            thv10(j,i) = thvx(j,i,kz) + hfxv(j,i)/(vonkar*ustr(j,i)) * &
                        6.0_rkx*log(m2p%za(j,i,kz)*d_r10)
          end if
        end if
      end do
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      if ( ifaholt == 1 ) then
        thv10(j,i) = max(thv10(j,i),m2p%tg(j,i))  ! gtb add to maximize
      else  if ( ifaholt  == 2 ) then
        thv10(j,i) = min(thv10(j,i),m2p%tg(j,i))  ! gtb add to minimize
      end if
      ! obklen compute obukhov length
      obklen(j,i) = comp_obklen(thv10(j,i),ustr(j,i),hfxv(j,i))
    end do
    !
    ! compute diffusivities and counter gradient terms
    ! for momentum, heat and moisture and the counter-gradient
    ! terms for heat and moisture.
    !
    ! reference : holtslag, de bruijn and pan - mwr - 8/90
    !
    ! input arguments :  j       longitudinal position index
    !                    ubx3d   u wind component
    !                    vbx3d   v wind component
    !                    thatm   potential temperature
    !                    thvx    virtual potential temperature
    !                    za      height of half sigma levels
    !                    f       coriolis parameter
    !                    shum    specific humidity
    !                    xhfx    sensible heat flux
    !                    xqfx    sfc kinematic moisture flux
    !                    thv10   virt. pot. temp. at 10m
    !                    hfxv    surface virtual heat flux
    !                    obklen  monin obukov length
    !                    ustr    friction velocity
    !
    ! input/output
    ! arguments :        therm   thermal temperature excess
    !
    ! output arguments : cgh     counter-gradient term for heat
    !                    cgs     counter-gradient star term
    !                    kvm     eddy diffusivity for momentum
    !                    kvh     eddy diffusivity for heat
    !                    kvq     eddy diffusivity for moisture
    !                   zpbl     boundary layer height
    ! ------------------------------------------------------------
    !
    !
    ! note: kmxpbl, max no. of pbl levels (set in slice)
    ! compute Bulk Richardson Number (BRN)
    !
    if ( idynamic == 3 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        zlv = m2p%za(j,i,kz)
        tlv = thv10(j,i)
        ulv = m2p%uxatm(j,i,kz)
        vlv = m2p%vxatm(j,i,kz)
        do k = kzm1, kmxpbl(j,i), -1
          zkv = m2p%za(j,i,k)
          tkv = thvx(j,i,k)
          vvk = (m2p%uxatm(j,i,k)-ulv)**2+(m2p%vxatm(j,i,k)-vlv)**2
          vvk = vvk + 1.0e-10_rkx
          ri(k,j,i) = egrav*(tkv-tlv)*(zkv-zlv)/(tlv*vvk)
          ri(k,j,i) = max(-5.0_rkx,min(10.0_rkx,ri(k,j,i)))
        end do
      end do
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        do k = kzm1, kmxpbl(j,i), -1
          ri(k,j,i) = egrav*(thvx(j,i,k)-thv10(j,i))*m2p%za(j,i,k) / &
                      (thv10(j,i)*vv(j,i,k))
          ri(k,j,i) = max(-5.0_rkx,min(10.0_rkx,ri(k,j,i)))
        end do
      end do
    end if

    ! looking for first guess bl top
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      p2m%zpbl(j,i) = m2p%za(j,i,kz)
      do k = kzm1, kmxpbl(j,i) + 1, -1
        ! bl height lies between this level and the last
        ! use linear interp. of rich. no. to height of ri=ricr
        if ( (ri(k,j,i)   <  ricr(j,i)) .and. &
             (ri(k-1,j,i) >= ricr(j,i)) ) then
          p2m%zpbl(j,i) = m2p%za(j,i,k)+(m2p%za(j,i,k-1)-m2p%za(j,i,k)) * &
              ((ricr(j,i)-ri(k,j,i))/(ri(k-1,j,i)-ri(k,j,i)))
        end if
      end do
    end do

    ! recompute richardson no. at lowest model level
    if ( idynamic == 3 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        if ( lunstb(j,i) ) then
          ! estimate of convective velocity scale
          xfmt = (d_one-(binm*p2m%zpbl(j,i)/obklen(j,i)))**onet
          wsc = ustr(j,i)*xfmt
          ! thermal temperature excess
          therm = fak * hfxv(j,i)/wsc
          zlv = m2p%za(j,i,kz)
          tlv = thv10(j,i) + therm
          ulv = m2p%uxatm(j,i,kz)
          vlv = m2p%vxatm(j,i,kz)
          !zlv = max(obklen(j,i),d_10)
          !tlv = thv10(j,i) + therm
          vvk = ulv**2 + vlv**2 + fak*ustr(j,i)**2 + 1.0e-10_rkx
          ri(kz,j,i) = -egrav*therm*zlv/(thv10(j,i)*vvk)
          ri(kz,j,i) = max(-5.0_rkx,min(10.0_rkx,ri(kz,j,i)))
          ! recompute richardson no. at other model levels
          do k = kzm1, kmxpbl(j,i), -1
            zkv = m2p%za(j,i,k)
            tkv = thvx(j,i,k)
            vvk = (m2p%uxatm(j,i,k)-ulv)**2+(m2p%vxatm(j,i,k)-vlv)**2
            vvk = vvk + 1.0e-10_rkx
            ri(k,j,i) = egrav*(tkv-tlv)*(zkv-zlv)/(thv10(j,i)*vvk)
            ri(k,j,i) = max(-5.0_rkx,min(10.0_rkx,ri(k,j,i)))
          end do
        end if
      end do
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        if ( lunstb(j,i) ) then
          ! estimate of convective velocity scale
          xfmt = (d_one-(binm*p2m%zpbl(j,i)/obklen(j,i)))**onet
          wsc = ustr(j,i)*xfmt
          ! thermal temperature excess
          therm = fak * hfxv(j,i)/wsc
          tlv = thv10(j,i) + therm
          ri(kz,j,i) = -egrav*therm*m2p%za(j,i,kz)/(thv10(j,i)*vv(j,i,kz))
          ri(kz,j,i) = max(-5.0_rkx,min(10.0_rkx,ri(kz,j,i)))
          ! recompute richardson no. at other model levels
          do k = kzm1, kmxpbl(j,i), -1
            tkv = thvx(j,i,k)
            ri(k,j,i) = egrav*(tkv-tlv)*m2p%za(j,i,k) / &
               (thv10(j,i)*vv(j,i,k))
            ri(k,j,i) = max(-5.0_rkx,min(10.0_rkx,ri(k,j,i)))
          end do
        end if
      end do
    end if

    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      if ( lunstb(j,i) ) then
        ! improve estimate of bl height under convective conditions
        ! using convective temperature excess (therm)
        do k = kz, kmxpbl(j,i) + 1, -1
          ! bl height lies between this level and the last
          ! use linear interp. of rich. no. to height of ri=ricr
          if ( (ri(k,j,i) < ricr(j,i)) .and. &
               (ri(k-1,j,i) >= ricr(j,i)) ) then
            p2m%zpbl(j,i) = m2p%za(j,i,k) + &
              (m2p%za(j,i,k-1)-m2p%za(j,i,k))* &
              ((ricr(j,i)-ri(k,j,i))/(ri(k-1,j,i)-ri(k,j,i)))
          end if
        end do
      end if
    end do

    ! Find the k of the level of the pbl
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      ! BL height is at least the mechanical mixing depth
      ! PBL height must be greater than some minimum mechanical mixing depth
      ! Several investigators have proposed minimum mechanical mixing depth
      ! relationships as a function of the local friction velocity, u*.  We
      ! make use of a linear relationship of the form h = c u* where c=700.
      ! The scaling arguments that give rise to this relationship most often
      ! represent the coefficient c as some constant over the local coriolis
      ! parameter.  Here we make use of the experimental results of Koracin
      ! and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
      ! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
      ! latitude value for f so that c = 0.07/f = 700.
      !phpblm = 700.0_rkx*ustr(j,i)
      phpblm = (0.07_rkx*ustr(j,i))/pfcor(j,i)
      if ( p2m%zpbl(j,i) < phpblm ) then
        p2m%zpbl(j,i) = max(phpblm,p2m%zpbl(j,i))
      end if
      do k = kz, kmxpbl(j,i), -1
        p2m%kpbl(j,i) = k
        if ( m2p%za(j,i,k) > p2m%zpbl(j,i) ) exit
      end do
    end do

    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      zpbl = p2m%zpbl(j,i)
      fak1 = ustr(j,i)*zpbl*vonkar
      if ( lunstb(j,i) ) then
        xfmt = (d_one-binm*zpbl/obklen(j,i))**onet
        xfht = sqrt(d_one-binh*zpbl/obklen(j,i))
        wsc = ustr(j,i)*xfmt
        fak2 = wsc*zpbl*vonkar
      else
        xfmt = d_zero
        xfht = d_zero
        wsc = d_zero
        fak2 = d_zero
      end if
      do k = kz, p2m%kpbl(j,i), -1
        zm = m2p%za(j,i,k)
        zp = m2p%za(j,i,k-1)
        if ( zm < zpbl ) then
          zp = min(zp,zpbl)
          z = (zm+zp)*d_half
          zh = z/zpbl
          zl = z/obklen(j,i)
          term = max(d_one-zh,0.0_rkx)
          zzh = zh*term**pink
          zzhnew = zzh*zhnew_fac
          zzhnew2 = zzh*zhnew_fac
          if ( lunstb(j,i) ) then
            if ( zh < sffrac ) then
              term = (d_one-betam*zl)**onet
              pblk = fak1*zzh*term
              pblk1 = fak1*zzhnew*term
              pblk2 = fak1*zzhnew2*term
              pr = term/sqrt(d_one-betah*zl)
            else
              ! Convective velocity scale
              wstr = (hfxv(j,i)*egrav*zpbl/thv10(j,i))**onet
              fak3 = fakn*wstr/wsc
              pblk = fak2*zzh
              pblk1 = fak2*zzhnew
              pblk2 = fak2*zzhnew2
              ! compute counter gradient term
              pr = (xfmt/xfht) + ccon*fak3/fak
              cgs(j,i,k) = fak3/(zpbl*wsc)
              cgh(j,i,k) = xhfx(j,i)*cgs(j,i,k)
           end if
          else
            if ( zl < d_one ) then
              pblk = fak1*zzh/(d_one+betas*zl)
              pblk1 = fak1*zzhnew/(d_one+betas*zl)
              pblk2 = fak1*zzhnew2/(d_one+betas*zl)
            else
              pblk = fak1*zzh/(betas+zl)
              pblk1 = fak1*zzhnew/(betas+zl)
              pblk2 = fak1*zzhnew2/(betas+zl)
            end if
            pr = 1.0_rkx
          end if
          ! compute eddy diffusivities
          kvm(j,i,k) = max(pblk,kvm(j,i,k))
          kvh(j,i,k) = max(pblk/pr,kvh(j,i,k))
          kvq(j,i,k) = max(pblk1,kvq(j,i,k))
          if ( ichem == 1 ) then
            kvc(j,i,k) = max(pblk2,kvc(j,i,k))
          end if
        end if
      end do
    end do

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
      akzz1(j,i,k) = rhohf(j,i,k-1)*kvm(j,i,k)/dza(j,i,k-1)
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      akzz2(j,i,k) = hydf(j,i,k)
    end do

    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      uvdrage(j,i) = m2p%uvdrag(j,i)
    end do

    call exchange_lb(akzz1,1,jci1,jci2,ici1,ici2,1,kz)
    call exchange_lb(akzz2,1,jci1,jci2,ici1,ici2,1,kz)
    call exchange_lb(uvdrage,1,jci1,jci2,ici1,ici2)

    if ( idynamic == 3 ) then
      !
      !   calculate coefficients at dot points for u
      !
      do concurrent ( j = jdii1:jdii2, i = ici1:ici2, k = 2:kz )
        betak(j,i,k) = 0.5_rkx * (akzz1(j-1,i,k)+akzz1(j,i,k))
      end do
      do concurrent ( j = jdii1:jdii2, i = ici1:ici2, k = 1:kz )
        alphak(j,i,k) = 0.5_rkx * (akzz2(j-1,i,k)+akzz2(j,i,k))
      end do
      !
      ! start now procedure for implicit diffusion calculations
      ! performed separately for wind (dot points)
      ! and temperature and water vapor (cross points)
      ! countergradient term is not included in the implicit diffusion
      ! scheme its effect is included as in the old explicit scheme
      ! calculations assume fluxes positive upward, so the sign in front
      ! of uflxsf and vflxsf has been changed in the various terms
      !
      ! wind components
      !
      ! first compute coefficients of the tridiagonal matrix
      !
      ! Atmosphere top
#ifdef STDPAR_FIXED
      do concurrent ( j = jdii1:jdii2, i = ici1:ici2 )
#else
      !$acc parallel loop collapse(2) gang vector
      do i = ici1, ici2
      do j = jdii1, jdii2
#endif
        coef1 = dt*alphak(j,i,1)*betak(j,i,2)
        coef2 = d_one + dt*alphak(j,i,1)*betak(j,i,2)
        coef3 = d_zero
        coefe(j,i,1) = coef1/coef2
        coeff1(j,i,1) = m2p%udatm(j,i,1)/coef2
        ! top to bottom
        !$acc loop seq
        do k = 2, kzm1
          coef1 = dt*alphak(j,i,k)*betak(j,i,k+1)
          coef2 = d_one+dt*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
          coef3 = dt*alphak(j,i,k)*betak(j,i,k)
          coefe(j,i,k) = coef1/(coef2-coef3*coefe(j,i,k-1))
          coeff1(j,i,k) = (m2p%udatm(j,i,k) + &
               coef3*coeff1(j,i,k-1))/(coef2-coef3*coefe(j,i,k-1))
        end do
        ! Nearest to surface
        coef1 = d_zero
        coef2 = d_one + dt*alphak(j,i,kz)*betak(j,i,kz)
        coef3 = dt*alphak(j,i,kz)*betak(j,i,kz)
        drgdot = 0.5_rkx * (uvdrage(j-1,i)+uvdrage(j,i))
        uflxsf = drgdot*m2p%udatm(j,i,kz)
        coefe(j,i,kz) = d_zero
        coeff1(j,i,kz) = (m2p%udatm(j,i,kz) - dt*alphak(j,i,kz)*uflxsf + &
             coef3*coeff1(j,i,kz-1))/(coef2-coef3*coefe(j,i,kz-1))
      !
      !   all coefficients have been computed, predict field and put it in
      !   temporary work space tpred
      !
        tpred1(j,i,kz) = coeff1(j,i,kz)
        !$acc loop seq
        do k = kzm1, 1, -1
          tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
        end do
#ifndef STDPAR_FIXED
      end do
#endif
      end do
      !
      !   calculate tendency due to vertical diffusion using temporary
      !   predicted field
      !
      do concurrent ( j = jdii1:jdii2, i = ici1:ici2, k = 1:kz )
        p2m%uten(j,i,k) = p2m%uten(j,i,k) + &
                          (tpred1(j,i,k)-m2p%udatm(j,i,k))*rdt
      end do
      !
      !   calculate coefficients at dot points for v wind
      !
      do concurrent ( j = jci1:jci2, i = idii1:idii2, k = 2:kz )
        betak(j,i,k) = 0.5_rkx * (akzz1(j,i-1,k)+akzz1(j,i,k))
      end do
      do concurrent ( j = jci1:jci2, i = idii1:idii2, k = 1:kz )
        alphak(j,i,k) = 0.5_rkx * (akzz2(j,i-1,k)+akzz2(j,i,k))
      end do
      !
      ! start now procedure for implicit diffusion calculations
      ! performed separately for wind (dot points)
      ! and temperature and water vapor (cross points)
      ! countergradient term is not included in the implicit diffusion
      ! scheme its effect is included as in the old explicit scheme
      ! calculations assume fluxes positive upward, so the sign in front
      ! of uflxsf and vflxsf has been changed in the various terms
      !
      ! wind components
      !
      ! first compute coefficients of the tridiagonal matrix
      !
      ! Atmosphere top
#ifdef STDPAR_FIXED
      do concurrent ( j = jci1:jci2, i = idii1:idii2 )
#else
      !$acc parallel loop collapse(2) gang vector
      do i = idii1, idii2
      do j = jci1, jci2
#endif
        coef1 = dt*alphak(j,i,1)*betak(j,i,2)
        coef2 = d_one + dt*alphak(j,i,1)*betak(j,i,2)
        coef3 = d_zero
        coefe(j,i,1) = coef1/coef2
        coeff2(j,i,1) = m2p%vdatm(j,i,1)/coef2
        ! top to bottom
        !$acc loop seq
        do k = 2, kzm1
          coef1 = dt*alphak(j,i,k)*betak(j,i,k+1)
          coef2 = d_one+dt*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
          coef3 = dt*alphak(j,i,k)*betak(j,i,k)
          coefe(j,i,k) = coef1/(coef2-coef3*coefe(j,i,k-1))
          coeff2(j,i,k) = (m2p%vdatm(j,i,k) + &
               coef3*coeff2(j,i,k-1))/(coef2-coef3*coefe(j,i,k-1))
        end do
        ! Nearest to surface
        coef1 = d_zero
        coef2 = d_one + dt*alphak(j,i,kz)*betak(j,i,kz)
        coef3 = dt*alphak(j,i,kz)*betak(j,i,kz)
        drgdot = 0.5_rkx * (uvdrage(j,i-1)+uvdrage(j,i))
        vflxsf = drgdot*m2p%vdatm(j,i,kz)
        coefe(j,i,kz) = d_zero
        coeff2(j,i,kz) = (m2p%vdatm(j,i,kz) - dt*alphak(j,i,kz)*vflxsf + &
             coef3*coeff2(j,i,kz-1))/(coef2-coef3*coefe(j,i,kz-1))
      !
      !   all coefficients have been computed, predict field and put it in
      !   temporary work space tpred
      !
        tpred2(j,i,kz) = coeff2(j,i,kz)
        !$acc loop seq
        do k = kzm1, 1, -1
          tpred2(j,i,k) = coefe(j,i,k)*tpred2(j,i,k+1) + coeff2(j,i,k)
        end do
#ifndef STDPAR_FIXED
      end do
#endif
      end do
      !
      !   calculate tendency due to vertical diffusion using temporary
      !   predicted field
      !
      do concurrent ( j = jci1:jci2, i = idii1:idii2, k = 1:kz )
        p2m%vten(j,i,k) = p2m%vten(j,i,k) + &
                          (tpred2(j,i,k)-m2p%vdatm(j,i,k))*rdt
      end do

    else
      !
      !   calculate coefficients at dot points for u and v wind
      !
      do concurrent ( j = jdii1:jdii2, i = idii1:idii2, k = 2:kz )
        betak(j,i,k) = (akzz1(j-1,i,k)+akzz1(j-1,i-1,k)+ &
                        akzz1(j,i,k)  +akzz1(j,i-1,k))*d_rfour
      end do
      do concurrent ( j = jdii1:jdii2, i = idii1:idii2, k = 1:kz )
        alphak(j,i,k) = (akzz2(j-1,i,k)+akzz2(j-1,i-1,k)+ &
                         akzz2(j,i,k)  +akzz2(j,i-1,k))*d_rfour
      end do
      !
      ! start now procedure for implicit diffusion calculations
      ! performed separately for wind (dot points)
      ! and temperature and water vapor (cross points)
      ! countergradient term is not included in the implicit diffusion
      ! scheme its effect is included as in the old explicit scheme
      ! calculations assume fluxes positive upward, so the sign in front
      ! of uflxsf and vflxsf has been changed in the various terms
      !
      ! wind components
      !
      ! first compute coefficients of the tridiagonal matrix
      !
      ! Atmosphere top
#ifdef STDPAR_FIXED
      do concurrent ( j = jdii1:jdii2, i = idii1:idii2 )
#else
      !$acc parallel loop collapse(2) gang vector
      do i = idii1, idii2
      do j = jdii1, jdii2
#endif
        coef1 = dt*alphak(j,i,1)*betak(j,i,2)
        coef2 = d_one + dt*alphak(j,i,1)*betak(j,i,2)
        coef3 = d_zero
        coefe(j,i,1) = coef1/coef2
        coeff1(j,i,1) = m2p%udatm(j,i,1)/coef2
        coeff2(j,i,1) = m2p%vdatm(j,i,1)/coef2
        ! top to bottom
        !$acc loop seq
        do k = 2, kzm1
          coef1 = dt*alphak(j,i,k)*betak(j,i,k+1)
          coef2 = d_one+dt*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
          coef3 = dt*alphak(j,i,k)*betak(j,i,k)
          coefe(j,i,k) = coef1/(coef2-coef3*coefe(j,i,k-1))
          coeff1(j,i,k) = (m2p%udatm(j,i,k) + &
               coef3*coeff1(j,i,k-1))/(coef2-coef3*coefe(j,i,k-1))
          coeff2(j,i,k) = (m2p%vdatm(j,i,k) + &
               coef3*coeff2(j,i,k-1))/(coef2-coef3*coefe(j,i,k-1))
        end do
        ! Nearest to surface
        coef1 = d_zero
        coef2 = d_one + dt*alphak(j,i,kz)*betak(j,i,kz)
        coef3 = dt*alphak(j,i,kz)*betak(j,i,kz)
        drgdot = (uvdrage(j-1,i-1)+uvdrage(j,i-1) + &
                  uvdrage(j-1,i)  +uvdrage(j,i))*d_rfour
        uflxsf = drgdot*m2p%udatm(j,i,kz)
        vflxsf = drgdot*m2p%vdatm(j,i,kz)
        coefe(j,i,kz) = d_zero
        coeff1(j,i,kz) = (m2p%udatm(j,i,kz) - dt*alphak(j,i,kz)*uflxsf + &
             coef3*coeff1(j,i,kz-1))/(coef2-coef3*coefe(j,i,kz-1))
        coeff2(j,i,kz) = (m2p%vdatm(j,i,kz) - dt*alphak(j,i,kz)*vflxsf + &
             coef3*coeff2(j,i,kz-1))/(coef2-coef3*coefe(j,i,kz-1))
      !
      !   all coefficients have been computed, predict field and put it in
      !   temporary work space tpred
      !
        tpred1(j,i,kz) = coeff1(j,i,kz)
        tpred2(j,i,kz) = coeff2(j,i,kz)
        !$acc loop seq
        do k = kzm1, 1, -1
          tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
          tpred2(j,i,k) = coefe(j,i,k)*tpred2(j,i,k+1) + coeff2(j,i,k)
        end do
#ifndef STDPAR_FIXED
      end do
#endif
      end do
      !
      !   calculate tendency due to vertical diffusion using temporary
      !   predicted field
      !
      do concurrent ( j = jdii1:jdii2, i = idii1:idii2, k = 1:kz )
        p2m%uten(j,i,k) = p2m%uten(j,i,k) + &
                          (tpred1(j,i,k)-m2p%udatm(j,i,k))*rdt*m2p%psdotb(j,i)
        p2m%vten(j,i,k) = p2m%vten(j,i,k) + &
                          (tpred2(j,i,k)-m2p%vdatm(j,i,k))*rdt*m2p%psdotb(j,i)
      end do
    end if
    !
    !   Common coefficients.
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      alphak(j,i,k) = hydf(j,i,k)
    end do
    !
    !   temperature
    !   calculate coefficients at cross points for temperature
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
      betak(j,i,k) = rhohf(j,i,k-1)*kvh(j,i,k)/dza(j,i,k-1)
    end do

#ifdef STDPAR_FIXED
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
#else
    !$acc parallel loop collapse(2) gang vector
    do i = ici1, ici2
    do j = jci1, jci2
#endif
      coef1 = dt*alphak(j,i,1)*betak(j,i,2)
      coef2 = d_one + dt*alphak(j,i,1)*betak(j,i,2)
      coef3 = d_zero
      coefe(j,i,1) = coef1/coef2
      coeff1(j,i,1) = m2p%thatm(j,i,1)/coef2
      !$acc loop seq
      do k = 2, kzm1
        coef1 = dt*alphak(j,i,k)*betak(j,i,k+1)
        coef2 = d_one+dt*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
        coef3 = dt*alphak(j,i,k)*betak(j,i,k)
        coefe(j,i,k) = coef1/(coef2-coef3*coefe(j,i,k-1))
        coeff1(j,i,k) = (m2p%thatm(j,i,k) + &
             coef3*coeff1(j,i,k-1))/(coef2-coef3*coefe(j,i,k-1))
      end do
      coef1 = d_zero
      coef2 = d_one + dt*alphak(j,i,kz)*betak(j,i,kz)
      coef3 = dt*alphak(j,i,kz)*betak(j,i,kz)
      coefe(j,i,kz) = d_zero
      coeff1(j,i,kz) = (m2p%thatm(j,i,kz) + &
           dt*alphak(j,i,kz)*rcpd*m2p%hfx(j,i) + &
           coef3*coeff1(j,i,kz-1))/(coef2-coef3*coefe(j,i,kz-1))
    !
    !   all coefficients have been computed, predict field and put it in
    !   temporary work space tpred
    !
      tpred1(j,i,kz) = coeff1(j,i,kz)
      !$acc loop seq
      do k = kzm1, 1, -1
        tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
      end do
#ifndef STDPAR_FIXED
    end do
#endif
    end do
    !
    !   calculate tendency due to vertical diffusion using temporary
    !   predicted field
    !
    if ( idynamic == 3 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        p2m%tten(j,i,k) = p2m%tten(j,i,k) + &
                 (tpred1(j,i,k)-m2p%thatm(j,i,k))*cfac(j,i,k)*rdt
      end do
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        p2m%tten(j,i,k) = p2m%tten(j,i,k) + &
                 (tpred1(j,i,k)-m2p%thatm(j,i,k))*cfac(j,i,k)*rdt*m2p%psb(j,i)
      end do
    end if
    !
    !   water vapor calculate coefficients at cross points for water vapor
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
      betak(j,i,k) = rhohf(j,i,k-1)*kvq(j,i,k)/dza(j,i,k-1)
    end do

#ifdef STDPAR_FIXED
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
#else
    !$acc parallel loop collapse(2) gang vector
    do i = ici1, ici2
    do j = jci1, jci2
#endif
      coef1 = dt*alphak(j,i,1)*betak(j,i,2)
      coef2 = d_one + dt*alphak(j,i,1)*betak(j,i,2)
      coef3 = d_zero
      coefe(j,i,1) = coef1/coef2
      coeff1(j,i,1) = m2p%qxatm(j,i,1,iqv)/coef2
      !$acc loop seq
      do k = 2, kzm1
        coef1 = dt*alphak(j,i,k)*betak(j,i,k+1)
        coef2 = d_one+dt*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
        coef3 = dt*alphak(j,i,k)*betak(j,i,k)
        coefe(j,i,k) = coef1/(coef2-coef3*coefe(j,i,k-1))
        coeff1(j,i,k) = (m2p%qxatm(j,i,k,iqv) +         &
             coef3*coeff1(j,i,k-1))/(coef2-coef3*coefe(j,i,k-1))
      end do
      coef1 = d_zero
      coef2 = d_one + dt*alphak(j,i,kz)*betak(j,i,kz)
      coef3 = dt*alphak(j,i,kz)*betak(j,i,kz)
      coefe(j,i,kz) = d_zero
      coeff1(j,i,kz) = (m2p%qxatm(j,i,kz,iqv) +  &
           dt*alphak(j,i,kz)*m2p%qfx(j,i) +  &
           coef3*coeff1(j,i,kz-1))/(coef2-coef3*coefe(j,i,kz-1))
    !
    !   all coefficients have been computed, predict field and put it in
    !   temporary work space tpred
    !
      tpred1(j,i,kz) = coeff1(j,i,kz)
      !$acc loop seq
      do k = kzm1, 1, -1
        tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
      end do
#ifndef STDPAR_FIXED
    end do
#endif
    end do
    !
    !   calculate tendency due to vertical diffusion using temporary
    !   predicted field
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      qtenv(j,i,k) = d_zero
    end do
    if ( idynamic == 3 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        qtenv(j,i,k) = (tpred1(j,i,k)-m2p%qxatm(j,i,k,iqv))*rdt
      end do
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        qtenv(j,i,k) = (tpred1(j,i,k)-m2p%qxatm(j,i,k,iqv))*rdt*m2p%psb(j,i)
      end do
    end if
    !
    !   calculate coefficients at cross points for cloud water
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
      betak(j,i,k) = rhohf(j,i,k-1)*kvq(j,i,k)/dza(j,i,k-1)
    end do

#ifdef STDPAR_FIXED
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
#else
    !$acc parallel loop collapse(2) gang vector
    do i = ici1, ici2
    do j = jci1, jci2
#endif
      coef1 = dt*alphak(j,i,1)*betak(j,i,2)
      coef2 = d_one + dt*alphak(j,i,1)*betak(j,i,2)
      coef3 = d_zero
      coefe(j,i,1) = coef1/coef2
      coeff1(j,i,1) = m2p%qxatm(j,i,1,iqc)/coef2
      !$acc loop seq
      do k = 2, kzm1
        coef1 = dt*alphak(j,i,k)*betak(j,i,k+1)
        coef2 = d_one+dt*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
        coef3 = dt*alphak(j,i,k)*betak(j,i,k)
        coefe(j,i,k) = coef1/(coef2-coef3*coefe(j,i,k-1))
        coeff1(j,i,k) = (m2p%qxatm(j,i,k,iqc) + &
             coef3*coeff1(j,i,k-1))/(coef2-coef3*coefe(j,i,k-1))
      end do
      coef1 = d_zero
      coef2 = d_one + dt*alphak(j,i,kz)*betak(j,i,kz)
      coef3 = dt*alphak(j,i,kz)*betak(j,i,kz)
      coefe(j,i,kz) = d_zero
      coeff1(j,i,kz) = (m2p%qxatm(j,i,kz,iqc) + &
           coef3*coeff1(j,i,kz-1))/(coef2-coef3*coefe(j,i,kz-1))
    !
    !   all coefficients have been computed, predict field and put it in
    !   temporary work space tpred
    !
      tpred1(j,i,kz) = coeff1(j,i,kz)
      !$acc loop seq
      do k = kzm1, 1, -1
        tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
      end do
#ifndef STDPAR_FIXED
    end do
#endif
    end do
    !
    !   calculate tendency due to vertical diffusion using temporary
    !   predicted field
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      qtenc(j,i,k) = d_zero
    end do
    if ( idynamic == 3 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        qtenc(j,i,k) = (tpred1(j,i,k)-m2p%qxatm(j,i,k,iqc))*rdt
      end do
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        qtenc(j,i,k) = (tpred1(j,i,k)-m2p%qxatm(j,i,k,iqc))*rdt*m2p%psb(j,i)
      end do
    end if

    if ( ipptls > 1 ) then
#ifdef STDPAR_FIXED
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
#else
      !$acc parallel loop collapse(2) gang vector
      do i = ici1, ici2
      do j = jci1, jci2
#endif
        coef1 = dt*alphak(j,i,1)*betak(j,i,2)
        coef2 = d_one + dt*alphak(j,i,1)*betak(j,i,2)
        coef3 = d_zero
        coefe(j,i,1) = coef1/coef2
        coeff1(j,i,1) = m2p%qxatm(j,i,1,iqi)/coef2
        !$acc loop seq
        do k = 2, kzm1
          coef1 = dt*alphak(j,i,k)*betak(j,i,k+1)
          coef2 = d_one+dt*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
          coef3 = dt*alphak(j,i,k)*betak(j,i,k)
          coefe(j,i,k) = coef1/(coef2-coef3*coefe(j,i,k-1))
          coeff1(j,i,k) = (m2p%qxatm(j,i,k,iqi) + &
               coef3*coeff1(j,i,k-1))/(coef2-coef3*coefe(j,i,k-1))
        end do
        coef1 = d_zero
        coef2 = d_one + dt*alphak(j,i,kz)*betak(j,i,kz)
        coef3 = dt*alphak(j,i,kz)*betak(j,i,kz)
        coefe(j,i,kz) = d_zero
        coeff1(j,i,kz) = (m2p%qxatm(j,i,kz,iqi) + &
             coef3*coeff1(j,i,kz-1))/(coef2-coef3*coefe(j,i,kz-1))
      !
      !   all coefficients have been computed, predict field and put it in
      !   temporary work space tpred
      !
        tpred1(j,i,kz) = coeff1(j,i,kz)
        !$acc loop seq
        do k = kzm1, 1, -1
          tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
        end do
#ifndef STDPAR_FIXED
      end do
#endif
      end do
      !
      !   calculate tendency due to vertical diffusion using temporary
      !   predicted field
      !
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        qteni(j,i,k) = d_zero
      end do
      if ( idynamic == 3 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          qteni(j,i,k) = (tpred1(j,i,k)-m2p%qxatm(j,i,k,iqi))*rdt
        end do
      else
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          qteni(j,i,k) = (tpred1(j,i,k)-m2p%qxatm(j,i,k,iqi))*rdt*m2p%psb(j,i)
        end do
      end if
    end if
    !
    !  now add countergradient term to temperature and water vapor equation
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      ttnp(j,i,1) = d_zero
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
      ttnp(j,i,k) = cfac(j,i,k)*hydf(j,i,k)*rhohf(j,i,k-1) * &
                        kvh(j,i,k)*cgh(j,i,k)
    end do
    if ( idynamic /= 3 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
        ttnp(j,i,k) = ttnp(j,i,k) * m2p%psb(j,i)
      end do
    end if
    !
    !   compute the tendencies:
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzm1 )
      p2m%tten(j,i,k) = p2m%tten(j,i,k)+(ttnp(j,i,k+1)-ttnp(j,i,k))
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      p2m%tten(j,i,kz) = p2m%tten(j,i,kz) - ttnp(j,i,kz)
    end do
    if ( idynamic == 3 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        ttnp(j,i,1) = d_zero
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
        ttnp(j,i,k) = hydf(j,i,k)*rhohf(j,i,k-1) * &
                          kvh(j,i,k)*cgs(j,i,k)*xqfx(j,i)
      end do
      !
      !   compute the tendencies:
      !
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzm1 )
        qtenv(j,i,k) = qtenv(j,i,k) + (ttnp(j,i,k+1)-ttnp(j,i,k))
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        qtenv(j,i,kz) = qtenv(j,i,kz) - ttnp(j,i,kz)
      end do
    end if

#ifdef RCEMIP
    !call force_water_conserve(qtenv,qtenc,qteni,m2p%qxatm,xqfx,m2p%psb)
#endif
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      p2m%qxten(j,i,k,iqv) = p2m%qxten(j,i,k,iqv) + qtenv(j,i,k)
      p2m%qxten(j,i,k,iqc) = p2m%qxten(j,i,k,iqc) + qtenc(j,i,k)
    end do
    if ( ipptls > 1 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        p2m%qxten(j,i,k,iqi) = p2m%qxten(j,i,k,iqi) + qteni(j,i,k)
      end do
    end if

    if ( ichem == 1 ) then
      !
      !     coef1, coef2, coef3 and coefe are the same as for water vapor
      !     and cloud water so they do not need to be recalculated
      !     recalculation of coef1,2,3  with tracer diffusivity kvc
      !
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
        betak(j,i,k) = rhohf(j,i,k-1)*kvc(j,i,k)/dza(j,i,k-1)
      end do

      do n = 1, ntr
#ifdef STDPAR_FIXED
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
#else
        !$acc parallel loop collapse(2) gang vector
        do i = ici1, ici2
        do j = jci1, jci2
#endif
          coef1 = dt*alphak(j,i,1)*betak(j,i,2)
          coef2 = d_one + dt*alphak(j,i,1)*betak(j,i,2)
          coef3 = d_zero
          coefe(j,i,1) = coef1/coef2
          coeff1(j,i,1) = m2p%chib(j,i,1,n)/coef2
          if ( abs(coeff1(j,i,1)) < dlowval ) coeff1(j,i,1) = d_zero
          if ( abs(coefe(j,i,1)) < dlowval ) coefe(j,i,1) = d_zero
          !$acc loop seq
          do k = 2, kzm1
            coef1 = dt*alphak(j,i,k)*betak(j,i,k+1)
            coef2 = d_one+dt*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
            coef3 = dt*alphak(j,i,k)*betak(j,i,k)
            coefe(j,i,k) = coef1/(coef2-coef3*coefe(j,i,k-1))
            coeff1(j,i,k) = (m2p%chib(j,i,k,n) + &
                 coef3*coeff1(j,i,k-1))/(coef2-coef3*coefe(j,i,k-1))
            if ( abs(coeff1(j,i,k)) < dlowval ) coeff1(j,i,k) = d_zero
            if ( abs(coefe(j,i,k)) < dlowval ) coefe(j,i,k) = d_zero
          end do
          coef1 = d_zero
          coef2 = d_one + dt*alphak(j,i,kz)*betak(j,i,kz)
          coef3 = dt*alphak(j,i,kz)*betak(j,i,kz)
          coefe(j,i,kz) = d_zero
          ! add dry deposition option1
          coeff1(j,i,kz) = (m2p%chib(j,i,kz,n)-dt*alphak(j,i,kz) * &
               m2p%chib(j,i,kz,n)*m2p%drydepv(j,i,n)*m2p%rhox2d(j,i) + &
               coef3*coeff1(j,i,kz-1))/(coef2-coef3*coefe(j,i,kz-1))
          if ( abs(coeff1(j,i,kz)) < dlowval ) coeff1(j,i,kz) = d_zero
          !
          !     all coefficients have been computed, predict field and put
          !     it in temporary work space tpred1
          !
          tpred1(j,i,kz) = coeff1(j,i,kz)
          !$acc loop seq
          do k = kzm1, 1, -1
            tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
          end do
#ifndef STDPAR_FIXED
        end do
#endif
        end do
        !
        !       calculate tendency due to vertical diffusion using temporary
        !       predicted field
        !       Dry deposition option 1 is included
        !
        if ( idynamic == 3 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            p2m%chiten(j,i,k,n) = p2m%chiten(j,i,k,n) +  &
                          (tpred1(j,i,k)-m2p%chib(j,i,k,n))*rdt
          end do
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            p2m%remdrd(j,i,n) = p2m%remdrd(j,i,n) + &
                m2p%chib(j,i,kz,n)*m2p%drydepv(j,i,n) * &
                dt*d_half*m2p%rhox2d(j,i)*hydf(j,i,kz)
          end do
        else
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            p2m%chiten(j,i,k,n) = p2m%chiten(j,i,k,n) +  &
                          (tpred1(j,i,k)-m2p%chib(j,i,k,n))*rdt*m2p%psb(j,i)
          end do
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            p2m%remdrd(j,i,n) = p2m%remdrd(j,i,n) + &
                m2p%chib(j,i,kz,n)*m2p%drydepv(j,i,n) * &
                dt*d_half*m2p%rhox2d(j,i)*hydf(j,i,kz)*m2p%psb(j,i)
          end do
        end if
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
    !@acc call nvtxEndRange
  end subroutine holtbl

#include <pfesat.inc>
#include <pfqsat.inc>

  pure real(rkx) function comp_obklen(thvs,ustar,bfs) result(obk)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: thvs, ustar, bfs
    obk = -(thvs*ustar**3) / (gvk*bfs+sign(1.0e-10_rkx,bfs))
  end function comp_obklen

  subroutine force_water_conserve(tendv,tendc,tendi,start,sflux,ps)
    implicit none
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: tendv
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: tendc
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: tendi
    real(rkx), dimension(:,:,:,:), pointer, contiguous, intent(in) :: start
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: sflux, ps
    integer(ik4) :: i, j, k
    real(rkx), dimension(kz) :: qi, qf
    real(rkx) :: sqtoti, sqtotf

#ifdef STDPAR_FIXED
    do concurrent ( j = jci1:jci2, i = ici1:ici2 ) local(qi,qf)
#else
    !$acc parallel loop collapse(2) gang vector private(qi,qf)
    do i = ici1, ici2
    do j = jci1, jci2
#endif
      qi(1:kz) = start(j,i,1:kz,iqv)
      qf(1:kz) = qi(1:kz) + tendv(j,i,1:kz) * dt
      if ( idynamic /= 3 ) then
        do k = 1, kz
          qi(k) = qi(k) / ps(j,i)
          qf(k) = qf(k) / ps(j,i)
        end do
      end if
      sqtoti = sum(qi)
      sqtotf = sum(qf)
      sqtoti = sqtoti + sflux(j,i)
      if ( abs(sqtotf-sqtoti) > minqq ) then
        k = maxloc(qi,1)
        tendv(j,i,k) = tendv(j,i,k) + (sqtoti-sqtotf) * rdt
      end if
#ifndef STDPAR_FIXED
    end do
#endif
    end do
#ifdef STDPAR_FIXED
    do concurrent ( j = jci1:jci2, i = ici1:ici2 ) local(qi,qf)
#else
    !$acc parallel loop collapse(2) gang vector private(qi,qf)
    do i = ici1, ici2
    do j = jci1, jci2
#endif
      qi(1:kz) = start(j,i,1:kz,iqc)
      qf(1:kz) = qi(1:kz) + tendc(j,i,1:kz) * dt
      if ( idynamic /= 3 ) then
        do k = 1, kz
          qi(k) = qi(k) / ps(j,i)
          qf(k) = qf(k) / ps(j,i)
        end do
      end if
      sqtoti = sum(qi)
      sqtotf = sum(qf)
      if ( abs(sqtotf-sqtoti) > minqq ) then
        k = maxloc(qi,1)
        tendc(j,i,k) = tendc(j,i,k) + (sqtoti-sqtotf) * rdt
      end if
#ifndef STDPAR_FIXED
    end do
#endif
    end do
    if ( ipptls > 1 ) then
#ifdef STDPAR_FIXED
      do concurrent ( j = jci1:jci2, i = ici1:ici2 ) local(qi,qf)
#else
      !$acc parallel loop collapse(2) gang vector private(qi,qf)
      do i = ici1, ici2
      do j = jci1, jci2
#endif
        qi(1:kz) = start(j,i,1:kz,iqi)
        qf(1:kz) = qi(1:kz) + tendi(j,i,1:kz) * dt
        if ( idynamic /= 3 ) then
          do k = 1, kz
            qi(k) = qi(k) / ps(j,i)
            qf(k) = qf(k) / ps(j,i)
          end do
        end if
        sqtoti = sum(qi)
        sqtotf = sum(qf)
        if ( abs(sqtotf-sqtoti) > minqq ) then
          k = maxloc(qi,1)
          tendi(j,i,k) = tendi(j,i,k) + (sqtoti-sqtotf) * rdt
        end if
#ifndef STDPAR_FIXED
      end do
#endif
      end do
    end if
  end subroutine force_water_conserve

end module mod_pbl_holtbl
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

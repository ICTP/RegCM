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
!

! fsolmon@ictp.it: this is the interface for calling rrtm

module mod_rrtmg_driver
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_rad_common
  use mod_rad_scenarios
  use mod_rad_tracer
  use mod_rad_aerosol
  use mod_rad_atmosphere
  use rrtmg_sw_rad
  use mcica_subcol_gen_sw
!  use rrtmg_lw_rad_nomcica
  use parrrsw
  use rrsw_wvn
  use parrrtm
  use rrtmg_lw_rad
  use mod_rad_outrad
  use mod_mpmessage
  use mod_runparams
  use mod_regcm_types

  implicit none

  private
!
  public :: allocate_mod_rad_rrtmg , rrtmg_driver

  real(rk8) , pointer , dimension(:) :: frsa , sabtp , clrst , solin , &
         clrss , firtp , frla , clrlt , clrls , empty1 , abv , sol ,  &
         sols , soll , solsd , solld , slwd , tsfc , psfc , asdir ,   &
         asdif , aldir , aldif , czen , dlat , xptrop , totcf ,       &
         totcl , totci

  real(rk8) , pointer , dimension(:,:) :: qrs , qrl , clwp_int , pint,   &
    rh , cld_int , empty2 , play , tlay , h2ovmr , o3vmr , co2vmr ,      &
    ch4vmr , n2ovmr , o2vmr , cfc11vmr , cfc12vmr , cfc22vmr,  ccl4vmr , &
    reicmcl , relqmcl , swhr , swhrc , ciwp , clwp , rei , rel,cldf ,    &
    lwhr , lwhrc , duflx_dt , duflxc_dt , ql1 , qi1

  real(rk8) , pointer , dimension(:,:) :: plev , tlev , swuflx , swdflx , &
    swuflxc , swdflxc , lwuflx , lwdflx , lwuflxc , lwdflxc ,             &
    swddiruviflx , swddifuviflx , swddirpirflx , swddifpirflx , swdvisflx

  real(rk8) , pointer , dimension(:,:,:) :: cldfmcl , taucmcl , ssacmcl , &
         asmcmcl , fsfcmcl , ciwpmcl , clwpmcl
  real(rk8) , pointer , dimension(:,:,:) :: cldfmcl_lw , taucmcl_lw , &
         ciwpmcl_lw , clwpmcl_lw , aermmr

  real(rk8) , pointer , dimension(:) :: aeradfo , aeradfos , &
    asaeradfo , asaeradfos
  real(rk8) , pointer , dimension(:) :: aerlwfo , aerlwfos , &
    asaerlwfo , asaerlwfos
  real(rk8) , pointer , dimension(:,:) :: fice , wcl , wci , gcl , gci , &
         fcl , fci , tauxcl , tauxci , h2ommr , n2ommr , ch4mmr ,       &
         cfc11mmr , cfc12mmr , deltaz , dzr
  real(rk8) , pointer , dimension(:,:,:) :: outtaucl , outtauci

  integer(ik4) , pointer , dimension(:) :: ioro

  ! spectral dependant quantities

  real(rk8) , pointer , dimension(:,:,:) :: tauaer , ssaaer , asmaer , ecaer
  ! tauc = in-cloud optical depth
  ! ssac = in-cloud single scattering albedo (non-delta scaled)
  ! asmc = in-cloud asymmetry parameter (non-delta scaled)
  ! fsfc = in-cloud forward scattering fraction (non-delta scaled)
  real(rk8) , pointer , dimension(:,:,:) :: tauc , ssac , asmc , fsfc
  real(rk8) , pointer , dimension(:,:) :: emis_surf
  real(rk8) , pointer , dimension(:,:,:) :: tauc_lw
  real(rk8) , pointer , dimension(:,:,:) :: tauaer_lw
  integer(ik4) :: npr , kth , ktf , kclimf , kclimh

  integer(ik4) :: permuteseed = 1 , mypid

  contains

  subroutine allocate_mod_rad_rrtmg
    implicit none
    integer(ik4) :: k

#if defined ( __PGI ) || defined ( __OPENCC__ ) || defined ( __INTEL_COMPILER )
    integer , external :: getpid
#endif
#if defined ( IBM )
    integer , external :: getpid_
#endif
    npr = (jci2-jci1+1)*(ici2-ici1+1)

    ! Define here the total number of vertical levels, including standard
    ! atmosphere hat replace kz by kth, kzp1 by ktf

    do k = 1 , n_prehlev
      kclimh = k
      if ( ptop*d_10 > stdplevh(k) ) exit
    end do
    kth = kz + n_prehlev - kclimh - 1
    do k = 1 , n_preflev
      kclimf = k
      if ( ptop*d_10 > stdplevf(k) ) exit
    end do
    ktf = kzp1 + n_preflev - kclimf - 1

    call getmem1d(frsa,1,npr,'rrtmg:frsa')
    call getmem1d(sabtp,1,npr,'rrtmg:sabtp')
    call getmem1d(clrst,1,npr,'rrtmg:clrst')
    call getmem1d(clrss,1,npr,'rrtmg:clrss')
    call getmem1d(firtp,1,npr,'rrtmg:firtp')
    call getmem1d(frla,1,npr,'rrtmg:frla')
    call getmem1d(clrlt,1,npr,'rrtmg:clrlt')
    call getmem1d(clrls,1,npr,'rrtmg:clrls')
    call getmem1d(empty1,1,npr,'rrtmg:empty1')
    call getmem1d(solin,1,npr,'rrtmg:solin')
    call getmem1d(sols,1,npr,'rrtmg:sols')
    call getmem1d(soll,1,npr,'rrtmg:soll')
    call getmem1d(solsd,1,npr,'rrtmg:solsd')
    call getmem1d(solld,1,npr,'rrtmg:solld')
    call getmem1d(slwd,1,npr,'rrtmg:slwd')
    call getmem2d(qrs,1,npr,1,kth,'rrtmg:qrs')
    call getmem2d(qrl,1,npr,1,kth,'rrtmg:qrl')
    call getmem2d(clwp_int,1,npr,1,kth,'rrtmg:clwp_int')
    call getmem2d(cld_int,1,npr,1,kth,'rrtmg:cld_int')
    call getmem2d(empty2,1,npr,1,kth,'rrtmg:empty2')
    call getmem1d(aeradfo,1,npr,'rrtmg:aeradfo')
    call getmem1d(aeradfos,1,npr,'rrtmg:aeradfos')
    call getmem1d(asaeradfo,1,npr,'rrtmg:aseradfo')
    call getmem1d(asaeradfos,1,npr,'rrtmg:aseradfos')
    call getmem1d(aerlwfo,1,npr,'rrtmg:aerlwfo')
    call getmem1d(aerlwfos,1,npr,'rrtmg:aerlwfos')
    call getmem1d(asaerlwfo,1,npr,'rrtmg:asaerlwfo')
    call getmem1d(asaerlwfos,1,npr,'rrtmg:asaerlwfos')
    call getmem1d(abv,1,npr,'rrtmg:sol')
    call getmem1d(sol,1,npr,'rrtmg:sol')
    call getmem1d(tsfc,1,npr,'rrtmg:tsfc')
    call getmem1d(psfc,1,npr,'rrtmg:psfc')
    call getmem1d(asdir,1,npr,'rrtmg:asdir')
    call getmem1d(asdif,1,npr,'rrtmg:asdif')
    call getmem1d(aldir,1,npr,'rrtmg:aldir')
    call getmem1d(aldif,1,npr,'rrtmg:aldif')
    call getmem1d(czen,1,npr,'rrtmg:czen')
    call getmem1d(dlat,1,npr,'rrtmg:dlat')
    call getmem1d(xptrop,1,npr,'rrtmg:xptrop')
    call getmem1d(ioro,1,npr,'rrtmg:ioro')
    call getmem1d(totcf,1,npr,'rrtmg:totcf')
    call getmem1d(totcl,1,npr,'rrtmg:totlf')
    call getmem1d(totci,1,npr,'rrtmg:totif')

    call getmem3d(aermmr,1,npr,1,kz,1,ntr,'rrtmg:aermmr')

    call getmem2d(pint,1,npr,1,kzp1,'rrtmg:pint')
    call getmem2d(rh,1,npr,1,kzp1,'rrtmg:rh')
    call getmem2d(play,1,npr,1,kth,'rrtmg:play')
    call getmem2d(tlay,1,npr,1,kth,'rrtmg:tlay')
    if ( ipptls == 2 ) then
      call getmem2d(ql1,1,npr,1,kth,'rrtmg:ql1')
      call getmem2d(qi1,1,npr,1,kth,'rrtmg:qi1')
    end if
    call getmem2d(h2ovmr,1,npr,1,kth,'rrtmg:h2ovmr')
    call getmem2d(o3vmr,1,npr,1,kth,'rrtmg:o3vmr')
    call getmem2d(co2vmr,1,npr,1,kth,'rrtmg:co2vmr')
    call getmem2d(ch4vmr,1,npr,1,kth,'rrtmg:ch4vmr')
    call getmem2d(n2ovmr,1,npr,1,kth,'rrtmg:n2ovmr')
    call getmem2d(o2vmr,1,npr,1,kth,'rrtmg:o2vmr')
    call getmem2d(cfc11vmr,1,npr,1,kth,'rrtmg:cfc11vmr')
    call getmem2d(cfc12vmr,1,npr,1,kth,'rrtmg:cfc12vmr')
    call getmem2d(cfc22vmr,1,npr,1,kth,'rrtmg:cfc22vmr')
    call getmem2d(ccl4vmr,1,npr,1,kth,'rrtmg:ccl4vmr')
    call getmem2d(reicmcl,1,npr,1,kth,'rrtmg:reicmcl')
    call getmem2d(relqmcl,1,npr,1,kth,'rrtmg:relqmcl')
    call getmem2d(swhr,1,npr,1,kth,'rrtmg:swhr')
    call getmem2d(swhrc,1,npr,1,kth,'rrtmg:swhrc')
    call getmem2d(ciwp,1,npr,1,kth,'rrtmg:ciwp')
    call getmem2d(clwp,1,npr,1,kth,'rrtmg:clwp')
    call getmem2d(rei,1,npr,1,kth,'rrtmg:rei')
    call getmem2d(rel,1,npr,1,kth,'rrtmg:rel')
    call getmem2d(cldf,1,npr,1,kth,'rrtmg:cldf')
    call getmem2d(lwhr,1,npr,1,kth,'rrtmg:lwhr')
    call getmem2d(lwhrc,1,npr,1,kth,'rrtmg:lwhrc')
    call getmem2d(duflx_dt,1,npr,1,kth,'rrtmg:duflx_dt')
    call getmem2d(duflxc_dt,1,npr,1,kth,'rrtmg:duflxc_dt')
    call getmem2d(plev,1,npr,1,ktf,'rrtmg:plev')
    call getmem2d(tlev,1,npr,1,ktf,'rrtmg:tlev')
    call getmem2d(swuflx,1,npr,1,ktf,'rrtmg:swuflx')
    call getmem2d(swdflx,1,npr,1,ktf,'rrtmg:swdflx')
    call getmem2d(swuflxc,1,npr,1,ktf,'rrtmg:swuflxc')
    call getmem2d(swdflxc,1,npr,1,ktf,'rrtmg:swdflxc')
    call getmem2d(lwuflx,1,npr,1,ktf,'rrtmg:lwuflx')
    call getmem2d(lwdflx,1,npr,1,ktf,'rrtmg:lwdflx')
    call getmem2d(lwuflxc,1,npr,1,ktf,'rrtmg:lwuflxc')
    call getmem2d(lwdflxc,1,npr,1,ktf,'rrtmg:lwdflxc')
    call getmem2d(swddiruviflx,1,npr,1,ktf,'rrtmg:swddiruviflx')
    call getmem2d(swddifuviflx,1,npr,1,ktf,'rrtmg:swddifuviflx')
    call getmem2d(swddirpirflx,1,npr,1,ktf,'rrtmg:swddirpirflx')
    call getmem2d(swddifpirflx,1,npr,1,ktf,'rrtmg:swddifpirflx')
    call getmem2d(swdvisflx,1,npr,1,ktf,'rrtmg:swdvisflx')
    call getmem3d(cldfmcl,1,ngptsw,1,npr,1,kth,'rrtmg:cldfmcl')
    call getmem3d(taucmcl,1,ngptsw,1,npr,1,kth,'rrtmg:taucmcl')
    call getmem3d(ssacmcl,1,ngptsw,1,npr,1,kth,'rrtmg:ssacmcl')
    call getmem3d(asmcmcl,1,ngptsw,1,npr,1,kth,'rrtmg:asmcmcl')
    call getmem3d(fsfcmcl,1,ngptsw,1,npr,1,kth,'rrtmg:fsfcmcl')
    call getmem3d(ciwpmcl,1,ngptsw,1,npr,1,kth,'rrtmg:ciwpmcl')
    call getmem3d(clwpmcl,1,ngptsw,1,npr,1,kth,'rrtmg:clwpmcl')
    call getmem3d(cldfmcl_lw,1,ngptlw,1,npr,1,kth,'rrtmg:cldfmcl_lw')
    call getmem3d(taucmcl_lw,1,ngptlw,1,npr,1,kth,'rrtmg:taucmcl_lw')
    call getmem3d(ciwpmcl_lw,1,ngptlw,1,npr,1,kth,'rrtmg:ciwpmcl_lw')
    call getmem3d(clwpmcl_lw,1,ngptlw,1,npr,1,kth,'rrtmg:clwpmcl_lw')
    call getmem3d(tauaer,1,npr,1,kth,1,nbndsw,'rrtmg:tauaer')
    call getmem3d(ssaaer,1,npr,1,kth,1,nbndsw,'rrtmg:ssaaer')
    call getmem3d(asmaer,1,npr,1,kth,1,nbndsw,'rrtmg:asmaer')
    call getmem3d(ecaer,1,npr,1,kth,1,nbndsw,'rrtmg:ecaer')
    call getmem3d(tauc,1,nbndsw,1,npr,1,kth,'rrtmg:tauc')
    call getmem3d(ssac,1,nbndsw,1,npr,1,kth,'rrtmg:ssac')
    call getmem3d(asmc,1,nbndsw,1,npr,1,kth,'rrtmg:asmc')
    call getmem3d(fsfc,1,nbndsw,1,npr,1,kth,'rrtmg:fsfc')
    call getmem2d(emis_surf,1,npr,1,nbndlw,'rrtmg:emis_surf')
    call getmem3d(tauaer_lw,1,npr,1,kth,1,nbndlw,'rrtmg:tauaer_lw')
    call getmem3d(tauc_lw,1,nbndlw,1,npr,1,kth,'rrtmg:tauc_lw')
    call getmem2d(fice,1,npr,1,kth,'rrtmg:fice')
    call getmem2d(wcl,1,npr,1,kth,'rrtmg:wcl')
    call getmem2d(wci,1,npr,1,kth,'rrtmg:wci')
    call getmem2d(gcl,1,npr,1,kth,'rrtmg:gcl')
    call getmem2d(gci,1,npr,1,kth,'rrtmg:gci')
    call getmem2d(fcl,1,npr,1,kth,'rrtmg:fcl')
    call getmem2d(fci,1,npr,1,kth,'rrtmg:fci')
    call getmem2d(tauxcl,1,npr,1,kth,'rrtmg:tauxcl')
    call getmem2d(tauxci,1,npr,1,kth,'rrtmg:tauxci')
    call getmem3d(outtaucl,1,npr,1,kz,1,4,'rrtmg:outtaucl')
    call getmem3d(outtauci,1,npr,1,kz,1,4,'rrtmg:outtauci')
    call getmem2d(h2ommr,1,npr,1,kth,'rrtmg:h2ommr')
    call getmem2d(n2ommr,1,npr,1,kth,'rrtmg:n2ommr')
    call getmem2d(ch4mmr,1,npr,1,kth,'rrtmg:ch4mmr')
    call getmem2d(cfc11mmr,1,npr,1,kth,'rrtmg:cfc11mmr')
    call getmem2d(cfc12mmr,1,npr,1,kth,'rrtmg:cfc12mmr')
    call getmem2d(deltaz,1,npr,1,kth,'rrtmg:deltaz')
    call getmem2d(dzr,1,npr,1,kth,'rrtmg:dzr')

#if defined ( IBM )
    mypid = getpid_()
#else
    mypid = getpid()
#endif

  end subroutine allocate_mod_rad_rrtmg
!
  subroutine rrtmg_driver(iyear,lout,m2r,r2a)
    implicit none
    type(mod_2_rad) , intent(in) :: m2r
    type(rad_2_mod) , intent(inout) :: r2a
    integer(ik4) , intent(in) :: iyear
    logical, intent(in) :: lout
    integer(ik4) :: dyofyr , iplon , k , kj , n , i , j
    real(rk8) :: adjes

!-----------------------------------------------------------------------

    iplon = 1 ! not effectively used

    ! adjustment of earthsun distance:
    ! if dyofyr = 0, we pass directly the ecc. factor in adjes,
    ! otherwise it is calculated in RRTM as a function of day of the year.
    dyofyr = 0
    adjes = eccf

    ! from water path and cloud radius / tauc_LW is not requested
    tauc_lw = dlowval
    call prep_dat_rrtm(m2r,iyear,inflgsw)

    !
    ! Call to the shortwave radiation code as soon one element of czen is > 0.
    !

    if ( maxval(m2r%coszrs) > dlowval) then
      ! generates cloud properties:
      permuteseed = permuteseed+mypid+ngptlw
      if ( permuteseed < 0 ) permuteseed = 2147483641+permuteseed
      call mcica_subcol_sw(iplon,npr,kth,icld,permuteseed,irng,play,   &
                           cldf,ciwp,clwp,rei,rel,tauc,ssac,asmc,fsfc, &
                           cldfmcl,ciwpmcl,clwpmcl,reicmcl,relqmcl,    &
                           taucmcl,ssacmcl,asmcmcl,fsfcmcl)
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          asdir(n) = m2r%aldirs(j,i)
          asdif(n) = m2r%aldifs(j,i)
          aldir(n) = m2r%aldirl(j,i)
          aldif(n) = m2r%aldifl(j,i)
          czen(n)  = m2r%coszrs(j,i)
          n = n + 1
        end do
      end do
      call rrtmg_sw(npr,kth,icld,idirect,play,plev,tlay,tlev,tsfc,   &
                    h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr, &
                    asdir,asdif,aldir,aldif,czen,adjes,      &
                    dyofyr,solcon,inflgsw,iceflgsw,liqflgsw, &
                    cldfmcl,taucmcl,ssacmcl,asmcmcl,fsfcmcl, &
                    ciwpmcl,clwpmcl,reicmcl,relqmcl,tauaer,  &
                    ssaaer,asmaer,ecaer,swuflx,swdflx,swhr,  &
                    swuflxc,swdflxc,swhrc,swddiruviflx,      &
                    swddifuviflx,swddirpirflx,swddifpirflx,  &
                    swdvisflx,aeradfo,aeradfos,asaeradfo,asaeradfos)

    else
      swuflx(:,:) = d_zero
      swdflx(:,:) = d_zero
      swuflxc(:,:) = d_zero
      swdflxc(:,:) = d_zero
      swhr(:,:) = d_zero
      swhrc(:,:) = d_zero
      swdvisflx(:,:) = d_zero
      swddiruviflx(:,:) = d_zero
      swddifuviflx(:,:) = d_zero
      swddirpirflx(:,:) = d_zero
      swddifpirflx(:,:) = d_zero
      aeradfo(:) = d_zero
      aeradfos(:) = d_zero
      asaeradfo(:) = d_zero
      asaeradfos(:) = d_zero
    end if ! end shortwave call

    ! LW call :
    idrv =0 
    if ( 1==1) then
    permuteseed = permuteseed+mypid+ngptsw
    if ( permuteseed < 0 ) permuteseed = 2147483641+permuteseed
    call mcica_subcol_lw(iplon,npr,kth,icld,permuteseed,irng,play,  &
                         cldf,ciwp,clwp,rei,rel,tauc_lw,cldfmcl_lw, &
                         ciwpmcl_lw,clwpmcl_lw,reicmcl,relqmcl,     &
                         taucmcl_lw)
   

    call rrtmg_lw(npr,kth,icld,idrv,idirect,play,plev,tlay,tlev,tsfc,  &
                  h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr,      &
                  cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,emis_surf, &
                  inflglw,iceflglw,liqflglw,cldfmcl_lw,         &
                  taucmcl_lw,ciwpmcl_lw,clwpmcl_lw,reicmcl,     &
                  relqmcl,tauaer_lw,lwuflx,lwdflx,lwhr,lwuflxc, &
                  lwdflxc,lwhrc,duflx_dt,duflxc_dt,aerlwfo,     &
                  aerlwfos,asaerlwfo,asaerlwfos)
!
!    print*, 'rrtmg', tlay(50,:)
!    print*, 'rrtmg', lwdflx(50,:)
   
    else 
!   call rrtmg_lw_nomcica( npr,kth,icld,idrv,play,plev,tlay,tlev,tsfc,  &
!                  h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr,      &
!                  cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,emis_surf, &
!                  inflglw,iceflglw,liqflglw,cldf,         &
!                  tauc,ciwp,clwp,rei,     &
!                  rel,tauaer_lw,lwuflx,lwdflx,lwhr,lwuflxc, &
!                  lwdflxc,lwhrc,duflx_dt,duflxc_dt )
!    print*, 'rrtmg2', tlay(50,:)
!    print*, 'rrtmg2', lwdflx(50,:)
     end if 
    ! Output and interface
    !
    ! EES  next 3 added, they are calculated in radcsw : but not used further
    ! fsnirt   - Near-IR flux absorbed at toa
    ! fsnrtc   - Clear sky near-IR flux absorbed at toa
    ! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
    ! fsds     - Flux Shortwave Downwelling Surface

    solin(:) = swdflx(:,kzp1)
    frsa(:)  = swdflx(:,1) - swuflx(:,1)
    sabtp(:) = swdflx(:,kzp1) - swuflx(:,kzp1)
    clrst(:) = swdflxc(:,kzp1) - swuflxc(:,kzp1)
    clrss(:) = swdflxc(:,1) - swuflxc(:,1)

    firtp(:) = -d_one * (lwdflx(:,kzp1) - lwuflx(:,kzp1))
    frla(:)  = -d_one * (lwdflx(:,1) - lwuflx(:,1))
    clrlt(:) = -d_one * (lwdflxc(:,kzp1) - lwuflxc(:,kzp1))
    clrls(:) = -d_one * (lwdflxc(:,1) - lwuflxc(:,1))

    ! coupling with BATS
    !  r2m%sabveg set to frsa (as in standard version : potential inconsistency
    !  if soil fraction is large)
    ! solar is normally the visible band only total incident surface flux
    abv(:) = frsa(:)
    sol(:) = swdvisflx(:,1)
    ! surface SW incident
    sols(:) =  swddiruviflx(:,1)
    solsd(:) =  swddifuviflx(:,1)
    soll(:) =  swddirpirflx(:,1)
    solld(:) =  swddifpirflx(:,1)
    ! LW incident
    slwd(:) = lwdflx(:,1)

    ! 3d heating rate back on regcm grid and converted to K.S-1
    do k = 1 , kz
      kj = kzp1-k
      qrs(:,kj) = swhr(:,k) / secpd
      qrl(:,kj) = lwhr(:,k) / secpd
      cld_int(:,kj) = cldf(:,k) !ouptut : these are in cloud diagnostics
      clwp_int(:,kj) = clwp(:,k)
      dzr(:,kj) = deltaz(:,k)
    end do

    ! Finally call radout for coupling to BATS/CLM/ATM and outputing fields
    ! in RAD files : these are diagnostics that are passed to radout but
    ! not used furthermore
    empty1 = dmissval
    empty2 = dmissval

    ! Calculate cloud parameters
    totcf(:) = d_one
    do k = 1 , kz
      totcf(:) = totcf*(d_one-cld_int(:,k))
    end do
    totcf(:) = d_one - totcf(:)
    totcf(:) = d_half * ( totcf(:) + maxval(cld_int(:,:),2) )
    if ( lsrfhack ) then
      where ( totcf > cftotmax )
        totcf = cftotmax
      end where
      where ( totcf < d_zero )
        totcf = d_zero
      end where
    end if

    totci(:) = sum(clwp_int*fice,2)*d_r1000
    totcl(:) = sum(clwp_int,2)*d_r1000

    call radout(lout,solin,sabtp,frsa,clrst,clrss,qrs,            &
                firtp,frla,clrlt,clrls,qrl,slwd,sols,soll,solsd,  &
                solld,totcf,totcl,totci,cld_int,clwp_int,abv,     &
                sol,aeradfo,aeradfos,aerlwfo,aerlwfos,tauxar3d,   &
                tauasc3d,gtota3d,dzr,outtaucl,outtauci,r2a, &
                 asaeradfo,asaeradfos,asaerlwfo,asaerlwfos)
  end subroutine rrtmg_driver

  subroutine prep_dat_rrtm(m2r,iyear,inflagsw)
    implicit none
    type(mod_2_rad) , intent(in) :: m2r
    integer(ik4) , intent(in) :: iyear , inflagsw
    real(rk8) :: ccvtem , clwtem , w1 , w2
    integer(ik4) :: i , j , k , kj , ncldm1 , ns , n , itr
    real(rk8) , parameter :: lowcld = 1.0D-30
    real(rk8) , parameter :: verynearone = 0.999999D0
    real(rk8) :: tmp1l , tmp2l , tmp3l , tmp1i , tmp2i , tmp3i
!
!   Set index for cloud particle properties based on the wavelength,
!   according to A. Slingo (1989) equations 1-3:
!   Use index 1 (0.25 to 0.69 micrometers) for visible
!   Use index 2 (0.69 - 1.19 micrometers) for near-infrared
!   Use index 3 (1.19 to 2.38 micrometers) for near-infrared
!   Use index 4 (2.38 to 4.00 micrometers) for near-infrared
!
!   Note that the minimum wavelength is encoded (with 0.001, 0.002,
!   0.003) in order to specify the index appropriate for the
!   near-infrared cloud absorption properties
!
!   the 14 RRTM SW bands intervall are given in wavenumber (wavenum1/2):
!   equivalence in micrometer are
!   3.07-3.84 / 2.50-3.07 / 2.15-2.5 / 1.94-2.15 / 1.62-1.94 /
!   1.29-1.62 / 1.24-1.29 / 0.778-1.24 / 0.625-0.778 / 0.441-0.625 /
!   0.344-0.441 / 0.263-0.344 / 0.200-0.263 / 3.846-12.195
!   A. Slingo's data for cloud particle radiative properties
!   (from 'A GCM Parameterization for the Shortwave Properties of Water
!   Clouds' JAS vol. 46 may 1989 pp 1419-1427)
!
!   abarl    - A coefficient for extinction optical depth
!   bbarl    - B coefficient for extinction optical depth
!   cbarl    - C coefficient for single particle scat albedo
!   dbarl    - D coefficient for single particle scat albedo
!   ebarl    - E coefficient for asymmetry parameter
!   fbarl    - F coefficient for asymmetry parameter
!
!   Caution... A. Slingo recommends no less than 4.0 micro-meters nor
!   greater than 20 micro-meters
!
!   ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
!
!   abari    - a coefficient for extinction optical depth
!   bbari    - b coefficient for extinction o
!
    integer, dimension (nbndsw) :: indsl
    real(rk8) , dimension(4) ::  abari , abarl , bbari , bbarl , cbari , &
                                cbarl , dbari , dbarl , ebari , ebarl , &
                                fbari , fbarl
    real(rk8) :: abarii , abarli , bbarii , bbarli , cbarii , cbarli , &
                dbarii , dbarli , ebarii , ebarli , fbarii , fbarli

    real(rk8) , parameter :: c287 = 0.287D+00
!
    data abarl / 2.817D-02 ,  2.682D-02 , 2.264D-02 , 1.281D-02/
    data bbarl / 1.305D+00 ,  1.346D+00 , 1.454D+00 , 1.641D+00/
    data cbarl /-5.620D-08 , -6.940D-06 , 4.640D-04 , 0.201D+00/
    data dbarl / 1.630D-07 ,  2.350D-05 , 1.240D-03 , 7.560D-03/
    data ebarl / 0.829D+00 ,  0.794D+00 , 0.754D+00 , 0.826D+00/
    data fbarl / 2.482D-03 ,  4.226D-03 , 6.560D-03 , 4.353D-03/
!
    data abari / 3.4480D-03 , 3.4480D-03 , 3.4480D-03 , 3.44800D-03/
    data bbari / 2.4310D+00 , 2.4310D+00 , 2.4310D+00 , 2.43100D+00/
    data cbari / 1.0000D-05 , 1.1000D-04 , 1.8610D-02 , 0.46658D+00/
    data dbari / 0.0000D+00 , 1.4050D-05 , 8.3280D-04 , 2.05000D-05/
    data ebari / 0.7661D+00 , 0.7730D+00 , 0.7940D+00 , 0.95950D+00/
    data fbari / 5.8510D-04 , 5.6650D-04 , 7.2670D-04 , 1.07600D-04/
    !
    ! define index pointing on appropriate parameter in slingo's table
    ! for eachRRTM SW band
    !
    data indsl /4,4,3,3,3,3,3,2,2,1,1,1,1,4 /
    ! CONVENTION : RRTMG driver takes layering form botom to TOA.
    ! regcm consider Top to bottom

    ! surface pressure and scaled pressure, from which level are computed
    ! RRTM SW takes pressure in mb,hpa
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        dlat(n) = dabs(m2r%xlat(j,i))
        xptrop(n) = m2r%ptrop(j,i) * d_r100
        n = n + 1
      end do
    end do

    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        psfc(n) = m2r%psatms(j,i) * d_r100
        n = n + 1
      end do
    end do

    do k = 1 , kz
      n = 1
      kj = kzp1-k
      do i = ici1 , ici2
        do j = jci1 , jci2
          play(n,kj) = m2r%phatms(j,i,k)*d_r100
          n = n + 1
        end do
      end do
    end do
    do k = kzp1 , kth
      play(:,k) = stdplevh(kclimh+k-kzp1)
    end do
    !
    ! convert pressures from Pa to mb and define interface pressures:
    !
    do k = 1 , kzp1
      n = 1
      kj = kzp1-k+1
      do i = ici1 , ici2
        do j = jci1 , jci2
          plev(n,kj) = m2r%pfatms(j,i,k)*d_r100
          n = n + 1
        end do
      end do
    end do
    do k = kzp1+1 , ktf
      plev(:,k) = stdplevf(kclimf+k-kzp1)
    end do
    !
    ! ground temperature
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        tsfc(n) = m2r%tg(j,i)
        n = n + 1
      end do
    end do
    !
    ! air temperatures
    !
    do k = 1 , kz
      kj = kzp1-k
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          tlay(n,kj) = m2r%tatms(j,i,k)
          n = n + 1
        end do
      end do
    end do
    do k = kzp1 , kth
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          tlay(n,k) = stdatm_val(calday,m2r%xlat(j,i),play(n,k),istdatm_tempk)
          n = n + 1
        end do
      end do
    end do
    !
    ! air temperature at the interface (same method as in vertical
    ! advection routine)
    !
    do n = 1 , npr
      tlev(n,1) = tsfc(n)
    end do
    do k = 2 , kz
      kj = kzp1-k+1
      do n = 1 , npr
        w1 =  (hsigma(kj) - sigma(kj)) / (hsigma(kj) - hsigma(kj-1))
        w2 =  (sigma(kj) - hsigma(kj-1) ) / (hsigma(kj) - hsigma(kj-1))
        if (k < kz-1) then
          tlev(n,k) =  w1*tlay(n,k-1) * (plev(n,k)/play(n,k-1))**c287 + &
                       w2*tlay(n,k)   * (plev(n,k)/play(n,k))**c287
        else
          ! linear interpolation for last upper tropospheric levels
          tlev(n,k) =  w1*tlay(n,k-1) + w2*tlay(n,k)
        end if
      end do
    end do
    do k = kzp1 , ktf
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          tlev(n,k) = stdatm_val(calday,m2r%xlat(j,i),plev(n,k),istdatm_tempk)
          n = n + 1
        end do
      end do
    end do
    if ( ipptls == 2 ) then
      do k = 1 , kz
        kj = kzp1 - k
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            ql1(n,kj) = m2r%qxatms(j,i,k,iqc)
            qi1(n,kj) = m2r%qxatms(j,i,k,iqi)
            n = n + 1
          end do
        end do
      end do
      do k = kzp1 , kth
        ql1(:,k) = minqx
        qi1(:,k) = minqx
      end do
    end if
    !
    ! h2o volume mixing ratio
    !
    do k = 1 , kz
      kj = kzp1 - k
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          h2ommr(n,kj) = dmax1(1.0D-7,m2r%qxatms(j,i,k,iqv))
          h2ovmr(n,kj) = h2ommr(n,kj) * ep2 
          n = n + 1
        end do
      end do
    end do
    do k = kzp1 , kth
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          h2ommr(n,k) = stdatm_val(calday,m2r%xlat(j,i), &
                                   play(n,k),istdatm_qdens) / &
                        stdatm_val(calday,m2r%xlat(j,i), &
                                   play(n,k),istdatm_airdn)
          h2ovmr(n,k) = h2ommr(n,k) * ep2
          n = n + 1
        end do
      end do
    end do
    !
    ! o3 volume mixing ratio (already on the right grid)
    !
    do k = 1 , kz
      n = 1
      kj = kzp1 - k
      do i = ici1 , ici2
        do j = jci1 , jci2
          o3vmr(n,k) = d_half*(o3prof(j,i,kj+1)+o3prof(j,i,kj)) * amo/amd
          n = n + 1
        end do
      end do
    end do
    do k = kzp1 , kth
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          o3vmr(n,k) = stdatm_val(calday,m2r%xlat(j,i), &
                                  play(n,k),istdatm_ozone) / &
                       stdatm_val(calday,m2r%xlat(j,i), &
                                  play(n,k),istdatm_airdn) * amo/amd
          n = n + 1
        end do
      end do
    end do
!
    ! other gas (n2o,ch4)
    !
    if ( iyear >= 1750 .and. iyear <= 2100 ) then
!     co2vmr = cgas(2,iyear)*1.0D-6
!     co2mmr = co2vmr*44.0D0/28.9644D0
      ch40 = cgas(3,iyear)*1.0D-9*0.55241D0
      n2o0 = cgas(4,iyear)*1.0D-9*1.51913D0
      cfc110 = cgas(5,iyear)*1.0D-12*4.69548D0
      cfc120 = cgas(6,iyear)*1.0D-12*4.14307D0
    else
      write(stderr,*) 'Loading gas scenario for simulation year: ', iyear
      call fatal(__FILE__,__LINE__,                                   &
            'CONCENTRATION VALUES OUTSIDE OF DATE RANGE (1750-2100)')
    end if

    call trcmix(1,npr,dlat,xptrop,play,n2ommr,ch4mmr,cfc11mmr,cfc12mmr)

    do k = 1 , kth
      do n = 1 , npr
        n2ovmr(n,k)   = n2ommr(n,k) * (44.D0/amd)
        ch4vmr(n,k)   = ch4mmr(n,k) * (16.D0/amd)
        co2vmr(n,k)   = cgas(2,iyear) * 1.0D-6
        o2vmr(n,k)    = 0.23143D0 * (32.D0/amd)
        cfc11vmr(n,k) = cfc11mmr(n,k) / (163.1278D0/amd)
        cfc12vmr(n,k) = cfc12mmr(n,k) / (175.1385D0/amd)
        !
        ! No data FOR NOW : IMPROVE !!
        !
        cfc22vmr(n,k) = d_zero
        ccl4vmr(n,k)  = d_zero
      end do
    end do
    !
    ! aerosols
    ! no stratospheric background for now
    !care : Tracers mixing ratios are on regcm grid
    !
    if ( ichem == 1 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              aermmr(n,k,itr) = m2r%chiatms(j,i,k,itr)
              n = n + 1
            end do
          end do
        end do
      end do
    end if
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          rh(n,k) = dmax1(dmin1(m2r%rhatms(j,i,k),0.99D0),d_zero)
          n = n + 1
        end do
      end do
    end do
    do k = 1 , kzp1
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          pint(n,k) = m2r%pfatms(j,i,k)
          n = n + 1
        end do
      end do
    end do
    call aeroppt(rh,aermmr,pint,1,npr)
    ! adapt reverse the vertical grid for RRTM
    do k = 1 , kz
      kj = kzp1 - k
      do n = 1 , npr
       tauaer(n,kj,:) = tauxar3d(n,k,:)
       ssaaer(n,kj,:) = tauasc3d(n,k,:)
       asmaer(n,kj,:) = gtota3d(n,k,:)
       ecaer  = 0.78D0! not used
       tauaer_lw(n,kj,:) = tauxar3d_lw(n,k,:)
      end do
    end do
    !
    ! deltaz
    !
    do k = 1 , kz
      n = 1
      kj = kzp1 - k
      do i = ici1 , ici2
        do j = jci1 , jci2
          deltaz(n,k) = m2r%deltaz(j,i,kj)
        end do
      end do
    end do
    do k = kzp1 , kth
      do n = 1 , npr
        deltaz(n,k) = rgas*tlay(n,k)*(plev(n,k) - &
                      plev(n,k+1))/(egrav*play(n,k))
      end do
    end do
    !
    ! cloud fraction and cloud liquid waterpath calculation:
    ! as in STANDARD SCHEME for now (getdat) : We need to improve this
    ! according to new cloud microphysics!
    ! fractional cloud cover (dependent on relative humidity)
    !
    ! qc   = gary's mods for clouds/radiation tie-in to exmois
    !
    do k = 1 , kz
      kj = kzp1 - k
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          ccvtem = d_zero   !cqc mod
          cldf(n,kj) = dmax1(m2r%cldfrc(j,i,k)*0.9999999D0,ccvtem)
          cldf(n,kj) = dmin1(cldf(n,kj),0.9999999D0)
          !
          ! convert liquid water content into liquid water path, i.e.
          ! multiply b deltaz
          !
          clwtem = m2r%cldlwc(j,i,kj) ! put on the right grid !
          !
          ! deltaz,clwp are on the right grid since plev and tlay are
          ! care pressure is on botom/toa grid
          !
          clwp(n,k) = clwtem * deltaz(n,k)
          if ( dabs(cldf(n,k)) < lowcld ) then
            cldf(n,k) = d_zero
            clwp(n,k) = d_zero
          end if
          n = n + 1
        end do
      end do
    end do
    !
    ! fabtest:  set cloud fractional cover at top model level = 0
    ! as in std scheme
    cldf(:,kz-1:kth) = d_zero
    clwp(:,kz-1:kth) = d_zero
    !
    ! set cloud fractional cover at bottom (ncld) model levels = 0
    !
    ncldm1 = ncld - 1
    do k = 1 , ncldm1
      do n = 1 , npr
        cldf(n,k) = d_zero
        clwp(n,k) = d_zero
      end do
    end do
    !
    ! maximum cloud fraction
    !----------------------------------------------------------------------
    do k = 1 , kz
      do n = 1 , npr
        if ( cldf(n,k) > 0.999D0 ) cldf(n,k) = 0.999D0
      end do
    end do
    !
    !
    ! CLOUD Properties:
    !
    ! cloud effective radius
    !   NB: orography types are specified in the following
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        ioro(n) = m2r%ldmsk(j,i)
        n = n + 1
      end do
    end do
    !
    call cldefr_rrtm(tlay,rel,rei,fice,play)
    !
    !  partition of total water path betwwen liquide and ice.
    !
    do k = 1 , kth
      do n = 1 , npr
        ciwp(n,k) =  clwp(n,k) *  fice(n,k)
        clwp(n,k) =  clwp(n,k) * (d_one - fice(n,k))
        ! now clwp is liquide only !
      end do
    end do
    !
    ! Cloud optical properties(tau,ssa,g,f) :
    ! 2 options :
    ! inflgsw == 0 : treated as the standard radiation scheme and
    !                passed to McICA/ RRTM
    ! inflgsw == 2 : direcly calculated within RRTMsw
    !
    ! initialise and  begin spectral loop
    !
    tauc  =  d_zero
    ssac  =  verynearone
    asmc  =  0.850D0
    fsfc  =  0.725D0

#if defined(CLM45) || defined(CLM)
    ! TS is the radiant temeprature
    emis_surf(:,:) = d_one
#else
    do k = 1 , nbndlw
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          emis_surf(n,k) = m2r%emiss(j,i)
          n = n + 1
        end do
      end do
    end do
#endif

    outtaucl(:,:,:) = d_zero
    outtauci(:,:,:) = d_zero

    if ( inflagsw == 0 ) then
      do ns = 1 , nbndsw
        !
        !     Set cloud extinction optical depth, single scatter albedo,
        !     asymmetry parameter, and forward scattered fraction:
        !
        abarli = abarl(indsl(ns))
        bbarli = bbarl(indsl(ns))
        cbarli = cbarl(indsl(ns))
        dbarli = dbarl(indsl(ns))
        ebarli = ebarl(indsl(ns))
        fbarli = fbarl(indsl(ns))
        abarii = abari(indsl(ns))
        bbarii = bbari(indsl(ns))
        cbarii = cbari(indsl(ns))
        dbarii = dbari(indsl(ns))
        ebarii = ebari(indsl(ns))
        fbarii = fbari(indsl(ns))
!
        do k = 1 , kz
          do n = 1 , npr
            if ( clwp(n,k) < dlowval .and. ciwp(n,k) < dlowval) cycle
            ! liquid
            tmp1l = abarli + bbarli/rel(n,k)
            tmp2l = d_one - cbarli - dbarli*rel(n,k)
            tmp3l = fbarli*rel(n,k)
            ! ice
            tmp1i = abarii + bbarii/rei(n,k)
            tmp2i = d_one - cbarii - dbarii*rei(n,k)
            tmp3i = fbarii*rei(n,k)
            !
            ! cloud optical properties extinction optical depth
            !              IN_CLOUD quantities !
            !
            tauxcl(n,k) = clwp(n,k)*tmp1l
            tauxci(n,k) = ciwp(n,k)*tmp1i
            outtaucl(n,k,indsl(ns)) = outtaucl(n,k,indsl(ns)) + tauxcl(n,k)
            outtauci(n,k,indsl(ns)) = outtauci(n,k,indsl(ns)) + tauxci(n,k)
            wcl(n,k) = dmin1(tmp2l,verynearone)
            gcl(n,k) = ebarli + tmp3l
            fcl(n,k) = gcl(n,k)*gcl(n,k)
            wci(n,k) = dmin1(tmp2i,verynearone)
            gci(n,k) = ebarii + tmp3i
            fci(n,k) = gci(n,k)*gci(n,k)
            tauc(ns,:,k) = tauxcl(n,k) + tauxci(n,k)
            ssac(ns,:,k) = (tauxcl(n,k) * wcl(n,k) + &
                            tauxci(n,k) * wci(n,k) ) / tauc (ns,:,k)
            asmc(ns,:,k) = (tauxcl(n,k) * gcl(n,k) + &
                            tauxci(n,k) * gci(n,k) ) / tauc (ns,:,k)
            fsfc(ns,:,k) = (tauxcl(n,k) * fcl(n,k) + &
                            tauxci(n,k) * fci(n,k) ) / tauc (ns,:,k)
          end do
        end do
      end do ! spectral loop
    end if  ! inflagsw
  end subroutine prep_dat_rrtm
!
! for now we use for RRTM the same param as in standard rad
!
  subroutine cldefr_rrtm(t,rel,rei,fice,pmid)
    implicit none
    real(rk8) , pointer , dimension(:,:) :: fice , pmid , rei , rel , t
    intent (in) pmid , t
    intent (inout) fice , rei , rel
    integer(ik4) :: k , n
    real(rk8) :: pnrml , rliq , weight
    ! reimax - maximum ice effective radius
    real(rk8) , parameter :: reimax = 30.0D0
    ! rirnge - range of ice radii (reimax - 10 microns)
    real(rk8) , parameter :: rirnge = 20.0D0
    ! pirnge - nrmlzd pres range for ice particle changes
    real(rk8) , parameter :: pirnge = 0.4D0
    ! picemn - normalized pressure below which rei=reimax
    real(rk8) , parameter :: picemn = 0.4D0
    ! Temperatures in K (263.16 , 243.16)
    real(rk8) , parameter :: minus10 = wattp-d_10
    real(rk8) , parameter :: minus30 = wattp-(d_three*d_10)

    do k = 1 , kth
      do n = 1 , npr
        !
        ! Define liquid drop size
        !
        if ( ioro(n) == 0 ) then
          !
          ! Effective liquid radius over ocean and sea ice
          !
          rliq = d_10
        else
          !
          ! Effective liquid radius over land
          !
          rliq = d_five+d_five* &
                  dmin1(d_one,dmax1(d_zero,(minus10-t(n,k))*0.05D0))
        end if
        rel(n,k) = rliq
        !
        ! Determine rei as function of normalized pressure
        !
        pnrml = pmid(n,k)/psfc(n)
        weight = dmax1(dmin1((pnrml-picemn)/pirnge,d_one),d_zero)
        rei(n,k) = reimax - rirnge*weight
        !
        ! Define fractional amount of cloud that is ice
        !
        if ( ipptls == 2 ) then
          if ( ql1(n,k)+qi1(n,k) > minqx ) then
            fice(n,k) = qi1(n,k) / (ql1(n,k)+qi1(n,k))
          else
            fice(n,k) = d_zero
          end if
        else
          if ( t(n,k) > minus10 ) then
            ! if warmer than -10 degrees C then water phase
            fice(n,k) = d_zero
          else if ( t(n,k) <= minus10 .and. t(n,k) >= minus30 ) then
            ! if colder than -10 degrees C but warmer than -30 C mixed phase
            fice(n,k) = (minus10-t(n,k))/20.0D0
          else
            ! if colder than -30 degrees C then ice phase
            fice(n,k) = d_one
          end if
          !
        end if
        ! Turn off ice radiative properties by setting fice = 0.0
        ! fice(n,k) = d_zero
        !
      end do
    end do
  end subroutine cldefr_rrtm

end module mod_rrtmg_driver
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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
  use mod_constants
  use mod_dynparam
  use mod_ipcc_scenario
  use mod_memutil
  use mod_rad_common
  use mod_rad_tracer
  use mod_rad_aerosol
  use mod_stdatm
  use mod_isaatm
  use rrtmg_sw_rad
  use mcica_subcol_gen_sw
  use parrrsw
  use rrsw_wvn
  use parrrtm
  use rrtmg_lw_rad
  use rrtmg_lw_rad_nomcica
  use rrtmg_sw_rad_nomcica
  use mod_rad_outrad
  use mod_mpmessage
  use mod_runparams
  use mod_mppparam , only : italk
  use mod_regcm_types

  implicit none

  private
!
  public :: allocate_mod_rad_rrtmg , rrtmg_driver

  real(rkx) , pointer , dimension(:) :: frsa , sabtp , clrst , solin ,  &
         solout , clrss , firtp , lwout , lwin , frla , clrlt , clrls , &
         abv , sol , sols , soll , solsd , solld , slwd , tsfc , psfc , &
         asdir , asdif , aldir , aldif , czen , dlat , xptrop , totcf , &
         totwv, totcl , totci

  real(rkx) , pointer , dimension(:,:) :: qrs , qrl , clwp_int , pint, &
    rh , cld_int , tlay , h2ovmr , o3vmr , co2vmrk , play , ch4vmr ,   &
    n2ovmr , o2vmr , cfc11vmr , cfc12vmr , cfc22vmr,  ccl4vmr , o3 ,   &
    reicmcl , relqmcl , swhr , swhrc , ciwp , clwp , rei , rel , cldf ,&
    lwhr , lwhrc , duflx_dt , duflxc_dt , ql1 , qi1

  real(rkx) , pointer , dimension(:,:) :: plev , tlev , swuflx , swdflx , &
    swuflxc , swdflxc , lwuflx , lwdflx , lwuflxc , lwdflxc ,             &
    swddiruviflx , swddifuviflx , swddirpirflx , swddifpirflx , swdvisflx

  real(rkx) , pointer , dimension(:,:,:) :: cldfmcl , taucmcl , ssacmcl , &
         asmcmcl , fsfcmcl , ciwpmcl , clwpmcl
  real(rkx) , pointer , dimension(:,:,:) :: cldfmcl_lw , taucmcl_lw , &
         ciwpmcl_lw , clwpmcl_lw

  real(rkx) , pointer , dimension(:) :: aeradfo , aeradfos , &
    asaeradfo , asaeradfos
  real(rkx) , pointer , dimension(:) :: aerlwfo , aerlwfos , &
    asaerlwfo , asaerlwfos
  real(rkx) , pointer , dimension(:,:) :: fice , wcl , wci , gcl , gci , &
         fcl , fci , tauxcl , tauxci , h2ommr , n2ommr , ch4mmr ,       &
         cfc11mmr , cfc12mmr , deltaz , dzr
  real(rkx) , pointer , dimension(:) :: topz
  real(rkx) , pointer , dimension(:,:,:) :: outtaucl , outtauci

  integer(ik4) , pointer , dimension(:) :: ioro

  ! spectral dependant quantities

  real(rkx) , pointer , dimension(:,:,:) :: tauaer , ssaaer , asmaer , ecaer
  ! tauc = in-cloud optical depth
  ! ssac = in-cloud single scattering albedo (non-delta scaled)
  ! asmc = in-cloud asymmetry parameter (non-delta scaled)
  ! fsfc = in-cloud forward scattering fraction (non-delta scaled)
  real(rkx) , pointer , dimension(:,:,:) :: tauc , ssac , asmc , fsfc
  real(rkx) , pointer , dimension(:,:) :: emis_surf
  real(rkx) , pointer , dimension(:,:,:) :: tauc_lw
  real(rkx) , pointer , dimension(:,:,:) :: tauaer_lw
  integer(ik4) :: npr

  integer(ik4) :: permuteseed = 1 , mypid

  logical , parameter :: luse_max_rnovl = .true.

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

    call getmem1d(frsa,1,npr,'rrtmg:frsa')
    call getmem1d(sabtp,1,npr,'rrtmg:sabtp')
    call getmem1d(clrst,1,npr,'rrtmg:clrst')
    call getmem1d(clrss,1,npr,'rrtmg:clrss')
    call getmem1d(firtp,1,npr,'rrtmg:firtp')
    call getmem1d(lwout,1,npr,'rrtmg:lwout')
    call getmem1d(lwin,1,npr,'rrtmg:lwin')
    call getmem1d(frla,1,npr,'rrtmg:frla')
    call getmem1d(clrlt,1,npr,'rrtmg:clrlt')
    call getmem1d(clrls,1,npr,'rrtmg:clrls')
    call getmem1d(solin,1,npr,'rrtmg:solin')
    call getmem1d(solout,1,npr,'rrtmg:solout')
    call getmem1d(sols,1,npr,'rrtmg:sols')
    call getmem1d(soll,1,npr,'rrtmg:soll')
    call getmem1d(solsd,1,npr,'rrtmg:solsd')
    call getmem1d(solld,1,npr,'rrtmg:solld')
    call getmem1d(slwd,1,npr,'rrtmg:slwd')
    if ( idiag == 1 ) then
      call getmem2d(qrs,1,npr,1,kth,'rrtmg:qrs')
      call getmem2d(qrl,1,npr,1,kth,'rrtmg:qrl')
      call getmem2d(o3,1,npr,1,kth,'rrtmg:o3')
    end if
    call getmem2d(clwp_int,1,npr,1,kz,'rrtmg:clwp_int')
    call getmem2d(cld_int,1,npr,1,kzp1,'rrtmg:cld_int')
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
    call getmem1d(totwv,1,npr,'rrtmg:totwv')
    call getmem1d(totcl,1,npr,'rrtmg:totlf')
    call getmem1d(totci,1,npr,'rrtmg:totif')

    if ( ichem == 1 .or. iclimaaer > 0 ) then
      call getmem2d(pint,1,npr,1,kzp1,'rrtmg:pint')
      call getmem2d(rh,1,npr,1,kzp1,'rrtmg:rh')
    end if
    call getmem2d(play,1,npr,1,kth,'rrtmg:play')
    call getmem2d(tlay,1,npr,1,kth,'rrtmg:tlay')
    if ( ipptls > 1 ) then
      call getmem2d(ql1,1,npr,1,kth,'rrtmg:ql1')
      call getmem2d(qi1,1,npr,1,kth,'rrtmg:qi1')
    end if
    call getmem2d(h2ovmr,1,npr,1,kth,'rrtmg:h2ovmr')
    call getmem2d(o3vmr,1,npr,1,kth,'rrtmg:o3vmr')
    call getmem2d(co2vmrk,1,npr,1,kth,'rrtmg:co2vmrk')
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
    call getmem1d(topz,1,npr,'rrtmg:topz')
    call getmem2d(dzr,1,npr,1,kth,'rrtmg:dzr')

    call allocate_tracers(1,npr)

#if defined ( IBM )
    mypid = getpid_()
#else
    mypid = getpid()
#endif
  end subroutine allocate_mod_rad_rrtmg

  subroutine rrtmg_driver(iyear,imonth,lout,m2r,r2m)
    implicit none
    type(mod_2_rad) , intent(in) :: m2r
    type(rad_2_mod) , intent(inout) :: r2m
    integer(ik4) , intent(in) :: iyear , imonth
    logical , intent(in) :: lout
    integer(ik4) :: k , kj , n , i , j , kmincld , kmaxcld ,ldirect
    logical :: lradfor
    real(rkx) :: adjes

    ! from water path and cloud radius / tauc_LW is not requested
    tauc_lw(:,:,:) = dlowval
    call prep_dat_rrtm(m2r,iyear,imonth)
    adjes = real(eccf,rkx)

    lradfor = ( rcmtimer%start( ) .or. syncro_radfor%will_act( ) )

    !
    ! Call to the shortwave radiation code as soon one element of czen is > 0.
    !
    swuflx(:,:) = d_zero
    swdflx(:,:) = d_zero
    swuflxc(:,:) = d_zero
    swdflxc(:,:) = d_zero
    swhr(:,:) = d_zero        !
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

    ! hanlde aerosol direct effect in function of ichem or iclimaaer

    if ( ichem == 1 .and. iaerosol ==1 ) then
      ldirect = idirect
    elseif ( iclimaaer > 0 ) then
      ldirect = 2
    end if

    if ( maxval(m2r%coszrs) > 1.0e-3_rkx ) then
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          asdir(n) = m2r%aldirs(j,i)
          asdif(n) = m2r%aldifs(j,i)
          aldir(n) = m2r%aldirl(j,i)
          aldif(n) = m2r%aldifl(j,i)
          czen(n)  = m2r%coszrs(j,i)
          if ( czen(n) < 1.0e-3_rkx ) czen(n) = 0.0_rkx
          n = n + 1
        end do
      end do
      if ( imcica == 1 ) then
        ! generates cloud properties:
        permuteseed = permuteseed+mypid+ngptlw
        if ( permuteseed < 0 ) permuteseed = 2147483641+permuteseed
        call mcica_subcol_sw(npr,kth,icld,permuteseed,irng,play,         &
                             cldf,ciwp,clwp,rei,rel,tauc,ssac,asmc,fsfc, &
                             cldfmcl,ciwpmcl,clwpmcl,reicmcl,relqmcl,    &
                             taucmcl,ssacmcl,asmcmcl,fsfcmcl)
        call rrtmg_sw(npr,kth,icld,lradfor,ldirect,play,plev,tlay,tlev, &
                      tsfc,h2ovmr,o3vmr,co2vmrk,ch4vmr,n2ovmr,  &
                      o2vmr,asdir,asdif,aldir,aldif,czen,adjes, &
                      0,solcon,inflgsw,iceflgsw,liqflgsw,       &
                      cldfmcl,taucmcl,ssacmcl,asmcmcl,fsfcmcl,  &
                      ciwpmcl,clwpmcl,reicmcl,relqmcl,tauaer,   &
                      ssaaer,asmaer,ecaer,swuflx,swdflx,swhr,   &
                      swuflxc,swdflxc,swhrc,swddiruviflx,       &
                      swddifuviflx,swddirpirflx,swddifpirflx,   &
                      swdvisflx,aeradfo,aeradfos,asaeradfo,asaeradfos)
      else
        call rrtmg_sw_nomcica(npr,kth,icld,ldirect,play,plev,tlay,tlev, &
                              tsfc,h2ovmr,o3vmr,co2vmrk,ch4vmr,  &
                              n2ovmr,o2vmr,asdir,asdif,aldir,    &
                              aldif,czen,adjes,0,solcon,inflgsw, &
                              iceflgsw,liqflgsw,cldf,tauc,       &
                              ssac,asmc,fsfc,ciwp,clwp,rei,rel,  &
                              tauaer,ssaaer,asmaer,ecaer,        &
                              swuflx,swdflx,swhr,swuflxc,        &
                              swdflxc,swhrc,swddiruviflx,        &
                              swddifuviflx,swddirpirflx,         &
                              swddifpirflx,swdvisflx,            &
                              aeradfo,aeradfos,                  &
                              asaeradfo,asaeradfos)
      end if
    end if ! end shortwave call

    ! LW call :
    if ( imcica == 1 ) then
      permuteseed = permuteseed+mypid+ngptsw
      if ( permuteseed < 0 ) permuteseed = 2147483641+permuteseed
      call mcica_subcol_lw(npr,kth,icld,permuteseed,irng,play,        &
                           cldf,ciwp,clwp,rei,rel,tauc_lw,cldfmcl_lw, &
                           ciwpmcl_lw,clwpmcl_lw,reicmcl,relqmcl,     &
                           taucmcl_lw)
      call rrtmg_lw(npr,kth,icld,0,lradfor,ldirect,play,plev,tlay,tlev, &
                    tsfc,h2ovmr,o3vmr,co2vmrk,ch4vmr,n2ovmr,o2vmr,      &
                    cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,emis_surf,       &
                    inflglw,iceflglw,liqflglw,cldfmcl_lw,               &
                    taucmcl_lw,ciwpmcl_lw,clwpmcl_lw,reicmcl,           &
                    relqmcl,tauaer_lw,lwuflx,lwdflx,lwhr,lwuflxc,       &
                    lwdflxc,lwhrc,duflx_dt,duflxc_dt,aerlwfo,           &
                    aerlwfos,asaerlwfo,asaerlwfos)
    else
      call rrtmg_lw_nomcica(npr,kth,icld,0,play,plev,tlay,tlev,tsfc,       &
                            h2ovmr,o3vmr,co2vmrk,ch4vmr,n2ovmr,o2vmr,      &
                            cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,emis_surf,  &
                            inflglw,iceflglw,liqflglw,cldf,tauc,ciwp,clwp, &
                            rei,rel,tauaer_lw,lwuflx,lwdflx,lwhr,lwuflxc,  &
                            lwdflxc,lwhrc,duflx_dt,duflxc_dt)
    end if
    ! Output and interface
    !
    ! EES  next 3 added, they are calculated in radcsw : but not used further
    ! fsnirt   - Near-IR flux absorbed at toa
    ! fsnrtc   - Clear sky near-IR flux absorbed at toa
    ! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
    ! fsds     - Flux Shortwave Downwelling Surface

    sabtp(:)  = swdflx(:,kth) - swuflx(:,kth)
    solin(:)  = swdflx(:,kth)
    solout(:) = swuflx(:,kth)
    frsa(:)   = swdflx(:,1) - swuflx(:,1)
    clrst(:)  = swdflxc(:,kth) - swuflxc(:,kth)
    clrss(:)  = swdflxc(:,1) - swuflxc(:,1)

    firtp(:)  = -d_one * (lwdflx(:,kth) - lwuflx(:,kth))
    lwout(:)  = -lwuflx(:,kth)
    lwin(:)   = -lwdflx(:,kth)
    frla(:)   = -d_one * (lwdflx(:,1) - lwuflx(:,1))
    clrlt(:)  = -d_one * (lwdflxc(:,kth) - lwuflxc(:,kth))
    clrls(:)  = -d_one * (lwdflxc(:,1) - lwuflxc(:,1))

    ! coupling with BATS
    !  r2m%sabveg set to frsa (as in standard version : potential inconsistency
    !  if soil fraction is large)
    ! solar is normally the visible band only total incident surface flux
    abv(:) = frsa(:) + frla(:)
    sol(:) = swdvisflx(:,1)
    ! surface SW incident
    sols(:)  = swddiruviflx(:,1)
    solsd(:) = swddifuviflx(:,1)
    soll(:)  = swddirpirflx(:,1)
    solld(:) = swddifpirflx(:,1)
    ! LW incident
    slwd(:)  = lwdflx(:,1)

    ! 3d heating rate back on regcm grid and converted to K.S-1
    if ( idiag == 1 ) then
      do k = 1 , kz
        kj = kzp1-k
        o3(:,kj) = o3vmr(:,k)
        qrs(:,kj) = swhr(:,k) / secpd
        qrl(:,kj) = lwhr(:,k) / secpd
      end do
    end if
    do k = 1 , kz
      kj = kzp1-k
      dzr(:,kj) = deltaz(:,k)
      clwp_int(:,kj) = clwp(:,k)
    end do

    kmaxcld = 1+ncld
    kmincld = kz-ncld
    cld_int(:,:) = d_zero
    do k = kmaxcld , kmincld
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( clwp_int(n,k) > d_zero ) then
            cld_int(n,k) = m2r%cldfrc(j,i,k-1)+m2r%cldfrc(j,i,k) - &
                          (m2r%cldfrc(j,i,k-1)*m2r%cldfrc(j,i,k))
            cld_int(n,k) = min(cld_int(n,k),cftotmax)
          end if
          n = n + 1
        end do
      end do
    end do

    ! Calculate cloud parameters
    do n = 1 , npr
      totcf(n) = d_one
      if ( luse_max_rnovl ) then
        do k = 2 , kzp1
          totcf(n) = totcf(n) * (d_one - max(cld_int(n,k-1),cld_int(n,k)))/ &
                                (d_one - cld_int(n,k-1))
        end do
      else
        do k = 1 , kz
          totcf(n) = totcf(n)*(d_one-cld_int(n,k))
        end do
      end if
      totcf(n) = d_one - totcf(n)
    end do
    !totcf(:) = d_half * ( totcf(:) + maxval(cld_int(:,:),2) )
    totwv(:) = d_zero
    totci(:) = d_zero
    totcl(:) = d_zero
    do k = 1 , kz
      kj = kzp1-k
      do n = 1 , npr
        totci(n) = totci(n) + &
           clwp_int(n,k)*cld_int(n,k)*fice(n,kj)*d_r1000
        totcl(n) = totcl(n) + &
           clwp_int(n,k)*cld_int(n,k)*(d_one-fice(n,kj))*d_r1000
        totwv(n) =  totwv(n) + &
          h2ommr(n,k)*(play(n,k)*d_100)/(rgas*tlay(n,k))*deltaz(n,k)
      end do
    end do

    ! Finally call radout for coupling to BATS/CLM/ATM and outputing fields
    ! in RAD files

    call radout(lout,solin,solout,frsa,clrst,clrss,qrs,lwout,       &
                frla,clrlt,clrls,qrl,slwd,sols,soll,solsd,          &
                solld,totcf,totwv,totcl,totci,cld_int,clwp_int,abv, &
                sol,aeradfo,aeradfos,aerlwfo,aerlwfos,tauxar3d,     &
                tauasc3d,gtota3d,dzr,o3,outtaucl,outtauci,          &
                asaeradfo,asaeradfos,asaerlwfo,asaerlwfos,r2m,m2r)

  end subroutine rrtmg_driver

  subroutine prep_dat_rrtm(m2r,iyear,imonth)
    implicit none
    type(mod_2_rad) , intent(in) :: m2r
    integer(ik4) , intent(in) :: iyear , imonth
    integer(ik4) :: i , j , k , kj , ns , n , itr , kmaxcld , kmincld
    real(rkx) , parameter :: verynearone = 0.999999_rkx
    real(rkx) :: tmp1l , tmp2l , tmp3l , tmp1i , tmp2i , tmp3i
    real(rkx) :: w1 , w2 , p1 , p2
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
    real(rkx) , dimension(4) ::  abari , abarl , bbari , bbarl , cbari , &
                                cbarl , dbari , dbarl , ebari , ebarl , &
                                fbari , fbarl
    real(rkx) :: abarii , abarli , bbarii , bbarli , cbarii , cbarli , &
                dbarii , dbarli , ebarii , ebarli , fbarii , fbarli

    data abarl / 2.817e-2_rkx ,  2.682e-2_rkx , 2.264e-2_rkx , 1.281e-2_rkx/
    data bbarl / 1.305e+0_rkx ,  1.346e+0_rkx , 1.454e+0_rkx , 1.641e+0_rkx/
    data cbarl /-5.620e-8_rkx , -6.940e-6_rkx , 4.640e-4_rkx , 0.201e+0_rkx/
    data dbarl / 1.630e-7_rkx ,  2.350e-5_rkx , 1.240e-3_rkx , 7.560e-3_rkx/
    data ebarl / 0.829e+0_rkx ,  0.794e+0_rkx , 0.754e+0_rkx , 0.826e+0_rkx/
    data fbarl / 2.482e-3_rkx ,  4.226e-3_rkx , 6.560e-3_rkx , 4.353e-3_rkx/
!
    data abari / 3.4480e-3_rkx , 3.4480e-3_rkx , 3.4480e-3_rkx , 3.44800e-3_rkx/
    data bbari / 2.4310e+0_rkx , 2.4310e+0_rkx , 2.4310e+0_rkx , 2.43100e+0_rkx/
    data cbari / 1.0000e-5_rkx , 1.1000e-4_rkx , 1.8610e-2_rkx , 0.46658e+0_rkx/
    data dbari / 0.0000e+0_rkx , 1.4050e-5_rkx , 8.3280e-4_rkx , 2.05000e-5_rkx/
    data ebari / 0.7661e+0_rkx , 0.7730e+0_rkx , 0.7940e+0_rkx , 0.95950e+0_rkx/
    data fbari / 5.8510e-4_rkx , 5.6650e-4_rkx , 7.2670e-4_rkx , 1.07600e-4_rkx/
    !
    ! define index pointing on appropriate parameter in slingo's table
    ! for eachRRTM SW band
    !
    data indsl /4,4,3,3,3,3,3,2,2,1,1,1,1,4 /
    ! CONVENTION : RRTMG driver takes layering from bottom to TOA.
    ! Regcm consider TOA to bottom

    ! surface pressure and scaled pressure, from which level are computed
    ! RRTM SW takes pressure in mb,hpa
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        dlat(n) = m2r%xlat(j,i)
        topz(n) = m2r%zq(j,i,1)-m2r%za(j,i,1)
        n = n + 1
      end do
    end do

    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        xptrop(n) = m2r%ptrop(j,i) * d_r100
        psfc(n) = m2r%psatms(j,i) * d_r100
        n = n + 1
      end do
    end do

    do k = 1 , kz
      n = 1
      kj = kzp1-k
      do i = ici1 , ici2
        do j = jci1 , jci2
          play(n,k) = m2r%phatms(j,i,kj)*d_r100
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
          plev(n,k) = m2r%pfatms(j,i,kj)*d_r100
          n = n + 1
        end do
      end do
    end do
    do k = kzp1+1 , ktf
      plev(:,k) = stdplevf(kclimf+k-kzp1-1)
    end do
    ! smooth transition from top to climato at kzp1
    play(:,kzp1) = d_half * (plev(:,kzp1) + plev(:,kzp2) )
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
          tlay(n,k) = m2r%tatms(j,i,kj)
          n = n + 1
        end do
      end do
    end do
    do k = kzp1 , kth
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          tlay(n,k) = max(tlay(n,kz), &
            stdatm_val(calday,dlat(n),play(n,k),istdatm_tempk))
          n = n + 1
        end do
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
          n = n + 1
        end do
      end do
    end do
    do k = kzp1 , kth
      deltaz(:,k) = (stdhlevf(kclimf+k-kzp1) - &
                     stdhlevf(kclimf+k-kzp1-1)) * d_1000
    end do
    !
    ! air temperature at the interface
    !
    tlev(:,1) = tsfc(:)
    if ( idynamic /= 3 ) then
      do k = 2 , kz
        kj = kzp1-k+1
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            p1 = (plev(n,k)/play(n,k-1))**c287
            p2 = (plev(n,k)/play(n,k))**c287
            w1 = (hsigma(kj) - sigma(kj)) / (hsigma(kj) - hsigma(kj-1))
            w2 = (sigma(kj) - hsigma(kj-1) ) / (hsigma(kj) - hsigma(kj-1))
            tlev(n,k) = w1 * tlay(n,k-1) * p1 + w2 * tlay(n,k) * p2
            n = n + 1
          end do
        end do
      end do
    else
      do k = 2 , kz
        kj = kzp1-k+1
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            tlev(n,k) = 0.5_rkx * (tlay(n,k-1) + tlay(n,k))
            n = n + 1
          end do
        end do
      end do
    end if
    do k = kzp1 , ktf
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          tlev(n,k) = max(tlev(n,kz), &
            stdatm_val(calday,dlat(n),plev(n,k),istdatm_tempk))
          n = n + 1
        end do
      end do
    end do

    if ( ipptls > 1 ) then
      do k = 1 , kz
        kj = kzp1 - k
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            ql1(n,k) = m2r%qxatms(j,i,kj,iqc)
            qi1(n,k) = m2r%qxatms(j,i,kj,iqi)
            n = n + 1
          end do
        end do
      end do
      do k = kzp1 , kth
        ql1(:,k) = d_zero
        qi1(:,k) = d_zero
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
          h2ommr(n,k) = max(1.0e-8_rkx,m2r%qxatms(j,i,kj,iqv))
          h2ovmr(n,k) = h2ommr(n,k) * rep2
          n = n + 1
        end do
      end do
    end do
    do k = kzp1 , kth
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          h2ommr(n,k) = 1.0e-10_rkx
          h2ovmr(n,k) = h2ommr(n,k) * rep2
          n = n + 1
        end do
      end do
    end do
    !
    ! o3 volume mixing ratio
    !
    do k = 1 , kz
      n = 1
      kj = kzp1 - k
      do i = ici1 , ici2
        do j = jci1 , jci2
          o3vmr(n,k) = d_half*(o3prof(j,i,kj)+o3prof(j,i,kj+1)) * amd/amo3
          n = n + 1
        end do
      end do
    end do
    do k = kzp1 , kth
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          o3vmr(n,k) = &
            stdatm_val(calday,dlat(n),play(n,k),istdatm_ozone) / &
            stdatm_val(calday,dlat(n),play(n,k),istdatm_airdn) * amd/amo3
          n = n + 1
        end do
      end do
    end do
    !
    ! cgas is in ppm , ppb , ppt
    !
    ! Transform in mass mixing ratios (g/g) for trcmix
    !
    do n = 1 , npr
      co2vmr(n) = ghgval(igh_co2,iyear,imonth,dlat(n))
      co2mmr(n) = co2vmr(n)*(amco2/amd)
      ch40(n) = ghgval(igh_ch4,iyear,imonth,dlat(n))*(amch4/amd)
      n2o0(n) = ghgval(igh_n2o,iyear,imonth,dlat(n))*(amn2o/amd)
      cfc110(n) = ghgval(igh_cfc11,iyear,imonth,dlat(n))*(amcfc11/amd)
      cfc120(n) = ghgval(igh_cfc12,iyear,imonth,dlat(n))*(amcfc12/amd)
    end do

    do k = 1 , kz
      co2vmrk(:,k) = co2vmr(:)
    end do

    call trcmix(1,npr,dlat,xptrop,play,n2ommr,ch4mmr,cfc11mmr,cfc12mmr)

    do k = 1 , kth
      do n = 1 , npr
        !
        ! Form mass mixing ratios to vomlume mixing ratios
        !
        o2vmr(n,k)    = 0.209460_rkx
        n2ovmr(n,k)   = n2ommr(n,k) * (amd/amn2o)
        ch4vmr(n,k)   = ch4mmr(n,k) * (amd/amch4)
        cfc11vmr(n,k) = cfc11mmr(n,k) * (amd/amcfc11)
        cfc12vmr(n,k) = cfc12mmr(n,k) * (amd/amcfc12)
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
        kj = kzp1-k+1
        do k = 1 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              aermmr(n,k,itr) = m2r%chiatms(j,i,kj,itr)
              n = n + 1
            end do
          end do
        end do
      end do
    end if
    if ( ichem == 1 .or. iclimaaer > 0 ) then
      do k = 1 , kz
        n = 1
        kj = kz-k+1
        do i = ici1 , ici2
          do j = jci1 , jci2
            rh(n,k) = m2r%rhatms(j,i,kj)
            n = n + 1
          end do
        end do
      end do
      do k = 1 , kzp1
        n = 1
        kj = kzp1-k+1
        do i = ici1 , ici2
          do j = jci1 , jci2
            pint(n,k) = m2r%pfatms(j,i,kj)
            n = n + 1
          end do
        end do
      end do
      call aeroppt(rh,pint,1,npr)
      ! adapt reverse the vertical grid for RRTM
      do i = 1 , nbndsw
        do k = 1 , kth
          kj = kth + 1 - k
          do n = 1 , npr
            ecaer(n,k,i)  = 0.78_rkx ! not used
            tauaer(n,k,i) = tauxar3d(n,kj,i)
            ssaaer(n,k,i) = tauasc3d(n,kj,i)
            asmaer(n,k,i) = gtota3d(n,kj,i)
          end do
        end do
      end do
      do i = 1 , nbndlw
        do k = 1 , kz
          kj = kzp1 - k
          do n = 1 , npr
            tauaer_lw(n,k,i) = tauxar3d_lw(n,kj,i)
          end do
        end do
      end do
    else
      ecaer(:,:,:)  = 0.78_rkx ! not used
      tauaer(:,:,:) = d_zero
      ssaaer(:,:,:) = d_zero
      asmaer(:,:,:) = d_zero
      tauaer_lw(:,:,:) = d_zero
    end if
    !
    ! cloud fraction and cloud liquid waterpath calculation:
    ! as in STANDARD SCHEME for now (getdat) : We need to improve this
    ! according to new cloud microphysics!
    ! fractional cloud cover (dependent on relative humidity)
    !
    ! qc   = gary's mods for clouds/radiation tie-in to exmois
    !
    kmaxcld = 1+ncld
    kmincld = kz-ncld
    cldf = d_zero
    clwp = d_zero
    do k = kmaxcld , kmincld
      kj = kzp1 - k
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          !
          ! convert liquid water content into liquid water path, i.e.
          ! multiply b deltaz
          ! deltaz,clwp are on the right grid since plev and tlay are
          ! care pressure is on bottom/toa grid
          !
          cldf(n,k) = min(m2r%cldfrc(j,i,kj),cftotmax)
          clwp(n,k) = m2r%cldlwc(j,i,kj) * deltaz(n,k)
          n = n + 1
        end do
      end do
    end do
    !
    ! fabtest:  set cloud fractional cover at top model level = 0
    ! as in std scheme
    cldf(:,kzp1:kth) = d_zero
    clwp(:,kzp1:kth) = d_zero
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
    call cldefr_rrtm(tlay,play,rel,rei,fice)
    !
    ! partition of total water path betwwen liquide and ice.
    ! now clwp is liquide only !
    !
    do k = 1 , kth
      do n = 1 , npr
        ciwp(n,k) =  clwp(n,k) * fice(n,k)
        clwp(n,k) =  clwp(n,k) * (d_one - fice(n,k))
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
    asmc  =  0.850_rkx
    fsfc  =  0.725_rkx

#if defined(CLM45) || defined(CLM)
    ! TS is the radiant temeprature
    emis_surf(:,:) = 1.0_rkx
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

    if ( inflgsw == 0 ) then
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
            ! if McICA approach is used (irng = 0, or 1) ,
            ! cloud optical depth (extinction) is in cloud quantity.
            ! if not, RRTM receive grid levl cloud optical depth and
            ! optical properties, as in the standard scheme
            if ( imcica == 1 ) then
              ! cloud optical properties extinction optical depth
              !              IN_CLOUD quantities !
              !
              tauxcl(n,k) = clwp(n,k)*tmp1l
              tauxci(n,k) = ciwp(n,k)*tmp1i
              outtaucl(n,k,indsl(ns)) = outtaucl(n,k,indsl(ns)) + tauxcl(n,k)
              outtauci(n,k,indsl(ns)) = outtauci(n,k,indsl(ns)) + tauxci(n,k)
            else
              !
              ! use here the same parametrisation as in standard scheme
              ! in mod_radiatiom. Take care: grid level quantity !
              ! ciwp and clwp are already calculated
              !
              tauxcl(n,k) = clwp(n,k)*tmp1l*cldf(n,k) / &
                          (d_one+(d_one-0.85_rkx)*(d_one-cldf(n,k))*      &
                          clwp(n,k)*tmp1l)
              outtaucl(n,k,indsl(ns)) = outtaucl(n,k,indsl(ns)) + tauxcl(n,k)
              tauxci(n,k) = ciwp(n,k)*tmp1i*cldf(n,k) /     &
                          (d_one+(d_one-0.78_rkx)*(d_one-cldf(n,k)) * &
                          ciwp(n,k)*tmp1i)
              outtauci(n,k,indsl(ns)) = outtauci(n,k,indsl(ns)) + tauxci(n,k)
            end if

            wcl(n,k) = min(tmp2l,verynearone)
            gcl(n,k) = ebarli + tmp3l
            fcl(n,k) = gcl(n,k)*gcl(n,k)
            wci(n,k) = min(tmp2i,verynearone)
            gci(n,k) = ebarii + tmp3i
            fci(n,k) = gci(n,k)*gci(n,k)
            tauc(ns,n,k) = tauxcl(n,k) + tauxci(n,k)
            ssac(ns,n,k) = (tauxcl(n,k) * wcl(n,k) + &
                            tauxci(n,k) * wci(n,k) ) / tauc (ns,n,k)
            asmc(ns,n,k) = (tauxcl(n,k) * gcl(n,k) + &
                            tauxci(n,k) * gci(n,k) ) / tauc (ns,n,k)
            fsfc(ns,n,k) = (tauxcl(n,k) * fcl(n,k) + &
                            tauxci(n,k) * fci(n,k) ) / tauc (ns,n,k)
          end do
        end do
      end do ! spectral loop
      !
      ! if mcica is not enabled, the column is either clear sky
      ! either overcast with effective cloud optical properties
      ! similar to std rad (accounting for cld fraction).
      ! Needs to set cldf to 0 or 1 for the the rrtm_nomcica code to work.
      !
      if ( imcica == 0 ) then
        do k = 1 , kz
          do n = 1 , npr
            if ( cldf(n,k) > dlowval ) then
              cldf(n,k) = d_one
            else
              cldf(n,k) = d_zero
            end if
          end do
        end do
     end if
   end if  ! inflagsw
  end subroutine prep_dat_rrtm
  !
  ! For now we use for RRTM the same param as in standard rad
  !
  subroutine cldefr_rrtm(t,pmid,rel,rei,fice)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: pmid , t
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: fice , rei , rel
    integer(ik4) :: k , n
    real(rkx) :: pnrml , weight
    ! real(rkx) :: tpara
    ! reimax - maximum ice effective radius
    real(rkx) , parameter :: reimax = 30.0_rkx
    ! rirnge - range of ice radii (reimax - 10 microns)
    real(rkx) , parameter :: rirnge = 20.0_rkx
    ! pirnge - nrmlzd pres range for ice particle changes
    real(rkx) , parameter :: pirnge = 0.4_rkx
    ! picemn - normalized pressure below which rei=reimax
    real(rkx) , parameter :: picemn = 0.4_rkx
    ! Temperatures in K (263.16 , 243.16)
    real(rkx) , parameter :: minus10 = wattp-10.0_rkx
    real(rkx) , parameter :: minus30 = wattp-30.0_rkx

    do k = 1 , kth
      do n = 1 , npr
        !
        ! Define liquid drop size
        ! rel : liquid effective drop size (microns)
        !
        ! Global distribution of Cloud Droplet Effective Radius from
        ! POLDER polarization measurements
        ! On average, droplets are 2 to 3 micron smaller overland than over
        ! the ocean, but the smaller droplets are found over highly
        ! polluted regions and in areas affected by smoke from biomass
        ! burning activity. Largest droplets are found in remote tropical
        ! oceans, away from polarized reflectance of liquid clouds.
        ! Toward the poles, for colder temperature, large droplets freeze
        ! and only the smaller ones are kept liquid.
        !
        ! tpara = min(d_one,max(d_zero,(minus10-tm1(n,k))*0.05_rkx))
        !
        if ( ioro(n) == 1 ) then
          ! Effective liquid radius over land
          ! rel(n,k) = 6.0_rkx + 5.0_rkx * tpara
          rel(n,k) = 8.50_rkx
        else
          ! Effective liquid radius over ocean and sea ice
          ! rel(n,k) = 7.0_rkx + 5.0_rkx * tpara
          rel(n,k) = 11.0_rkx
        end if
        !
        ! Determine rei as function of normalized pressure
        !
        pnrml = pmid(n,k)/psfc(n)
        weight = max(min((pnrml-picemn)/pirnge,d_one),d_zero)
        rei(n,k) = reimax - rirnge*weight
        !
        ! Define fractional amount of cloud that is ice
        !
        if ( ipptls > 1 ) then
          if ( qi1(n,k) > minqq ) then
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
            fice(n,k) = (minus10-t(n,k))/20.0_rkx
          else
            ! if colder than -30 degrees C then ice phase
            fice(n,k) = d_one
          end if
        end if
        ! Turn off ice radiative properties by setting fice = 0.0
        ! fice(n,k) = d_zero
      end do
    end do
  end subroutine cldefr_rrtm

end module mod_rrtmg_driver
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

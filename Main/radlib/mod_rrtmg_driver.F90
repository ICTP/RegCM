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
  use mod_date
  use mod_stdatm
  use mod_isaatm
  use mcica_subcol_gen_sw
  use parrrsw
  use rrsw_wvn
  use parrrtm
  use rrtmg_sw_rad
  use rrtmg_lw_rad
  use rrtmg_lw_rad_nomcica
  use rrtmg_sw_rad_nomcica
  use mod_rad_outrad
  use mod_mpmessage
  use mod_runparams
  use mod_regcm_types

  implicit none

  private
!
  public :: allocate_mod_rad_rrtmg, rrtmg_driver

  real(rk8), pointer, contiguous, dimension(:) :: frsa => null( )
  real(rk8), pointer, contiguous, dimension(:) :: sabtp => null( )
  real(rk8), pointer, contiguous, dimension(:) :: clrst => null( )
  real(rk8), pointer, contiguous, dimension(:) :: solin => null( )
  real(rk8), pointer, contiguous, dimension(:) :: solout => null( )
  real(rk8), pointer, contiguous, dimension(:) :: clrss => null( )
  real(rk8), pointer, contiguous, dimension(:) :: firtp => null( )
  real(rk8), pointer, contiguous, dimension(:) :: lwout => null( )
  real(rk8), pointer, contiguous, dimension(:) :: lwin => null( )
  real(rk8), pointer, contiguous, dimension(:) :: frla => null( )
  real(rk8), pointer, contiguous, dimension(:) :: clrlt => null( )
  real(rk8), pointer, contiguous, dimension(:) :: clrls => null( )
  real(rk8), pointer, contiguous, dimension(:) :: abv => null( )
  real(rk8), pointer, contiguous, dimension(:) :: sol => null( )
  real(rk8), pointer, contiguous, dimension(:) :: sols => null( )
  real(rk8), pointer, contiguous, dimension(:) :: soll => null( )
  real(rk8), pointer, contiguous, dimension(:) :: solsd => null( )
  real(rk8), pointer, contiguous, dimension(:) :: solld => null( )
  real(rk8), pointer, contiguous, dimension(:) :: slwd => null( )
  real(rk8), pointer, contiguous, dimension(:) :: tsfc => null( )
  real(rk8), pointer, contiguous, dimension(:) :: psfc => null( )
  real(rk8), pointer, contiguous, dimension(:) :: asdir => null( )
  real(rk8), pointer, contiguous, dimension(:) :: asdif => null( )
  real(rk8), pointer, contiguous, dimension(:) :: aldir => null( )
  real(rk8), pointer, contiguous, dimension(:) :: aldif => null( )
  real(rk8), pointer, contiguous, dimension(:) :: czen => null( )
  real(rk8), pointer, contiguous, dimension(:) :: dlat => null( )
  real(rk8), pointer, contiguous, dimension(:) :: xptrop => null( )
  real(rk8), pointer, contiguous, dimension(:) :: totcf => null( )
  real(rk8), pointer, contiguous, dimension(:) :: totwv => null( )
  real(rk8), pointer, contiguous, dimension(:) :: totcl => null( )
  real(rk8), pointer, contiguous, dimension(:) :: totci => null( )

  real(rk8), pointer, contiguous, dimension(:,:) :: qrs => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: qrl => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: clwp_int => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: pint => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: rh => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: cld_int => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: tlay => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: h2ovmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: o3vmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: co2vmrk => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: play => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: ch4vmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: n2ovmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: o2vmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: cfc11vmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: cfc12vmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: cfc22vmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: ccl4vmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: o3 => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: reicmcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: relqmcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swhr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swhrc => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: ciwp => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: clwp => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: rei => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: rel => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: cldf => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: lwhr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: lwhrc => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: duflx_dt => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: duflxc_dt => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: ql1 => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: qi1 => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: alpha => null( )

  real(rk8), pointer, contiguous, dimension(:,:) :: plev => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: tlev => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swuflx => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swdflx => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swuflxc => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swdflxc => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: lwuflx => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: lwdflx => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: lwuflxc => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: lwdflxc => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swddiruviflx => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swddifuviflx => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swddirpirflx => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swddifpirflx => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: swdvisflx => null( )

  real(rk8), pointer, contiguous, dimension(:,:,:) :: cldfmcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: taucmcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: ssacmcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: asmcmcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: fsfcmcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: ciwpmcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: clwpmcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: cldfmcl_lw => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: taucmcl_lw => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: ciwpmcl_lw => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: clwpmcl_lw => null( )

  real(rk8), pointer, contiguous, dimension(:) :: aeradfo => null( )
  real(rk8), pointer, contiguous, dimension(:) :: aeradfos => null( )
  real(rk8), pointer, contiguous, dimension(:) :: asaeradfo => null( )
  real(rk8), pointer, contiguous, dimension(:) :: asaeradfos => null( )
  real(rk8), pointer, contiguous, dimension(:) :: aerlwfo => null( )
  real(rk8), pointer, contiguous, dimension(:) :: aerlwfos => null( )
  real(rk8), pointer, contiguous, dimension(:) :: asaerlwfo => null( )
  real(rk8), pointer, contiguous, dimension(:) :: asaerlwfos => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: fice => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: wcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: wci => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: gcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: gci => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: fcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: fci => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: tauxcl => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: tauxci => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: h2ommr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: n2ommr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: ch4mmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: cfc11mmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: cfc12mmr => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: deltaz => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: dzr => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: outtaucl => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: outtauci => null( )

  integer(ik4), pointer, contiguous, dimension(:) :: ioro => null( )

  ! spectral dependant quantities

  real(rk8), pointer, contiguous, dimension(:,:,:) :: tauaer => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: ssaaer => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: asmaer => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: ecaer => null( )
  ! tauc = in-cloud optical depth
  ! ssac = in-cloud single scattering albedo (non-delta scaled)
  ! asmc = in-cloud asymmetry parameter (non-delta scaled)
  ! fsfc = in-cloud forward scattering fraction (non-delta scaled)
  real(rk8), pointer, contiguous, dimension(:,:,:) :: tauc => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: ssac => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: asmc => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: fsfc => null( )
  real(rk8), pointer, contiguous, dimension(:,:) :: emis_surf => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: tauc_lw => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: tauaer_lw => null( )

  real(rk8), pointer, contiguous, dimension(:) :: cfc110 => null( )
  real(rk8), pointer, contiguous, dimension(:) :: cfc120 => null( )
  real(rk8), pointer, contiguous, dimension(:) :: ch40 => null( )
  real(rk8), pointer, contiguous, dimension(:) :: co2mmr => null( )
  real(rk8), pointer, contiguous, dimension(:) :: co2vmr => null( )
  real(rk8), pointer, contiguous, dimension(:) :: n2o0 => null( )

  integer(ik4) :: npr, npj

  integer(ik4) :: permuteseed = 37_ik4

  logical, parameter :: luse_max_rnovl = .true.

  contains

  subroutine allocate_mod_rad_rrtmg
    implicit none
    npj = (jci2-jci1+1)
    npr = npj*(ici2-ici1+1)
    call getmem(frsa,1,npr,'rrtmg:frsa')
    call getmem(sabtp,1,npr,'rrtmg:sabtp')
    call getmem(clrst,1,npr,'rrtmg:clrst')
    call getmem(clrss,1,npr,'rrtmg:clrss')
    call getmem(firtp,1,npr,'rrtmg:firtp')
    call getmem(lwout,1,npr,'rrtmg:lwout')
    call getmem(lwin,1,npr,'rrtmg:lwin')
    call getmem(frla,1,npr,'rrtmg:frla')
    call getmem(clrlt,1,npr,'rrtmg:clrlt')
    call getmem(clrls,1,npr,'rrtmg:clrls')
    call getmem(solin,1,npr,'rrtmg:solin')
    call getmem(solout,1,npr,'rrtmg:solout')
    call getmem(sols,1,npr,'rrtmg:sols')
    call getmem(soll,1,npr,'rrtmg:soll')
    call getmem(solsd,1,npr,'rrtmg:solsd')
    call getmem(solld,1,npr,'rrtmg:solld')
    call getmem(slwd,1,npr,'rrtmg:slwd')
    call getmem(qrs,1,npr,1,kz,'rrtmg:qrs')
    call getmem(qrl,1,npr,1,kz,'rrtmg:qrl')
    call getmem(o3,1,npr,1,kz,'rrtmg:o3')
    call getmem(clwp_int,1,npr,1,kz,'rrtmg:clwp_int')
    call getmem(cld_int,1,npr,1,kzp1,'rrtmg:cld_int')
    call getmem(aeradfo,1,npr,'rrtmg:aeradfo')
    call getmem(aeradfos,1,npr,'rrtmg:aeradfos')
    call getmem(asaeradfo,1,npr,'rrtmg:aseradfo')
    call getmem(asaeradfos,1,npr,'rrtmg:aseradfos')
    call getmem(aerlwfo,1,npr,'rrtmg:aerlwfo')
    call getmem(aerlwfos,1,npr,'rrtmg:aerlwfos')
    call getmem(asaerlwfo,1,npr,'rrtmg:asaerlwfo')
    call getmem(asaerlwfos,1,npr,'rrtmg:asaerlwfos')
    call getmem(abv,1,npr,'rrtmg:sol')
    call getmem(sol,1,npr,'rrtmg:sol')
    call getmem(tsfc,1,npr,'rrtmg:tsfc')
    call getmem(psfc,1,npr,'rrtmg:psfc')
    call getmem(asdir,1,npr,'rrtmg:asdir')
    call getmem(asdif,1,npr,'rrtmg:asdif')
    call getmem(aldir,1,npr,'rrtmg:aldir')
    call getmem(aldif,1,npr,'rrtmg:aldif')
    call getmem(czen,1,npr,'rrtmg:czen')
    call getmem(dlat,1,npr,'rrtmg:dlat')
    call getmem(xptrop,1,npr,'rrtmg:xptrop')
    call getmem(ioro,1,npr,'rrtmg:ioro')
    call getmem(totcf,1,npr,'rrtmg:totcf')
    call getmem(totwv,1,npr,'rrtmg:totwv')
    call getmem(totcl,1,npr,'rrtmg:totlf')
    call getmem(totci,1,npr,'rrtmg:totif')
    call getmem(co2mmr,1,npr,'rrtmg:co2mmr')
    call getmem(co2vmr,1,npr,'rrtmg:co2vmr')
    call getmem(n2o0,1,npr,'rrtmg:n2o0')
    call getmem(ch40,1,npr,'rrtmg:ch40')
    call getmem(cfc110,1,npr,'rrtmg:cfc110')
    call getmem(cfc120,1,npr,'rrtmg:cfc120')

    if ( ichem == 1 .or. iclimaaer > 0 ) then
      call getmem(pint,1,npr,1,kzp1,'rrtmg:pint')
      call getmem(rh,1,npr,1,kz,'rrtmg:rh')
    end if
    call getmem(play,1,npr,1,kth,'rrtmg:play')
    call getmem(tlay,1,npr,1,kth,'rrtmg:tlay')
    call getmem(alpha,1,npr,1,kth,'rrtmg:alpha')
    if ( ipptls > 1 ) then
      call getmem(ql1,1,npr,1,kth,'rrtmg:ql1')
      call getmem(qi1,1,npr,1,kth,'rrtmg:qi1')
    end if
    call getmem(h2ovmr,1,npr,1,kth,'rrtmg:h2ovmr')
    call getmem(o3vmr,1,npr,1,kth,'rrtmg:o3vmr')
    call getmem(co2vmrk,1,npr,1,kth,'rrtmg:co2vmrk')
    call getmem(ch4vmr,1,npr,1,kth,'rrtmg:ch4vmr')
    call getmem(n2ovmr,1,npr,1,kth,'rrtmg:n2ovmr')
    call getmem(o2vmr,1,npr,1,kth,'rrtmg:o2vmr')
    call getmem(cfc11vmr,1,npr,1,kth,'rrtmg:cfc11vmr')
    call getmem(cfc12vmr,1,npr,1,kth,'rrtmg:cfc12vmr')
    call getmem(cfc22vmr,1,npr,1,kth,'rrtmg:cfc22vmr')
    call getmem(ccl4vmr,1,npr,1,kth,'rrtmg:ccl4vmr')
    call getmem(reicmcl,1,npr,1,kth,'rrtmg:reicmcl')
    call getmem(relqmcl,1,npr,1,kth,'rrtmg:relqmcl')
    call getmem(swhr,1,npr,1,kth,'rrtmg:swhr')
    call getmem(swhrc,1,npr,1,kth,'rrtmg:swhrc')
    call getmem(ciwp,1,npr,1,kth,'rrtmg:ciwp')
    call getmem(clwp,1,npr,1,kth,'rrtmg:clwp')
    call getmem(rei,1,npr,1,kth,'rrtmg:rei')
    call getmem(rel,1,npr,1,kth,'rrtmg:rel')
    call getmem(cldf,1,npr,1,kth,'rrtmg:cldf')
    call getmem(lwhr,1,npr,1,kth,'rrtmg:lwhr')
    call getmem(lwhrc,1,npr,1,kth,'rrtmg:lwhrc')
    call getmem(duflx_dt,1,npr,1,kth,'rrtmg:duflx_dt')
    call getmem(duflxc_dt,1,npr,1,kth,'rrtmg:duflxc_dt')
    call getmem(plev,1,npr,1,ktf,'rrtmg:plev')
    call getmem(tlev,1,npr,1,ktf,'rrtmg:tlev')
    call getmem(swuflx,1,npr,1,ktf,'rrtmg:swuflx')
    call getmem(swdflx,1,npr,1,ktf,'rrtmg:swdflx')
    call getmem(swuflxc,1,npr,1,ktf,'rrtmg:swuflxc')
    call getmem(swdflxc,1,npr,1,ktf,'rrtmg:swdflxc')
    call getmem(lwuflx,1,npr,1,ktf,'rrtmg:lwuflx')
    call getmem(lwdflx,1,npr,1,ktf,'rrtmg:lwdflx')
    call getmem(lwuflxc,1,npr,1,ktf,'rrtmg:lwuflxc')
    call getmem(lwdflxc,1,npr,1,ktf,'rrtmg:lwdflxc')
    call getmem(swddiruviflx,1,npr,1,ktf,'rrtmg:swddiruviflx')
    call getmem(swddifuviflx,1,npr,1,ktf,'rrtmg:swddifuviflx')
    call getmem(swddirpirflx,1,npr,1,ktf,'rrtmg:swddirpirflx')
    call getmem(swddifpirflx,1,npr,1,ktf,'rrtmg:swddifpirflx')
    call getmem(swdvisflx,1,npr,1,ktf,'rrtmg:swdvisflx')
    call getmem(cldfmcl,1,ngptsw,1,npr,1,kth,'rrtmg:cldfmcl')
    call getmem(taucmcl,1,ngptsw,1,npr,1,kth,'rrtmg:taucmcl')
    call getmem(ssacmcl,1,ngptsw,1,npr,1,kth,'rrtmg:ssacmcl')
    call getmem(asmcmcl,1,ngptsw,1,npr,1,kth,'rrtmg:asmcmcl')
    call getmem(fsfcmcl,1,ngptsw,1,npr,1,kth,'rrtmg:fsfcmcl')
    call getmem(ciwpmcl,1,ngptsw,1,npr,1,kth,'rrtmg:ciwpmcl')
    call getmem(clwpmcl,1,ngptsw,1,npr,1,kth,'rrtmg:clwpmcl')
    call getmem(cldfmcl_lw,1,ngptlw,1,npr,1,kth,'rrtmg:cldfmcl_lw')
    call getmem(taucmcl_lw,1,ngptlw,1,npr,1,kth,'rrtmg:taucmcl_lw')
    call getmem(ciwpmcl_lw,1,ngptlw,1,npr,1,kth,'rrtmg:ciwpmcl_lw')
    call getmem(clwpmcl_lw,1,ngptlw,1,npr,1,kth,'rrtmg:clwpmcl_lw')
    call getmem(tauaer,1,npr,1,kth,1,nbndsw,'rrtmg:tauaer')
    call getmem(ssaaer,1,npr,1,kth,1,nbndsw,'rrtmg:ssaaer')
    call getmem(asmaer,1,npr,1,kth,1,nbndsw,'rrtmg:asmaer')
    call getmem(ecaer,1,npr,1,kth,1,nbndsw,'rrtmg:ecaer')
    call getmem(tauc,1,nbndsw,1,npr,1,kth,'rrtmg:tauc')
    call getmem(ssac,1,nbndsw,1,npr,1,kth,'rrtmg:ssac')
    call getmem(asmc,1,nbndsw,1,npr,1,kth,'rrtmg:asmc')
    call getmem(fsfc,1,nbndsw,1,npr,1,kth,'rrtmg:fsfc')
    call getmem(emis_surf,1,npr,1,nbndlw,'rrtmg:emis_surf')
    call getmem(tauaer_lw,1,npr,1,kth,1,nbndlw,'rrtmg:tauaer_lw')
    call getmem(tauc_lw,1,nbndlw,1,npr,1,kth,'rrtmg:tauc_lw')
    call getmem(fice,1,npr,1,kth,'rrtmg:fice')
    call getmem(wcl,1,npr,1,kth,'rrtmg:wcl')
    call getmem(wci,1,npr,1,kth,'rrtmg:wci')
    call getmem(gcl,1,npr,1,kth,'rrtmg:gcl')
    call getmem(gci,1,npr,1,kth,'rrtmg:gci')
    call getmem(fcl,1,npr,1,kth,'rrtmg:fcl')
    call getmem(fci,1,npr,1,kth,'rrtmg:fci')
    call getmem(tauxcl,1,npr,1,kth,'rrtmg:tauxcl')
    call getmem(tauxci,1,npr,1,kth,'rrtmg:tauxci')
    call getmem(outtaucl,1,npr,1,kz,1,4,'rrtmg:outtaucl')
    call getmem(outtauci,1,npr,1,kz,1,4,'rrtmg:outtauci')
    call getmem(h2ommr,1,npr,1,kth,'rrtmg:h2ommr')
    call getmem(n2ommr,1,npr,1,kth,'rrtmg:n2ommr')
    call getmem(ch4mmr,1,npr,1,kth,'rrtmg:ch4mmr')
    call getmem(cfc11mmr,1,npr,1,kth,'rrtmg:cfc11mmr')
    call getmem(cfc12mmr,1,npr,1,kth,'rrtmg:cfc12mmr')
    call getmem(deltaz,1,npr,1,kth,'rrtmg:deltaz')
    call getmem(dzr,1,npr,1,kth,'rrtmg:dzr')

  end subroutine allocate_mod_rad_rrtmg

  subroutine rrtmg_driver(iyear,imonth,iday,lout,m2r,r2m)
    implicit none
    type(mod_2_rad), intent(in) :: m2r
    type(rad_2_mod), intent(inout) :: r2m
    integer(ik4), intent(in) :: iyear, imonth, iday
    logical, intent(in) :: lout
    integer(ik4) :: k, kj, n, i, j, ldirect
    integer(ik4) :: idcor           ! Decorrelation length type
    integer(ik4) :: juldat          ! Julian day of year
    real(rk8) :: decorr_con         ! decorrelation length, constant (m)
    logical :: lradfor
    real(rk8) :: adjes
    integer(ik4), parameter :: idrv = 0
    integer(ik4), parameter :: iaer = 10
    integer(ik4), parameter :: isolvar = -1

    ! from water path and cloud radius / tauc_LW is not requested
    tauc_lw(:,:,:) = 1.0e-10_rk8
    call prep_dat_rrtm(m2r,iyear,imonth)
    adjes = real(eccf,rk8)

    lradfor = ( rcmtimer%start( ) .or. syncro_radfor%will_act( ) )
    !
    ! Call to the shortwave radiation code as soon one element of czen is > 0.
    !
    swuflx(:,:) = 0.0_rk8
    swdflx(:,:) = 0.0_rk8
    swuflxc(:,:) = 0.0_rk8
    swdflxc(:,:) = 0.0_rk8
    swhr(:,:) = 0.0_rk8        !
    swhrc(:,:) = 0.0_rk8
    swdvisflx(:,:) = 0.0_rk8
    swddiruviflx(:,:) = 0.0_rk8
    swddifuviflx(:,:) = 0.0_rk8
    swddirpirflx(:,:) = 0.0_rk8
    swddifpirflx(:,:) = 0.0_rk8
    aeradfo(:) = 0.0_rk8
    aeradfos(:) = 0.0_rk8
    asaeradfo(:) = 0.0_rk8
    asaeradfos(:) = 0.0_rk8

    ! hanlde aerosol direct effect in function of ichem or iclimaaer
    !
    ldirect = 0
    if ( ichem == 1 .and. iaerosol == 1 ) then
      ldirect = idirect
    else if ( iclimaaer > 0 ) then
      ldirect = 2
    end if

    juldat = int(julianday(iyear, imonth, iday),ik4)
    alpha(:,:) = 0.0_rk8
    if ( icld == 4 .or. icld == 5 ) then
      idcor = 1
      decorr_con = 5000.0_rk8
      call get_alpha(npr,kth,icld,idcor,decorr_con,deltaz, &
                     dlat,juldat,cldf,alpha)
    end if

    do i = ici1, ici2
      do j = jci1, jci2
        n = (j-jci1+1)+(i-ici1)*npj
        asdir(n) = m2r%aldirs(j,i)
        asdif(n) = m2r%aldifs(j,i)
        aldir(n) = m2r%aldirl(j,i)
        aldif(n) = m2r%aldifl(j,i)
        czen(n)  = m2r%coszrs(j,i)
        if ( czen(n) < 1.0e-3_rk8 ) czen(n) = 0.0_rk8
      end do
    end do

    if ( any(m2r%coszrs > 1.0e-3_rk8) ) then
      if ( imcica == 1 ) then
        permuteseed = permuteseed + ngptlw*ngptsw*kz*nicross*njcross
        do while ( permuteseed < 0 )
          permuteseed = 2147483641+permuteseed
        end do
        ! generates cloud properties:
        call mcica_subcol_sw(npr,kth,icld,permuteseed,irng,play,    &
                             cldf,ciwp,clwp,rei,rel,tauc,ssac,asmc, &
                             fsfc,alpha,cldfmcl,ciwpmcl,clwpmcl,    &
                             reicmcl,relqmcl,taucmcl,ssacmcl,       &
                             asmcmcl,fsfcmcl)
        call rrtmg_sw(npr,kth,icld,iaer,lradfor,ldirect,play,plev,  &
                      tlay,tlev,tsfc,h2ovmr,o3vmr,co2vmrk,ch4vmr,   &
                      n2ovmr,o2vmr,asdir,asdif,aldir,aldif,czen,    &
                      adjes,-1,solcon,isolvar,inflgsw,iceflgsw,     &
                      liqflgsw,cldfmcl,taucmcl,ssacmcl,asmcmcl,     &
                      fsfcmcl,ciwpmcl,clwpmcl,reicmcl,relqmcl,      &
                      tauaer,ssaaer,asmaer,ecaer,swuflx,swdflx,     &
                      swhr,swuflxc,swdflxc,swhrc,swddiruviflx,      &
                      swddifuviflx,swddirpirflx,swddifpirflx,       &
                      swdvisflx,aeradfo,aeradfos,asaeradfo,asaeradfos)
      else
        call rrtmg_sw_nomcica(npr,kth,icld,iaer,lradfor,ldirect,play,   &
                              plev,tlay,tlev,tsfc,h2ovmr,o3vmr,         &
                              co2vmrk,ch4vmr,n2ovmr,o2vmr,asdir,        &
                              asdif,aldir,aldif,czen,adjes,-1,          &
                              solcon,isolvar,inflgsw,iceflgsw,liqflgsw, &
                              cldf,tauc,ssac,asmc,fsfc,ciwp,clwp,       &
                              rei,rel,tauaer,ssaaer,asmaer,ecaer,       &
                              swuflx,swdflx,swhr,swuflxc,swdflxc,       &
                              swhrc,swddiruviflx,swddifuviflx,          &
                              swddirpirflx,swddifpirflx,swdvisflx,      &
                              aeradfo,aeradfos,asaeradfo,asaeradfos)
      end if
    end if ! end shortwave call

    ! LW call :
    if ( imcica == 1 ) then
      permuteseed = permuteseed + ngptlw*ngptsw*kz*nicross*njcross
      do while ( permuteseed < 0 )
        permuteseed = 2147483641+permuteseed
      end do
      call mcica_subcol_lw(npr,kth,icld,permuteseed,irng,play,   &
                           cldf,ciwp,clwp,rei,rel,tauc_lw,alpha, &
                           cldfmcl_lw,ciwpmcl_lw,clwpmcl_lw,     &
                           reicmcl,relqmcl,taucmcl_lw)
      call rrtmg_lw(npr,kth,icld,idrv,lradfor,ldirect,play,plev,tlay,tlev, &
                    tsfc,h2ovmr,o3vmr,co2vmrk,ch4vmr,n2ovmr,o2vmr,         &
                    cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,emis_surf,          &
                    inflglw,iceflglw,liqflglw,cldfmcl_lw,                  &
                    taucmcl_lw,ciwpmcl_lw,clwpmcl_lw,reicmcl,              &
                    relqmcl,tauaer_lw,lwuflx,lwdflx,lwhr,lwuflxc,          &
                    lwdflxc,lwhrc,aerlwfo,aerlwfos,asaerlwfo,asaerlwfos,   &
                    duflx_dt,duflxc_dt)
    else
      call rrtmg_lw_nomcica(npr,kth,icld,idrv,lradfor,ldirect,play,plev, &
                            tlay,tlev,tsfc,h2ovmr,o3vmr,co2vmrk,ch4vmr,  &
                            n2ovmr,o2vmr,cfc11vmr,cfc12vmr,cfc22vmr,     &
                            ccl4vmr,emis_surf,inflglw,iceflglw,liqflglw, &
                            cldf,tauc,ciwp,clwp,rei,rel,tauaer_lw,lwuflx,&
                            lwdflx,lwhr,lwuflxc,lwdflxc,lwhrc,aerlwfo,   &
                            aerlwfos,asaerlwfo,asaerlwfos,duflx_dt,duflxc_dt)
    end if
    ! Output and interface
    !
    ! EES  next 3 added, they are calculated in radcsw : but not used further
    ! fsnirt   - Near-IR flux absorbed at toa
    ! fsnrtc   - Clear sky near-IR flux absorbed at toa
    ! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
    ! fsds     - Flux Shortwave Downwelling Surface

    do n = 1, npr
      sabtp(n)  = swdflx(n,kth) - swuflx(n,kth)
      solin(n)  = swdflx(n,kth)
      solout(n) = swuflx(n,kth)
      frsa(n)   = swdflx(n,1) - swuflx(n,1)
      clrst(n)  = swdflxc(n,kth) - swuflxc(n,kth)
      clrss(n)  = swdflxc(n,1) - swuflxc(n,1)

      firtp(n)  = -1.0_rk8 * (lwdflx(n,kth) - lwuflx(n,kth))
      lwout(n)  = lwuflx(n,kth)
      lwin(n)   = -lwdflx(n,kth)
      frla(n)   = -1.0_rk8 * (lwdflx(n,1) - lwuflx(n,1))
      clrlt(n)  = -1.0_rk8 * (lwdflxc(n,kth) - lwuflxc(n,kth))
      clrls(n)  = -1.0_rk8 * (lwdflxc(n,1) - lwuflxc(n,1))
    end do

    ! coupling with BATS
    !  r2m%sabveg set to frsa (as in standard version: potential inconsistency
    !  if soil fraction is large)
    ! solar is normally the visible band only total incident surface flux
    do n = 1, npr
      abv(n) = frsa(n) + frla(n)
      sol(n) = swdvisflx(n,1)
      ! surface SW incident
      sols(n)  = swddiruviflx(n,1)
      solsd(n) = swddifuviflx(n,1)
      soll(n)  = swddirpirflx(n,1)
      solld(n) = swddifpirflx(n,1)
      ! LW incident
      slwd(n)  = lwdflx(n,1)
    end do

    ! 3d heating rate back on regcm grid and converted to K.S-1
    ! used to calculate heatinge rate tendency
    !
    do k = 1, kz
      kj = kzp1 -k
      do n = 1, npr
        qrs(n,kj) = swhr(n,k) / secpd
        qrl(n,kj) = lwhr(n,k) / secpd
      end do
    end do

    do k = 1, kz
      kj = kzp1-k
      do n = 1, npr
        o3(n,kj) = o3vmr(n,k)
        dzr(n,kj) = deltaz(n,k)
        clwp_int(n,kj) = clwp(n,k)
      end do
    end do

    cld_int(:,:) = 0.0_rk8
    do k = 1, kz
      do i = ici1, ici2
        do j = jci1, jci2
          n = (j-jci1+1)+(i-ici1)*npj
          if ( clwp_int(n,k) > 0.0_rk8 ) then
            cld_int(n,k) = m2r%cldfrc(j,i,k-1)+m2r%cldfrc(j,i,k) - &
                          (m2r%cldfrc(j,i,k-1)*m2r%cldfrc(j,i,k))
            cld_int(n,k) = min(cld_int(n,k),cftotmax)
          end if
        end do
      end do
    end do

    ! Calculate cloud parameters
    do n = 1, npr
      totcf(n) = 1.0_rk8
      if ( luse_max_rnovl ) then
        do k = 2, kzp1
          totcf(n) = totcf(n) * &
            (1.0001_rk8 - max(cld_int(n,k-1),cld_int(n,k)))/ &
            (1.0001_rk8 - cld_int(n,k-1))
        end do
      else
        do k = 1, kz
          totcf(n) = totcf(n)*(1.0_rk8-cld_int(n,k))
        end do
      end if
      totcf(n) = 1.0_rk8 - totcf(n)
    end do
    !totcf(:) = 0.5_rk8 * ( totcf(:) + maxval(cld_int(:,:),2) )
    totwv(:) = 0.0_rk8
    totci(:) = 0.0_rk8
    totcl(:) = 0.0_rk8
    do k = 1, kz
      kj = kzp1-k
      do n = 1, npr
        totci(n) = totci(n) + &
           clwp_int(n,k)*cld_int(n,k)*fice(n,kj)*0.001_rk8
        totcl(n) = totcl(n) + &
           clwp_int(n,k)*cld_int(n,k)*(1.0_rk8-fice(n,kj))*0.001_rk8
        totwv(n) =  totwv(n) + &
          h2ommr(n,k)*(play(n,k)*100.0_rk8)/(rgas*tlay(n,k))*deltaz(n,k)
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
    type(mod_2_rad), intent(in) :: m2r
    integer(ik4), intent(in) :: iyear, imonth
    integer(ik4) :: i, j, k, kj, ns, n, itr
    real(rk8), parameter :: verynearone = 0.999999_rk8
    real(rk8) :: tmp1l, tmp2l, tmp3l, tmp1i, tmp2i, tmp3i
    real(rk8) :: w1, w2, p1, p2
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
    real(rk8) :: abarii, abarli, bbarii, bbarli, cbarii, cbarli, &
                dbarii, dbarli, ebarii, ebarli, fbarii, fbarli

    real(rk8), dimension(4), parameter :: abarl = &
         [  2.817e-2_rk8,  2.682e-2_rk8, 2.264e-2_rk8, 1.281e-2_rk8 ]
    real(rk8), dimension(4), parameter :: bbarl = &
         [  1.305e+0_rk8,  1.346e+0_rk8, 1.454e+0_rk8, 1.641e+0_rk8 ]
    real(rk8), dimension(4), parameter :: cbarl = &
         [ -5.620e-8_rk8, -6.940e-6_rk8, 4.640e-4_rk8, 0.201e+0_rk8 ]
    real(rk8), dimension(4), parameter :: dbarl = &
         [  1.630e-7_rk8,  2.350e-5_rk8, 1.240e-3_rk8, 7.560e-3_rk8 ]
    real(rk8), dimension(4), parameter :: ebarl = &
         [  0.829e+0_rk8,  0.794e+0_rk8, 0.754e+0_rk8, 0.826e+0_rk8 ]
    real(rk8), dimension(4), parameter :: fbarl = &
         [  2.482e-3_rk8,  4.226e-3_rk8, 6.560e-3_rk8, 4.353e-3_rk8 ]
!
    real(rk8), dimension(4), parameter :: abari = &
         [  3.4480e-3_rk8, 3.4480e-3_rk8, 3.4480e-3_rk8, 3.44800e-3_rk8 ]
    real(rk8), dimension(4), parameter :: bbari = &
         [  2.4310e+0_rk8, 2.4310e+0_rk8, 2.4310e+0_rk8, 2.43100e+0_rk8 ]
    real(rk8), dimension(4), parameter :: cbari = &
         [  1.0000e-5_rk8, 1.1000e-4_rk8, 1.8610e-2_rk8, 0.46658e+0_rk8 ]
    real(rk8), dimension(4), parameter :: dbari = &
         [  0.0000e+0_rk8, 1.4050e-5_rk8, 8.3280e-4_rk8, 2.05000e-5_rk8 ]
    real(rk8), dimension(4), parameter :: ebari = &
         [  0.7661e+0_rk8, 0.7730e+0_rk8, 0.7940e+0_rk8, 0.95950e+0_rk8 ]
    real(rk8), dimension(4), parameter :: fbari = &
         [  5.8510e-4_rk8, 5.6650e-4_rk8, 7.2670e-4_rk8, 1.07600e-4_rk8 ]
    !
    ! define index pointing on appropriate parameter in slingo's table
    ! for eachRRTM SW band
    !
    integer, dimension (nbndsw), parameter :: indsl = &
          [ 4, 4, 3, 3, 3, 3, 3, 2, 2, 1, 1, 1, 1, 4 ]
    ! CONVENTION : RRTMG driver takes layering from bottom to TOA.
    ! Regcm consider TOA to bottom

    ! surface pressure and scaled pressure, from which level are computed
    ! RRTM SW takes pressure in mb,hpa
    do i = ici1, ici2
      do j = jci1, jci2
        n = (j-jci1+1)+(i-ici1)*npj
        dlat(n) = m2r%xlat(j,i)
        ioro(n) = m2r%ldmsk(j,i)
      end do
    end do
    !
    ! The pressures in the RRTMG are expressed in hPa: must divide by 100.0
    !
    do i = ici1, ici2
      do j = jci1, jci2
        n = (j-jci1+1)+(i-ici1)*npj
        xptrop(n) = m2r%ptrop(j,i) * 0.01_rk8
        psfc(n) = m2r%psatms(j,i) * 0.01_rk8
      end do
    end do

    do k = 1, kz
      kj = kzp1-k
      do i = ici1, ici2
        do j = jci1, jci2
          n = (j-jci1+1)+(i-ici1)*npj
          play(n,k) = m2r%phatms(j,i,kj)*0.01_rk8
        end do
      end do
    end do
    !
    ! convert pressures from Pa to mb and define interface pressures:
    !
    do k = 1, kzp1
      kj = kzp2-k
      do i = ici1, ici2
        do j = jci1, jci2
          n = (j-jci1+1)+(i-ici1)*npj
          plev(n,k) = m2r%pfatms(j,i,kj)*0.01_rk8
        end do
      end do
    end do
    do k = kzp2, ktf
      do n = 1, npr
        plev(n,k) = stdplevf(kclimf+k-kzp2)
      end do
    end do
    ! smooth transition from top to climate at kzp1
    do k = kzp1, kth
      do n = 1, npr
        play(n,k) = stdplevh(kclimh+k-kzp1)
      end do
    end do
    !
    ! ground temperature
    !
    do i = ici1, ici2
      do j = jci1, jci2
        n = (j-jci1+1)+(i-ici1)*npj
        tsfc(n) = m2r%tg(j,i)
      end do
    end do
    !
    ! air temperatures
    !
    do k = 1, kz
      kj = kzp1-k
      do i = ici1, ici2
        do j = jci1, jci2
          n = (j-jci1+1)+(i-ici1)*npj
          tlay(n,k) = m2r%tatms(j,i,kj)
        end do
      end do
    end do
    do k = kzp1, kth
      do i = ici1, ici2
        do j = jci1, jci2
          n = (j-jci1+1)+(i-ici1)*npj
#ifdef RCEMIP
          tlay(n,k) = max(stdatm_val(real(calday,rk8),real(dayspy,rk8), &
            0.0_rk8,play(n,k),istdatm_tempk),tlay(n,kz))
#else
          tlay(n,k) = max(stdatm_val(real(calday,rk8),real(dayspy,rk8), &
            dlat(n),play(n,k),istdatm_tempk),tlay(n,kz))
#endif
        end do
      end do
    end do
    !
    ! deltaz
    !
    do k = 1, kz
      kj = kzp1 - k
      do i = ici1, ici2
        do j = jci1, jci2
          n = (j-jci1+1)+(i-ici1)*npj
          deltaz(n,k) = m2r%deltaz(j,i,kj)
        end do
      end do
    end do
    do k = kzp1, kth
      do n = 1, npr
        deltaz(n,k) = (stdhlevf(kclimf+k-kzp1) - &
                       stdhlevf(kclimf+k-kzp1-1)) * 1000.0_rk8
      end do
    end do
    !
    ! air temperature at the interface
    !
    do n = 1, npr
      tlev(n,1) = tsfc(n)
    end do
    if ( idynamic /= 3 ) then
      do k = 2, kz
        kj = kzp2-k
        do i = ici1, ici2
          do j = jci1, jci2
            n = (j-jci1+1)+(i-ici1)*npj
            p1 = (plev(n,k)/play(n,k-1))**c287
            p2 = (plev(n,k)/play(n,k))**c287
            w1 = (hsigma(kj) - sigma(kj)) / (hsigma(kj) - hsigma(kj-1))
            w2 = (sigma(kj) - hsigma(kj-1) ) / (hsigma(kj) - hsigma(kj-1))
            tlev(n,k) = w1 * tlay(n,k-1) * p1 + w2 * tlay(n,k) * p2
          end do
        end do
      end do
    else
      do k = 2, kz
        do n = 1, npr
          tlev(n,k) = 0.5_rk8 * (tlay(n,k-1) + tlay(n,k))
        end do
      end do
    end if
    do k = kzp1, ktf
      do n = 1, npr
#ifdef RCEMIP
        tlev(n,k) = max(stdatm_val(real(calday,rk8),real(dayspy,rk8), &
          0.0_rk8,plev(n,k),istdatm_tempk),tlev(n,kz))
#else
        tlev(n,k) = max(stdatm_val(real(calday,rk8),real(dayspy,rk8), &
          dlat(n),plev(n,k),istdatm_tempk),tlev(n,kz))
#endif
      end do
    end do

    if ( ipptls > 1 ) then
      do k = 1, kz
        kj = kzp1 - k
        do i = ici1, ici2
          do j = jci1, jci2
            n = (j-jci1+1)+(i-ici1)*npj
            ql1(n,k) = m2r%qxatms(j,i,kj,iqc)
            qi1(n,k) = m2r%qxatms(j,i,kj,iqi)
          end do
        end do
      end do
      do k = kzp1, kth
        do n = 1, npr
          ql1(n,k) = 0.0_rk8
          qi1(n,k) = 0.0_rk8
        end do
      end do
    end if
    !
    ! h2o volume mixing ratio
    !
    do k = 1, kz
      kj = kzp1 - k
      do i = ici1, ici2
        do j = jci1, jci2
          n = (j-jci1+1)+(i-ici1)*npj
          h2ommr(n,k) = m2r%qxatms(j,i,kj,iqv)
          h2ovmr(n,k) = h2ommr(n,k) * rep2
        end do
      end do
    end do
    do k = kzp1, kth
      do n = 1, npr
#ifdef RCEMIP
        h2ommr(n,k) = 1.0e-14_rk8
#else
        h2ommr(n,k) = stdatm_val(real(calday,rk8),real(dayspy,rk8), &
                         dlat(n),play(n,k),istdatm_qdens) / &
                      stdatm_val(real(calday,rk8),real(dayspy,rk8), &
                         dlat(n),play(n,k),istdatm_airdn) * amd/amw
#endif
        h2ovmr(n,k) = h2ommr(n,k) * rep2
      end do
    end do
    !
    ! o3 volume mixing ratio
    !
    do k = 1, kz
      kj = kzp1 - k
      do i = ici1, ici2
        do j = jci1, jci2
          n = (j-jci1+1)+(i-ici1)*npj
          o3vmr(n,k) = 0.5_rk8*(o3prof(j,i,kj)+o3prof(j,i,kj+1)) * amd/amo3
        end do
      end do
    end do
    do k = kzp1, kth
      do n = 1, npr
#ifdef RCEMIP
        o3vmr(n,k) = &
          stdatm_val(real(calday,rk8),real(dayspy,rk8), &
                     0.0_rk8,play(n,k),istdatm_ozone) / &
          stdatm_val(real(calday,rk8),real(dayspy,rk8), &
                     0.0_rk8,play(n,k),istdatm_airdn) * amd/amo3
#else
        o3vmr(n,k) = &
          stdatm_val(real(calday,rk8),real(dayspy,rk8), &
                     dlat(n),play(n,k),istdatm_ozone) / &
          stdatm_val(real(calday,rk8),real(dayspy,rk8), &
                     dlat(n),play(n,k),istdatm_airdn) * amd/amo3
#endif
      end do
    end do
    !
    ! Transform in mass mixing ratios (g/g)
    !
#ifdef RCEMIP
    do n = 1, npr
      co2vmr(n) = ghgval(igh_co2,iyear,imonth,0.0_rk8)
      co2mmr(n) = co2vmr(n)*(amco2/amd)
      ch40(n) = ghgval(igh_ch4,iyear,imonth,0.0_rk8)*(amch4/amd)
      n2o0(n) = ghgval(igh_n2o,iyear,imonth,0.0_rk8)*(amn2o/amd)
      cfc110(n) = ghgval(igh_cfc11,iyear,imonth,0.0_rk8)*(amcfc11/amd)
      cfc120(n) = ghgval(igh_cfc12,iyear,imonth,0.0_rk8)*(amcfc12/amd)
    end do
#else
    do n = 1, npr
      co2vmr(n) = ghgval(igh_co2,iyear,imonth,dlat(n))
      co2mmr(n) = co2vmr(n)*(amco2/amd)
      ch40(n) = ghgval(igh_ch4,iyear,imonth,dlat(n))*(amch4/amd)
      n2o0(n) = ghgval(igh_n2o,iyear,imonth,dlat(n))*(amn2o/amd)
      cfc110(n) = ghgval(igh_cfc11,iyear,imonth,dlat(n))*(amcfc11/amd)
      cfc120(n) = ghgval(igh_cfc12,iyear,imonth,dlat(n))*(amcfc12/amd)
    end do
#endif

    do k = 1, kth
      do n = 1, npr
        co2vmrk(n,k) = co2vmr(n)
      end do
    end do

    call trcmix(1,npr,dlat,xptrop,play,n2o0,ch40,cfc110,cfc120, &
                n2ommr,ch4mmr,cfc11mmr,cfc12mmr)

    do k = 1, kth
      do n = 1, npr
        !
        ! Form mass mixing ratios to vomlume mixing ratios
        !
        o2vmr(n,k)    = 0.209460_rk8
        n2ovmr(n,k)   = n2ommr(n,k) * (amd/amn2o)
        ch4vmr(n,k)   = ch4mmr(n,k) * (amd/amch4)
        cfc11vmr(n,k) = cfc11mmr(n,k) * (amd/amcfc11)
        cfc12vmr(n,k) = cfc12mmr(n,k) * (amd/amcfc12)
        !
        ! No data FOR NOW : IMPROVE !!
        !
        cfc22vmr(n,k) = 0.0_rk8
        ccl4vmr(n,k)  = 0.0_rk8
      end do
    end do
    !
    ! aerosols
    ! no stratospheric background for now
    !care : aermmr,rh and pint must be passed to aeroppt consitent with
    !       the regcm vertical grid i.e.top down
    !       the tauxar,tauasc,gtota calculated in aeroppt are then
    !       switched on RRTM grid.
    if ( ichem == 1 ) then
      do itr = 1, ntr
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              n = (j-jci1+1)+(i-ici1)*npj
              aermmr(n,k,itr) = max(m2r%chiatms(j,i,k,itr),0.0_rkx)
            end do
          end do
        end do
      end do
    end if
    if ( ichem == 1 .or. iclimaaer > 0 ) then
      do k = 1, kz
        do i = ici1, ici2
          do j = jci1, jci2
            n = (j-jci1+1)+(i-ici1)*npj
            rh(n,k) = m2r%rhatms(j,i,k)
          end do
        end do
      end do
      do k = 1, kzp1
        do i = ici1, ici2
          do j = jci1, jci2
            n = (j-jci1+1)+(i-ici1)*npj
            pint(n,k) = m2r%pfatms(j,i,k)
          end do
        end do
      end do
      call aeroppt(rh,pint,1,npr)
      ! adapt reverse the vertical grid for RRTM
      do i = 1, nbndsw
        do k = 1, kth
          kj = kth + 1 - k
          do n = 1, npr
            ecaer(n,k,i)  = 0.78_rk8 ! not used
            tauaer(n,k,i) = tauxar3d(n,kj,i)
            ssaaer(n,k,i) = tauasc3d(n,kj,i)
            asmaer(n,k,i) = gtota3d(n,kj,i)
          end do
        end do
      end do
      do i = 1, nbndlw
        do k = 1, kth
          kj = kth + 1 - k
          do n = 1, npr
            tauaer_lw(n,k,i) = tauxar3d_lw(n,kj,i)
          end do
        end do
      end do
    else
      ecaer(:,:,:)  = 0.78_rk8 ! not used
      tauaer(:,:,:) = 0.0_rk8
      ssaaer(:,:,:) = 0.0_rk8
      asmaer(:,:,:) = 0.0_rk8
      tauaer_lw(:,:,:) = 0.0_rk8
    end if
    !
    ! cloud fraction and cloud liquid waterpath calculation:
    ! as in STANDARD SCHEME for now (getdat) : We need to improve this
    ! according to new cloud microphysics!
    ! fractional cloud cover (dependent on relative humidity)
    !
    ! qc   = gary's mods for clouds/radiation tie-in to exmois
    !
    cldf = 0.0_rk8
    clwp = 0.0_rk8
    do k = 1, kz
      kj = kzp1 - k
      do i = ici1, ici2
        do j = jci1, jci2
          !
          ! convert liquid water content into liquid water path, i.e.
          ! multiply b deltaz
          ! deltaz,clwp are on the right grid since plev and tlay are
          ! care pressure is on bottom/toa grid
          !
          n = (j-jci1+1)+(i-ici1)*npj
          cldf(n,k) = min(m2r%cldfrc(j,i,kj),real(cftotmax,rkx))
          clwp(n,k) = m2r%cldlwc(j,i,kj) * deltaz(n,k)
          if ( clwp(n,k) > 0.0_rk8 ) then
            cldf(n,k) = min(m2r%cldfrc(j,i,kj),real(cftotmax,rkx))
          else
            cldf(n,k) = 0.0_rk8
          end if
        end do
      end do
    end do
    !
    ! fabtest:  set cloud fractional cover at top model level = 0
    ! as in std scheme
    do k = kzp1, kth
      do n = 1, npr
        cldf(n,k) = 0.0_rk8
        clwp(n,k) = 0.0_rk8
      end do
    end do
    !
    ! CLOUD Properties:
    !
    ! cloud effective radius
    !
    call cldefr_rrtm(tlay,play,rel,rei,fice)
    !
    ! partition of total water path betwwen liquide and ice.
    ! now clwp is liquide only !
    !
    do k = 1, kth
      do n = 1, npr
        ciwp(n,k) =  clwp(n,k) * fice(n,k)
        clwp(n,k) =  clwp(n,k) * (1.0_rk8 - fice(n,k))
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
    tauc  =  1.0e-10_rk8
    ssac  =  verynearone
    asmc  =  0.850_rk8
    fsfc  =  0.725_rk8

    do k = 1, nbndlw
      do i = ici1, ici2
        do j = jci1, jci2
          n = (j-jci1+1)+(i-ici1)*npj
          emis_surf(n,k) = m2r%emiss(j,i)
        end do
      end do
    end do

    outtaucl(:,:,:) = 0.0_rk8
    outtauci(:,:,:) = 0.0_rk8

    if ( inflgsw == 0 ) then
      do ns = 1, nbndsw
        !
        ! Set cloud extinction optical depth, single scatter albedo,
        ! asymmetry parameter, and forward scattered fraction:
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

        do k = 1, kz
          do n = 1, npr
            if ( clwp(n,k) < dlowval .and. ciwp(n,k) < dlowval) cycle
            ! liquid
            tmp1l = abarli + bbarli/rel(n,k)
            tmp2l = 1.0_rk8 - cbarli - dbarli*rel(n,k)
            tmp3l = fbarli*rel(n,k)
            ! ice
            tmp1i = abarii + bbarii/rei(n,k)
            tmp2i = 1.0_rk8 - cbarii - dbarii*rei(n,k)
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
                          (1.0_rk8+(1.0_rk8-0.85_rk8)*(1.0_rk8-cldf(n,k))*  &
                          clwp(n,k)*tmp1l)
              outtaucl(n,k,indsl(ns)) = outtaucl(n,k,indsl(ns)) + tauxcl(n,k)
              tauxci(n,k) = ciwp(n,k)*tmp1i*cldf(n,k) /     &
                          (1.0_rk8+(1.0_rk8-0.78_rk8)*(1.0_rk8-cldf(n,k)) * &
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
        do k = 1, kz
          do n = 1, npr
            if ( cldf(n,k) > dlowval ) then
              cldf(n,k) = 1.0_rk8
            else
              cldf(n,k) = 0.0_rk8
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
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: pmid, t
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: fice, rei, rel
    integer(ik4) :: k, n
    real(rk8) :: pnrml, weight
    ! real(rk8) :: tpara
    ! reimax - maximum ice effective radius
    real(rk8), parameter :: reimax = 30.0_rk8
    ! rirnge - range of ice radii (reimax - 10 microns)
    real(rk8), parameter :: rirnge = 20.0_rk8
    ! pirnge - nrmlzd pres range for ice particle changes
    real(rk8), parameter :: pirnge = 0.4_rk8
    ! picemn - normalized pressure below which rei=reimax
    real(rk8), parameter :: picemn = 0.4_rk8
    ! Temperatures in K (263.16, 243.16)
    real(rk8), parameter :: minus10 = wattp-10.0_rk8
    real(rk8), parameter :: minus30 = wattp-30.0_rk8

    do k = 1, kth
      do n = 1, npr
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
        ! tpara = min(1.0_rk8,max(0.0_rk8,(minus10-tm1(n,k))*0.05_rk8))
        !
        if ( ioro(n) == 1 ) then
          ! Effective liquid radius over land
          ! rel(n,k) = 6.0_rk8 + 5.0_rk8 * tpara
          rel(n,k) = 8.50_rk8
        else
          ! Effective liquid radius over ocean and sea ice
          ! rel(n,k) = 7.0_rk8 + 5.0_rk8 * tpara
          rel(n,k) = 11.0_rk8
        end if
        !
        ! Determine rei as function of normalized pressure
        !
        pnrml = pmid(n,k)/psfc(n)
        weight = max(min((pnrml-picemn)/pirnge,1.0_rk8),0.0_rk8)
        rei(n,k) = reimax - rirnge*weight
        !
        ! Define fractional amount of cloud that is ice
        !
      end do
    end do
    if ( ipptls > 1 ) then
      do k = 1, kth
        do n = 1, npr
          if ( qi1(n,k) > minqq ) then
            fice(n,k) = qi1(n,k) / (ql1(n,k)+qi1(n,k))
          else
            fice(n,k) = 0.0_rk8
          end if
        end do
      end do
    else
      do k = 1, kth
        do n = 1, npr
          if ( t(n,k) > minus10 ) then
            ! if warmer than -10 degrees C then water phase
            fice(n,k) = 0.0_rk8
          else if ( t(n,k) <= minus10 .and. t(n,k) >= minus30 ) then
            ! if colder than -10 degrees C but warmer than -30 C mixed phase
            fice(n,k) = (minus10-t(n,k))/20.0_rk8
          else
            ! if colder than -30 degrees C then ice phase
            fice(n,k) = 1.0_rk8
          end if
          ! Turn off ice radiative properties by setting fice = 0.0
          ! fice(n,k) = 0.0_rk8
        end do
      end do
    end if
  end subroutine cldefr_rrtm

end module mod_rrtmg_driver
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

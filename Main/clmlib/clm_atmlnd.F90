#include <misc.h>
#include <preproc.h>

module clm_atmlnd

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_atmlnd
!
! !DESCRIPTION:
! Handle atm2lnd, lnd2atm mapping/downscaling/upscaling/data
!
! !USES:
  use clm_varpar  , only : numrad, ndst   !ndst = number of dust bins.
                                          !only used # ifdef DUST
  use clm_varcon  , only : rair, grav, cpair
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan
  use spmdMod     , only : masterproc
  use abortutils,   only : endrun

  use clm_varpar  , only : nvoc
  use clm_drydep  , only : n_drydep
!
! !PUBLIC TYPES:
  implicit none

!----------------------------------------------------
! atmosphere -> land variables structure
!----------------------------------------------------
type atm2lnd_type
#if (defined OFFLINE)
  real(r8), pointer :: flfall(:)       !frac of liquid water in falling precip
#endif
  real(r8), pointer :: forc_t(:)       !atmospheric temperature (Kelvin)
  real(r8), pointer :: forc_u(:)       !atm wind speed, east direction (m/s)
  real(r8), pointer :: forc_v(:)       !atm wind speed, north direction (m/s)
  real(r8), pointer :: forc_wind(:)    !atmospheric wind speed   
  real(r8), pointer :: forc_q(:)       !atmospheric specific humidity (kg/kg)
  real(r8), pointer :: forc_hgt(:)     !atmospheric reference height (m)
  real(r8), pointer :: forc_hgt_u(:)   !obs height of wind [m] (new)
  real(r8), pointer :: forc_hgt_t(:)   !obs height of temperature [m] (new)
  real(r8), pointer :: forc_hgt_q(:)   !obs height of humidity [m] (new)
  real(r8), pointer :: forc_pbot(:)    !atmospheric pressure (Pa)
  real(r8), pointer :: forc_th(:)      !atm potential temperature (Kelvin)
  real(r8), pointer :: forc_vp(:)      !atmospheric vapor pressure (Pa)
  real(r8), pointer :: forc_rho(:)     !density (kg/m**3)
  real(r8), pointer :: forc_psrf(:)    !surface pressure (Pa)
  real(r8), pointer :: forc_pco2(:)    !CO2 partial pressure (Pa)
  real(r8), pointer :: forc_lwrad(:)   !downwrd IR longwave radiation (W/m**2)
  real(r8), pointer :: forc_solad(:,:) !direct beam radiation (numrad)
                                       !(vis=forc_sols , nir=forc_soll )
  real(r8), pointer :: forc_solai(:,:) !diffuse radiation (numrad)
                                       !(vis=forc_solsd, nir=forc_solld)
  real(r8), pointer :: forc_solar(:)   !incident solar radiation
  real(r8), pointer :: forc_rain(:)    !rain rate [mm/s]
  real(r8), pointer :: forc_snow(:)    !snow rate [mm/s]
  real(r8), pointer :: forc_ndep(:)    !nitrogen deposition rate (gN/m2/s)
  ! 4/14/05: PET
  ! Adding isotope code
  real(r8), pointer :: forc_pc13o2(:)  !C13O2 partial pressure (Pa)
  real(r8), pointer :: forc_po2(:)     !O2 partial pressure (Pa)
end type atm2lnd_type

!----------------------------------------------------
! land -> atmosphere variables structure
!----------------------------------------------------
type lnd2atm_type
  real(r8), pointer :: t_rad(:)        !radiative temperature (Kelvin)
  real(r8), pointer :: t_ref2m(:)      !2m surface air temperature (Kelvin)
  real(r8), pointer :: q_ref2m(:)      !2m surface specific humidity (kg/kg)
  real(r8), pointer :: h2osno(:)       !snow water (mm H2O)
  real(r8), pointer :: albd(:,:)       !(numrad) surface albedo (direct)
  real(r8), pointer :: albi(:,:)       !(numrad) surface albedo (diffuse)
  real(r8), pointer :: taux(:)         !wind stress: e-w (kg/m/s**2)
  real(r8), pointer :: tauy(:)         !wind stress: n-s (kg/m/s**2)
  real(r8), pointer :: eflx_lh_tot(:)  !total latent HF (W/m**2)  [+ to atm]
  real(r8), pointer :: eflx_sh_tot(:)  !total sensible HF (W/m**2) [+ to atm]
  real(r8), pointer :: eflx_lwrad_out(:) !IR (longwave) radiation (W/m**2)
  real(r8), pointer :: qflx_evap_tot(:)!qflx_evap(_soi + _veg) + qflx_tran_veg
  real(r8), pointer :: fsa(:)          !solar rad absorbed (total) (W/m**2)
  real(r8), pointer :: nee(:)          !net CO2 flux (kg C/m**2/s) [+ to atm]
!abt rcm below
  real(r8), pointer :: smtot(:)          !total soil moisture
  real(r8), pointer :: sm1m(:)           !soil moisture of first 1m
  real(r8), pointer :: tlef(:)           !leaf temperature [Kelvin]
  real(r8), pointer :: u10(:)            !10 meter wind
  real(r8), pointer :: tgrnd(:)          !ground temperature
  real(r8), pointer :: qflx_infl(:)      !infiltration [mm H2O/s]
  real(r8), pointer :: qflx_surf(:)      !surface runoff [mm H2O/s]
  real(r8), pointer :: qflx_drain(:)     !sub-surface runoff [mm H2O/s]
  real(r8), pointer :: uvdrag(:)         !surface drag stress
  real(r8), pointer :: sm10cm(:)         !soil moisture of top 10cm
  real(r8), pointer :: frac_sno(:)       !fraction ground covered by snow
  real(r8), pointer :: frac_veg_nosno(:) !fraction veg with no snow
  real(r8), pointer :: vdrydep(:,:)      !dry deposition velocity
#if (defined VOC)
  real(r8), pointer :: voc2rcm(:,:) !transfer to chem scheme
#endif
!abt rcm above
#if (defined DUST || defined  PROGSSLT )
  real(r8), pointer :: ram1(:)         !aerodynamical resistance (s/m)
  real(r8), pointer :: fv(:)           !friction velocity (m/s) (for dust model)
#endif
#if (defined DUST  )
  real(r8), pointer :: flxdst(:,:)       !dust flux (size bins)
#endif
end type lnd2atm_type

  type(atm2lnd_type),public,target :: atm_a2l      ! a2l fields on atm grid
  type(lnd2atm_type),public,target :: atm_l2a      ! l2a fields on atm grid

  type(atm2lnd_type),public,target :: clm_a2l      ! a2l fields on clm grid
  type(lnd2atm_type),public,target :: clm_l2a      ! l2a fields on clm grid

! !PUBLIC MEMBER FUNCTIONS:
  public :: init_atm2lnd_type
  public :: init_lnd2atm_type
  public :: clm_mapa2l
  public :: clm_mapl2a
  public :: clm_map2gcell
!abt rcm below
  public :: clm2rcm
!abt rcm above
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein and Tony Craig, 2006-01-10
! 2008.5.8    A Tawfik Revised to work with RegCM
!
!
! !PRIVATE MEMBER FUNCTIONS:

!EOP
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_atm2lnd_type
!
! !INTERFACE:
  subroutine init_atm2lnd_type(beg, end, a2l)
!
! !DESCRIPTION:
! Initialize atmospheric variables required by the land
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (atm2lnd_type), intent(inout):: a2l
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Modified by T Craig, 11/01/05 for finemesh project
!
!EOP
!
! !LOCAL VARIABLES:
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

#if (defined OFFLINE)
  allocate(a2l%flfall(beg:end))
#endif
  allocate(a2l%forc_t(beg:end))
  allocate(a2l%forc_u(beg:end))
  allocate(a2l%forc_v(beg:end))
  allocate(a2l%forc_wind(beg:end))
  allocate(a2l%forc_q(beg:end))
  allocate(a2l%forc_hgt(beg:end))
  allocate(a2l%forc_hgt_u(beg:end))
  allocate(a2l%forc_hgt_t(beg:end))
  allocate(a2l%forc_hgt_q(beg:end))
  allocate(a2l%forc_pbot(beg:end))
  allocate(a2l%forc_th(beg:end))
  allocate(a2l%forc_vp(beg:end))
  allocate(a2l%forc_rho(beg:end))
  allocate(a2l%forc_psrf(beg:end))
  allocate(a2l%forc_pco2(beg:end))
  allocate(a2l%forc_lwrad(beg:end))
  allocate(a2l%forc_solad(beg:end,numrad))
  allocate(a2l%forc_solai(beg:end,numrad))
  allocate(a2l%forc_solar(beg:end))
  allocate(a2l%forc_rain(beg:end))
  allocate(a2l%forc_snow(beg:end))
  allocate(a2l%forc_ndep(beg:end))
  ! 4/14/05: PET
  ! Adding isotope code
  allocate(a2l%forc_pc13o2(beg:end))
  allocate(a2l%forc_po2(beg:end))

! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

#if (defined OFFLINE)
  a2l%flfall(beg:end) = ival
#endif
  a2l%forc_t(beg:end) = ival
  a2l%forc_u(beg:end) = ival
  a2l%forc_v(beg:end) = ival
  a2l%forc_wind(beg:end) = ival
  a2l%forc_q(beg:end) = ival
  a2l%forc_hgt(beg:end) = ival
  a2l%forc_hgt_u(beg:end) = ival
  a2l%forc_hgt_t(beg:end) = ival
  a2l%forc_hgt_q(beg:end) = ival
  a2l%forc_pbot(beg:end) = ival
  a2l%forc_th(beg:end) = ival
  a2l%forc_vp(beg:end) = ival
  a2l%forc_rho(beg:end) = ival
  a2l%forc_psrf(beg:end) = ival
  a2l%forc_pco2(beg:end) = ival
  a2l%forc_lwrad(beg:end) = ival
  a2l%forc_solad(beg:end,1:numrad) = ival
  a2l%forc_solai(beg:end,1:numrad) = ival
  a2l%forc_solar(beg:end) = ival
  a2l%forc_rain(beg:end) = ival
  a2l%forc_snow(beg:end) = ival
  a2l%forc_ndep(beg:end) = ival
  ! 4/14/05: PET
  ! Adding isotope code
  a2l%forc_pc13o2(beg:end) = ival
  a2l%forc_po2(beg:end) = ival

end subroutine init_atm2lnd_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2atm_type
!
! !INTERFACE:
  subroutine init_lnd2atm_type(beg, end, l2a)
!
! !DESCRIPTION:
! Initialize land variables required by the atmosphere
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (lnd2atm_type), intent(inout):: l2a
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Modified by T Craig, 11/01/05 for finemesh project
!
!EOP
!
! !LOCAL VARIABLES:
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(l2a%t_rad(beg:end))
  allocate(l2a%t_ref2m(beg:end))
  allocate(l2a%q_ref2m(beg:end))
  allocate(l2a%h2osno(beg:end))
  allocate(l2a%albd(beg:end,1:numrad))
  allocate(l2a%albi(beg:end,1:numrad))
  allocate(l2a%taux(beg:end))
  allocate(l2a%tauy(beg:end))
  allocate(l2a%eflx_lwrad_out(beg:end))
  allocate(l2a%eflx_sh_tot(beg:end))
  allocate(l2a%eflx_lh_tot(beg:end))
  allocate(l2a%qflx_evap_tot(beg:end))
  allocate(l2a%fsa(beg:end))
  allocate(l2a%nee(beg:end))
!abt rcm below
  allocate(l2a%smtot(beg:end))
  allocate(l2a%sm1m(beg:end))
  allocate(l2a%u10(beg:end))
  allocate(l2a%tlef(beg:end))
  allocate(l2a%sm10cm(beg:end))
  allocate(l2a%uvdrag(beg:end))
  allocate(l2a%qflx_infl(beg:end))
  allocate(l2a%qflx_surf(beg:end))
  allocate(l2a%qflx_drain(beg:end))
  allocate(l2a%tgrnd(beg:end))
  allocate(l2a%frac_sno(beg:end))
  allocate(l2a%frac_veg_nosno(beg:end))
#if (defined VOC)
  allocate(l2a%voc2rcm(beg:end,nvoc))
#endif
  allocate(l2a%vdrydep(beg:end,n_drydep))
!abt rcm above


#if (defined DUST || defined  PROGSSLT )
  allocate(l2a%ram1(beg:end))
  allocate(l2a%fv(beg:end))
#endif
#if (defined DUST )
  allocate(l2a%flxdst(beg:end,1:ndst))
#endif

! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  l2a%t_rad(beg:end) = ival
  l2a%t_ref2m(beg:end) = ival
  l2a%q_ref2m(beg:end) = ival
  l2a%h2osno(beg:end) = ival
  l2a%albd(beg:end,1:numrad) = ival
  l2a%albi(beg:end,1:numrad) = ival
  l2a%taux(beg:end) = ival
  l2a%tauy(beg:end) = ival
  l2a%eflx_lwrad_out(beg:end) = ival
  l2a%eflx_sh_tot(beg:end) = ival
  l2a%eflx_lh_tot(beg:end) = ival
  l2a%qflx_evap_tot(beg:end) = ival
  l2a%fsa(beg:end) = ival
  l2a%nee(beg:end) = ival
!abt rcm below
  l2a%smtot(beg:end)      = ival
  l2a%sm1m(beg:end)       = ival
  l2a%u10(beg:end)        = ival
  l2a%tlef(beg:end)       = ival
  l2a%sm10cm(beg:end)     = ival
  l2a%uvdrag(beg:end)     = ival
  l2a%qflx_infl(beg:end)  = ival
  l2a%qflx_surf(beg:end)  = ival
  l2a%qflx_drain(beg:end) = ival
  l2a%tgrnd(beg:end)      = ival
  l2a%frac_sno(beg:end)   = ival
  l2a%frac_veg_nosno(beg:end)      = ival
#if (defined VOC)
  l2a%voc2rcm(beg:end,1:nvoc) = ival
#endif
  l2a%vdrydep(beg:end,:) = ival
!abt rcm above
#if (defined DUST || defined  PROGSSLT )
  l2a%ram1(beg:end) = ival
  l2a%fv(beg:end) = ival
#endif
#if (defined DUST )
  l2a%flxdst(beg:end,1:ndst) = ival
#endif
end subroutine init_lnd2atm_type

! abt below rcm subrountine  
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_mapr2l
!
! !INTERFACE:
  subroutine clm_mapr2l(a2l_src, a2l_dst)
!
! !DESCRIPTION:
! Maps atm2lnd fields from external grid to clm grid
!
! !USES:
  use decompMod, only : get_proc_bounds, get_proc_bounds_atm
  use areaMod  , only : map_maparrayl, map1dl_a2l, map1dl_l2a, map_setptrs
  use decompMod, only : ldecomp,adecomp
  use domainMod, only : ldomain,adomain
  use QSatMod,   only : QSat
!
! !ARGUMENTS:
  implicit none
  type(atm2lnd_type), intent(in)  :: a2l_src
  type(atm2lnd_type), intent(out) :: a2l_dst
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n                     ! loop counter
  integer :: ix                    ! field index
  integer :: nflds                 ! number of fields to be mapped
  integer :: nradflds              ! size of 2nd dim in arrays
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  real(r8),pointer :: asrc(:,:)    ! temporary source data
  real(r8),pointer :: adst(:,:)    ! temporary dest data
  integer          :: nmap         ! size of map
  integer          :: mo           ! size of map
  integer, pointer :: src(:)       ! map src index
  integer, pointer :: dst(:)       ! map dst index
  real(r8),pointer :: wts(:)       ! map wts values
  integer :: ns                    !source (atm) indexes
  integer :: nd                    !destination (lnd) indexes
  ! temporaries for topo downscaling:
  real(r8):: hsurf_a,hsurf_l,Hbot,Hsrf,lapse
  real(r8):: zbot_a, tbot_a, pbot_a, thbot_a, qbot_a, qs_a, es_a
  real(r8):: zbot_l, tbot_l, pbot_l, thbot_l, qbot_l, qs_l, es_l
  real(r8):: tsrf_l, psrf_l, egcm_l, rhos_l
  real(r8):: dum1,dum2,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8
  real(r8),allocatable :: qsum(:)
  logical :: first_call = .true.
!------------------------------------------------------------------------------

  if (first_call .and. masterproc) then
    write(6,*) 'clm_mapa2l subroutine'
  endif

  nradflds = size(a2l_src%forc_solad,dim=2)
  if (nradflds /= numrad) then
    write(6,*) 'clm_mapa2l ERROR: nradflds ne numrad ',nradflds,numrad
    call endrun()
  endif

  !--- allocate temporaries
  call get_proc_bounds_atm(begg_s, endg_s)
  call get_proc_bounds    (begg_d, endg_d)

!  nflds = 21+2*numrad

!  allocate(asrc(begg_s:endg_s,nflds))
!  allocate(adst(begg_d:endg_d,nflds))

!  ix = 0
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_t(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_u(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_v(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_wind(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_q(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt_u(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt_t(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt_q(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_pbot(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_th(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_vp(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_rho(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_psrf(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_pco2(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_lwrad(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_solar(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_rain(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_snow(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_pc13o2(:)  
!  ix=ix+1; asrc(:,ix) = a2l_src%forc_po2(:)  
!  do n = 1,numrad
!     ix=ix+1; asrc(:,ix) = a2l_src%forc_solad(:,n)  
!     ix=ix+1; asrc(:,ix) = a2l_src%forc_solai(:,n)  
!  enddo
!-forc_ndep is not recd from atm,don't know why it's in a2l (TCFIX) ---
!-forc_ndep cannot be updated here, array will be trashed and CN will fail ---
!  asrc(:,xx) = a2l_src%forc_ndep(:)  

#if (defined OFFLINE)
  call map_maparrayl(begg_s, endg_s, begg_d, endg_d, 1, a2l_src%flfall, a2l_dst%flfall, map1dl_a2l)
#endif
!  call map_maparrayl(begg_s, endg_s, begg_d, endg_d, nflds, asrc, adst, map1dl_a2l)

  a2l_dst%forc_t(:)     =     a2l_src%forc_t(:)
  a2l_dst%forc_u(:)     =     a2l_src%forc_u(:)
  a2l_dst%forc_v(:)     =     a2l_src%forc_v(:)
  a2l_dst%forc_wind(:)  =     a2l_src%forc_wind(:)
  a2l_dst%forc_q(:)     =     a2l_src%forc_q(:)
  a2l_dst%forc_hgt(:)   =     a2l_src%forc_hgt(:)
  a2l_dst%forc_hgt_u(:) =     a2l_src%forc_hgt_u(:)
  a2l_dst%forc_hgt_t(:) =     a2l_src%forc_hgt_t(:)
  a2l_dst%forc_hgt_q(:) =     a2l_src%forc_hgt_q(:)
  a2l_dst%forc_pbot(:)  =   a2l_src%forc_pbot(:)
  a2l_dst%forc_th(:)    =   a2l_src%forc_th(:)
  a2l_dst%forc_vp(:)    =   a2l_src%forc_vp(:)
  a2l_dst%forc_rho(:)   =   a2l_src%forc_rho(:)
  a2l_dst%forc_psrf(:)  =   a2l_src%forc_psrf(:)
  a2l_dst%forc_pco2(:)  =   a2l_src%forc_pco2(:)
  a2l_dst%forc_lwrad(:) =   a2l_src%forc_lwrad(:)
  a2l_dst%forc_solar(:) =   a2l_src%forc_solar(:)
  a2l_dst%forc_rain(:)  =   a2l_src%forc_rain(:)
  a2l_dst%forc_snow(:)  =   a2l_src%forc_snow(:)
  a2l_dst%forc_pc13o2(:)=   a2l_src%forc_pc13o2(:)
  a2l_dst%forc_po2(:)   =   a2l_src%forc_po2(:)
  do n = 1,numrad
     a2l_dst%forc_solad(:,n)  = a2l_src%forc_solad(:,n)
     a2l_dst%forc_solai(:,n)  = a2l_src%forc_solai(:,n)
  enddo

!  deallocate(asrc)
!  deallocate(adst)

  if (first_call.and.masterproc) then
    write(6,*) 'clm_mapr2l mapping complete'
  endif

!-topographic downscaling
!-only call this if there is more than 1 land cell / atm cell somewhere
  call map_setptrs(map1dl_l2a,dstmo=mo)
  if (mo > 1) then

  if (first_call.and.masterproc) then
    write(6,*) 'clm_mapr2l downscaling ON'
  endif

  call map_setptrs(map1dl_a2l,nwts=nmap,src=src,dst=dst,dstmo=mo)
  if (mo /= 1) then
     write(6,*)' clm_mapr2l ERROR: map1dl_a2l mo not 1 ',mo
     call endrun()
  endif

  lapse   = 0.0065_r8                  ! hardwired in multiple places in cam

  do n = 1,nmap
    ns = src(n)
    nd = dst(n)

    hsurf_a = adomain%topo(ns)        ! atm elevation
    hsurf_l = ldomain%ntop(nd)        ! lnd elevation

    if (abs(hsurf_a - hsurf_l) .gt. 0.1_r8) then

       tbot_a = a2l_src%forc_t(ns)        ! atm temp at bot
       thbot_a= a2l_src%forc_th(ns)       ! atm pot temp at bot
       pbot_a = a2l_src%forc_pbot(ns)     ! atm press at bot
       qbot_a = a2l_src%forc_q(ns)        ! atm sp humidity at bot
       zbot_a = a2l_src%forc_hgt(ns)      ! atm ref height

       zbot_l = zbot_a
       tbot_l = tbot_a-lapse*(hsurf_l-hsurf_a)          ! lnd temp for topo

       Hbot   = rair*0.5_r8*(tbot_a+tbot_l)/grav        ! scale ht at avg temp
       pbot_l = pbot_a*exp(-(hsurf_l-hsurf_a)/Hbot)     ! lnd press for topo
       thbot_l= tbot_l*exp((zbot_l/Hbot)*(rair/cpair))  ! pot temp calc

       tsrf_l = tbot_l-lapse*(-zbot_l)                  ! lnd temp at surface
       Hsrf   = rair*0.5_r8*(tbot_l+tsrf_l)/grav        ! scale ht at avg temp
       psrf_l = pbot_l*exp(-(zbot_l)/Hsrf)              ! lnd press for topo

       call Qsat(tbot_a,pbot_a,es_a,dum1,qs_a,dum2)
       call Qsat(tbot_l,pbot_l,es_l,dum1,qs_l,dum2)
       qbot_l = qbot_a*(qs_l/qs_a)

       a2l_dst%forc_hgt(nd)  = zbot_l
       a2l_dst%forc_t(nd)    = tbot_l
       a2l_dst%forc_pbot(nd) = pbot_l
       a2l_dst%forc_th(nd)   = thbot_l
       a2l_dst%forc_q(nd)    = qbot_l
       a2l_dst%forc_vp(nd)   = es_l
       a2l_dst%forc_psrf(nd) = psrf_l

    endif
  enddo

  allocate(qsum(begg_s:endg_s))
  qsum = 0.0_r8
  call map_setptrs(map1dl_l2a,nwts=nmap,src=src,dst=dst,wts=wts)
  do n = 1,nmap
    ns = dst(n)
    nd = src(n)
    qsum(ns) = qsum(ns) + wts(n)* a2l_dst%forc_q(nd)
  enddo

!abt below call map_setptrs(map1dl_a2l,nwts=nmap,src=src,dst=dst)
!  do n = 1,nmap
!    ns = src(n)
!    nd = dst(n)

!    qbot_a = a2l_src%forc_q(ns)        ! atm specific humidity
!    qbot_l = a2l_dst%forc_q(nd)        ! lnd specific humidity
!    pbot_l = a2l_dst%forc_pbot(nd)  
!    tbot_l = a2l_dst%forc_t(nd)  

!    qbot_l = qbot_l - (qsum(ns) - qbot_a)        ! normalize
!    egcm_l = qbot_l*pbot_l/(0.622+0.378*qbot_l)
!    rhos_l = (pbot_l-0.378*egcm_l) / (rair*tbot_l)

!    a2l_dst%forc_q(nd)    = qbot_l
!    a2l_dst%forc_rho(nd)  = rhos_l

!abt above  enddo

  deallocate(qsum)

! --- check ---
  call map_setptrs(map1dl_l2a,nwts=nmap,src=src,dst=dst,wts=wts)
  do ns = begg_s,endg_s
    sum1 = 0.0_r8
    sum2 = 0.0_r8
    sum3 = 0.0_r8
    sum4 = 0.0_r8
    sum5 = 0.0_r8
    sum6 = 0.0_r8
    do n = 1,nmap
      if (dst(n) == ns) then
        nd = src(n)
        sum1 = sum1 + ldomain%ntop(nd)   * wts(n)
        sum2 = sum2 + a2l_dst%forc_t(nd)    * wts(n)
        sum3 = sum3 + a2l_dst%forc_q(nd)    * wts(n)
        sum4 = sum4 + a2l_dst%forc_hgt(nd)  * wts(n)
        sum5 = sum5 + a2l_dst%forc_pbot(nd) * wts(n)
        sum6 = sum6 + a2l_dst%forc_th(nd)   * wts(n)
      endif
    enddo
    if   ((abs(sum1 - adomain%topo(ns))   > 1.0e-8) &
      .or.(abs(sum2 - a2l_src%forc_t(ns))    > 1.0e-3) &
      .or.(abs(sum3 - a2l_src%forc_q(ns))    > 1.0e-8) &
      .or.(abs(sum4 - a2l_src%forc_hgt(ns))  > 1.0e-6) &
!      .or.(abs(sum5 - a2l_src%forc_pbot(ns)) > 1.0e-6) &
!      .or.(abs(sum6 - a2l_src%forc_th(ns))   > 1.0e-6) &
       ) then
      write(6,*) 'clm_map2l check ERROR topo ',sum1,adomain%topo(ns)
      write(6,*) 'clm_map2l check ERROR t    ',sum2,a2l_src%forc_t(ns)
      write(6,*) 'clm_map2l check ERROR q    ',sum3,a2l_src%forc_q(ns)
      write(6,*) 'clm_map2l check ERROR hgt  ',sum4,a2l_src%forc_hgt(ns)
      write(6,*) 'clm_map2l check ERROR pbot ',sum5,a2l_src%forc_pbot(ns)
      write(6,*) 'clm_map2l check ERROR th   ',sum6,a2l_src%forc_th(ns)
!      call endrun()
    endif
  enddo

  endif   ! mx_ovr > 1

  first_call = .false.

end subroutine clm_mapr2l




!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_mapa2l
!
! !INTERFACE:
  subroutine clm_mapa2l(a2l_src, a2l_dst)
!
! !DESCRIPTION:
! Maps atm2lnd fields from external grid to clm grid
!
! !USES:
  use decompMod, only : get_proc_bounds, get_proc_bounds_atm
  use areaMod  , only : map_maparrayl, map1dl_a2l, map1dl_l2a, map_setptrs
  use decompMod, only : ldecomp,adecomp
  use domainMod, only : ldomain,adomain
  use QSatMod,   only : QSat
!
! !ARGUMENTS:
  implicit none
  type(atm2lnd_type), intent(in)  :: a2l_src
  type(atm2lnd_type), intent(out) :: a2l_dst
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n                     ! loop counter
  integer :: ix                    ! field index
  integer :: nflds                 ! number of fields to be mapped
  integer :: nradflds              ! size of 2nd dim in arrays
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  real(r8),pointer :: asrc(:,:)    ! temporary source data
  real(r8),pointer :: adst(:,:)    ! temporary dest data
  integer          :: nmap         ! size of map
  integer          :: mo           ! size of map
  integer, pointer :: src(:)       ! map src index
  integer, pointer :: dst(:)       ! map dst index
  real(r8),pointer :: wts(:)       ! map wts values
  integer :: ns                    !source (atm) indexes
  integer :: nd                    !destination (lnd) indexes
  ! temporaries for topo downscaling:
  real(r8):: hsurf_a,hsurf_l,Hbot,Hsrf,lapse
  real(r8):: zbot_a, tbot_a, pbot_a, thbot_a, qbot_a, qs_a, es_a
  real(r8):: zbot_l, tbot_l, pbot_l, thbot_l, qbot_l, qs_l, es_l
  real(r8):: tsrf_l, psrf_l, egcm_l, rhos_l
  real(r8):: dum1,dum2,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8
  real(r8),allocatable :: qsum(:)
  logical :: first_call = .true.
!------------------------------------------------------------------------------

  if (first_call .and. masterproc) then
    write(6,*) 'clm_mapa2l subroutine'
  endif

  nradflds = size(a2l_src%forc_solad,dim=2)
  if (nradflds /= numrad) then
    write(6,*) 'clm_mapa2l ERROR: nradflds ne numrad ',nradflds,numrad
    call endrun()
  endif

  !--- allocate temporaries
  call get_proc_bounds_atm(begg_s, endg_s)
  call get_proc_bounds    (begg_d, endg_d)

  nflds = 21+2*numrad

  allocate(asrc(begg_s:endg_s,nflds))
  allocate(adst(begg_d:endg_d,nflds))

  ix = 0
  ix=ix+1; asrc(:,ix) = a2l_src%forc_t(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_u(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_v(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_wind(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_q(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt_u(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt_t(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt_q(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_pbot(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_th(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_vp(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_rho(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_psrf(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_pco2(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_lwrad(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_solar(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_rain(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_snow(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_pc13o2(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_po2(:)  
  do n = 1,numrad
     ix=ix+1; asrc(:,ix) = a2l_src%forc_solad(:,n)  
     ix=ix+1; asrc(:,ix) = a2l_src%forc_solai(:,n)  
  enddo
!-forc_ndep is not recd from atm,don't know why it's in a2l (TCFIX) ---
!-forc_ndep cannot be updated here, array will be trashed and CN will fail ---
!  asrc(:,xx) = a2l_src%forc_ndep(:)  

#if (defined OFFLINE)
  call map_maparrayl(begg_s, endg_s, begg_d, endg_d, 1, a2l_src%flfall, a2l_dst%flfall, map1dl_a2l)
#endif
  call map_maparrayl(begg_s, endg_s, begg_d, endg_d, nflds, asrc, adst, map1dl_a2l)

  ix = 0
  ix=ix+1; a2l_dst%forc_t(:)     =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_u(:)     =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_v(:)     =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_wind(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_q(:)     =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_hgt(:)   =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_hgt_u(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_hgt_t(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_hgt_q(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_pbot(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_th(:)    =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_vp(:)    =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_rho(:)   =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_psrf(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_pco2(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_lwrad(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_solar(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_rain(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_snow(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_pc13o2(:)=   adst(:,ix)
  ix=ix+1; a2l_dst%forc_po2(:)   =   adst(:,ix)
  do n = 1,numrad
     ix=ix+1; a2l_dst%forc_solad(:,n)  = adst(:,ix)
     ix=ix+1; a2l_dst%forc_solai(:,n)  = adst(:,ix)
  enddo

  deallocate(asrc)
  deallocate(adst)

  if (first_call.and.masterproc) then
    write(6,*) 'clm_mapa2l mapping complete'
  endif

!-topographic downscaling
!-only call this if there is more than 1 land cell / atm cell somewhere
  call map_setptrs(map1dl_l2a,dstmo=mo)
  if (mo > 1) then

  if (first_call.and.masterproc) then
    write(6,*) 'clm_mapa2l downscaling ON'
  endif

  call map_setptrs(map1dl_a2l,nwts=nmap,src=src,dst=dst,dstmo=mo)
  if (mo /= 1) then
     write(6,*)' clm_mapa2l ERROR: map1dl_a2l mo not 1 ',mo
     call endrun()
  endif

  lapse   = 0.0065_r8                  ! hardwired in multiple places in cam

  do n = 1,nmap
    ns = src(n)
    nd = dst(n)

    hsurf_a = adomain%topo(ns)        ! atm elevation
    hsurf_l = ldomain%ntop(nd)        ! lnd elevation

    if (abs(hsurf_a - hsurf_l) .gt. 0.1_r8) then

       tbot_a = a2l_src%forc_t(ns)        ! atm temp at bot
       thbot_a= a2l_src%forc_th(ns)       ! atm pot temp at bot
       pbot_a = a2l_src%forc_pbot(ns)     ! atm press at bot
       qbot_a = a2l_src%forc_q(ns)        ! atm sp humidity at bot
       zbot_a = a2l_src%forc_hgt(ns)      ! atm ref height

       zbot_l = zbot_a
       tbot_l = tbot_a-lapse*(hsurf_l-hsurf_a)          ! lnd temp for topo

       Hbot   = rair*0.5_r8*(tbot_a+tbot_l)/grav        ! scale ht at avg temp
       pbot_l = pbot_a*exp(-(hsurf_l-hsurf_a)/Hbot)     ! lnd press for topo
       thbot_l= tbot_l*exp((zbot_l/Hbot)*(rair/cpair))  ! pot temp calc

       tsrf_l = tbot_l-lapse*(-zbot_l)                  ! lnd temp at surface
       Hsrf   = rair*0.5_r8*(tbot_l+tsrf_l)/grav        ! scale ht at avg temp
       psrf_l = pbot_l*exp(-(zbot_l)/Hsrf)              ! lnd press for topo

       call Qsat(tbot_a,pbot_a,es_a,dum1,qs_a,dum2)
       call Qsat(tbot_l,pbot_l,es_l,dum1,qs_l,dum2)
       qbot_l = qbot_a*(qs_l/qs_a)

       a2l_dst%forc_hgt(nd)  = zbot_l
       a2l_dst%forc_t(nd)    = tbot_l
       a2l_dst%forc_pbot(nd) = pbot_l
       a2l_dst%forc_th(nd)   = thbot_l
       a2l_dst%forc_q(nd)    = qbot_l
       a2l_dst%forc_vp(nd)   = es_l
       a2l_dst%forc_psrf(nd) = psrf_l

    endif
  enddo

  allocate(qsum(begg_s:endg_s))
  qsum = 0.0_r8
  call map_setptrs(map1dl_l2a,nwts=nmap,src=src,dst=dst,wts=wts)
  do n = 1,nmap
    ns = dst(n)
    nd = src(n)
    qsum(ns) = qsum(ns) + wts(n)* a2l_dst%forc_q(nd)
  enddo

  call map_setptrs(map1dl_a2l,nwts=nmap,src=src,dst=dst)
  do n = 1,nmap
    ns = src(n)
    nd = dst(n)

    qbot_a = a2l_src%forc_q(ns)        ! atm specific humidity
    qbot_l = a2l_dst%forc_q(nd)        ! lnd specific humidity
    pbot_l = a2l_dst%forc_pbot(nd)  
    tbot_l = a2l_dst%forc_t(nd)  

    qbot_l = qbot_l - (qsum(ns) - qbot_a)        ! normalize
    egcm_l = qbot_l*pbot_l/(0.622+0.378*qbot_l)
    rhos_l = (pbot_l-0.378*egcm_l) / (rair*tbot_l)

    a2l_dst%forc_q(nd)    = qbot_l
    a2l_dst%forc_rho(nd)  = rhos_l

  enddo

  deallocate(qsum)

! --- check ---
  call map_setptrs(map1dl_l2a,nwts=nmap,src=src,dst=dst,wts=wts)
  do ns = begg_s,endg_s
    sum1 = 0.0_r8
    sum2 = 0.0_r8
    sum3 = 0.0_r8
    sum4 = 0.0_r8
    sum5 = 0.0_r8
    sum6 = 0.0_r8
    do n = 1,nmap
      if (dst(n) == ns) then
        nd = src(n)
        sum1 = sum1 + ldomain%ntop(nd)   * wts(n)
        sum2 = sum2 + a2l_dst%forc_t(nd)    * wts(n)
        sum3 = sum3 + a2l_dst%forc_q(nd)    * wts(n)
        sum4 = sum4 + a2l_dst%forc_hgt(nd)  * wts(n)
        sum5 = sum5 + a2l_dst%forc_pbot(nd) * wts(n)
        sum6 = sum6 + a2l_dst%forc_th(nd)   * wts(n)
      endif
    enddo
    if   ((abs(sum1 - adomain%topo(ns))   > 1.0e-8) &
      .or.(abs(sum2 - a2l_src%forc_t(ns))    > 1.0e-3) &
      .or.(abs(sum3 - a2l_src%forc_q(ns))    > 1.0e-8) &
      .or.(abs(sum4 - a2l_src%forc_hgt(ns))  > 1.0e-6) &
!      .or.(abs(sum5 - a2l_src%forc_pbot(ns)) > 1.0e-6) &
!      .or.(abs(sum6 - a2l_src%forc_th(ns))   > 1.0e-6) &
       ) then
      write(6,*) 'clm_map2l check ERROR topo ',sum1,adomain%topo(ns)
      write(6,*) 'clm_map2l check ERROR t    ',sum2,a2l_src%forc_t(ns)
      write(6,*) 'clm_map2l check ERROR q    ',sum3,a2l_src%forc_q(ns)
      write(6,*) 'clm_map2l check ERROR hgt  ',sum4,a2l_src%forc_hgt(ns)
      write(6,*) 'clm_map2l check ERROR pbot ',sum5,a2l_src%forc_pbot(ns)
      write(6,*) 'clm_map2l check ERROR th   ',sum6,a2l_src%forc_th(ns)
!      call endrun()
    endif
  enddo

  endif   ! mx_ovr > 1

  first_call = .false.

end subroutine clm_mapa2l

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_mapl2a
!
! !INTERFACE:
  subroutine clm_mapl2a(l2a_src, l2a_dst)
!
! !DESCRIPTION:
! Maps lnd2atm fields from clm grid to external grid
!
! !USES:
  use decompMod, only : get_proc_bounds, get_proc_bounds_atm
  use areaMod  , only : map_maparrayl, map1dl_l2a
!
! !ARGUMENTS:
  implicit none
  type(lnd2atm_type), intent(in)  :: l2a_src
  type(lnd2atm_type), intent(out) :: l2a_dst
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n                     ! loop counter
  integer :: ix                    ! field index
  integer :: nflds                 ! number of fields to be mapped
  integer :: nradflds              ! size of 2nd dim in arrays
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  real(r8),pointer :: asrc(:,:)    ! temporary source data
  real(r8),pointer :: adst(:,:)    ! temporary dest data
#if (defined DUST )
  integer :: m                     ! loop counter
#endif
!------------------------------------------------------------------------------

  nradflds = size(l2a_src%albd,dim=2)
  if (nradflds /= numrad) then
    write(6,*) 'clm_mapl2a ERROR: nradflds ne numrad ',nradflds,numrad
    call endrun()
  endif

  !--- allocate temporaries
  call get_proc_bounds    (begg_s, endg_s)
  call get_proc_bounds_atm(begg_d, endg_d)

  nflds = 12+2*numrad

#if (defined DUST || defined  PROGSSLT )  
  ! add on fv, ram1 
  nflds = nflds + 2 
#endif

#if (defined DUST  ) 
  ! add on the number of dust bins (for flxdust)
  nflds = nflds + ndst
#endif


  allocate(asrc(begg_s:endg_s,nflds))
  allocate(adst(begg_d:endg_d,nflds))

  ix = 0
  ix=ix+1; asrc(:,ix) = l2a_src%t_rad(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%t_ref2m(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%q_ref2m(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%h2osno(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%taux(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%tauy(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%eflx_lh_tot(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%eflx_sh_tot(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%eflx_lwrad_out(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%qflx_evap_tot(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%fsa(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%nee(:)  
  do n = 1,numrad
     ix=ix+1; asrc(:,ix) = l2a_src%albd(:,n)  
     ix=ix+1; asrc(:,ix) = l2a_src%albi(:,n)  
  enddo
#if (defined DUST || defined  PROGSSLT )
  ix=ix+1; asrc(:,ix) = l2a_src%ram1(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%fv(:)
#endif
#if (defined DUST )
  do m = 1,ndst  ! dust bins
     ix=ix+1; asrc(:,ix) = l2a_src%flxdst(:,m)  
  end do !m
#endif

  call map_maparrayl(begg_s, endg_s, begg_d, endg_d, nflds, asrc, adst, map1dl_l2a)

  ix = 0
  ix=ix+1; l2a_dst%t_rad(:)          = adst(:,ix)
  ix=ix+1; l2a_dst%t_ref2m(:)        = adst(:,ix)
  ix=ix+1; l2a_dst%q_ref2m(:)        = adst(:,ix)
  ix=ix+1; l2a_dst%h2osno(:)         = adst(:,ix)
  ix=ix+1; l2a_dst%taux(:)           = adst(:,ix)
  ix=ix+1; l2a_dst%tauy(:)           = adst(:,ix)
  ix=ix+1; l2a_dst%eflx_lh_tot(:)    = adst(:,ix)
  ix=ix+1; l2a_dst%eflx_sh_tot(:)    = adst(:,ix)
  ix=ix+1; l2a_dst%eflx_lwrad_out(:) = adst(:,ix)
  ix=ix+1; l2a_dst%qflx_evap_tot(:)  = adst(:,ix)
  ix=ix+1; l2a_dst%fsa(:)            = adst(:,ix)
  ix=ix+1; l2a_dst%nee(:)            = adst(:,ix)
  do n = 1,numrad
     ix=ix+1; l2a_dst%albd(:,n)      = adst(:,ix)
     ix=ix+1; l2a_dst%albi(:,n)      = adst(:,ix)
  enddo

#if (defined DUST || defined  PROGSSLT )
  ix=ix+1; l2a_dst%ram1(:)           = adst(:,ix)
  ix=ix+1; l2a_dst%fv(:)             = adst(:,ix)
#endif
#if (defined DUST  )
  do m = 1,ndst  ! dust bins
     ix=ix+1; l2a_dst%flxdst(:,m)    = adst(:,ix)
  end do !m
#endif

  deallocate(asrc)
  deallocate(adst)

end subroutine clm_mapl2a

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_map2gcell
!
! !INTERFACE: subroutine clm_map2gcell(init)
  subroutine clm_map2gcell(init)
!
! !DESCRIPTION:
! Compute l2a component of gridcell derived type
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use subgridAveMod
  use decompMod   , only : get_proc_bounds
  use clm_varcon  , only : sb
  use clm_varpar  , only : numrad
!abt added below
  use clm_varcon  , only : denh2o
  use clm_varpar  , only : nlevsoi
  use mod_dynparam
  use clm_varsur  , only : cgaschem
!abt added above
!
! !ARGUMENTS:
  implicit none
  save
  logical, optional, intent(in) :: init  ! if true=>only set a subset of arguments
!
! !REVISION HISTORY:
! Mariana Vertenstein: created 03/10-25
! 2008.5.8    A Tawfik Revised to work with RegCM
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: begp, endp      ! per-proc beginning and ending pft indices
  integer :: begc, endc      ! per-proc beginning and ending column indices
  integer :: begl, endl      ! per-proc beginning and ending landunit indices
  integer :: begg, endg      ! per-proc gridcell ending gridcell indices
!
! !USES:
!
! !REVISION HISTORY:
! 03-04-27 : Created by Mariana Vertenstein
! 03-08-25 : Updated to vector data structure (Mariana Vertenstein)
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: g,lev,c                    ! indices
  type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
  type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
  type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
  type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!abt rcm below
  real(r8), pointer :: wtgcell(:)        ! weight of columns relative to gridcells
  integer , pointer :: cgrid(:)          ! gridcell of corresponding column
  real(r8), pointer :: rootfr(:,:)       ! fraction of roots in each soil layer
  integer , pointer :: pcolumn(:)        ! column index of corresponding pft
  real(r8),allocatable :: temp_sm10cm(:) ! temporary array for 10cm soil moisture
  real(r8),allocatable :: temp_sm1m(:)   ! temporary array for 1m soil moisture
  real(r8),allocatable :: temp_smtot(:)  ! temporary array for total soil moisture
  real(r8),allocatable :: temp_root(:,:) ! temporary array for root fraction (nlevsoi)
  integer :: p
!abt rcm above
#if (defined DUST)
  type(pft_dflux_type),pointer :: pdf   ! local pointer to derived subtype
  integer n
#endif

!------------------------------------------------------------------------

  ! Set pointers into derived type

  gptr => clm3%g
  lptr => clm3%g%l
  cptr => clm3%g%l%c
  pptr => clm3%g%l%c%p
!abt rcm below
  wtgcell   => clm3%g%l%c%wtgcell
  cgrid     => clm3%g%l%c%gridcell
  rootfr    => clm3%g%l%c%p%pps%rootfr
  pcolumn   => clm3%g%l%c%p%column
!abt rcm above

  ! Determine processor bounds

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)


  ! Compute gridcell averages. 

  if (present(init)) then

     call c2g(begc, endc, begl, endl, begg, endg, &
          cptr%cws%h2osno, clm_l2a%h2osno, &
          c2l_scale_type= 'unity', l2g_scale_type='unity')
!dir$ concurrent
!cdir nodep
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albd, clm_l2a%albd,&
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albi, clm_l2a%albi,&
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out,&
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
!dir$ concurrent
!cdir nodep
     do g = begg,endg
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
     end do

  else

     call c2g(begc, endc, begl, endl, begg, endg, cptr%cws%h2osno, clm_l2a%h2osno,&
          c2l_scale_type= 'unity', l2g_scale_type='unity')
!dir$ concurrent
!cdir nodep

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albd, clm_l2a%albd, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albi, clm_l2a%albi, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pes%t_ref2m, clm_l2a%t_ref2m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pes%q_ref2m, clm_l2a%q_ref2m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pmf%taux, clm_l2a%taux, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pmf%tauy, clm_l2a%tauy, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_lh_tot, clm_l2a%eflx_lh_tot, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_sh_tot, clm_l2a%eflx_sh_tot, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pwf%qflx_evap_tot, clm_l2a%qflx_evap_tot, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%fsa, clm_l2a%fsa, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
                  
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
                  
#if (defined CN)
     call c2g(begc, endc, begl, endl, begg, endg, &
          cptr%ccf%nee, clm_l2a%nee, &
          c2l_scale_type= 'unity', l2g_scale_type='unity')
#elif (defined CASA)
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pps%co2flux, clm_l2a%nee, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
#else
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pcf%fco2, clm_l2a%nee, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
     ! Note that fco2 in is umolC/m2/sec so units need to be changed to gC/m2/sec
     do g = begg,endg
        clm_l2a%nee(g) = clm_l2a%nee(g)*12.011e-6_r8
     end do
#endif

#if (defined DUST || defined  PROGSSLT )
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pps%fv, clm_l2a%fv, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pps%ram1, clm_l2a%ram1, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
#endif

#if (defined DUST )
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, ndst, &
           pptr%pdf%flx_mss_vrt_dst, clm_l2a%flxdst, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
#endif

     ! Convert from gC/m2/s to kgC/m2/s
!dir$ concurrent
!cdir nodep
     do g = begg,endg
        clm_l2a%nee(g) = clm_l2a%nee(g)*1.0e-3_r8
     end do

!dir$ concurrent
!cdir nodep
     do g = begg,endg
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
     end do

!!!!!! abt rcm below for grid-average for additional variables to be passed to atmo!!!!

     allocate(temp_sm10cm(begg:endg))
     allocate(temp_smtot(begg:endg))
     allocate(temp_sm1m(begg:endg))
     allocate(temp_root(begg:endg,1:nlevsoi))
     temp_sm1m(:)   = 0._r8
     temp_sm10cm(:) = 0._r8
     temp_smtot(:)  = 0._r8
     temp_root(:,:) = 0._r8

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           nlevsoi,rootfr, temp_root, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     do c = begc,endc
        g = cgrid(c)
!! averaging for smtot and sm1m are done without the p2g or c2g subrountines
       do lev = 1,nlevsoi
        if(cptr%cws%h2osoi_liq(c,lev) /= spval) then
          if(lev.le.3) temp_sm10cm(g) = temp_sm10cm(g) + cptr%cws%h2osoi_liq(c,lev) &
                                        * (1000./denh2o) * wtgcell(c) 
          temp_smtot(g)  = temp_smtot(g) + cptr%cws%h2osoi_liq(c,lev) * &
                           (1000./denh2o) * wtgcell(c) 
          if(temp_root(g,lev) > 0._r8) temp_sm1m(g)   = temp_sm1m(g) +  &
                            cptr%cws%h2osoi_liq(c,lev) * (1000./denh2o) *wtgcell(c)
        endif       
       enddo
     enddo

     clm_l2a%smtot(:)  = temp_smtot(:)
     clm_l2a%sm10cm(:) = temp_sm10cm(:)
     clm_l2a%sm1m(:)   = temp_sm1m(:)

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pps%u10, clm_l2a%u10, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pes%t_veg, clm_l2a%tlef, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call c2g(begc, endc, begl, endl, begg, endg, cptr%ces%t_grnd, &
          clm_l2a%tgrnd,c2l_scale_type= 'unity', l2g_scale_type='unity')

     call c2g(begc, endc, begl, endl, begg, endg, cptr%cps%frac_sno, &
          clm_l2a%frac_sno,c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          dble(pptr%pps%frac_veg_nosno), clm_l2a%frac_veg_nosno, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%uvdrag, clm_l2a%uvdrag, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call c2g(begc, endc, begl, endl, begg, endg, cptr%cwf%qflx_infl, &
          clm_l2a%qflx_infl,c2l_scale_type= 'unity', l2g_scale_type='unity')

     call c2g(begc, endc, begl, endl, begg, endg, cptr%cwf%qflx_surf, &
          clm_l2a%qflx_surf,c2l_scale_type= 'unity', l2g_scale_type='unity')

     call c2g(begc, endc, begl, endl, begg, endg, cptr%cwf%qflx_drain, &
          clm_l2a%qflx_drain,c2l_scale_type= 'unity', l2g_scale_type='unity')

#if (defined VOC)
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, nvoc, &
           pptr%pvf%vocflx, clm_l2a%voc2rcm, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
     ! convert from ug/m2/h to kg/m2/s
     clm_l2a%voc2rcm(:,:) = clm_l2a%voc2rcm(:,:)*1.e-9_r8/(3600._r8)


     where(clm_l2a%voc2rcm.ne.clm_l2a%voc2rcm)
        clm_l2a%voc2rcm = 0._r8
     endwhere
#endif

     if( cgaschem == 1 ) then
        !For dry deposition velocity
        call p2g(begp, endp, begc, endc, begl, endl, begg, endg, n_drydep, &
           pptr%pdd%drydepvel, clm_l2a%vdrydep, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
     end if


     deallocate(temp_sm10cm)
     deallocate(temp_sm1m)
     deallocate(temp_smtot)
     deallocate(temp_root)
     call clm2rcm(clm_l2a)

!!!! abt above up for grid averaging !!!!!!!!!!!!!!!

  end if  !if statement for initialize or not

end subroutine clm_map2gcell


! abt rcm below

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm2rcm
!
! !INTERFACE:
  subroutine clm2rcm(clm_l2a)
!
! !DESCRIPTION:
! Maps lnd2atm fields from clm grid to external grid
!
! !USES:
!abt rcm below
  use spmdMod 
  use decompMod        , only : ldecomp, get_proc_clumps, get_clump_bounds
  use decompMod        , only : get_proc_global
  use areaMod          , only : map_maparrayl, map1dl_l2a
  use domainMod        , only : adomain,ldomain
  use clm_varsur
  use clm_varpar       , only : lsmlon,lsmlat
  use clm_drydep       , only : c2r_depout
  use mod_clm
  use mod_dynparam
!
! !ARGUMENTS:
  implicit none
!rcm above
  type(lnd2atm_type), intent(in)  :: clm_l2a
!  type(lnd2atm_type), intent(out) :: l2a_dst
!
! !REVISION HISTORY:
! 2009.1.8    A Tawfik Revised to work with RegCM
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n,ni,nj,nt,nn         ! loop counter
  integer :: nlon,nlat,ntot,ierr   ! abt grid values
  integer :: nclumps               ! field index
  integer :: nflds,nc              ! number of fields to be mapped
  integer :: nradflds              ! size of 2nd dim in arrays
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  integer :: begg,endg,begc,endc   ! beg and end of grid and column
  integer :: begl,endl,begp,endp
  real(r8),allocatable :: c2r_all(:)    ! used to capture all c2r vars 
  integer ,allocatable :: displace(:)   ! used for gathering
  real(r8),allocatable :: dry_local(:)  ! used to capture all c2r vars
  integer :: numg,numl,numc,nump        ! proc totals
  integer :: nout                       ! number of vars going to regcm
  integer :: iv                         ! index run
!------------------------------------------------------------------------------

    ! Determine clump bounds for this processor

    nclumps = get_proc_clumps()
    if(nclumps.gt.1) then
      write(6,*) "ERROR: not formatted for more than 1 clump per processor"
      call endrun()
    endif

    ! Loop over clumps on this processor

    do nc = 1,nclumps

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
     call get_proc_global(numg,numl,numc,nump) 

     nt = endg-begg+1
     nout = 20
     if( ichem == 1 ) then  !Aerosol scheme on
#if (defined VOC)
       if(cgaschem == 1)    nout = nout + 1
#endif
       if(caerosol == 1)    nout = nout + 2
     end if
     if ( .not. allocated(c2r_all)    ) allocate(c2r_all(nt*nout))
     if ( .not. allocated(c2r_allout) ) allocate(c2r_allout(numg*nout))
     if ( .not. allocated(displace)   ) allocate(displace(npes))

     !** Used for dry deposition
     !** for chemistry on only
     if( cgaschem == 1 ) then
        if(.not.allocated(dry_local))  allocate(dry_local(nt*n_drydep))
        if(.not.allocated(c2r_depout)) allocate(c2r_depout(numg * n_drydep))
     end if


!!!!!! below is transfering variables from CLM to RegCM !!!
!!!!!!!!!!!!!!!!!!! by c2r* variables !!!!!!!!!!!!!!!!!!!!!
!!!!!! rest of the transfer from clm to regcm is in !!!!!!!
!!!!!! interf_clm(para/ser).F for inout = 2 !!!!!!!!!!!!!!!
!!!!!! below merely stores all output into one array for !!
!!!!!! quicker gathering                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     nn = 0   
     do n  = begg,endg   
 
        nn = nn + 1
        c2r_all(nn)       = clm_l2a%tgrnd(n)
        c2r_all(nn+nt)    = clm_l2a%h2osno(n)
        c2r_all(nn+nt*2)  = clm_l2a%eflx_sh_tot(n)
        c2r_all(nn+3*nt)  = clm_l2a%qflx_evap_tot(n)
        c2r_all(nn+4*nt)  = clm_l2a%uvdrag(n)
        c2r_all(nn+5*nt)  = clm_l2a%albd(n,1)
        c2r_all(nn+6*nt)  = clm_l2a%albd(n,2)
        c2r_all(nn+7*nt)  = clm_l2a%albi(n,1)
        c2r_all(nn+8*nt)  = clm_l2a%albi(n,2)
        c2r_all(nn+9*nt)  = clm_l2a%t_rad(n)
        c2r_all(nn+10*nt) = clm_l2a%t_ref2m(n)
        c2r_all(nn+11*nt) = clm_l2a%q_ref2m(n)
        c2r_all(nn+12*nt) = clm_l2a%u10(n)
        c2r_all(nn+13*nt) = clm_l2a%tlef(n)
        c2r_all(nn+14*nt) = clm_l2a%sm10cm(n)
        c2r_all(nn+15*nt) = clm_l2a%sm1m(n)
        c2r_all(nn+16*nt) = clm_l2a%smtot(n)
        c2r_all(nn+17*nt) = clm_l2a%qflx_infl(n)*(24._r8 * 3600._r8)
        c2r_all(nn+18*nt) = clm_l2a%qflx_surf(n)*(24._r8 * 3600._r8)
        c2r_all(nn+19*nt) = clm_l2a%qflx_drain(n)*(24._r8 * 3600._r8)

        if(    ichem == 1 ) then
        if( caerosol == 1 ) then
           c2r_all(nn+20*nt) = clm_l2a%frac_sno(n)
           c2r_all(nn+21*nt) = clm_l2a%frac_veg_nosno(n)
        end if
        end if

#if (defined VOC)
        if( cgaschem == 1 ) then
          if( caerosol == 1 ) then  !Aerosol scheme on
            c2r_all(nn+22*nt)   = clm_l2a%voc2rcm(n,1)
          else
            c2r_all(nn+20*nt)   = clm_l2a%voc2rcm(n,1)
          end if
        end if
#endif

	!***** Dry deposition
	!** for chemistry on only
	if( cgaschem == 1 ) then
          do iv=1,n_drydep
            dry_local(nn + (nt*(iv-1))) =  clm_l2a%vdrydep(n,iv)
          end do
        end if

        
    enddo !gridcell
    enddo !nclumps


    !****** Gather CLM variables to the RegCM grid
    do nc = 1,npes
      if(nc.eq.1) then
        displace(nc) = 0
      else
        displace(nc) = displace(nc-1) + c2rngc(nc-1)*nout
      endif
    enddo  

    if(nn .ne. nt) then
      write(6,*) "ERROR: nn not equal to nt"
      call endrun()
    endif

    call mpi_allgatherv(c2r_all(1),nn*nout,MPI_REAL8,c2r_allout, &
                        c2rngc*nout,displace,MPI_REAL8,mpicom,ierr)     
    


    !****** Gather Dry Dep CLM variables to the RegCM grid
    !** for chemistry on only
    if(cgaschem == 1) then
      do nc = 1,npes
        if(nc.eq.1) then
          displace(nc) = 0
        else
          displace(nc) = displace(nc-1) + c2rngc(nc-1)*n_drydep
        endif
      enddo
      call mpi_allgatherv(dry_local(1),nn*n_drydep,MPI_REAL8, &
                          c2r_depout, c2rngc*n_drydep, displace, &
                          MPI_REAL8,mpicom,ierr)
      deallocate(dry_local)
    end if


    deallocate(c2r_all,displace)        



end subroutine clm2rcm


!!!!!! abt rcm above



!------------------------------------------------------------------------
!------------------------------------------------------------------------
end module clm_atmlnd


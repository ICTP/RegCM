module seq_flds_indices

  use seq_flds_mod  
  use shr_string_mod
  use shr_sys_mod, only: shr_sys_abort

  implicit none
  save
  public

  ! atm -> drv

  integer :: index_a2x_Sa_z            ! bottom atm level height
  integer :: index_a2x_Sa_u            ! bottom atm level zon wind
  integer :: index_a2x_Sa_v            ! bottom atm level mer wind
  integer :: index_a2x_Sa_tbot         ! bottom atm level temp
  integer :: index_a2x_Sa_ptem         ! bottom atm level pot temp
  integer :: index_a2x_Sa_shum         ! bottom atm level spec hum
  integer :: index_a2x_Sa_dens         ! bottom atm level air den
  integer :: index_a2x_Sa_pbot         ! bottom atm level pressure
  integer :: index_a2x_Sa_pslv         ! sea level atm pressure
  integer :: index_a2x_Faxa_lwdn       ! downward lw heat flux
  integer :: index_a2x_Faxa_rainc      ! prec: liquid "convective"
  integer :: index_a2x_Faxa_rainl      ! prec: liquid "large scale"
  integer :: index_a2x_Faxa_snowc      ! prec: frozen "convective"
  integer :: index_a2x_Faxa_snowl      ! prec: frozen "large scale"
  integer :: index_a2x_Faxa_swndr      ! sw: nir direct  downward
  integer :: index_a2x_Faxa_swvdr      ! sw: vis direct  downward
  integer :: index_a2x_Faxa_swndf      ! sw: nir diffuse downward
  integer :: index_a2x_Faxa_swvdf      ! sw: vis diffuse downward
  integer :: index_a2x_Faxa_swnet      ! sw: net
  integer :: index_a2x_Sa_co2prog      ! bottom atm level prognostic co2
  integer :: index_a2x_Sa_co2diag      ! bottom atm level diagnostic co2
  integer :: nflds_a2x

  ! drv -> atm

  integer :: index_x2a_Sx_t            ! surface temperature             
  integer :: index_x2a_So_t            ! sea surface temperature         
  integer :: index_x2a_Sx_lfrac        ! surface land fraction           
  integer :: index_x2a_Sx_ifrac        ! surface ice fraction            
  integer :: index_x2a_Sx_ofrac        ! surface ocn fraction            
  integer :: index_x2a_Sx_tref         ! 2m reference temperature        
  integer :: index_x2a_Sx_qref         ! 2m reference specific humidity  
  integer :: index_x2a_Sx_avsdr        ! albedo, visible, direct         
  integer :: index_x2a_Sx_anidr        ! albedo, near-ir, direct         
  integer :: index_x2a_Sx_avsdf        ! albedo, visible, diffuse        
  integer :: index_x2a_Sx_anidf        ! albedo, near-ir, diffuse        
  integer :: index_x2a_Sx_snowh        ! surface snow depth              
  integer :: index_x2a_Sx_fv           ! friction velocity
  integer :: index_x2a_Sx_ram1         ! aerodynamical resistance
  integer :: index_x2a_Faxx_taux       ! wind stress, zonal              
  integer :: index_x2a_Faxx_tauy       ! wind stress, meridional         
  integer :: index_x2a_Faxx_lat        ! latent          heat flux       
  integer :: index_x2a_Faxx_sen        ! sensible        heat flux       
  integer :: index_x2a_Faxx_lwup       ! upward longwave heat flux       
  integer :: index_x2a_Faxx_evap       ! evaporation    water flux       
  integer :: index_x2a_Faxx_flxdst1    ! dust flux size bin 1    
  integer :: index_x2a_Faxx_flxdst2    ! dust flux size bin 2    
  integer :: index_x2a_Faxx_flxdst3    ! dust flux size bin 3    
  integer :: index_x2a_Faxx_flxdst4    ! dust flux size bin 4
  integer :: index_x2a_So_ustar	
  integer :: index_x2a_So_re
  integer :: index_x2a_So_ssq
  !
  integer :: index_x2a_Si_avsdr        ! albedo, visible, direct         
  integer :: index_x2a_Si_anidr        ! albedo, near-ir, direct         
  integer :: index_x2a_Si_avsdf        ! albedo, visible, diffuse        
  integer :: index_x2a_Si_anidf        ! albedo, near-ir, diffuse        
  integer :: index_x2a_So_avsdr        ! albedo, visible, direct   ! (flux module)
  integer :: index_x2a_So_anidr        ! albedo, near-ir, direct   ! (flux module)
  integer :: index_x2a_So_avsdf        ! albedo, visible, diffuse  ! (flux module)
  integer :: index_x2a_So_anidf        ! albedo, near-ir, diffuse  ! (flux module)
  integer :: index_x2a_Faii_lat        ! latent          heat flux       
  integer :: index_x2a_Faii_sen        ! sensible        heat flux       
  integer :: index_x2a_Faii_lwup       ! upward longwave heat flux       
  integer :: index_x2a_Faox_lat        ! latent          heat flux ! (flux module)
  integer :: index_x2a_Faox_sen        ! sensible        heat flux ! (flux module)      
  integer :: index_x2a_Faox_lwup       ! upward longwave heat flux       
  !
  integer :: nflds_x2a

  ! ice -> drv

  integer :: index_i2x_Si_ifrac        ! fractional ice coverage wrt ocean
  integer :: index_i2x_Si_sicthk       ! sea ice thickness (needed only for cam/som)
  integer :: index_i2x_Si_t            ! temperature                     
  integer :: index_i2x_Si_tref         ! 2m reference temperature        
  integer :: index_i2x_Si_qref         ! 2m reference specific humidity  
  integer :: index_i2x_Si_avsdr        ! albedo: visible, direct         
  integer :: index_i2x_Si_avsdf        ! albedo: near ir, direct         
  integer :: index_i2x_Si_anidr        ! albedo: visible, diffuse        
  integer :: index_i2x_Si_anidf        ! albedo: near ir, diffuse        
  integer :: index_i2x_Faii_lwup       ! upward longwave heat flux  
  integer :: index_i2x_Faii_lat        ! latent          heat flux  
  integer :: index_i2x_Faii_sen        ! sensible        heat flux      
  integer :: index_i2x_Faii_evap       ! evaporation    water flux      
  integer :: index_i2x_Faii_taux       ! wind stress, zonal            
  integer :: index_i2x_Faii_tauy       ! wind stress, meridional       
  integer :: index_i2x_Faii_swnet      ! sw: net
  integer :: index_i2x_Fioi_swpen      ! sw: net penetrating ice
  integer :: index_i2x_Fioi_melth      ! heat  flux from melting ice (<0)
  integer :: index_i2x_Fioi_meltw      ! water flux from melting ice
  integer :: index_i2x_Fioi_salt       ! salt  flux from meting  ice
  integer :: index_i2x_Fioi_taux       ! ice/ocn stress, zonal
  integer :: index_i2x_Fioi_tauy       ! ice/ocn stress, zonal
  integer :: nflds_i2x

  ! drv -> ice

  integer :: index_x2i_So_t            ! sea surface temperature
  integer :: index_x2i_Sa_z            ! bottom atm level height
  integer :: index_x2i_Sa_u            ! bottom atm level zon wind
  integer :: index_x2i_Sa_v            ! bottom atm level mer wind
  integer :: index_x2i_Sa_tbot         ! bottom atm level temp
  integer :: index_x2i_Sa_pbot         ! bottom atm level pressure
  integer :: index_x2i_Sa_ptem         ! bottom atm level pot temp
  integer :: index_x2i_Sa_shum         ! bottom atm level spec hum
  integer :: index_x2i_Sa_dens         ! bottom atm level air den
  integer :: index_x2i_Faxa_lwdn       ! downward lw heat flux
  integer :: index_x2i_Faxa_rain       ! prec: liquid 
  integer :: index_x2i_Faxa_snow       ! prec: frozen 
  integer :: index_x2i_Faxa_swndr      ! sw: nir direct  downward
  integer :: index_x2i_Faxa_swvdr      ! sw: vis direct  downward
  integer :: index_x2i_Faxa_swndf      ! sw: nir diffuse downward
  integer :: index_x2i_Faxa_swvdf      ! sw: vis diffuse downward
  integer :: index_x2i_Faxa_swnet      ! sw: net
  integer :: index_x2i_Fioo_q          ! ocn freeze or melt heat  
  integer :: nflds_x2i

  ! hub atm/ocn flufxes and states (computed by flux module)	

  integer :: index_xao_So_tref    
  integer :: index_xao_So_qref    
  integer :: index_xao_So_avsdr   
  integer :: index_xao_So_avsdf   
  integer :: index_xao_So_anidr   
  integer :: index_xao_So_anidf   
  integer :: index_xao_Sx_duu10n 
  integer :: index_xao_Faox_taux  
  integer :: index_xao_Faox_tauy   
  integer :: index_xao_Faox_lat   
  integer :: index_xao_Faox_sen   
  integer :: index_xao_Faox_evap  
  integer :: index_xao_Faox_lwup
  integer :: index_xao_Faox_swnet 
  integer :: index_xao_So_ustar         ! optional
  integer :: index_xao_So_re            ! optional
  integer :: index_xao_So_ssq           ! optional
  integer :: nflds_xao

  ! ocn -> drv

  integer :: index_o2x_So_t      
  integer :: index_o2x_So_u
  integer :: index_o2x_So_v
  integer :: index_o2x_So_s
  integer :: index_o2x_So_dhdx
  integer :: index_o2x_So_dhdy
  integer :: index_o2x_Fioo_q
  integer :: nflds_o2x	

  ! drv -> ocn

  integer :: index_x2o_Si_ifrac        ! fractional ice wrt ocean
  integer :: index_x2o_Si_sicthk       ! needed only for cam-som
  integer :: index_x2o_Sx_duu10n
  integer :: index_x2o_Sa_pslv
  integer :: index_x2o_Foxx_taux  
  integer :: index_x2o_Foxx_tauy  
  integer :: index_x2o_Foxx_swnet 
  integer :: index_x2o_Foxx_sen   
  integer :: index_x2o_Foxx_lat   
  integer :: index_x2o_Foxx_lwdn  
  integer :: index_x2o_Foxx_lwup
  integer :: index_x2o_Foxx_melth 
  integer :: index_x2o_Foxx_salt
  integer :: index_x2o_Foxx_prec
  integer :: index_x2o_Foxx_snow
  integer :: index_x2o_Foxx_rain
  integer :: index_x2o_Foxx_evap
  integer :: index_x2o_Foxx_meltw 
  integer :: index_x2o_Forr_roff
  integer :: nflds_x2o

  ! lnd -> drv

  integer :: index_l2x_Sl_landfrac     ! land fraction
  integer :: index_l2x_Sl_t            ! temperature
  integer :: index_l2x_Sl_tref         ! 2m reference temperature
  integer :: index_l2x_Sl_qref         ! 2m reference specific humidity
  integer :: index_l2x_Sl_avsdr        ! albedo: direct , visible
  integer :: index_l2x_Sl_anidr        ! albedo: direct , near-ir
  integer :: index_l2x_Sl_avsdf        ! albedo: diffuse, visible
  integer :: index_l2x_Sl_anidf        ! albedo: diffuse, near-ir
  integer :: index_l2x_Sl_snowh        ! snow height
  integer :: index_l2x_Fall_taux       ! wind stress, zonal
  integer :: index_l2x_Fall_tauy       ! wind stress, meridional
  integer :: index_l2x_Fall_lat        ! latent          heat flux
  integer :: index_l2x_Fall_sen        ! sensible        heat flux
  integer :: index_l2x_Fall_lwup       ! upward longwave heat flux
  integer :: index_l2x_Fall_evap       ! evaporation     water flux
  integer :: index_l2x_Fall_swnet      ! heat flux       shortwave net       
  integer :: index_l2x_Fall_nee        ! co2 flux **For testing set to 0
  integer :: index_l2x_Sl_fv           ! friction velocity  
  integer :: index_l2x_Sl_ram1         ! aerodynamical resistance
  integer :: index_l2x_Fall_flxdst1    ! dust flux size bin 1    
  integer :: index_l2x_Fall_flxdst2    ! dust flux size bin 2    
  integer :: index_l2x_Fall_flxdst3    ! dust flux size bin 3    
  integer :: index_l2x_Fall_flxdst4    ! dust flux size bin 4
  integer :: nflds_l2x

  ! roff to driver (part of land for now)

  integer :: index_r2x_Forr_roff       ! runoff to ocean

  ! atm -> drv

  integer :: index_x2l_Sa_z            ! bottom atm level height
  integer :: index_x2l_Sa_u            ! bottom atm level zon wind
  integer :: index_x2l_Sa_v            ! bottom atm level mer wind
  integer :: index_x2l_Sa_ptem         ! bottom atm level pot temp
  integer :: index_x2l_Sa_shum         ! bottom atm level spec hum
  integer :: index_x2l_Sa_pbot         ! bottom atm level pressure
  integer :: index_x2l_Sa_tbot         ! bottom atm level temp
  integer :: index_x2l_Faxa_lwdn       ! downward lw heat flux
  integer :: index_x2l_Faxa_rainc      ! prec: liquid "convective"
  integer :: index_x2l_Faxa_rainl      ! prec: liquid "large scale"
  integer :: index_x2l_Faxa_snowc      ! prec: frozen "convective"
  integer :: index_x2l_Faxa_snowl      ! prec: frozen "large scale"
  integer :: index_x2l_Faxa_swndr      ! sw: nir direct  downward
  integer :: index_x2l_Faxa_swvdr      ! sw: vis direct  downward
  integer :: index_x2l_Faxa_swndf      ! sw: nir diffuse downward
  integer :: index_x2l_Faxa_swvdf      ! sw: vis diffuse downward
  integer :: index_x2l_Faxa_swnet      ! sw: net
  integer :: index_x2l_Sa_co2prog      ! bottom atm level prognostic co2
  integer :: index_x2l_Sa_co2diag      ! bottom atm level diagnostic co2
  integer :: nflds_x2l

contains

  subroutine seq_flds_indices_set

    !-------------------------------------------------------------
    ! atm -> drv
    !-------------------------------------------------------------

    index_a2x_Sa_z          = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_z')
    if (index_a2x_Sa_z == 0) call shr_sys_abort('index_a2x_Sa_z is zero')

    index_a2x_Sa_u          = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_u')
    if (index_a2x_Sa_u == 0) call shr_sys_abort('index_a2x_Sa_u is zero')

    index_a2x_Sa_v          = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_v')
    if (index_a2x_Sa_v == 0) call shr_sys_abort('index_a2x_Sa_v is zero')

    index_a2x_Sa_tbot       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_tbot')
    if (index_a2x_Sa_tbot == 0) call shr_sys_abort('index_a2x_Sa_tbot is zero')

    index_a2x_Sa_ptem       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_ptem')
    if (index_a2x_Sa_ptem == 0) call shr_sys_abort('index_a2x_Sa_ptem is zero')

    index_a2x_Sa_pbot       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_pbot')
    if (index_a2x_Sa_pbot == 0) call shr_sys_abort('index_a2x_Sa_pbot is zero')

    index_a2x_Sa_pslv       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_pslv')
    if (index_a2x_Sa_pslv == 0) call shr_sys_abort('index_a2x_Sa_pslv is zero')

    index_a2x_Sa_shum       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_shum')
    if (index_a2x_Sa_shum == 0) call shr_sys_abort('index_a2x_Sa_shum is zero')

    index_a2x_Sa_dens       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_dens')
    if (index_a2x_Sa_dens == 0) call shr_sys_abort('index_a2x_Sa_dens is zero')

    index_a2x_Faxa_swnet    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swnet')
    if (index_a2x_Faxa_swnet == 0) call shr_sys_abort('index_a2x_Faxa_swnet is zero')

    index_a2x_Faxa_lwdn     = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_lwdn')
    if (index_a2x_Faxa_lwdn == 0) call shr_sys_abort('index_a2x_Faxa_lwdn is zero')

    index_a2x_Faxa_rainc    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_rainc')
    if (index_a2x_Faxa_rainc == 0) call shr_sys_abort('index_a2x_Faxa_rainc is zero')

    index_a2x_Faxa_rainl    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_rainl')
    if (index_a2x_Faxa_rainl == 0) call shr_sys_abort('index_a2x_Faxa_rainl is zero')

    index_a2x_Faxa_snowc    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_snowc')
    if (index_a2x_Faxa_snowc == 0) call shr_sys_abort('index_a2x_Faxa_snowc is zero')

    index_a2x_Faxa_snowl    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_snowl')
    if (index_a2x_Faxa_snowl == 0) call shr_sys_abort('index_a2x_Faxa_snowl is zero')

    index_a2x_Faxa_swndr    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swndr')
    if (index_a2x_Faxa_swndr == 0) call shr_sys_abort('index_a2x_Faxa_swndr is zero')

    index_a2x_Faxa_swvdr    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swvdr')
    if (index_a2x_Faxa_swvdr == 0) call shr_sys_abort('index_a2x_Faxa_swvdr is zero')

    index_a2x_Faxa_swndf    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swndf')
    if (index_a2x_Faxa_swndf == 0) call shr_sys_abort('index_a2x_Faxa_swndf is zero')

    index_a2x_Faxa_swvdf    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swvdf')
    if (index_a2x_Faxa_swvdf == 0) call shr_sys_abort('index_a2x_Faxa_swvdf is zero')

    index_a2x_Sa_co2prog    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_co2prog')
    index_a2x_Sa_co2diag    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_co2diag')

    nflds_a2x = shr_string_listGetNum(seq_flds_a2x_fields)

    !-------------------------------------------------------------
    ! drv -> atm
    !-------------------------------------------------------------

    index_x2a_Sx_avsdr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_avsdr')
    if (index_x2a_Sx_avsdr == 0) call shr_sys_abort('index_x2a_Sx_avsdr is zero')

    index_x2a_Sx_anidr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_anidr')
    if (index_x2a_Sx_anidr == 0) call shr_sys_abort('index_x2a_Sx_anidr is zero')

    index_x2a_Sx_avsdf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_avsdf')
    if (index_x2a_Sx_avsdf == 0) call shr_sys_abort('index_x2a_Sx_avsdf is zero')

    index_x2a_Sx_anidf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_anidf')
    if (index_x2a_Sx_anidf == 0) call shr_sys_abort('index_x2a_Sx_anidf is zero')

    index_x2a_Sx_t          = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_t')
    if (index_x2a_Sx_t == 0) call shr_sys_abort('index_x2a_Sx_t is zero')

    index_x2a_So_t          = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_t')
    if (index_x2a_So_t == 0) call shr_sys_abort('index_x2a_So_t is zero')

    index_x2a_Sx_snowh      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_snowh')
    if (index_x2a_Sx_snowh == 0) call shr_sys_abort('index_x2a_Sx_snowh is zero')

    index_x2a_Sx_tref       = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_tref')
    if (index_x2a_Sx_tref == 0) call shr_sys_abort('index_x2a_Sx_tref is zero')

    index_x2a_Sx_qref       = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_qref')
    if (index_x2a_Sx_qref == 0) call shr_sys_abort('index_x2a_Sx_qref is zero')

    index_x2a_Sx_ifrac      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_ifrac')
    if (index_x2a_Sx_ifrac == 0) call shr_sys_abort('index_x2a_Sx_ifrac is zero')

    index_x2a_Sx_ofrac      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_ofrac')
    if (index_x2a_Sx_ofrac == 0) call shr_sys_abort('index_x2a_Sx_ofrac is zero')

    index_x2a_Sx_lfrac      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_lfrac')
    if (index_x2a_Sx_lfrac == 0) call shr_sys_abort('index_x2a_Sx_lfrac is zero')

    index_x2a_Faxx_taux     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_taux')
    if (index_x2a_Faxx_taux == 0) call shr_sys_abort('index_x2a_Faxx_taux is zero')

    index_x2a_Faxx_tauy     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_tauy')
    if (index_x2a_Faxx_tauy == 0) call shr_sys_abort('index_x2a_Faxx_tauy is zero')

    index_x2a_Faxx_lat      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_lat')
    if (index_x2a_Faxx_lat == 0) call shr_sys_abort('index_x2a_Faxx_lat is zero')

    index_x2a_Faxx_sen      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_sen')
    if (index_x2a_Faxx_sen == 0) call shr_sys_abort('index_x2a_Faxx_sen is zero')

    index_x2a_Faxx_lwup     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_lwup')
    if (index_x2a_Faxx_lwup == 0) call shr_sys_abort('index_x2a_Faxx_lwup is zero')

    index_x2a_Faxx_evap     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_evap')
    if (index_x2a_Faxx_evap == 0) call shr_sys_abort('index_x2a_Faxx_evap is zero')

    ! fields needed to calculate water isotopes to ocean evaporation processes

    index_x2a_So_ustar     = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_ustar')
    if (index_x2a_So_ustar == 0) call shr_sys_abort('index_x2a_So_ustar is zero')
    
    index_x2a_So_re     = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_re')
    if (index_x2a_So_re == 0) call shr_sys_abort('index_x2a_So_re is zero')

    index_x2a_So_ssq     = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_ssq')
    if (index_x2a_So_ssq == 0) call shr_sys_abort('index_x2a_So_ssq is zero')

    ! optional fields

    index_x2a_Si_avsdr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Si_avsdr')
    index_x2a_si_anidr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Si_anidr')
    index_x2a_si_avsdf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Si_avsdf')
    index_x2a_si_anidf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Si_anidf')
    index_x2a_So_avsdr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_avsdr')
    index_x2a_So_anidr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_anidr')
    index_x2a_So_avsdf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_avsdf')
    index_x2a_So_anidf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_anidf')
    index_x2a_Faox_lat      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faox_lat')
    index_x2a_Faox_sen      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faox_sen')
    index_x2a_Faox_lwup     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faox_lwup')
    index_x2a_Faii_lat      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faii_lat')
    index_x2a_Faii_sen      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faii_sen')
    index_x2a_Faii_lwup     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faii_lwup')
    index_x2a_Sx_fv         = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_fv')
    index_x2a_Sx_ram1       = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_ram1')
    index_x2a_Faxx_flxdst1  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_flxdst1')
    index_x2a_Faxx_flxdst2  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_flxdst2')
    index_x2a_Faxx_flxdst3  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_flxdst3')
    index_x2a_Faxx_flxdst4  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_flxdst4')

    nflds_x2a = shr_string_listGetNum(seq_flds_x2a_fields)

    !-------------------------------------------------------------
    ! ocn -> drv
    !-------------------------------------------------------------

    index_o2x_So_t          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_t')
    if (index_o2x_So_t == 0) call shr_sys_abort('index_o2x_So_t is zero')

    index_o2x_So_u          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_u')
    if (index_o2x_So_u == 0) call shr_sys_abort('index_o2x_So_u is zero')

    index_o2x_So_v          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_v')
    if (index_o2x_So_v == 0) call shr_sys_abort('index_o2x_So_v is zero')

    index_o2x_So_s          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_s')
    if (index_o2x_So_s == 0) call shr_sys_abort('index_o2x_So_s is zero')

    index_o2x_So_dhdx          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_dhdx')
    if (index_o2x_So_dhdx == 0) call shr_sys_abort('index_o2x_So_dhdx is zero')

    index_o2x_So_dhdy          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_dhdy')
    if (index_o2x_So_dhdy == 0) call shr_sys_abort('index_o2x_So_dhdy is zero')

    index_o2x_Fioo_q        = shr_string_listGetIndexF(seq_flds_o2x_fields,'Fioo_q')
    if (index_o2x_Fioo_q == 0) call shr_sys_abort('index_o2x_Fioo_q is zero')

    nflds_o2x = shr_string_listGetNum(seq_flds_o2x_fields)

    !-------------------------------------------------------------
    ! hub atm/ocn fluxes/states
    !-------------------------------------------------------------

    index_xao_So_tref       = shr_string_listGetIndexF(seq_flds_xao_fields,'So_tref')
    if (index_xao_So_tref == 0) call shr_sys_abort('index_xao_So_tref is zero')

    index_xao_So_qref       = shr_string_listGetIndexF(seq_flds_xao_fields,'So_qref')
    if (index_xao_So_qref == 0) call shr_sys_abort('index_xao_So_qref is zero')

    index_xao_So_avsdr      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_avsdr')
    if (index_xao_So_avsdr == 0) call shr_sys_abort('index_xao_So_avsdr is zero')

    index_xao_So_anidr      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_anidr')
    if (index_xao_So_anidr == 0) call shr_sys_abort('index_xao_So_anidr is zero')

    index_xao_So_avsdf      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_avsdf')
    if (index_xao_So_avsdf == 0) call shr_sys_abort('index_xao_So_avsdf is zero')

    index_xao_So_anidf      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_anidf')
    if (index_xao_So_anidf == 0) call shr_sys_abort('index_xao_So_anidf is zero')

    index_xao_Sx_duu10n      = shr_string_listGetIndexF(seq_flds_xao_fields,'Sx_duu10n')
    if (index_xao_Sx_duu10n == 0) call shr_sys_abort('index_xao_Sx_duu10n is zero')

    index_xao_So_ustar      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_ustar')
    if (index_xao_So_ustar == 0) call shr_sys_abort('index_xao_So_ustar is zero')

    index_xao_So_re      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_re')
    if (index_xao_So_re == 0) call shr_sys_abort('index_xao_So_re is zero')

    index_xao_So_ssq      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_ssq')
    if (index_xao_So_ssq == 0) call shr_sys_abort('index_xao_So_ssq is zero')

    index_xao_Faox_taux     = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_taux')
    if (index_xao_Faox_taux == 0) call shr_sys_abort('index_xao_Faox_taux is zero')

    index_xao_Faox_tauy     = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_tauy')
    if (index_xao_Faox_tauy == 0) call shr_sys_abort('index_xao_Faox_tauy is zero')

    index_xao_Faox_lat      = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_lat')
    if (index_xao_Faox_lat == 0) call shr_sys_abort('index_xao_Faox_lat is zero')

    index_xao_Faox_sen      = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_sen')
    if (index_xao_Faox_sen == 0) call shr_sys_abort('index_xao_Faox_sen is zero')

    index_xao_Faox_evap     = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_evap')
    if (index_xao_Faox_evap == 0) call shr_sys_abort('index_xao_Faox_evap is zero')

    index_xao_Faox_lwup     = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_lwup')
    if (index_xao_Faox_lwup == 0) call shr_sys_abort('index_xao_Faox_lwup is zero')

    index_xao_Faox_swnet    = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_swnet')
    if (index_xao_Faox_swnet == 0) call shr_sys_abort('index_xao_Faox_swnet is zero')

    nflds_xao = shr_string_listGetNum(seq_flds_xao_fields)

    !-------------------------------------------------------------
    ! drv -> ocn
    !-------------------------------------------------------------

    index_x2o_Si_ifrac      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Si_ifrac')
    if (index_x2o_Si_ifrac == 0) call shr_sys_abort('index_x2o_Si_ifrac is zero')

    index_x2o_Si_sicthk     = shr_string_listGetIndexF(seq_flds_x2o_fields,'Si_sicthk')
    if (index_x2o_Si_sicthk == 0) call shr_sys_abort('index_x2o_Si_sicthk is zero')

    index_x2o_Sa_pslv      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Sa_pslv')
    if (index_x2o_Sa_pslv == 0) call shr_sys_abort('index_x2o_Sa_pslv is zero')

    index_x2o_Sx_duu10n      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Sx_duu10n')
    if (index_x2o_Sx_duu10n == 0) call shr_sys_abort('index_x2o_Sx_duu10n is zero')

    index_x2o_Foxx_tauy     = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_tauy')
    if (index_x2o_Foxx_tauy == 0) call shr_sys_abort('index_x2o_Foxx_tauy is zero')

    index_x2o_Foxx_taux     = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_taux')
    if (index_x2o_Foxx_taux == 0) call shr_sys_abort('index_x2o_Foxx_taux is zero')

    index_x2o_Foxx_swnet    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_swnet')
    if (index_x2o_Foxx_swnet == 0) call shr_sys_abort('index_x2o_Foxx_swnet is zero')

    index_x2o_Foxx_lat      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_lat')
    if (index_x2o_Foxx_lat == 0) call shr_sys_abort('index_x2o_Foxx_lat is zero')

    index_x2o_Foxx_sen      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_sen')
    if (index_x2o_Foxx_sen == 0) call shr_sys_abort('index_x2o_Foxx_sen is zero')

    index_x2o_Foxx_lwdn     = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_lwdn')
    if (index_x2o_Foxx_lwdn == 0) call shr_sys_abort('index_x2o_Foxx_lwdn is zero')

    index_x2o_Foxx_lwup     = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_lwup')
    if (index_x2o_Foxx_lwup == 0) call shr_sys_abort('index_x2o_Foxx_lwup is zero')

    index_x2o_Foxx_melth    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_melth')   
    if (index_x2o_Foxx_melth == 0) call shr_sys_abort('index_x2o_Foxx_melth is zero')

    index_x2o_Foxx_salt    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_salt')   
    if (index_x2o_Foxx_salt == 0) call shr_sys_abort('index_x2o_Foxx_salt is zero')

    index_x2o_Foxx_prec    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_prec')   
    if (index_x2o_Foxx_prec == 0) call shr_sys_abort('index_x2o_Foxx_prec is zero')

    index_x2o_Foxx_snow    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_snow')   
    if (index_x2o_Foxx_snow == 0) call shr_sys_abort('index_x2o_Foxx_snow is zero')

    index_x2o_Foxx_rain    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_rain')   
    if (index_x2o_Foxx_rain == 0) call shr_sys_abort('index_x2o_Foxx_rain is zero')

    index_x2o_Foxx_evap    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_evap')
    if (index_x2o_Foxx_evap == 0) call shr_sys_abort('index_x2o_Foxx_evap is zero')

    index_x2o_Foxx_meltw    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_meltw')
    if (index_x2o_Foxx_meltw == 0) call shr_sys_abort('index_x2o_Foxx_meltw is zero')

    index_x2o_Forr_roff    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Forr_roff')
    if (index_x2o_Forr_roff == 0) call shr_sys_abort('index_x2o_Forr_roff is zero')

    nflds_x2o = shr_string_listGetNum(seq_flds_x2o_fields)
    
    !-------------------------------------------------------------
    ! ice -> drv
    !-------------------------------------------------------------

    index_i2x_Si_t          = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_t')
    if (index_i2x_Si_t == 0) call shr_sys_abort('index_i2x_Si_t is zero')

    index_i2x_Si_tref       = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_tref')
    if (index_i2x_Si_tref == 0) call shr_sys_abort('index_i2x_Si_tref is zero')

    index_i2x_Si_qref       = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_qref')
    if (index_i2x_Si_qref == 0) call shr_sys_abort('index_i2x_Si_qref is zero')

    index_i2x_Si_ifrac      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_ifrac')
    if (index_i2x_Si_ifrac == 0) call shr_sys_abort('index_i2x_Si_ifrac is zero')

    index_i2x_Si_avsdr      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_avsdr')
    if (index_i2x_Si_avsdr == 0) call shr_sys_abort('index_i2x_Si_avsdr is zero')

    index_i2x_Si_anidr      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_anidr')
    if (index_i2x_Si_anidr == 0) call shr_sys_abort('index_i2x_Si_anidr is zero')

    index_i2x_Si_avsdf      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_avsdf')
    if (index_i2x_Si_avsdf == 0) call shr_sys_abort('index_i2x_Si_avsdf is zero')

    index_i2x_Si_anidf      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_anidf')
    if (index_i2x_Si_anidf == 0) call shr_sys_abort('index_i2x_Si_anidf is zero')

    index_i2x_Si_sicthk     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_sicthk')
    if (index_i2x_Si_sicthk == 0) call shr_sys_abort('index_i2x_Si_sicthk is zero')

    index_i2x_Faii_taux     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_taux')
    if (index_i2x_Faii_taux == 0) call shr_sys_abort('index_i2x_Faii_taux is zero')

    index_i2x_Faii_tauy     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_tauy')
    if (index_i2x_Faii_tauy == 0) call shr_sys_abort('index_i2x_Faii_tauy is zero')

    index_i2x_Faii_lat      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_lat')
    if (index_i2x_Faii_lat == 0) call shr_sys_abort('index_i2x_Faii_lat is zero')

    index_i2x_Faii_sen      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_sen')
    if (index_i2x_Faii_sen == 0) call shr_sys_abort('index_i2x_Faii_sen is zero')

    index_i2x_Faii_lwup     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_lwup')
    if (index_i2x_Faii_lwup == 0) call shr_sys_abort('index_i2x_Faii_lwup is zero')

    index_i2x_Faii_evap     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_evap')
    if (index_i2x_Faii_evap == 0) call shr_sys_abort('index_i2x_Faii_evap is zero')

    index_i2x_Faii_swnet    = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_swnet')
    if (index_i2x_Faii_swnet == 0) call shr_sys_abort('index_i2x_Faii_swnet is zero')

    index_i2x_Fioi_swpen    = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_swpen')
    if (index_i2x_Fioi_swpen == 0) call shr_sys_abort('index_i2x_Fioi_swpen is zero')

    index_i2x_Fioi_melth    = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_melth')
    if (index_i2x_Fioi_melth == 0) call shr_sys_abort('index_i2x_Fioi_melth is zero')

    index_i2x_Fioi_meltw    = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_meltw')
    if (index_i2x_Fioi_meltw == 0) call shr_sys_abort('index_i2x_Fioi_meltw is zero')

    index_i2x_Fioi_salt     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_salt')
    if (index_i2x_Fioi_salt == 0) call shr_sys_abort('index_i2x_Fioi_salt is zero')

    index_i2x_Fioi_taux     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_taux')
    if (index_i2x_Fioi_taux == 0) call shr_sys_abort('index_i2x_Fioi_taux is zero')

    index_i2x_Fioi_tauy     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_tauy')
    if (index_i2x_Fioi_tauy == 0) call shr_sys_abort('index_i2x_Fioi_tauy is zero')

    nflds_i2x               = shr_string_listGetNum(seq_flds_i2x_fields)

    !-------------------------------------------------------------
    ! drv -> ice
    !-------------------------------------------------------------

    index_x2i_So_t          = shr_string_listGetIndexF(seq_flds_x2i_fields,'So_t')
    if (index_x2i_So_t == 0) call shr_sys_abort('index_x2i_So_t is zero')

    index_x2i_Sa_z          = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_z')
    if (index_x2i_Sa_z == 0) call shr_sys_abort('index_x2i_Sa_z is zero')

    index_x2i_Sa_u          = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_u')
    if (index_x2i_Sa_u == 0) call shr_sys_abort('index_x2i_Sa_u is zero')

    index_x2i_Sa_v          = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_v')
    if (index_x2i_Sa_v == 0) call shr_sys_abort('index_x2i_Sa_v is zero')

    index_x2i_Sa_tbot       = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_tbot')
    if (index_x2i_Sa_tbot == 0) call shr_sys_abort('index_x2i_Sa_tbot is zero')

    index_x2i_Sa_ptem       = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_ptem')
    if (index_x2i_Sa_ptem == 0) call shr_sys_abort('index_x2i_Sa_ptem is zero')

    index_x2i_Sa_pbot       = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_pbot')
    if (index_x2i_Sa_pbot == 0) call shr_sys_abort('index_x2i_Sa_pbot is zero')

    index_x2i_Sa_shum       = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_shum')
    if (index_x2i_Sa_shum == 0) call shr_sys_abort('index_x2i_Sa_shum is zero')

    index_x2i_Sa_dens       = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_dens')
    if (index_x2i_Sa_dens == 0) call shr_sys_abort('index_x2i_Sa_dens is zero')

    index_x2i_Faxa_lwdn     = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_lwdn')
    if (index_x2i_Faxa_lwdn == 0) call shr_sys_abort('index_x2i_Faxa_lwdn is zero')

    index_x2i_Faxa_rain     = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_rain')
    if (index_x2i_Faxa_rain == 0) call shr_sys_abort('index_x2i_Faxa_rain is zero')

    index_x2i_Faxa_snow     = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_snow')
    if (index_x2i_Faxa_snow == 0) call shr_sys_abort('index_x2i_Faxa_snow is zero')

    index_x2i_Faxa_swndr    = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_swndr')
    if (index_x2i_Faxa_swndr == 0) call shr_sys_abort('index_x2i_Faxa_swndr is zero')

    index_x2i_Faxa_swvdr    = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_swvdr')
    if (index_x2i_Faxa_swvdr == 0) call shr_sys_abort('index_x2i_Faxa_swvdr is zero')

    index_x2i_Faxa_swndf    = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_swndf')
    if (index_x2i_Faxa_swndf == 0) call shr_sys_abort('index_x2i_Faxa_swndf is zero')

    index_x2i_Faxa_swvdf    = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_swvdf')
    if (index_x2i_Faxa_swvdf == 0) call shr_sys_abort('index_x2i_Faxa_swvdf is zero')

    index_x2i_Fioo_q        = shr_string_listGetIndexF(seq_flds_x2i_fields,'Fioo_q')
    if (index_x2i_Fioo_q == 0) call shr_sys_abort('index_x2i_Fioo_q is zero')

    nflds_x2i               = shr_string_listGetNum(seq_flds_x2i_fields)

    !-------------------------------------------------------------
    ! lnd -> drv 
    !-------------------------------------------------------------

    index_l2x_Sl_landfrac   = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_landfrac')
    if (index_l2x_Sl_landfrac == 0) call shr_sys_abort('index_l2x_Sl_landfrac is zero')

    index_l2x_Sl_t          = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_t')
    if (index_l2x_Sl_t == 0) call shr_sys_abort('index_l2x_Sl_t is zero')

    index_l2x_Sl_snowh      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_snowh')
    if (index_l2x_Sl_snowh == 0) call shr_sys_abort('index_l2x_Sl_snowh is zero')

    index_l2x_Sl_avsdr      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_avsdr')
    if (index_l2x_Sl_avsdr == 0) call shr_sys_abort('index_l2x_Sl_avsdr is zero')

    index_l2x_Sl_anidr      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_anidr')
    if (index_l2x_Sl_anidr == 0) call shr_sys_abort('index_l2x_Sl_anidr is zero')

    index_l2x_Sl_avsdf      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_avsdf')
    if (index_l2x_Sl_avsdf == 0) call shr_sys_abort('index_l2x_Sl_avsdf is zero')

    index_l2x_Sl_anidf      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_anidf')
    if (index_l2x_Sl_anidf == 0) call shr_sys_abort('index_l2x_Sl_anidf is zero')

    index_l2x_Sl_tref       = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_tref')
    if (index_l2x_Sl_tref == 0) call shr_sys_abort('index_l2x_Sl_tref is zero')

    index_l2x_Sl_qref       = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_qref')
    if (index_l2x_Sl_qref == 0) call shr_sys_abort('index_l2x_Sl_qref is zero')

    index_l2x_Fall_taux     = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_taux')
    if (index_l2x_Fall_taux == 0) call shr_sys_abort('index_l2x_Fall_taux is zero')

    index_l2x_Fall_tauy     = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_tauy')
    if (index_l2x_Fall_tauy == 0) call shr_sys_abort('index_l2x_Fall_tauy is zero')

    index_l2x_Fall_lat      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_lat')
    if (index_l2x_Fall_lat == 0) call shr_sys_abort('index_l2x_Fall_lat is zero')

    index_l2x_Fall_sen      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_sen')
    if (index_l2x_Fall_sen == 0) call shr_sys_abort('index_l2x_Fall_sen is zero')

    index_l2x_Fall_lwup     = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_lwup')
    if (index_l2x_Fall_lwup == 0) call shr_sys_abort('index_l2x_Fall_lwup is zero')

    index_l2x_Fall_evap     = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_evap')
    if (index_l2x_Fall_evap == 0) call shr_sys_abort('index_l2x_Fall_evap is zero')

    index_l2x_Fall_swnet    = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_swnet')
    if (index_l2x_Fall_swnet == 0) call shr_sys_abort('index_l2x_Fall_swnet is zero')

    index_l2x_Fall_nee      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_nee')

    index_l2x_Sl_ram1      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_ram1')
    index_l2x_Sl_fv        = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_fv')
    index_l2x_Fall_flxdst1 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxdst1')
    index_l2x_Fall_flxdst2 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxdst2')
    index_l2x_Fall_flxdst3 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxdst3')
    index_l2x_Fall_flxdst4 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxdst4')

    nflds_l2x              = shr_string_listGetNum(seq_flds_l2x_fields)

    !-------------------------------------------------------------
    ! runoff
    !-------------------------------------------------------------

    index_r2x_Forr_roff  = shr_string_listGetIndexF(seq_flds_r2x_fields,'Forr_roff')
    if (index_r2x_Forr_roff == 0) call shr_sys_abort('index_r2x_Forr_roff is zero')

    !-------------------------------------------------------------
    ! lnd -> drv
    !-------------------------------------------------------------

    index_x2l_Sa_z          = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_z')
    if (index_x2l_Sa_z == 0) call shr_sys_abort('index_x2l_Sa_z is zero')

    index_x2l_Sa_u          = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_u')
    if (index_x2l_Sa_u == 0) call shr_sys_abort('index_x2l_Sa_u is zero')

    index_x2l_Sa_v          = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_v')
    if (index_x2l_Sa_v == 0) call shr_sys_abort('index_x2l_Sa_v is zero')

    index_x2l_Sa_ptem       = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_ptem')
    if (index_x2l_Sa_ptem == 0) call shr_sys_abort('index_x2l_Sa_ptem is zero')

    index_x2l_Sa_pbot       = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_pbot')
    if (index_x2l_Sa_pbot == 0) call shr_sys_abort('index_x2l_Sa_pbot is zero')

    index_x2l_Sa_tbot       = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_tbot')
    if (index_x2l_Sa_tbot == 0) call shr_sys_abort('index_x2l_Sa_tbot is zero')

    index_x2l_Sa_shum       = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_shum')
    if (index_x2l_Sa_shum == 0) call shr_sys_abort('index_x2l_Sa_shum is zero')

    index_x2l_Faxa_lwdn     = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_lwdn')
    if (index_x2l_Faxa_lwdn == 0) call shr_sys_abort('index_x2l_Faxa_lwdn is zero')

    index_x2l_Faxa_rainc    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_rainc')
    if (index_x2l_Faxa_rainc == 0) call shr_sys_abort('index_x2l_Faxa_rainc is zero')

    index_x2l_Faxa_rainl    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_rainl')
    if (index_x2l_Faxa_rainl == 0) call shr_sys_abort('index_x2l_Faxa_rainl is zero')

    index_x2l_Faxa_snowc    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_snowc')
    if (index_x2l_Faxa_snowc == 0) call shr_sys_abort('index_x2l_Faxa_snowc is zero')

    index_x2l_Faxa_snowl    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_snowl')
    if (index_x2l_Faxa_snowl == 0) call shr_sys_abort('index_x2l_Faxa_snowl is zero')

    index_x2l_Faxa_swndr    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_swndr')
    if (index_x2l_Faxa_swndr == 0) call shr_sys_abort('index_x2l_Faxa_swndr is zero')

    index_x2l_Faxa_swvdr    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_swvdr')
    if (index_x2l_Faxa_swvdr == 0) call shr_sys_abort('index_x2l_Faxa_swvdr is zero')

    index_x2l_Faxa_swndf    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_swndf')
    if (index_x2l_Faxa_swndf == 0) call shr_sys_abort('index_x2l_Faxa_swndf is zero')

    index_x2l_Faxa_swvdf    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_swvdf')
    if (index_x2l_Faxa_swvdf == 0) call shr_sys_abort('index_x2l_Faxa_swvdf is zero')

    index_x2l_Sa_co2prog    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_co2prog')
    index_x2l_Sa_co2diag    = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_co2diag')

    nflds_x2l = shr_string_listGetNum(seq_flds_x2l_fields)

end subroutine seq_flds_indices_set

end module seq_flds_indices

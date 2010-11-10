module seq_flds_mod

   !----------------------------------------------------------------------------
   ! atm fields: atm->drv
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_a2x_states = &
         'Sa_z'        &    ! bottom atm level height         DEF
      //':Sa_u'        &    ! bottom atm level zon wind       DEF
      //':Sa_v'        &    ! bottom atm level mer wind       DEF
      //':Sa_tbot'     &    ! bottom atm level temp           DEF
      //':Sa_ptem'     &    ! bottom atm level pot temp       DEF
      //':Sa_shum'     &    ! bottom atm level spec hum       DEF
      //':Sa_dens'     &    ! bottom atm level air den        DEF
      //':Sa_pbot'     &    ! bottom atm level pressurea      DEF
      //':Sa_pslv'          ! sea level atm pressure          DEF

   ! Fluxes
   character(*), parameter :: seq_flds_a2x_fluxes = &
         'Faxa_lwdn'   &    ! downward lw heat flux           DEF
      //':Faxa_rainc'  &    ! prec: liquid "convective"       DEF
      //':Faxa_rainl'  &    ! prec: liquid "large scale"      DEF
      //':Faxa_snowc'  &    ! prec: frozen "convective"       DEF
      //':Faxa_snowl'  &    ! prec: frozen "large scale"      DEF
      //':Faxa_swndr'  &    ! sw: nir direct  downward        DEF
      //':Faxa_swvdr'  &    ! sw: vis direct  downward        DEF
      //':Faxa_swndf'  &    ! nir diffuse downward            DEF
      //':Faxa_swvdf'  &    ! sw: vis diffuse downward        DEF
      //':Faxa_swnet'       ! sw: net                         DEF


   character(*), parameter :: seq_flds_a2x_fields = &
      trim(seq_flds_a2x_states)//":"//trim(seq_flds_a2x_fluxes)

   !----------------------------------------------------------------------------
   ! atm fields: drv->atm
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   ! 
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_x2a_states = &
         'Sx_tref'     &    ! 2m reference temperature        DEF
      //':Sx_qref'     &    ! 2m reference specific humidity  DEF
      //':Sx_avsdr'    &    ! albedo, visible, direct         DEF
      //':Sx_anidr'    &    ! albedo, near-ir, direct         DEF
      //':Sx_avsdf'    &    ! albedo, visible, diffuse        DEF
      //':Sx_anidf'    &    ! albedo, near-ir, diffuse        DEF
      //':Sx_t'        &    ! surface temperature             DEF
      //':So_t'        &    ! sea surface temperature         DEF
      //':Sx_snowh'    &    ! surface snow depth              DEF
      //':Sx_lfrac'    &    ! surface land fraction           DEF
      //':Sx_ifrac'    &    ! surface ice fraction            DEF
      //':Sx_ofrac'    &    ! surface ocn fraction            DEF
      //':So_ustar'    &    ! needed for isoptope calc        DEF
      //':So_re'       &    ! needed for isoptope calc        DEF
      //':So_ssq'      &    ! needed for isoptope calc        DEF
      //':Sx_fv'       &    !
      //':Sx_ram1'          !

   ! Fluxes
   character(*), parameter :: seq_flds_x2a_fluxes = &
         'Faxx_taux'   &    ! wind stress, zonal              DEF
      //':Faxx_tauy'   &    ! wind stress, meridional         DEF
      //':Faxx_lat'    &    ! latent          heat flux       DEF
      //':Faxx_sen'    &    ! sensible        heat flux       DEF
      //':Faxx_lwup'   &    ! upward longwave heat flux       DEF
      //':Faxx_evap'   &    ! evaporation    water flux       DEF
      //':Faxx_flxdst1' &   ! dust flux bin 1
      //':Faxx_flxdst2' &   ! dust flux bin 2
      //':Faxx_flxdst3' &   ! dust flux bin 3
      //':Faxx_flxdst4'     ! dust flux bin 4

   character(*), parameter :: seq_flds_x2a_fields = &
      trim(seq_flds_x2a_states)//":"//trim(seq_flds_x2a_fluxes)


   !----------------------------------------------------------------------------
   ! ice fields: ice->drv
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   ! 
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_i2x_states = &
         'Si_t'        &    ! temperature                     DEF
      //':Si_tref'     &    ! 2m reference temperature        DEF
      //':Si_qref'     &    ! 2m reference specific humidity  DEF
      //':Si_ifrac'    &    ! fractional ice cov wrt ocean    DEF
      //':Si_avsdr'    &    ! albedo: visible, direct         DEF
      //':Si_anidr'    &    ! albedo: near ir, direct         DEF
      //':Si_avsdf'    &    ! albedo: visible, diffuse        DEF
      //':Si_anidf'    &    ! albedo: near ir, diffuse        DEF
      //':Si_sicthk'        ! sea ice thickness (m)           DEF (needed only for cam-som)

   ! Fluxes
   character(*), parameter :: seq_flds_i2x_fluxes = &
         'Faii_taux'   &    ! wind stress, zonal              DEF
      //':Faii_tauy'   &    ! wind stress, meridional         DEF
      //':Faii_lat'    &    ! latent          heat flux       DEF
      //':Faii_sen'    &    ! sensible        heat flux       DEF
      //':Faii_lwup'   &    ! upward longwave heat flux       DEF
      //':Faii_evap'   &    ! evaporation    water flux       DEF
      //':Faii_swnet'  &    ! shortwave: net absorbed         DEF
      //':Fioi_swpen'  &    ! net SW penetrating ice          DEF
      //':Fioi_melth'  &    ! heat  flux from melting ice     DEF
      //':Fioi_meltw'  &    ! water flux from melting ice     DEF
      //':Fioi_salt'   &    ! salt  flux from melting ice     DEF
      //':Fioi_taux'   &    ! ice/ocn stress, zonal           DEF
      //':Fioi_tauy'        ! ice/ocn stress, meridional      DEF

   character(*), parameter :: seq_flds_i2x_fields = &
      trim(seq_flds_i2x_states)//":"//trim(seq_flds_i2x_fluxes)


   !----------------------------------------------------------------------------
   ! ice fields: drv->ice
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_x2i_states = &
         'So_t'        &    ! ocean layer temperature         DEF
      //':Sa_z'        &    ! atm bottom layer height         DEF
      //':Sa_u'        &    ! atm u velocity                  DEF
      //':Sa_v'        &    ! atm v velocity                  DEF
      //':Sa_ptem'     &    ! atm potential temp              DEF
      //':Sa_tbot'     &    ! atm bottom temp                 DEF
      //':Sa_pbot'     &    ! atm bottom pressure             DEF
      //':Sa_shum'     &    ! atm specfic humidity            DEF
      //':Sa_dens'          ! atm air density                 DEF

   ! Fluxes
   character(*), parameter :: seq_flds_x2i_fluxes = &
         'Fioo_q'      &    ! ocn freeze or melt heat         DEF
      //':Faxa_swndr'  &    ! atm sw near-ir, direct          DEF
      //':Faxa_swvdr'  &    ! atm sw visable, direct          DEF
      //':Faxa_swndf'  &    ! atm sw near-ir, diffuse         DEF
      //':Faxa_swvdf'  &    ! atm sw visable, diffuse         DEF
      //':Faxa_lwdn'   &    ! long-wave down                  DEF
      //':Faxa_rain'   &    ! prec: liquid                    DEF
      //':Faxa_snow'        ! prec: frozen                    DEF

   character(*), parameter :: seq_flds_x2i_fields = &
      trim(seq_flds_x2i_states)//":"//trim(seq_flds_x2i_fluxes)

   !----------------------------------------------------------------------------
   ! lnd fields: lnd->drv 
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_l2x_states = &
         'Sl_t'        &    ! temperature                     DEF
      //':Sl_tref'     &    ! 2m reference temperature        DEF
      //':Sl_qref'     &    ! 2m reference specific humidity  DEF
      //':Sl_avsdr'    &    ! albedo: direct , visible        DEF
      //':Sl_anidr'    &    ! albedo: direct , near-ir        DEF
      //':Sl_avsdf'    &    ! albedo: diffuse, visible        DEF
      //':Sl_anidf'    &    ! albedo: diffuse, near-ir        DEF
      //':Sl_snowh'    &    ! snow height                     DEF
      //':Sl_landfrac' &    ! fractional land                 DEF
      //':Sl_fv'       &    ! friction velocity  
      //':Sl_ram1'          ! aerodynamical resistance

   ! Fluxes
   character(*), parameter :: seq_flds_l2x_fluxes = &
         'Fall_taux'   &    ! wind stress, zonal              DEF
      //':Fall_tauy'   &    ! wind stress, meridional         DEF
      //':Fall_lat'    &    ! latent          heat flux       DEF
      //':Fall_sen'    &    ! sensible        heat flux       DEF
      //':Fall_lwup'   &    ! upward longwave heat flux       DEF
      //':Fall_evap'   &    ! evaporation    water flux       DEF
      //':Fall_swnet'  &    ! shortwave: net absorbed         DEF
      //':Fall_flxdst1' &    ! dust flux bin 1
      //':Fall_flxdst2' &    ! dust flux bin 2
      //':Fall_flxdst3' &    ! dust flux bin 3
      //':Fall_flxdst4'      ! dust flux bin 4

   character(*), parameter :: seq_flds_l2x_fields = &
      trim(seq_flds_l2x_states)//":"//trim(seq_flds_l2x_fluxes)


   !----------------------------------------------------------------------------
   ! lnd fields: drv->lnd
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_x2l_states = &
         'Sa_z'        &    ! bottom atm level height         DEF
      //':Sa_u'        &    ! bottom atm level zon wind       DEF
      //':Sa_v'        &    ! bottom atm level mer wind       DEF
      //':Sa_tbot'     &    ! bottom atm level temp           DEF
      //':Sa_ptem'     &    ! bottom atm level pot temp       DEF
      //':Sa_shum'     &    ! bottom atm level spec hum       DEF
      //':Sa_pbot'          ! bottom atm level pressure       DEF

   ! Fluxes
   character(*), parameter :: seq_flds_x2l_fluxes = &
         'Faxa_lwdn'   &    ! downward longwave heat flux     DEF
      //':Faxa_rainc'  &    ! precip: liquid, convective      DEF
      //':Faxa_rainl'  &    ! precip: liquid, large-scale     DEF
      //':Faxa_snowc'  &    ! precip: frozen, convective      DEF
      //':Faxa_snowl'  &    ! precip: frozen, large-scale     DEF
      //':Faxa_swndr'  &    ! shortwave: nir direct  down     DEF
      //':Faxa_swvdr'  &    ! shortwave: vis direct  down     DEF
      //':Faxa_swndf'  &    ! shortwave: nir diffuse down     DEF
      //':Faxa_swvdf'       ! shortwave: vis diffuse down     DEF

   character(*), parameter :: seq_flds_x2l_fields = &
      trim(seq_flds_x2l_states)//":"//trim(seq_flds_x2l_fluxes)

   !----------------------------------------------------------------------------
   ! ocn fields: ocn->drv
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_o2x_states = &
         'So_t'        &    ! temperature                     DEF
      //':So_u'        &    ! velocity, zonal                 DEF
      //':So_v'        &    ! velocity, meridional            DEF
      //':So_s'        &    ! salinity                        DEF
      //':So_dhdx'     &    ! surface slope, zonal            DEF
      //':So_dhdy'          ! surface slope, meridional       DEF

   ! Fluxes
   character(*), parameter :: seq_flds_o2x_fluxes = &
         'Fioo_q'           ! heat of fusion (q>0) melt pot (q<0)  DEF 

   character(*), parameter :: seq_flds_o2x_fields = &
      trim(seq_flds_o2x_states)//":"//trim(seq_flds_o2x_fluxes)

   !----------------------------------------------------------------------------
   ! ocn fields: drv->ocn
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_x2o_states = &
         'Si_ifrac'    &    ! state: ice fraction wrt ocean   DEF
      //':Si_sicthk'   &    ! state: sea ice thickness (m)    DEF (needed only for cam-som)
      //':Sa_pslv'     &    ! state: sea level pressure       DEF
      //':Sx_duu10n'        ! state: 10m wind speed squared   DEF 

   ! Fluxes
   character(*), parameter :: seq_flds_x2o_fluxes = &
         'Foxx_taux'   &    ! wind stress: zonal              DEF
      //':Foxx_tauy'   &    ! wind stress: meridional         DEF
      //':Foxx_swnet'  &    ! heat flux: shortwave net        DEF
      //':Foxx_lat'    &    ! heat flux: latent               DEF
      //':Foxx_sen'    &    ! heat flux: sensible             DEF
      //':Foxx_lwdn'   &    ! heat flux: long-wave down       DEF
      //':Foxx_lwup'   &    ! heat flux: long-wave up         DEF
      //':Foxx_melth'  &    ! heat flux: melt                 DEF
      //':Foxx_salt'   &    ! salt flux                       DEF
      //':Foxx_prec'   &    ! water flux: rain+snow           DEF
      //':Foxx_snow'   &    ! water flux: snow                DEF
      //':Foxx_rain'   &    ! water flux: rain                DEF
      //':Foxx_evap'   &    ! water flux: evap                DEF
      //':Foxx_meltw'  &    ! water flux: melt                DEF
      //':Forr_roff'        ! water flux: runoff              DEF

   character(*), parameter :: seq_flds_x2o_fields = &
      trim(seq_flds_x2o_states)//":"//trim(seq_flds_x2o_fluxes)

   !----------------------------------------------------------------------------
   ! hub computed fields: atm/ocn states/fluxes
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_xao_states = &
         'So_tref'     &    ! 2m reference temperature        DEF 
      //':So_qref'     &    ! 2m reference specific humidity  DEF 
      //':So_avsdr'    &    ! albedo: visible, direct         DEF 
      //':So_anidr'    &    ! albedo: near ir, direct         DEF 
      //':So_avsdf'    &    ! albedo: visible, diffuse        DEF 
      //':So_anidf'    &    ! albedo: near ir, diffuse        DEF 
      //':Sx_duu10n'   &    ! diag 10m wind speed squared     DEF
      //':So_ustar'    &    ! ustar                           DEF
      //':So_ssq'      &    ! surface saturation spec. hum.   DEF
      //':So_re'            ! sqrt of exch. coeff (tracers)   DEF
   
   ! Fluxes
   character(*), parameter :: seq_flds_xao_fluxes = &
         'Faox_taux'   &    ! wind stress, zonal              DEF 
      //':Faox_tauy'   &    ! wind stress, meridional         DEF 
      //':Faox_lat'    &    ! latent          heat flux       DEF 
      //':Faox_sen'    &    ! sensible        heat flux       DEF 
      //':Faox_evap'   &    ! evaporation    water flux       DEF 
      //':Faox_lwup'   &    ! upward longwave heat flux       DEF 
      //':Faox_swnet'       ! shortwave: net absorbed         DEF 

   character(*), parameter :: seq_flds_xao_fields = &
      trim(seq_flds_xao_states)//":"//trim(seq_flds_xao_fluxes)

   !----------------------------------------------------------------------------
   ! run-off field 
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: seq_flds_r2x_states = &
         ''

   ! Fluxes
   character(*), parameter :: seq_flds_r2x_fluxes = &
         'Forr_roff'        ! runoff to ocean                 DEF

   !  character(*), parameter :: seq_flds_r2x_fields = &
   !     trim(seq_flds_r2x_states)//":"//trim(seq_flds_r2x_fluxes)

   character(*), parameter :: seq_flds_r2x_fields = &
      trim(seq_flds_r2x_fluxes)

   !----------------------------------------------------------------------------
   ! component names
   !----------------------------------------------------------------------------

   character(32),parameter :: seq_flds_atmname='atm'
   character(32),parameter :: seq_flds_ocnname='ocn'
   character(32),parameter :: seq_flds_icename='ice'
   character(32),parameter :: seq_flds_lndname='lnd'
   character(32),parameter :: seq_flds_rtmname='roff'
   character(32),parameter :: seq_flds_cplname='cpl'

end module seq_flds_mod


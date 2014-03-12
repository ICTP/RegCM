Module mod_clm_drydepvelocity                                              

  !-----------------------------------------------------------------------  
  !  
  ! Purpose:  
  ! Deposition velocity (m/s) 
  !  
  ! Method:  
  ! This code simulates dry deposition velocities using the Wesely scheme.   
  ! Details of this method can be found in:  
  ! 
  ! M.L Wesely. Parameterization of surface resistances to gaseous dry deposition 
  ! in regional-scale numericl models. 1989. Atmospheric Environment vol.23 No.6  
  ! pp. 1293-1304. 
  ! 
  ! In Wesely (1998) "the magnitude of the dry deposition velocity can be found 
  ! as: 
  ! 
  !  |vd|=(ra+rb+rc)^-1 
  ! 
  ! where ra is the aerodynamic resistance (common to all gases) between a  
  ! specific height and the surface, rb is the quasilaminar sublayer resistance 
  ! (whose only dependence on the porperties of the gas of interest is its 
  ! molecular diffusivity in air), and rc is the bulk surface resistance". 
  ! 
  ! In this subroutine both ra and rb are calculated elsewhere in CLM.  Thus ra  
  ! and rb were "globalized" in order to gain access to them for the calculation. 
  ! "ram1" is the CLM variable used for ra.  ram1 was globalized in the following 
  ! subroutines; BareGroundFluxes.F90, Biogeophysics_lake.F90, CanopyFluxes.F90, 
  ! and clmtype.F90.  
  ! 
  ! "rb" is the CLM variable used for rb in the Wesely equation above.  rb was  
  ! globalized in the following subroutines; clmtype.F90     
  ! 
  ! In Wesely (1989) rc is estimated for five seasonal categories and 11 landuse 
  ! types.  For each season and landuse type, Wesely compiled data into a  
  ! look-up-table for several parameters used to calculate rc. In this subroutine 
  ! the same values are used as found in wesely's look-up-tables, the only  
  ! difference is that this subroutine uses a CLM generated LAI to select values 
  ! from the look-up-table instead of seasonality.  Inaddition, Wesely(1989)  
  ! land use types are "mapped" into CLM plant function types (PFT). 
  ! 
  ! Subroutine written to operate at the patch level. 
  ! 
  ! Output: 
  ! 
  ! vd(n_species) !Dry deposition velocity [m s-1] for each molecule or species 
  !  
  ! Author: Beth Holland and  James Sulzman 
  ! 
  ! Modified: Francis Vitt -- 30 Mar 2007
  !----------------------------------------------------------------------- 

  use mod_realkinds
  use mod_mpmessage
  use mod_stdio
  use mod_clm_type 
  use mod_clm_atmlnd,           only : clm_a2l
  use mod_clm_drydep,       only : n_drydep, drydep_list
  use mod_clm_drydep,       only : drydep_method, DD_XLND
  use mod_clm_drydep,       only : index_o3=>o3_ndx, index_o3a=>o3a_ndx, index_so2=>so2_ndx, index_h2=>h2_ndx, &
                                   index_co=>co_ndx, index_ch4=>ch4_ndx, index_pan=>pan_ndx, &
                                   index_xpan=>xpan_ndx
  implicit none 
  save 

  private

  save

  public :: depvel_compute

CONTAINS 

  !----------------------------------------------------------------------- 
  ! computes the dry deposition velocity of tracers
  !----------------------------------------------------------------------- 
  subroutine depvel_compute( lbp , ubp ) 
    use mod_clm_varcon     , only :  tmelt => tfrz
    use mod_clm_drydep    , only :  seq_drydep_setHCoeff, mapping, drat, foxd, &
                                    rcls, h2_a, h2_b, h2_c, ri, rac, rclo, rlu, &
                                    rgss, rgso
    use mod_clm_varcon        , only : istsoil, istice, istslak, istdlak, &
            istwet, isturb
    use mod_clm_pftvarcon         , only : noveg, ndllf_evr_tmp_tree, ndllf_evr_brl_tree,   &
                                   ndllf_dcd_brl_tree,        nbrdlf_evr_trp_tree,  &
                                   nbrdlf_evr_tmp_tree,       nbrdlf_dcd_trp_tree,  &
                                   nbrdlf_dcd_tmp_tree,       nbrdlf_dcd_brl_tree,  &
                                   nbrdlf_evr_shrub,          nbrdlf_dcd_tmp_shrub, &
                                   nbrdlf_dcd_brl_shrub,      nc3_arctic_grass,     &
                                   nc3_nonarctic_grass,       nc4_grass, nc3crop,   &
                                   nc3irrig,       npcropmin, npcropmax

    implicit none 

    !----Arguments-----------------------------------------------------

    integer, intent(in) :: lbp, ubp                    ! pft bounds

    ! ------------------------ local variables ------------------------ 
    ! local pointers to implicit in arguments 
    logical , pointer :: pactive(:)       ! true=>do computations on this pft (see reweightMod for details)
    integer , pointer :: plandunit(:)     !pft's landunit index
    integer , pointer :: ivt(:)           !landunit type
    integer , pointer :: pgridcell(:)     !pft's gridcell index
    real(rk8), pointer :: elai(:)          !one-sided leaf area index with burying by snow 
    real(rk8), pointer :: forc_t(:)        !atmospheric temperature (Kelvin) 
    real(rk8), pointer :: forc_q(:)        !atmospheric specific humidity (kg/kg) 
    real(rk8), pointer :: forc_psrf(:)     !surface pressure (Pa) 
    real(rk8), pointer :: latdeg(:)        !latitude (degrees) 
    real(rk8), pointer :: londeg(:)        !longitude (degrees) 
    real(rk8), pointer :: forc_rain(:)     !rain rate [mm/s] 
    real(rk8), pointer :: forc_solad(:,:)  !direct beam radiation (visible only) 
    real(rk8), pointer :: forc_solai(:,:)  !direct beam radiation (visible only) 
    real(rk8), pointer :: ram1(:)          !aerodynamical resistance 
    real(rk8), pointer :: vds(:)           !aerodynamical resistance 
    real(rk8), pointer :: rssun(:)         !stomatal resistance 
    real(rk8), pointer :: rssha(:)         !shaded stomatal resistance (s/m) 
    real(rk8), pointer :: fsun(:)          !sunlit fraction of canopy 
    real(rk8), pointer :: rb1(:)           !leaf boundary layer resistance [s/m]
    real(rk8), pointer :: annlai(:,:)      !12 months of monthly lai from input data set 
    real(rk8), pointer :: mlaidiff(:)      !difference in lai between month one and month two 
    real(rk8), pointer :: velocity(:,:)
    real(rk8), pointer :: snow_depth(:)        ! snow height (m)

    integer, pointer :: pcolumn(:)        ! column index associated with each pft
    integer :: c
    integer , pointer :: itypelun(:)      ! landunit type

    real(rk8), pointer :: h2osoi_vol(:,:)    ! volumetric soil water (0<=h2osoi_vol<=watsat)
    real(rk8) :: soilw, var_soilw, fact_h2, dv_soil_h2

    ! new local variables  
    integer :: pi,g, l
    integer :: ispec 
    integer :: length 
    integer :: wesveg              !wesely vegegation index  
    integer :: clmveg              !clm veg index from ivegtype 
    integer :: i 
    integer :: index_season        !seasonal index based on LAI.  This indexs wesely data tables 
    integer :: nstep               !current step 
    integer :: indexp 

    real(rk8) :: pg          ! surface pressure 
    real(rk8) :: tc          ! temperature in celsius  
    real(rk8) :: rs          ! constant for calculating rsmx  
    real(rk8) :: es          ! saturation vapor pressur 
    real(rk8) :: ws          ! saturation mixing ratio 
    real(rk8) :: rmx         ! resistance by vegetation 
    real(rk8) :: qs          ! saturation specific humidity 
    real(rk8) :: dewm        ! multiplier for rs when dew occurs 
    real(rk8) :: crs         ! multiplier to calculate crs 
    real(rk8) :: rdc         ! part of lower canopy resistance 
    real(rk8) :: rain        ! rain fall 
    real(rk8) :: spec_hum    ! specific humidity 
    real(rk8) :: solar_flux  ! solar radiation(direct beam) W/m2 
    real(rk8) :: lat         ! latitude in degrees 
    real(rk8) :: lon         ! longitude in degrees 
    real(rk8) :: sfc_temp    ! surface temp 
    real(rk8) :: minlai      ! minimum of monthly lai 
    real(rk8) :: maxlai      ! maximum of monthly lai 
    real(rk8) :: rds

    logical  :: has_dew
    logical  :: has_rain
    real(rk8), parameter :: rain_threshold      = 1.D-7  ! of the order of 1cm/day expressed in m/s

    ! local arrays: dependent on species only 
    ! 

    real(rk8), dimension(n_drydep) :: rsmx !vegetative resistance (plant mesophyll) 
    real(rk8), dimension(n_drydep) :: rclx  !lower canopy resistance  
    real(rk8), dimension(n_drydep) :: rlux  !vegetative resistance (upper canopy) 
    real(rk8), dimension(n_drydep) :: rgsx  !gournd resistance 
    real(rk8), dimension(n_drydep) :: heff   
    real(rk8) :: rc    !combined surface resistance 
    real(rk8) :: cts   !correction to flu rcl and rgs for frost 
    real(rk8) :: rlux_o3

    ! constants 
    real(rk8), parameter :: slope = 0.D0      ! Used to calculate  rdc in (lower canopy resistance) 
    integer, parameter :: wveg_unset = -1     ! Unset Wesley vegetation type

    character(len=32), parameter :: subname = "depvel_compute"

    !-------------------------------------------------------------------------------------
    ! jfl : mods for PAN
    !-------------------------------------------------------------------------------------
    real(rk8)                  :: dv_pan
    real(rk8) :: c0_pan(11) = (/ 0.000D0, 0.006D0, 0.002D0, 0.009D0, 0.015D0, &
                                0.006D0, 0.000D0, 0.000D0, 0.000D0, 0.002D0, 0.002D0 /)
    real(rk8) :: k_pan (11) = (/ 0.000D0, 0.010D0, 0.005D0, 0.004D0, 0.003D0, &
                                0.005D0, 0.000D0, 0.000D0, 0.000D0, 0.075D0, 0.002D0 /)
    !----------------------------------------------------------------------- 
    if ( n_drydep == 0 .or. drydep_method /= DD_XLND ) return

    ! local pointers to original implicit out arrays 

    ! Assign local pointers to derived subtypes components (column-level) 
    forc_t     => clm_a2l%forc_t
    forc_q     => clm_a2l%forc_q
    forc_psrf  => clm_a2l%forc_pbot
    forc_rain  => clm_a2l%forc_rain 

    latdeg     => clm3%g%latdeg
    londeg     => clm3%g%londeg
    pactive    => clm3%g%l%c%p%active
    ivt        => clm3%g%l%c%p%itype
    elai       => clm3%g%l%c%p%pps%elai 
    ram1       => clm3%g%l%c%p%pps%ram1 
    vds        => clm3%g%l%c%p%pps%vds
    fsun       => clm3%g%l%c%p%pps%fsun 
    rssun      => clm3%g%l%c%p%pps%rssun 
    rssha      => clm3%g%l%c%p%pps%rssha 
    rb1        => clm3%g%l%c%p%pps%rb1 
    mlaidiff   => clm3%g%l%c%p%pps%mlaidiff 
    annlai     => clm3%g%l%c%p%pps%annlai    

    forc_solai => clm_a2l%forc_solai 
    forc_solad => clm_a2l%forc_solad

    pgridcell  => clm3%g%l%c%p%gridcell
    plandunit  => clm3%g%l%c%p%landunit

    pcolumn    => clm3%g%l%c%p%column
    itypelun   => clm3%g%l%itype

    h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol

    velocity   => clm3%g%l%c%p%pdd%drydepvel ! cm/sec

    snow_depth        => clm3%g%l%c%cps%snow_depth

    ! Assign local pointers to original implicit out arrays 
    !_________________________________________________________________ 
    ! 
    ! Begin loop through pfts 
    pft_loop: do pi = lbp,ubp
       l = plandunit(pi)

       active: if (pactive(pi)) then

          c = pcolumn(pi)
          g = pgridcell(pi)
          !solar_flux = forc_lwrad  !rename CLM variables to fit with Dry Dep variables 

          pg         = forc_psrf(g)  
          spec_hum   = forc_q(g)
          rain       = forc_rain(g) 
          sfc_temp   = forc_t(g) 
          lat        = latdeg(g) 
          lon        = londeg(g) 
          solar_flux = forc_solad(g,1) 
          clmveg     = ivt(pi) 
          soilw = h2osoi_vol(c,1)

          !  print *,'bb',pi,cps%npfts,lat,lon,clmveg 
          !map CLM veg type into Wesely veg type  
          wesveg = wveg_unset 
          if (clmveg == noveg                               ) wesveg = 8 
          if (clmveg == ndllf_evr_tmp_tree                  ) wesveg = 5 
          if (clmveg == ndllf_evr_brl_tree                  ) wesveg = 5 
          if (clmveg == ndllf_dcd_brl_tree                  ) wesveg = 5 
          if (clmveg == nbrdlf_evr_trp_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_evr_tmp_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_dcd_trp_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_dcd_tmp_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_dcd_brl_tree                 ) wesveg = 4 
          if (clmveg == nbrdlf_evr_shrub                    ) wesveg = 11 
          if (clmveg == nbrdlf_dcd_tmp_shrub                ) wesveg = 11 
          if (clmveg == nbrdlf_dcd_brl_shrub                ) wesveg = 11 
          if (clmveg == nc3_arctic_grass                    ) wesveg = 3 
          if (clmveg == nc3_nonarctic_grass                 ) wesveg = 3 
          if (clmveg == nc4_grass                           ) wesveg = 3 
          if (clmveg == nc3crop                             ) wesveg = 2 
          if (clmveg == nc3irrig                            ) wesveg = 2 
          if (clmveg >= npcropmin .and. clmveg <= npcropmax ) wesveg = 2 
          if (wesveg == wveg_unset )then
             write(stderr,*) 'clmveg = ', clmveg, 'itypelun = ', itypelun(l)
             call fatal(__FILE__,__LINE__, &
               subname//': Not able to determine Wesley vegetation type')
          end if

          ! creat seasonality index used to index wesely data tables from LAI,  Bascially 
          !if elai is between max lai from input data and half that max the index_season=1 


          !mail1j and mlai2j are the two monthly lai values pulled from a CLM input data set 
          !/fs/cgd/csm/inputdata/lnd/clm2/rawdata/mksrf_lai.nc.  lai for dates in the middle  
          !of the month are interpolated using using these values and stored in the variable  
          !elai (done elsewhere).  If the difference between mlai1j and mlai2j is greater 
          !than zero it is assumed to be fall and less than zero it is assumed to be spring. 

          !wesely seasonal "index_season" 
          ! 1 - midsummer with lush vegetation 
          ! 2 - Autumn with unharvested cropland 
          ! 3 - Late autumn after frost, no snow 
          ! 4 - Winter, snow on ground and subfreezing 
          ! 5 - Transitional spring with partially green short annuals 


          !mlaidiff=jan-feb 
          minlai=minval(annlai(:,pi)) 
          maxlai=maxval(annlai(:,pi)) 

          index_season = -1

          if ( itypelun(l) /= istsoil )then
             if ( itypelun(l) == istice ) then
                wesveg       = 8
                index_season = 4
             elseif ( itypelun(l) == istdlak .or. itypelun(l) == istslak ) then
                wesveg       = 7
                index_season = 4
             elseif ( itypelun(l) == istwet ) then
                wesveg       = 9
                index_season = 2
             elseif ( itypelun(l) == isturb ) then
                wesveg       = 1
                index_season = 2
             end if
          else if ( snow_depth(c) > 0 ) then
             index_season = 4
          else if(elai(pi).gt.0.5D0*maxlai) then  
             index_season = 1  
          endif

          if (index_season<0) then 
             if (elai(pi).lt.(minlai+0.05*(maxlai-minlai))) then  
                index_season = 3 
             endif
          endif

          if (index_season<0) then 
             if (mlaidiff(pi).gt.0.0D0) then 
                index_season = 2 
             elseif (mlaidiff(pi).lt.0.0D0) then 
                index_season = 5 
             elseif (mlaidiff(pi).eq.0.0D0) then 
                index_season = 3 
             endif
          endif

          if (index_season<0) then 
             call fatal(__FILE__,__LINE__,&
               subname//': not able to determine season')
          endif

          ! saturation specific humidity 
          ! 
          es = 611D0*exp(5414.77D0*((1.D0/tmelt)-(1.D0/sfc_temp))) 
          ws = .622D0*es/(pg-es) 
          qs = ws/(1.D0+ws) 

          has_dew = .false.
          if( qs <= spec_hum ) then
             has_dew = .true.
          end if
          if( sfc_temp < tmelt ) then
             has_dew = .false.
          end if

          has_rain = rain > rain_threshold

          if ( has_dew .or. has_rain ) then
             dewm = 3.D0
          else
             dewm = 1.D0
          end if

          ! 
          ! constant in determining rs 
          ! 
          crs = 1.D36

          tc = sfc_temp - tmelt 
          if(sfc_temp.gt.tmelt.and.sfc_temp.lt.313.15D0) then 
             crs = (1.D0+(200.D0/(solar_flux+.1D0))**2) * (400.D0/(tc*(40.D0-tc))) 
          endif
          ! 
          ! rdc (lower canopy res) 
          ! 
          rdc=100.D0*(1.D0+1000.D0/(solar_flux+10.D0))/(1.D0+1000.D0*slope) 

          ! surface resistance : depends on both land type and species 
          ! land types are computed seperately, then resistance is computed as average of values 
          ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1 
          ! 
          ! compute rsmx = 1/(rs+rm) : multiply by 3 if surface is wet 
          ! 

          !******************************************************* 
          call seq_drydep_setHCoeff( sfc_temp, heff(:n_drydep) )
          !********************************************************* 

          species_loop1: do ispec=1, n_drydep
             if(mapping(ispec).le.0) cycle 

             if(ispec.eq.index_o3.or.ispec.eq.index_o3a.or.ispec.eq.index_so2) then 
                rmx=0.D0
             else 
                rmx=1.D0/((heff(ispec)/3000.D0)+(100.D0*foxd(ispec))) 
             endif

             ! correction for frost 
             cts = 1000.D0*exp( -tc - 4.D0 )                 ! correction for frost
             rgsx(ispec) = cts + 1.D0/((heff(ispec)/(1.D5*rgss(index_season,wesveg))) + & 
                                        (foxd(ispec)/rgso(index_season,wesveg))) 

             !-------------------------------------------------------------------------------------
             ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
             !-------------------------------------------------------------------------------------
             if( ispec == index_h2 .or. ispec == index_co .or. ispec == index_ch4 ) then

                if( ispec == index_co ) then
                   fact_h2 = 1.0D0
                elseif ( ispec == index_h2 ) then
                   fact_h2 = 0.5D0
                elseif ( ispec == index_ch4 ) then
                   fact_h2 = 50.0D0
                end if
                !-------------------------------------------------------------------------------------
                ! no deposition on snow, ice, desert, and water
                !-------------------------------------------------------------------------------------
                if( wesveg == 1 .or. wesveg == 7 .or. wesveg == 8 .or. index_season == 4 ) then
                   rgsx(ispec) = 1.D36
                else
                   var_soilw = max( .1D0,min( soilw,.3D0 ) )
                   if( wesveg == 3 ) then
                      var_soilw = log( var_soilw )
                   end if
                   dv_soil_h2 = h2_c(wesveg) + var_soilw*(h2_b(wesveg) + var_soilw*h2_a(wesveg))
                   if( dv_soil_h2 > 0.D0 ) then
                      rgsx(ispec) = fact_h2/(dv_soil_h2*1.D-4)
                   end if
                end if
             end if

             rs=ri(index_season,wesveg)*crs 
             if(wesveg.eq.7) then ! over water
                rclx(ispec)=1.D36
                rsmx(ispec)=1.D36
                rlux(ispec)=1.D36
             else 
                ! ??? fvitt   rs=(fsun(pi)*rssun(pi))+(rssha(pi)*(1.-fsun(pi))) -- made the same as mo_drydep
                rsmx(ispec) = (dewm*rs*drat(ispec)+rmx) 
             endif

             !-------------------------------------------------------------------------------------
             ! jfl : special case for PAN
             !-------------------------------------------------------------------------------------
             if( ispec == index_pan .or. ispec == index_xpan ) then
                dv_pan =  c0_pan(wesveg) * (1.D0 - exp( -k_pan(wesveg)*(dewm*rs*drat(ispec))*1.D-2 ))
                if( dv_pan > 0.D0 .and. index_season /= 4 ) then
                   rsmx(ispec) = ( 1.D0/dv_pan )
                end if
             end if

             rclx(ispec) = cts + 1.D0/((heff(ispec)/(1.D5*rcls(index_season,wesveg))) + & 
                           (foxd(ispec)/rclo(index_season,wesveg))) 
             rlux(ispec) = cts + rlu(index_season,wesveg)/(1.D-5*heff(ispec)+foxd(ispec)) 

          end do species_loop1
          
          ! 
          ! no effect over water
          ! 
          no_water: if(wesveg.ne.7) then 
             ! 
             ! no effect if sfc_temp < O C 
             ! 
             non_freezing: if(sfc_temp.gt.tmelt) then

                if( has_dew ) then 
                   rlux_o3 = 1.D0/((1.D0/3000.D0)+(1.D0/(3.D0*rlu(index_season,wesveg)))) 
                   if (index_o3 > 0) then
                      rlux(index_o3) = rlux_o3
                   endif
                   if (index_o3a > 0) then
                      rlux(index_o3a) = rlux_o3
                   endif
                endif

                if(has_rain) then 
                   rlux_o3 = 1.D0/((1.D0/1000.D0)+(1.D0/(3.D0*rlu(index_season,wesveg)))) 
                   if (index_o3 > 0) then
                      rlux(index_o3) = rlux_o3
                   endif
                   if (index_o3a > 0) then
                      rlux(index_o3a) = rlux_o3
                   endif
                endif

                if ( index_o3 > 0 ) then
                   rclx(index_o3) = cts + rclo(index_season,wesveg)
                   rlux(index_o3) = cts + rlux(index_o3)
                end if
                if ( index_o3a > 0 ) then
                   rclx(index_o3a) = cts + rclo(index_season,wesveg)
                   rlux(index_o3a) = cts + rlux(index_o3a)
                end if

                species_loop2: do ispec=1,n_drydep 
                   if(mapping(ispec).le.0) cycle 
                   if(ispec.ne.index_o3.and.ispec.ne.index_o3a.and.ispec.ne.index_so2) then 

                      if( has_dew ) then
                         rlux(ispec)=1.D0/((1.D0/(3.D0*rlux(ispec)))+ & 
                              (1.D-7*heff(ispec))+(foxd(ispec)/rlux_o3)) 
                      endif

                   elseif(ispec.eq.index_so2) then 

                      if( has_dew ) then
                         rlux(ispec) = 100.D0
                      endif

                      if(has_rain) then 
                         rlux(ispec) = 1.D0/((1.D0/5000.D0)+(1.D0/(3.D0*rlu(index_season,wesveg)))) 
                      endif

                      rclx(ispec) = cts + rcls(index_season,wesveg)
                      rlux(ispec) = cts + rlux(ispec)

                      if( has_dew .or. has_rain ) then
                         rlux(ispec)=50.D0
                      endif

                   endif
                end do species_loop2

             endif non_freezing

          endif no_water

          rds = 1.D0/vds(pi)

          species_loop3: do ispec=1,n_drydep 
             if(mapping(ispec).le.0) cycle 

             ! 
             ! compute rc 
             ! 
             rc = 1.D0/((1.D0/rsmx(ispec))+(1.D0/rlux(ispec)) + & 
                        (1.D0/(rdc+rclx(ispec)))+(1.D0/(rac(index_season,wesveg)+rgsx(ispec)))) 
             rc = max( 10.D0, rc)
             !
             ! assume no surface resistance for SO2 over water
             !
             if ( drydep_list(ispec) == 'SO2' .and. wesveg == 7 ) then
               rc = 0.D0
             end if

             select case( drydep_list(ispec) )
             case ( 'SO4' )
                velocity(pi,ispec) = (1.D0/(ram1(pi)+rds))*100.D0
             case ( 'NH4','NH4NO3','XNH4NO3' )
                velocity(pi,ispec) = (1.D0/(ram1(pi)+0.5D0*rds))*100.D0
             case ( 'Pb' )
                velocity(pi,ispec) = 0.2D0
             case ( 'CB1', 'CB2', 'OC1', 'OC2', 'SOAM', 'SOAI', 'SOAT', 'SOAB', 'SOAX' )
                velocity(pi,ispec) = 0.10D0
             case default
                velocity(pi,ispec) = (1.D0/(ram1(pi)+rb1(pi)+rc))*100.D0
             end select
          end do species_loop3

       endif active

    end do pft_loop

  end subroutine depvel_compute

end Module mod_clm_drydepvelocity                                              

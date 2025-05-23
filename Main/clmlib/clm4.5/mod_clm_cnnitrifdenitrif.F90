module mod_clm_cnnitrifdenitrif
#ifdef CN
#ifdef NITRIF_DENITRIF
  !
  ! Calculate nitrification and denitrification rates
  !
  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_clm_varcon, only: secspday, tfrz

  implicit none

  save

  private

  public:: nitrif_denitrif
  ! stop nitrification and denitrification in frozen soils
  logical, public :: no_frozen_nitrif_denitrif = .false.

  contains
  !
  !  calculate nitrification and denitrification rates
  !
  subroutine nitrif_denitrif(lbc, ubc, num_soilc, filter_soilc)
    use mod_clm_type
    use mod_clm_varpar     , only : nlevgrnd,nlevdecomp
    use mod_clm_varcon, only: rpi, denh2o, dzsoi, zisoi, grav
#ifdef LCH4
    use mod_clm_varcon, only: d_con_g, d_con_w
    use mod_clm_ch4varcon,        only : organic_max
#ifdef CENTURY_DECOMP
    use mod_clm_cndecompcascadecentury, only: anoxia_wtsat
#else
    use mod_clm_cndecompcascadebgc, only: anoxia_wtsat
#endif
#endif
    use mod_clm_varcon, only : spval
    implicit none
    integer(ik4), intent(in) :: lbc, ubc   ! column bounds
    integer(ik4), intent(in) :: num_soilc  ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns

    integer(ik4) :: c, fc, reflev, j

    ! column level
    ! total soil respiration rate (g C / m3 / s)
    real(rk8) :: soil_hr_vr(lbc:ubc,1:nlevdecomp)
    real(rk8) :: g_per_m3__to__ug_per_gsoil
    real(rk8) :: g_per_m3_sec__to__ug_per_gsoil_day

    ! put on clm structure for diagnostic purposes
    ! (ugN / g soil) soil nitrate concentration
    real(rk8), pointer, contiguous :: smin_no3_massdens_vr(:,:)
    real(rk8), pointer, contiguous :: k_nitr_t_vr(:,:)
    real(rk8), pointer, contiguous :: k_nitr_ph_vr(:,:)
    real(rk8), pointer, contiguous :: k_nitr_h2o_vr(:,:)
    real(rk8), pointer, contiguous :: k_nitr_vr(:,:)
    real(rk8), pointer, contiguous :: wfps_vr(:,:)
    real(rk8), pointer, contiguous :: fmax_denit_carbonsubstrate_vr(:,:)
    real(rk8), pointer, contiguous :: fmax_denit_nitrate_vr(:,:)
    real(rk8), pointer, contiguous :: f_denit_base_vr(:,:)

    real(rk8) :: k_nitr_max ! maximum nitrification rate constant (1/s)
    real(rk8) :: mu, sigma
    real(rk8) :: t
#ifdef LCH4
    real(rk8), pointer, contiguous :: pH(:)
#else
    real(rk8) :: pH(lbc:ubc)
#endif
    real(rk8), pointer, contiguous :: phr_vr(:,:)   ! potential hr (not N-limited)
    real(rk8), pointer, contiguous :: w_scalar(:,:) ! soil water scalar for decomp
    real(rk8), pointer, contiguous :: t_scalar(:,:) ! temperature scalar for decomp
    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
    real(rk8), pointer, contiguous :: h2osoi_vol(:,:)
    ! liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous :: h2osoi_liq(:,:)
    ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(rk8), pointer, contiguous :: watsat(:,:)
    ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous :: t_soisno(:,:)
    ! (gN/m3) soil mineral NH4 pool
    real(rk8), pointer, contiguous :: smin_nh4_vr(:,:)
    real(rk8), pointer, contiguous :: smin_no3_vr(:,:)   ! (gN/m3) soil mineral NO3 pool
    ! bulk density of dry soil material [kg/m3]
    real(rk8), pointer, contiguous :: bd(:,:)
    ! layer thickness (m)  (-nlevsno+1:nlevgrnd)
    real(rk8), pointer, contiguous :: dz(:,:)
    ! maximumn monthly-mean soil temperature
    real(rk8), pointer, contiguous :: tmean_monthly_max_vr(:,:)
    ! monthly-mean soil temperature
    real(rk8), pointer, contiguous :: tmean_monthly_vr(:,:)
    ! (gN/m3/s) potential soil nitrification flux
    real(rk8), pointer, contiguous :: pot_f_nit_vr(:,:)
    ! (gN/m3/s) potential soil denitrification flux
    real(rk8), pointer, contiguous :: pot_f_denit_vr(:,:)
    ! volumetric soil water at field capacity (nlevsoi)
    real(rk8), pointer, contiguous :: watfc(:,:)
    ! Clapp and Hornberger "b" (nlevgrnd)
    real(rk8), pointer, contiguous :: bsw(:,:)
    ! ratio of N2 to N2O production by denitrification [gN/gN]
    real(rk8), pointer, contiguous :: n2_n2o_ratio_denit_vr(:,:)

    !debug-- put these in clmtype for outing to hist files
    !diffusivity (unitless fraction of total diffusivity)
    real(rk8), pointer, contiguous :: diffus(:,:)
    real(rk8), pointer, contiguous :: ratio_k1(:,:)
    real(rk8), pointer, contiguous :: ratio_no3_co2(:,:)
    real(rk8), pointer, contiguous :: soil_co2_prod(:,:)  ! (ug C / g soil / day)
    real(rk8), pointer, contiguous :: fr_WFPS(:,:)
    ! (kg soil / m3) bulk density of soil (including water)
    real(rk8), pointer, contiguous :: soil_bulkdensity(:,:)

    real(rk8) :: co2diff_con(2)  ! diffusion constants for CO2
    real(rk8) :: eps
    real(rk8) :: f_a

    ! (J/m^2), Arah and Vinten 1995
    real(rk8) :: surface_tension_water = 73.e-3_rk8

    real(rk8) :: rij_kro_a = 1.5e-10_rk8       !  Arah and Vinten 1995
    real(rk8) :: rij_kro_alpha = 1.26_rk8      !  Arah and Vinten 1995
    real(rk8) :: rij_kro_beta = 0.6_rk8        !  Arah and Vinten 1995
    real(rk8) :: rij_kro_gamma = 0.6_rk8       !  Arah and Vinten 1995
    real(rk8) :: rij_kro_delta = 0.85_rk8      !  Arah and Vinten 1995

    real(rk8) ::  rho_w  = 1.e3_rk8              ! (kg/m3)
    real(rk8) :: r_max
    real(rk8) :: r_min(lbc:ubc,1:nlevdecomp)
    real(rk8) :: ratio_diffusivity_water_gas(lbc:ubc,1:nlevdecomp)
    real(rk8) :: om_frac
    ! scalar values in sat portion for averaging
    real(rk8) :: anaerobic_frac_sat, r_psi_sat, r_min_sat

    real(rk8), pointer, contiguous :: r_psi(:,:)
    real(rk8), pointer, contiguous :: anaerobic_frac(:,:)

    ! soil water potential in each soil layer (MPa)
    real(rk8), pointer, contiguous :: soilpsi(:,:)
    ! O2 consumption during decomposition in each soil layer
    ! (nlevsoi) (mol/m3/s)
    real(rk8), pointer, contiguous :: o2_decomp_depth_unsat(:,:)
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8), pointer, contiguous :: conc_o2_unsat(:,:)
    ! O2 consumption during decomposition in each soil layer
    ! (nlevsoi) (mol/m3/s)
    real(rk8), pointer, contiguous :: o2_decomp_depth_sat(:,:)
    ! O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(rk8), pointer, contiguous :: conc_o2_sat(:,:)
    ! fractional inundated area in soil column (excluding dedicated
    ! wetland columns)
    real(rk8), pointer, contiguous :: finundated(:)
    real(rk8), pointer, contiguous :: sucsat(:,:) ! minimum soil suction (mm)

    ! column 3D org (kg/m3 organic matter) (nlevgrnd)
    real(rk8), pointer, contiguous :: cellorg(:,:)
    character(len=32) :: subname='nitrif_denitrif' ! subroutine name
    real(rk8) :: arg

    !-- implicit in
#ifdef LCH4
    pH                       => clm3%g%l%c%cps%pH
#else
    pH(lbc:ubc) = 6.5_rkx !!! set all soils with the same pH
#endif
    phr_vr                   => clm3%g%l%c%ccf%phr_vr
    w_scalar                 => clm3%g%l%c%ccf%w_scalar
    t_scalar                 => clm3%g%l%c%ccf%t_scalar
    h2osoi_vol               => clm3%g%l%c%cws%h2osoi_vol
    h2osoi_liq               => clm3%g%l%c%cws%h2osoi_liq
    watsat                   => clm3%g%l%c%cps%watsat
    t_soisno                 => clm3%g%l%c%ces%t_soisno
    smin_nh4_vr              => clm3%g%l%c%cns%smin_nh4_vr
    smin_no3_vr              => clm3%g%l%c%cns%smin_no3_vr
    bd                       => clm3%g%l%c%cps%bd
    dz                       => clm3%g%l%c%cps%dz
    watfc                    => clm3%g%l%c%cps%watfc
    bsw                      => clm3%g%l%c%cps%bsw

    soilpsi                  => clm3%g%l%c%cps%soilpsi
#ifdef LCH4
    o2_decomp_depth_unsat    => clm3%g%l%c%cch4%o2_decomp_depth_unsat
    conc_o2_unsat            => clm3%g%l%c%cch4%conc_o2_unsat
    o2_decomp_depth_sat      => clm3%g%l%c%cch4%o2_decomp_depth_sat
    conc_o2_sat              => clm3%g%l%c%cch4%conc_o2_sat
    finundated               => clm3%g%l%c%cws%finundated
#endif
    sucsat                   => clm3%g%l%c%cps%sucsat

    r_psi                   => clm3%g%l%c%cnf%r_psi
    anaerobic_frac          => clm3%g%l%c%cnf%anaerobic_frac

    ! ! subsets of the n flux calcs (for diagnostic/debugging purposes)
    smin_no3_massdens_vr    => clm3%g%l%c%cnf%smin_no3_massdens_vr
    k_nitr_t_vr             => clm3%g%l%c%cnf%k_nitr_t_vr
    k_nitr_ph_vr            => clm3%g%l%c%cnf%k_nitr_ph_vr
    k_nitr_h2o_vr           => clm3%g%l%c%cnf%k_nitr_h2o_vr
    k_nitr_vr               => clm3%g%l%c%cnf%k_nitr_vr
    wfps_vr                 => clm3%g%l%c%cnf%wfps_vr
    fmax_denit_carbonsubstrate_vr => &
            clm3%g%l%c%cnf%fmax_denit_carbonsubstrate_vr
    fmax_denit_nitrate_vr   => clm3%g%l%c%cnf%fmax_denit_nitrate_vr
    f_denit_base_vr         => clm3%g%l%c%cnf%f_denit_base_vr
    diffus                  => clm3%g%l%c%cnf%diffus
    ratio_k1                => clm3%g%l%c%cnf%ratio_k1
    ratio_no3_co2           => clm3%g%l%c%cnf%ratio_no3_co2
    soil_co2_prod           => clm3%g%l%c%cnf%soil_co2_prod
    fr_WFPS                 => clm3%g%l%c%cnf%fr_WFPS
    soil_bulkdensity        => clm3%g%l%c%cnf%soil_bulkdensity

    cellorg   => clm3%g%l%c%cps%cellorg

    !-- implicit out
    pot_f_nit_vr                 => clm3%g%l%c%cnf%pot_f_nit_vr
    pot_f_denit_vr               => clm3%g%l%c%cnf%pot_f_denit_vr
    n2_n2o_ratio_denit_vr        => clm3%g%l%c%cnf%n2_n2o_ratio_denit_vr

    k_nitr_max =  0.1_rk8 / secspday   ! [1/sec] 10%/day  Parton et al., 2001

    co2diff_con(1) =   0.1325_rk8
    co2diff_con(2) =   0.0009_rk8

    do j = 1, nlevdecomp
      do fc = 1, num_soilc
        c = filter_soilc(fc)

        !---------------- calculate soil anoxia state
        ! calculate gas diffusivity of soil at field capacity here
        ! use expression from methane code, but neglect OM for now
        f_a = 1._rk8 - watfc(c,j) / watsat(c,j)
        ! Air-filled fraction of total soil volume
        eps =  watsat(c,j)-watfc(c,j)

        ! use diffusivity calculation including peat
#ifdef LCH4
        if (organic_max > 0._rk8) then
          om_frac = min(cellorg(c,j)/organic_max, 1._rk8)
          ! Use first power, not square as in iniTimeConst
        else
          om_frac = 1._rk8
        end if
        diffus (c,j) = (d_con_g(2,1) + d_con_g(2,2)*t_soisno(c,j)) * 1.e-4_rk8 * &
              (om_frac * f_a**(10._rk8/3._rk8) / watsat(c,j)**2 + &
              (1._rk8-om_frac) * eps**2 * f_a**(3._rk8 / bsw(c,j)) )

        ! calculate anoxic fraction of soils
        ! use rijtema and kroess model after Riley et al., 2000
        ! caclulated r_psi as a function of psi
        r_min(c,j) = 2.0_rk8 * surface_tension_water / &
               (rho_w * grav * abs(soilpsi(c,j)))
        r_max = 2.0_rk8 * surface_tension_water / (rho_w * grav * 0.1_rk8)
        r_psi(c,j) = sqrt(r_min(c,j) * r_max)
        ratio_diffusivity_water_gas(c,j) = (d_con_g(2,1) + &
               d_con_g(2,2)*t_soisno(c,j) ) * 1.e-4_rk8 / &
              ((d_con_w(2,1) + d_con_w(2,2)*t_soisno(c,j) + &
              d_con_w(2,3)*t_soisno(c,j)**2) * 1.e-9_rk8)

        if ( o2_decomp_depth_unsat(c,j) /= spval .and. &
             conc_o2_unsat(c,j) /= spval .and. &
             o2_decomp_depth_unsat(c,j) > 0._rk8 ) then
          arg = -rij_kro_a * r_psi(c,j)**(-rij_kro_alpha) *        &
            o2_decomp_depth_unsat(c,j)**(-rij_kro_beta) *          &
            conc_o2_unsat(c,j)**rij_kro_gamma * (h2osoi_vol(c,j) + &
            ratio_diffusivity_water_gas(c,j) * watsat(c,j))**rij_kro_delta
          if ( arg < 25._rk8 ) then
            anaerobic_frac(c,j) = exp(-arg)
          else
            anaerobic_frac(c,j) = 0._rk8
          end if
        else
          anaerobic_frac(c,j) = 0._rk8
        end if

        if (anoxia_wtsat) then
         ! Average saturated fraction values into anaerobic_frac(c,j).
          r_min_sat = 2._rk8 * surface_tension_water / &
                 (rho_w * grav * abs(grav * 1.e-6_rk8 * sucsat(c,j)))
          r_psi_sat = sqrt(r_min_sat * r_max)
          if (o2_decomp_depth_sat(c,j) /= spval .and. &
             conc_o2_sat(c,j) /= spval .and. &
             o2_decomp_depth_sat(c,j) > 0._rk8) then
            anaerobic_frac_sat = exp(-rij_kro_a * &
                   r_psi_sat**(-rij_kro_alpha) * &
                   o2_decomp_depth_sat(c,j)**(-rij_kro_beta) * &
                   conc_o2_sat(c,j)**rij_kro_gamma * &
                   (watsat(c,j) + ratio_diffusivity_water_gas(c,j) * &
                   watsat(c,j))**rij_kro_delta)
          else
            anaerobic_frac_sat = 0._rk8
          end if
          anaerobic_frac(c,j) = (1._rk8 - finundated(c))*anaerobic_frac(c,j) + &
                 finundated(c)*anaerobic_frac_sat
        end if
#else
        ! NITRIF_DENITRIF requires Methane model to be active,
        ! otherwise diffusivity will be zeroed out here. EBK CDK 10/18/2011
        anaerobic_frac(c,j) = 0._rk8
        diffus (c,j) = 0._rk8
        !call fatal(__FILE__,__LINE__, &
        ! trim(subname)// &
        !   ' ERROR: NITRIF_DENITRIF requires Methane model to be active' )
#endif
        !---------------- nitrification
        ! follows CENTURY nitrification scheme (Parton et al., (2001, 1996))

        ! assume nitrification temp function equal to the HR scalar
        k_nitr_t_vr(c,j) = min(t_scalar(c,j), 1._rk8)

        ! ph function from Parton et al., (2001, 1996)
        k_nitr_ph_vr(c,j) = 0.56 + atan(rpi * 0.45 * (-5.+ pH(c)))/rpi

        ! moisture function-- assume the same moisture function as
        ! limits heterotrophic respiration
        ! Parton et al. base their nitrification- soil moisture rate
        ! constants based on heterotrophic rates-- can we do the same?
        k_nitr_h2o_vr(c,j) = w_scalar(c,j)

        ! nitrification constant is a set scalar * temp, moisture,
        ! and ph scalars
        k_nitr_vr(c,j) = k_nitr_max * k_nitr_t_vr(c,j) * &
               k_nitr_h2o_vr(c,j) * k_nitr_ph_vr(c,j)

        ! first-order decay of ammonium pool with scalar defined above
        pot_f_nit_vr(c,j) = max(smin_nh4_vr(c,j) * k_nitr_vr(c,j), 0._rk8)

        ! limit to oxic fraction of soils
        pot_f_nit_vr(c,j)  = pot_f_nit_vr(c,j) * (1._rk8 - anaerobic_frac(c,j))

        ! limit to non-frozen soil layers
        if ( t_soisno(c,j) <= tfrz .and. no_frozen_nitrif_denitrif) then
          pot_f_nit_vr(c,j) = 0._rk8
        end if

        !---------------- denitrification
        ! first some input variables an unit conversions
        soil_hr_vr(c,j) = phr_vr(c,j)

        ! CENTURY papers give denitrification in units of per gram soil;
        ! need to convert from volumetric to mass-based units here
        soil_bulkdensity(c,j) = bd(c,j) + h2osoi_liq(c,j)/dz(c,j)

        g_per_m3__to__ug_per_gsoil = 1.e3_rk8 / soil_bulkdensity(c,j)

        g_per_m3_sec__to__ug_per_gsoil_day = &
               g_per_m3__to__ug_per_gsoil * secspday

        smin_no3_massdens_vr(c,j) = max(smin_no3_vr(c,j), 0._rk8) * &
               g_per_m3__to__ug_per_gsoil

        soil_co2_prod(c,j) = (soil_hr_vr(c,j) * &
               (g_per_m3_sec__to__ug_per_gsoil_day))

        !! maximum potential denitrification rates based on heterotrophic
        !! respiration rates or nitrate concentrations,
        !! from (del Grosso et al., 2000)
        fmax_denit_carbonsubstrate_vr(c,j) = (0.1_rk8 * &
               (soil_co2_prod(c,j)**1.3_rk8)) &
              / g_per_m3_sec__to__ug_per_gsoil_day
        !
        fmax_denit_nitrate_vr(c,j) = (1.15_rk8 * &
               smin_no3_massdens_vr(c,j)**0.57_rk8)  &
              / g_per_m3_sec__to__ug_per_gsoil_day

        ! find limiting denitrification rate
        f_denit_base_vr(c,j) = max(min(fmax_denit_carbonsubstrate_vr(c,j), &
               fmax_denit_nitrate_vr(c,j)),0._rk8)

        ! limit to non-frozen soil layers
        if ( t_soisno(c,j) <= tfrz .and. no_frozen_nitrif_denitrif ) then
          f_denit_base_vr(c,j) = 0._rk8
        end if

        ! limit to anoxic fraction of soils
        pot_f_denit_vr(c,j) = f_denit_base_vr(c,j) * anaerobic_frac(c,j)

        ! now calculate the ratio of N2O to N2 from denitrifictaion,
        ! following Del Grosso et al., 2000
        ! diffusivity constant (figure 6b)
        ratio_k1(c,j) = max(1.7_rk8, 38.4_rk8 - 350._rk8 * diffus(c,j))

        ! ratio function (figure 7c)
        if ( soil_co2_prod(c,j) > 0 ) then
          ratio_no3_co2(c,j) = smin_no3_massdens_vr(c,j) / soil_co2_prod(c,j)
        else
          ! fucntion saturates at large no3/co2 ratios, so set as some
          ! nominally large number
          ratio_no3_co2(c,j) = 100._rk8
        end if

        ! total water limitation function (Del Grosso et al., 2000, figure 7a)
        wfps_vr(c,j) = max(min(h2osoi_vol(c,j)/watsat(c, j), 1._rk8), 0._rk8) * &
               100._rk8
        fr_WFPS(c,j) = max(0.1_rk8, 0.015_rk8 * wfps_vr(c,j) - 0.32_rk8)
#ifdef LCH4
        if (anoxia_wtsat) then
          fr_WFPS(c,j) = fr_WFPS(c,j)*(1._rk8 - finundated(c)) + &
                 finundated(c)*1.18_rk8
        end if
#endif
        ! final ratio expression
        arg = 0.8 * ratio_no3_co2(c,j)
        if ( arg < 25.0_rk8 ) then
          n2_n2o_ratio_denit_vr(c,j) = max(0.16_rk8*ratio_k1(c,j), &
              ratio_k1(c,j) * exp(-0.8 * ratio_no3_co2(c,j))) * fr_WFPS(c,j)
        else
          n2_n2o_ratio_denit_vr(c,j) = 0.16_rk8*ratio_k1(c,j)
        end if
      end do
    end do
  end subroutine nitrif_denitrif
#endif
#endif

end module mod_clm_cnnitrifdenitrif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

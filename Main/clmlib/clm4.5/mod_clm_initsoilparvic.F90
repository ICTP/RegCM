module mod_clm_initsoilparvic
  !
  ! Performs mapping between VIC and CLM layers
  !
#if (defined VICHYDRO)
  use mod_realkinds
  use mod_intkinds
  use mod_clm_type
  use mod_clm_varcon , only : denh2o , denice , pondmx
  use mod_clm_varpar , only : nlevsoi , nlayer , nlayert , nlevgrnd

  implicit none

  private

  save

  public :: initSoilParVIC      ! map clm soil parameters to vic parameters

  contains
  !
  ! This subroutine converts default CLM soil properties to VIC
  ! parameters to be used for runoff simulations
  ! added by M. Huang
  !
  subroutine initSoilParVIC(c, claycol, sandcol, om_fraccol)
    implicit none
    integer(ik4) , intent(in)  :: c ! column bounds
    ! read in - soil texture: percent sand
    real(rk8) , pointer :: sandcol(:,:)
    ! read in - soil texture: percent clay
    real(rk8) , pointer :: claycol(:,:)
    ! read in - organic matter: kg/m3
    real(rk8) , pointer :: om_fraccol(:,:)
    ! porosity of organic soil
    real(rk8) :: om_watsat = 0.9_rk8
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(rk8) :: om_hksat = 0.1_rk8
    ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(rk8) :: om_tkm = 0.25_rk8
    ! saturated suction for organic matter (Letts, 2000)
    real(rk8) :: om_sucsat = 10.3_rk8
    ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(rk8) :: om_csol = 2.5_rk8
    ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(rk8) :: om_tkd = 0.05_rk8
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(rk8) :: om_b = 2.7_rk8
    ! soil expt for VIC
    real(rk8) :: om_expt = 3._rk8+2._rk8*2.7_rk8
    ! organic matter (kg/m3) where soil is assumed to act like peat
    real(rk8) :: organic_max = 130._rk8
    ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
    real(rk8) :: csol_bedrock = 2.0e6_rk8
    ! percolation threshold
    real(rk8) :: pc = 0.5_rk8
    ! percolation exponent
    real(rk8) :: pcbeta = 0.139_rk8
    real(rk8) :: xksat       ! maximum hydraulic conductivity of soil [mm/s]
    real(rk8) :: perc_frac   ! "percolating" fraction of organic soil
    real(rk8) :: perc_norm   ! normalize to 1 when 100% organic soil
    real(rk8) :: uncon_hksat ! series conductivity of mineral/organic soil
    real(rk8) :: uncon_frac  ! fraction of "unconnected" soil

    real(rk8), pointer :: dz(:,:)  !layer depth (m)
    real(rk8), pointer :: zi(:,:)  !interface level below a "z" level (m)
    real(rk8), pointer :: z(:,:)   !layer thickness (m)
    !fraction of VIC layers in CLM layers
    real(rk8), pointer :: vic_clm_fract(:,:,:)
    ! sum of node fractions in each VIC layer
    real(rk8) :: temp_sum_frac
    !temporary, weighted averaged sand% for VIC layers
    real(rk8) :: sandvic(1:nlayert)
    !temporary, weighted averaged clay% for VIC layers
    real(rk8) :: clayvic(1:nlayert)
    !temporary, weighted averaged organic matter fract for VIC layers
    real(rk8) :: om_fracvic(1:nlayert)

    !b infiltration parameter
    real(rk8), pointer :: b_infil(:)
    !fracton of Dsmax where non-linear baseflow begins
    real(rk8), pointer :: ds(:)
    !max. velocity of baseflow (mm/day)
    real(rk8), pointer :: dsmax(:)
    !fraction of maximum soil moisutre where non-liear base flow occurs
    real(rk8), pointer :: Wsvic(:)
    !baseflow exponent (Qb)
    real(rk8), pointer :: c_param(:)
    !pore-size distribution related paramter(Q12)
    real(rk8), pointer :: expt(:,:)
    !Saturated hydrologic conductivity (mm/s)
    real(rk8), pointer :: ksat(:,:)
    !soil moisture dissusion parameter
    real(rk8), pointer :: phi_s(:,:)
    !layer depth of upper layer(m)
    real(rk8), pointer :: depth(:,:)
    !soil porosity
    real(rk8), pointer :: porosity(:,:)
    !maximum soil moisture (ice + liq)
    real(rk8), pointer :: max_moist(:,:)

    ! other local variables

    integer(ik4) :: i , j

    ! Assign local pointers to derived subtypes components (column-level)

    dz         => clm3%g%l%c%cps%dz
    zi         => clm3%g%l%c%cps%zi
    z          => clm3%g%l%c%cps%z
    depth      => clm3%g%l%c%cps%depth
    vic_clm_fract => clm3%g%l%c%cps%vic_clm_fract

    ! VIC soil parameters
    b_infil        => clm3%g%l%c%cps%b_infil
    dsmax          => clm3%g%l%c%cps%dsmax
    ds             => clm3%g%l%c%cps%ds
    Wsvic          => clm3%g%l%c%cps%Wsvic
    c_param        => clm3%g%l%c%cps%c_param
    expt           => clm3%g%l%c%cps%expt
    ksat           => clm3%g%l%c%cps%ksat
    phi_s          => clm3%g%l%c%cps%phi_s
    porosity       => clm3%g%l%c%cps%porosity
    max_moist      => clm3%g%l%c%cps%max_moist

    !  map parameters between VIC layers and CLM layers

    c_param(c) = 2._rk8
    ! map the CLM layers to VIC layers
    ! There might have better way to do this process

    do i = 1 , nlayer
      sandvic(i) = 0._rk8
      clayvic(i) = 0._rk8
      om_fracvic(i) = 0._rk8
      temp_sum_frac = 0._rk8
      do j = 1 , nlevsoi
        sandvic(i) = sandvic(i) + sandcol(c,j) * vic_clm_fract(c,i,j)
        clayvic(i) = clayvic(i) + claycol(c,j) * vic_clm_fract(c,i,j)
        om_fracvic(i) = om_fracvic(i) + om_fraccol(c,j) * vic_clm_fract(c,i,j)
        temp_sum_frac =temp_sum_frac + vic_clm_fract(c,i,j)
      end do
      !average soil properties, M.Huang, 08/11/2010
      sandvic(i) = sandvic(i)/temp_sum_frac
      clayvic(i) = clayvic(i)/temp_sum_frac
      om_fracvic(i) = om_fracvic(i)/temp_sum_frac
      !make sure sand, clay and om fractions are between 0 and 100%
      sandvic(i) = min(100._rk8, sandvic(i))
      clayvic(i) = min(100._rk8, clayvic(i))
      om_fracvic(i) = min(100._rk8, om_fracvic(i))
      sandvic(i) = max(0._rk8, sandvic(i))
      clayvic(i) = max(0._rk8, clayvic(i))
      om_fracvic(i) = max(0._rk8, om_fracvic(i))
      !calculate other parameters based on teh percentages
      porosity(c, i) = 0.489_rk8 - 0.00126_rk8*sandvic(i)
      expt(c, i) = 3._rk8+ 2._rk8*(2.91_rk8 + 0.159_rk8*clayvic(i))
      xksat = 0.0070556 *( 10.**(-0.884+0.0153*sandvic(i)) )
      !consider organic matter, M.Huang
      expt(c, i) = (1._rk8-om_fracvic(i))*expt(c, i) + om_fracvic(i)*om_expt
      porosity(c,i) = (1._rk8 - om_fracvic(i))*porosity(c,i) + &
              om_watsat*om_fracvic(i)
      ! perc_frac is zero unless perf_frac greater than percolation threshold
      if (om_fracvic(i) > pc) then
        perc_norm = (1._rk8 - pc)**(-pcbeta)
        perc_frac = perc_norm*(om_fracvic(i) - pc)**pcbeta
      else
        perc_frac = 0._rk8
      endif
      ! uncon_frac is fraction of mineral soil plus fraction of
      ! "nonpercolating" organic soil
      uncon_frac=(1._rk8-om_fracvic(i))+(1._rk8-perc_frac)*om_fracvic(i)
      ! uncon_hksat is series addition of mineral/organic conductivites
      if (om_fracvic(i) < 1._rk8) then
        uncon_hksat=uncon_frac/((1._rk8-om_fracvic(i))/xksat &
                    +((1._rk8-perc_frac)*om_fracvic(i))/om_hksat)
      else
        uncon_hksat = 0._rk8
      end if
      ksat(c,i)  = uncon_frac*uncon_hksat + (perc_frac*om_fracvic(i))*om_hksat
      max_moist(c,i) = porosity(c,i)*depth(c,i)*1000._rk8 !in mm!

      phi_s(c,i)=-(exp((1.54_rk8 - 0.0095_rk8*sandvic(i) + &
               0.0063_rk8*(100.0_rk8-sandvic(i)-clayvic(i)))*log(10.0_rk8))*9.8e-5_rk8)
    end do
  end subroutine initSoilParVIC
#endif

end module mod_clm_initsoilparvic
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

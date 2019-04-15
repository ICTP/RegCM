module mod_clm_cnmresp
#ifdef CN
  !
  ! Module holding maintenance respiration routines for coupled carbon
  ! nitrogen code.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_clm_varpar , only : nlevgrnd
  use mod_clm_varcon , only : tfrz

  implicit none

  save

  private

  public :: CNMResp

  contains

  subroutine CNMResp(lbc, ubc, num_soilc, filter_soilc, num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_pftvarcon , only : npcropmin
    use mod_clm_subgridave , only : p2c
    use mod_clm_varctl , only : q10_maintenance
    implicit none
    integer(ik4), intent(in) :: lbc, ubc  ! column-index bounds
    integer(ik4), intent(in) :: num_soilc ! number of soil points in col filter
    integer(ik4), intent(in) :: filter_soilc(:) ! column filter for soil points
    integer(ik4), intent(in) :: num_soilp  ! number of soil points in pft filter
    integer(ik4), intent(in) :: filter_soilp(:) ! pft filter for soil points

    ! column level
    ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8), pointer :: t_soisno(:,:)
    ! pft level
    ! 2 m height surface air temperature (Kelvin)
    real(rk8), pointer :: t_ref2m(:)
    real(rk8), pointer :: q10m(:)
    real(rk8), pointer :: col_q10m(:)
    real(rk8), pointer :: leafn(:)      ! (gN/m2) leaf N
    real(rk8), pointer :: frootn(:)     ! (gN/m2) fine root N
    real(rk8), pointer :: livestemn(:)  ! (gN/m2) live stem N
    real(rk8), pointer :: livecrootn(:) ! (gN/m2) live coarse root N
    real(rk8), pointer :: grainn(:)     ! (kgN/m2) grain N
    ! fraction of roots in each soil layer  (nlevgrnd)
    real(rk8), pointer :: rootfr(:,:)
    integer(ik4) , pointer :: ivt(:)       ! pft vegetation type
    integer(ik4) , pointer :: pcolumn(:)   ! index into column level quantities
    integer(ik4) , pointer :: plandunit(:) ! index into land level quantities
    integer(ik4) , pointer :: clandunit(:) ! index into land level quantities
    integer(ik4) , pointer :: itypelun(:)  ! landunit type
    ! ecophysiological constants
    ! binary flag for woody lifeform (1=woody, 0=not woody)
    real(rk8), pointer :: woody(:)
    logical , pointer :: croplive(:) ! Flag, true if planted, not harvested

    ! pft level
    real(rk8), pointer :: leaf_mr(:)
    real(rk8), pointer :: froot_mr(:)
    real(rk8), pointer :: livestem_mr(:)
    real(rk8), pointer :: livecroot_mr(:)
    real(rk8), pointer :: grain_mr(:)
    ! sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8), pointer :: lmrsun(:)
    ! shaded leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rk8), pointer :: lmrsha(:)
    real(rk8), pointer :: laisun(:)   ! sunlit projected leaf area index
    real(rk8), pointer :: laisha(:)   ! shaded projected leaf area index
    ! fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4) , pointer :: frac_veg_nosno(:)

    integer(ik4) :: c,p,j  ! indices
    integer(ik4) :: fp     ! soil filter pft index
    integer(ik4) :: fc     ! soil filter column index
    real(rk8):: br   ! base rate (gC/gN/s)
    real(rk8):: q10  ! temperature dependence
    real(rk8):: tc   ! temperature correction, 2m air temp (unitless)
    real(rk8):: tcc   ! converting from Kelvin to Celisus
    ! temperature correction by soil layer (unitless)
    real(rk8):: tcsoi(lbc:ubc,nlevgrnd)
    type(pft_estate_type), pointer :: peisos
    type(column_estate_type), pointer :: ceisos
    peisos => clm3%g%l%c%p%pes
    ceisos => clm3%g%l%c%ces

    ! Assign local pointers to derived type arrays
    t_soisno       => clm3%g%l%c%ces%t_soisno
    t_ref2m        => clm3%g%l%c%p%pes%t_ref2m
    col_q10m       => ceisos%pes_a%q10m
    q10m           => peisos%q10m
    leafn          => clm3%g%l%c%p%pns%leafn
    frootn         => clm3%g%l%c%p%pns%frootn
    livestemn      => clm3%g%l%c%p%pns%livestemn
    livecrootn     => clm3%g%l%c%p%pns%livecrootn
    grainn         => clm3%g%l%c%p%pns%grainn
    rootfr         => clm3%g%l%c%p%pps%rootfr
    leaf_mr        => clm3%g%l%c%p%pcf%leaf_mr
    froot_mr       => clm3%g%l%c%p%pcf%froot_mr
    livestem_mr    => clm3%g%l%c%p%pcf%livestem_mr
    livecroot_mr   => clm3%g%l%c%p%pcf%livecroot_mr
    grain_mr       => clm3%g%l%c%p%pcf%grain_mr
    lmrsun         => clm3%g%l%c%p%pcf%lmrsun
    lmrsha         => clm3%g%l%c%p%pcf%lmrsha
    laisun         => clm3%g%l%c%p%pps%laisun
    laisha         => clm3%g%l%c%p%pps%laisha
    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    ivt            => clm3%g%l%c%p%itype
    pcolumn        => clm3%g%l%c%p%column
    plandunit      => clm3%g%l%c%p%landunit
    clandunit      => clm3%g%l%c%landunit
    itypelun       => clm3%g%l%itype
    woody          => pftcon%woody
    croplive       => clm3%g%l%c%p%pps%croplive

    ! base rate for maintenance respiration is from:
    ! M. Ryan, 1991. Effects of climate change on plant respiration.
    ! Ecological Applications, 1(2), 157-167.
    ! Original expression is br = 0.0106 molC/(molN h)
    ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
    br = 2.525e-6_rk8
    ! Peter Thornton: 3/13/09
    ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to
    ! 1.5 as part of the tuning to improve seasonal cycle of atmospheric
    ! CO2 concentration in global simulatoins
    !q10 = 1.5_rk8


    ! pft loop for leaves and live wood
    do fp = 1 , num_soilp
      p = filter_soilp(fp)

      ! calculate maintenance respiration fluxes in
      ! gC/m2/s for each of the live plant tissues.
      ! Leaf and live wood MR

      ! samy : calculation of gridded Q10 to take account for
      ! spatial heterogenity
      ! Choosing between 2 models for Q10 : 0 == linear; 1 == polynomial
      ! linear model Ref :
      ! http://www.ntsg.umt.edu/sites/ntsg.umt.edu/files/modis/MOD17UsersGuide2015_v3.pdf
      ! Polynomial model Ref : Modelling temperature acclimation effects
      ! on the carbon dynamics of forest ecosystems in the conterminous
      ! United states; MIN Chen and QIANLAI Zhuang
      !(2013); Tellus B: Chemical and Physical Meteorology
      tcc = t_ref2m(p) - tfrz ! to convert it to Celsius

      if ( q10_maintenance == 1 ) then
         q10m(p) = 2.35665_rk8 - (0.05308_rk8*tcc) + &
           (0.00238_rk8*tcc*tcc) - (0.00004_rk8*tcc*tcc*tcc)
      else
         q10m(p) = 3.22_rk8 - (0.046_rk8*tcc)
      end if

      ! Set this to expected range
      q10m(p) = max(1.0_rkx,min(3.0_rkx,q10m(p)))

      tc = q10m(p)**((t_ref2m(p)-tfrz - 20.0_rk8)/10.0_rk8)
      if ( frac_veg_nosno(p) == 1 ) then
        leaf_mr(p) = lmrsun(p) * laisun(p) * 12.011e-6_rk8 + &
                     lmrsha(p) * laisha(p) * 12.011e-6_rk8
      else
        leaf_mr(p) = 0.0_rk8
      end if

      if ( abs(woody(ivt(p))-1._rk8) < epsilon(1.0) ) then
        livestem_mr(p) = livestemn(p)*br*tc
        livecroot_mr(p) = livecrootn(p)*br*tc
      else if (ivt(p) >= npcropmin) then
        livestem_mr(p) = livestemn(p)*br*tc
        grain_mr(p) = grainn(p)*br*tc
      end if
    end do

   call p2c(num_soilc,filter_soilc,q10m,col_q10m)

   ! column loop to calculate temperature factors in each soil layer
    do j = 1 , nlevgrnd
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        ! calculate temperature corrections for each soil layer, for use in
        ! estimating fine root maintenance respiration with depth
        ! samy : to calculate q10m at column level for soil temperature
        ! correction
        tcsoi(c,j) = col_q10m(c)**((t_soisno(c,j)-tfrz - 20.0_rk8)/10.0_rk8)
      end do
    end do

    ! soil and pft loop for fine root
    do j = 1 , nlevgrnd
      do fp = 1 , num_soilp
        p = filter_soilp(fp)
        c = pcolumn(p)

        ! Fine root MR
        ! rootfr(j) sums to 1.0 over all soil layers, and
        ! describes the fraction of root mass that is in each
        ! layer.  This is used with the layer temperature correction
        ! to estimate the total fine root maintenance respiration as a
        ! function of temperature and N content.

        froot_mr(p) = froot_mr(p) + frootn(p)*br*tcsoi(c,j)*rootfr(p,j)
      end do
    end do
  end subroutine CNMResp

#endif

end module mod_clm_cnmresp
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

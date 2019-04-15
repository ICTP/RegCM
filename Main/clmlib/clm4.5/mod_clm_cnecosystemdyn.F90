module mod_clm_cnecosystemdyn
#ifdef CN
  !
  ! Ecosystem dynamics: phenology, vegetation
  !
  use mod_intkinds
  use mod_realkinds
  use mod_clm_type
  use mod_clm_varctl , only: use_c13, use_c14
  use mod_clm_cnallocation , only : CNAllocationInit
  use mod_clm_cnphenology , only : CNPhenologyInit
  use mod_clm_cnfire , only : CNFireInit , CNFireArea, CNFireFluxes
  use mod_clm_cnc14decay , only : C14_init_BombSpike
  use mod_clm_surfrd , only : crop_prog
  use mod_clm_cnsetvalue , only : CNZeroFluxes
  use mod_clm_cnndynamics , only : CNNDeposition , CNNFixation
  use mod_clm_cnndynamics , only : CNNLeaching, CNNFert, CNSoyfix
  use mod_clm_cnmresp , only : CNMResp
  use mod_clm_cndecomp , only : CNDecompAlloc
  use mod_clm_cnphenology , only : CNPhenology
  use mod_clm_cngresp , only : CNGResp
  use mod_clm_cncstateupdate1 , only : CStateUpdate1,CStateUpdate0
  use mod_clm_cnnstateupdate1 , only : NStateUpdate1
  use mod_clm_cngapmortality , only : CNGapMortality
  use mod_clm_cncstateupdate2 , only : CStateUpdate2, CStateUpdate2h
  use mod_clm_cnnstateupdate2 , only : NStateUpdate2, NStateUpdate2h
  use mod_clm_cncstateupdate3 , only : CStateUpdate3
  use mod_clm_cnnstateupdate3 , only : NStateUpdate3
  use mod_clm_cnprecisioncontrol , only : CNPrecisionControl
  use mod_clm_cnvegstructupdate , only : CNVegStructUpdate
  use mod_clm_cnannualupdate , only : CNAnnualUpdate
  use mod_clm_cnsummary , only : CSummary, NSummary
  use mod_clm_cncisoflux , only : CIsoFlux1, CIsoFlux2, CIsoFlux2h, CIsoFlux3
  use mod_clm_cnc14decay , only : C14Decay, C14BombSpike
  use mod_clm_pftdyn , only : CNHarvest
  use mod_clm_cnwoodproducts , only : CNWoodProducts
  use mod_clm_cnsoillittverttransp , only : CNSoilLittVertTransp

  implicit none

  private

  save

  public :: CNEcosystemDynInit   ! Ecosystem dynamics initialization
  public :: CNEcosystemDyn       ! Ecosystem dynamics: phenology, vegetation

  contains
  !
  ! Initialzation of the CN Ecosystem dynamics.
  !
  subroutine CNEcosystemDynInit(lbg, ubg, lbc, ubc, lbp, ubp )
    implicit none
    integer(ik4) , intent(in) :: lbg , ubg        ! gridcell bounds
    integer(ik4) , intent(in) :: lbc , ubc        ! column bounds
    integer(ik4) , intent(in) :: lbp , ubp        ! pft bounds
    call CNAllocationInit ( lbc, ubc, lbp, ubp )
    call CNPhenologyInit  ( lbp, ubp )
    call CNFireInit       ( lbg, ubg )
    if ( use_c14 ) call C14_init_BombSpike()
  end subroutine CNEcosystemDynInit
  !
  ! The core CN code is executed here. Calculates fluxes for maintenance
  ! respiration, decomposition, allocation, phenology, and growth respiration.
  ! These routines happen on the radiation time step so that canopy structure
  ! stays synchronized with albedo calculations.
  !
  subroutine CNEcosystemDyn(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                     num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb)
    implicit none
    integer(ik4), intent(in) :: lbc, ubc ! column bounds
    integer(ik4), intent(in) :: lbp, ubp ! pft bounds
    integer(ik4), intent(in) :: num_soilc ! number of soil columns in filter
    ! filter for soil columns
    integer(ik4), intent(in) :: filter_soilc(ubc-lbc+1)
    integer(ik4), intent(in) :: num_soilp ! number of soil pfts in filter
    ! filter for soil pfts
    integer(ik4), intent(in) :: filter_soilp(ubp-lbp+1)
    ! number of prog. crop pfts in filter
    integer(ik4), intent(in) :: num_pcropp
    ! filter for prognostic crop pfts
    integer(ik4), intent(in) :: filter_pcropp(ubp-lbp+1)
    ! true = surface albedo calculation time step
    logical, intent(in) :: doalb

    ! Call the main CN routines
    call CNZeroFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)

    call CNNDeposition(lbc, ubc)

    call CNNFixation(num_soilc,filter_soilc)

    if ( crop_prog ) then
      call CNNFert(num_soilc,filter_soilc)
      call CNSoyfix(num_soilc, filter_soilc, num_soilp, filter_soilp)
    end if

    call CNMResp(lbc, ubc, num_soilc, filter_soilc, num_soilp, filter_soilp)

    call CNDecompAlloc(lbp, ubp, lbc, ubc, num_soilc, filter_soilc, &
                          num_soilp, filter_soilp )

    ! CNphenology needs to be called after CNdecompAlloc, because it
    ! depends on current time-step fluxes to new growth on the last
    ! litterfall timestep in deciduous systems

    call CNPhenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                     num_pcropp, filter_pcropp, doalb)

    call CNGResp(num_soilp, filter_soilp)

    call CStateUpdate0(num_soilp, filter_soilp, 'bulk')

    if ( use_c13 ) then
      call CStateUpdate0(num_soilp, filter_soilp, 'c13')
      call CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')
    end if

    if ( use_c14 ) then
      call CStateUpdate0(num_soilp, filter_soilp, 'c14')
      call CIsoFlux1(num_soilc,filter_soilc,num_soilp,filter_soilp,'c14')
    end if

    call CStateUpdate1(num_soilc,filter_soilc,num_soilp,filter_soilp,'bulk')

    if ( use_c13 ) then
      call CStateUpdate1(num_soilc,filter_soilc,num_soilp,filter_soilp,'c13')
    end if
    if ( use_c14 ) then
      call CStateUpdate1(num_soilc,filter_soilc,num_soilp,filter_soilp,'c14')
    end if

    call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)

    call CNSoilLittVertTransp(lbc, ubc, num_soilc, filter_soilc)

    call CNGapMortality(num_soilc, filter_soilc, num_soilp, filter_soilp)

    if ( use_c13 ) then
      call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')
    end if

    if ( use_c14 ) then
      call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
    end if

    call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')

    if ( use_c13 ) then
      call CStateUpdate2(num_soilc,filter_soilc,num_soilp,filter_soilp,'c13')
    end if

    if ( use_c14 ) then
      call CStateUpdate2(num_soilc,filter_soilc,num_soilp,filter_soilp,'c14')
    end if

    call NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)

#ifdef DYNPFT
    call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp)
#endif

    if ( use_c13 ) then
      call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')
    end if

    if ( use_c14 ) then
      call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
    end if

    call CStateUpdate2h(num_soilc,filter_soilc,num_soilp,filter_soilp,'bulk')

    if ( use_c13 ) then
      call CStateUpdate2h(num_soilc,filter_soilc,num_soilp,filter_soilp,'c13')
    end if

    if ( use_c14 ) then
      call CStateUpdate2h(num_soilc,filter_soilc,num_soilp,filter_soilp,'c14')
    end if

    call NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)

    call CNWoodProducts(num_soilc, filter_soilc)

    call CNFireArea(num_soilc, filter_soilc,num_soilp, filter_soilp)

    call CNFireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)

    call CNNLeaching(lbc, ubc, num_soilc, filter_soilc)

    if ( use_c13 ) then
      call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')
    end if

    if ( use_c14 ) then
      call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
    end if

    call CStateUpdate3(num_soilc,filter_soilc,num_soilp,filter_soilp,'bulk')

    if ( use_c13 ) then
      call CStateUpdate3(num_soilc,filter_soilc,num_soilp,filter_soilp,'c13')
    end if

    if ( use_c14 ) then
      call CStateUpdate3(num_soilc,filter_soilc,num_soilp,filter_soilp,'c14')
      call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp)
      call C14BombSpike(num_soilp, filter_soilp)
    end if

    call NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)

    call CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)

    if (doalb) then
      call CNVegStructUpdate(num_soilp, filter_soilp)
    end if

    call CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')

    if ( use_c13 ) then
      call CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')
    end if

    if ( use_c14 ) then
      call CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
    end if

    call NSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)

  end subroutine CNEcosystemDyn
#endif

end module mod_clm_cnecosystemdyn
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

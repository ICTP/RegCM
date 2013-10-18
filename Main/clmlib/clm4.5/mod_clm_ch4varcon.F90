module mod_clm_ch4varcon

#ifdef LCH4
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ch4varcon
!
! !DESCRIPTION:
! Module containing CH4 parameters and logical switches and routine to read constants from CLM namelist.
!
! !USES:
  use mod_stdio
  use mod_realkinds
  use mod_mpmessage
  use mod_dynparam
  use mod_mppparam
  use mod_clm_varctl  , only : NLFileName_in
!
! !PUBLIC TYPES:
  implicit none
  save

!
! Methane Model Parameters
!
  real(rk8) :: q10ch4base = 295.D0  ! Rough estimate from comparison between Walter and previous CLM-CH4 data
  ! Uses Michigan bog data from Shannon & White
  ! This is the temperature at which the effective f_ch4 actually equals the constant f_ch4.

  real(rk8) :: q10ch4 = 1.33D0 ! Production Q10
  ! Note that this is the additional Q10 for methane production ABOVE the soil decomposition temperature relationship.
  ! Corresponds to a methanogenesis Q10 of 2 when SOM HR has Q10 of 1.5.
  ! (No assumption is made in CH4 code about SOM HR temperature relationship.)
  ! Note that this formulation should be improved by making anaerobic decomposition in general a stronger function of
  ! temperature, not just methane production.

  real(rk8) :: vmax_ch4_oxid = 45.D-6 * 1000.D0 / 3600.D0 ! [mol/m3-w/s];
  ! 45 uM/h from Walter and Heimann for the Mich. site (2000)
  ! oxidation rate constant (Walter and Heimann 2000)

  real(rk8) :: k_m = 5.D-6 * 1000.D0 ! [mol/m3-w]
  ! Michaelis-Menten oxidation rate constant for CH4 concentration (Walter and Heimann 2000)

  real(rk8) :: q10_ch4oxid = 1.9D0 ! Segers, 1998
  ! Q10 oxidation constant (Walter and Heimann 2000)

  real(rk8) :: smp_crit = -2.4D5 ! mm. From Schnell & King, 1996.
  ! Critical soil moisture potential to reduce oxidation (mm) due to dessication of methanotrophs above the water table.
  ! To turn off limitation, set to very large negative value.

  real(rk8) :: aereoxid = -1.D0 ! fraction of methane flux entering aerenchyma rhizosphere that will be
  ! oxidized rather than emitted.  From Wania.
  ! Note, this has been replaced by prognostic O2 diffusion into aerenchyma and is set to -1 by default.
  ! Set to value between 0 & 1 (inclusive) for sensitivity tests.  In particular, to calculate the 
  ! prognostic fraction analogous to the Wania parameter, set to 0. and compare to a normal run.

  real(rk8) :: mino2lim = 0.2D0 ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate
  ! for soil decomposition (or diagnostic O2-limitation / seasonal inundation factor)

  real(rk8) :: rootlitfrac = 0.50D0 ! Fraction of soil organic matter associated with roots
  ! Used to partition the production between rootfr and the top 5 layers (ifndef VERTSOILC)

  real(rk8) :: scale_factor_aere = 1.0D0 ! scale factor on the aerenchyma area for sensitivity tests

  real(rk8) :: vgc_max = 0.15D0 ! ratio of saturation pressure triggering ebullition

  real(rk8) :: organic_max  = 130.D0 ! organic matter content (kg/m3) where soil is assumed to act like peat
  ! for diffusion. Very large values will lead to all soil being treated as mineral. Negative values will lead   ! to all soil being treated as peat.

  real(rk8) :: satpow = 2.D0 ! exponent on watsat for saturated soil solute diffusion
  ! (2 = Buckingham / Moldrup; 4/3 = Millington-Quirk)

  real(rk8) :: cnscalefactor = 1.D0 ! scale factor on CN decomposition for assigning methane flux
                                    ! This should equal 1 except for sensitivity studies.

  real(rk8) :: f_ch4 = 0.2D0
                          ! originally 25% / (100% + 25%) from Wania. This is the ratio of CH4 production to total C
                          ! mineralization.
                          ! fraction of total decomposition that comes off as CH4 rather than CO2
                          ! Effective value will depend on temperature, redox, & pH but cannot exceed 50%
                          ! based on stoichiometry.
                          ! Note this is a crude parameterization: values in the field vary widely.

  logical :: transpirationloss = .true. ! switch for activating CH4 loss from transpiration
                                      ! Transpiration loss assumes that the methane concentration in dissolved soil
                                      ! water remains constant through the plant and is released when the water evaporates
                                      ! from the stomata.
                                      ! Currently hard-wired to true; impact is < 1 Tg CH4/yr

  real(rk8) :: k_m_o2 = 20.D-6 * 1000.D0 ! [mol/m3-w]
  ! Michaelis-Menten oxidation rate constant for O2 concentration (Segers 1998, Lidstrom and Somers 1984)

  real(rk8) :: nongrassporosratio = 0.33D0 ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport
                                           ! Some values in Colmer 2003

  logical :: allowlakeprod = .false. ! Switch to allow production under lakes based on soil carbon dataset
                                     ! (Methane can be produced, and CO2 produced from methane oxidation,
                                     ! which will slowly reduce the available carbon stock, if ! replenishlakec, but no other biogeochem is done.)
                                     ! Note: switching this off turns off ALL lake methane biogeochem. However, 0 values
                                     ! will still be averaged into the concentration _sat history fields.

  real(rk8) :: lake_decomp_fact = 9.D-11 ! Base decomposition rate (1/s) at 25C
                       ! Equates to about a 200 year lifetime.

  logical :: usephfact = .false. ! Switch to use pH factor in methane production

  real(rk8) :: k_m_unsat = 5.D-6 * 1000.D0 / 10.D0 ! [mol/m3-w]
  ! Michaelis-Menten oxidation rate constant for CH4 concentration: literature suggests that methanotrophs
  ! in upland areas have higher affinity for methane in order to access the low ambient concentrations above the water
  ! table. (See Bender & Conrad, 1992, etc.)

  real(rk8) :: vmax_oxid_unsat = 45.D-6 * 1000.D0 / 3600.D0 / 10.D0 ! [mol/m3-w/s]
  ! Literature suggests that while k_m is lower, vmax is also lower in upland areas.

  logical :: replenishlakec = .true. ! Switch for keeping carbon storage under lakes constant
                                      ! so that lakes do not affect the carbon balance
                                      ! Good for long term rather than transient warming experiments
               ! NOTE SWITCHING THIS OFF ASSUMES TRANSIENT CARBON SUPPLY FROM LAKES; COUPLED MODEL WILL NOT CONSERVE CARBON
               ! IN THIS MODE.

  real(rk8) :: scale_factor_gasdiff = 1.0D0 ! For sensitivity tests; convection would allow this to be > 1

  real(rk8) :: scale_factor_liqdiff = 1.0D0 ! For sensitivity tests; convection would allow this to be > 1

  real(rk8) :: redoxlag = 30.0D0 ! Number of days to lag in the calculation of finundated_lag, which will
                                 ! be used to assess the availability of alternative electron acceptors in recently
                                 ! inundated area, reducing production.
                                 ! Set to 0 to turn off this feature.

  ! New namelists added 6/12/11

  logical :: fin_use_fsat = .false. ! Use fsat rather than the inversion to Prigent satellite inundation obs. (applied to
                                    ! CLM water table depth and surface runoff) to calculated finundated which is
                                    ! used in methane code and potentially soil code
                                    !!!! Attn EK: Set this to true when Sean Swenson's prognostic, tested
                                       ! fsat is integrated. (CLM4 fsat is bad for these purposes.)

  real(rk8) :: unsat_aere_ratio = 0.05D0 / 0.3D0 ! Ratio to multiply upland vegetation aerenchyma porosity by compared to
                                        ! inundated systems. Note: porosity will be kept at above a minimum residual value
                                        ! porosmin set in subroutine ch4_aere.

  logical :: usefrootc = .false.    ! Use CLMCN fine root C rather than ann NPP & LAI based parameterization to
                                    ! calculate tiller C for aerenchyma area calculation.
                                    ! The NPP & LAI param. was based on Wania for Arctic sedges and may not be
                                    ! appropriate for woody PFTs, although nongrassporosratio above partly adjusts
                                    ! for this.  However, using fine root C reduces the aerenchyma area by a large
                                    ! factor.

  logical :: ch4offline = .true.    ! true --> Methane is not passed between the land & atmosphere.
                                    ! NEM is not added to NEE flux to atm. to correct for methane production,
                                    ! and ambient CH4 is set to constant 2009 value.

  logical :: ch4rmcnlim = .false.   ! Remove the N and low moisture limitations on SOM HR when calculating
                                    ! methanogenesis.
                                    ! Note: this option has not been extensively tested.
                                    ! Currently hardwired off.

  logical :: anoxicmicrosites = .false. ! Use Arah & Stephen 1998 expression to allow production above the water table
                                        ! Currently hardwired off; expression is crude.

  logical :: ch4frzout = .false.    ! Exclude CH4 from frozen fraction of soil pore H2O, to simulate "freeze-out" pulse
                                    ! as in Mastepanov 2008.
                                    ! Causes slight increase in emissions in the fall and decrease in the spring.
                                    ! Currently hardwired off; small impact.

  real(rk8) :: redoxlag_vertical = 0.D0   ! time lag (days) to inhibit production for newly unsaturated layers
                                          ! when decreasing WT depth for unsat. zone. See the description for redoxlag.

  real(rk8) :: atmch4 = 1.7D-6    ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model
                                    ! (mol/mol)
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: ch4conrd ! Read and initialize CH4 constants
!
! !REVISION HISTORY:
! Created by Zack Subin
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4conrd
!
! !INTERFACE:
  subroutine ch4conrd ()
!
! !DESCRIPTION:
! Read and initialize CH4 constants
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initch4
!
! !REVISION HISTORY:
! Created by Zack Subin
! Modified 6/12/2011 to use namelist rather than ASCII file.
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,n                ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'ch4conrd'  ! subroutine name

!-----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! Driver
    namelist /ch4par_in/ &
         ch4offline, fin_use_fsat, replenishlakec, allowlakeprod, atmch4

    ! Production
    namelist /ch4par_in/ &
         q10ch4base, q10ch4, rootlitfrac, f_ch4, cnscalefactor, &
         redoxlag, & ! ch4rmcnlim, anoxicmicrosites,
         mino2lim, lake_decomp_fact, usephfact, redoxlag_vertical

    ! Oxidation
    namelist /ch4par_in/ &
         vmax_ch4_oxid, k_m, q10_ch4oxid, smp_crit, k_m_o2, k_m_unsat, vmax_oxid_unsat

    ! Aerenchyma
    namelist /ch4par_in/ &
         aereoxid, scale_factor_aere,  & !transpirationloss,
         nongrassporosratio, unsat_aere_ratio, &
         usefrootc

    ! Ebullition
    namelist /ch4par_in/ &
         vgc_max

    ! Transport
    namelist /ch4par_in/ &
         organic_max, satpow, scale_factor_gasdiff, scale_factor_liqdiff !, ch4frzout


       ! ----------------------------------------------------------------------
       ! Read namelist from standard input.
       ! ----------------------------------------------------------------------

    if (myid == iocpu) then

       write(stdout,*) 'Attempting to read CH4 parameters .....'
       unitn = file_getunit( )
       write(stdout,*) 'Read in ch4par_in namelist from: ', trim(NLFilename_in)
       open( unitn, file=trim(NLFilename_in), status='old' )
       read(unitn, ch4par_in, iostat=ierr)
       if (ierr /= 0) then
         call fatal(__FILE__,__LINE__, &
              subname//' error in reading in ch4par_in namelist' )
       end if
       call file_freeunit( unitn )

    end if ! masterproc

    call bcast(q10ch4base)
    call bcast(q10ch4)
    call bcast(vmax_ch4_oxid)
    call bcast(k_m)
    call bcast(q10_ch4oxid)
    call bcast(smp_crit)
    call bcast(aereoxid)
    call bcast(mino2lim)
    call bcast(rootlitfrac)
    call bcast(scale_factor_aere)
    call bcast(vgc_max)
    call bcast(organic_max)
    call bcast(satpow)
    call bcast(cnscalefactor)
    call bcast(f_ch4)
    ! call bcast(transpirationloss)
    call bcast(k_m_o2)
    call bcast(nongrassporosratio)
    call bcast(allowlakeprod)
    call bcast(lake_decomp_fact)
    call bcast(usephfact)
    call bcast(k_m_unsat)
    call bcast(vmax_oxid_unsat)
    call bcast(replenishlakec)
    call bcast(scale_factor_gasdiff)
    call bcast(scale_factor_liqdiff)
    call bcast(redoxlag)
    call bcast(fin_use_fsat)
    call bcast(unsat_aere_ratio)
    call bcast(usefrootc)
    call bcast(ch4offline)
    ! call bcast(ch4rmcnlim)
    ! call bcast(anoxicmicrosites)
    ! call bcast(ch4frzout)
    call bcast(redoxlag_vertical)
    call bcast(atmch4)

    if (myid == iocpu) then
       write(stdout,*) 'Successfully read CH4 namelist'
       write(stdout,*)' '
       write(stdout,*)'q10ch4base = ', q10ch4base
       write(stdout,*)'q10ch4 = ', q10ch4
       write(stdout,*)'vmax_ch4_oxid = ', vmax_ch4_oxid
       write(stdout,*)'k_m = ', k_m
       write(stdout,*)'q10_ch4oxid = ', q10_ch4oxid
       write(stdout,*)'smp_crit = ', smp_crit
       write(stdout,*)'aereoxid = ', aereoxid
       write(stdout,*)'mino2lim = ', mino2lim
       write(stdout,*)'rootlitfrac = ', rootlitfrac
       write(stdout,*)'scale_factor_aere = ', scale_factor_aere
       write(stdout,*)'vgc_max = ', vgc_max
       write(stdout,*)'organic_max = ', organic_max
       write(stdout,*)'satpow = ', satpow
       write(stdout,*)'cnscalefactor = ', cnscalefactor
       write(stdout,*)'f_ch4 = ', f_ch4
       !write(stdout,*)'transpirationloss = ', transpirationloss
       write(stdout,*)'k_m_o2 = ', k_m_o2
       write(stdout,*)'nongrassporosratio = ', nongrassporosratio
       write(stdout,*)'allowlakeprod = ', allowlakeprod
       write(stdout,*)'lake_decomp_fact = ', lake_decomp_fact
       write(stdout,*)'usephfact = ', usephfact
       write(stdout,*)'k_m_unsat = ', k_m_unsat
       write(stdout,*)'vmax_oxid_unsat = ', vmax_oxid_unsat
       write(stdout,*)'replenishlakec = ', replenishlakec
       write(stdout,*)'scale_factor_gasdiff = ', scale_factor_gasdiff
       write(stdout,*)'scale_factor_liqdiff = ', scale_factor_liqdiff
       write(stdout,*)'redoxlag = ', redoxlag
       write(stdout,*)'fin_use_fsat = ', fin_use_fsat
       write(stdout,*)'unsat_aere_ratio = ', unsat_aere_ratio
       write(stdout,*)'usefrootc = ', usefrootc
       write(stdout,*)'ch4offline = ', ch4offline
       !write(stdout,*)'ch4rmcnlim = ', ch4rmcnlim
       !write(stdout,*)'anoxicmicrosites = ', anoxicmicrosites
       !write(stdout,*)'ch4frzout = ', ch4frzout
       write(stdout,*)'redoxlag_vertical = ', redoxlag_vertical
       write(stdout,*)'atmch4 = ', atmch4

       if (ch4offline) write(stdout,*)'CH4 Model will be running offline and not affect fluxes to atmosphere.'
       if (aereoxid >= 0.D0) write(stdout,*) 'Fixed aerenchyma oxidation has been selected.'
       if (.not. allowlakeprod) &
         write(stdout,*) 'Lake production has been disabled.  Lakes will '//&
         'not factor into CH4 BGC.  "Sat" history fields will not average '//&
         'over lakes except for concentrations, which will average zero '//&
         'from lakes.'
       if (.not. replenishlakec .and. .not. ch4offline) &
         write(stdout,*)'LAKE SOIL CARBON WILL NOT BE REPLENISHED BUT '//&
         'INSTEAD WILL BE TRANSIENTLY RELEASED: COUPLED MODEL WILL NOT '//&
         'CONSERVE CARBON IN THIS MODE!'
       write(stdout,*)'Successfully initialized CH4 parameters from namelist.'
       write(stdout,*)

    end if

  end subroutine ch4conrd

#endif
! defined LCH4

end module mod_clm_ch4varcon

#include <misc.h>
#include <preproc.h>

module VOCEmissionMod

#if (defined VOC)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: VOCEmissionMod
!
! !DESCRIPTION:
! Volatile organic compound emission
!
! !USES:
  use abortutils, only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: VOCEmission
!abt add below
  private :: gamma_calc
!abt above
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: VOCEmission
!
! !INTERFACE:
  subroutine VOCEmission (lbp, ubp, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Volatile organic compound emission
! This code simulates volatile organic compound emissions
! following the algorithm presented in Guenther, A., 1999: Modeling
! Biogenic Volatile Organic Compound Emissions to the Atmosphere. In
! Reactive Hydrocarbons in the Atmosphere, Ch. 3
! This model relies on the assumption that 90% of isoprene and monoterpene
! emissions originate from canopy foliage:
!    E = epsilon * gamma * density * delta
! The factor delta (longterm activity factor) applies to isoprene emission
! from deciduous plants only. We neglect this factor at the present time.
! This factor is discussed in Guenther (1997).
! Subroutine written to operate at the patch level.
! IN FINAL IMPLEMENTATION, REMEMBER:
! 1. may wish to call this routine only as freq. as rad. calculations
! 2. may wish to place epsilon values directly in pft-physiology file
! Output: vocflx(nvoc) !VOC flux [ug m-2 h-1]
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varpar   , only : nvoc,nlevsoi
    use clm_atmlnd   , only : clm_a2l
    use shr_const_mod, only : SHR_CONST_RGAS
!abt added below
    use clm_varvoc
    use clm_time_manager, only : get_step_size, get_curr_calday, get_curr_date
    use clm_time_manager, only : get_nstep
    use clm_varctl      , only : nsrest
!abt added above
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(num_nolakep) ! pft filter for non-lake points
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
! 2/1/02, Peter Thornton: migration to new data structure
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pgridcell(:)     ! gridcell index of corresponding pft
    integer , pointer :: ivt(:)           ! pft vegetation type for current
    real(r8), pointer :: t_veg(:)         ! pft vegetation temperature (Kelvin)
    real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
    real(r8), pointer :: elai(:)          ! one-sided leaf area index with burying by snow
    real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (visible only)
    real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation     (visible only)
    real(r8), pointer :: sla(:)           ! ecophys constant - specific leaf area [m2 leaf g-1 carbon]
!
!!!abt rcm below
    real(r8), pointer :: epsilon_iso(:)   ! emission factor map for isoprene
    real(r8), pointer :: epsilon_apin(:)  ! emission factor map for a-pinene
    real(r8), pointer :: epsilon_bpin(:)  ! emission factor map for b-pinene
    real(r8), pointer :: epsilon_mbo(:)   ! emission factor map for methylbutenol
    real(r8), pointer :: epsilon_myrc(:)  ! emission factor map for myrcene
    real(r8), pointer :: epsilon_sabi(:)  ! emission factor map for sabinene
    real(r8), pointer :: epsilon_limo(:)  ! emission factor map for limonene
    real(r8), pointer :: epsilon_acar(:)  ! emission factor map for a-3carene
    real(r8), pointer :: epsilon_ocim(:)  ! emission factor map for ocimene
    real(r8), pointer :: epsilon_omtp(:)  ! emission factor map for other monoterps
    real(r8), pointer :: epsilon_farn(:)  ! emission factor map for farnicene
    real(r8), pointer :: epsilon_bcar(:)  ! emission factor map for b-caryophyllene
    real(r8), pointer :: epsilon_osqt(:)  ! emission factor map for other sesquiterps
    real(r8), pointer :: epsilon_meoh(:)  ! emission factor map for methanol
    real(r8), pointer :: epsilon_acto(:)  ! emission factor map for acetone
    real(r8), pointer :: epsilon_meth(:)  ! emission factor map for methane
    real(r8), pointer :: epsilon_no(:)    ! emission factor map for no,n2o,nh3
    real(r8), pointer :: epsilon_acta(:)  ! emission factor map for acetaldehyde etc
    real(r8), pointer :: epsilon_form(:)  ! emission factor map for formaldehyde formic acid
    real(r8), pointer :: epsilon_co(:)    ! emission factor map for carbonmoxide
!!!abt rcm above
!
! local pointers to original implicit out arrays
!
    real(r8), pointer :: vocflx(:,:)      ! VOC flux [ug m-2 h-1]
    real(r8), pointer :: vocflx_tot(:)    ! VOC flux [ug m-2 h-1]
    real(r8), pointer :: vocflx_1(:)      ! VOC flux(1) [ug m-2 h-1]
    real(r8), pointer :: vocflx_2(:)      ! VOC flux(2) [ug m-2 h-1]
    real(r8), pointer :: vocflx_3(:)      ! VOC flux(3) [ug m-2 h-1]
    real(r8), pointer :: vocflx_4(:)      ! VOC flux(4) [ug m-2 h-1]
    real(r8), pointer :: vocflx_5(:)      ! VOC flux(5) [ug m-2 h-1]
!abt added below
    real(r8), pointer :: vocflx_6(:)      ! VOC flux(6) [ug m-2 h-1]
    real(r8), pointer :: vocflx_7(:)      ! VOC flux(7) [ug m-2 h-1]
    real(r8), pointer :: vocflx_8(:)      ! VOC flux(8) [ug m-2 h-1]
    real(r8), pointer :: vocflx_9(:)      ! VOC flux(9) [ug m-2 h-1]
    real(r8), pointer :: vocflx_10(:)     ! VOC flux(10) [ug m-2 h-1]
    real(r8), pointer :: vocflx_11(:)     ! VOC flux(11) [ug m-2 h-1]
    real(r8), pointer :: vocflx_12(:)     ! VOC flux(12) [ug m-2 h-1]
    real(r8), pointer :: vocflx_13(:)     ! VOC flux(13) [ug m-2 h-1]
    real(r8), pointer :: vocflx_14(:)     ! VOC flux(14) [ug m-2 h-1]
    real(r8), pointer :: vocflx_15(:)     ! VOC flux(15) [ug m-2 h-1]
    real(r8), pointer :: vocflx_16(:)     ! VOC flux(16) [ug m-2 h-1]
    real(r8), pointer :: vocflx_17(:)     ! VOC flux(17) [ug m-2 h-1]
    real(r8), pointer :: vocflx_18(:)     ! VOC flux(18) [ug m-2 h-1]
    real(r8), pointer :: vocflx_19(:)     ! VOC flux(19) [ug m-2 h-1]
    real(r8), pointer :: vocflx_20(:)     ! VOC flux(20) [ug m-2 h-1]
!abt added above
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: fp,p,g,n         ! indices
    real(r8) :: epsilon(lbp:ubp) ! emission factor [ug g-1 h-1]
    real(r8) :: gamma(lbp:ubp)   ! activity factor (instantaneous light and temp. condns)
    real(r8) :: density          ! source density factor [g dry wgt foliar mass/m2 ground]
    real(r8) :: cl               ! temporary
    real(r8) :: ct               ! temporary
    real(r8) :: par              ! temporary
    real(r8) :: reciprod         ! temporary
!abt below added
    real(r8) :: gamma_ce(lbp:ubp,nvoc)         ! activation factor for LAI,temp,humidity,& light in canopy (abt)
    real(r8) :: gamma_age(lbp:ubp,nvoc)        ! activation factor for leaf age (abt)
    real(r8) :: gamma_sm(lbp:ubp)              ! activation factor for soil moisture (abt)
!abt above added
! Constants

    real(r8), parameter :: R   = SHR_CONST_RGAS*0.001_r8 ! univ. gas constant [J K-1 mol-1]
    real(r8), parameter :: alpha = 0.0027_r8 ! empirical coefficient
    real(r8), parameter :: cl1 = 1.066_r8    ! empirical coefficient
    real(r8), parameter :: ct1 = 95000.0_r8  ! empirical coefficient [J mol-1]
    real(r8), parameter :: ct2 = 230000.0_r8 ! empirical coefficient [J mol-1]
    real(r8), parameter :: ct3 = 0.961_r8    ! empirical coefficient
    real(r8), parameter :: tm  = 314.0_r8    ! empirical coefficient [K]
    real(r8), parameter :: tstd = 303.0_r8   ! std temperature [K]
    real(r8), parameter :: bet = 0.09_r8     ! beta empirical coefficient [K-1]
    real(r8), parameter :: scale_factor = 5._r8/7._r8  ! J-F Larmaque's empirical factor

! These are the values from version of genesis-ibis / 1000.
! With DGVM defined, use LPJ's sla [m2 leaf g-1 carbon]
! Divide by 2 in the equation to get dry weight foliar mass from grams carbon

    real(r8) :: hardwire_sla(0:16)
    real(r8) :: slarea(lbp:ubp)           ! Specific leaf areas [m2 leaf g-1 carbon]
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    forc_solad => clm_a2l%forc_solad
    forc_solai => clm_a2l%forc_solai

    epsilon_iso  => clm3%g%gem%epsilon_iso     !abt
    epsilon_bpin => clm3%g%gem%epsilon_bpin    !abt
    epsilon_apin => clm3%g%gem%epsilon_apin    !abt
    epsilon_mbo  => clm3%g%gem%epsilon_mbo     !abt
    epsilon_sabi => clm3%g%gem%epsilon_sabi    !abt
    epsilon_limo => clm3%g%gem%epsilon_limo    !abt
    epsilon_myrc => clm3%g%gem%epsilon_myrc    !abt
    epsilon_acar => clm3%g%gem%epsilon_acar    !abt
    epsilon_ocim => clm3%g%gem%epsilon_ocim    !abt
    epsilon_omtp => clm3%g%gem%epsilon_omtp    !abt
    epsilon_farn => clm3%g%gem%epsilon_farn    !abt
    epsilon_bcar => clm3%g%gem%epsilon_bcar    !abt
    epsilon_osqt => clm3%g%gem%epsilon_osqt    !abt
    epsilon_meoh => clm3%g%gem%epsilon_meoh    !abt
    epsilon_acto => clm3%g%gem%epsilon_acto    !abt
    epsilon_meth => clm3%g%gem%epsilon_meth    !abt
    epsilon_no   => clm3%g%gem%epsilon_no      !abt
    epsilon_acta => clm3%g%gem%epsilon_acta    !abt
    epsilon_form => clm3%g%gem%epsilon_form    !abt
    epsilon_co   => clm3%g%gem%epsilon_co      !abt

    ! Assign local pointers to derived subtypes components (pft-level)

    pgridcell  => clm3%g%l%c%p%gridcell
    ivt        => clm3%g%l%c%p%itype
    t_veg      => clm3%g%l%c%p%pes%t_veg
    fsun       => clm3%g%l%c%p%pps%fsun
    elai       => clm3%g%l%c%p%pps%elai
    vocflx     => clm3%g%l%c%p%pvf%vocflx
    vocflx_tot => clm3%g%l%c%p%pvf%vocflx_tot
    vocflx_1   => clm3%g%l%c%p%pvf%vocflx_1
    vocflx_2   => clm3%g%l%c%p%pvf%vocflx_2
    vocflx_3   => clm3%g%l%c%p%pvf%vocflx_3
    vocflx_4   => clm3%g%l%c%p%pvf%vocflx_4
    vocflx_5   => clm3%g%l%c%p%pvf%vocflx_5
    vocflx_6   => clm3%g%l%c%p%pvf%vocflx_6     !abt
    vocflx_7   => clm3%g%l%c%p%pvf%vocflx_7     !abt
    vocflx_8   => clm3%g%l%c%p%pvf%vocflx_8     !abt
    vocflx_9   => clm3%g%l%c%p%pvf%vocflx_9     !abt
    vocflx_10  => clm3%g%l%c%p%pvf%vocflx_10    !abt
    vocflx_11  => clm3%g%l%c%p%pvf%vocflx_11    !abt
    vocflx_12  => clm3%g%l%c%p%pvf%vocflx_12    !abt
    vocflx_13  => clm3%g%l%c%p%pvf%vocflx_13    !abt
    vocflx_14  => clm3%g%l%c%p%pvf%vocflx_14    !abt
    vocflx_15  => clm3%g%l%c%p%pvf%vocflx_15    !abt
    vocflx_16  => clm3%g%l%c%p%pvf%vocflx_16    !abt
    vocflx_17  => clm3%g%l%c%p%pvf%vocflx_17    !abt
    vocflx_18  => clm3%g%l%c%p%pvf%vocflx_18    !abt
    vocflx_19  => clm3%g%l%c%p%pvf%vocflx_19    !abt
    vocflx_20  => clm3%g%l%c%p%pvf%vocflx_20    !abt
    sla        => pftcon%sla

#if (!defined DGVM)
    hardwire_sla( 0) = 0._r8
    hardwire_sla( 1) = 0.0125_r8 !needleleaf
    hardwire_sla( 2) = 0.0125_r8 !Gordon Bonan suggests NET = 0.0076
    hardwire_sla( 3) = 0.0125_r8 !Gordon Bonan suggests NDT = 0.0200
    hardwire_sla( 4) = 0.0250_r8 !broadleaf
    hardwire_sla( 5) = 0.0250_r8 !Gordon Bonan suggests BET = 0.0178
    hardwire_sla( 6) = 0.0250_r8 !Gordon Bonan suggests BDT = 0.0274
    hardwire_sla( 7) = 0.0250_r8
    hardwire_sla( 8) = 0.0250_r8
    hardwire_sla( 9) = 0.0250_r8
    hardwire_sla(10) = 0.0250_r8
    hardwire_sla(11) = 0.0250_r8
    hardwire_sla(12) = 0.0200_r8 !grass
    hardwire_sla(13) = 0.0200_r8
    hardwire_sla(14) = 0.0200_r8
    hardwire_sla(15) = 0.0200_r8
    hardwire_sla(16) = 0.0200_r8 !numpft = 16
#endif

!!!abt added below
!!!! To calculate indices used in accumulating PPFD and !!!!
!!!!      vegetation temperature during gamma calc      !!!!         

  ! Get indices to be used for time accumulation variables
    c24  = c24 + 1
    if(c24 > n24) c24 = 1
    c240 = c240 + 1
    if(c240 > n240) c240 = 1

  ! calculate the gamma factors for each compound
    call gamma_calc (lbp,ubp,num_nolakep,filter_nolakep,gamma_ce,gamma_age,gamma_sm)







!!!abt added above



    ! Determine specific leaf array
!dir$ concurrent
!cdir nodep
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       g = pgridcell(p)

#if (!defined DGVM)
       slarea(p) = hardwire_sla(ivt(p))
#else
       slarea(p) = sla(ivt(p))
#endif
       vocflx_tot(p) = 0._r8
    end do


    ! Begin loop through voc species

    do n = 1, nvoc

!dir$ concurrent
!cdir nodep
          do fp = 1,num_nolakep
             p = filter_nolakep(fp)
             g = pgridcell(p)


             if(n.eq.1) then
             ! isoprenes: (values received from C. Wiedinmyer)

               epsilon(p) = 0._r8
               if(r2cefmap(n) == 1) then
                  epsilon(p) = epsilon_iso(g)
               else
                  if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 2000._r8
                  else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 13000._r8
                  else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 11000._r8
                  else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 400._r8
                  end if
               endif

!dir$ concurrent
!cdir nodep

             ! Myrcene: (values received from C. Wiedinmyer)
             elseif(n.eq.2) then    

               epsilon(p) = 0._r8
!             epsilon(p) = 1.0_r8                 !Guenther (personal communication)

               if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_myrc(g)
               else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 75._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 20._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 22._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 0.3_r8
                 end if
               endif 

!dir$ concurrent
!cdir nodep

            ! Sabinene: (values received from C. Wiedinmyer)
            elseif(n.eq.3) then

!             epsilon(p) = 1.0_r8                 !Guenther (personal communication)
              epsilon(p) = 0._r8                 !Guenther (personal communication)

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_sabi(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 70._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 45._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 50._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 0.7_r8
                 end if
              endif


!dir$ concurrent
!cdir nodep


            ! Limonene: (values received from C. Wiedinmyer)
             elseif(n.eq.4) then

!             epsilon(p) = 1.0_r8                 !Guenther (personal communication)
               epsilon(p) = 0._r8                 

               if(r2cefmap(n) == 1) then
                  epsilon(p) = epsilon_limo(g)
               else
                  if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                     epsilon(p) = 100._r8
                  else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                     epsilon(p) = 45._r8
                  else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                     epsilon(p) = 52._r8
                  else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                     epsilon(p) = 0.7_r8
                  end if
               endif 

!dir$ concurrent
!cdir nodep

             ! CO
             elseif(n.eq.5) then

!             epsilon(p) = 0.3_r8                 !Guenther (personal communication)
               epsilon(p) = 0._r8  

               if(r2cefmap(n) == 1) then
                  epsilon(p) = epsilon_co(g)    
               else
                  epsilon(p) = 1000._r8      !Wiedinmyer (personal communication)
               endif

!dir$ concurrent
!cdir nodep


            ! methylbutenol: (values received from C. Wiedinmyer)
            elseif(n.eq.6) then

              ! epsilon:  for methylbutenol
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_mbo(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 100._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 0.1_r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 1._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 0.01_r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep

             ! a-pinene: (values received from C. Wiedinmyer)
             elseif(n.eq.7) then

               ! epsilon:  for a-pinene
               epsilon(p) = 0._r8

               if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_apin(g)
               else
!abt changed to increase EMISSIONS
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                  epsilon(p) = 450._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 180._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 200._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 2._r8
                 end if
               endif	 

!dir$ concurrent
!cdir nodep


            ! b-pinene: (values received from C. Wiedinmyer)
            elseif(n.eq.8) then

              ! epsilon: for b-pinene
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_bpin(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 300._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 90._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 100._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 1.5_r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep

            ! ocimene: (values received from C. Wiedinmyer)
            elseif(n.eq.9) then

              ! epsilon: for ocimene
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_ocim(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 90._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 60._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 85._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 1._r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep

            ! a-3carene: (values received from C. Wiedinmyer)
            elseif(n.eq.10) then

              ! epsilon: for a-3carene
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_acar(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 18._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 160._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 25._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 0.3_r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep


            ! other monoterpenes: (values received from C. Wiedinmyer)
            elseif(n.eq.11) then

              ! epsilon: for other monoterpenes
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_omtp(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 90._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 180._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 110._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 4.8_r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep

            ! farnicene: (values received from C. Wiedinmyer)
            elseif(n.eq.12) then

              ! epsilon: for farnicene
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_farn(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 35._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 30._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 30._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 0.5_r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep

            ! b-caryophyllene: (values received from C. Wiedinmyer)
            elseif(n.eq.13) then

              ! epsilon: for b-caryophyllene
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_bcar(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 30._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 60._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 45._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 0.9_r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep


            ! other sesquiterpenes: (values received from C. Wiedinmyer)
            elseif(n.eq.14) then

              ! epsilon: for other sesquiterpenes
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_osqt(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 75._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 110._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 85._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 1.4_r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep

            ! methanol: (values received from C. Wiedinmyer)
            elseif(n.eq.15) then

              ! epsilon: for methanol
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_meoh(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 800._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 800._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 800._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 800._r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep


            ! acetone: (values received from C. Wiedinmyer)
            elseif(n.eq.16) then

              ! epsilon: for acetone
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_acto(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 240._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 240._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 240._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 80._r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep


            ! methane: (values received from C. Wiedinmyer)
            elseif(n.eq.17) then

              ! epsilon: for methane
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_meth(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 30._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 30._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 30._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 30._r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep


            ! no2/n2o/nh3: (values received from C. Wiedinmyer)
            elseif(n.eq.18) then

              ! epsilon: for no2/n2o/nh3
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_no(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 5._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 6._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 30._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 70._r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep


            ! acetaldehyde/ethanol: (values received from C. Wiedinmyer)
            elseif(n.eq.19) then

              ! epsilon: for acetaldehyde/ethanol
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_acta(g)
              else
                 if (ivt(p) >= 1 .and. ivt(p) <= 3) then         !needleleaf
                    epsilon(p) = 240._r8
                 else if (ivt(p) >= 4 .and. ivt(p) <= 8) then    !broadleaf
                    epsilon(p) = 240._r8
                 else if (ivt(p) >= 9 .and. ivt(p) <= 11) then   !shrub
                    epsilon(p) = 240._r8
                 else if (ivt(p) >= 12 .and. ivt(p) <= 16) then  !herbaceous
                    epsilon(p) = 80._r8
                 end if
              endif			 

!dir$ concurrent
!cdir nodep


            ! formic acid/formaldehyde/acetic acid: (values received from C. Wiedinmyer)
            elseif(n.eq.20) then

              ! epsilon: for formic acid/formaldehyde/acetic acid
              epsilon(p) = 0._r8

              if(r2cefmap(n) == 1) then
                 epsilon(p) = epsilon_form(g)
              else
                 epsilon(p) = 70._r8
              endif			 

            else
              write(*,*) "ERROR: voc number exceeds available" ; call endrun()
            endif  !end storing of epsilon values

!dir$ concurrent
!cdir nodep

             ! gamma: Activity factor. Units [dimensionless]

             ! isoprenes:
             ! scale total incident par by fraction of sunlit leaves (added on 1/2002)
             ! multiply w/m2 by 4.6 to get umol/m2/s for par (added 8/14/02)
             ! got this value from subr. Stomata

!abt original below
!             reciprod = 1._r8 / (R * t_veg(p) * tstd)
!             ct = exp(ct1 * (t_veg(p) - tstd) * reciprod) / &
!                 (ct3 + exp(ct2 * (t_veg(p) - tm) * reciprod))

!             par = (forc_solad(g,1) + fsun(p) * forc_solai(g,1)) * 4.6_r8
!             cl = alpha * cl1 * par * (1._r8 + alpha * alpha * par * par)**(-0.5_r8)
!             gamma(p) = cl * ct !gamma = 1 under std temp & light condns

!             par = ((1._r8 - fsun(p)) * forc_solai(g,1)) * 4.6_r8
!             cl = alpha * cl1 * par * (1._r8 + alpha * alpha * par * par)**(-0.5_r8)
!             gamma(p) = gamma(p) + cl * ct !gamma(sun) + gamma(sha)
!abt original above

             gamma(p) = gamma_ce(p,n) * gamma_age(p,n) * gamma_sm(p)

!dir$ concurrent
!cdir nodep


          if(r2cefmap(n) <= 2) then
             ! density here refers to the production/loss rate within canopy (dimensionless)  !abt
             ! currently all density values will be set to 1 except isoprene
             if(n.eq.1) then
               density = 0.96_r8   !***NOTE: density here refers to amount lost or gained within the canopy
             else                  !         not the source density factor described above
               density = 1._r8
             endif                 

          else
             ! density: Source density factor [g dry weight foliar mass m-2 ground]
             if (ivt(p) > 0) then
               density = elai(p) / (slarea(p) * 0.5_r8)
             else
               density = 0._r8
             end if
          endif

          ! calculate the voc flux


          if(n.eq.1) then
               vocflx(p,n) = epsilon(p) * gamma(p) * density * scale_factor
          else
               vocflx(p,n) = epsilon(p) * gamma(p) * density
          end if

          vocflx_tot(p) = vocflx_tot(p) + vocflx(p,n)

          if(n.eq.nvoc) then
             vocflx_1(p) = vocflx(p,1)
             vocflx_2(p) = vocflx(p,2)
             vocflx_3(p) = vocflx(p,3)
             vocflx_4(p) = vocflx(p,4)
             vocflx_5(p) = vocflx(p,5)
             vocflx_6(p) = vocflx(p,6)    !abt added
             vocflx_7(p) = vocflx(p,7)    !abt added
             vocflx_8(p) = vocflx(p,8)    !abt added
             vocflx_9(p) = vocflx(p,9)    !abt added
             vocflx_10(p) = vocflx(p,10)  !abt added
             vocflx_11(p) = vocflx(p,11)  !abt added
             vocflx_12(p) = vocflx(p,12)  !abt added
             vocflx_13(p) = vocflx(p,13)  !abt added
             vocflx_14(p) = vocflx(p,14)  !abt added
             vocflx_15(p) = vocflx(p,15)  !abt added
             vocflx_16(p) = vocflx(p,16)  !abt added
             vocflx_17(p) = vocflx(p,17)  !abt added
             vocflx_18(p) = vocflx(p,18)  !abt added
             vocflx_19(p) = vocflx(p,19)  !abt added
             vocflx_20(p) = vocflx(p,20)  !abt added
          endif


       end do   ! end pft loop


    end do   ! end voc species loop

    ! Calculate total voc flux and individual components for history output

!!dir$ concurrent
!!cdir nodep
!    do fp = 1,num_nolakep
!       p = filter_nolakep(fp)
!       vocflx_tot(p) = 0._r8
!    end do
!    do n = 1, nvoc
!!dir$ concurrent
!!cdir nodep
!       do fp = 1,num_nolakep
!          p = filter_nolakep(fp)
!          vocflx_tot(p) = vocflx_tot(p) + vocflx(p,n)
!       end do
!    end do
!!dir$ concurrent
!!cdir nodep
!    do fp = 1,num_nolakep
!       p = filter_nolakep(fp)
!       vocflx_1(p) = vocflx(p,1)
!       vocflx_2(p) = vocflx(p,2)
!       vocflx_3(p) = vocflx(p,3)
!       vocflx_4(p) = vocflx(p,4)
!       vocflx_5(p) = vocflx(p,5)
!       vocflx_6(p) = vocflx(p,6)    !abt added
!       vocflx_7(p) = vocflx(p,7)    !abt added
!       vocflx_8(p) = vocflx(p,8)    !abt added
!       vocflx_9(p) = vocflx(p,9)    !abt added
!       vocflx_10(p) = vocflx(p,10)  !abt added
!       vocflx_11(p) = vocflx(p,11)  !abt added
!       vocflx_12(p) = vocflx(p,12)  !abt added
!       vocflx_13(p) = vocflx(p,13)  !abt added
!       vocflx_14(p) = vocflx(p,14)  !abt added
!       vocflx_15(p) = vocflx(p,15)  !abt added
!       vocflx_16(p) = vocflx(p,16)  !abt added
!       vocflx_17(p) = vocflx(p,17)  !abt added
!       vocflx_18(p) = vocflx(p,18)  !abt added
!       vocflx_19(p) = vocflx(p,19)  !abt added
!       vocflx_20(p) = vocflx(p,20)  !abt added
!    end do

  end subroutine VOCEmission
  
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gamma_calc
!
! !INTERFACE:
  subroutine gamma_calc (lbp,ubp,num_nolakep,filter_nolakep,gamma_ce,gamma_age,gamma_sm)
!
! !DESCRIPTION:
! Calculates the gamma factors described in Guenther et al 2006
! Uses constant values for empirical values such as Epot,Topt, and alpha
! until there is enough past model information to calculate these 
! values directly.
! EXAMPLE: Topt = 312.5 K from Guenther 1999 until average leaf temperature
!          from the past 240hrs,T240, is calculated.  At which time this eq 
!          will be used:
!          Topt = 313 + (0.6 · (T240œôø²297))   (Guenther 2006)
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varpar   , only : nvoc,nlevsoi
    use clm_atmlnd   , only : clm_a2l
    use shr_const_mod, only : SHR_CONST_RGAS
    use clm_time_manager, only : get_step_size, get_curr_calday, get_curr_date
    use clm_time_manager, only : get_nstep
    use clm_varvoc
    use clm_varctl   , only : nsrest
    use spmdMod      , only : iam
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num_nolakep                    ! indices
    integer , intent(in) :: lbp,ubp                        ! indices
    integer , intent(in) :: filter_nolakep(num_nolakep)    ! indices
    real(r8), intent(inout) :: gamma_ce(lbp:ubp,nvoc)      ! activation factor for LAI,temp,humidity,& light in canopy
    real(r8), intent(inout) :: gamma_age(lbp:ubp,nvoc)     ! activation factor for leaf age
    real(r8), intent(inout) :: gamma_sm(lbp:ubp)           ! activation factor for soil moisture
!
! !CALLED FROM: VOCEmission
!
! !REVISION HISTORY:
! Author: Ahmed Tawfik
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: h2osoi_vol(:,:)  ! volumetric soil water (0<=h2osoi_vol<=watsat [m3/m3]
    real(r8), pointer :: sucsat(:,:)      ! minimum soil suction (mm) (nlevsoi)
    real(r8), pointer :: bsw(:,:)         ! Clapp and Hornberger "b" (nlevsoi)
    real(r8), pointer :: rootfr(:,:)      ! fraction of roots in each soil layer  
    real(r8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity) (nlevsoi)
    integer , pointer :: pcolumn(:)       ! column index of corresponding pft
    integer , pointer :: pgridcell(:)     ! gridcell index of corresponding pft
!!!    real(r8), pointer :: t_veg(:)         ! pft vegetation temperature (Kelvin)
    real(r8), pointer :: t_ref2m(:)       ! 2meter temperature (Kelvin)  *** buggin
    real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
    real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (visible only)
    real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation     (visible only)
    real(r8), pointer :: t_sum24(:)       ! pft vegetation temperature (Kelvin)
    real(r8), pointer :: t_sum240(:)      ! pft vegetation temperature (Kelvin)
    real(r8), pointer :: fsd24(:)         ! pft direct radiation averaged avg 24hrs
    real(r8), pointer :: fsd240(:)        ! pft direct radiation averaged avg 240hrs
    real(r8), pointer :: fsi24(:)         ! pft indirect radiation averaged avg 24hrs
    real(r8), pointer :: fsi240(:)        ! pft indirect radiation averaged avg 240hrs
    real(r8), pointer :: monlai(:,:)      ! monthly lai
    real(r8), pointer :: fsun24(:)        ! pft vegetation temperature (Kelvin)
    real(r8), pointer :: fsun240(:)       ! pft vegetation temperature (Kelvin)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: fp,c          ! indices
    integer  :: p,g,n         ! indices
    real(r8) :: cl            ! temporary
    real(r8) :: ct            ! temporary
    real(r8) :: par           ! temporary
    real(r8) :: x             ! temporary
    real(r8) :: upper,lower   !numerator and denomenator for ct calc
    real(r8) :: calday        !calendar day at Greenwich (1.00 -> 365.99)
    real(r8) :: calday1st     !calendar day for 1st timestep at Greenwich (1.00 -> 365.99)
    integer  :: kda           !day (1 -> 31)
    integer  :: kmo           !month (1 -> 12)
    integer  :: kyr           !year (0 -> ...)
    integer  :: pda           !day (1 -> 31) for previous timestep
    integer  :: pmo           !month (1 -> 12) for previous timestep
    integer  :: pyr           !year (0 -> ...) for previous timestep
    integer  :: ksec          !current seconds of current date (0 -> 86400)
    integer  :: mcdate        !current date in integer format [yyyymmdd]
    integer  :: mcsec         !current time of day [seconds]
    real(r8) :: alpha_sh      ! empirical coefficient (shade)
    real(r8) :: Cp_sh         ! empirical coefficient (shade)
    real(r8) :: alpha_su      ! empirical coefficient (in the sun)
    real(r8) :: Cp_su         ! empirical coefficient (in the sun)
    real(r8) :: Eopt          ! empirical coefficient
    real(r8) :: Topt          ! empirical coefficient	
    real(r8) :: t_24,t_240    ! average temp over 24h and 240h
    real(r8) :: p_240sh,p_240su         ! average PPFD over 24h and 240h
    real(r8) :: p_24su,p_24sh           ! average PPFD over 24h and 240h 
    real(r8) :: p_24,p_240
    real(r8) :: fsun_24,fsun_240
    integer  :: counter,lev
    integer  :: nstep,dtime             ! step size and timestep number in CLM
!    real(r8) :: Aold,Agro,Amat,Anew   ! relative emission rate to each canopy fraction (old,grow,mature,new)
    real(r8) :: Fold,Fgro,Fmat,Fnew     ! canopy fraction (old,grow,mature,new)
    real(r8) :: gamma_pt,gamma_sun,gamma_shade
    real(r8) :: gamma_lai,gamma_t,gamma_p
    real(r8) :: wilt(nlevsoi),wilt1(nlevsoi),sm_sum(nlevsoi) ! wilting point,wilt point + del_vol and
	                                                     ! sum of gamma_sm over all soil levels
	
! Constants
    real(r8), parameter :: R    = SHR_CONST_RGAS*0.001_r8 ! univ. gas constant [J K-1 mol-1]
    real(r8), parameter :: ct1  = 95.0_r8     ! empirical coefficient [J mol-1]
    real(r8), parameter :: ct2  = 230.0_r8    ! empirical coefficient [J mol-1]
    real(r8), parameter :: ct3  = 0.961_r8    ! empirical coefficient
    real(r8), parameter :: tm   = 314.0_r8    ! empirical coefficient [K]
    real(r8), parameter :: tstd = 303.0_r8    ! std temperature [K]
    real(r8), parameter :: bet  = 0.09_r8     ! beta empirical coefficient [K-1]
    real(r8), parameter :: Cce  = 0.4_r8     ! sets emission activity to unity at std conditions
    real(r8), parameter :: Psun = 200.0_r8    ! std PPFD for sun portion [umol/m2/s]
    real(r8), parameter :: Psha = 50.0_r8     ! std PPFD for shaded portion [umol/m2/s]
    real(r8), parameter :: del_vol = 0.06     ! empirical paramter

    logical :: WEIGHTED

!-----------------------------------------------------------------------


  ! Get time variables
    dtime    = get_step_size()
    nstep    = get_nstep()
    calday1st = get_curr_calday(offset=-int(dtime*nstep))
    calday   = get_curr_calday()
    call get_curr_date(kyr, kmo, kda, mcsec)
    call get_curr_date(pyr, pmo, pda, mcsec,offset=-int(dtime))

  ! Assign local pointers to derived type members (gridcell-level)
    forc_solad => clm_a2l%forc_solad
    forc_solai => clm_a2l%forc_solai

  ! Assign local pointers to derived subtypes components (column-level)
    bsw        => clm3%g%l%c%cps%bsw
    watsat     => clm3%g%l%c%cps%watsat
    sucsat     => clm3%g%l%c%cps%sucsat
    h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol
  
  ! Assign local pointers to derived subtypes components (pft-level)
    pgridcell  => clm3%g%l%c%p%gridcell
!!    t_veg      => clm3%g%l%c%p%pes%t_veg
    t_ref2m    => clm3%g%l%c%p%pes%t_ref2m  !** buggin use 2m Temp instead of veg temp ** !!!
    fsun       => clm3%g%l%c%p%pps%fsun
    rootfr     => clm3%g%l%c%p%pps%rootfr	
    pcolumn    => clm3%g%l%c%p%column
    t_sum24    => clm3%g%l%c%p%pva%t_sum24
    t_sum240   => clm3%g%l%c%p%pva%t_sum240
    fsd24      => clm3%g%l%c%p%pva%fsd24
    fsd240     => clm3%g%l%c%p%pva%fsd240
    fsi24      => clm3%g%l%c%p%pva%fsi24
    fsi240     => clm3%g%l%c%p%pva%fsi240
    fsun24     => clm3%g%l%c%p%pva%fsun24
    fsun240    => clm3%g%l%c%p%pva%fsun240
    monlai     => clm3%g%l%c%p%pva%monlai

  ! Set some constants (note these are only constant when there isn't enough info)    
    alpha_su = 0.001_r8 
    Cp_su    = 1.21_r8  
    alpha_sh = 0.001_r8 
    Cp_sh    = 1.21_r8  	  
    Eopt     = 2.26_r8    
    Topt     = 317._r8   

  ! Initialize
    gamma_ce(:,:)  = 0._r8
    gamma_age(:,:) = 0._r8
    gamma_sm(:)    = 0._r8
    t_24      = 0._r8
    t_240     = 0._r8
    p_24su    = 0._r8
    p_240su   = 0._r8
    p_24sh    = 0._r8
    p_240sh   = 0._r8
    fsun_24   = 0._r8
    fsun_240  = 0._r8

  do fp = 1,num_nolakep           !pft loop
      p = filter_nolakep(fp)
      g = pgridcell(p)
      c = pcolumn(p)




  ! calculate gamma_ce:  gamma_ce = gamma_pt * LAI * Cce (immediately below if statements
  ! are for calculating temp & light effects on emissions over last 1 and 10 days)

    if(nstep .ge. n240+n_start) then
      t_240     = t_sum240(p) 
      t_24      = t_sum24(p)  
      p_240sh   = fsi240(p) 
      p_240su   = fsd240(p)
      p_24sh    = fsi24(p) 
      p_24su    = fsd24(p) 
      fsun_24   = fsun24(p) 
      fsun_240  = fsun240(p) 

      p_240su  = (p_240su + fsun_240  * p_240sh ) * 4.6_r8	
      p_24su   = (p_24su  + fsun_24   * p_24sh  ) * 4.6_r8	
      p_24sh   = ((1._r8  - fsun_24)  * p_24sh  ) * 4.6_r8
      p_240sh  = ((1._r8  - fsun_240) * p_240sh ) * 4.6_r8

      alpha_sh = 0.004_r8 - 0.0005_r8*log(p_240sh)
      alpha_su = 0.004_r8 - 0.0005_r8*log(p_240su)
      Cp_su    = 0.0468_r8*exp(0.0005_r8*(p_24su - Psun))*(p_240su)**(0.6_r8)
      Cp_sh    = 0.0468_r8*exp(0.0005_r8*(p_24sh - Psha))*(p_240sh)**(0.6_r8)
      Topt     = 313._r8 + (0.6_r8*(t_240 - 297._r8))
      Eopt     = 2.034_r8*exp(0.05_r8*(t_24 - 297._r8))*exp(0.05_r8*(t_240 - 297._r8))
    endif




   ! calculate gamma_sm: gamma_sm varies depending on the soil moisture and wilting point
   ! wilting point calculation from Muller 2008

     sm_sum(:) = 0._r8
     counter   = 0
     do lev = 1,nlevsoi
	 wilt(lev)  = 0.5_r8*watsat(c,lev)*(200._r8/sucsat(c,lev))**(-1._r8/bsw(c,lev))
         wilt1(lev) = del_vol + wilt(lev)

         sm_sum(2)  = (h2osoi_vol(c,lev) - wilt(lev)) / del_vol
         sm_sum(1)  = rootfr(p,lev) * max(0._r8, min(1._r8, sm_sum(2)) ) + sm_sum(1)
     enddo

     if(sum(rootfr(p,:)) > 0._r8) then
          gamma_sm(p) = sm_sum(1)
     else
          gamma_sm(p) = 0._r8
     end if
     gamma_sm_out(p) = gamma_sm(p)
  



  ! calculate gamma_age: gamma_age = FnewAnew + FoldAold + FmatAmat + FgroAgro
  ! eqs below hold assuming that timestep in days is less than or equal to the number of days
  ! between budbreak and induction of emission (Guenther 2006)
    if(pmo == kmo .or. monlai(p,1) == monlai(p,2)) then
      Fnew = 0.0_r8
      Fgro = 0.1_r8
      Fmat = 0.8_r8
      Fold = 0.1_r8
    elseif(pmo /= kmo) then
      if(monlai(p,1) > monlai(p,2)) then
        Fnew = 0.0_r8
	Fgro = 0.0_r8
	Fold = (monlai(p,1)-monlai(p,2))/monlai(p,1)
        Fmat = 1._r8 - Fold
      elseif(monlai(p,1) < monlai(p,2)) then 
        Fnew = 1._r8 - (monlai(p,1)/monlai(p,2))
        Fmat = monlai(p,1)/monlai(p,2)
        Fgro = 1._r8 - Fnew - Fmat
	Fold = 0.0_r8
      endif
    endif


  do  n = 1,nvoc                  !begin voc loop

  !Calculate light and temperature dependency factor
    if(n.eq.1) then  !isoprene  (Guenther 2006)
!      x     = ((1._r8 /Topt) - (1._r8/t_veg(p))) / 0.00831_r8
      x     = ((1._r8 /Topt) - (1._r8/t_ref2m(p))) / 0.00831_r8  !*** buggin USE 2m Temp **!!!!
      upper = exp(ct1*x)*ct2*Eopt
      lower = ct2 - ct1*(1._r8 - exp(ct2*x))
      ct    = upper/lower
      par = (forc_solad(g,1) + fsun(p) * forc_solai(g,1)) * 4.6_r8
      cl  = alpha_su * Cp_su * par * (1._r8 + alpha_su * alpha_su * par * par)**(-0.5_r8)
      gamma_sun = cl * ct * fsun(p) !gamma = 1 under std temp & light condns
      par = ((1._r8 - fsun(p)) * forc_solai(g,1)) * 4.6_r8
      cl = alpha_sh * Cp_sh * par * (1._r8 + alpha_sh * alpha_sh * par * par)**(-0.5_r8)
      gamma_shade = ct * cl * (1 - fsun(p))
      gamma_pt      = gamma_sun + gamma_shade !gamma(sun) + gamma(sha)
      gamma_ce(p,n) = Cce * monlai(p,2) * gamma_pt

      gamma_t_out(p)  = ct             !for output
      gamma_p_out(p)  = gamma_pt/ct    !for output
      Cpsun_out(p)    = Cp_su
      Cpsh_out(p)     = Cp_sh
      alphasu_out(p)  = alpha_su
      alphash_out(p)  = alpha_sh
      Topt_out(p)     = Topt
      Eopt_out(p)     = Eopt
      t24_out(p)      = t_24
      t240_out(p)     = t_240
      p24su_out(p)    = p_24su
      p24sh_out(p)    = p_24sh
      p240su_out(p)   = p_240su
      p240sh_out(p)   = p_240sh


    else ! non-isoprene (Guenther 2006 and MEGANv2.03)
!!!      gamma_t   = exp(beta_V(n) * (t_veg(p) - tstd))
      gamma_t   = exp(beta_V(n) * (t_ref2m(p) - tstd))     !*** buggin USE 2m Temp ** !!!
      par = (forc_solad(g,1) + fsun(p) * forc_solai(g,1)) * 4.6_r8
      cl  = alpha_su * Cp_su * par * (1._r8 + alpha_su * alpha_su * par * par)**(-0.5_r8)
      gamma_sun = cl !gamma = 1 under std temp & light condns
      par = ((1._r8 - fsun(p)) * forc_solai(g,1)) * 4.6_r8
      cl = alpha_sh * Cp_sh * par * (1._r8 + alpha_sh * alpha_sh * par * par)**(-0.5_r8)
      gamma_shade = cl
      gamma_p = gamma_sun + gamma_shade !gamma(sun) + gamma(shade)
      gamma_lai = (0.49*monlai(p,2))/(1+0.2*monlai(p,2)**2)**(0.5_r8)     
      gamma_ce(p,n) = gamma_lai*gamma_t*((1-LDF(n))+ LDF(n)*gamma_p)
    endif

  ! Calculate Leaf Age Factor    
    gamma_age(p,n) = Fnew*Anew(n) + Fold*Aold(n) + Fmat*Amat(n) + Fgro*Agro(n)
    if(n.eq.1) gamma_age_out(p) = gamma_age(p,n)

  enddo    !end nvoc loop
  enddo    !end pft loop  





  end subroutine gamma_calc


#endif

end module VOCEmissionMod

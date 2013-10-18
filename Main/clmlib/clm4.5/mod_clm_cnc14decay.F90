
module mod_clm_cnc14decay
#if (defined CN)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: C14FluxMod
!
! !DESCRIPTION:
! Module for 14-carbon flux variable update, non-mortality fluxes.
!
! !USES:
    use mod_realkinds
    use mod_nchelper
    use mod_mpmessage
    use mod_stdio
    use mod_dynparam
    use mod_mppparam
    use mod_clm_varpar   , only: ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
    implicit none
    save
    private
!
! !PUBLIC MEMBER FUNCTIONS:
    public:: C14Decay
    public:: C14BombSpike
    public:: C14_init_BombSpike

! !PUBLIC TYPES:
    logical, public :: use_c14_bombspike = .false.         ! do we use time-varying atmospheric C14?
    character(len=256), public :: atm_c14_filename = ' '   ! file name of C14 input data
    
! !PRIVATE TYPES:
    real(rk8), allocatable, private :: atm_c14file_time(:)
    real(rk8), allocatable, private :: atm_delta_c14(:)

! !REVISION HISTORY:
! 2/15/2011 by C. Koven
!
!EOP
!-----------------------------------------------------------------------

contains



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C14Decay
!
! !INTERFACE:
subroutine C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, calculate the radioactive decay of C14
!
! !USES:
   use mod_clm_type
   use mod_clm_time_manager, only: get_step_size, get_days_per_year
   use mod_clm_varcon, only: secspday
   use mod_clm_varctl, only : spinup_state
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays

   real(rk8), pointer :: decomp_cpools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   real(rk8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(rk8), pointer :: xsmrpool(:)           ! (gC/m2) execss maint resp C pool
   real(rk8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(rk8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(rk8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(rk8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(rk8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(rk8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(rk8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(rk8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(rk8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(rk8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(rk8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(rk8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(rk8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(rk8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(rk8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(rk8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(rk8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(rk8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(rk8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(rk8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(rk8), pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   real(rk8), pointer :: seedc(:)

! !OTHER LOCAL VARIABLES:
   integer :: fp,j,l,p,fc,c,i
   real(rk8) :: dt                             ! radiation time step (seconds)
   real(rk8) :: half_life
   real(rk8) :: decay_const
   real(rk8), pointer :: spinup_factor(:)      ! factor for AD spinup associated with each pool
   real(rk8) :: days_per_year                  ! days per year
   real(rk8) :: spinup_term               ! spinup accelerated decomposition factor, used to accelerate transport as well


    ! assign local pointers at the column level
    decomp_cpools_vr                   => clm3%g%l%c%cc14s%decomp_cpools_vr

    ! ! assign local pointers at the column level
    ! new pointers for dynamic landcover
    seedc                          => clm3%g%l%c%cc14s%seedc

   ! assign local pointers at the pft level
    cpool                          => clm3%g%l%c%p%pc14s%cpool
    xsmrpool                       => clm3%g%l%c%p%pc14s%xsmrpool
    deadcrootc                     => clm3%g%l%c%p%pc14s%deadcrootc
    deadcrootc_storage             => clm3%g%l%c%p%pc14s%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%l%c%p%pc14s%deadcrootc_xfer
    deadstemc                      => clm3%g%l%c%p%pc14s%deadstemc
    deadstemc_storage              => clm3%g%l%c%p%pc14s%deadstemc_storage
    deadstemc_xfer                 => clm3%g%l%c%p%pc14s%deadstemc_xfer
    frootc                         => clm3%g%l%c%p%pc14s%frootc
    frootc_storage                 => clm3%g%l%c%p%pc14s%frootc_storage
    frootc_xfer                    => clm3%g%l%c%p%pc14s%frootc_xfer
    gresp_storage                  => clm3%g%l%c%p%pc14s%gresp_storage
    gresp_xfer                     => clm3%g%l%c%p%pc14s%gresp_xfer
    leafc                          => clm3%g%l%c%p%pc14s%leafc
    leafc_storage                  => clm3%g%l%c%p%pc14s%leafc_storage
    leafc_xfer                     => clm3%g%l%c%p%pc14s%leafc_xfer
    livecrootc                     => clm3%g%l%c%p%pc14s%livecrootc
    livecrootc_storage             => clm3%g%l%c%p%pc14s%livecrootc_storage
    livecrootc_xfer                => clm3%g%l%c%p%pc14s%livecrootc_xfer
    livestemc                      => clm3%g%l%c%p%pc14s%livestemc
    livestemc_storage              => clm3%g%l%c%p%pc14s%livestemc_storage
    livestemc_xfer                 => clm3%g%l%c%p%pc14s%livestemc_xfer
    pft_ctrunc                     => clm3%g%l%c%p%pc14s%pft_ctrunc
    spinup_factor                  => decomp_cascade_con%spinup_factor


    ! set time steps
    dt = real( get_step_size(), rk8 )
    days_per_year = get_days_per_year()

    half_life = 5568.D0 * secspday * days_per_year  !! libby half-life value, for comparison against ages calculated with this value
    ! half_life = 5730.D0 * secspday * days_per_year  !! recent half-life value
    decay_const = - log(0.5D0) / half_life


    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       seedc(c) = seedc(c) *  (1.D0 - decay_const * dt)
    end do ! end of columns loop

    do l = 1, ndecomp_pools
       if ( spinup_state .eq. 1) then
          ! speed up radioactive decay by the same factor as decomposition so tat SOM ages prematurely in all respects
          spinup_term = spinup_factor(l) 
       else
          spinup_term = 1.
       endif
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             decomp_cpools_vr(c,j,l) = decomp_cpools_vr(c,j,l) * (1.D0 - decay_const * spinup_term * dt)
          end do
       end do
    end do ! end of columns loop

 
    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       cpool(p)              = cpool(p)               * (1.D0 - decay_const * dt)
       xsmrpool(p)           = xsmrpool(p)            * (1.D0 - decay_const * dt)
       deadcrootc(p)         = deadcrootc(p)          * (1.D0 - decay_const * dt)
       deadcrootc_storage(p) = deadcrootc_storage(p)  * (1.D0 - decay_const * dt)
       deadcrootc_xfer(p)    = deadcrootc_xfer(p)     * (1.D0 - decay_const * dt)
       deadstemc(p)          = deadstemc(p)           * (1.D0 - decay_const * dt)
       deadstemc_storage(p)  = deadstemc_storage(p)   * (1.D0 - decay_const * dt)
       deadstemc_xfer(p)     = deadstemc_xfer(p)      * (1.D0 - decay_const * dt)
       frootc(p)             = frootc(p)              * (1.D0 - decay_const * dt)
       frootc_storage(p)     = frootc_storage(p)      * (1.D0 - decay_const * dt)
       frootc_xfer(p)        = frootc_xfer(p)         * (1.D0 - decay_const * dt)
       gresp_storage(p)      = gresp_storage(p)       * (1.D0 - decay_const * dt)
       gresp_xfer(p)         = gresp_xfer(p)          * (1.D0 - decay_const * dt)
       leafc(p)              = leafc(p)               * (1.D0 - decay_const * dt)
       leafc_storage(p)      = leafc_storage(p)       * (1.D0 - decay_const * dt)
       leafc_xfer(p)         = leafc_xfer(p)          * (1.D0 - decay_const * dt)
       livecrootc(p)         = livecrootc(p)          * (1.D0 - decay_const * dt)
       livecrootc_storage(p) = livecrootc_storage(p)  * (1.D0 - decay_const * dt)
       livecrootc_xfer(p)    = livecrootc_xfer(p)     * (1.D0 - decay_const * dt)
       livestemc(p)          = livestemc(p)           * (1.D0 - decay_const * dt)
       livestemc_storage(p)  = livestemc_storage(p)   * (1.D0 - decay_const * dt)
       livestemc_xfer(p)     = livestemc_xfer(p)      * (1.D0 - decay_const * dt)
       pft_ctrunc(p)         = pft_ctrunc(p)          * (1.D0 - decay_const * dt)
    end do

 end subroutine C14Decay


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C14BombSpike
!
! !INTERFACE:
subroutine C14BombSpike(num_soilp, filter_soilp)
!
! !DESCRIPTION:
! for transient pulse simulation, impose a simplified bomb spike
!
! !USES:
   use mod_clm_type
   use mod_clm_time_manager, only: get_curr_date,get_days_per_year
   use mod_clm_varcon  , only : c14ratio,secspday
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts


! !OTHER LOCAL VARIABLES:
   integer :: yr, mon, day, tod, offset
   real(rk8) :: dateyear
   real(rk8), pointer :: rc14_atm(:)             !C14O2/C12O2 in atmosphere
   real(rk8) :: delc14o2_atm 
   real(rk8) :: days_per_year                    ! days per year
   integer :: fp, p, nt
   integer :: ind_below
   integer :: ntim_atm_ts
   real(rk8) :: twt_1, twt_2                     ! weighting fractions for interpolating

   rc14_atm       => clm3%g%l%c%p%pepv%rc14_atm

   if ( use_c14_bombspike ) then

      ! get current date
      call get_curr_date(yr, mon, day, tod, offset)
      days_per_year = get_days_per_year()
      dateyear = real(yr) + real(mon)/12.D0 + real(day)/days_per_year + real(tod)/(secspday*days_per_year)

      ! find points in atm timeseries to interpolate between
      ntim_atm_ts = size(atm_c14file_time)
      ind_below = 0
      do nt = 1, ntim_atm_ts
         if (dateyear .ge. atm_c14file_time(nt) ) then
            ind_below = ind_below+1
         endif
      end do

      ! interpolate between nearest two points in atm c14 timeseries
      if (ind_below .eq. 0 ) then 
         delc14o2_atm = atm_delta_c14(1)
      elseif (ind_below .eq. ntim_atm_ts ) then
         delc14o2_atm = atm_delta_c14(ntim_atm_ts)
      else
         twt_2 = min(1.D0, max(0.D0,(dateyear-atm_c14file_time(ind_below))/ &
               (atm_c14file_time(ind_below+1)-atm_c14file_time(ind_below))))
         twt_1 = 1.D0 - twt_2
         delc14o2_atm = atm_delta_c14(ind_below) * twt_1 +  atm_delta_c14(ind_below+1) * twt_2
      endif

      ! change delta units to ratio, put on pft loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)   
         rc14_atm(p) = (delc14o2_atm * 1.D-3 + 1.D0) * c14ratio
      end do
      
   else
      ! for constant 14c concentration
      ! pft loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)   
         rc14_atm(p) = c14ratio
      end do
   endif

end subroutine C14BombSpike


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C14_init_BombSpike
!
! !INTERFACE:
subroutine C14_init_BombSpike()
!
! !DESCRIPTION:
! read netcdf file containing a timeseries of atmospheric delta C14 values; save in module-level array 
!
! !USES:
!
! !ARGUMENTS:
   implicit none

! !OTHER LOCAL VARIABLES:
   integer :: ncid            ! netcdf id
   integer :: ntim            ! number of input data time samples
   integer :: t
   
   if ( use_c14_bombspike ) then
      
      if ( myid == italk ) then
         write(stdout, *) 'C14_init_BombSpike: preparing to open file:'
         write(stdout, *) trim(atm_c14_filename)
      endif

      call openfile_withname(atm_c14_filename,ncid)
      call ncd_inqdim(ncid,'time',ntim)
      
      !! allocate arrays based on size of netcdf timeseries
      allocate(atm_c14file_time(ntim))
      allocate(atm_delta_c14(ntim))

      call read_var1d_static(ncid,'time',ntim,atm_c14file_time)
      call read_var1d_static(ncid,'atm_delta_c14',ntim,atm_delta_c14)
      
      call closefile(ncid)

      ! check to make sure that time dimension is well behaved
      do t = 2, ntim
         if ( atm_c14file_time(t) - atm_c14file_time(t-1) .le. 0.D0 ) then
            write(stderr, *) 'C14_init_BombSpike: error.  time axis must be monotonically increasing'
            call fatal(__FILE__,__LINE__,'clm now stopping')
         endif
      end do
   endif
   
end subroutine C14_init_BombSpike

#endif

end module mod_clm_cnc14decay

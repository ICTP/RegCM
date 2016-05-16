module mod_clm_soillittverttransp
#ifdef CN
  !
  ! calculate vertical mixing of all decomposing C and N pools
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_stdio
  use mod_runparams , only : dtsrf
  use mod_clm_type
  use mod_clm_varctl , only : use_c13 , use_c14 , spinup_state
  use mod_clm_varpar , only : nlevdecomp , ndecomp_pools , nlevdecomp_full
  use mod_clm_varcon , only : secspday , zsoi , dzsoi_decomp , zisoi
  use mod_clm_tridiagonal , only : Tridiagonal

  implicit none

  private

  save

  public :: CNSoilLittVertTransp

  ! [m^2/sec] = 1 cm^2 / yr
  real(rk8), public :: som_diffus = 1e-4_rk8 / (secspday * 365._rk8)
  real(rk8), public :: som_adv_flux =  0._rk8
  ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr
  real(rk8), public :: cryoturb_diffusion_k =  5e-4_rk8 / (secspday * 365._rk8)
  ! (m) this is the maximum depth of cryoturbation
  real(rk8), public :: max_depth_cryoturb = 3._rk8
  ! (m) maximum active layer thickness for cryoturbation to occur
  real(rk8), public :: max_altdepth_cryoturbation = 2._rk8

  contains
  !
  ! Calculate vertical mixing of soil and litter pools.
  ! Also reconcile sources and sinks of these pools
  ! calculated in the CStateUpdate1 and NStateUpdate1 subroutines.
  ! Advection-diffusion code based on algorithm in Patankar (1980)
  ! Initial code by C. Koven and W. Riley
  !

  subroutine CNSoilLittVertTransp(lbc, ubc, num_soilc, filter_soilc)
    implicit none
    integer(ik4) , intent(in) :: num_soilc ! number of soil columns in filter
    integer(ik4) , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer(ik4) , intent(in) :: lbc , ubc ! column-index bounds
    ! SOM advective flux (m/s)
    real(rk8), pointer :: som_adv_coef(:,:)
    ! SOM diffusivity due to bio/cryo-turbation (m2/s)
    real(rk8), pointer :: som_diffus_coef(:,:)
    ! pointer for concentration state variable being transported
    real(rk8), pointer :: conc_ptr(:,:,:)
    ! pointer for source term
    real(rk8), pointer :: source(:,:,:)
    ! TRUE => pool is a cwd pool
    logical, pointer  :: is_cwd(:)
#ifdef VERTSOILC
    ! difference between nodes
    real(rk8) :: dz_node(1:nlevdecomp+1)
    ! diffusivity (m2/s)  (includes spinup correction, if any)
    real(rk8) :: diffus (lbc:ubc,1:nlevdecomp+1)
    ! Weights for calculating harmonic mean of diffusivity
    real(rk8) :: w_m1, w_p1
    ! Harmonic mean of diffusivity
    real(rk8) :: d_m1, d_p1
    ! diffusivity/delta_z for next j (set to zero for no diffusion)
    real(rk8) :: d_p1_zp1(lbc:ubc,1:nlevdecomp+1)
    ! diffusivity/delta_z for previous j (set to zero for no diffusion)
    real(rk8) :: d_m1_zm1 (lbc:ubc,1:nlevdecomp+1)
    ! advective flux (m/s)  (includes spinup correction, if any)
    real(rk8) :: adv_flux(lbc:ubc,1:nlevdecomp+1)
    ! "a" vector for tridiagonal matrix
    real(rk8) :: a_tri(lbc:ubc,0:nlevdecomp+1)
    ! "b" vector for tridiagonal matrix
    real(rk8) :: b_tri(lbc:ubc,0:nlevdecomp+1)
    ! "c" vector for tridiagonal matrix
    real(rk8) :: c_tri(lbc:ubc,0:nlevdecomp+1)
    ! "r" vector for tridiagonal solution
    real(rk8) :: r_tri(lbc:ubc,0:nlevdecomp+1)
    real(rk8) :: a_p_0
    real(rk8) :: deficit
    real(rk8) :: conc_trcr(lbc:ubc,0:nlevdecomp+1)
    real(rk8) :: f_p1(lbc:ubc,1:nlevdecomp+1)  ! water flux for next j
    real(rk8) :: f_m1(lbc:ubc,1:nlevdecomp+1)  ! water flux for previous j
    real(rk8) :: pe_p1(lbc:ubc,1:nlevdecomp+1) ! Peclet # for next j
    real(rk8) :: pe_m1(lbc:ubc,1:nlevdecomp+1) ! Peclet # for previous j
    integer(ik4) :: jtop(lbc:ubc)     ! top level at each column
    integer(ik4) :: s ! indices
#endif
    ! pointer to store the vertical tendency
    ! (gain/loss due to vertical transport)
    real(rk8), pointer :: trcr_tendency_ptr(:,:,:)
    ! maximum annual depth of thaw
    real(rk8), pointer :: altmax(:)
    ! prior year maximum annual depth of thaw
    real(rk8), pointer :: altmax_lastyear(:)

    integer(ik4) :: ntype
    integer(ik4) :: i_type , fc , c , j , l ! indices
    ! spinup accelerated decomposition factor, used to accelerate
    ! transport as well
    real(rk8), pointer :: spinup_factor(:)
    ! spinup accelerated decomposition factor, used to accelerate
    ! transport as well
    real(rk8) :: spinup_term
    real(rk8) :: epsilon  ! small number

    is_cwd             => decomp_cascade_con%is_cwd
    spinup_factor      => decomp_cascade_con%spinup_factor
    altmax             => clm3%g%l%c%cps%altmax
    altmax_lastyear    => clm3%g%l%c%cps%altmax_lastyear
    som_adv_coef       => clm3%g%l%c%cps%som_adv_coef
    som_diffus_coef    => clm3%g%l%c%cps%som_diffus_coef

    ntype = 2
    if ( use_c13 ) then
      ntype = ntype+1
    end if
    if ( use_c14 ) then
      ntype = ntype+1
    end if
    spinup_term = 1._rk8
    epsilon = 1.e-30

#ifdef VERTSOILC
    !------ first get diffusivity / advection terms -------!
    ! use different mixing rates for bioturbation and cryoturbation, with
    ! fixed bioturbation and cryoturbation set to a maximum depth
    do fc = 1, num_soilc
      c = filter_soilc (fc)
      if  ( (max(altmax(c),altmax_lastyear(c)) <= max_altdepth_cryoturbation) &
              .and. ( max(altmax(c), altmax_lastyear(c)) > 0._rk8) ) then
        ! use mixing profile modified slightly from Koven et al. (2009):
        !   constant through active layer, linear decrease from base of
        !   active layer to zero at a fixed depth
        do j = 1,nlevdecomp+1
          if ( zisoi(j) < max(altmax(c), altmax_lastyear(c)) ) then
            som_diffus_coef(c,j) = cryoturb_diffusion_k
            som_adv_coef(c,j) = 0._rk8
          else
            ! go linearly to zero between ALT and max_depth_cryoturb
            som_diffus_coef(c,j) = max(cryoturb_diffusion_k * &
              (1._rk8-( zisoi(j)-max(altmax(c),altmax_lastyear(c)) ) / &
              (max_depth_cryoturb-max(altmax(c),altmax_lastyear(c)))), 0._rk8)
            som_adv_coef(c,j) = 0._rk8
          end if
        end do
      else if (  max(altmax(c), altmax_lastyear(c)) > 0._rk8 ) then
        ! constant advection, constant diffusion
        do j = 1,nlevdecomp+1
          som_adv_coef(c,j) = som_adv_flux
          som_diffus_coef(c,j) = som_diffus
        end do
      else
        ! completely frozen soils--no mixing
        do j = 1,nlevdecomp+1
           som_adv_coef(c,j) = 0._rk8
           som_diffus_coef(c,j) = 0._rk8
        end do
      end if
    end do

    ! Set the distance between the node and the one ABOVE it
    dz_node(1) = zsoi(1)
    do j = 2,nlevdecomp+1
      dz_node(j)= zsoi(j) - zsoi(j-1)
    end do

#endif

    !------ loop over litter/som types
    do i_type = 1, ntype
      select case (i_type)
        case (1)  ! C
          conc_ptr => clm3%g%l%c%ccs%decomp_cpools_vr
          source    => clm3%g%l%c%ccf%decomp_cpools_sourcesink
          trcr_tendency_ptr => clm3%g%l%c%ccf%decomp_cpools_transport_tendency
        case (2)  ! N
          conc_ptr => clm3%g%l%c%cns%decomp_npools_vr
          source    => clm3%g%l%c%cnf%decomp_npools_sourcesink
          trcr_tendency_ptr => clm3%g%l%c%cnf%decomp_npools_transport_tendency
        case (3)
          if ( use_c13 ) then
             ! C13
             conc_ptr => clm3%g%l%c%cc13s%decomp_cpools_vr
             source    => clm3%g%l%c%cc13f%decomp_cpools_sourcesink
             trcr_tendency_ptr => &
                     clm3%g%l%c%cc13f%decomp_cpools_transport_tendency
          else
             ! C14
             conc_ptr => clm3%g%l%c%cc14s%decomp_cpools_vr
             source    => clm3%g%l%c%cc14f%decomp_cpools_sourcesink
             trcr_tendency_ptr => &
                     clm3%g%l%c%cc14f%decomp_cpools_transport_tendency
          end if
        case (4)
          if ( use_c14 .and. use_c13 ) then
            ! C14
            conc_ptr => clm3%g%l%c%cc14s%decomp_cpools_vr
            source    => clm3%g%l%c%cc14f%decomp_cpools_sourcesink
            trcr_tendency_ptr => &
                    clm3%g%l%c%cc14f%decomp_cpools_transport_tendency
          else
            write(stderr,*) &
                    'error.  ncase = 4, but c13 and c14 not both enabled.'
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
      end select

#ifdef VERTSOILC
      do s = 1, ndecomp_pools
        if ( spinup_state == 1 ) then
          ! increase transport (both advection and diffusion)
          ! by the same factor as accelerated decomposition for a given pool
          spinup_term = spinup_factor(s)
        else
          spinup_term = 1.
        end if
        if ( .not. is_cwd(s) ) then
          do j = 1,nlevdecomp+1
            do fc = 1, num_soilc
              c = filter_soilc (fc)
              if ( abs(som_adv_coef(c,j)) * spinup_term < epsilon ) then
                adv_flux(c,j) = epsilon
              else
                adv_flux(c,j) = som_adv_coef(c,j) * spinup_term
              end if
              if ( abs(som_diffus_coef(c,j)) * spinup_term < epsilon ) then
                diffus(c,j) = epsilon
              else
                diffus(c,j) = som_diffus_coef(c,j) * spinup_term
              end if
            end do
          end do
          ! Set Pe (Peclet #) and D/dz throughout column
          do fc = 1, num_soilc ! dummy terms here
            c = filter_soilc (fc)
            conc_trcr(c,0) = 0._rk8
            conc_trcr(c,nlevdecomp+1) = 0._rk8
          end do
          do j = 1,nlevdecomp+1
            do fc = 1, num_soilc
              c = filter_soilc (fc)

              conc_trcr(c,j) = conc_ptr(c,j,s)
              ! dz_tracer below is the difference between gridcell edges
              !  (dzsoi_decomp)
              ! dz_node_tracer is difference between cell centers

              ! Calculate the D and F terms in the Patankar algorithm
              if (j == 1) then
                d_m1_zm1(c,j) = 0._rk8
                w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
                if ( diffus(c,j+1) > 0._rk8 .and. diffus(c,j) > 0._rk8) then
                  d_p1 = 1._rk8 / ((1._rk8 - w_p1) / diffus(c,j) + &
                          w_p1 / diffus(c,j+1)) ! Harmonic mean of diffus
                else
                  d_p1 = 0._rk8
                end if
                d_p1_zp1(c,j) = d_p1 / dz_node(j+1)
                f_m1(c,j) = adv_flux(c,j)  ! Include infiltration here
                f_p1(c,j) = adv_flux(c,j+1)
                pe_m1(c,j) = 0._rk8
                pe_p1(c,j) = f_p1(c,j) / d_p1_zp1(c,j) ! Peclet #
              else if (j == nlevdecomp+1) then
                ! At the bottom, assume no gradient in d_z
                ! (i.e., they're the same)
                w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
                if ( diffus(c,j) > 0._rk8 .and. diffus(c,j-1) > 0._rk8) then
                  d_m1 = 1._rk8 / ((1._rk8 - w_m1) / diffus(c,j) + &
                          w_m1 / diffus(c,j-1)) ! Harmonic mean of diffus
                else
                  d_m1 = 0._rk8
                end if
                d_m1_zm1(c,j) = d_m1 / dz_node(j)
                d_p1_zp1(c,j) = d_m1_zm1(c,j) ! Set to be the same
                f_m1(c,j) = adv_flux(c,j)
                !f_p1(c,j) = adv_flux(c,j+1)
                f_p1(c,j) = 0._rk8
                pe_m1(c,j) = f_m1(c,j) / d_m1_zm1(c,j) ! Peclet #
                pe_p1(c,j) = f_p1(c,j) / d_p1_zp1(c,j) ! Peclet #
              else
                ! Use distance from j-1 node to interface with j divided
                ! by distance between nodes
                w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
                if ( diffus(c,j-1) > 0._rk8 .and. diffus(c,j) > 0._rk8) then
                  d_m1 = 1._rk8 / ((1._rk8 - w_m1) / diffus(c,j) + &
                          w_m1 / diffus(c,j-1)) ! Harmonic mean of diffus
                else
                  d_m1 = 0._rk8
                end if
                w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
                if ( diffus(c,j+1) > 0._rk8 .and. diffus(c,j) > 0._rk8) then
                  d_p1 = 1._rk8 / ((1._rk8 - w_p1) / diffus(c,j) + &
                          w_p1 / diffus(c,j+1)) ! Harmonic mean of diffus
                else
                  ! Arithmetic mean of diffus
                  d_p1 = (1._rk8 - w_m1) * diffus(c,j) + w_p1 * diffus(c,j+1)
                end if
                d_m1_zm1(c,j) = d_m1 / dz_node(j)
                d_p1_zp1(c,j) = d_p1 / dz_node(j+1)
                f_m1(c,j) = adv_flux(c,j)
                f_p1(c,j) = adv_flux(c,j+1)
                pe_m1(c,j) = f_m1(c,j) / d_m1_zm1(c,j) ! Peclet #
                pe_p1(c,j) = f_p1(c,j) / d_p1_zp1(c,j) ! Peclet #
              end if
            end do ! fc
          end do ! j; nlevdecomp
          ! Calculate the tridiagonal coefficients
          do j = 0,nlevdecomp +1
            do fc = 1, num_soilc
              c = filter_soilc (fc)
              ! g = cgridcell(c)

              if (j > 0 .and. j < nlevdecomp+1) then
                a_p_0 =  dzsoi_decomp(j) / dtsrf
              end if

              if (j == 0) then ! top layer (atmosphere)
                a_tri(c,j) = 0._rk8
                b_tri(c,j) = 1._rk8
                c_tri(c,j) = -1._rk8
                r_tri(c,j) = 0._rk8
              else if (j == 1) then
                a_tri(c,j) = -(d_m1_zm1(c,j) * aaa(pe_m1(c,j)) + &
                        max( f_m1(c,j), 0._rk8)) ! Eqn 5.47 Patankar
                c_tri(c,j) = -(d_p1_zp1(c,j) * aaa(pe_p1(c,j)) + &
                        max(-f_p1(c,j), 0._rk8))
                b_tri(c,j) = -a_tri(c,j) - c_tri(c,j) + a_p_0
                r_tri(c,j) = source(c,j,s) * dzsoi_decomp(j) /dtsrf + &
                        (a_p_0 - adv_flux(c,j)) * conc_trcr(c,j)
              else if (j < nlevdecomp+1) then
                a_tri(c,j) = -(d_m1_zm1(c,j) * aaa(pe_m1(c,j)) + &
                        max( f_m1(c,j), 0._rk8)) ! Eqn 5.47 Patankar
                c_tri(c,j) = -(d_p1_zp1(c,j) * aaa(pe_p1(c,j)) + &
                        max(-f_p1(c,j), 0._rk8))
                b_tri(c,j) = -a_tri(c,j) - c_tri(c,j) + a_p_0
                r_tri(c,j) = source(c,j,s) * dzsoi_decomp(j) /dtsrf + &
                        a_p_0 * conc_trcr(c,j)
              else ! j==nlevdecomp+1; 0 concentration gradient at bottom
                a_tri(c,j) = -1._rk8
                b_tri(c,j) = 1._rk8
                c_tri(c,j) = 0._rk8
                r_tri(c,j) = 0._rk8
              end if
            end do ! fc; column
          end do ! j; nlevdecomp

          do fc = 1, num_soilc
            c = filter_soilc (fc)
            jtop(c) = 0
          end do
          ! subtract initial concentration and source terms for
          ! tendency calculation
          do fc = 1 , num_soilc
            c = filter_soilc (fc)
            do j = 1 , nlevdecomp
              trcr_tendency_ptr(c,j,s) = 0.-(conc_trcr(c,j) + source(c,j,s))
            end do
          end do
          ! Solve for the concentration profile for this time step
          call Tridiagonal(lbc, ubc, 0, nlevdecomp+1, jtop, num_soilc, &
                  filter_soilc, a_tri, b_tri, c_tri, r_tri, &
                  conc_trcr(lbc:ubc,0:nlevdecomp+1))
          ! add post-transport concentration to calculate tendency term
          do fc = 1 , num_soilc
            c = filter_soilc (fc)
            do j = 1 , nlevdecomp
              trcr_tendency_ptr(c,j,s) = trcr_tendency_ptr(c,j,s) + &
                      conc_trcr(c,j)
              trcr_tendency_ptr(c,j,s) = trcr_tendency_ptr(c,j,s) / dtsrf
            end do
          end do
        else
          ! for CWD pools, just add
          do j = 1 , nlevdecomp
            do fc = 1 , num_soilc
              c = filter_soilc (fc)
              conc_trcr(c,j) = conc_ptr(c,j,s) + source(c,j,s)
            end do
          end do
        end if ! not CWD
        do j = 1 , nlevdecomp
          do fc = 1 , num_soilc
            c = filter_soilc (fc)
            conc_ptr(c,j,s) = conc_trcr(c,j)
          end do
        end do
      end do ! s (pool loop)
#else
      ! for single level case, no transport
      ! just update the fluxes calculated in the StateUpdate1 subroutines
      do l = 1 , ndecomp_pools
        do j = 1 , nlevdecomp
          do fc = 1 , num_soilc
            c = filter_soilc (fc)
            conc_ptr(c,j,l) = conc_ptr(c,j,l) + source(c,j,l)
            trcr_tendency_ptr(c,j,l) = 0._rk8
          end do
        end do
      end do
#endif
    end do  ! i_type

    contains

      ! A function from Patankar, Table 5.2, pg 95
      pure real(rk8) function aaa(pe)
        implicit none
        real(rk8) , intent(in) :: pe
        aaa = max (0._rk8, (1._rk8 - 0.1_rk8 * abs(pe))**5)
      end function aaa

  end subroutine CNSoilLittVertTransp

#endif

end module mod_clm_soillittverttransp
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

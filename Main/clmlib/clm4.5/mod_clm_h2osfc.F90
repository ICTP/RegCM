module mod_clm_h2osfc
  !
  ! Calculate surface water hydrology
  !
  implicit none

  private

  save

  ! Calculate fraction of land surface that is wet
  public :: FracH2oSfc

  contains
  !
  ! Determine fraction of land surfaces which are submerged
  ! based on surface microtopography and surface water storage.
  !
  subroutine FracH2oSfc(lbc,ubc,num_h2osfc,filter_h2osfc,frac_h2osfc,no_update)
    use mod_intkinds
    use mod_realkinds
    use mod_clm_type
    use mod_clm_varcon , only : rpi , istsoil, istcrop
    implicit none
    integer(ik4) , intent(in) :: lbc, ubc   ! column bounds
    ! number of column points in column filter
    integer(ik4) , intent(in) :: num_h2osfc
    integer(ik4) , intent(in) :: filter_h2osfc(ubc-lbc+1)  ! column filter
    ! fractional surface water (mm)
    real(rk8), intent(inout) :: frac_h2osfc(lbc:ubc)
    ! flag to make calculation w/o updating variables
    integer(ik4) , intent(in), optional :: no_update

    real(rk8), pointer :: h2osfc(:)         ! surface water (mm)

    ! microtopography pdf sigma (m)
    real(rk8), pointer :: micro_sigma(:)
    ! fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: frac_sno(:)
    ! eff. fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: frac_sno_eff(:)
    integer(ik4) , pointer :: snl(:)      ! minus number of snow layers
    real(rk8), pointer :: h2osno(:)       ! snow water (mm H2O)
    real(rk8), pointer :: h2osoi_liq(:,:) ! liquid water (col,lyr) [kg/m2]
    real(rk8), pointer :: topo_slope(:)   ! topographic slope
    real(rk8), pointer :: topo_ndx(:)     ! topographic slope
    integer(ik4) , pointer :: ltype(:)    ! landunit type
    integer(ik4) , pointer :: clandunit(:)  ! columns's landunit

    integer(ik4)  :: c,f,l     ! indices
    real(rk8):: d,fd,dfdd      ! temporary variable for frac_h2oscs iteration
    real(rk8):: sigma          ! microtopography pdf sigma in mm
    real(rk8):: min_h2osfc

    ! Assign local pointers to derived subtypes components (column-level)

    h2osoi_liq          => clm3%g%l%c%cws%h2osoi_liq
    h2osfc              => clm3%g%l%c%cws%h2osfc
    micro_sigma         => clm3%g%l%c%cps%micro_sigma
    topo_slope          => clm3%g%l%c%cps%topo_slope
    topo_ndx            => clm3%g%l%c%cps%topo_ndx
    ltype               => clm3%g%l%itype
    clandunit           => clm3%g%l%c%landunit

    frac_sno            => clm3%g%l%c%cps%frac_sno
    frac_sno_eff        => clm3%g%l%c%cps%frac_sno_eff
    snl                 => clm3%g%l%c%cps%snl
    h2osno              => clm3%g%l%c%cws%h2osno

    ! arbitrary lower limit on h2osfc for safer numerics...
    min_h2osfc = 1.D-8

    do f = 1 , num_h2osfc
      c = filter_h2osfc(f)
      l = clandunit(c)
      ! h2osfc only calculated for soil vegetated land units
      if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then

        !  Use newton-raphson method to iteratively determine frac_h20sfc
        !  based on amount of surface water storage (h2osfc) and
        !  microtopography variability (micro_sigma)
        if ( h2osfc(c) > min_h2osfc ) then
          ! a cutoff is needed for numerical reasons...
          ! (nonconvergence after 5 iterations)
          d = 0.0D0
          sigma = 1.0D3 * micro_sigma(c) ! convert to mm
          do l = 1 , 10
            fd = 0.5D0*d*(1.0D0+erf(d/(sigma*sqrt(2.0D0)))) + &
                    sigma/sqrt(2.0D0*rpi)*exp(-d**2/(2.0D0*sigma**2)) - &
                    h2osfc(c)
            dfdd = 0.5D0*(1.0D0+erf(d/(sigma*sqrt(2.0D0))))
            d = d - fd/dfdd
          end do
          !--  update the submerged areal fraction using the new d value
          frac_h2osfc(c) = 0.5D0*(1.0D0+erf(d/(sigma*sqrt(2.0D0))))
        else
          frac_h2osfc(c) = 0.D0
          h2osoi_liq(c,1) = h2osoi_liq(c,1) + h2osfc(c)
          h2osfc(c) = 0.D0
        end if

        if ( .not. present(no_update) ) then

          ! adjust fh2o, fsno when sum is greater than zero
          ! energy balance error when h2osno > 0 and snl = 0
          if ( frac_sno(c) > (1.D0 - frac_h2osfc(c)) .and. h2osno(c) > 0 ) then
            if ( frac_h2osfc(c) > 0.01D0 ) then
              frac_h2osfc(c) = max(1.0D0 - frac_sno(c),0.01D0)
              frac_sno(c) = 1.0D0 - frac_h2osfc(c)
            else
              frac_sno(c) = 1.0D0 - frac_h2osfc(c)
            end if
            frac_sno_eff(c)=frac_sno(c)
          end if
        end if ! end of no_update construct
      else !if landunit not istsoil/istcrop, set frac_h2osfc to zero
        frac_h2osfc(c) = 0.D0
      end if
    end do
  end subroutine FracH2oSfc

end module mod_clm_h2osfc

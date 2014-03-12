module mod_clm_activelayer
  !
  ! Module holding routines for calculation of active layer dynamics
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_runparams
  use mod_date
  use mod_clm_type , only : clm3
  use mod_clm_varpar , only : nlevgrnd
  use mod_clm_varcon , only : zsoi
  
  implicit none

  private

  save

  public:: alt_calc
  
  contains

  ! Define active layer thickness similarly to frost_table, except set as
  ! deepest thawed layer and define on nlevgrnd
  ! also update annual maxima, and keep track of prior year for
  ! rooting memory
  subroutine alt_calc(lbc,ubc,num_soilc,filter_soilc)
    implicit none
    integer(ik4) , intent(in) :: lbc , ubc  ! column bounds
    integer(ik4) , intent(in) :: num_soilc  ! number of soil columns in filter
    ! filter for soil columns
    integer(ik4) , intent(in) , dimension(:) :: filter_soilc
    !
    ! local pointers to implicit in scalars
    !
    ! column level
    ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8) , pointer , dimension(:,:) :: t_soisno
    ! current depth of thaw
    real(rk8) , pointer , dimension(:) :: alt
    ! maximum annual depth of thaw
    real(rk8) , pointer , dimension(:) :: altmax
    ! prior year maximum annual depth of thaw
    real(rk8) , pointer , dimension(:) :: altmax_lastyear
    ! current depth of thaw
    integer(ik4) , pointer , dimension(:) :: alt_indx
    ! maximum annual depth of thaw
    integer(ik4) , pointer , dimension(:) :: altmax_indx
    ! prior year maximum annual depth of thaw
    integer(ik4) , pointer , dimension(:) :: altmax_lastyear_indx
    ! gridcell latitude (radians)
    real(rk8) , pointer , dimension(:) :: lat
    ! gridcell index of column
    integer(ik4) , pointer , dimension(:) :: cgridcell
 
    integer(ik4) :: c , j , fc , g  ! counters
    integer(ik4) :: alt_ind         ! index of base of activel layer
    integer(ik4) :: year            ! year (0, ...) for nstep+1
    integer(ik4) :: mon             ! month (1, ..., 12) for nstep+1
    integer(ik4) :: day             ! day of month (1, ..., 31) for nstep+1
    integer(ik4) :: sec             ! seconds into current date for nstep+1
    integer(ik4) :: dtime           ! time step length in seconds
    integer(ik4) :: k_frz           ! index of first nonfrozen soil layer
    real(rk8) :: t1 , t2 , z1 , z2  ! temporary variables
    ! used to break loop when first unfrozen layer reached
    logical :: found_thawlayer

    ! Assign local pointers to derived type arrays
    t_soisno              => clm3%g%l%c%ces%t_soisno
    alt                   => clm3%g%l%c%cps%alt
    altmax                => clm3%g%l%c%cps%altmax
    altmax_lastyear       => clm3%g%l%c%cps%altmax_lastyear
    alt_indx              => clm3%g%l%c%cps%alt_indx
    altmax_indx           => clm3%g%l%c%cps%altmax_indx
    altmax_lastyear_indx  => clm3%g%l%c%cps%altmax_lastyear_indx

    ! Assign local pointers to derived subtypes components
    ! (gridcell-level and mapping)
    lat             => clm3%g%lat
    cgridcell       => clm3%g%l%c%gridcell

    ! on a set annual timestep, update annual maxima
    ! make this 1 January for NH columns, 1 July for SH columns
    if ( date_is(idatex,1,1) .and. time_is(idatex,1,dtsrf) ) then
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        g = cgridcell(c)
        if ( lat(g) > d_zero ) then 
          altmax_lastyear(c) = altmax(c)
          altmax_lastyear_indx(c) = altmax_indx(c)
          altmax(c) = d_zero
          altmax_indx(c) = 0
        end if
      end do
    end if
    if ( date_is(idatex,7,1) .and. time_is(idatex,1,dtsrf) ) then
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = cgridcell(c)
        if ( lat(g) <= d_zero ) then 
          altmax_lastyear(c) = altmax(c)
          altmax_lastyear_indx(c) = altmax_indx(c)
          altmax(c) = d_zero
          altmax_indx(c) = 0
        end if
      end do
    end if

    do fc = 1 , num_soilc
      c = filter_soilc(fc)
      ! calculate alt for a given timestep
      ! start from base of soil and search upwards for first thawed layer.
      ! note that this will put talik in with active layer
      ! a different way of doing this could be to keep track of how long a
      ! given layer has ben frozen for, and define ALT as the first layer
      ! that has been frozen for less than 2 years.
      if ( t_soisno(c,nlevgrnd) > tzero ) then
        alt(c) = zsoi(nlevgrnd)
        alt_indx(c) = nlevgrnd
      else
        k_frz = 0
        found_thawlayer = .false.
        do j = nlevgrnd-1 , 1 , -1
          if ( ( t_soisno(c,j) > tzero ) .and. .not. found_thawlayer ) then
            k_frz = j
            found_thawlayer = .true.
          end if
        end do
        if ( k_frz > 0 ) then
          ! define active layer as the depth at which the linearly
          ! interpolated temperature line intersects with zero
          z1 = zsoi(k_frz)
          z2 = zsoi(k_frz+1)
          t1 = t_soisno(c,k_frz)
          t2 = t_soisno(c,k_frz+1)
          alt(c) = z1 + (t1-tzero)*(z2-z1)/(t1-t2)
          alt_indx(c) = k_frz
        else
          alt(c) = 0.0D0
          alt_indx(c) = 0
        end if
      end if

      ! if appropriate, update maximum annual active layer thickness
      if ( alt(c) > altmax(c) ) then
        altmax(c) = alt(c)
        altmax_indx(c) = alt_indx(c)
      end if
    end do
  end subroutine alt_calc

end module mod_clm_activelayer

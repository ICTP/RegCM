module mod_clm_vicmap
  !
  ! Performs  the mapping from CLM layers to VIC layers
  ! Specifically, 10 (or 23 when more_vertlayers == .true.)
  ! CLM hydrologically active soil layers are mapped to three VIC layers
  ! by assigning the first nlvic(1) layers to VIC layer 1
  !              the next nlvic(2) layers  to VIC alyer 2
  !              and the remaining to VIC layer 3
  !
  use mod_intkinds
  use mod_realkinds

  implicit none

  private

  save

#if (defined VICHYDRO)
  public  :: initCLMVICMap   ! map layer/node fractions
  public  :: CLMVICMap       ! map from VIC to CLM layers

  contains
  !
  ! This subroutine calculates mapping between CLM and VIC layers
  ! added by AWang
  ! modified by M.Huang for CLM4
  !
  subroutine initCLMVICMap(c)
    use mod_clm_type
    use mod_clm_varcon , only : denh2o, denice, pondmx
    use mod_clm_varpar , only : nlevsoi, nlayer, nlayert, nlevgrnd
    implicit none
    integer(ik4) , intent(in)  :: c

    real(rk8), pointer :: dz(:,:)       ! layer depth (m)
    real(rk8), pointer :: zi(:,:)       ! interface level below a "z" level (m)
    real(rk8), pointer :: z(:,:)        ! layer thickness (m)
    real(rk8), pointer :: depth(:,:)    ! layer depth of VIC (m)
    ! fraction of VIC layers in clm layers
    real(rk8), pointer :: vic_clm_fract(:,:,:)
    ! sum of fraction for each layer
    real(rk8) :: sum_frac(1:nlayer)
    real(rk8) :: deltal(1:nlayer+1) ! temporary
    real(rk8) :: zsum               ! temporary
    real(rk8) :: lsum               ! temporary
    real(rk8) :: temp               ! temporary

    integer(ik4) :: i, j, fc

    ! note: in CLM h2osoil_liq unit is kg/m2, in VIC moist is mm
    ! h2osoi_ice is actually water equavlent ice content.
    ! Assign local pointers to derived subtypes components (column-level)

    dz         => clm3%g%l%c%cps%dz
    zi         => clm3%g%l%c%cps%zi
    z          => clm3%g%l%c%cps%z
    depth      => clm3%g%l%c%cps%depth
    vic_clm_fract => clm3%g%l%c%cps%vic_clm_fract

    !  set fraction of VIC layer in each CLM layer

    lsum = 0._rk8
    do i = 1, nlayer
      deltal(i) = depth(c,i)
    end do
    do i = 1, nlayer
      zsum = 0._rk8
      sum_frac(i) = 0._rk8
      do j = 1, nlevsoi
        if ( (zsum < lsum) .and. (zsum + dz(c,j) >= lsum ) )  then
          temp = linear_interp(lsum, zsum, zsum + dz(c,j), 0._rk8, 1._rk8)
          vic_clm_fract(c,i,j) = 1._rk8 - temp
          if ( lsum + deltal(i) < zsum + dz(c,j) ) then
            temp = linear_interp(lsum + deltal(i), zsum, &
                    zsum + dz(c,j), 1._rk8, 0._rk8)
            vic_clm_fract(c,i,j) = vic_clm_fract(c,i,j) - temp
          end if
        else if ( (zsum < lsum + deltal(i)) .and. &
                  (zsum + dz(c,j) >= lsum + deltal(i)) ) then
          temp = linear_interp(lsum + deltal(i), zsum, &
                  zsum + dz(c,j), 0._rk8, 1._rk8)
          vic_clm_fract(c,i,j) = temp
          if ( zsum <= lsum ) then
            temp = linear_interp(lsum, zsum, zsum + dz(c,j), 0._rk8, 1._rk8)
            vic_clm_fract(c,i,j) = vic_clm_fract(c,i,j) - temp
          end if
        else if ( (zsum >= lsum .and. &
                  zsum + dz(c,j) <= lsum + deltal(i)) )  then
          vic_clm_fract(c,i,j) = 1._rk8
        else
          vic_clm_fract(c,i,j) = 0._rk8
        end if
        zsum = zsum + dz(c,j)
        sum_frac(i) = sum_frac(i) + vic_clm_fract(c,i,j)
      end do                           ! end CLM layer calculation
      lsum = lsum + deltal(i)
    end do                             ! end VIC layer calcultion
    contains

    ! This subroutine provides linear interpolation
    pure real(rk8) function linear_interp(x,x0,x1,y0,y1) result (y)
      implicit none
      real(rk8) , intent(in) :: x , x0 , y0 , x1 , y1
      y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
    end function linear_interp

  end subroutine initCLMVICMap
  !
  ! mapping from VIC to CLM layers, M.Huang
  !
  subroutine CLMVICMap(lbc, ubc, numf, filter)
    use mod_clm_type
    use mod_clm_varcon , only : denh2o, denice, pondmx, watmin
    use mod_clm_varpar , only : nlevsoi, nlayer, nlayert, nlevgrnd
    implicit none
    integer(ik4) , intent(in)  :: lbc, ubc ! column bounds
    ! number of column soil points in column filter
    integer(ik4) , intent(in)  :: numf
    ! column filter for soil points
    integer(ik4) , intent(in)  :: filter(ubc-lbc+1)

    real(rk8), pointer :: dz(:,:)  !layer depth (m)
    real(rk8), pointer :: zi(:,:)  !interface level below a "z" level (m)
    real(rk8), pointer :: z(:,:)   !layer thickness (m)
    real(rk8), pointer :: h2osoi_liq(:,:)  !liquid water (kg/m2)
    real(rk8), pointer :: h2osoi_ice(:,:)  !ice lens (kg/m2)
    real(rk8), pointer :: moist(:,:)       !liquid water (mm)
    real(rk8), pointer :: ice(:,:)         !ice lens (mm)
    real(rk8), pointer :: depth(:,:)       !layer depth of upper layer (m)
    real(rk8), pointer :: max_moist(:,:)   !max layer moist + ice (mm)
    !volumetric soil moisture for VIC soil layers
    real(rk8), pointer :: moist_vol(:,:)
    !volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
    real(rk8), pointer :: h2osoi_vol(:,:)
    !soil porisity (1-bulk_density/soil_density)
    real(rk8), pointer :: porosity(:,:)
    !fraction of VIC layers in each CLM layer
    real(rk8), pointer :: vic_clm_fract(:,:,:)

    real(rk8) :: ice0(1:nlayer)            ! last step ice lens (mm)  (new)
    real(rk8) :: moist0(1:nlayer)          ! last step soil water (mm)  (new)
    integer(ik4)  :: i, j, c, fc

    ! note: in CLM3 h2osoil_liq unit is kg/m2, in VIC moist is mm
    ! h2osoi_ice is actually water equavlent ice content.
    ! Assign local pointers to derived subtypes components (column-level)
    dz            => clm3%g%l%c%cps%dz
    zi            => clm3%g%l%c%cps%zi
    z             => clm3%g%l%c%cps%z
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice
    moist         => clm3%g%l%c%cws%moist
    ice           => clm3%g%l%c%cws%ice
    h2osoi_vol    => clm3%g%l%c%cws%h2osoi_vol
    moist_vol     => clm3%g%l%c%cws%moist_vol
    porosity      => clm3%g%l%c%cps%porosity
    depth         => clm3%g%l%c%cps%depth
    max_moist     => clm3%g%l%c%cps%max_moist
    vic_clm_fract => clm3%g%l%c%cps%vic_clm_fract

    ! map CLM to VIC
    do fc = 1 , numf
      c = filter(fc)
      do i = 1 , nlayer
        ice0(i) = ice(c,i)
        moist0(i) = moist(c,i)
        ice(c,i) = 0._rk8
        moist(c,i) = 0._rk8
        do j = 1, nlevsoi
          ice(c,i) = ice(c,i) + h2osoi_ice(c,j) * vic_clm_fract(c,i,j)
          moist(c,i) = moist(c,i) + h2osoi_liq(c,j) * vic_clm_fract(c,i,j)
        end do
        ice(c,i) = min((moist0(i) + ice0(i)), ice(c,i))
        ice(c,i) = max(0._rk8, ice(c,i))
        moist(c,i) =max(watmin, moist(c,i))
        moist(c,i) =min(max_moist(c,i)-ice(c,i), moist(c,i))
        moist_vol(c,i) = moist(c,i)/(depth(c,i)*denice) &
                       + ice(c,i)/(depth(c,i)*denh2o)
        moist_vol(c,i) = min(porosity(c,i), moist_vol(c,i))
        moist_vol(c,i) = max(0.01_rk8, moist_vol(c,i))
      end do

      ! hydrologic inactive layers
      ice(c, nlayer+1:nlayert) = h2osoi_ice(c, nlevsoi+1:nlevgrnd)
      moist(c, nlayer+1:nlayert) = h2osoi_liq(c, nlevsoi+1:nlevgrnd)
      moist_vol(c, nlayer+1:nlayert) = h2osoi_vol(c, nlevsoi+1:nlevgrnd)
    end do
  end subroutine CLMVICMap

#endif

end module mod_clm_vicmap
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

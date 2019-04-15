module mod_clm_dynland
  !
  ! Compute heat and water content to track conservation wrt dynamic land use
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_clm_type
  use mod_clm_decomp , only : get_proc_bounds
  use mod_clm_varcon , only : istsoil , istice , istwet , istdlak , &
                      istslak , isturb , istcrop , icol_road_perv , &
                      icol_road_imperv , icol_roof , icol_sunwall , &
                      icol_shadewall , cpice , cpliq , denh2o
  use mod_clm_varpar , only : nlevsno , nlevgrnd , nlevurb , nlevlak

  implicit none

  private

  save

  public :: dynland_hwcontent

  contains
  !
  ! Compute grid-level heat and water content
  !
  subroutine dynland_hwcontent(begg,endg,gcell_liq,gcell_ice,gcell_heat)
    implicit none

    ! proc beg & end gridcell indices
    integer(ik4) , intent(in)  :: begg , endg
    real(rk8) , intent(out) :: gcell_liq(begg:endg)
    real(rk8) , intent(out) :: gcell_ice(begg:endg)
    real(rk8) , intent(out) :: gcell_heat(begg:endg)

    integer(ik4)  :: li , lf         ! loop initial/final indicies
    integer(ik4)  :: ci , cf         ! loop initial/final indicies
    integer(ik4)  :: pi , pf         ! loop initial/final indicies

    ! loop indicies (grid,lunit,column,pft,vertical level)
    integer(ik4)  :: g , l , c , p , k

    real(rk8) :: wtgcell       ! weight relative to grid cell
    real(rk8) :: wtcol         ! weight relative to column
    real(rk8) :: liq           ! sum of liquid water at column level
    real(rk8) :: ice           ! sum of frozen water at column level
    real(rk8) :: heat          ! sum of heat content at column level
    real(rk8) :: cv            ! heat capacity [J/(m^2 K)]

    ! true=>do computations on this pft (see reweightMod for details)
    logical ,pointer :: pactive(:)

    integer(ik4) ,pointer :: ltype(:)   ! landunit type index
    integer(ik4) ,pointer :: ctype(:)   ! column   type index
    integer(ik4) ,pointer :: ptype(:)   ! pft      type index

    ! number of impervious road layers
    integer(ik4),  pointer :: nlev_improad(:)

    ! thermal conductivity of urban wall
    real(rk8), pointer :: cv_wall(:,:)
    ! thermal conductivity of urban roof
    real(rk8), pointer :: cv_roof(:,:)
    ! thermal conductivity of urban impervious road
    real(rk8), pointer :: cv_improad(:,:)

    integer(ik4) , pointer :: snl(:)       ! number of snow layers
    real(rk8), pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)
    real(rk8), pointer :: h2osno(:)        ! snow water (mm H2O)
    real(rk8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(rk8), pointer :: h2osoi_ice(:,:)  ! frozen water (kg/m2)
    ! volumetric soil water at saturation (porosity)
    real(rk8), pointer :: watsat(:,:)
    ! heat capacity, soil solids (J/m**3/Kelvin)
    real(rk8), pointer :: csol(:,:)
    real(rk8), pointer :: dz(:,:)  ! layer depth (m)

    type(gridcell_type) , pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type) , pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type) , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type) , pointer :: pptr  ! pointer to pft derived subtype

    ! Note: this routine does not compute heat or water content of lakes.

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    pactive => clm3%g%l%c%p%active

    ltype => clm3%g%l%itype
    ctype => clm3%g%l%c%itype
    ptype => clm3%g%l%c%p%itype

    nlev_improad => clm3%g%l%lps%nlev_improad
    cv_wall      => clm3%g%l%lps%cv_wall
    cv_roof      => clm3%g%l%lps%cv_roof
    cv_improad   => clm3%g%l%lps%cv_improad

    snl          => clm3%g%l%c%cps%snl
    watsat       => clm3%g%l%c%cps%watsat
    csol         => clm3%g%l%c%cps%csol
    dz           => clm3%g%l%c%cps%dz
    t_soisno     => clm3%g%l%c%ces%t_soisno
    h2osoi_liq   => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice   => clm3%g%l%c%cws%h2osoi_ice
    h2osno       => clm3%g%l%c%cws%h2osno

    ! Get relevant sizes

    do g = begg , endg ! loop over grid cells

      gcell_liq  (g) = 0.0_rk8   ! sum for one grid cell
      gcell_ice  (g) = 0.0_rk8   ! sum for one grid cell
      gcell_heat (g) = 0.0_rk8   ! sum for one grid cell

      li = gptr%luni(g)
      lf = gptr%lunf(g)
      do l = li , lf   ! loop over land units

        ci = lptr%coli(l)
        cf = lptr%colf(l)
        do c = ci,cf   ! loop over columns

          liq   = 0.0_rk8 ! sum for one column
          ice   = 0.0_rk8
          heat  = 0.0_rk8

          !--- water & ice, above ground only ---
          if ( (ltype(l) == istsoil .or. ltype(l) == istcrop         )  &
          .or. (ltype(l) == istwet                                   )  &
          .or. (ltype(l) == istice                                   )  &
          .or. (ltype(l) == isturb .and. ctype(c) == icol_roof       )  &
          .or. (ltype(l) == isturb .and. ctype(c) == icol_road_imperv)  &
          .or. (ltype(l) == istdlak                                  )  &
          .or. (ltype(l) == isturb .and. ctype(c) == icol_road_perv  )) then

            if ( snl(c) < 0 ) then
              do k = snl(c)+1 , 0 ! loop over snow layers
                liq   = liq   + clm3%g%l%c%cws%h2osoi_liq(c,k)
                ice   = ice   + clm3%g%l%c%cws%h2osoi_ice(c,k)
              end do
            else                 ! no snow layers exist
              ice = ice + cptr%cws%h2osno(c)
            end if
          end if

          !--- water & ice, below ground only ---
          if ( (ltype(l) == istsoil .or. ltype(l) == istcrop         )  &
          .or. (ltype(l) == istwet                                   )  &
          .or. (ltype(l) == istice                                   )  &
          .or. (ltype(l) == istdlak                                  )  &
          .or. (ltype(l) == isturb .and. ctype(c) == icol_road_perv  )) then
            do k = 1 , nlevgrnd
              liq   = liq   + cptr%cws%h2osoi_liq(c,k)
              ice   = ice   + cptr%cws%h2osoi_ice(c,k)
            end do
          end if

          !--- water & ice, below ground, for lakes ---
          if ( ltype(l) == istdlak ) then
            do k = 1 , nlevlak
              liq = liq + (1.0_rk8 - &
                      cptr%cws%lake_icefrac(c,k))*cptr%cps%dz_lake(c,k)*denh2o
              ice = ice + &
                      cptr%cws%lake_icefrac(c,k)*cptr%cps%dz_lake(c,k)*denh2o
              ! lake layers do not change thickness when freezing,
              ! so denh2o should be used (thermal properties are
              ! appropriately adjusted; see SLakeTemperatureMod)
            end do
          end if

          !--- water in aquifer ---
          if ( (ltype(l) == istsoil .or. ltype(l) == istcrop         )  &
          .or. (ltype(l) == istwet                                   )  &
          .or. (ltype(l) == istice                                   )  &
          .or. (ltype(l) == isturb .and. ctype(c) == icol_road_perv  )) then
            liq = liq + cptr%cws%wa(c)
          end if

          !--- water in canopy (at pft level) ---
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
            ! note: soil specified at LU level
            pi = cptr%pfti(c)
            pf = cptr%pftf(c)
            do p = pi , pf ! loop over pfts
              if (pactive(p)) then
                wtcol = pptr%wtcol(p)
                liq = liq + pptr%pws%h2ocan(p) * wtcol
              end if
            end do
          end if

          if ( (ltype(l) /= istslak) ) then
            ! in new lake code, all lakes are "istdlak", but
            ! have variable depth

            !--- heat content, below ground only ---
            if (nlevurb > 0) then
              do k = 1,nlevurb
                if (ctype(c)==icol_sunwall .or. ctype(c)==icol_shadewall) then
                  cv = cv_wall(l,k) * dz(c,k)
                  heat = heat + cv*t_soisno(c,k) / 1.e6_rk8
                else if (ctype(c) == icol_roof) then
                  cv = cv_roof(l,k) * dz(c,k)
                  heat = heat + cv*t_soisno(c,k) / 1.e6_rk8
                end if
              end do
            end if
            do k = 1,nlevgrnd
              if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall &
              .and. ctype(c) /= icol_roof) then
                if (ctype(c) == icol_road_imperv &
                .and. k >= 1 .and. k <= nlev_improad(l)) then
                  cv = cv_improad(l,k) * dz(c,k)
                else
                  cv = (h2osoi_ice(c,k)*cpice + h2osoi_liq(c,k)*cpliq)
                end if
                heat = heat + cv*t_soisno(c,k) / 1.e6_rk8
              end if
            end do

            !--- heat content, below ground in lake water, for lakes ---
            do k = 1 , nlevlak
              if (ltype(l) == istdlak) then
                cv = denh2o*cptr%cps%dz_lake(c,k) * &
                        ( cptr%cws%lake_icefrac(c,k)*cpice + &
                          (1.0_rk8 - cptr%cws%lake_icefrac(c,k))*cpliq )
                heat = heat + cv*cptr%ces%t_lake(c,k) / 1.e6_rk8
              end if
            end do

            !--- heat content, above ground only ---
            if ( snl(c) < 0 ) then
              do k = snl(c)+1,0 ! loop over snow layers
                cv = cpliq*h2osoi_liq(c,k) + cpice*h2osoi_ice(c,k)
                heat = heat + cv*t_soisno(c,k) / 1.e6_rk8
              end do
            else if ( h2osno(c) > 0.0_rk8 .and. ltype(l) /= istdlak) then
              ! the heat capacity (not latent heat) of snow without
              ! snow layers is currently ignored in SLakeTemperature,
              ! so it should be ignored here
              k = 1
              cv = cpice*h2osno(c)
              heat = heat + cv*t_soisno(c,k) / 1.e6_rk8
            end if

          end if

          !- scale x/m^2 column-level values into x/m^2 gridcell-level values
          wtgcell = cptr%wtgcell(c)
          gcell_liq  (g) = gcell_liq  (g) + liq   * wtgcell
          gcell_ice  (g) = gcell_ice  (g) + ice   * wtgcell
          gcell_heat (g) = gcell_heat (g) + heat  * wtgcell

        end do ! column loop
      end do ! landunit loop
    end do ! grid cell loop
  end subroutine dynland_hwcontent

end module mod_clm_dynland
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

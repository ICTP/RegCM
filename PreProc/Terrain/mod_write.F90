!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_write

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_ncstream_types
  use mod_ncstream
  use mod_memutil
  implicit none

  private

  public :: setup_outvars, write_domain, lrmoist, lrts, domain_datafile

  logical :: lrmoist, lrts

  integer(ik4) :: nvar2d
  integer(ik4) :: nvar3d
  integer(ik4) :: nvar_bnds  ! number of CF cell-boundary coordinate variables
  integer(ik4) :: idlnd ! The position of landuse in the v2dvar_base
  integer(ik4) :: idtxt ! The position of texture in the v2dvar_base
  ! Indices for centre-point coordinate variables in v2dvar_base,
  ! used to attach 'bounds' attributes after the variables are defined
  integer(ik4) :: idxlat, idxlon, iddlat, iddlon
  integer(ik4) :: idulat, idulon, idvlat, idvlon

#ifdef SINGLE_PRECISION_REAL
  type(ncvariable2d_real), save, dimension(:), allocatable :: v2dvar_base
  type(ncvariable3d_real), save, dimension(:), allocatable :: v3dvar_base
  type(ncvariable3d_real), save, dimension(:), allocatable :: v3dvar_bnds
  type(ncvariable2d_real), save :: v2dvar_lake
#else
  type(ncvariable2d_double), save, dimension(:), allocatable :: v2dvar_base
  type(ncvariable3d_double), save, dimension(:), allocatable :: v3dvar_base
  type(ncvariable3d_double), save, dimension(:), allocatable :: v3dvar_bnds
  type(ncvariable2d_double), save :: v2dvar_lake
#endif

  type domain_datafile
    character (len=256) :: fname
    logical :: lsub, lndfudge, texfudge, lakfudge
    integer(ik4) :: ntype
    real(rkx) :: ts0
    real(rkx), dimension(:), pointer, contiguous :: sigma
    real(rkx), dimension(:,:), pointer, contiguous :: xlat
    real(rkx), dimension(:,:), pointer, contiguous :: xlon
    real(rkx), dimension(:,:,:), pointer, contiguous :: xlat_bnds
    real(rkx), dimension(:,:,:), pointer, contiguous :: xlon_bnds
    real(rkx), dimension(:,:), pointer, contiguous :: dlat
    real(rkx), dimension(:,:), pointer, contiguous :: dlon
    real(rkx), dimension(:,:,:), pointer, contiguous :: dlat_bnds
    real(rkx), dimension(:,:,:), pointer, contiguous :: dlon_bnds
    real(rkx), dimension(:,:), pointer, contiguous :: ulat
    real(rkx), dimension(:,:), pointer, contiguous :: ulon
    real(rkx), dimension(:,:,:), pointer, contiguous :: ulat_bnds
    real(rkx), dimension(:,:,:), pointer, contiguous :: ulon_bnds
    real(rkx), dimension(:,:), pointer, contiguous :: vlat
    real(rkx), dimension(:,:), pointer, contiguous :: vlon
    real(rkx), dimension(:,:,:), pointer, contiguous :: vlat_bnds
    real(rkx), dimension(:,:,:), pointer, contiguous :: vlon_bnds
    real(rkx), dimension(:,:), pointer, contiguous :: xmap
    real(rkx), dimension(:,:), pointer, contiguous :: dmap
    real(rkx), dimension(:,:), pointer, contiguous :: umap
    real(rkx), dimension(:,:), pointer, contiguous :: vmap
    real(rkx), dimension(:,:), pointer, contiguous :: coriol
    real(rkx), dimension(:,:), pointer, contiguous :: mask
    real(rkx), dimension(:,:), pointer, contiguous :: htgrid
    real(rkx), dimension(:,:), pointer, contiguous :: topou
    real(rkx), dimension(:,:), pointer, contiguous :: topov
    real(rkx), dimension(:,:), pointer, contiguous :: lndout
    real(rkx), dimension(:,:), pointer, contiguous :: snowam
    real(rkx), dimension(:,:), pointer, contiguous :: smoist
    real(rkx), dimension(:,:), pointer, contiguous :: dpth
    real(rkx), dimension(:,:), pointer, contiguous :: texout
    real(rkx), dimension(:,:,:), pointer, contiguous :: rmoist
    real(rkx), dimension(:,:,:), pointer, contiguous :: rts
    real(rkx), dimension(:,:,:), pointer, contiguous :: frac_tex
    real(rkx), dimension(:,:), pointer, contiguous :: ps0
    real(rkx), dimension(:,:,:), pointer, contiguous :: pr0
    real(rkx), dimension(:,:,:), pointer, contiguous :: t0
    real(rkx), dimension(:,:,:), pointer, contiguous :: rho0
    real(rkx), dimension(:,:,:), pointer, contiguous :: z0
    real(rkx), dimension(:,:,:), pointer, contiguous :: zeta
    real(rkx), dimension(:,:,:), pointer, contiguous :: zetau
    real(rkx), dimension(:,:,:), pointer, contiguous :: zetav
    real(rkx), dimension(:,:,:), pointer, contiguous :: fmz
  end type domain_datafile

  character(len=512) :: landuse_legend =                     &
               '1  => Crop/mixed farming'//char(10)//        &
               '2  => Short grass'//char(10)//               &
               '3  => Evergreen needleleaf tree'//char(10)// &
               '4  => Deciduous needleleaf tree'//char(10)// &
               '5  => Deciduous broadleaf tree'//char(10)//  &
               '6  => Evergreen broadleaf tree'//char(10)//  &
               '7  => Tall grass'//char(10)//                &
               '8  => Desert'//char(10)//                    &
               '9  => Tundra'//char(10)//                    &
               '10 => Irrigated Crop'//char(10)//            &
               '11 => Semi-desert'//char(10)//               &
               '12 => Ice cap/glacier'//char(10)//           &
               '13 => Bog or marsh'//char(10)//              &
               '14 => Inland water'//char(10)//              &
               '15 => Ocean'//char(10)//                     &
               '16 => Evergreen shrub'//char(10)//           &
               '17 => Deciduous shrub'//char(10)//           &
               '18 => Mixed Woodland'//char(10)//            &
               '19 => Forest/Field mosaic'//char(10)//       &
               '20 => Water and Land mixture'//char(10)//    &
               '21 => Urban'//char(10)//                     &
               '22 => Sub-Urban'

  character(len=256) :: texture_legend =                     &
                '1  => Sand'//char(10)//                     &
                '2  => Loamy Sand'//char(10)//               &
                '3  => Sandy Loam'//char(10)//               &
                '4  => Silt Loam'//char(10)//                &
                '5  => Silt'//char(10)//                     &
                '6  => Loam'//char(10)//                     &
                '7  => Sandy Clay Loam'//char(10)//          &
                '8  => Silty Clay Loam'//char(10)//          &
                '9  => Clay Loam'//char(10)//                &
                '10 => Sandy Clay'//char(10)//               &
                '11 => Silty Clay'//char(10)//               &
                '12 => Clay'//char(10)//                     &
                '13 => OM'//char(10)//                       &
                '14 => Water'//char(10)//                    &
                '15 => Bedrock'//char(10)//                  &
                '16 => Other'//char(10)//                    &
                '17 => No data'

  contains

  subroutine setup_outvars
    implicit none

    lrmoist = .false.
    lrts = .false.

    ! Cross-grid bounds (xlat_bnds, xlon_bnds) and dot/B-grid bounds
    ! (dlat_bnds, dlon_bnds) are always present.  For idynamic==3 we also
    ! have u- and v-grid bounds, giving 8 bound variables total.
    if ( idynamic == 1 ) then
      nvar2d = 14
      nvar3d = 3
      nvar_bnds = 4    ! xlat_bnds, xlon_bnds, dlat_bnds, dlon_bnds
    else if ( idynamic == 2 ) then
      nvar2d = 15
      nvar3d = 7
      nvar_bnds = 4
    else if ( idynamic == 3 ) then
      nvar2d = 21
      nvar3d = 7
      nvar_bnds = 8    ! + ulat_bnds, ulon_bnds, vlat_bnds, vlon_bnds
    end if
    allocate(v2dvar_base(nvar2d))
    allocate(v3dvar_base(nvar3d))
    allocate(v3dvar_bnds(nvar_bnds))

    ! Record positions of centre-point coordinate variables so we can
    ! attach 'bounds' attributes to them after outstream_addvar.
    idxlon = 1
    idxlat = 2
    iddlon = 3
    iddlat = 4

    v2dvar_base(1)%vname = 'xlon'
    v2dvar_base(1)%vunit = 'degrees_east'
    v2dvar_base(1)%long_name = 'Longitude on Cross Points'
    v2dvar_base(1)%standard_name = 'longitude'
    v2dvar_base(2)%vname = 'xlat'
    v2dvar_base(2)%vunit = 'degrees_north'
    v2dvar_base(2)%standard_name = 'latitude'
    v2dvar_base(2)%long_name = 'Latitude on Cross Points'
    v2dvar_base(3)%vname = 'dlon'
    v2dvar_base(3)%vunit = 'degrees_east'
    v2dvar_base(3)%standard_name = 'longitude'
    v2dvar_base(3)%long_name = 'Longitude on Dot Points'
    v2dvar_base(4)%vname = 'dlat'
    v2dvar_base(4)%vunit = 'degrees_north'
    v2dvar_base(4)%standard_name = 'latitude'
    v2dvar_base(4)%long_name = 'Latitude on Dot Points'
    v2dvar_base(5)%vname = 'coriol'
    v2dvar_base(5)%vunit = 's-1'
    v2dvar_base(5)%long_name = 'Coriolis Parameter'
    v2dvar_base(5)%standard_name = 'coriolis_parameter'
    v2dvar_base(6)%vname = 'mask'
    v2dvar_base(6)%vunit = '1'
    v2dvar_base(6)%long_name = 'Land Mask'
    v2dvar_base(6)%standard_name = 'land_binary_mask'
    v2dvar_base(7)%vname = 'topo'
    v2dvar_base(7)%vunit = 'm'
    v2dvar_base(7)%long_name = 'Surface Model Elevation'
    v2dvar_base(7)%standard_name = 'surface_altitude'
    v2dvar_base(8)%vname = 'landuse'
    v2dvar_base(8)%vunit = '1'
    v2dvar_base(8)%long_name = 'Landuse category as defined in BATS1E'
    v2dvar_base(8)%standard_name = 'land_type'
    idlnd = 8
    v2dvar_base(9)%vname = 'snowam'
    v2dvar_base(9)%vunit = 'mm'
    v2dvar_base(9)%long_name = 'Snow initial LWE in mm'
    v2dvar_base(9)%standard_name = 'snowfall_amount'
    v2dvar_base(9)%lfillvalue = .true.
    v2dvar_base(10)%vname = 'smoist'
    v2dvar_base(10)%vunit = '1'
    v2dvar_base(10)%long_name = 'Soil Moisture'
    v2dvar_base(10)%standard_name = 'volume_fraction_of_water_in_soil'
    v2dvar_base(10)%lfillvalue = .true.
    v2dvar_base(11)%vname = 'texture'
    v2dvar_base(11)%vunit = '1'
    v2dvar_base(11)%long_name = 'Texture dominant category'
    v2dvar_base(11)%standard_name = 'soil_type'
    idtxt = 11
    v2dvar_base(12)%vname = 'areacella'
    v2dvar_base(12)%vunit = 'm2'
    v2dvar_base(12)%long_name = 'Atmosphere >Grid-Cell Area'
    v2dvar_base(12)%standard_name = 'cell_area'

    v2dvar_lake%vname = 'dhlake'
    v2dvar_lake%vunit = 'm'
    v2dvar_lake%long_name = 'Depth below MSL'
    v2dvar_lake%standard_name = 'depth'

    v3dvar_base(1)%vname = 'rmoist'
    v3dvar_base(1)%vunit = 'kg m-2'
    v3dvar_base(1)%long_name = 'Soil Moisture'
    v3dvar_base(1)%standard_name = 'volume_fraction_of_water_in_soil'
    v3dvar_base(1)%lfillvalue = .true.
    v3dvar_base(1)%axis = 'xys'
    v3dvar_base(2)%vname = 'rts'
    v3dvar_base(2)%vunit = 'K'
    v3dvar_base(2)%long_name = 'Soil Temperature'
    v3dvar_base(2)%standard_name = 'soil_temperature'
    v3dvar_base(2)%lfillvalue = .true.
    v3dvar_base(2)%axis = 'xys'
    v3dvar_base(3)%vname = 'texture_fraction'
    v3dvar_base(3)%vunit = '1'
    v3dvar_base(3)%long_name = 'Texture category fraction'
    v3dvar_base(3)%standard_name = 'soil_type_fraction'
    v3dvar_base(3)%axis = 'xyT'
    v3dvar_base(3)%lfillvalue = .true.

    if ( idynamic == 3 ) then
      v2dvar_base(13)%vname = 'ulon'
      v2dvar_base(13)%vunit = 'degrees_east'
      v2dvar_base(13)%long_name = 'Longitude on U Points'
      v2dvar_base(13)%standard_name = 'longitude'
      v2dvar_base(14)%vname = 'ulat'
      v2dvar_base(14)%vunit = 'degrees_north'
      v2dvar_base(14)%standard_name = 'latitude'
      v2dvar_base(14)%long_name = 'Latitude on U Points'
      v2dvar_base(15)%vname = 'vlon'
      v2dvar_base(15)%vunit = 'degrees_east'
      v2dvar_base(15)%standard_name = 'longitude'
      v2dvar_base(15)%long_name = 'Longitude on V Points'
      v2dvar_base(16)%vname = 'vlat'
      v2dvar_base(16)%vunit = 'degrees_north'
      v2dvar_base(16)%standard_name = 'latitude'
      v2dvar_base(16)%long_name = 'Latitude on V Points'
      v2dvar_base(17)%vname = 'xmap'
      v2dvar_base(17)%vunit = '1'
      v2dvar_base(17)%long_name = 'Map Factor on Cross Points'
      v2dvar_base(17)%standard_name = 'map_factor'
      v2dvar_base(18)%vname = 'umap'
      v2dvar_base(18)%vunit = '1'
      v2dvar_base(18)%long_name = 'Map Factor on U Points'
      v2dvar_base(18)%standard_name = 'map_factor'
      v2dvar_base(19)%vname = 'vmap'
      v2dvar_base(19)%vunit = '1'
      v2dvar_base(19)%long_name = 'Map Factor on V Points'
      v2dvar_base(19)%standard_name = 'map_factor'
      v2dvar_base(20)%vname = 'topou'
      v2dvar_base(20)%vunit = 'm'
      v2dvar_base(20)%long_name = 'Surface Model Elevation'
      v2dvar_base(20)%standard_name = 'surface_altitude'
      v2dvar_base(21)%vname = 'topov'
      v2dvar_base(21)%vunit = 'm'
      v2dvar_base(21)%long_name = 'Surface Model Elevation'
      v2dvar_base(21)%standard_name = 'surface_altitude'
    else
      v2dvar_base(13)%vname = 'xmap'
      v2dvar_base(13)%vunit = '1'
      v2dvar_base(13)%long_name = 'Map Factor on Cross Points'
      v2dvar_base(13)%standard_name = 'map_factor'
      v2dvar_base(14)%vname = 'dmap'
      v2dvar_base(14)%vunit = '1'
      v2dvar_base(14)%long_name = 'Map Factor on Dot Points'
      v2dvar_base(14)%standard_name = 'map_factor'
      if ( idynamic == 2 ) then
        v2dvar_base(15)%vname = 'ps0'
        v2dvar_base(15)%vunit = 'Pa'
        v2dvar_base(15)%long_name = 'Reference State Surface Pressure'
        v2dvar_base(15)%standard_name = 'air_pressure'
      end if
    end if

    if ( idynamic == 3 ) then
      v3dvar_base(4)%vname = 'zeta'
      v3dvar_base(4)%vunit = 'm'
      v3dvar_base(4)%long_name = 'Elevation above ground'
      v3dvar_base(4)%standard_name = 'heigth'
      v3dvar_base(4)%axis = 'xyz'
      v3dvar_base(5)%vname = 'fmz'
      v3dvar_base(5)%vunit = ''
      v3dvar_base(5)%long_name = 'Vertical factor'
      v3dvar_base(5)%standard_name = ''
      v3dvar_base(5)%axis = 'xyz'
      v3dvar_base(6)%vname = 'zetau'
      v3dvar_base(6)%vunit = 'm'
      v3dvar_base(6)%long_name = 'Elevation above ground'
      v3dvar_base(6)%standard_name = 'heigth'
      v3dvar_base(6)%axis = 'xyz'
      v3dvar_base(7)%vname = 'zetav'
      v3dvar_base(7)%vunit = 'm'
      v3dvar_base(7)%long_name = 'Elevation above ground'
      v3dvar_base(7)%standard_name = 'heigth'
      v3dvar_base(7)%axis = 'xyz'
    else if ( idynamic == 2 ) then
      v3dvar_base(4)%vname = 'pr0'
      v3dvar_base(4)%vunit = 'Pa'
      v3dvar_base(4)%long_name = 'Reference State Pressure'
      v3dvar_base(4)%standard_name = 'air_pressure'
      v3dvar_base(4)%axis = 'xyz'
      v3dvar_base(5)%vname = 't0'
      v3dvar_base(5)%vunit = 'K'
      v3dvar_base(5)%long_name = 'Reference State Temperature'
      v3dvar_base(5)%standard_name = 'air_temperature'
      v3dvar_base(5)%axis = 'xyz'
      v3dvar_base(6)%vname = 'rho0'
      v3dvar_base(6)%vunit = 'kg m-3'
      v3dvar_base(6)%long_name = 'Reference State Density'
      v3dvar_base(6)%standard_name = 'air_density'
      v3dvar_base(6)%axis = 'xyz'
      v3dvar_base(7)%vname = 'z0'
      v3dvar_base(7)%vunit = 'm'
      v3dvar_base(7)%long_name = 'Reference State elevation'
      v3dvar_base(7)%standard_name = 'heigth'
      v3dvar_base(7)%axis = 'xyz'
    end if

    ! ----------------------------------------------------------------
    ! CF cell-boundary coordinate variables  shape: (nv=4, jx, iy)
    ! axis code 'nxy' : n=nv dimension first (CF convention), then x, y
    ! ----------------------------------------------------------------
    v3dvar_bnds(1)%vname = 'xlon_bnds'
    v3dvar_bnds(1)%vunit = 'degrees_east'
    v3dvar_bnds(1)%long_name = 'Cell corner longitudes on Cross Points'
    v3dvar_bnds(1)%standard_name = 'longitude'
    v3dvar_bnds(1)%axis = 'nxy'
    v3dvar_bnds(2)%vname = 'xlat_bnds'
    v3dvar_bnds(2)%vunit = 'degrees_north'
    v3dvar_bnds(2)%long_name = 'Cell corner latitudes on Cross Points'
    v3dvar_bnds(2)%standard_name = 'latitude'
    v3dvar_bnds(2)%axis = 'nxy'
    v3dvar_bnds(3)%vname = 'dlon_bnds'
    v3dvar_bnds(3)%vunit = 'degrees_east'
    v3dvar_bnds(3)%long_name = 'Cell corner longitudes on Dot Points'
    v3dvar_bnds(3)%standard_name = 'longitude'
    v3dvar_bnds(3)%axis = 'nxy'
    v3dvar_bnds(4)%vname = 'dlat_bnds'
    v3dvar_bnds(4)%vunit = 'degrees_north'
    v3dvar_bnds(4)%long_name = 'Cell corner latitudes on Dot Points'
    v3dvar_bnds(4)%standard_name = 'latitude'
    v3dvar_bnds(4)%axis = 'nxy'

    if ( idynamic == 3 ) then
      ! For idynamic==3 the u/v coords live at positions 13/14 and 15/16
      idulon = 13
      idulat = 14
      idvlon = 15
      idvlat = 16
      v3dvar_bnds(5)%vname = 'ulon_bnds'
      v3dvar_bnds(5)%vunit = 'degrees_east'
      v3dvar_bnds(5)%long_name = 'Cell corner longitudes on U Points'
      v3dvar_bnds(5)%standard_name = 'longitude'
      v3dvar_bnds(5)%axis = 'nxy'
      v3dvar_bnds(6)%vname = 'ulat_bnds'
      v3dvar_bnds(6)%vunit = 'degrees_north'
      v3dvar_bnds(6)%long_name = 'Cell corner latitudes on U Points'
      v3dvar_bnds(6)%standard_name = 'latitude'
      v3dvar_bnds(6)%axis = 'nxy'
      v3dvar_bnds(7)%vname = 'vlon_bnds'
      v3dvar_bnds(7)%vunit = 'degrees_east'
      v3dvar_bnds(7)%long_name = 'Cell corner longitudes on V Points'
      v3dvar_bnds(7)%standard_name = 'longitude'
      v3dvar_bnds(7)%axis = 'nxy'
      v3dvar_bnds(8)%vname = 'vlat_bnds'
      v3dvar_bnds(8)%vunit = 'degrees_north'
      v3dvar_bnds(8)%long_name = 'Cell corner latitudes on V Points'
      v3dvar_bnds(8)%standard_name = 'latitude'
      v3dvar_bnds(8)%axis = 'nxy'
    else
      idulat = 0
      idulon = 0
      idvlat = 0
      idvlon = 0
    end if

  end subroutine setup_outvars

  subroutine write_domain(df)
    implicit none
    type(domain_datafile), intent(in) :: df
    type(nc_output_stream) :: ncout
    type(ncoutstream_params) :: opar
    real(rkx) :: dx
    integer(ik4) :: ivar

    allocate(v2dvar_base(12)%rval(size(df%xlat,1), size(df%xlat,2)))

    opar%fname = df%fname
    opar%pname = 'terrain'
    opar%l_bound = .true.
    opar%l_subgrid = df%lsub
    opar%l_full_sigma = .true.
    call outstream_setup(ncout,opar)

    call outstream_addatt(ncout, &
      ncattribute_logical('preliminary_resampling',lresamp))
    call outstream_addatt(ncout, &
          ncattribute_real8('radius_interpolation',roidem))
    call outstream_addatt(ncout, &
      ncattribute_logical('boundary_smoothing',smthbdy))
    if ( smthbdy ) then
      call outstream_addatt(ncout, &
        ncattribute_integer('smoothing_npoints',nsmthbdy))
    end if
    call outstream_addatt(ncout, &
      ncattribute_real8('minimum_h2o_pct_for_water',h2opct))
    call outstream_addatt(ncout, &
      ncattribute_integer('smoothing_level',ismthlev))
    call outstream_addatt(ncout, &
      ncattribute_logical('h2o_hgt_over_water',h2ohgt))
    call outstream_addatt(ncout, &
      ncattribute_integer('intermediate_resolution',df%ntype))
    call outstream_addatt(ncout, &
      ncattribute_logical('landuse_fudging',df%lndfudge))
    call outstream_addatt(ncout, &
      ncattribute_logical('texture_fudging',df%texfudge))
    if ( idynamic == 2 ) then
      call outstream_addatt(ncout, &
        ncattribute_real8('base_state_surface_temperature',df%ts0))
    else if ( idynamic == 3 ) then
      call outstream_addatt(ncout, &
        ncattribute_real8('top_pressure_stdatm',ptop))
    end if

    if ( lakedpth ) then
      call outstream_addatt(ncout, &
        ncattribute_logical('lake_fudging',df%lakfudge))
    end if

    call outstream_addatt(ncout, &
       ncattribute_logical('initialized_soil_moisture',lrmoist))

    call outstream_addatt(ncout, &
       ncattribute_logical('initialized_soil_temperature',lrts))

    call outstream_addatt(ncout, &
       ncattribute_logical('initialized_surface_moisture',lsmoist))

    ! ---- Define all 2D centre-point variables ----
    do ivar = 1, nvar2d
      v2dvar_base(ivar)%j1 = -1
      v2dvar_base(ivar)%j2 = -1
      v2dvar_base(ivar)%i1 = -1
      v2dvar_base(ivar)%i2 = -1
      call outstream_addvar(ncout,v2dvar_base(ivar))
    end do
    call outstream_addvaratt(ncout,v2dvar_base(idlnd), &
      ncattribute_string('legend',landuse_legend))
    call outstream_addvaratt(ncout,v2dvar_base(idtxt), &
      ncattribute_string('legend',texture_legend))

    ! ---- Define all 3D atmospheric variables ----
    do ivar = 1, nvar3d
      v3dvar_base(ivar)%j1 = -1
      v3dvar_base(ivar)%j2 = -1
      v3dvar_base(ivar)%i1 = -1
      v3dvar_base(ivar)%i2 = -1
      v3dvar_base(ivar)%k1 = -1
      v3dvar_base(ivar)%k2 = -1
      call outstream_addvar(ncout,v3dvar_base(ivar))
    end do

    ! ---- Define CF cell-boundary coordinate variables ----
    do ivar = 1, nvar_bnds
      v3dvar_bnds(ivar)%j1 = -1
      v3dvar_bnds(ivar)%j2 = -1
      v3dvar_bnds(ivar)%i1 = -1
      v3dvar_bnds(ivar)%i2 = -1
      v3dvar_bnds(ivar)%k1 = -1
      v3dvar_bnds(ivar)%k2 = -1
      call outstream_addvar(ncout,v3dvar_bnds(ivar))
    end do

    ! ---- Attach CF 'bounds' attributes to centre-point coordinate vars ----
    call outstream_addvaratt(ncout,v2dvar_base(idxlon), &
      ncattribute_string('bounds','xlon_bnds'))
    call outstream_addvaratt(ncout,v2dvar_base(idxlat), &
      ncattribute_string('bounds','xlat_bnds'))
    call outstream_addvaratt(ncout,v2dvar_base(iddlon), &
      ncattribute_string('bounds','dlon_bnds'))
    call outstream_addvaratt(ncout,v2dvar_base(iddlat), &
      ncattribute_string('bounds','dlat_bnds'))
    if ( idynamic == 3 ) then
      call outstream_addvaratt(ncout,v2dvar_base(idulon), &
        ncattribute_string('bounds','ulon_bnds'))
      call outstream_addvaratt(ncout,v2dvar_base(idulat), &
        ncattribute_string('bounds','ulat_bnds'))
      call outstream_addvaratt(ncout,v2dvar_base(idvlon), &
        ncattribute_string('bounds','vlon_bnds'))
      call outstream_addvaratt(ncout,v2dvar_base(idvlat), &
        ncattribute_string('bounds','vlat_bnds'))
    end if

    ! ---- Assign data pointers for 2D centre-point variables ----
    v2dvar_base(1)%rval => df%xlon
    v2dvar_base(2)%rval => df%xlat
    v2dvar_base(3)%rval => df%dlon
    v2dvar_base(4)%rval => df%dlat
    v2dvar_base(5)%rval => df%coriol
    v2dvar_base(6)%rval => df%mask
    v2dvar_base(7)%rval => df%htgrid
    v2dvar_base(8)%rval => df%lndout
    v2dvar_base(9)%rval => df%snowam
    v2dvar_base(10)%rval => df%smoist
    v2dvar_base(11)%rval => df%texout

    ! ---- Assign data pointers for 3D atmospheric variables ----
    v3dvar_base(1)%rval => df%rmoist
    v3dvar_base(2)%rval => df%rts
    v3dvar_base(3)%rval => df%frac_tex

    if ( lakedpth ) then
      v2dvar_lake%j1 = -1
      v2dvar_lake%j2 = -1
      v2dvar_lake%i1 = -1
      v2dvar_lake%i2 = -1
      call outstream_addvar(ncout,v2dvar_lake)
      v2dvar_lake%rval => df%dpth
    end if

    if ( idynamic == 2 ) then
      v3dvar_base(4)%rval => df%pr0
      v3dvar_base(5)%rval => df%t0
      v3dvar_base(6)%rval => df%rho0
      v3dvar_base(7)%rval => df%z0
    end if

    if ( idynamic == 3 ) then
      v2dvar_base(13)%rval => df%ulon
      v2dvar_base(14)%rval => df%ulat
      v2dvar_base(15)%rval => df%vlon
      v2dvar_base(16)%rval => df%vlat
      v2dvar_base(17)%rval => df%xmap
      v2dvar_base(18)%rval => df%umap
      v2dvar_base(19)%rval => df%vmap
      v2dvar_base(20)%rval => df%topou
      v2dvar_base(21)%rval => df%topov
      v3dvar_base(4)%rval => df%zeta
      v3dvar_base(5)%rval => df%fmz
      v3dvar_base(6)%rval => df%zetau
      v3dvar_base(7)%rval => df%zetav
    else
      v2dvar_base(13)%rval => df%xmap
      v2dvar_base(14)%rval => df%dmap
      if ( idynamic == 2 ) then
        v2dvar_base(15)%rval => df%ps0
      end if
    end if

    ! ---- Assign data pointers for cell-boundary variables ----
    v3dvar_bnds(1)%rval => df%xlon_bnds
    v3dvar_bnds(2)%rval => df%xlat_bnds
    v3dvar_bnds(3)%rval => df%dlon_bnds
    v3dvar_bnds(4)%rval => df%dlat_bnds
    if ( idynamic == 3 ) then
      v3dvar_bnds(5)%rval => df%ulon_bnds
      v3dvar_bnds(6)%rval => df%ulat_bnds
      v3dvar_bnds(7)%rval => df%vlon_bnds
      v3dvar_bnds(8)%rval => df%vlat_bnds
    end if

    dx = ds * 1000.0
    if ( iproj == 'ROTLLR' ) then
      v2dvar_base(12)%rval = (dx*dx)/df%xmap
    else
      v2dvar_base(12)%rval = (dx/df%xmap)**2
    end if

    call outstream_enable(ncout,df%sigma)

    ! ---- Write all variables ----
    do ivar = 1, nvar2d
      call outstream_writevar(ncout,v2dvar_base(ivar))
    end do

    do ivar = 1, nvar3d
      call outstream_writevar(ncout,v3dvar_base(ivar))
    end do

    do ivar = 1, nvar_bnds
      call outstream_writevar(ncout,v3dvar_bnds(ivar))
    end do

    if ( lakedpth ) then
      call outstream_writevar(ncout,v2dvar_lake)
    end if
    call outstream_dispose(ncout)
    deallocate(v2dvar_base(12)%rval)
  end subroutine write_domain

end module mod_write
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

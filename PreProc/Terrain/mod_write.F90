!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_write

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_ncstream_types
  use mod_ncstream
  use mod_memutil

  private

  public :: setup_outvars , write_domain , lrmoist

  logical :: lrmoist

  integer(ik4) :: nvar2d
  integer(ik4) :: nvar3d
  integer(ik4) :: idlnd ! The position of landuse in the v2dvar_base
  integer(ik4) :: idtxt ! The position of texture in the v2dvar_base

#ifdef SINGLE_PRECISION_REAL
  type(ncvariable2d_real) , save, dimension(:), allocatable :: v2dvar_base
  type(ncvariable3d_real) , save, dimension(:), allocatable :: v3dvar_base
  type(ncvariable2d_real) , save :: v2dvar_lake
#else
  type(ncvariable2d_double) , save, dimension(:), allocatable :: v2dvar_base
  type(ncvariable3d_double) , save, dimension(:), allocatable :: v3dvar_base
  type(ncvariable2d_double) , save :: v2dvar_lake
#endif

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

    if ( idynamic == 1 ) then
      nvar2d = 14
      nvar3d = 2
    else if ( idynamic == 2 ) then
      nvar2d = 15
      nvar3d = 6
    else if ( idynamic == 3 ) then
      nvar2d = 19
      nvar3d = 4
    end if
    allocate(v2dvar_base(nvar2d))
    allocate(v3dvar_base(nvar3d))
    v2dvar_base(1)%vname = 'xlon'
    v2dvar_base(1)%vunit = 'degrees_east'
    v2dvar_base(1)%long_name = 'Longitude on Cross Points'
    v2dvar_base(1)%standard_name = 'longitude'
    v2dvar_base(2)%vname = 'xlat'
    v2dvar_base(2)%vunit = 'degrees_north'
    v2dvar_base(2)%long_name = 'Latitude on Cross Points'
    v2dvar_base(2)%standard_name = 'latitude'
    v2dvar_base(3)%vname = 'dlon'
    v2dvar_base(3)%vunit = 'degrees_east'
    v2dvar_base(3)%long_name = 'Longitude on Dot Points'
    v2dvar_base(3)%standard_name = 'longitude'
    v2dvar_base(4)%vname = 'dlat'
    v2dvar_base(4)%vunit = 'degrees_north'
    v2dvar_base(4)%long_name = 'latitude'
    v2dvar_base(4)%standard_name = 'Latitude on Dot Points'
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
    v3dvar_base(2)%vname = 'texture_fraction'
    v3dvar_base(2)%vunit = '1'
    v3dvar_base(2)%long_name = 'Texture category fraction'
    v3dvar_base(2)%standard_name = 'soil_type_fraction'
    v3dvar_base(2)%axis = 'xyT'
    v3dvar_base(2)%lfillvalue = .true.

    if ( idynamic == 3 ) then
      v2dvar_base(13)%vname = 'ulon'
      v2dvar_base(13)%vunit = 'degrees_east'
      v2dvar_base(13)%long_name = 'Longitude on U Points'
      v2dvar_base(13)%standard_name = 'longitude'
      v2dvar_base(14)%vname = 'ulat'
      v2dvar_base(14)%vunit = 'degrees_north'
      v2dvar_base(14)%long_name = 'latitude'
      v2dvar_base(14)%standard_name = 'Latitude on U Points'
      v2dvar_base(15)%vname = 'vlon'
      v2dvar_base(15)%vunit = 'degrees_east'
      v2dvar_base(15)%long_name = 'Longitude on V Points'
      v2dvar_base(15)%standard_name = 'longitude'
      v2dvar_base(16)%vname = 'vlat'
      v2dvar_base(16)%vunit = 'degrees_north'
      v2dvar_base(16)%long_name = 'latitude'
      v2dvar_base(16)%standard_name = 'Latitude on V Points'
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
      v3dvar_base(3)%vname = 'zeta'
      v3dvar_base(3)%vunit = 'm'
      v3dvar_base(3)%long_name = 'Elevation above ground'
      v3dvar_base(3)%standard_name = 'heigth'
      v3dvar_base(3)%axis = 'xyz'
      v3dvar_base(4)%vname = 'fmz'
      v3dvar_base(4)%vunit = ''
      v3dvar_base(4)%long_name = 'Vertical factor'
      v3dvar_base(4)%standard_name = ''
      v3dvar_base(4)%axis = 'xyz'
    else if ( idynamic == 2 ) then
      v3dvar_base(3)%vname = 'pr0'
      v3dvar_base(3)%vunit = 'Pa'
      v3dvar_base(3)%long_name = 'Reference State Pressure'
      v3dvar_base(3)%standard_name = 'air_pressure'
      v3dvar_base(3)%axis = 'xyz'
      v3dvar_base(4)%vname = 't0'
      v3dvar_base(4)%vunit = 'K'
      v3dvar_base(4)%long_name = 'Reference State Temperature'
      v3dvar_base(4)%standard_name = 'air_temperature'
      v3dvar_base(4)%axis = 'xyz'
      v3dvar_base(5)%vname = 'rho0'
      v3dvar_base(5)%vunit = 'kg m-3'
      v3dvar_base(5)%long_name = 'Reference State Density'
      v3dvar_base(5)%standard_name = 'air_density'
      v3dvar_base(5)%axis = 'xyz'
      v3dvar_base(6)%vname = 'z0'
      v3dvar_base(6)%vunit = 'm'
      v3dvar_base(6)%long_name = 'Reference State elevation'
      v3dvar_base(6)%standard_name = 'heigth'
      v3dvar_base(6)%axis = 'xyz'
    end if

  end subroutine setup_outvars

  subroutine write_domain(fname,lsub,lndfudge,texfudge,lakfudge,ntype,sigma, &
                          xlat,xlon,dlat,dlon,ulat,ulon,vlat,vlon,xmap,dmap, &
                          umap,vmap,coriol,mask,htgrid,lndout,snowam,smoist, &
                          rmoist,dpth,texout,frac_tex,ps0,pr0,t0,rho0,z0,ts0,&
                          zeta,fmz)
    implicit none
    character (len=*) , intent(in) :: fname
    logical , intent(in) :: lsub , lndfudge , texfudge , lakfudge
    integer(ik4) , intent(in) :: ntype
    real(rkx) , dimension(:) , pointer , intent(in) :: sigma
    real(rkx) , dimension(:,:) , pointer , intent(in) :: xlat , xlon
    real(rkx) , dimension(:,:) , pointer , intent(in) :: dlat , dlon
    real(rkx) , dimension(:,:) , pointer , intent(in) :: ulat , ulon
    real(rkx) , dimension(:,:) , pointer , intent(in) :: vlat , vlon
    real(rkx) , dimension(:,:) , pointer , intent(in) :: xmap , dmap , coriol
    real(rkx) , dimension(:,:) , pointer , intent(in) :: umap , vmap
    real(rkx) , dimension(:,:) , pointer , intent(in) :: mask
    real(rkx) , dimension(:,:) , pointer , intent(in) :: htgrid , lndout
    real(rkx) , dimension(:,:) , pointer , intent(in) :: snowam
    real(rkx) , dimension(:,:) , pointer , intent(in) :: smoist
    real(rkx) , dimension(:,:) , pointer , intent(in) :: dpth
    real(rkx) , dimension(:,:) , pointer , intent(in) :: texout
    real(rkx) , dimension(:,:) , pointer , intent(in) :: ps0
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: rmoist
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: frac_tex
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: pr0
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: t0
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: rho0
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: z0
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: zeta
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: fmz
    real(rkx) , intent(in) :: ts0

    type(nc_output_stream) :: ncout
    type(ncoutstream_params) :: opar
    real(rkx) :: dx
    integer(ik4) :: ivar

    allocate(v2dvar_base(12)%rval(size(xlat,1), size(xlat,2)))

    opar%fname = fname
    opar%pname = 'terrain'
    opar%l_bound = .true.
    opar%l_subgrid = lsub
    opar%l_full_sigma = .true.
    call outstream_setup(ncout,opar)

    call outstream_addatt(ncout, &
      ncattribute_logical('preliminary_resampling',lresamp))
    call outstream_addatt(ncout, &
          ncattribute_real8('radius_interpolation',roidem))
    call outstream_addatt(ncout, &
      ncattribute_logical('boundary_smoothing',smthbdy))
    call outstream_addatt(ncout, &
      ncattribute_real8('minimum_h2o_pct_for_water',h2opct))
    call outstream_addatt(ncout, &
      ncattribute_integer('smoothing_level',ismthlev))
    call outstream_addatt(ncout, &
      ncattribute_logical('h2o_hgt_over_water',h2ohgt))
    call outstream_addatt(ncout, &
      ncattribute_integer('intermediate_resolution',ntype))
    call outstream_addatt(ncout, &
      ncattribute_logical('landuse_fudging',lndfudge))
    call outstream_addatt(ncout, &
      ncattribute_logical('texture_fudging',texfudge))
    if ( idynamic == 2 ) then
      call outstream_addatt(ncout, &
        ncattribute_real8('base_state_surface_temperature',ts0))
    else if ( idynamic == 3 ) then
      call outstream_addatt(ncout, &
        ncattribute_real8('top_pressure_stdatm',ptop))
    end if

    if ( lakedpth ) then
      call outstream_addatt(ncout, &
        ncattribute_logical('lake_fudging',lakfudge))
    end if

    call outstream_addatt(ncout, &
       ncattribute_logical('initialized_soil_moisture',lrmoist))

    do ivar = 1 , nvar2d
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
    v2dvar_base(1)%rval => xlon
    v2dvar_base(2)%rval => xlat
    v2dvar_base(3)%rval => dlon
    v2dvar_base(4)%rval => dlat
    v2dvar_base(5)%rval => coriol
    v2dvar_base(6)%rval => mask
    v2dvar_base(7)%rval => htgrid
    v2dvar_base(8)%rval => lndout
    v2dvar_base(9)%rval => snowam
    v2dvar_base(10)%rval => smoist
    v2dvar_base(11)%rval => texout

    do ivar = 1 , nvar3d
      v3dvar_base(ivar)%j1 = -1
      v3dvar_base(ivar)%j2 = -1
      v3dvar_base(ivar)%i1 = -1
      v3dvar_base(ivar)%i2 = -1
      v3dvar_base(ivar)%k1 = -1
      v3dvar_base(ivar)%k2 = -1
      call outstream_addvar(ncout,v3dvar_base(ivar))
    end do

    v3dvar_base(1)%rval => rmoist
    v3dvar_base(2)%rval => frac_tex

    if ( lakedpth ) then
      v2dvar_lake%j1 = -1
      v2dvar_lake%j2 = -1
      v2dvar_lake%i1 = -1
      v2dvar_lake%i2 = -1
      call outstream_addvar(ncout,v2dvar_lake)
      v2dvar_lake%rval => dpth
    end if

    if ( idynamic == 2 ) then
      v3dvar_base(3)%rval => pr0
      v3dvar_base(4)%rval => t0
      v3dvar_base(5)%rval => rho0
      v3dvar_base(6)%rval => z0
    end if

    if ( idynamic == 3 ) then
      v2dvar_base(13)%rval => ulon
      v2dvar_base(14)%rval => ulat
      v2dvar_base(15)%rval => vlon
      v2dvar_base(16)%rval => vlat
      v2dvar_base(17)%rval => xmap
      v2dvar_base(18)%rval => umap
      v2dvar_base(19)%rval => vmap
      v3dvar_base(3)%rval => zeta
      v3dvar_base(4)%rval => fmz
    else
      v2dvar_base(13)%rval => xmap
      v2dvar_base(14)%rval => dmap
      if ( idynamic == 2 ) then
        v2dvar_base(15)%rval => ps0
      end if
    end if

    dx = ds * 1000.0
    if ( iproj == 'ROTLLR' ) then
      v2dvar_base(12)%rval = (dx*dx)/xmap
    else
      v2dvar_base(12)%rval = (dx/xmap)**2
    end if

    call outstream_enable(ncout,sigma)

    do ivar = 1 , nvar2d
      call outstream_writevar(ncout,v2dvar_base(ivar))
    end do

    do ivar = 1 , nvar3d
      call outstream_writevar(ncout,v3dvar_base(ivar))
    end do

    if ( lakedpth ) then
      call outstream_writevar(ncout,v2dvar_lake)
    end if
    call outstream_dispose(ncout)
    deallocate(v2dvar_base(12)%rval)
  end subroutine write_domain

end module mod_write
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

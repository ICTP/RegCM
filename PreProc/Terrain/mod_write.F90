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

  private

  public :: setup_outvars , write_domain

  integer(ik4) , parameter :: nvar2d = 11
  type(ncvariable2d_real) , save, dimension(nvar2d) :: v2dvar_base
  integer(ik4) :: idlnd ! The position of landuse in the v2dvar_base
  type(ncvariable2d_real) , save :: v2dvar_lake
  type(ncvariable2d_real) , save :: v2dvar_texture
  type(ncvariable3d_real) , save :: v3dvar_texture

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
!
  subroutine setup_outvars
    implicit none

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
    v2dvar_base(5)%vname = 'xmap'
    v2dvar_base(5)%vunit = '1'
    v2dvar_base(5)%long_name = 'Map Factor on Cross Points'
    v2dvar_base(5)%standard_name = 'map_factor'
    v2dvar_base(6)%vname = 'dmap'
    v2dvar_base(6)%vunit = '1'
    v2dvar_base(6)%long_name = 'Map Factor on Dot Points'
    v2dvar_base(6)%standard_name = 'map_factor'
    v2dvar_base(7)%vname = 'coriol'
    v2dvar_base(7)%vunit = 's-1'
    v2dvar_base(7)%long_name = 'Coriolis Parameter'
    v2dvar_base(7)%standard_name = 'coriolis_parameter'
    v2dvar_base(8)%vname = 'mask'
    v2dvar_base(8)%vunit = '1'
    v2dvar_base(8)%long_name = 'Land Mask'
    v2dvar_base(8)%standard_name = 'land_binary_mask'
    v2dvar_base(9)%vname = 'topo'
    v2dvar_base(9)%vunit = 'm'
    v2dvar_base(9)%long_name = 'Surface Model Elevation'
    v2dvar_base(9)%standard_name = 'surface_altitude'
    v2dvar_base(10)%vname = 'landuse'
    v2dvar_base(10)%vunit = '1'
    v2dvar_base(10)%long_name = 'Landuse category as defined in BATS1E'
    v2dvar_base(10)%standard_name = 'land_type'
    idlnd = 10
    v2dvar_base(11)%vname = 'snowam'
    v2dvar_base(11)%vunit = 'mm'
    v2dvar_base(11)%long_name = 'Snow initial LWE in mm'
    v2dvar_base(11)%standard_name = 'snowfall_amount'
    v2dvar_base(11)%lfillvalue = .true.

    v2dvar_lake%vname = 'dhlake'
    v2dvar_lake%vunit = 'm'
    v2dvar_lake%long_name = 'Depth below MSL'
    v2dvar_lake%standard_name = 'depth'

    v2dvar_texture%vname = 'texture'
    v2dvar_texture%vunit = '1'
    v2dvar_texture%long_name = 'Texture dominant category'
    v2dvar_texture%standard_name = 'soil_type'
    v3dvar_texture%vname = 'texture_fraction'
    v3dvar_texture%vunit = '1'
    v3dvar_texture%long_name = 'Texture category fraction'
    v3dvar_texture%standard_name = 'soil_type_fraction'
    v3dvar_texture%axis = 'xyT'
    v3dvar_texture%lfillvalue = .true.
  end subroutine setup_outvars

  subroutine write_domain(fname,lsub,lndfudge,texfudge,lakfudge,ntype,sigma, &
                          xlat,xlon,dlat,dlon,xmap,dmap,coriol,mask,    &
                          htgrid,lndout,snowam,dpth,texout,frac_tex)
    implicit none
    character (len=*) , intent(in) :: fname
    logical , intent(in) :: lsub , lndfudge , texfudge , lakfudge
    integer(ik4) , intent(in) :: ntype
    real(rk8) , dimension(:) , pointer , intent(in) :: sigma
    real(rk8) , dimension(:,:) , pointer , intent(in) :: xlat , xlon
    real(rk8) , dimension(:,:) , pointer , intent(in) :: dlat , dlon
    real(rk8) , dimension(:,:) , pointer , intent(in) :: xmap , dmap , coriol
    real(rk8) , dimension(:,:) , pointer , intent(in) :: mask
    real(rk8) , dimension(:,:) , pointer , intent(in) :: htgrid , lndout
    real(rk8) , dimension(:,:) , pointer , intent(in) :: snowam
    real(rk8) , dimension(:,:) , pointer , intent(in) :: dpth
    real(rk8) , dimension(:,:) , pointer , intent(in) :: texout
    real(rk8) , dimension(:,:,:) , pointer , intent(in) :: frac_tex

    type(nc_output_stream) :: ncout
    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar

    opar%fname = fname
    opar%pname = 'terrain'
    opar%l_bound = .true.
    opar%l_subgrid = lsub
    opar%l_full_sigma = .true.
    call outstream_setup(ncout,opar)

    call outstream_addatt(ncout, &
      ncattribute_logical('boundary_smoothing',smthbdy))
    call outstream_addatt(ncout, &
      ncattribute_real8('minimum_h2o_pct_for_water',h2opct))
    call outstream_addatt(ncout, &
      ncattribute_integer('smoothing_level',ismthlev))
    call outstream_addatt(ncout, &
      ncattribute_logical('h2o_hgt_over_water',h2ohgt))
    call outstream_addatt(ncout, &
      ncattribute_logical('landuse_fudging',lndfudge))
    call outstream_addatt(ncout, &
      ncattribute_integer('intermediate_resolution',ntype))

    if ( lakedpth ) then
      call outstream_addatt(ncout, &
        ncattribute_logical('lake_fudging',lakfudge))
    end if
    if ( ltexture ) then
      call outstream_addatt(ncout, &
        ncattribute_logical('texture_fudging',texfudge))
    end if

    do ivar = 1 , nvar2d
      v2dvar_base(ivar)%j1 = -1
      v2dvar_base(ivar)%j2 = -1
      v2dvar_base(ivar)%i1 = -1
      v2dvar_base(ivar)%i2 = -1
      call outstream_addvar(ncout,v2dvar_base(ivar))
    end do
    call outstream_addvaratt(ncout,v2dvar_base(idlnd), &
      ncattribute_string('legend',landuse_legend))
    v2dvar_base(1)%rval => xlon
    v2dvar_base(2)%rval => xlat
    v2dvar_base(3)%rval => dlon
    v2dvar_base(4)%rval => dlat
    v2dvar_base(5)%rval => xmap
    v2dvar_base(6)%rval => dmap
    v2dvar_base(7)%rval => coriol
    v2dvar_base(8)%rval => mask
    v2dvar_base(9)%rval => htgrid
    v2dvar_base(10)%rval => lndout
    v2dvar_base(11)%rval => snowam

    if ( lakedpth ) then
      v2dvar_lake%j1 = -1
      v2dvar_lake%j2 = -1
      v2dvar_lake%i1 = -1
      v2dvar_lake%i2 = -1
      call outstream_addvar(ncout,v2dvar_lake)
      v2dvar_lake%rval => dpth
    end if

    if ( ltexture ) then
      v2dvar_texture%j1 = -1
      v2dvar_texture%j2 = -1
      v2dvar_texture%i1 = -1
      v2dvar_texture%i2 = -1
      v3dvar_texture%j1 = -1
      v3dvar_texture%j2 = -1
      v3dvar_texture%i1 = -1
      v3dvar_texture%i2 = -1
      call outstream_addvar(ncout,v2dvar_texture)
      call outstream_addvaratt(ncout,v2dvar_texture, &
        ncattribute_string('legend',texture_legend))
      call outstream_addvar(ncout,v3dvar_texture)
      v2dvar_texture%rval => texout
      v3dvar_texture%rval => frac_tex
    end if

    call outstream_enable(ncout,sigma)

    do ivar = 1 , nvar2d
      call outstream_writevar(ncout,v2dvar_base(ivar))
    end do
    if ( lakedpth ) then
      call outstream_writevar(ncout,v2dvar_lake)
    end if
    if ( ltexture ) then
      call outstream_writevar(ncout,v2dvar_texture)
      call outstream_writevar(ncout,v3dvar_texture)
    end if

    call outstream_dispose(ncout)

  end subroutine write_domain

end module mod_write

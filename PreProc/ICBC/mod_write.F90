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
  use mod_grid
  use mod_memutil
  use mod_message
  use mod_ncstream_types
  use mod_ncstream
  use mod_nhinterp
  use mod_vectutil
  use mod_stdio
  use mod_zita

  private

  real(rkx) , pointer , dimension(:) :: pcoord
  real(rkx) , pointer , dimension(:,:) :: ps4 , ts4 , wtop4 , psd0 , topod
  real(rkx) , pointer , dimension(:,:) :: pd4
  real(rkx) , pointer , dimension(:,:) :: pr , ssr , strd , clt
  real(rkx) , pointer , dimension(:,:,:) :: q4 , qc4 , qi4
  real(rkx) , pointer , dimension(:,:,:) :: t4 , u4 , v4 , z4
  real(rkx) , pointer , dimension(:,:,:) :: pp4 , ww4 , tv4 , tvd4
  real(rkx) , pointer , dimension(:,:,:) :: zud4 , zvd4

  public :: ps4 , pd4 , ts4 , q4 , t4 , u4 , v4 , pp4 , ww4 , zud4 , zvd4
  public :: qc4 , qi4 , z4 , pr , ssr , strd , clt
  public :: init_output , close_output , dispose_output , newfile , writef
  public :: init_houtput , newhfile , writehf
  public :: init_outpgw , newpgwfile , writepgwf

  type(nc_output_stream) , save :: ncout
  integer(ik4) :: nvar3d
  integer(ik4) :: nvar2d
  type(ncvariable2d_mixed) , allocatable , save , dimension(:) :: v2dvar_icbc
  type(ncvariable3d_mixed) , allocatable , save , dimension(:) :: v3dvar_icbc
  logical :: qli_present = .false.

  contains

  subroutine init_outpgw(plevs)
    implicit none
    real(rkx) , pointer , dimension(:) , intent(in) :: plevs
    integer :: nplevs , ierr

    nplevs = size(plevs)
    call getmem1d(pcoord,1,nplevs,'mod_write:pcoord')
    pcoord = plevs
    call getmem2d(ps4,1,jx,1,iy,'mod_write:ps4')
    call getmem2d(ts4,1,jx,1,iy,'mod_write:ts4')
    call getmem3d(z4,1,jx,1,iy,1,nplevs,'mod_write:z4')
    call getmem3d(q4,1,jx,1,iy,1,nplevs,'mod_write:q4')
    call getmem3d(t4,1,jx,1,iy,1,nplevs,'mod_write:t4')
    call getmem3d(u4,1,jx,1,iy,1,nplevs,'mod_write:u4')
    call getmem3d(v4,1,jx,1,iy,1,nplevs,'mod_write:v4')
    if ( idynamic == 3 ) then
      nvar2d = 10
    else
      nvar2d = 8
    end if
    nvar3d = 5
    allocate(v2dvar_icbc(nvar2d), v3dvar_icbc(nvar3d), stat=ierr)
    if ( ierr /= 0 ) then
      write(stderr,*) 'Allocation error in init_output'
      call die('icbc','Allocation error',1)
    end if
    v2dvar_icbc(1)%vname = 'xlon'
    v2dvar_icbc(1)%vunit = 'degrees_east'
    v2dvar_icbc(1)%long_name = 'Longitude on Cross Points'
    v2dvar_icbc(1)%standard_name = 'longitude'
    v2dvar_icbc(2)%vname = 'xlat'
    v2dvar_icbc(2)%vunit = 'degrees_north'
    v2dvar_icbc(2)%long_name = 'Latitude on Cross Points'
    v2dvar_icbc(2)%standard_name = 'latitude'
    v2dvar_icbc(3)%vname = 'mask'
    v2dvar_icbc(3)%vunit = '1'
    v2dvar_icbc(3)%long_name = 'Land Mask'
    v2dvar_icbc(3)%standard_name = 'land_binary_mask'
    v2dvar_icbc(4)%vname = 'topo'
    v2dvar_icbc(4)%vunit = 'm'
    v2dvar_icbc(4)%long_name = 'Surface Model Elevation'
    v2dvar_icbc(4)%standard_name = 'surface_altitude'
    v2dvar_icbc(5)%vname = 'ps'
    v2dvar_icbc(5)%vunit = 'hPa'
    v2dvar_icbc(5)%long_name = 'Surface pressure'
    v2dvar_icbc(5)%standard_name = 'surface_air_pressure'
    v2dvar_icbc(5)%lrecords = .true.
    v2dvar_icbc(6)%vname = 'ts'
    v2dvar_icbc(6)%vunit = 'K'
    v2dvar_icbc(6)%long_name = 'Surface Temperature'
    v2dvar_icbc(6)%standard_name = 'surface_temperature'
    v2dvar_icbc(6)%lrecords = .true.
    v3dvar_icbc(1)%vname = 't'
    v3dvar_icbc(1)%vunit = 'K'
    v3dvar_icbc(1)%long_name = 'Temperature'
    v3dvar_icbc(1)%standard_name = 'air_temperature'
    v3dvar_icbc(1)%lrecords = .true.
    v3dvar_icbc(2)%vname = 'qv'
    v3dvar_icbc(2)%vunit = 'kg kg-1'
    v3dvar_icbc(2)%long_name = 'Water vapor mixing ratio'
    v3dvar_icbc(2)%standard_name = 'humidity_mixing_ratio'
    v3dvar_icbc(2)%lrecords = .true.
    v3dvar_icbc(3)%vname = 'u'
    v3dvar_icbc(3)%vunit = 'm s-1'
    v3dvar_icbc(3)%long_name = 'Zonal component (westerly) of wind'
    v3dvar_icbc(3)%standard_name = 'grid_eastward_wind'
    v3dvar_icbc(3)%lrecords = .true.
    v3dvar_icbc(4)%vname = 'v'
    v3dvar_icbc(4)%vunit = 'm s-1'
    v3dvar_icbc(4)%long_name = 'Meridional component (southerly) of wind'
    v3dvar_icbc(4)%standard_name = 'grid_northward_wind'
    v3dvar_icbc(4)%lrecords = .true.
    v3dvar_icbc(5)%vname = 'z'
    v3dvar_icbc(5)%vunit = 'm'
    v3dvar_icbc(5)%long_name = 'Geopotential Height'
    v3dvar_icbc(5)%standard_name = 'geopotential_height'
    v3dvar_icbc(5)%lrecords = .true.
    if ( idynamic == 3 ) then
      v2dvar_icbc(7)%vname = 'ulon'
      v2dvar_icbc(7)%vunit = 'degrees_east'
      v2dvar_icbc(7)%long_name = 'Longitude on U Points'
      v2dvar_icbc(7)%standard_name = 'longitude'
      v2dvar_icbc(8)%vname = 'ulat'
      v2dvar_icbc(8)%vunit = 'degrees_north'
      v2dvar_icbc(8)%long_name = 'Latitude on U Points'
      v2dvar_icbc(8)%standard_name = 'latitude'
      v2dvar_icbc(9)%vname = 'vlon'
      v2dvar_icbc(9)%vunit = 'degrees_east'
      v2dvar_icbc(9)%long_name = 'Longitude on V Points'
      v2dvar_icbc(9)%standard_name = 'longitude'
      v2dvar_icbc(10)%vname = 'vlat'
      v2dvar_icbc(10)%vunit = 'degrees_north'
      v2dvar_icbc(10)%long_name = 'Latitude on V Points'
      v2dvar_icbc(10)%standard_name = 'latitude'
    else
      v2dvar_icbc(7)%vname = 'dlon'
      v2dvar_icbc(7)%vunit = 'degrees_east'
      v2dvar_icbc(7)%long_name = 'Longitude on Dot Points'
      v2dvar_icbc(7)%standard_name = 'longitude'
      v2dvar_icbc(8)%vname = 'dlat'
      v2dvar_icbc(8)%vunit = 'degrees_north'
      v2dvar_icbc(8)%long_name = 'Latitude on Dot Points'
      v2dvar_icbc(8)%standard_name = 'latitude'
    end if
  end subroutine init_outpgw

  subroutine init_houtput
    implicit none
    integer(ik4) :: ierr
    call getmem2d(pr,1,jx,1,iy,'mod_write:pr')
    call getmem2d(ssr,1,jx,1,iy,'mod_write:ssr')
    call getmem2d(strd,1,jx,1,iy,'mod_write:strd')
    call getmem2d(clt,1,jx,1,iy,'mod_write:clt')
    nvar2d = 8
    allocate(v2dvar_icbc(nvar2d), stat=ierr)
    if ( ierr /= 0 ) then
      write(stderr,*) 'Allocation error in init_houtput'
      call die('icbc','Allocation error',1)
    end if
    v2dvar_icbc(1)%vname = 'xlon'
    v2dvar_icbc(1)%vunit = 'degrees_east'
    v2dvar_icbc(1)%long_name = 'Longitude on Cross Points'
    v2dvar_icbc(1)%standard_name = 'longitude'
    v2dvar_icbc(2)%vname = 'xlat'
    v2dvar_icbc(2)%vunit = 'degrees_north'
    v2dvar_icbc(2)%long_name = 'Latitude on Cross Points'
    v2dvar_icbc(2)%standard_name = 'latitude'
    v2dvar_icbc(3)%vname = 'mask'
    v2dvar_icbc(3)%vunit = '1'
    v2dvar_icbc(3)%long_name = 'Land Mask'
    v2dvar_icbc(3)%standard_name = 'land_binary_mask'
    v2dvar_icbc(4)%vname = 'topo'
    v2dvar_icbc(4)%vunit = 'm'
    v2dvar_icbc(4)%long_name = 'Surface Model Elevation'
    v2dvar_icbc(4)%standard_name = 'surface_altitude'
    v2dvar_icbc(5)%vname = 'pr'
    v2dvar_icbc(5)%vunit = 'kg m-2 s-1'
    v2dvar_icbc(5)%long_name = 'Precipitation flux'
    v2dvar_icbc(5)%standard_name = 'precipitation_flux'
    v2dvar_icbc(5)%lrecords = .true.
    v2dvar_icbc(6)%vname = 'ssr'
    v2dvar_icbc(6)%vunit = 'W m-2'
    v2dvar_icbc(6)%long_name = 'Surface Downwelling Shortwave Flux'
    v2dvar_icbc(6)%standard_name = 'surface_downwelling_shortwave_flux_in_air'
    v2dvar_icbc(6)%lrecords = .true.
    v2dvar_icbc(7)%vname = 'strd'
    v2dvar_icbc(7)%vunit = 'W m-2'
    v2dvar_icbc(7)%long_name = 'Surface Downwelling Longwave Flux'
    v2dvar_icbc(7)%standard_name = 'surface_downwelling_longwave_flux_in_air'
    v2dvar_icbc(7)%lrecords = .true.
    v2dvar_icbc(8)%vname = 'clt'
    v2dvar_icbc(8)%vunit = '1'
    v2dvar_icbc(8)%long_name = 'Total Cloud Fraction'
    v2dvar_icbc(8)%standard_name = 'cloud_area_fraction'
    v2dvar_icbc(8)%lrecords = .true.
  end subroutine init_houtput

  subroutine init_output
    implicit none
    integer(ik4) :: ierr , i3i
    if ( dattyp == 'FNEST' .or. dattyp == 'IFSXX' ) then
      qli_present = .true.
    else
      qli_present = .false.
    end if
    call getmem2d(ps4,1,jx,1,iy,'mod_write:ps4')
    call getmem2d(ts4,1,jx,1,iy,'mod_write:ts4')
    call getmem3d(q4,1,jx,1,iy,1,kz,'mod_write:q4')
    if ( qli_present ) then
      call getmem3d(qc4,1,jx,1,iy,1,kz,'mod_write:qc4')
      call getmem3d(qi4,1,jx,1,iy,1,kz,'mod_write:qi4')
    end if
    call getmem3d(t4,1,jx,1,iy,1,kz,'mod_write:t4')
    call getmem3d(u4,1,jx,1,iy,1,kz,'mod_write:u4')
    call getmem3d(v4,1,jx,1,iy,1,kz,'mod_write:v4')
    if ( idynamic == 1 ) then
      nvar2d = 8
      nvar3d = 4
      call getmem2d(pd4,1,jx,1,iy,'mod_write:pd4')
    else if ( idynamic == 2 ) then
      nvar2d = 11
      nvar3d = 8
      call getmem2d(pd4,1,jx,1,iy,'mod_write:pd4')
      call getmem2d(wtop4,1,jx,1,iy,'mod_write:wtop4')
      call getmem2d(psd0,1,jx,1,iy,'mod_write:psd0')
      call getmem2d(topod,1,jx,1,iy,'mod_write:topod')
      call getmem3d(pp4,1,jx,1,iy,1,kz,'mod_write:pp4')
      call getmem3d(ww4,1,jx,1,iy,1,kz,'mod_write:ww4')
      call getmem3d(tv4,1,jx,1,iy,1,kz,'mod_write:tv4')
      call getmem3d(tvd4,1,jx,1,iy,1,kz,'mod_write:tvd4')
    else if ( idynamic == 3 ) then
      nvar2d = 10
      nvar3d = 4
      call getmem2d(psd0,1,jx,1,iy,'mod_write:psd0')
      call getmem3d(zud4,1,jx,1,iy,1,kz,'mod_write:zud4')
      call getmem3d(zvd4,1,jx,1,iy,1,kz,'mod_write:zvd4')
    end if
    if ( qli_present ) then
      nvar3d = nvar3d + 2
    end if
    allocate(v2dvar_icbc(nvar2d), v3dvar_icbc(nvar3d), stat=ierr)
    if ( ierr /= 0 ) then
      write(stderr,*) 'Allocation error in init_output'
      call die('icbc','Allocation error',1)
    end if
    v2dvar_icbc(1)%vname = 'xlon'
    v2dvar_icbc(1)%vunit = 'degrees_east'
    v2dvar_icbc(1)%long_name = 'Longitude on Cross Points'
    v2dvar_icbc(1)%standard_name = 'longitude'
    v2dvar_icbc(2)%vname = 'xlat'
    v2dvar_icbc(2)%vunit = 'degrees_north'
    v2dvar_icbc(2)%long_name = 'Latitude on Cross Points'
    v2dvar_icbc(2)%standard_name = 'latitude'
    v2dvar_icbc(3)%vname = 'mask'
    v2dvar_icbc(3)%vunit = '1'
    v2dvar_icbc(3)%long_name = 'Land Mask'
    v2dvar_icbc(3)%standard_name = 'land_binary_mask'
    v2dvar_icbc(4)%vname = 'topo'
    v2dvar_icbc(4)%vunit = 'm'
    v2dvar_icbc(4)%long_name = 'Surface Model Elevation'
    v2dvar_icbc(4)%standard_name = 'surface_altitude'
    v2dvar_icbc(5)%vname = 'ps'
    v2dvar_icbc(5)%vunit = 'hPa'
    v2dvar_icbc(5)%long_name = 'Surface pressure'
    v2dvar_icbc(5)%standard_name = 'surface_air_pressure'
    v2dvar_icbc(5)%lrecords = .true.
    v2dvar_icbc(6)%vname = 'ts'
    v2dvar_icbc(6)%vunit = 'K'
    v2dvar_icbc(6)%long_name = 'Surface Temperature'
    v2dvar_icbc(6)%standard_name = 'surface_temperature'
    v2dvar_icbc(6)%lrecords = .true.
    v3dvar_icbc(1)%vname = 't'
    v3dvar_icbc(1)%vunit = 'K'
    v3dvar_icbc(1)%long_name = 'Temperature'
    v3dvar_icbc(1)%standard_name = 'air_temperature'
    v3dvar_icbc(1)%lrecords = .true.
    v3dvar_icbc(2)%vname = 'qv'
    v3dvar_icbc(2)%vunit = 'kg kg-1'
    v3dvar_icbc(2)%long_name = 'Water vapor mixing ratio'
    v3dvar_icbc(2)%standard_name = 'humidity_mixing_ratio'
    v3dvar_icbc(2)%lrecords = .true.
    v3dvar_icbc(3)%vname = 'u'
    v3dvar_icbc(3)%vunit = 'm s-1'
    v3dvar_icbc(3)%long_name = 'Zonal component (westerly) of wind'
    v3dvar_icbc(3)%standard_name = 'grid_eastward_wind'
    v3dvar_icbc(3)%lrecords = .true.
    v3dvar_icbc(4)%vname = 'v'
    v3dvar_icbc(4)%vunit = 'm s-1'
    v3dvar_icbc(4)%long_name = 'Meridional component (southerly) of wind'
    v3dvar_icbc(4)%standard_name = 'grid_northward_wind'
    v3dvar_icbc(4)%lrecords = .true.
    if ( qli_present ) then
      v3dvar_icbc(5)%vname = 'qc'
      v3dvar_icbc(5)%vunit = 'kg kg-1'
      v3dvar_icbc(5)%long_name = 'Mass Fraction of Cloud Liquid Water'
      v3dvar_icbc(5)%standard_name = &
                   'mass_fraction_of_cloud_liquid_water_in_air'
      v3dvar_icbc(5)%lrecords = .true.
      v3dvar_icbc(6)%vname = 'qi'
      v3dvar_icbc(6)%vunit = 'kg kg-1'
      v3dvar_icbc(6)%long_name = 'Mass Fraction of Cloud Ice'
      v3dvar_icbc(6)%standard_name = 'mass_fraction_of_cloud_ice_in_air'
      v3dvar_icbc(6)%lrecords = .true.
      i3i = 7
    else
      i3i = 5
    end if
    if ( idynamic == 2 ) then
      v2dvar_icbc(9)%vname = 'wtop'
      v2dvar_icbc(9)%vunit = 'm s-1'
      v2dvar_icbc(9)%long_name = 'Model top vertical velocity'
      v2dvar_icbc(9)%standard_name = 'upward_air_velocity'
      v2dvar_icbc(9)%lrecords = .true.
      v2dvar_icbc(10)%vname = 'xmap'
      v2dvar_icbc(10)%vunit = '1'
      v2dvar_icbc(10)%long_name = 'Map Factor on Cross Points'
      v2dvar_icbc(10)%standard_name = 'map_factor'
      v2dvar_icbc(10)%lrecords = .false.
      v2dvar_icbc(11)%vname = 'ps0'
      v2dvar_icbc(11)%vunit = 'Pa'
      v2dvar_icbc(11)%long_name = 'Surface Reference pressure'
      v2dvar_icbc(11)%standard_name = 'air_presure'
      v2dvar_icbc(11)%lrecords = .false.
      v3dvar_icbc(i3i)%vname = 'p0'
      v3dvar_icbc(i3i)%vunit = 'Pa'
      v3dvar_icbc(i3i)%long_name = 'Reference atmospheric pressure'
      v3dvar_icbc(i3i)%standard_name = 'air_pressure'
      v3dvar_icbc(i3i)%lrecords = .false.
      v3dvar_icbc(i3i+1)%vname = 't0'
      v3dvar_icbc(i3i+1)%vunit = 'K'
      v3dvar_icbc(i3i+1)%long_name = 'Reference atmospheric temperature'
      v3dvar_icbc(i3i+1)%standard_name = 'air_temperature'
      v3dvar_icbc(i3i+1)%lrecords = .false.
      v3dvar_icbc(i3i+2)%vname = 'w'
      v3dvar_icbc(i3i+2)%vunit = 'm s-1'
      v3dvar_icbc(i3i+2)%long_name = 'Vertical wind'
      v3dvar_icbc(i3i+2)%standard_name = 'upward_air_velocity'
      v3dvar_icbc(i3i+2)%lrecords = .true.
      v3dvar_icbc(i3i+3)%vname = 'pp'
      v3dvar_icbc(i3i+3)%vunit = 'Pa'
      v3dvar_icbc(i3i+3)%long_name = 'Pressure perturbation'
      v3dvar_icbc(i3i+3)%standard_name = &
        'difference_of_air_pressure_from_model_reference'
      v3dvar_icbc(i3i+3)%lrecords = .true.
    end if
    if ( idynamic == 3 ) then
      v2dvar_icbc(7)%vname = 'ulon'
      v2dvar_icbc(7)%vunit = 'degrees_east'
      v2dvar_icbc(7)%long_name = 'Longitude on U Points'
      v2dvar_icbc(7)%standard_name = 'longitude'
      v2dvar_icbc(8)%vname = 'ulat'
      v2dvar_icbc(8)%vunit = 'degrees_north'
      v2dvar_icbc(8)%long_name = 'Latitude on U Points'
      v2dvar_icbc(8)%standard_name = 'latitude'
      v2dvar_icbc(9)%vname = 'vlon'
      v2dvar_icbc(9)%vunit = 'degrees_east'
      v2dvar_icbc(9)%long_name = 'Longitude on V Points'
      v2dvar_icbc(9)%standard_name = 'longitude'
      v2dvar_icbc(10)%vname = 'vlat'
      v2dvar_icbc(10)%vunit = 'degrees_north'
      v2dvar_icbc(10)%long_name = 'Latitude on V Points'
      v2dvar_icbc(10)%standard_name = 'latitude'
    else
      v2dvar_icbc(7)%vname = 'dlon'
      v2dvar_icbc(7)%vunit = 'degrees_east'
      v2dvar_icbc(7)%long_name = 'Longitude on Dot Points'
      v2dvar_icbc(7)%standard_name = 'longitude'
      v2dvar_icbc(8)%vname = 'dlat'
      v2dvar_icbc(8)%vunit = 'degrees_north'
      v2dvar_icbc(8)%long_name = 'Latitude on Dot Points'
      v2dvar_icbc(8)%standard_name = 'latitude'
    end if
    if ( idynamic == 2 ) then
      call crs2dot(psd0,ps0,jx,iy,i_band,i_crm)
      call crs2dot(topod,topogm,jx,iy,i_band,i_crm)
    end if
  end subroutine init_output

  subroutine close_output
    implicit none
    call outstream_dispose(ncout)
  end subroutine close_output

  subroutine dispose_output
    implicit none
    if ( allocated(v2dvar_icbc) ) deallocate(v2dvar_icbc)
    if ( allocated(v3dvar_icbc) ) deallocate(v3dvar_icbc)
  end subroutine dispose_output

  subroutine newhfile(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1

    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar
    character(len=256) :: ofname

    call outstream_dispose(ncout)
    write (ofname,'(a,a,a,a,a,a)') trim(dirglob), pthsep, &
      trim(domname), '_SFBC.', trim(tochar10(idate1)), '.nc'
    opar%fname = ofname
    opar%pname = 'clmbc'
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncout,opar)
    call outstream_addatt(ncout,ncattribute_string('global_atm_source','ERA5'))
    v2dvar_icbc(1)%rval => xlon
    v2dvar_icbc(2)%rval => xlat
    v2dvar_icbc(3)%rval => mask
    v2dvar_icbc(4)%rval => topogm
    v2dvar_icbc(5)%rval => pr
    v2dvar_icbc(6)%rval => ssr
    v2dvar_icbc(7)%rval => strd
    v2dvar_icbc(8)%rval => clt
    do ivar = 1 , nvar2d
      call outstream_addvar(ncout,v2dvar_icbc(ivar))
    end do
    call outstream_enable(ncout,sigmaf)
    do ivar = 1 , 4
      call outstream_writevar(ncout,v2dvar_icbc(ivar))
    end do
  end subroutine newhfile

  subroutine newpgwfile(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar
    character(len=256) :: ofname

    call outstream_dispose(ncout)
    write (ofname,'(a,a,a,a)') trim(dirglob), pthsep, &
      trim(domname), '_PGWBC.nc'

    opar%pname = 'PGWBC'
    opar%fname = ofname
    opar%zero_date = idate1
    opar%l_bound = .true.
    opar%l_plev = .true.
    call outstream_setup(ncout,opar)

    v2dvar_icbc(1)%rval => xlon
    v2dvar_icbc(2)%rval => xlat
    v2dvar_icbc(3)%rval => mask
    v2dvar_icbc(4)%rval => topogm
    v2dvar_icbc(5)%rval => ps4
    v2dvar_icbc(6)%rval => ts4
    v3dvar_icbc(1)%rval => t4
    v3dvar_icbc(2)%rval => q4
    v3dvar_icbc(3)%rval => u4
    v3dvar_icbc(4)%rval => v4
    v3dvar_icbc(5)%rval => z4
    if ( idynamic == 3 ) then
      v2dvar_icbc(7)%rval => ulon
      v2dvar_icbc(8)%rval => ulat
      v2dvar_icbc(9)%rval => vlon
      v2dvar_icbc(10)%rval => vlat
    else
      v2dvar_icbc(7)%rval => dlon
      v2dvar_icbc(8)%rval => dlat
    end if
    do ivar = 1 , nvar2d
      call outstream_addvar(ncout,v2dvar_icbc(ivar))
    end do
    do ivar = 1 , nvar3d
      call outstream_addvar(ncout,v3dvar_icbc(ivar))
    end do
    call outstream_enable(ncout,pcoord)
    do ivar = 1 , nvar2d
      if ( .not. v2dvar_icbc(ivar)%lrecords ) then
        call outstream_writevar(ncout,v2dvar_icbc(ivar))
      end if
    end do
  end subroutine newpgwfile

  subroutine newfile(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1

    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar , i3i
    character(len=256) :: ofname

    call outstream_dispose(ncout)
    write (ofname,'(a,a,a,a,a,a,a)') trim(dirglob), pthsep, &
      trim(domname), '_ICBC', '.', trim(tochar10(idate1)), '.nc'

    opar%pname = 'icbc'
    opar%fname = ofname
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncout,opar)

    call outstream_addatt(ncout, &
                 ncattribute_string('global_atm_source',dattyp))
    v2dvar_icbc(1)%rval => xlon
    v2dvar_icbc(2)%rval => xlat
    v2dvar_icbc(3)%rval => mask
    v2dvar_icbc(4)%rval => topogm
    v2dvar_icbc(5)%rval => ps4
    v2dvar_icbc(6)%rval => ts4
    v3dvar_icbc(1)%rval => t4
    v3dvar_icbc(2)%rval => q4
    v3dvar_icbc(3)%rval => u4
    v3dvar_icbc(4)%rval => v4
    if ( qli_present ) then
      v3dvar_icbc(5)%rval => qc4
      v3dvar_icbc(6)%rval => qi4
      i3i = 7
    else
      i3i = 5
    end if
    if ( idynamic == 2 ) then
      v2dvar_icbc(9)%rval => wtop4
      v2dvar_icbc(10)%rval => msfx
      v2dvar_icbc(11)%rval => ps0
      v3dvar_icbc(i3i)%rval => pr0
      v3dvar_icbc(i3i+1)%rval => t0
      v3dvar_icbc(i3i+2)%rval => ww4
      v3dvar_icbc(i3i+3)%rval => pp4
    end if
    if ( idynamic == 3 ) then
      v2dvar_icbc(7)%rval => ulon
      v2dvar_icbc(8)%rval => ulat
      v2dvar_icbc(9)%rval => vlon
      v2dvar_icbc(10)%rval => vlat
    else
      v2dvar_icbc(7)%rval => dlon
      v2dvar_icbc(8)%rval => dlat
    end if
    do ivar = 1 , nvar2d
      call outstream_addvar(ncout,v2dvar_icbc(ivar))
    end do
    do ivar = 1 , nvar3d
      call outstream_addvar(ncout,v3dvar_icbc(ivar))
    end do
    call outstream_enable(ncout,sigmah)
    do ivar = 1 , nvar2d
      if ( .not. v2dvar_icbc(ivar)%lrecords ) then
        call outstream_writevar(ncout,v2dvar_icbc(ivar))
      end if
    end do
    if ( idynamic == 2 ) then
      do ivar = 1 , nvar3d
        if ( .not. v3dvar_icbc(ivar)%lrecords ) then
          call outstream_writevar(ncout,v3dvar_icbc(ivar))
        end if
      end do
    end if
  end subroutine newfile

  subroutine writehf(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ivar

    call outstream_addrec(ncout,idate)
    do ivar = 5 , nvar2d
      call outstream_writevar(ncout,v2dvar_icbc(ivar))
    end do
  end subroutine writehf

  subroutine writepgwf(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ivar
    call outstream_addrec(ncout,idate)
    do ivar = 1, nvar2d
      if ( v2dvar_icbc(ivar)%lrecords ) then
        call outstream_writevar(ncout,v2dvar_icbc(ivar))
      end if
    end do
    do ivar = 1 , nvar3d
      call outstream_writevar(ncout,v3dvar_icbc(ivar))
    end do
  end subroutine writepgwf

  subroutine writef(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ivar , k
    real(rkx) :: dx

    if ( idynamic == 1 ) then
      ps4 = (ps4+ptop)*d_10
    end if

    if ( idynamic == 2 ) then
      dx = ds * d_1000
      call meandiv(u4,v4,pd4,msfd,sigmah,dsigma,jx,iy,kz,dx,jx-1,iy-1)
      tv4 = t4 * (d_one + ep1 * q4)
      do k = 1 , kz
        call crs2dot(tvd4(:,:,k),tv4(:,:,k),jx,iy,i_band,i_crm)
      end do
      ! Compute nonhydrostatic vertical velocity (w) on full sigma levels.
      call nhw(1,iy,1,jx,kz,sigmaf,dsigma,u4,v4,tv4, &
               ps4,pd4,ps0,msfx,ww4,wtop4,dx,i_band,i_crm)
      call nhinterp(1,iy,1,jx,kz,sigmah,sigmaf,u4,tvd4,pd4,psd0,1)
      call nhinterp(1,iy,1,jx,kz,sigmah,sigmaf,v4,tvd4,pd4,psd0,1)
      call nhinterp(1,iy,1,jx,kz,sigmah,sigmaf,t4,tv4,ps4,ps0,1)
      call nhinterp(1,iy,1,jx,kz,sigmah,sigmaf,q4,tv4,ps4,ps0,2)
      if ( qli_present ) then
        call nhinterp(1,iy,1,jx,kz,sigmah,sigmaf,qc4,tv4,ps4,ps0,1)
        call nhinterp(1,iy,1,jx,kz,sigmah,sigmaf,qi4,tv4,ps4,ps0,1)
      end if
      ! Recompute virtual temperature on non hydrostatic sigma.
      tv4 = t4 * (d_one + ep1 * q4)
      ! Compute the nonhydrostatic perturbation pressure field (pp).
      call nhpp(1,iy,1,jx,kz,sigmaf,t4,pr0,t0,tv4,ps4,ps0,pp4)
      ps4 = (ps4+ptop)*d_10
    end if

    if ( idynamic == 3 ) then
      ! Remember in this case ptop is zero!
      ps4 = ps4*d_10
    end if

    call outstream_addrec(ncout,idate)
    do ivar = 1, nvar2d
      if ( v2dvar_icbc(ivar)%lrecords ) then
        call outstream_writevar(ncout,v2dvar_icbc(ivar))
      end if
    end do
    do ivar = 1 , nvar3d
      call outstream_writevar(ncout,v3dvar_icbc(ivar))
    end do
  end subroutine writef

end module mod_write
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

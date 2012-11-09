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

  private

  real(rk8) , pointer , dimension(:,:) :: ps4 , ts4
  real(rk8) , pointer , dimension(:,:,:) :: h4 , q4
  real(rk8) , pointer , dimension(:,:,:) :: t4 , u4 , v4

  public :: ps4 , ts4 , h4 , q4 , t4 , u4 , v4
  public :: init_output , close_output , newfile , writef

  type(nc_output_stream) , save :: ncout
  integer(ik4) , parameter :: nvar2d = 6
  integer(ik4) , parameter :: nvar2d_static = 4
  integer(ik4) , parameter :: nvar3d = 4
  type(ncvariable2d_real) , save , dimension(nvar2d) :: v2dvar_icbc
  type(ncvariable3d_real) , save , dimension(nvar3d) :: v3dvar_icbc

  contains

  subroutine init_output
  implicit none
    call getmem2d(ps4,1,jx,1,iy,'mod_write:ps4')
    call getmem2d(ts4,1,jx,1,iy,'mod_write:ts4')
    call getmem3d(h4,1,jx,1,iy,1,kz,'mod_write:h4')
    call getmem3d(q4,1,jx,1,iy,1,kz,'mod_write:q4')
    call getmem3d(t4,1,jx,1,iy,1,kz,'mod_write:t4')
    call getmem3d(u4,1,jx,1,iy,1,kz,'mod_write:u4')
    call getmem3d(v4,1,jx,1,iy,1,kz,'mod_write:v4')
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
    v3dvar_icbc(3)%standard_name = 'eastward_wind'
    v3dvar_icbc(3)%lrecords = .true.
    v3dvar_icbc(4)%vname = 'v'
    v3dvar_icbc(4)%vunit = 'm s-1'
    v3dvar_icbc(4)%long_name = 'Meridional component (southerly) of wind'
    v3dvar_icbc(4)%standard_name = 'northward_wind'
    v3dvar_icbc(4)%lrecords = .true.
  end subroutine init_output

  subroutine close_output
    implicit none
    call outstream_dispose(ncout)
  end subroutine close_output

  subroutine newfile(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1

    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar
    character(len=256) :: ofname

    call outstream_dispose(ncout)
    write (ofname,'(a,a,a,a,i10,a)') trim(dirglob), pthsep, &
      trim(domname), '_ICBC.', toint10(idate1), '.nc'
    opar%fname = ofname
    opar%pname = 'icbc'
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncout,opar)
    call outstream_addatt(ncout,ncattribute_string('global_atm_source',dattyp))
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
    do ivar = 1 , nvar2d
      call outstream_addvar(ncout,v2dvar_icbc(ivar))
    end do
    do ivar = 1 , nvar3d
      call outstream_addvar(ncout,v3dvar_icbc(ivar))
    end do
    call outstream_enable(ncout,sigma2)
    do ivar = 1 , nvar2d_static
      call outstream_writevar(ncout,v2dvar_icbc(ivar))
    end do
  end subroutine newfile

  subroutine writef(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ivar

    ps4 = (ps4+ptop)*d_10
    call outstream_addrec(ncout,idate)
    do ivar = nvar2d_static , nvar2d
      call outstream_writevar(ncout,v2dvar_icbc(ivar))
    end do
    do ivar = 1 , nvar3d
      call outstream_writevar(ncout,v3dvar_icbc(ivar))
    end do
  end subroutine writef
!
end module mod_write

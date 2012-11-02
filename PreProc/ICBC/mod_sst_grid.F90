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

module mod_sst_grid

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_message
  use mod_nchelper
  use mod_ncstream
  use mod_domain
  use netcdf

  private

  real(rk8) , public , pointer , dimension(:,:) :: sstmm , icemm
  real(rk8) , public , pointer , dimension(:,:) :: xlat , xlon
  real(rk8) , pointer , dimension(:,:) :: topo , mask
  real(rk8) , pointer , dimension(:) :: sigma

  type(nc_output_stream) :: ncout
  integer(ik4) , parameter :: nvar2d = 4
  type(ncvariable2d_real) , dimension(nvar2d) :: v2dvar_base
  type(ncvariable2d_real) :: v2dvar_sst

  public :: init_grid , read_domain_info , setup_outvars , open_sstfile , &
            close_sstfile , writerec

  contains

  subroutine init_grid
    implicit none
    call getmem2d(sstmm,1,jx,1,iy,'mod_sst_grid:sstmm')
    call getmem2d(icemm,1,jx,1,iy,'mod_sst_grid:icemm')
    call getmem2d(xlat,1,jx,1,iy,'mod_sst_grid:xlat')
    call getmem2d(xlon,1,jx,1,iy,'mod_sst_grid:xlon')
    call getmem2d(topo,1,jx,1,iy,'mod_sst_grid:topo')
    call getmem2d(mask,1,jx,1,iy,'mod_sst_grid:mask')
    call getmem1d(sigma,1,kzp1,'mod_sst_grid:sigma')
  end subroutine init_grid

  subroutine read_domain_info(terfile)
    implicit none
    character(len=256) , intent(in) :: terfile
    integer(ik4) :: istatus , incin
    call openfile_withname(terfile,incin)
    call read_domain(incin,sigma,xlat,xlon,ht=topo,mask=mask)
    istatus = nf90_close(incin)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error closing file '//trim(terfile))
  end subroutine read_domain_info

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
    v2dvar_base(3)%vname = 'mask'
    v2dvar_base(3)%vunit = '1'
    v2dvar_base(3)%long_name = 'Land Mask'
    v2dvar_base(3)%standard_name = 'land_binary_mask'
    v2dvar_base(4)%vname = 'topo'
    v2dvar_base(4)%vunit = 'm'
    v2dvar_base(4)%long_name = 'Surface Model Elevation'
    v2dvar_base(4)%standard_name = 'surface_altitude'

    v2dvar_sst%vname = 'sst'
    v2dvar_sst%vunit = 'K'
    v2dvar_sst%long_name = 'Sea Surface Temperature'
    v2dvar_sst%standard_name = 'sea_surface_temperature'
    v2dvar_sst%lfillvalue = .true.
    v2dvar_sst%rmissval = -9999.0
    v2dvar_sst%lrecords = .true.
  end subroutine setup_outvars

  subroutine open_sstfile(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    character(len=256) :: sstname
    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar
    sstname = trim(dirglob)//pthsep//trim(domname)//'_SST.nc'
    opar%fname = sstname
    opar%pname = 'sst'
    opar%zero_date = idate1
    opar%l_bound = .true.
    opar%l_full_sigma = .true.
    call outstream_setup(ncout,opar)
    call outstream_addatt(ncout,ncattribute_string('sst_source',ssttyp))
    v2dvar_base(1)%rval => xlon
    v2dvar_base(2)%rval => xlat
    v2dvar_base(3)%rval => mask
    v2dvar_base(4)%rval => topo
    do ivar = 1 , nvar2d
      call outstream_addvar(ncout,v2dvar_base(ivar))
    end do
    call outstream_addvar(ncout,v2dvar_sst)
    v2dvar_sst%rval => sstmm
    call outstream_enable(ncout,sigma)
    do ivar = 1 , nvar2d
      call outstream_writevar(ncout,v2dvar_base(ivar))
    end do
  end subroutine open_sstfile

  subroutine close_sstfile
    implicit none
    call outstream_dispose(ncout)
  end subroutine close_sstfile

  subroutine writerec(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    call outstream_addrec(ncout,idate)
    call outstream_writevar(ncout,v2dvar_sst)
  end subroutine writerec

end module mod_sst_grid

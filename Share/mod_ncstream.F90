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

module mod_ncstream

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_constants
  use mod_memutil
  use mod_dynparam
  use mod_message
  use netcdf

  private

  integer(ik4) , parameter :: maxdims = 16
  integer(ik4) , parameter :: maxname = 32
  integer(ik4) , parameter :: maxunit = 36
  integer(ik4) , parameter :: maxstring = 256
  integer(ik4) , parameter :: maxpath = 256

  integer(ik4) , parameter :: jx_dim          = 1
  integer(ik4) , parameter :: iy_dim          = 2
  integer(ik4) , parameter :: kz_dim          = 3
  integer(ik4) , parameter :: time_dim        = 4
  integer(ik4) , parameter :: time_bound_dim  = 5
  integer(ik4) , parameter :: texture_dim     = 6
  integer(ik4) , parameter :: h2m_level_dim   = 7
  integer(ik4) , parameter :: h10m_level_dim  = 8
  integer(ik4) , parameter :: soil_layer_dim  = 9
  integer(ik4) , parameter :: water_depth_dim = 10

  integer(ik4) :: ncstat

  logical , parameter :: nocopy = .false.

  type ncstream_params
    ! The name of the output file
    character(len=maxstring) :: fname = 'regcm_generic.nc'
    ! The name of the program writing the file
    character(len=maxname) :: pname = 'ncstream test'
    ! To enable Parallel I/O. Must use hdf5 library
    integer(ik4) :: mpi_comm = -1
    integer(ik4) :: mpi_info = -1
    integer(ik4) :: mpi_iotype = -1
    ! If boundary values are part of the grids
    ! Note that in the model configuration as of now, the output is on
    ! CROSS points INTERNAL
    ! Instead for the ICBC, the output is on EXTERNAL and the CROSS grid
    ! has an extra unused line and row.
    ! JX = number of DOT zonal points WITH boundary
    ! IY = number of DOT meridional points WITH boundary
    logical :: l_bound = .false.
    ! If this is a subgrid output file.
    logical :: l_subgrid = .false.
    ! If the vertical coordinate is on full sigma levels
    logical :: l_full_sigma = .false.
    ! Initial time for this run
    real(rk8) :: zero_time = 0.0D0
    ! If parallel I/O, the processor patch indexes on the global grid
    integer(ik4) :: global_jstart , global_jend
    integer(ik4) :: global_istart , global_iend
  end type ncstream_params

  type internal_obuffer
    logical :: lhas1dint = .false.
    logical :: lhas2dint = .false.
    logical :: lhas3dint = .false.
    logical :: lhas4dint = .false.
    logical :: lhas1dreal = .false.
    logical :: lhas2dreal = .false.
    logical :: lhas3dreal = .false.
    logical :: lhas4dreal = .false.
    integer(ik4) , dimension(1) :: max1d_int = 0
    integer(ik4) , dimension(2) :: max2d_int = 0
    integer(ik4) , dimension(3) :: max3d_int = 0
    integer(ik4) , dimension(4) :: max4d_int = 0
    integer(ik4) , dimension(1) :: max1d_real = 0
    integer(ik4) , dimension(2) :: max2d_real = 0
    integer(ik4) , dimension(3) :: max3d_real = 0
    integer(ik4) , dimension(4) :: max4d_real = 0
    integer(ik4) , dimension(:) , allocatable :: intbuff
    real(rk4) , dimension(:) , allocatable :: realbuff
  end type internal_obuffer

  type obuff_p
    type(internal_obuffer) , pointer :: xb => null()
  end type obuff_p

  type ncstream
    character(len=maxpath) :: filename
    character(len=maxname) :: progname
    integer(ik4) :: id = -1
    !
    ! Defaults is output from model, i.e.:
    !   -) normal grid, not subgrid
    !   -) no boundary on cross points
    !   -) half sigma level on vertical
    !
    logical :: l_bound = .false.
    logical :: l_subgrid = .false.
    logical :: l_full_sigma = .false.
    real(rk8) :: zero_time = 0.0D0
    !
    ! Work flags
    !
    logical :: l_sync = .false.
    logical :: l_parallel = .false.
    logical :: l_hasrec = .false.
    logical :: l_hastbound = .false.
    logical :: l_hasgrid = .false.
    logical :: l_enabled = .false.
    !
    ! Dimension identifiers for 'coded' dimensions
    !
    integer(ik4) , dimension(maxdims) :: id_dims
    integer(ik4) , dimension(maxdims) :: len_dims
    integer(ik4) :: irec = 0
    ! Implemented up to 4d var with records
    integer(ik4) , dimension(5) :: istart , icount
    ! If using parallel I/O, those are the patch window
    integer(ik4) , dimension(2) :: jparbound
    integer(ik4) , dimension(2) :: iparbound
    integer(ik4) :: global_nj , global_ni , parsize
  end type ncstream
!
  type ncglobal_attribute_standard
    character(len=maxname) :: aname
  end type ncglobal_attribute_standard

  type, extends(ncglobal_attribute_standard) :: ncattribute_string
    character(len=maxstring) :: theval
  end type ncattribute_string

  type, extends(ncglobal_attribute_standard) :: ncattribute_logical
    logical :: theval
  end type ncattribute_logical

  type, extends(ncglobal_attribute_standard) :: ncattribute_integer
    integer(ik4) :: theval
  end type ncattribute_integer

  type, extends(ncglobal_attribute_standard) :: ncattribute_real4
    real(rk4) :: theval
  end type ncattribute_real4

  type, extends(ncglobal_attribute_standard) :: ncattribute_real8
    real(rk8) :: theval
  end type ncattribute_real8

  type, extends(ncglobal_attribute_standard) :: ncattribute_real4_array
    real(rk4) , dimension(maxdims) :: theval
    integer(ik4) :: numval = maxdims
  end type ncattribute_real4_array

  type, extends(ncglobal_attribute_standard) :: ncattribute_real8_array
    real(rk8) , dimension(maxdims) :: theval
    integer(ik4) :: numval = maxdims
  end type ncattribute_real8_array
!
  type ncvariable_standard
    integer(ik4) :: id = -1
    integer(ik4) :: nctype = -1
    character(len=maxname) :: vname
    character(len=maxunit) :: vunit
    character(len=maxstring) :: long_name
    character(len=maxstring) :: standard_name
    character(len=maxstring) :: cell_method = "time: point"
    integer(ik4) :: totsize = 0
    logical :: lgridded = .false.
    logical :: laddmethod = .true.
    logical :: lfillvalue = .false.
  end type ncvariable_standard

  type, extends(ncvariable_standard) :: ncvariable_0d
    logical :: lrecords = .false.
  end type ncvariable_0d

  type, extends(ncvariable_0d) :: ncvariable0d_real
    real(rk4) , dimension(1) :: rval = 0.0
  end type ncvariable0d_real

  type, extends(ncvariable_0d) :: ncvariable0d_integer
    integer(ik4) , dimension(1) :: ival = 0
  end type ncvariable0d_integer

  type, extends(ncvariable_0d) :: ncvariable0d_char
    character, dimension(1) :: cval = 'x'
  end type ncvariable0d_char

  type, extends(ncvariable_standard) :: ncvariable_1d
    logical :: lrecords = .false.
    character(len=1) :: axis = 'x'
    integer(ik4) , dimension(1) :: nval = 0
  end type ncvariable_1d

  type, extends(ncvariable_1d) :: ncvariable1d_real
    real(rk8) , dimension(:) , pointer :: rval => null()
  end type ncvariable1d_real

  type, extends(ncvariable_1d) :: ncvariable1d_integer
    integer(ik4) , dimension(:) , pointer :: ival => null()
  end type ncvariable1d_integer

  type, extends(ncvariable_standard) :: ncvariable_2d
    logical :: lrecords = .false.
    character(len=2) :: axis = 'xy'
    integer(ik4) , dimension(2) :: nval = 0
    integer(ik4) :: i1 = -1 , i2 = -1
    integer(ik4) :: j1 = -1 , j2 = -1
  end type ncvariable_2d

  type, extends(ncvariable_2d) :: ncvariable2d_real
    real(rk8) , dimension(:,:) , pointer :: rval => null()
  end type ncvariable2d_real

  type, extends(ncvariable_2d) :: ncvariable2d_integer
    integer(ik4) , dimension(:,:) , pointer :: ival => null()
  end type ncvariable2d_integer

  type, extends(ncvariable_standard) :: ncvariable_3d
    logical :: lrecords = .false.
    character(len=3) :: axis = 'xyz'
    integer(ik4) , dimension(3) :: nval = 0
    integer(ik4) :: i1 = -1 , i2 = -1
    integer(ik4) :: j1 = -1 , j2 = -1
  end type ncvariable_3d

  type, extends(ncvariable_3d) :: ncvariable3d_real
    real(rk8) , dimension(:,:,:) , pointer :: rval => null()
  end type ncvariable3d_real

  type, extends(ncvariable_3d) :: ncvariable3d_integer
    integer(rk4) , dimension(:,:,:) , pointer :: ival => null()
  end type ncvariable3d_integer

  type, extends(ncvariable_standard) :: ncvariable_4d
    logical :: lrecords = .false.
    character(len=4) :: axis = 'xyzd'
    integer(ik4) , dimension(4) :: nval = 0
    integer(ik4) :: i1 = -1 , i2 = -1
    integer(ik4) :: j1 = -1 , j2 = -1
  end type ncvariable_4d

  type, extends(ncvariable_4d) :: ncvariable4d_real
    real(rk8) , dimension(:,:,:,:) , pointer :: rval => null()
  end type ncvariable4d_real

  type, extends(ncvariable_4d) :: ncvariable4d_integer
    integer(rk4) , dimension(:,:,:,:) , pointer :: ival => null()
  end type ncvariable4d_integer

  integer(ik4) , dimension(maxdims) :: id_dim
  integer(ik4) , dimension(maxdims) :: len_dim

  type ncstream_p
    type(ncstream) , pointer :: xs => null()
  end type ncstream_p

  type basic_variables
    type(ncvariable0d_real) :: time_var
    type(ncvariable1d_real) :: tbound_var
    type(ncvariable0d_real) :: ptop_var
    type(ncvariable1d_real) :: sigma_var
    type(ncvariable1d_real) :: jx_var
    type(ncvariable1d_real) :: iy_var
    type(ncvariable0d_char) :: map_var
  end type basic_variables

  type basic_variables_p
    type(basic_variables) , pointer :: xv => null()
  end type basic_variables_p

  type nc_output_stream
    type(ncstream_p) :: ncp
    type(obuff_p) :: obp
    type(basic_variables_p) :: svp
  end type nc_output_stream

  public :: ncstream_params
  public :: nc_output_stream
  public :: ncvariable0d_char
  public :: ncvariable0d_real , ncvariable0d_integer
  public :: ncvariable1d_real , ncvariable1d_integer
  public :: ncvariable2d_real , ncvariable2d_integer
  public :: ncvariable3d_real , ncvariable3d_integer
  public :: ncvariable4d_real , ncvariable4d_integer

  public :: ncattribute_string
  public :: ncattribute_logical , ncattribute_integer
  public :: ncattribute_real4 , ncattribute_real8
  public :: ncattribute_real4_array , ncattribute_real8_array

  public :: outstream_setup
  public :: outstream_enable , outstream_dispose
  public :: outstream_addvar , outstream_addatt
  public :: outstream_addvaratt
  public :: outstream_writevar
  public :: outstream_addrec


  contains

    subroutine outstream_setup(ncout,params)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      type(ncstream_params) , intent(in) :: params
      integer(ik4) :: iomode = nf90_clobber
      type(ncstream) , pointer :: stream

      if ( associated(ncout%ncp%xs) ) call outstream_dispose(ncout)
      ! Allocate all space
      allocate(ncout%ncp%xs)
      allocate(ncout%obp%xb)
      allocate(ncout%svp%xv)
      stream => ncout%ncp%xs
      stream%filename = params%fname
#ifdef NETCDF4_HDF5
      if ( params%mpi_comm /= -1 ) then
        if ( params%mpi_iotype /= -1 ) then
          iomode = ior(ior(nf90_netcdf4,params%mpi_iotype),nf90_clobber)
        else
          iomode = ior(ior(nf90_netcdf4,nf90_mpiio),nf90_clobber)
        end if
        ncstat = nf90_create(stream%filename,iomode, &
          stream%id,comm=params%mpi_comm,info=params%mpi_info)
        stream%l_parallel = .true.
      else
        iomode = ior(ior(nf90_classic_model,nf90_clobber),nf90_netcdf4)
        ncstat = nf90_create(stream%filename,iomode,stream%id)
      end if
#else
      if ( params%mpi_comm /= -1 ) then
        !iomode = ior(nf90_pnetcdf,nf90_clobber)
        !ncstat = nf90_create(stream%filename,iomode, &
        !  params%mpi_comm,params%mpi_info,stream%id)
        !stream%l_parallel = .true.
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        call die('nc_stream','Parallel netcdf with Pnetcdf crash',1)
      else
        ncstat = nf90_create(stream%filename,iomode,stream%id)
      end if
#endif
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) nf90_strerror(ncstat)
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        call die('nc_stream','Cannot create file '//trim(stream%filename),1)
      end if
      if ( stream%l_parallel ) then
        stream%jparbound(1) = params%global_jstart
        stream%jparbound(2) = params%global_jend
        stream%iparbound(1) = params%global_istart
        stream%iparbound(2) = params%global_iend
        stream%global_ni = params%global_iend-params%global_istart+1
        stream%global_nj = params%global_jend-params%global_jstart+1
        stream%parsize = stream%global_ni*stream%global_nj
#ifdef DEBUG
        write(stdout,*) 'Parallel I/O enabled.'
        write(stdout,*) 'Processor ', myid, ' window: ', &
          stream%global_nj,' x ' , stream%global_ni, ' at ', &
          stream%jparbound(1) , ',', stream%iparbound(1)
#endif
      end if
      stream%progname     = params%pname
      stream%zero_time    = params%zero_time
      stream%l_bound      = params%l_bound
      stream%l_subgrid    = params%l_subgrid
      stream%l_full_sigma = params%l_full_sigma
      stream%id_dims(:)   = -1
      stream%len_dims(:)  = 0
      call add_common_global_params(ncout)
    end subroutine outstream_setup

    subroutine outstream_dispose(ncout)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      type(ncstream) , pointer :: stream

      if ( .not. associated(ncout%ncp%xs) ) return
      stream => ncout%ncp%xs
      if ( stream%id > 0 ) then
        ncstat = nf90_close(stream%id)
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          write(stderr,*) nf90_strerror(ncstat)
          call die('nc_stream','Cannot close file '//trim(stream%filename),1)
        end if
      end if
      call deallocate_buffers(ncout%obp%xb)
      deallocate(ncout%ncp%xs)
      deallocate(ncout%obp%xb)
      deallocate(ncout%svp%xv)
    end subroutine outstream_dispose

    subroutine deallocate_buffers(xbf)
      implicit none
      type(internal_obuffer) , pointer :: xbf
      if ( allocated(xbf%intbuff) )  deallocate(xbf%intbuff)
      if ( allocated(xbf%realbuff) ) deallocate(xbf%realbuff)
    end subroutine deallocate_buffers

    subroutine outstream_enable(ncout,sigma)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      real(rk8) , dimension(:) , pointer , intent(in) :: sigma

      type(ncstream) , pointer :: stream
      type(internal_obuffer) , pointer :: buffer
      type(basic_variables) , pointer :: stvar
      integer(ik4) :: maxnum_int , maxnum_real , i
      real(rk8) , dimension(2) :: trlat
      real(rk8) :: xds

      if ( .not. associated(ncout%ncp%xs) ) return
      stream => ncout%ncp%xs
      buffer => ncout%obp%xb
      stvar  => ncout%svp%xv
      if ( stream%l_enabled ) return

      if ( stream%l_hasrec ) then
        stvar%time_var = ncvariable0d_real(vname='time', &
            vunit='hours since 1949-12-01 00:00:00 UTC', &
            long_name='time',standard_name='time',       &
            lrecords = .true.)
        call outstream_addvar(ncout,stvar%time_var)
        call add_attribute(stream, &
          ncattribute_string('calendar',calendar), &
            stvar%time_var%id,stvar%time_var%vname)
        if ( stream%l_hastbound ) then
          call add_attribute(stream, &
            ncattribute_string('bounds','time_bnds'),stvar%time_var%id, &
            stvar%time_var%vname)
          stvar%tbound_var =  ncvariable1d_real(vname='time_bnds', &
                  vunit='hours since 1949-12-01 00:00:00 UTC',     &
                  axis='b',long_name = '', standard_name = '',     &
                  lrecords = .true.)
          call outstream_addvar(ncout,stvar%tbound_var)
          call add_attribute(stream, &
            ncattribute_string('calendar',calendar), &
              stvar%tbound_var%id,stvar%tbound_var%vname)
        end if
      end if
      if ( stream%l_hasgrid ) then
        stvar%map_var =  ncvariable0d_char(vname='rcm_map', &
          vunit='',long_name = '', standard_name = '',     &
          lrecords = .false.)
        call outstream_addvar(ncout,stvar%map_var)
        select case (iproj)
          case('LAMCON')
            call add_attribute(stream,ncattribute_string('grid_mapping_name', &
              'lambert_conformal_conic'),stvar%map_var%id,stvar%map_var%vname)
            trlat(1) = truelatl
            trlat(2) = truelath
            call add_attribute(stream, &
              ncattribute_real8_array('standard_parallel', &
              trlat,2),stvar%map_var%id,stvar%map_var%vname)
          case('POLSTR')
            call add_attribute(stream,ncattribute_string('grid_mapping_name', &
              'stereographic'),stvar%map_var%id,stvar%map_var%vname)
          case('NORMER')
            call add_attribute(stream,ncattribute_string('grid_mapping_name', &
              'mercator'),stvar%map_var%id,stvar%map_var%vname)
          case('ROTMER')
            call add_attribute(stream,ncattribute_string('grid_mapping_name', &
             'rotated_latitude_longitude'),stvar%map_var%id,stvar%map_var%vname)
            call add_attribute(stream, &
              ncattribute_real8('grid_north_pole_latitude', &
              plat),stvar%map_var%id,stvar%map_var%vname)
            call add_attribute(stream, &
              ncattribute_real8('grid_north_pole_longitude', &
              plon),stvar%map_var%id,stvar%map_var%vname)
        end select
        call add_attribute(stream, &
          ncattribute_real8('longitude_of_central_meridian', &
          clon),stvar%map_var%id,stvar%map_var%vname)
        call add_attribute(stream, &
          ncattribute_real8('latitude_of_projection_origin', &
          clat),stvar%map_var%id,stvar%map_var%vname)
      end if
      ncstat = nf90_enddef(stream%id)
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        write(stderr,*) nf90_strerror(ncstat)
        call die('nc_stream','Cannot enable file '//trim(stream%filename),1)
      end if
      !
      ! Allocate buffer space shared by all vars
      !
      maxnum_int  = buffer%max1d_int(1)
      maxnum_real = buffer%max1d_real(1)

      if ( buffer%lhas2dint ) then
        maxnum_int = max(product(buffer%max2d_int),maxnum_int)
      end if
      if ( buffer%lhas2dreal ) then
        maxnum_real = max(product(buffer%max2d_real),maxnum_real)
      end if
      if ( buffer%lhas3dint ) then
        maxnum_int = max(product(buffer%max3d_int),maxnum_int)
      end if
      if ( buffer%lhas3dreal ) then
        maxnum_real = max(product(buffer%max3d_real),maxnum_real)
      end if
      if ( buffer%lhas4dint ) then
        maxnum_int = max(product(buffer%max4d_int),maxnum_int)
      end if
      if ( buffer%lhas4dreal ) then
        maxnum_real = max(product(buffer%max4d_real),maxnum_real)
      end if
      if ( maxnum_int > 0 ) allocate(buffer%intbuff(maxnum_int))
      if ( maxnum_real > 0 ) allocate(buffer%realbuff(maxnum_real))
      stream%l_enabled = .true.
#ifdef DEBUG
      write(stdout,*) 'Enabled netCDF output stream ',trim(stream%filename)
      if ( allocated(buffer%intbuff) ) &
        write(stdout,*) 'Total buffer integer size :', &
          size(buffer%intbuff)
      if ( allocated(buffer%realbuff) ) &
        write(stdout,*) 'Total buffer float size   :', &
          size(buffer%realbuff)
#endif
      ! Put "basic" information in the file
      if ( stream%l_subgrid ) then
        xds = ds/dble(nsg)
      else
        xds = ds
      end if
      buffer%realbuff(1) = -real((dble(jx-1)/d_two) * xds * d_1000)
      do i = 2 , jx
        buffer%realbuff(i) = real(dble(buffer%realbuff(i-1))+xds)
      end do
      call outstream_writevar(ncout,stvar%jx_var,nocopy)
      buffer%realbuff(1) = -real((dble(iy-1)/d_two) * xds * d_1000)
      do i = 2 , iy
        buffer%realbuff(i) = real(dble(buffer%realbuff(i-1))+xds)
      end do
      call outstream_writevar(ncout,stvar%iy_var,nocopy)
      buffer%realbuff(1:size(sigma)) = real(sigma)
      call outstream_writevar(ncout,stvar%sigma_var,nocopy)
      stvar%ptop_var%rval(1) = real(ptop*10.0D0)
      call outstream_writevar(ncout,stvar%ptop_var)
    end subroutine outstream_enable

    subroutine outstream_sync(stream)
      implicit none
      type(ncstream) , pointer , intent(in) :: stream
      if ( .not. stream%l_enabled ) return
      if ( stream%l_sync ) then
        ncstat = nf90_sync(stream%id)
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          write(stderr,*) nf90_strerror(ncstat)
          call die('nc_stream', &
                   'Cannot sync file '//trim(stream%filename), 1)
        end if
      end if
    end subroutine outstream_sync

    subroutine outstream_addrec(ncout,val)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      real(rk4) , intent(in) :: val
      type(ncstream) , pointer :: stream
      type(basic_variables) , pointer :: stvar
      type(internal_obuffer) , pointer :: buffer
      if ( .not. associated(ncout%ncp%xs) ) return
      stream => ncout%ncp%xs
      stvar  => ncout%svp%xv
      buffer => ncout%obp%xb
      if ( .not. stream%l_enabled ) return
      stream%irec = stream%irec+1
      stvar%time_var%rval = val
      call outstream_writevar(ncout,stvar%time_var)
      if ( stream%l_hastbound ) then
        buffer%realbuff(1) = real(stream%zero_time)
        buffer%realbuff(2) = val
        call outstream_writevar(ncout,stvar%tbound_var,nocopy)
        stream%zero_time = val
      end if
    end subroutine outstream_addrec

    subroutine add_dimension(stream,dname)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      character(len=*) , intent(in) :: dname
      character(len=maxname) :: the_name
      integer(ik4) :: pdim , num
      if ( stream%l_enabled ) return
      if ( stream%id < 0 ) return
      select case (dname)
        case ('JX','jx')
          if ( stream%l_bound ) then
            ! this is the number of dot points WITH bondary
            num = jx
          else
            ! this is the number of cross points WITHOUT bondary
            num = jx - 3
          end if
          ! In subgrid , multiply for the number of points
          if ( stream%l_subgrid ) num = num*nsg
          the_name = 'jx'
          pdim = jx_dim
        case ('IY','iy')
          if ( stream%l_bound ) then
            ! this is the number of dot points WITH bondary
            num = iy
          else
            ! this is the number of cross points WITHOUT bondary
            num = iy - 3
          end if
          ! In subgrid , multiply for the number of points
          if ( stream%l_subgrid ) num = num*nsg
          the_name = 'iy'
          pdim = iy_dim
        case ('KZ','kz')
          if ( stream%l_full_sigma ) then
            ! this is the number of FULL sigma levels
            num = kz + 1
          else
            ! this is the number of HALF sigma levels
            num = kz
          end if
          the_name = 'kz'
          pdim = kz_dim
        case ('TIME','time','Time','T')
          num = nf90_unlimited
          the_name = 'time'
          pdim = time_dim
        case ('TBOUND','time_bounds','Time_bounds','TimeBounds','TB')
          num = 2
          the_name = 'time_bounds'
          pdim = time_bound_dim
        case ('NTEX','ntex','ITEX','itex','texture','textures')
          num = ntex
          the_name = 'ntex'
          pdim = texture_dim
        case ('2M','2m','lev2m','2mlev','lev_2m')
          num = 1
          the_name = 'm2'
          pdim = h2m_level_dim
        case ('10M','10m','lev10m','10mlev','lev_10m')
          num = 1
          the_name = 'm10'
          pdim = h10m_level_dim
        case ('NSOIL','nsoil','n_soil_layer','nlay','layers')
          num = 2
          the_name = 'soil_layer'
          pdim = soil_layer_dim
        case ('DPTH','dpth','NDEPTH','ndepth','depth','DEPTH')
          num = ndpmax
          the_name = 'depth'
          pdim = water_depth_dim
        case default
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', &
            'Cannot add dimension '//trim(dname)//' to file '// &
            trim(stream%filename)//': Undefined in add_dimension', 1)
      end select
      ncstat = nf90_def_dim(stream%id,the_name,num,stream%id_dims(pdim))
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        write(stderr,*) nf90_strerror(ncstat)
        call die('nc_stream', &
          'Cannot add dimension '//trim(the_name)//' to file '// &
          trim(stream%filename), 1)
      end if
      stream%len_dims(pdim) = num
    end subroutine add_dimension

    subroutine outstream_addatt(ncout,att)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      class(ncglobal_attribute_standard) , intent(in) :: att
      if ( .not. associated(ncout%ncp%xs) ) return
      call add_attribute(ncout%ncp%xs,att)
    end subroutine outstream_addatt

    subroutine outstream_addvaratt(ncout,var,att)
      implicit none
      type(nc_output_stream) , intent(in) :: ncout
      class(ncvariable_standard) , intent(in) :: var
      class(ncglobal_attribute_standard) , intent(in) :: att
      if ( .not. associated(ncout%ncp%xs) ) return
      call add_attribute(ncout%ncp%xs,att,var%id,var%vname)
    end subroutine outstream_addvaratt

    subroutine add_attribute(stream,att,iloc,vname)
      implicit none
      type(ncstream) , pointer , intent(in) :: stream
      class(ncglobal_attribute_standard) , intent(in) :: att
      integer(ik4) , optional :: iloc
      character(len=4) :: cdum
      character(len=*) , optional :: vname
      integer(ik4) :: iv
      if ( stream%l_enabled ) return
      if ( stream%id < 0 ) return
      if ( present(iloc) ) then
        iv = iloc
      else
        iv = nf90_global
      end if
      select type(att)
        class is (ncattribute_string)
          ncstat = nf90_put_att(stream%id,iv,att%aname,att%theval)
        class is (ncattribute_logical)
          call cdumlogical(cdum,att%theval)
          ncstat = nf90_put_att(stream%id,iv,att%aname,cdum)
        class is (ncattribute_integer)
          ncstat = nf90_put_att(stream%id,iv,att%aname,att%theval)
        class is (ncattribute_real4)
          ncstat = nf90_put_att(stream%id,iv,att%aname,att%theval)
        class is (ncattribute_real8)
          ncstat = nf90_put_att(stream%id,iv,att%aname,att%theval)
        class is (ncattribute_real4_array)
          ncstat = nf90_put_att(stream%id,iv, &
            att%aname,att%theval(1:att%numval))
        class is (ncattribute_real8_array)
          ncstat = nf90_put_att(stream%id,iv, &
            att%aname,att%theval(1:att%numval))
        class default
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', &
            'Cannot add attribute '//trim(att%aname)// &
            ' in file '//trim(stream%filename), 1)
      end select
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) nf90_strerror(ncstat)
        if ( present(iloc) ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', &
            'Cannot add attribute '//trim(att%aname)// &
            'to variable '//vname//' in file '//trim(stream%filename), 1)
        else
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', &
            'Cannot add global attribute '//trim(att%aname)// &
            ' in file '//trim(stream%filename), 1)
        end if
      end if
    end subroutine add_attribute

    subroutine add_varatts(stream,var)
      implicit none
      type(ncstream) , pointer , intent(in) :: stream
      class(ncvariable_standard) , intent(in) :: var
      character(len=16) :: coords_cross = 'xlat xlon'
      character(len=16) :: coords_dot   = 'dlat dlon'
      if ( len_trim(var%long_name) > 0 ) &
        call add_attribute(stream, &
          ncattribute_string('long_name',var%long_name),var%id,var%vname)
      if ( len_trim(var%standard_name) > 0 ) &
        call add_attribute(stream, &
         ncattribute_string('standard_name',var%standard_name),var%id,var%vname)
      if ( len_trim(var%vunit) > 0 ) &
        call add_attribute(stream, &
          ncattribute_string('units',var%vunit),var%id,var%vname)
      if ( var%lgridded ) then
        if ( stream%l_bound .and. &
            (var%vname == 'u' .or. var%vname == 'v') ) then
          call add_attribute(stream, &
            ncattribute_string('coordinates',coords_dot),var%id,var%vname)
        else
          call add_attribute(stream, &
            ncattribute_string('coordinates',coords_cross),var%id,var%vname)
        end if
        call add_attribute(stream, &
          ncattribute_string('grid_mapping','rcm_map'),var%id,var%vname)
        stream%l_hasgrid = .true.
      end if
      if ( var%laddmethod .and. var%vname(1:4) /= 'time' ) then
        call add_attribute(stream, &
          ncattribute_string('cell_methods',var%cell_method),var%id,var%vname)
        if ( (verify(var%cell_method(1:9),'time: mea') == 0 .or.    &
              verify(var%cell_method(1:9),'time: max') == 0 .or. &
              verify(var%cell_method(1:9),'time: min') == 0 .or. &
              verify(var%cell_method(1:9),'time: sum') == 0) .and.   &
              .not. stream%l_hastbound ) then
          stream%l_hastbound = .true.
        end if
      end if
      if ( var%lfillvalue ) then
        call add_attribute(stream, &
          ncattribute_real4('_FillValue',smissval),var%id,var%vname)
      end if
    end subroutine add_varatts

    subroutine add_variable(ncout,var,ndims)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      class(ncvariable_standard) , intent(inout) :: var
      integer(ik4) , intent(in) :: ndims
      type(ncstream) , pointer :: stream
      if ( .not. associated(ncout%ncp%xs) ) return
      stream => ncout%ncp%xs
      if ( stream%l_enabled ) return
      if ( ndims == 0 ) then
        ncstat = nf90_def_var(stream%id,var%vname,var%nctype,var%id)
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          write(stderr,*) nf90_strerror(ncstat)
          call die('nc_stream', &
            'Cannot add variable '//trim(var%vname)//' to file '// &
            trim(stream%filename), 1)
        end if
      else
        ncstat = nf90_def_var(stream%id,var%vname,var%nctype, &
                              id_dim(1:ndims),var%id)
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          write(stderr,*) nf90_strerror(ncstat)
          call die('nc_stream', &
            'Cannot define variable '//trim(var%vname)// &
            ' in file '//trim(stream%filename), 1)
        end if
#if defined(NETCDF4_HDF5) && defined (NETCDF4_COMPRESS)
        if ( ndims > 3 ) then
          ncstat = nf90_def_var_deflate(stream%id,var%id,1,1,9)
          if ( ncstat /= nf90_noerr ) then
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            write(stderr,*) nf90_strerror(ncstat)
            call die('nc_stream', &
              'Cannot set compression on variable '//trim(var%vname)// &
              ' in file '//trim(stream%filename), 1)
          end if
        end if
#endif
      end if
      call add_varatts(stream,var)
    end subroutine add_variable

    subroutine outstream_addvar(ncout,var)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      class(ncvariable_standard) :: var
      type(ncstream) , pointer :: stream
      type(internal_obuffer) , pointer :: buffer
      if ( .not. associated(ncout%ncp%xs) ) return
      stream => ncout%ncp%xs
      buffer => ncout%obp%xb
      if ( stream%l_enabled ) return
      if ( stream%id < 0 ) return
      select type(var)
        class is (ncvariable0d_char)
          var%nctype = nf90_char
        class is (ncvariable0d_real)
          var%nctype = nf90_float
        class is (ncvariable1d_real)
          var%nctype = nf90_float
        class is (ncvariable2d_real)
          var%nctype = nf90_float
        class is (ncvariable3d_real)
          var%nctype = nf90_float
        class is (ncvariable4d_real)
          var%nctype = nf90_float
        class is (ncvariable0d_integer)
          var%nctype = nf90_int
        class is (ncvariable1d_integer)
          var%nctype = nf90_int
        class is (ncvariable2d_integer)
          var%nctype = nf90_int
        class is (ncvariable3d_integer)
          var%nctype = nf90_int
        class is (ncvariable4d_integer)
          var%nctype = nf90_int
        class default
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', &
                   'Cannot add variable '//trim(var%vname)// &
                   ' in file '//trim(stream%filename), 1)
      end select
      select type(var)
        class is (ncvariable_0d)
          if ( var%lrecords ) then
            stream%l_hasrec = .true.
            call dimlist(stream,'t')
            var%lgridded = .false.
            call add_variable(ncout,var,1)
          else
            var%lgridded = .false.
            var%laddmethod = .false.
            call add_variable(ncout,var,0)
          end if
        class is (ncvariable_1d)
          if ( var%lrecords ) then
            stream%l_hasrec = .true.
            call dimlist(stream,var%axis//'t')
            var%lgridded = .false.
            call add_variable(ncout,var,2)
          else
            call dimlist(stream,var%axis)
            var%lgridded = .false.
            var%laddmethod = .false.
            call add_variable(ncout,var,1)
          end if
          var%nval(1) = len_dim(1)
          var%totsize = product(var%nval)
        class is (ncvariable_2d)
          if ( var%lrecords ) then
            stream%l_hasrec = .true.
            if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
                 var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
              var%lgridded = .true.
            end if
            call dimlist(stream,var%axis//'t')
            call add_variable(ncout,var,3)
          else
            call dimlist(stream,var%axis)
            if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
                 var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
              var%lgridded = .true.
            end if
            var%laddmethod = .false.
            call add_variable(ncout,var,2)
          end if
          var%nval(1) = len_dim(1)
          var%nval(2) = len_dim(2)
          if ( var%i1 < 1 ) then
          end if
          if ( var%j1 < 1 .and. var%j2 < 1 ) then
            var%j1 = 1
            var%j2 = len_dim(1)
          end if
          if ( var%i1 < 1 .and. var%i2 < 1 ) then
            var%i1 = 1
            var%i2 = len_dim(1)
          end if
          if ( var%j2-var%j1+1 /= var%nval(1) .or. &
               var%i2-var%i1+1 /= var%nval(2) ) then
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            write(stderr,*) 'Size of variable do not match span requested'
            call die('nc_stream', &
                   'Cannot add variable '//trim(var%vname)// &
                   ' in file '//trim(stream%filename), 1)
          end if
          var%totsize = product(var%nval)
        class is (ncvariable_3d)
          if ( var%lrecords ) then
            stream%l_hasrec = .true.
            if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
                 var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
              var%lgridded = .true.
            end if
            call dimlist(stream,var%axis//'t')
            call add_variable(ncout,var,4)
          else
            call dimlist(stream,var%axis)
            if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
                 var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
              var%lgridded = .true.
            end if
            var%laddmethod = .false.
            call add_variable(ncout,var,3)
          end if
          var%nval(1) = len_dim(1)
          var%nval(2) = len_dim(2)
          if ( var%i1 < 1 ) then
          end if
          if ( var%j1 < 1 .and. var%j2 < 1 ) then
            var%j1 = 1
            var%j2 = len_dim(1)
          end if
          if ( var%i1 < 1 .and. var%i2 < 1 ) then
            var%i1 = 1
            var%i2 = len_dim(1)
          end if
          if ( var%j2-var%j1+1 /= var%nval(1) .or. &
               var%i2-var%i1+1 /= var%nval(2) ) then
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            write(stderr,*) 'Size of variable do not match span requested'
            call die('nc_stream', &
                   'Cannot add variable '//trim(var%vname)// &
                   ' in file '//trim(stream%filename), 1)
          end if
          var%nval(3) = len_dim(3)
          var%totsize = product(var%nval)
        class is (ncvariable_4d)
          if ( var%lrecords ) then
            stream%l_hasrec = .true.
            if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
                 var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
              var%lgridded = .true.
            end if
            call dimlist(stream,var%axis//'t')
            call add_variable(ncout,var,5)
          else
            call dimlist(stream,var%axis)
            if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
                 var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
              var%lgridded = .true.
            end if
            var%laddmethod = .false.
            call add_variable(ncout,var,4)
          end if
          var%nval(1) = len_dim(1)
          var%nval(2) = len_dim(2)
          if ( var%i1 < 1 ) then
          end if
          if ( var%j1 < 1 .and. var%j2 < 1 ) then
            var%j1 = 1
            var%j2 = len_dim(1)
          end if
          if ( var%i1 < 1 .and. var%i2 < 1 ) then
            var%i1 = 1
            var%i2 = len_dim(1)
          end if
          if ( var%j2-var%j1+1 /= var%nval(1) .or. &
               var%i2-var%i1+1 /= var%nval(2) ) then
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            write(stderr,*) 'Size of variable do not match span requested'
            call die('nc_stream', &
                   'Cannot add variable '//trim(var%vname)// &
                   ' in file '//trim(stream%filename), 1)
          end if
          var%nval(3) = len_dim(3)
          var%nval(4) = len_dim(4)
          var%totsize = product(var%nval)
        class default
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', &
                   'Cannot add variable '//trim(var%vname)// &
                   ' in file '//trim(stream%filename), 1)
      end select
      select type(var)
        class is (ncvariable1d_real)
          buffer%lhas1dreal = .true.
          buffer%max1d_real(1) = max(buffer%max1d_real(1),len_dim(1))
        class is (ncvariable2d_real) 
          buffer%lhas2dreal = .true.
          if ( stream%l_parallel .and. var%lgridded ) then
            buffer%max2d_real(1) = max(buffer%max2d_real(1),stream%global_nj)
            buffer%max2d_real(2) = max(buffer%max2d_real(2),stream%global_ni)
          else
            buffer%max2d_real(1) = max(buffer%max2d_real(1),len_dim(1))
            buffer%max2d_real(2) = max(buffer%max2d_real(2),len_dim(2))
          end if
        class is (ncvariable3d_real)
          buffer%lhas3dreal = .true.
          if ( stream%l_parallel .and. var%lgridded ) then
            buffer%max3d_real(1) = max(buffer%max3d_real(1),stream%global_nj)
            buffer%max3d_real(2) = max(buffer%max3d_real(2),stream%global_ni)
          else
            buffer%max3d_real(1) = max(buffer%max3d_real(1),len_dim(1))
            buffer%max3d_real(2) = max(buffer%max3d_real(2),len_dim(2))
          end if
          buffer%max3d_real(3) = max(buffer%max3d_real(3),len_dim(3))
        class is (ncvariable4d_real)
          buffer%lhas4dreal = .true.
          if ( stream%l_parallel .and. var%lgridded ) then
            buffer%max4d_real(1) = max(buffer%max4d_real(1),stream%global_nj)
            buffer%max4d_real(2) = max(buffer%max4d_real(2),stream%global_ni)
          else
            buffer%max4d_real(1) = max(buffer%max4d_real(1),len_dim(1))
            buffer%max4d_real(2) = max(buffer%max4d_real(2),len_dim(2))
          end if
          buffer%max4d_real(3) = max(buffer%max4d_real(3),len_dim(3))
          buffer%max4d_real(4) = max(buffer%max4d_real(4),len_dim(4))
        class is (ncvariable1d_integer)
          buffer%lhas1dint = .true.
          buffer%max1d_int(1) =  max(buffer%max1d_int(1),len_dim(1))
        class is (ncvariable2d_integer)
          buffer%lhas2dint = .true.
          if ( stream%l_parallel .and. var%lgridded ) then
            buffer%max2d_int(1) = max(buffer%max2d_int(1),len_dim(1))
            buffer%max2d_int(2) = max(buffer%max2d_int(2),len_dim(2))
          else
            buffer%max2d_int(1) = max(buffer%max2d_int(1),stream%global_nj)
            buffer%max2d_int(2) = max(buffer%max2d_int(2),stream%global_ni)
          end if
        class is (ncvariable3d_integer)
          buffer%lhas3dint = .true.
          if ( stream%l_parallel .and. var%lgridded ) then
            buffer%max3d_int(1) = max(buffer%max3d_int(1),len_dim(1))
            buffer%max3d_int(2) = max(buffer%max3d_int(2),len_dim(2))
          else
            buffer%max3d_int(1) = max(buffer%max3d_int(1),stream%global_nj)
            buffer%max3d_int(2) = max(buffer%max3d_int(2),stream%global_ni)
          end if
          buffer%max3d_int(3) = max(buffer%max3d_int(3),len_dim(3))
        class is (ncvariable4d_integer)
          buffer%lhas4dint = .true.
          if ( stream%l_parallel .and. var%lgridded ) then
            buffer%max4d_int(1) = max(buffer%max4d_int(1),len_dim(1))
            buffer%max4d_int(2) = max(buffer%max4d_int(2),len_dim(2))
          else
            buffer%max4d_int(1) = max(buffer%max4d_int(1),stream%global_nj)
            buffer%max4d_int(2) = max(buffer%max4d_int(2),stream%global_ni)
          end if
          buffer%max4d_int(3) = max(buffer%max4d_int(3),len_dim(3))
          buffer%max4d_int(4) = max(buffer%max4d_int(4),len_dim(4))
        class default
          continue
      end select
    end subroutine outstream_addvar

    subroutine dimlist(stream,code)
      implicit none
      type(ncstream) , intent(inout) , pointer :: stream
      character(len=*) , intent(in) :: code
      character(len=maxdims) :: safecode
      integer(ik4) :: ic

      do ic = 1 , maxdims
        safecode(ic:ic) = ' '
      end do
      safecode = code
      ic = 1
      do while ( safecode(ic:ic) /= ' ' )
        select case (safecode(ic:ic))
          case ('x')
            if ( stream%id_dims(jx_dim) < 0 ) then
              call add_dimension(stream,'jx')
            end if
            id_dim(ic) = stream%id_dims(jx_dim)
            len_dim(ic) = stream%len_dims(jx_dim)
          case ('y')
            if ( stream%id_dims(iy_dim) < 0 ) then
              call add_dimension(stream,'iy')
            end if
            id_dim(ic) = stream%id_dims(iy_dim)
            len_dim(ic) = stream%len_dims(iy_dim)
          case ('z')
            if ( stream%id_dims(kz_dim) < 0 ) then
              call add_dimension(stream,'kz')
            end if
            id_dim(ic) = stream%id_dims(kz_dim)
            len_dim(ic) = stream%len_dims(kz_dim)
          case ('t')
            if ( stream%id_dims(time_dim) < 0 ) then
              call add_dimension(stream,'time')
            end if
            id_dim(ic) = stream%id_dims(time_dim)
            len_dim(ic) = stream%len_dims(time_dim)
          case ('b')
            if ( stream%id_dims(time_bound_dim) < 0 ) then
              call add_dimension(stream,'time_bounds')
            end if
            id_dim(ic) = stream%id_dims(time_bound_dim)
            len_dim(ic) = stream%len_dims(time_bound_dim)
          case ('T')
            if ( stream%id_dims(texture_dim) < 0 ) then
              call add_dimension(stream,'ntex')
            end if
            id_dim(ic) = stream%id_dims(texture_dim)
            len_dim(ic) = stream%len_dims(texture_dim)
          case ('2')
            if ( stream%id_dims(h2m_level_dim) < 0 ) then
              call add_dimension(stream,'2m')
            end if
            id_dim(ic) = stream%id_dims(h2m_level_dim)
            len_dim(ic) = stream%len_dims(h2m_level_dim)
          case ('w')
            if ( stream%id_dims(h10m_level_dim) < 0 ) then
              call add_dimension(stream,'h10')
            end if
            id_dim(ic) = stream%id_dims(h10m_level_dim)
            len_dim(ic) = stream%len_dims(h10m_level_dim)
          case ('s')
            if ( stream%id_dims(soil_layer_dim) < 0 ) then
              call add_dimension(stream,'soil_layer')
            end if
            id_dim(ic) = stream%id_dims(soil_layer_dim)
            len_dim(ic) = stream%len_dims(soil_layer_dim)
          case ('d')
            if ( stream%id_dims(water_depth_dim) < 0 ) then
              call add_dimension(stream,'depth')
            end if
            id_dim(ic) = stream%id_dims(water_depth_dim)
            len_dim(ic) = stream%len_dims(water_depth_dim)
          case default
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            write(stderr,*) 'Not in list. Known dimension codes: xyztbT2wsd'
            call die('nc_stream', &
              'Cannot select dimension '//safecode(ic:ic)//' on file '// &
              trim(stream%filename), 1)
        end select
        ic = ic + 1
      end do
    end subroutine dimlist

    subroutine outstream_writevar(ncout,var,lcopy)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      class(ncvariable_standard) , intent(inout) :: var
      logical , intent(in) , optional :: lcopy
      type(ncstream) , pointer :: stream
      type(internal_obuffer) , pointer :: buffer
      logical :: docopy
      docopy = .true.
      if ( present(lcopy) ) docopy = lcopy
      stream => ncout%ncp%xs
      if ( .not. stream%l_enabled ) return
      if ( stream%id < 0 ) return
      buffer => ncout%obp%xb
      select type(var)
        class is (ncvariable0d_real)
          if ( var%lrecords ) then
            stream%istart(1) = stream%irec
            stream%icount(1) = 1
            ncstat = nf90_put_var(stream%id,var%id,var%rval, &
              stream%istart(1:1),stream%icount(1:1))
          else
            ncstat = nf90_put_var(stream%id,var%id,var%rval(1))
          end if
        class is (ncvariable1d_real)
          if ( docopy ) then
            if ( .not. associated(var%rval) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%realbuff(1:var%nval(1)) = real(var%rval(1:var%nval(1)))
          end if
          if ( var%lrecords ) then
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = stream%irec
            stream%icount(2) = 1
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%realbuff,stream%istart(1:2),stream%icount(1:2))
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%realbuff,stream%istart(1:1),stream%icount(1:1))
          end if
        class is (ncvariable2d_real)
          if ( docopy ) then
            if ( .not. associated(var%rval) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%realbuff(1:var%totsize) = &
              real(reshape(var%rval(var%j1:var%j2,var%i1:var%i2), &
              (/var%totsize/)))
          end if
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
          end if
          if ( var%lrecords ) then
            stream%istart(3) = stream%irec
            stream%icount(3) = 1
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%realbuff,stream%istart(1:3),stream%icount(1:3))
          else
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%realbuff,stream%istart(1:2),stream%icount(1:2))
          end if
        class is (ncvariable3d_real)
          if ( docopy ) then
            if ( .not. associated(var%rval) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%realbuff(1:var%totsize) = &
              real(reshape(var%rval(var%j1:var%j2,var%i1:var%i2,:), &
              (/var%totsize/)))
          end if
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
          end if
          stream%istart(3) = 1
          stream%icount(3) = var%nval(3)
          if ( var%lrecords ) then
            stream%istart(4) = stream%irec
            stream%icount(4) = 1
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%realbuff,stream%istart(1:4),stream%icount(1:4))
          else
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%realbuff,stream%istart(1:3),stream%icount(1:3))
          end if
        class is (ncvariable4d_real)
          if ( docopy ) then
            if ( .not. associated(var%rval) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%realbuff(1:var%totsize) = &
              real(reshape(var%rval(var%j1:var%j2,var%i1:var%i2,:,:), &
              (/var%totsize/)))
          end if
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
          end if
          stream%istart(3) = 1
          stream%icount(3) = var%nval(3)
          stream%istart(4) = 1
          stream%icount(4) = var%nval(4)
          if ( var%lrecords ) then
            stream%istart(5) = stream%irec
            stream%icount(5) = 1
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%realbuff,stream%istart(1:5),stream%icount(1:5))
          else
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%realbuff,stream%istart(1:4),stream%icount(1:4))
          end if
        class is (ncvariable0d_integer)
          if ( var%lrecords ) then
            stream%istart(1) = stream%irec
            stream%icount(1) = 1
            ncstat = nf90_put_var(stream%id,var%id,var%ival, &
              stream%istart(1:1),stream%icount(1:1))
          else
            ncstat = nf90_put_var(stream%id,var%id,var%ival(1))
          end if
        class is (ncvariable1d_integer)
          if ( docopy ) then
            if ( .not. associated(var%ival) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%intbuff(1:var%nval(1)) = var%ival(1:var%nval(1))
          end if
          if ( var%lrecords ) then
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = stream%irec
            stream%icount(2) = 1
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%intbuff,stream%istart(1:2),stream%icount(1:2))
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%intbuff,stream%istart(1:1),stream%icount(1:1))
          end if
        class is (ncvariable2d_integer)
          if ( docopy ) then
            if ( .not. associated(var%ival) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%intbuff(1:var%totsize) = &
              reshape(var%ival(var%j1:var%j2,var%i1:var%i2),(/var%totsize/))
          end if
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
          end if
          if ( var%lrecords ) then
            stream%istart(3) = stream%irec
            stream%icount(3) = 1
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%intbuff,stream%istart(1:3),stream%icount(1:3))
          else
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%intbuff,stream%istart(1:2),stream%icount(1:2))
          end if
        class is (ncvariable3d_integer)
          if ( docopy ) then
            if ( .not. associated(var%ival) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%intbuff(1:var%totsize) = &
              reshape(var%ival(var%j1:var%j2,var%i1:var%i2,:),(/var%totsize/))
          end if
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
          end if
          stream%istart(3) = 1
          stream%icount(3) = var%nval(3)
          if ( var%lrecords ) then
            stream%istart(4) = stream%irec
            stream%icount(4) = 1
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%intbuff,stream%istart(1:4),stream%icount(1:4))
          else
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%intbuff,stream%istart(1:3),stream%icount(1:3))
          end if
        class is (ncvariable4d_integer)
          if ( docopy ) then
            if ( .not. associated(var%ival) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%intbuff(1:var%totsize) = &
              reshape(var%ival(var%j1:var%j2,var%i1:var%i2,:,:),(/var%totsize/))
          end if
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
          end if
          stream%istart(3) = 1
          stream%icount(3) = var%nval(3)
          stream%istart(4) = 1
          stream%icount(4) = var%nval(4)
          if ( var%lrecords ) then
            stream%istart(5) = stream%irec
            stream%icount(5) = 1
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%intbuff,stream%istart(1:5),stream%icount(1:5))
          else
            ncstat = nf90_put_var(stream%id,var%id, &
              buffer%intbuff,stream%istart(1:4),stream%icount(1:4))
          end if
        class default
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(stream%filename), 1)
      end select
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        write(stderr,*) nf90_strerror(ncstat)
        call die('nc_stream','Cannot write variable '//trim(var%vname)// &
          ' in file '//trim(stream%filename), 1)
      end if
      call outstream_sync(stream)
    end subroutine outstream_writevar

    subroutine cdumlogical(cdum,yesno)
      implicit none
      character(len=*) , intent(out) :: cdum
      logical , intent(in) :: yesno
      if (yesno) then
        write(cdum,'(a)') 'Yes'
      else
        write(cdum,'(a)') 'No'
      end if
    end subroutine cdumlogical

    subroutine add_common_global_params(ncout)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      type(ncstream) , pointer :: stream
      type(basic_variables) , pointer :: stvar
      character(maxstring) :: history
      real(rk8) , dimension(2) :: trlat
      integer(ik4) , dimension(8) :: tvals

      stream => ncout%ncp%xs
      stvar  => ncout%svp%xv
      call add_attribute(stream, &
        ncattribute_string('title','ICTP Regional Climatic model V4'))
      call add_attribute(stream,ncattribute_string('institution','ICTP'))
      call add_attribute(stream, &
        ncattribute_string('source','RegCM Model output file'))
      call add_attribute(stream,ncattribute_string('Conventions','CF-1.4'))
      call add_attribute(stream,ncattribute_string('references', &
        'http://gforge.ictp.it/gf/project/regcm'))
      call add_attribute(stream,ncattribute_string('model_revision',SVN_REV))
      call date_and_time(values=tvals)
      write(history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)') &
        tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,          &
        tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,                &
        ' : Created by RegCM '//trim(stream%progname)//' program'
      call add_attribute(stream,ncattribute_string('history',history))
      call add_attribute(stream,ncattribute_string('experiment',domname))
      call add_attribute(stream,ncattribute_string('projection',iproj))
      if ( stream%l_subgrid ) then
        call add_attribute(stream,ncattribute_real8('grid_size_in_meters', &
          (ds*1000.0)/dble(nsg)))
        call add_attribute(stream,ncattribute_string('model_subgrid','Yes'))
      else
        call add_attribute(stream,ncattribute_real8('grid_size_in_meters', &
          (ds*1000.0)))
      end if
      call add_attribute(stream, &
        ncattribute_real8('latitude_of_projection_origin',clat))
      call add_attribute(stream, &
        ncattribute_real8('longitude_of_projection_origin',clon))
      if ( iproj == 'ROTMER' ) then
        call add_attribute(stream, &
          ncattribute_real8('grid_north_pole_latitude',plat))
        call add_attribute(stream, &
          ncattribute_real8('grid_north_pole_longitude',plon))
      else if ( iproj == 'LAMCON' ) then
        trlat(1) = truelatl
        trlat(2) = truelath
        call add_attribute(stream, &
          ncattribute_real8_array('standard_parallel',trlat,2))
      end if
      call add_attribute(stream,ncattribute_real8('grid_factor',xcone))
      stvar%jx_var = ncvariable1d_real(vname='jx',vunit='m', &
        long_name = 'x-coordinate in Cartesian system', &
        standard_name = 'projection_x_coordinate',      &
        axis = 'x',lrecords = .false.)
      stvar%iy_var = ncvariable1d_real(vname='iy',vunit='m', &
        long_name = 'y-coordinate in Cartesian system', &
        standard_name = 'projection_y_coordinate',      &
        axis = 'y',lrecords = .false.)
      if ( stream%l_full_sigma ) then
        stvar%sigma_var =  ncvariable1d_real(vname='sigma',vunit='1', &
          long_name = "Sigma at full model layers",             &
          standard_name = 'atmosphere_sigma_coordinate' ,       &
          axis = 'z',lrecords = .false.)
      else
        stvar%sigma_var =  ncvariable1d_real(vname='sigma',vunit='1', &
          long_name = "Sigma at half model layers",             &
          standard_name = 'atmosphere_sigma_coordinate' ,       &
          axis = 'z',lrecords = .false.)
      end if
      stvar%ptop_var =  ncvariable0d_real(vname='ptop',vunit='hPa', &
        long_name = "Pressure at model top",                  &
        standard_name = 'air_pressure',lrecords = .false.)
      call outstream_addvar(ncout,stvar%jx_var)
      call outstream_addvar(ncout,stvar%iy_var)
      call outstream_addvar(ncout,stvar%sigma_var)
      call outstream_addvar(ncout,stvar%ptop_var)
      call add_attribute(stream, &
        ncattribute_string('axis','X'),stvar%jx_var%id,stvar%jx_var%vname)
      call add_attribute(stream, &
        ncattribute_string('axis','Y'),stvar%iy_var%id,stvar%iy_var%vname)
      call add_attribute(stream, &
        ncattribute_string('axis','Z'),stvar%sigma_var%id,stvar%sigma_var%vname)
      call add_attribute(stream, &
        ncattribute_string('positive','down'), &
          stvar%sigma_var%id,stvar%sigma_var%vname)
      call add_attribute(stream, &
        ncattribute_string('formula_terms','sigma: sigma ps: ps ptop: ptop'), &
          stvar%sigma_var%id,stvar%sigma_var%vname)
    end subroutine add_common_global_params

end module mod_ncstream

#ifdef TESTNCSTREAM
subroutine myabort
  call abort
end subroutine myabort

!
! Example program to use this module and write a netcdf file
!
program test
  use mod_ncstream
  use mod_dynparam

  type(nc_output_stream) :: ncout
  type(ncstream_params) :: opar

  type(ncvariable0d_integer) :: var0dint
  type(ncvariable1d_real) :: var1dreal
  type(ncvariable2d_real) :: var2dreal
  type(ncvariable3d_real) :: var3dreal

  real(rk8) , dimension(18) :: sigma
  real(rk8) , dimension(36,36) :: 2dvar
  real(rk8) , dimension(36,36,18) :: 3dvar
  integer(ik4) :: k

  data sigma /0.025000000372529, 0.0750000011175871, 0.129999998956919, &
              0.195000000298023, 0.270000003278255, 0.349999994039536, &
              0.429999992251396, 0.510000005364418, 0.590000003576279, &
              0.669999986886978, 0.744999974966049, 0.809999972581863, &
              0.864999979734421, 0.909999996423721, 0.944999992847443, &
              0.969999998807907, 0.985000014305115, 0.995000004768372 /
  jx = 36
  iy = 36
  kz = 18
  ds = 50.0D0
  domname = 'domname'
  calendar = 'gregorian'
  iproj = 'LAMCON'
  xcone = 0.71
  clat = 43.0
  clon = 15.0
  ptop = 5.0
  truelatl = 30.0
  truelath = 60.0

  var0dint%vname = 'int0dvar'
  var0dint%vunit = 'units_int0dvar'
  var0dint%long_name = 'longname int0dvar'
  var0dint%standard_name = 'stdname_int0dvar'
  var0dint%lrecords = .true.

  var1dreal%vname = 'real1dvar'
  var1dreal%vunit = 'units_real1dvar'
  var1dreal%long_name = 'longname real1dvar'
  var1dreal%standard_name = 'stdname_real1dvar'
  var1dreal%lrecords = .false.
  var1dreal%lgridded = .false.
  var1dreal%axis = 'z'

  var2dreal%vname = 'real2dvar'
  var2dreal%vunit = 'units_real2dvar'
  var2dreal%long_name = 'longname real2dvar'
  var2dreal%standard_name = 'stdname_real2dvar'
  var2dreal%lrecords = .false.
  var2dreal%axis = 'xy'

  var3dreal%vname = 'real3dvar'
  var3dreal%vunit = 'units_real3dvar'
  var3dreal%long_name = 'longname real3dvar'
  var3dreal%standard_name = 'stdname_real3dvar'
  var3dreal%lrecords = .true.
  var3dreal%lfillvalue = .true.
  var3dreal%axis = 'xyz'
  var3dreal%cell_method = 'time: mean'

  ! Setup an output stream
  opar%zero_time = 256188.0
  opar%l_bound = .true.
  call outstream_setup(ncout,opar)

  ! Add variables with different dimensions. Can be in a loop !!
  call outstream_addvar(ncout,var0dint)
  call outstream_addvar(ncout,var1dreal)
  call outstream_addvar(ncout,var2dreal)
  call outstream_addvar(ncout,var3dreal)

  ! Add some more attributes
  call outstream_addatt(ncout,ncattribute_string('boundary_smoothing','No'))
  call outstream_addatt(ncout,ncattribute_string('lake_fudging','No'))
  call outstream_addatt(ncout, &
    ncattribute_real8('minimum_h2o_pct_for_water',50.0D0))

  ! Enable the output stream for write
  call outstream_enable(ncout,sigma)

  var1dreal%rval => sigma
  2dvar(1:jx,1:iy) = 1.0D0
  2dvar(jx/2,iy/2) = 2.0D0
  var2dreal%rval => 2dvar

  ! Write some static variables
  call outstream_writevar(ncout,var1dreal)
  call outstream_writevar(ncout,var2dreal)

  var0dint%ival(1) = 12
  3dvar(1:jx,1:iy,:) = 1.0D0
  do k = 1 , kz
    3dvar(jx/4,iy/4,k) = dble(k)
  end do

  ! Write variables in the current record step
  call outstream_addrec(ncout,256212.0)
  call outstream_writevar(ncout,var0dint)
  call outstream_writevar(ncout,var3dreal)
  var3dreal%rval => 3dvar

  var0dint%ival(1) = 13
  3dvar(1:jx,1:iy,:) = 1.0D0
  do k = 1 , kz
    3dvar(jx/4,iy/4,k) = dble(k)*1.5D0
  end do

  ! Add a new record
  call outstream_addrec(ncout,256236.0)
  call outstream_writevar(ncout,var0dint)
  call outstream_writevar(ncout,var3dreal)

  ! Finally, close the file and cleanup all
  call outstream_dispose(ncout)

end program test
#endif

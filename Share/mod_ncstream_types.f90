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

module mod_ncstream_types

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_date

  public

  integer(ik4) , parameter :: ncmaxdims = 16
  integer(ik4) , private , parameter :: maxname = 32
  integer(ik4) , private , parameter :: maxunit = 36
  integer(ik4) , private , parameter :: maxstring = 256
  integer(ik4) , private , parameter :: maxpath = 256

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

  type ncinstream_params
    ! The name of the input file
    character(len=maxstring) :: fname
    ! To enable Parallel I/O. Must use hdf5 library
    integer(ik4) :: mpi_comm = -1
    integer(ik4) :: mpi_info = -1
    ! If parallel I/O, the processor patch indexes on the global grid
    integer(ik4) :: global_jstart , global_jend
    integer(ik4) :: global_istart , global_iend
  end type ncinstream_params

  type ncoutstream_params
    ! The name of the output file
    character(len=maxstring) :: fname = 'regcm_generic.nc'
    ! The name of the program writing the file
    character(len=maxname) :: pname = 'ncoutstream test'
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
    ! Sync the output each timestep
    logical :: l_sync = .false.
    ! Initial time for this run
    type(rcm_time_and_date) :: zero_date = rcm_time_and_date(1,18231,0)
    ! If parallel I/O, the processor patch indexes on the global grid
    integer(ik4) :: global_jstart , global_jend
    integer(ik4) :: global_istart , global_iend
  end type ncoutstream_params

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

  type internal_ibuffer
    integer(ik4) , dimension(:) , allocatable :: intbuff
    real(rk4) , dimension(:) , allocatable :: realbuff
  end type internal_ibuffer

  type obuff_p
    type(internal_obuffer) , pointer :: xb => null()
  end type obuff_p

  type ibuff_p
    type(internal_ibuffer) , pointer :: xb => null()
  end type ibuff_p

  type ncinstream
    character(len=maxpath) :: filename
    logical :: l_parallel = .false.
    integer(ik4) :: id = -1
    integer(ik4) :: timeid = -1
    integer(ik4) :: nrec = -1
    integer(ik4) , dimension(2) :: jparbound
    integer(ik4) , dimension(2) :: iparbound
    integer(ik4) :: global_nj , global_ni , parsize
    real(rk8) , dimension(2) :: xtime = (/dmissval,dmissval/)
    type(rcm_time_and_date) :: refdate
    character(len=maxunit) :: tunit , tcal
    real(rk8) :: deltat = 1.0D0
    integer(ik4) , dimension(5) :: istart , icount , istride
    integer(ik4) :: ndims
    integer(ik4) , allocatable , dimension(:) :: len_dims
  end type ncinstream

  type ncoutstream
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
    integer(ik4) , dimension(ncmaxdims) :: id_dims
    integer(ik4) , dimension(ncmaxdims) :: len_dims
    integer(ik4) :: irec = 0
    ! Implemented up to 4d var with records
    integer(ik4) , dimension(5) :: istart , icount
    ! If using parallel I/O, those are the patch window
    integer(ik4) , dimension(2) :: jparbound
    integer(ik4) , dimension(2) :: iparbound
    integer(ik4) :: global_nj , global_ni , parsize
  end type ncoutstream
!
  type ncattribute_standard
    character(len=maxname) :: aname
  end type ncattribute_standard

  type, extends(ncattribute_standard) :: ncattribute_string
    character(len=maxstring) :: theval
  end type ncattribute_string

  type, extends(ncattribute_standard) :: ncattribute_logical
    logical :: theval
  end type ncattribute_logical

  type, extends(ncattribute_standard) :: ncattribute_integer
    integer(ik4) :: theval
  end type ncattribute_integer

  type, extends(ncattribute_standard) :: ncattribute_real4
    real(rk4) :: theval
  end type ncattribute_real4

  type, extends(ncattribute_standard) :: ncattribute_real8
    real(rk8) :: theval
  end type ncattribute_real8

  type, extends(ncattribute_standard) :: ncattribute_real4_array
    real(rk4) , pointer , dimension(:) :: theval
    integer(ik4) :: numval = 0
  end type ncattribute_real4_array

  type, extends(ncattribute_standard) :: ncattribute_real8_array
    real(rk8) , pointer , dimension(:) :: theval
    integer(ik4) :: numval = 0
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
    real(rk4) :: rmissval = smissval
    integer(ik4) :: imissval = -9999
    integer(ik4) :: ndims = 0
    integer(ik4) , dimension(5) :: idims = (/-1,-1,-1,-1,-1/)
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
    logical :: is_slice = .false.
    real(rk8) , dimension(:,:) , pointer :: rval => null()
    real(rk8) , dimension(:,:,:) , pointer :: rval_slice => null()
  end type ncvariable2d_real

  type, extends(ncvariable_2d) :: ncvariable2d_integer
    logical :: is_slice = .false.
    integer(ik4) , dimension(:,:) , pointer :: ival => null()
    integer(ik4) , dimension(:,:,:) , pointer :: ival_slice => null()
  end type ncvariable2d_integer

  type, extends(ncvariable_standard) :: ncvariable_3d
    logical :: lrecords = .false.
    character(len=3) :: axis = 'xyz'
    integer(ik4) , dimension(3) :: nval = 0
    integer(ik4) :: i1 = -1 , i2 = -1
    integer(ik4) :: j1 = -1 , j2 = -1
    integer(ik4) :: k1 = -1 , k2 = -1
    logical :: is_level = .false.
  end type ncvariable_3d

  type, extends(ncvariable_3d) :: ncvariable3d_real
    logical :: is_slice = .false.
    real(rk8) , dimension(:,:) , pointer :: rval_level => null()
    real(rk8) , dimension(:,:,:) , pointer :: rval => null()
    real(rk8) , dimension(:,:,:,:) , pointer :: rval_slice => null()
  end type ncvariable3d_real

  type, extends(ncvariable_3d) :: ncvariable3d_integer
    logical :: is_slice = .false.
    integer(rk4) , dimension(:,:) , pointer :: ival_level => null()
    integer(rk4) , dimension(:,:,:) , pointer :: ival => null()
    integer(rk4) , dimension(:,:,:,:) , pointer :: ival_slice => null()
  end type ncvariable3d_integer

  type, extends(ncvariable_standard) :: ncvariable_4d
    logical :: lrecords = .false.
    character(len=4) :: axis = 'xyzd'
    integer(ik4) , dimension(4) :: nval = 0
    integer(ik4) :: i1 = -1 , i2 = -1
    integer(ik4) :: j1 = -1 , j2 = -1
    integer(ik4) :: k1 = -1 , k2 = -1
    integer(ik4) :: n1 = -1 , n2 = -1
  end type ncvariable_4d

  type, extends(ncvariable_4d) :: ncvariable4d_real
    real(rk8) , dimension(:,:,:,:) , pointer :: rval => null()
  end type ncvariable4d_real

  type, extends(ncvariable_4d) :: ncvariable4d_integer
    integer(rk4) , dimension(:,:,:,:) , pointer :: ival => null()
  end type ncvariable4d_integer

  type ncoutstream_p
    type(ncoutstream) , pointer :: xs => null()
  end type ncoutstream_p

  type ncinstream_p
    type(ncinstream) , pointer :: xs => null()
  end type ncinstream_p

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
    type(ncoutstream_p) :: ncp
    type(obuff_p) :: obp
    type(basic_variables_p) :: svp
  end type nc_output_stream

  type nc_input_stream
    type(ncinstream_p) :: ncp
    type(ibuff_p) :: ibp
  end type nc_input_stream

  type nc_variable_p
    class(ncvariable_standard) , pointer :: vp => null()
  end type nc_variable_p

  type nc_varlist
    type(nc_variable_p) , dimension(:) , pointer :: vlist
  end type nc_varlist

end module mod_ncstream_types

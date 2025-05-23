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

module mod_ncstream_types

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_date
#ifdef PNETCDF
  use mpi, only: mpi_offset_kind
#endif

  implicit none

  public

  integer(ik4), parameter :: ncmaxdims = 16
  integer(ik4), private, parameter :: maxname = 64
  integer(ik4), private, parameter :: maxunit = 36
  integer(ik4), private, parameter :: maxattarr = 8
  integer(ik4), private, parameter :: maxstring = 512
  integer(ik4), private, parameter :: maxpath = maxstring

  integer(ik4), parameter :: jx_dim          = 1
  integer(ik4), parameter :: iy_dim          = 2
  integer(ik4), parameter :: kz_dim          = 3
  integer(ik4), parameter :: time_dim        = 4
  integer(ik4), parameter :: time_bound_dim  = 5
  integer(ik4), parameter :: texture_dim     = 6
  integer(ik4), parameter :: h2m_level_dim   = 7
  integer(ik4), parameter :: h10m_level_dim  = 8
  integer(ik4), parameter :: h50m_level_dim  = 9
  integer(ik4), parameter :: h100m_level_dim = 10
  integer(ik4), parameter :: h150m_level_dim = 11
  integer(ik4), parameter :: soil_layer_dim  = 12
  integer(ik4), parameter :: water_depth_dim = 13
  integer(ik4), parameter :: months_dim      = 14
  integer(ik4), parameter :: spectral_dim    = 15
  integer(ik4), parameter :: spectral_b_dim  = 16

  type ncinstream_params
    ! The name of the input file
    character(len=maxstring) :: fname
    ! To enable Parallel I/O. Must use hdf5 library
    integer(ik4) :: mpi_comm = -1
    integer(ik4) :: mpi_info = -1
    ! If parallel I/O, the processor patch indexes on the global grid
    integer(ik4) :: global_jstart, global_jend
    integer(ik4) :: global_istart, global_iend
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
    logical :: l_keep = .false.
    ! If boundary values are part of the grids
    ! Note that in the model configuration as of now, the output is on
    ! CROSS points INTERNAL
    ! Instead for the ICBC, the output is on EXTERNAL and the CROSS grid
    ! has an extra unused line and row.
    ! JX = number of DOT zonal points WITH boundary
    ! IY = number of DOT meridional points WITH boundary
    logical :: l_bound = .false.
    ! BAND output has no East or West boundary.
    logical :: l_band  = .false.
    ! CRM output has no boundaries
    logical :: l_crm   = .false.
    ! If this is a subgrid output file.
    logical :: l_subgrid = .false.
    ! If the vertical coordinate is on full sigma levels
    logical :: l_full_sigma = .false.
    ! If the vertical coordinate is pressure levels
    logical :: l_plev = .false.
    ! Sync the output each timestep
    logical :: l_sync = .false.
    ! Has spectral intervals dimensions (radfile)
    logical :: l_specint = .false.
    ! Initial time for this run
    type(rcm_time_and_date) :: zero_date = rcm_time_and_date(1,18231,0)
    ! If parallel I/O, the processor patch indexes on the global grid
    integer(ik4) :: global_jstart, global_jend
    integer(ik4) :: global_istart, global_iend
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
    logical :: lhas1ddouble = .false.
    logical :: lhas2ddouble = .false.
    logical :: lhas3ddouble = .false.
    logical :: lhas4ddouble = .false.
    integer(ik4), dimension(1) :: max1d_int = 0
    integer(ik4), dimension(2) :: max2d_int = 0
    integer(ik4), dimension(3) :: max3d_int = 0
    integer(ik4), dimension(4) :: max4d_int = 0
    integer(ik4), dimension(1) :: max1d_real = 0
    integer(ik4), dimension(2) :: max2d_real = 0
    integer(ik4), dimension(3) :: max3d_real = 0
    integer(ik4), dimension(4) :: max4d_real = 0
    integer(ik4), dimension(1) :: max1d_double = 0
    integer(ik4), dimension(2) :: max2d_double = 0
    integer(ik4), dimension(3) :: max3d_double = 0
    integer(ik4), dimension(4) :: max4d_double = 0
    integer(ik4), dimension(:), allocatable :: intbuff
    real(rk4), dimension(:), allocatable :: realbuff
    real(rk8), dimension(:), allocatable :: doublebuff
  end type internal_obuffer

  type internal_ibuffer
    integer(ik4), dimension(:), allocatable :: intbuff
    real(rk4), dimension(:), allocatable :: realbuff
    real(rk8), dimension(:), allocatable :: doublebuff
  end type internal_ibuffer

  type obuff_p
    type(internal_obuffer), pointer :: xb => null()
  end type obuff_p

  type ibuff_p
    type(internal_ibuffer), pointer :: xb => null()
  end type ibuff_p

  type ncinstream
    character(len=maxpath) :: filename
    logical :: l_parallel = .false.
    integer(ik4) :: id = -1
    integer(ik4) :: timeid = -1
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) :: nrec = -1
#else
    integer(ik4) :: nrec = -1
#endif
    integer(ik4), dimension(2) :: jparbound
    integer(ik4), dimension(2) :: iparbound
    integer(ik4) :: global_nj, global_ni, parsize
    real(rk8), dimension(2) :: xtime = [dmissval,dmissval]
    type(rcm_time_and_date) :: refdate
    character(len=maxunit) :: tunit, tcal
    real(rk8) :: deltat = 1.0_rk8
#ifdef PNETCDF
    integer(kind=mpi_offset_kind), dimension(5) :: istart, icount, istride
#else
    integer(ik4), dimension(5) :: istart, icount, istride
#endif
    integer(ik4) :: ndims
#ifdef PNETCDF
    integer(kind=mpi_offset_kind), allocatable, dimension(:) :: len_dims
#else
    integer(ik4), allocatable, dimension(:) :: len_dims
#endif
  end type ncinstream

  type ncoutstream
    character(len=maxpath) :: filename
    character(len=maxname) :: progname
    integer(ik4) :: id = -1
    !
    ! Defaults is output from model, i.e.:
    !   -) normal grid, not subgrid
    !   -) BAND output, no E/W boundary
    !   -) CRM output no boundaries E/W and N/S
    !   -) no boundary on cross points
    !   -) half sigma level on vertical
    !
    logical :: l_keep = .false.
    logical :: l_bound = .false.
    logical :: l_band  = .false.
    logical :: l_crm   = .false.
    logical :: l_subgrid = .false.
    logical :: l_full_sigma = .false.
    logical :: l_plev = .false.
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
    logical :: l_has2mlev = .false.
    logical :: l_has10mlev = .false.
    logical :: l_has50mlev = .false.
    logical :: l_has100mlev = .false.
    logical :: l_has150mlev = .false.
    logical :: l_hassoillev = .false.
    logical :: l_hasspectral = .false.
    !
    ! Dimension identifiers for 'coded' dimensions
    !
    integer(ik4), dimension(ncmaxdims) :: id_dims
    integer(ik4), dimension(ncmaxdims) :: len_dims
    integer(ik4) :: irec = 0
    ! Implemented up to 4d var with records
#ifdef PNETCDF
    integer(kind=mpi_offset_kind), dimension(5) :: istart, icount
#else
    integer(ik4), dimension(5) :: istart, icount
#endif
    ! If using parallel I/O, those are the patch window
    integer(ik4), dimension(2) :: jparbound
    integer(ik4), dimension(2) :: iparbound
    integer(ik4) :: global_nj, global_ni, parsize
  end type ncoutstream
!
  type ncattribute_standard
    character(len=maxname) :: aname
  end type ncattribute_standard

  type, extends(ncattribute_standard) :: ncattribute_string
    character(len=2*maxstring) :: theval
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
    real(rk4), dimension(maxattarr) :: theval
    integer(ik4) :: numval = 0
  end type ncattribute_real4_array

  type, extends(ncattribute_standard) :: ncattribute_real8_array
    real(rk8), dimension(maxattarr) :: theval
    integer(ik4) :: numval = 0
  end type ncattribute_real8_array

  type ncvariable_standard
    integer(ik4) :: id = -1
    integer(ik4) :: nctype = -1
    character(len=maxname) :: vname
    character(len=maxunit) :: vunit
    character(len=maxstring) :: long_name
    character(len=maxstring) :: standard_name
    character(len=maxstring) :: notes = " "
    character(len=maxstring) :: cell_method = "time: point"
    integer(ik4) :: totsize = 0
    logical :: lgridded = .false.
    logical :: laddmethod = .true.
    logical :: lfillvalue = .false.
    real(rk4) :: rmissval = smissval
    integer(ik4) :: imissval = -9999
    integer(ik4) :: ndims = 0
    integer(ik4), dimension(5) :: idims = -1
    logical :: lrecords = .false.
  end type ncvariable_standard

  type, extends(ncvariable_standard) :: ncvariable_0d
  end type ncvariable_0d

  type, extends(ncvariable_0d) :: ncvariable0d_real
    real(rk4), dimension(1) :: rval = 0.0
  end type ncvariable0d_real

  type, extends(ncvariable_0d) :: ncvariable0d_double
    real(rk8), dimension(1) :: rval = 0.0
  end type ncvariable0d_double

  type, extends(ncvariable_0d) :: ncvariable0d_mixed
    logical :: is_mixed = .true.
    real(rkx), dimension(1) :: rval = 0.0
  end type ncvariable0d_mixed

  type, extends(ncvariable_0d) :: ncvariable0d_integer
    integer(ik4), dimension(1) :: ival = 0
  end type ncvariable0d_integer

  type, extends(ncvariable_0d) :: ncvariable0d_char
    character, dimension(1) :: cval = 'x'
  end type ncvariable0d_char

  type, extends(ncvariable_standard) :: ncvariable_1d
    character(len=1) :: axis = 'x'
    integer(ik4), dimension(1) :: nval = 0
  end type ncvariable_1d

  type, extends(ncvariable_1d) :: ncvariable1d_real
    real(rk4), dimension(:), pointer, contiguous :: rval => null()
  end type ncvariable1d_real

  type, extends(ncvariable_1d) :: ncvariable1d_double
    real(rk8), dimension(:), pointer, contiguous :: rval => null()
  end type ncvariable1d_double

  type, extends(ncvariable_1d) :: ncvariable1d_mixed
    logical :: is_mixed = .true.
    real(rkx), dimension(:), pointer, contiguous :: rval => null()
  end type ncvariable1d_mixed

  type, extends(ncvariable_1d) :: ncvariable1d_integer
    integer(ik4), dimension(:), pointer, contiguous :: ival => null()
  end type ncvariable1d_integer

  type, extends(ncvariable_standard) :: ncvariable_2d
    character(len=2) :: axis = 'xy'
    integer(ik4), dimension(2) :: nval = 0
    integer(ik4) :: i1 = -1, i2 = -1
    integer(ik4) :: j1 = -1, j2 = -1
  end type ncvariable_2d

  type, extends(ncvariable_2d) :: ncvariable2d_real
    logical :: is_slice = .false.
    real(rk4), dimension(:,:), pointer, contiguous :: rval => null()
    real(rk4), dimension(:,:,:), pointer, contiguous :: rval_slice => null()
  end type ncvariable2d_real

  type, extends(ncvariable_2d) :: ncvariable2d_double
    logical :: is_slice = .false.
    real(rk8), dimension(:,:), pointer, contiguous :: rval => null()
    real(rk8), dimension(:,:,:), pointer, contiguous :: rval_slice => null()
  end type ncvariable2d_double

  type, extends(ncvariable_2d) :: ncvariable2d_mixed
    logical :: is_slice = .false.
    logical :: is_mixed = .true.
    real(rkx), dimension(:,:), pointer, contiguous :: rval => null()
    real(rkx), dimension(:,:,:), pointer, contiguous :: rval_slice => null()
  end type ncvariable2d_mixed

  type, extends(ncvariable_2d) :: ncvariable2d_integer
    logical :: is_slice = .false.
    integer(ik4), dimension(:,:), pointer, contiguous :: ival => null()
    integer(ik4), dimension(:,:,:), pointer, contiguous :: ival_slice => null()
  end type ncvariable2d_integer

  type, extends(ncvariable_standard) :: ncvariable_3d
    character(len=3) :: axis = 'xyz'
    integer(ik4), dimension(3) :: nval = 0
    integer(ik4) :: i1 = -1, i2 = -1
    integer(ik4) :: j1 = -1, j2 = -1
    integer(ik4) :: k1 = -1, k2 = -1
  end type ncvariable_3d

  type, extends(ncvariable_3d) :: ncvariable3d_real
    logical :: is_slice = .false.
    real(rk4), dimension(:,:), pointer, contiguous :: rval_level => null()
    real(rk4), dimension(:,:,:), pointer, contiguous :: rval => null()
    real(rk4), dimension(:,:,:,:), pointer, contiguous :: rval_slice => null()
  end type ncvariable3d_real

  type, extends(ncvariable_3d) :: ncvariable3d_double
    logical :: is_slice = .false.
    real(rk8), dimension(:,:), pointer, contiguous :: rval_level => null()
    real(rk8), dimension(:,:,:), pointer, contiguous :: rval => null()
    real(rk8), dimension(:,:,:,:), pointer, contiguous :: rval_slice => null()
  end type ncvariable3d_double

  type, extends(ncvariable_3d) :: ncvariable3d_mixed
    logical :: is_slice = .false.
    logical :: is_mixed = .true.
    real(rkx), dimension(:,:), pointer, contiguous :: rval_level => null()
    real(rkx), dimension(:,:,:), pointer, contiguous :: rval => null()
    real(rkx), dimension(:,:,:,:), pointer, contiguous :: rval_slice => null()
  end type ncvariable3d_mixed

  type, extends(ncvariable_3d) :: ncvariable3d_integer
    logical :: is_slice = .false.
    integer(ik4), dimension(:,:), pointer, contiguous :: ival_level => null()
    integer(ik4), dimension(:,:,:), pointer, contiguous :: ival => null()
    integer(ik4), dimension(:,:,:,:), pointer, contiguous :: ival_slice => null()
  end type ncvariable3d_integer

  type, extends(ncvariable_standard) :: ncvariable_4d
    character(len=4) :: axis = 'xyzd'
    integer(ik4), dimension(4) :: nval = 0
    integer(ik4) :: i1 = -1, i2 = -1
    integer(ik4) :: j1 = -1, j2 = -1
    integer(ik4) :: k1 = -1, k2 = -1
    integer(ik4) :: n1 = -1, n2 = -1
  end type ncvariable_4d

  type, extends(ncvariable_4d) :: ncvariable4d_real
    real(rk4), dimension(:,:,:,:), pointer, contiguous :: rval => null()
  end type ncvariable4d_real

  type, extends(ncvariable_4d) :: ncvariable4d_double
    real(rk8), dimension(:,:,:,:), pointer, contiguous :: rval => null()
  end type ncvariable4d_double

  type, extends(ncvariable_4d) :: ncvariable4d_mixed
    logical :: is_mixed = .true.
    real(rkx), dimension(:,:,:,:), pointer, contiguous :: rval => null()
  end type ncvariable4d_mixed

  type, extends(ncvariable_4d) :: ncvariable4d_integer
    integer(ik4), dimension(:,:,:,:), pointer, contiguous :: ival => null()
  end type ncvariable4d_integer

  type ncoutstream_p
    type(ncoutstream), pointer :: xs => null()
  end type ncoutstream_p

  type ncinstream_p
    type(ncinstream), pointer :: xs => null()
  end type ncinstream_p

  type basic_variables
    type(ncvariable0d_double) :: time_var
    type(ncvariable1d_double) :: tbound_var
    type(ncvariable0d_double) :: ptop_var
    type(ncvariable1d_double) :: sigma_var
    type(ncvariable1d_double) :: ak_var
    type(ncvariable1d_double) :: bk_var
    type(ncvariable1d_double) :: jx_var
    type(ncvariable1d_double) :: iy_var
    type(ncvariable0d_char) :: map_var
    type(ncvariable1d_double) :: lev2m_var
    type(ncvariable1d_double) :: lev10m_var
    type(ncvariable1d_double) :: lev50m_var
    type(ncvariable1d_double) :: lev100m_var
    type(ncvariable1d_double) :: lev150m_var
    type(ncvariable1d_double) :: levsoil_var
    type(ncvariable2d_double) :: levsoilbound_var
    type(ncvariable2d_double) :: spectral_var
  end type basic_variables

  type basic_variables_p
    type(basic_variables), pointer :: xv => null()
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
    class(ncvariable_standard), pointer :: vp => null()
  end type nc_variable_p

  type nc_varlist
    type(nc_variable_p), dimension(:), pointer :: vlist
  end type nc_varlist

end module mod_ncstream_types
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

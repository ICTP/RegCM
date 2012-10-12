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
  integer(ik4) , parameter :: maxunit = 16
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

  integer(ik4) :: i_nc_status
  character(len=maxstring) :: cstr

  type internal_iobuffer
    logical :: lhas1dint = .false.
    logical :: lhas2dint = .false.
    logical :: lhas3dint = .false.
    logical :: lhas4dint = .false.
    logical :: lhas1dreal = .false.
    logical :: lhas2dreal = .false.
    logical :: lhas3dreal = .false.
    logical :: lhas4dreal = .false.
    integer(ik4) , dimension(1) :: max1d_int
    integer(ik4) , dimension(2) :: max2d_int
    integer(ik4) , dimension(3) :: max3d_int
    integer(ik4) , dimension(4) :: max4d_int
    integer(ik4) , dimension(1) :: max1d_real
    integer(ik4) , dimension(2) :: max2d_real
    integer(ik4) , dimension(3) :: max3d_real
    integer(ik4) , dimension(4) :: max4d_real
    integer(ik4) , dimension(:) , allocatable :: intbuff
    integer(ik4) , dimension(:,:) , allocatable :: ibuf2d
    integer(ik4) , dimension(:,:,:) , allocatable :: ibuf3d
    integer(ik4) , dimension(:,:,:,:) , allocatable :: ibuf4d
    real(rk4) , dimension(:) , allocatable :: realbuff
    real(rk4) , dimension(:,:) , allocatable :: rbuf2d
    real(rk4) , dimension(:,:,:) , allocatable :: rbuf3d
    real(rk4) , dimension(:,:,:,:) , allocatable :: rbuf4d
  end type internal_iobuffer

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
    !
    ! Work flags
    !
    logical :: l_sync = .false.
    logical :: l_enabled = .false.
    !
    ! Dimension identifiers for 'coded' dimensions
    !
    integer(ik4) , dimension(maxdims) :: id_dims
    integer(ik4) , dimension(maxdims) :: len_dims
    type(internal_iobuffer) :: buffer
    integer(ik4) :: irec
    ! Implemented up to 4d var with records
    integer(ik4) , dimension(5) :: istart , icount
  end type ncstream
!
  type ncglobal_attribute_standard
    character(len=maxname) :: aname
  end type ncglobal_attribute_standard

  type, extends(ncglobal_attribute_standard) :: ncattribute_string
    character(len=maxstring) :: theval
  end type ncattribute_string

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
    integer(ik4) :: id
    integer(ik4) :: nctype
    character(len=maxname) :: vname
    character(len=maxunit) :: vunit
    character(len=maxstring) :: long_name
    character(len=maxstring) :: standard_name
    character(len=maxstring) :: cell_method = "time: point"
    integer(ik4) :: totsize = 0
    logical :: lgridded = .true.
    logical :: laddmethod = .true.
    logical :: lfillvalue = .false.
    type(ncstream) , pointer :: stream
  end type ncvariable_standard

  type, extends(ncvariable_standard) :: ncvariable_0d
    logical :: lrecords = .false.
  end type ncvariable_0d

  type, extends(ncvariable_0d) :: ncvariable0d_real
    real(rk4) , dimension(1) :: rval
  end type ncvariable0d_real

  type, extends(ncvariable_0d) :: ncvariable0d_integer
    integer(ik4) , dimension(1) :: ival
  end type ncvariable0d_integer

  type, extends(ncvariable_standard) :: ncvariable_1d
    logical :: lrecords = .false.
    character(len=1) :: axis = ' '
    integer(ik4) , dimension(1) :: nval
  end type ncvariable_1d

  type, extends(ncvariable_1d) :: ncvariable1d_real
    real(rk4) , dimension(:) , pointer :: rval
  end type ncvariable1d_real

  type, extends(ncvariable_1d) :: ncvariable1d_integer
    integer(ik4) , dimension(:) , pointer :: ival
  end type ncvariable1d_integer

  type, extends(ncvariable_standard) :: ncvariable_2d
    logical :: lrecords = .false.
    character(len=2) :: axis = '  '
    integer(ik4) , dimension(2) :: nval
  end type ncvariable_2d

  type, extends(ncvariable_2d) :: ncvariable2d_real
    real(rk4) , dimension(:,:) , pointer :: rval
  end type ncvariable2d_real

  type, extends(ncvariable_2d) :: ncvariable2d_integer
    integer(ik4) , dimension(:,:) , pointer :: ival
  end type ncvariable2d_integer

  type, extends(ncvariable_standard) :: ncvariable_3d
    logical :: lrecords = .false.
    character(len=3) :: axis = '   '
    integer(ik4) , dimension(3) :: nval
  end type ncvariable_3d

  type, extends(ncvariable_3d) :: ncvariable3d_real
    real(rk4) , dimension(:,:,:) , pointer :: rval
  end type ncvariable3d_real

  type, extends(ncvariable_3d) :: ncvariable3d_integer
    integer(rk4) , dimension(:,:,:) , pointer :: ival
  end type ncvariable3d_integer

  type, extends(ncvariable_standard) :: ncvariable_4d
    logical :: lrecords = .false.
    character(len=4) :: axis = '    '
    integer(ik4) , dimension(4) :: nval
  end type ncvariable_4d

  type, extends(ncvariable_4d) :: ncvariable4d_real
    real(rk4) , dimension(:,:,:,:) , pointer :: rval
  end type ncvariable4d_real

  type, extends(ncvariable_4d) :: ncvariable4d_integer
    integer(rk4) , dimension(:,:,:,:) , pointer :: ival
  end type ncvariable4d_integer

  integer(ik4) , dimension(maxdims) :: id_dim
  integer(ik4) , dimension(maxdims) :: len_dim

  type ncstream_p
    type(ncstream) , pointer :: xs => null()
  end type ncstream_p

  interface stream_addvar
    module procedure stream_addvar0d_real
    module procedure stream_addvar0d_integer
    module procedure stream_addvar1d_real
    module procedure stream_addvar1d_integer
    module procedure stream_addvar2d_real
    module procedure stream_addvar2d_integer
    module procedure stream_addvar3d_real
    module procedure stream_addvar3d_integer
    module procedure stream_addvar4d_real
    module procedure stream_addvar4d_integer
  end interface stream_addvar

  interface stream_writevar
    module procedure stream_writevar0d_real
    module procedure stream_writevar0d_integer
    module procedure stream_writevar1d_real
    module procedure stream_writevar1d_integer
    module procedure stream_writevar2d_real
    module procedure stream_writevar2d_integer
    module procedure stream_writevar3d_real
    module procedure stream_writevar3d_integer
    module procedure stream_writevar4d_real
    module procedure stream_writevar4d_integer
  end interface stream_writevar

  public :: ncstream_p
  public :: ncvariable0d_real , ncvariable0d_integer
  public :: ncvariable1d_real , ncvariable1d_integer
  public :: ncvariable2d_real , ncvariable2d_integer
  public :: ncvariable3d_real , ncvariable3d_integer
  public :: ncvariable4d_real , ncvariable4d_integer

  public :: stream_setup , stream_enable , stream_dispose
  public :: stream_addvar
  public :: stream_writevar
  public :: stream_addrec

  contains

    subroutine stream_setup(ncp,fname,pname,l_bound,l_subgrid,l_full_sigma)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      character(len=*) , intent(in) :: fname
      character(len=*) , intent(in) :: pname
      logical , intent(in) , optional :: l_bound
      logical , intent(in) , optional :: l_subgrid
      logical , intent(in) , optional :: l_full_sigma
      type(ncstream) , pointer :: stream

      if ( associated(ncp%xs) ) call stream_dispose(ncp)
      allocate(ncp%xs)
      stream => ncp%xs
#ifdef NETCDF4_HDF5
      i_nc_status = nf90_create(fname, &
        ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model),stream%id)
#else
      i_nc_status = nf90_create(fname,nf90_clobber,stream%id)
#endif
      if ( i_nc_status /= nf90_noerr ) then
        write (stderr,*) nf90_strerror(i_nc_status)
        call die('nc_stream','Cannot create file '//trim(fname),1)
      end if
      stream%filename = fname
      stream%progname = pname
      if ( present(l_bound) )      stream%l_bound      = l_bound
      if ( present(l_subgrid) )    stream%l_subgrid    = l_subgrid
      if ( present(l_full_sigma) ) stream%l_full_sigma = l_full_sigma
      stream%id_dims(:) = -1
      stream%len_dims(:) = 0
      stream%buffer%max1d_int(:) = 0
      stream%buffer%max2d_int(:) = 0
      stream%buffer%max3d_int(:) = 0
      stream%buffer%max4d_int(:) = 0
      stream%buffer%max1d_real(:) = 0
      stream%buffer%max2d_real(:) = 0
      stream%buffer%max3d_real(:) = 0
      stream%buffer%max4d_real(:) = 0
      stream%irec = 1
    end subroutine stream_setup

    subroutine stream_dispose(ncp)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncstream) , pointer :: stream

      if ( .not. associated(ncp%xs) ) return
      stream => ncp%xs
      if ( stream%id > 0 ) then
        i_nc_status = nf90_close(stream%id)
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream','Cannot close file '//trim(stream%filename),1)
        end if
      end if
      call deallocate_buffers(stream%buffer)
      deallocate(ncp%xs)
    end subroutine stream_dispose

    subroutine deallocate_buffers(buffer)
      implicit none
      type(internal_iobuffer) , intent(inout) :: buffer
      if ( allocated(buffer%intbuff) )  deallocate(buffer%intbuff)
      if ( allocated(buffer%realbuff) ) deallocate(buffer%realbuff)
      if ( allocated(buffer%ibuf2d) )   deallocate(buffer%ibuf2d)
      if ( allocated(buffer%rbuf2d) )   deallocate(buffer%rbuf2d)
      if ( allocated(buffer%ibuf3d) )   deallocate(buffer%ibuf3d)
      if ( allocated(buffer%rbuf3d) )   deallocate(buffer%rbuf3d)
      if ( allocated(buffer%ibuf4d) )   deallocate(buffer%ibuf4d)
      if ( allocated(buffer%rbuf4d) )   deallocate(buffer%rbuf4d)
    end subroutine deallocate_buffers

    subroutine stream_enable(ncp)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncstream) , pointer :: stream
      integer(ik4) :: maxnum_int , maxnum_real
      if ( .not. associated(ncp%xs) ) return
      stream => ncp%xs
      if ( stream%l_enabled ) return

      call add_common_global_params(ncp)

      i_nc_status = nf90_enddef(stream%id)
      if ( i_nc_status /= nf90_noerr ) then
        write (stderr,*) nf90_strerror(i_nc_status)
        call die('nc_stream','Cannot enable file '//trim(stream%filename),1)
      end if
      !
      ! Allocate buffer space shared by all vars
      !
      maxnum_int  = stream%buffer%max1d_int(1)
      maxnum_real = stream%buffer%max1d_real(1)

      if ( stream%buffer%lhas2dint ) then
        allocate(stream%buffer%rbuf2d(stream%buffer%max2d_int(1), &
                                      stream%buffer%max2d_int(2)))
        maxnum_int = max(product(stream%buffer%max2d_int),maxnum_int)
      end if
      if ( stream%buffer%lhas2dreal ) then
        allocate(stream%buffer%rbuf2d(stream%buffer%max2d_real(1), &
                                      stream%buffer%max2d_real(2)))
        maxnum_real = max(product(stream%buffer%max2d_real),maxnum_real)
      end if
      if ( stream%buffer%lhas3dint ) then
        allocate(stream%buffer%rbuf3d(stream%buffer%max3d_int(1), &
                                      stream%buffer%max3d_int(2), &
                                      stream%buffer%max3d_int(3)))
        maxnum_int = max(product(stream%buffer%max3d_int),maxnum_int)
      end if
      if ( stream%buffer%lhas3dreal ) then
        allocate(stream%buffer%rbuf3d(stream%buffer%max3d_real(1), &
                                      stream%buffer%max3d_real(2), &
                                      stream%buffer%max3d_real(3)))
        maxnum_real = max(product(stream%buffer%max3d_real),maxnum_real)
      end if
      if ( stream%buffer%lhas4dint ) then
        allocate(stream%buffer%ibuf4d(stream%buffer%max4d_int(1), &
                                      stream%buffer%max4d_int(2), &
                                      stream%buffer%max4d_int(3), &
                                      stream%buffer%max4d_int(4)))
        maxnum_int = max(product(stream%buffer%max4d_int),maxnum_int)
      end if
      if ( stream%buffer%lhas4dreal ) then
        allocate(stream%buffer%rbuf4d(stream%buffer%max4d_real(1), &
                                      stream%buffer%max4d_real(2), &
                                      stream%buffer%max4d_real(3), &
                                      stream%buffer%max4d_real(4)))
        maxnum_real = max(product(stream%buffer%max4d_real),maxnum_real)
      end if
      if ( maxnum_int > 0 ) allocate(stream%buffer%intbuff(maxnum_int))
      if ( maxnum_real > 0 ) allocate(stream%buffer%realbuff(maxnum_real))
      stream%l_enabled = .true.
#ifdef DEBUG
      write(stdout,*) 'Enabled netCDF stream ',trim(stream%filename)
      if ( allocated(stream%buffer%intbuff) ) &
        write(stdout,*) 'Total buffer integer size :', &
          size(stream%buffer%intbuff)
      if ( allocated(stream%buffer%realbuff) ) &
        write(stdout,*) 'Total buffer float size   :', &
          size(stream%buffer%realbuff)
      if ( allocated(stream%buffer%ibuf2d) ) &
        write(stdout,*) 'Total buffer int  2d size :',size(stream%buffer%ibuf2d)
      if ( allocated(stream%buffer%rbuf2d) ) &
        write(stdout,*) 'Total buffer real 2d size :',size(stream%buffer%rbuf2d)
      if ( allocated(stream%buffer%ibuf3d) ) &
        write(stdout,*) 'Total buffer int  3d size :',size(stream%buffer%ibuf3d)
      if ( allocated(stream%buffer%rbuf3d) ) &
        write(stdout,*) 'Total buffer real 3d size :',size(stream%buffer%rbuf3d)
      if ( allocated(stream%buffer%ibuf4d) ) &
        write(stdout,*) 'Total buffer int  4d size :',size(stream%buffer%ibuf4d)
      if ( allocated(stream%buffer%rbuf4d) ) &
        write(stdout,*) 'Total buffer real 4d size :',size(stream%buffer%rbuf4d)
#endif
    end subroutine stream_enable

    subroutine stream_sync(stream)
      implicit none
      type(ncstream) , pointer , intent(in) :: stream
      if ( .not. stream%l_enabled ) return
      if ( stream%l_sync ) then
        i_nc_status = nf90_sync(stream%id)
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
                   'Cannot sync file '//trim(stream%filename), 1)
        end if
      end if
    end subroutine stream_sync

    subroutine stream_addrec(ncp)
      implicit none
      type(ncstream_p) , intent(in) :: ncp
      type(ncstream) , pointer :: stream
      if ( .not. associated(ncp%xs) ) return
      stream => ncp%xs
      if ( .not. stream%l_enabled ) return
      stream%irec = stream%irec+1
    end subroutine stream_addrec

    subroutine add_dimension(ncp,dname)
      implicit none
      type(ncstream_p) , intent(in) :: ncp
      character(len=*) , intent(in) :: dname
      type(ncstream) , pointer :: stream
      character(len=maxname) :: the_name
      integer(ik4) :: pdim , num
      if ( .not. associated(ncp%xs) ) return
      stream => ncp%xs
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
          call die('nc_stream', &
            'Cannot add dimension '//trim(dname)//' to file '// &
            trim(stream%filename)//': Undefined in add_dimension', 1)
      end select
      i_nc_status = nf90_def_dim(stream%id,the_name,num,stream%id_dims(pdim))
      if ( i_nc_status /= nf90_noerr ) then
        write (stderr,*) nf90_strerror(i_nc_status)
        call die('nc_stream', &
          'Cannot add dimension '//trim(the_name)//' to file '// &
          trim(stream%filename), 1)
      end if
      stream%len_dims(pdim) = num
    end subroutine add_dimension

    subroutine add_variable(ncp,var,ndims)
      implicit none
      type(ncstream_p) , intent(in) :: ncp
      class(ncvariable_standard) , intent(inout) :: var
      integer(ik4) , intent(in) :: ndims
      type(ncstream) , pointer :: stream
      if ( .not. associated(ncp%xs) ) return
      stream => ncp%xs
      if ( stream%l_enabled ) return
      if ( ndims == 0 ) then
        i_nc_status = nf90_def_var(stream%id,var%vname,var%nctype,var%id)
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot add variable '//trim(var%vname)//' to file '// &
            trim(stream%filename), 1)
        end if
      else
        i_nc_status = nf90_def_var(stream%id,var%vname,var%nctype, &
                              id_dim(1:ndims),var%id)
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot set compression on variable '//trim(var%vname)// &
            ' in file '//trim(stream%filename), 1)
        end if
#ifdef NETCDF4_HDF5
        i_nc_status = nf90_def_var_deflate(stream%id,var%id,1,1,9)
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot set compression on variable '//trim(var%vname)// &
            ' in file '//trim(stream%filename), 1)
        end if
#endif
      end if
      var%stream => stream
      call add_varatts(var)
    end subroutine add_variable

    subroutine add_varatts(var)
      implicit none
      class(ncvariable_standard) :: var
      character(len=16) :: coords_cross = 'xlat xlon'
      character(len=16) :: coords_dot   = 'dlat dlon'
      i_nc_status = nf90_put_att(var%stream%id,var%id,'long_name',var%long_name)
      if ( i_nc_status /= nf90_noerr ) then
        write (stderr,*) nf90_strerror(i_nc_status)
        call die('nc_stream', &
          'Cannot add attribute long_name to variable '//trim(var%vname)// &
          ' in file '//trim(var%stream%filename), 1)
      end if
      i_nc_status = nf90_put_att(var%stream%id,var%id,'standard_name', &
                                 var%standard_name)
      if ( i_nc_status /= nf90_noerr ) then
        write (stderr,*) nf90_strerror(i_nc_status)
        call die('nc_stream', &
          'Cannot add attribute standard_name to variable '//trim(var%vname)// &
          ' in file '//trim(var%stream%filename), 1)
      end if
      i_nc_status = nf90_put_att(var%stream%id,var%id,'units',var%vunit)
      if ( i_nc_status /= nf90_noerr ) then
        write (stderr,*) nf90_strerror(i_nc_status)
        call die('nc_stream', &
          'Cannot add attribute units to variable '//trim(var%vname)// &
          ' in file '//trim(var%stream%filename), 1)
      end if
      if ( var%lgridded ) then
        if ( var%stream%l_bound .and. &
            (var%vname == 'u' .or. var%vname == 'v') ) then
          i_nc_status = nf90_put_att(var%stream%id,var%id, &
                                     'coordinates',coords_dot)
          if ( i_nc_status /= nf90_noerr ) then
            write (stderr,*) nf90_strerror(i_nc_status)
            call die('nc_stream', &
              'Cannot add attribute coordinates to variable '// &
              trim(var%vname)//' in file '//trim(var%stream%filename), 1)
          end if
        else
          i_nc_status = nf90_put_att(var%stream%id,var%id, &
                                     'coordinates',coords_cross)
          if ( i_nc_status /= nf90_noerr ) then
            write (stderr,*) nf90_strerror(i_nc_status)
            call die('nc_stream', &
              'Cannot add attribute coordinates to variable '// &
              trim(var%vname)//' in file '//trim(var%stream%filename), 1)
          end if
        end if
        i_nc_status = nf90_put_att(var%stream%id,var%id, &
                                   'grid_mapping','rcm_map')
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot add attribute grid_mapping to variable '// &
            trim(var%vname)//' in file '//trim(var%stream%filename), 1)
        end if
      end if
      if ( var%laddmethod .and. var%vname(1:4) /= 'time' ) then
        i_nc_status = nf90_put_att(var%stream%id,var%id,'cell_methods', &
                                   var%cell_method)
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot add attribute cell_methods to variable '// &
            trim(var%vname)//' in file '//trim(var%stream%filename), 1)
        end if
      end if
      if ( var%lfillvalue ) then
        i_nc_status = nf90_put_att(var%stream%id,var%id,'_FillValue', &
                                   smissval)
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot add attribute _FillValue to variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
    end subroutine add_varatts

    subroutine add_globalatt(ncp,att)
      implicit none
      type(ncstream_p) , intent(in) :: ncp
      class(ncglobal_attribute_standard) :: att
      type(ncstream) , pointer :: stream
      if ( .not. associated(ncp%xs) ) return
      stream => ncp%xs
      if ( stream%l_enabled ) return
      if ( stream%id < 0 ) return
      select type(att)
        class is (ncattribute_string)
          i_nc_status = nf90_put_att(stream%id,nf90_global,att%aname,att%theval)
        class is (ncattribute_integer)
          i_nc_status = nf90_put_att(stream%id,nf90_global,att%aname,att%theval)
        class is (ncattribute_real4)
          i_nc_status = nf90_put_att(stream%id,nf90_global,att%aname,att%theval)
        class is (ncattribute_real8)
          i_nc_status = nf90_put_att(stream%id,nf90_global,att%aname,att%theval)
        class is (ncattribute_real4_array)
          i_nc_status = nf90_put_att(stream%id,nf90_global, &
            att%aname,att%theval(1:att%numval))
        class is (ncattribute_real8_array)
          i_nc_status = nf90_put_att(stream%id,nf90_global, &
            att%aname,att%theval(1:att%numval))
        class default
          call die('nc_stream', &
                   'Cannot add global attribute '//trim(att%aname)// &
                   ' in file '//trim(stream%filename), 1)
      end select
      if ( i_nc_status /= nf90_noerr ) then
        write (stderr,*) nf90_strerror(i_nc_status)
        call die('nc_stream', &
                 'Cannot add global attribute '//trim(att%aname)// &
                 ' in file '//trim(stream%filename), 1)
      end if
    end subroutine add_globalatt

    subroutine stream_addvar0d_norec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_0d) , intent(inout) :: var
      var%lgridded = .false.
      var%laddmethod = .false.
      call add_variable(ncp,var,0)
    end subroutine stream_addvar0d_norec

    subroutine stream_addvar1d_norec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_1d) , intent(inout) :: var
      call dimlist(ncp,var%axis)
      var%lgridded = .false.
      var%laddmethod = .false.
      call add_variable(ncp,var,1)
    end subroutine stream_addvar1d_norec

    subroutine stream_addvar2d_norec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_2d) , intent(inout) :: var
      call dimlist(ncp,var%axis)
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      var%laddmethod = .false.
      call add_variable(ncp,var,2)
    end subroutine stream_addvar2d_norec

    subroutine stream_addvar3d_norec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_3d) , intent(inout) :: var
      call dimlist(ncp,var%axis)
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      var%laddmethod = .false.
      call add_variable(ncp,var,3)
    end subroutine stream_addvar3d_norec

    subroutine stream_addvar4d_norec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_4d) , intent(inout) :: var
      call dimlist(ncp,var%axis)
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      var%laddmethod = .false.
      call add_variable(ncp,var,4)
    end subroutine stream_addvar4d_norec

    subroutine stream_addvar0d_rec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_0d) , intent(inout) :: var
      call dimlist(ncp,'t')
      var%lgridded = .false.
      call add_variable(ncp,var,1)
    end subroutine stream_addvar0d_rec

    subroutine stream_addvar1d_rec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_1d) , intent(inout) :: var
      call dimlist(ncp,var%axis//'t')
      var%lgridded = .false.
      call add_variable(ncp,var,2)
    end subroutine stream_addvar1d_rec

    subroutine stream_addvar2d_rec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_2d) , intent(inout) :: var
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      call dimlist(ncp,var%axis//'t')
      call add_variable(ncp,var,3)
    end subroutine stream_addvar2d_rec

    subroutine stream_addvar3d_rec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_3d) , intent(inout) :: var
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      call dimlist(ncp,var%axis//'t')
      call add_variable(ncp,var,4)
    end subroutine stream_addvar3d_rec

    subroutine stream_addvar4d_rec(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      class(ncvariable_4d) , intent(inout) :: var
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      call dimlist(ncp,var%axis//'t')
      call add_variable(ncp,var,5)
    end subroutine stream_addvar4d_rec

    subroutine stream_addvar0d_real(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable0d_real) , intent(inout) :: var
      var%nctype = nf90_float
      if ( var%lrecords ) then
        call stream_addvar0d_rec(ncp,var)
      else
        call stream_addvar0d_norec(ncp,var)
      end if
    end subroutine stream_addvar0d_real

    subroutine stream_addvar0d_integer(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable0d_integer) , intent(inout) :: var
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar0d_rec(ncp,var)
      else
        call stream_addvar0d_norec(ncp,var)
      end if
    end subroutine stream_addvar0d_integer

    subroutine stream_addvar1d_real(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable1d_real) , intent(inout) :: var
      type(ncstream) , pointer :: stream
      var%nctype = nf90_float
      if ( var%lrecords ) then
        call stream_addvar1d_rec(ncp,var)
      else
        call stream_addvar1d_norec(ncp,var)
      end if
      var%nval(1) = len_dim(1)
      var%totsize = product(var%nval)
      stream => ncp%xs
      stream%buffer%lhas1dreal = .true.
      stream%buffer%max1d_real(1) = max(stream%buffer%max1d_real(1),len_dim(1))
    end subroutine stream_addvar1d_real

    subroutine stream_addvar1d_integer(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable1d_integer) , intent(inout) :: var
      type(ncstream) , pointer :: stream
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar1d_rec(ncp,var)
      else
        call stream_addvar1d_norec(ncp,var)
      end if
      var%nval(1) = len_dim(1)
      var%totsize = product(var%nval)
      stream => ncp%xs
      stream%buffer%lhas1dint = .true.
      stream%buffer%max1d_int(1) = max(stream%buffer%max1d_int(1),len_dim(1))
    end subroutine stream_addvar1d_integer

    subroutine stream_addvar2d_real(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable2d_real) , intent(inout) :: var
      type(ncstream) , pointer :: stream
      var%nctype = nf90_real
      if ( var%lrecords ) then
        call stream_addvar2d_rec(ncp,var)
      else
        call stream_addvar2d_norec(ncp,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%totsize = product(var%nval)
      stream => ncp%xs
      stream%buffer%lhas2dreal = .true.
      stream%buffer%max2d_real(1) = max(stream%buffer%max2d_real(1),len_dim(1))
      stream%buffer%max2d_real(2) = max(stream%buffer%max2d_real(2),len_dim(2))
    end subroutine stream_addvar2d_real

    subroutine stream_addvar2d_integer(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable2d_integer) , intent(inout) :: var
      type(ncstream) , pointer :: stream
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar2d_rec(ncp,var)
      else
        call stream_addvar2d_norec(ncp,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%totsize = product(var%nval)
      stream => ncp%xs
      stream%buffer%lhas2dint = .true.
      stream%buffer%max2d_int(1) = max(stream%buffer%max2d_int(1),len_dim(1))
      stream%buffer%max2d_int(2) = max(stream%buffer%max2d_int(2),len_dim(2))
    end subroutine stream_addvar2d_integer

    subroutine stream_addvar3d_real(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable3d_real) , intent(inout) :: var
      type(ncstream) , pointer :: stream
      var%nctype = nf90_float
      if ( var%lrecords ) then
        call stream_addvar3d_rec(ncp,var)
      else
        call stream_addvar3d_norec(ncp,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%nval(3) = len_dim(3)
      var%totsize = product(var%nval)
      stream => ncp%xs
      stream%buffer%lhas3dreal = .true.
      stream%buffer%max3d_real(1) = max(stream%buffer%max3d_real(1),len_dim(1))
      stream%buffer%max3d_real(2) = max(stream%buffer%max3d_real(2),len_dim(2))
      stream%buffer%max3d_real(3) = max(stream%buffer%max3d_real(3),len_dim(3))
    end subroutine stream_addvar3d_real

    subroutine stream_addvar3d_integer(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable3d_integer) , intent(inout) :: var
      type(ncstream) , pointer :: stream
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar3d_rec(ncp,var)
      else
        call stream_addvar3d_norec(ncp,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%nval(3) = len_dim(3)
      var%totsize = product(var%nval)
      stream => ncp%xs
      stream%buffer%lhas3dint = .true.
      stream%buffer%max3d_int(1) = max(stream%buffer%max3d_int(1),len_dim(1))
      stream%buffer%max3d_int(2) = max(stream%buffer%max3d_int(2),len_dim(2))
      stream%buffer%max3d_int(3) = max(stream%buffer%max3d_int(3),len_dim(3))
    end subroutine stream_addvar3d_integer

    subroutine stream_addvar4d_real(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable4d_real) , intent(inout) :: var
      type(ncstream) , pointer :: stream
      var%nctype = nf90_float
      if ( var%lrecords ) then
        call stream_addvar4d_rec(ncp,var)
      else
        call stream_addvar4d_norec(ncp,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%nval(3) = len_dim(3)
      var%nval(4) = len_dim(4)
      var%totsize = product(var%nval)
      stream => ncp%xs
      stream%buffer%lhas4dreal = .true.
      stream%buffer%max4d_real(1) = max(stream%buffer%max4d_real(1),len_dim(1))
      stream%buffer%max4d_real(2) = max(stream%buffer%max4d_real(2),len_dim(2))
      stream%buffer%max4d_real(3) = max(stream%buffer%max4d_real(3),len_dim(3))
      stream%buffer%max4d_real(4) = max(stream%buffer%max4d_real(4),len_dim(4))
    end subroutine stream_addvar4d_real

    subroutine stream_addvar4d_integer(ncp,var)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      type(ncvariable4d_integer) , intent(inout) :: var
      type(ncstream) , pointer :: stream
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar4d_rec(ncp,var)
      else
        call stream_addvar4d_norec(ncp,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%nval(3) = len_dim(3)
      var%nval(4) = len_dim(4)
      var%totsize = product(var%nval)
      stream => ncp%xs
      stream%buffer%lhas4dint = .true.
      stream%buffer%max4d_int(1) = max(stream%buffer%max4d_int(1),len_dim(1))
      stream%buffer%max4d_int(2) = max(stream%buffer%max4d_int(2),len_dim(2))
      stream%buffer%max4d_int(3) = max(stream%buffer%max4d_int(3),len_dim(3))
      stream%buffer%max4d_int(4) = max(stream%buffer%max4d_int(4),len_dim(4))
    end subroutine stream_addvar4d_integer

    subroutine dimlist(ncp,code)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      character(len=*) , intent(in) :: code
      type(ncstream) , pointer :: stream
      character(len=maxdims) :: safecode
      integer(ik4) :: ic

      if ( .not. associated(ncp%xs) ) return
      stream => ncp%xs
      do ic = 1 , maxdims
        safecode(ic:ic) = ' '
      end do
      safecode = code
      ic = 1
      do while ( safecode(ic:ic) /= ' ' )
        select case (safecode(ic:ic))
          case ('x')
            if ( stream%id_dims(jx_dim) < 0 ) then
              call add_dimension(ncp,'jx')
            end if
            id_dim(ic) = stream%id_dims(jx_dim)
            len_dim(ic) = stream%len_dims(jx_dim)
          case ('y')
            if ( stream%id_dims(iy_dim) < 0 ) then
              call add_dimension(ncp,'iy')
            end if
            id_dim(ic) = stream%id_dims(iy_dim)
            len_dim(ic) = stream%len_dims(iy_dim)
          case ('z')
            if ( stream%id_dims(kz_dim) < 0 ) then
              call add_dimension(ncp,'kz')
            end if
            id_dim(ic) = stream%id_dims(kz_dim)
            len_dim(ic) = stream%len_dims(kz_dim)
          case ('t')
            if ( stream%id_dims(time_dim) < 0 ) then
              call add_dimension(ncp,'time')
            end if
            id_dim(ic) = stream%id_dims(time_dim)
            len_dim(ic) = stream%len_dims(time_dim)
          case ('b')
            if ( stream%id_dims(time_bound_dim) < 0 ) then
              call add_dimension(ncp,'time_bound')
            end if
            id_dim(ic) = stream%id_dims(time_bound_dim)
            len_dim(ic) = stream%len_dims(time_bound_dim)
          case ('T')
            if ( stream%id_dims(texture_dim) < 0 ) then
              call add_dimension(ncp,'ntex')
            end if
            id_dim(ic) = stream%id_dims(texture_dim)
            len_dim(ic) = stream%len_dims(texture_dim)
          case ('2')
            if ( stream%id_dims(h2m_level_dim) < 0 ) then
              call add_dimension(ncp,'2m')
            end if
            id_dim(ic) = stream%id_dims(h2m_level_dim)
            len_dim(ic) = stream%len_dims(h2m_level_dim)
          case ('w')
            if ( stream%id_dims(h10m_level_dim) < 0 ) then
              call add_dimension(ncp,'h10')
            end if
            id_dim(ic) = stream%id_dims(h10m_level_dim)
            len_dim(ic) = stream%len_dims(h10m_level_dim)
          case ('s')
            if ( stream%id_dims(soil_layer_dim) < 0 ) then
              call add_dimension(ncp,'soil_layer')
            end if
            id_dim(ic) = stream%id_dims(soil_layer_dim)
            len_dim(ic) = stream%len_dims(soil_layer_dim)
          case ('d')
            if ( stream%id_dims(water_depth_dim) < 0 ) then
              call add_dimension(ncp,'depth')
            end if
            id_dim(ic) = stream%id_dims(water_depth_dim)
            len_dim(ic) = stream%len_dims(water_depth_dim)
          case default
            write(stderr,*) 'Not in list. Known dimension codes: xyztbT2wsd'
            call die('nc_stream', &
              'Cannot select dimension '//safecode(ic:ic)//' on file '// &
              trim(stream%filename), 1)
        end select
        ic = ic + 1
      end do
    end subroutine dimlist

    subroutine stream_writevar0d_real(var)
      implicit none
      type(ncvariable0d_real) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      if ( var%lrecords ) then
        var%stream%istart(1) = var%stream%irec
        var%stream%icount(1) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id,var%rval, &
          var%stream%istart(1:1),var%stream%icount(1:1))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        i_nc_status = nf90_put_var(var%stream%id,var%id,var%rval(1))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar0d_real

    subroutine stream_writevar0d_integer(var)
      implicit none
      type(ncvariable0d_integer) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      if ( var%lrecords ) then
        var%stream%istart(1) = var%stream%irec
        var%stream%icount(1) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id,var%ival, &
          var%stream%istart(1:1),var%stream%icount(1:1))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        i_nc_status = nf90_put_var(var%stream%id,var%id,var%ival(1))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar0d_integer

    subroutine stream_writevar1d_real(var)
      implicit none
      type(ncvariable1d_real) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = var%stream%irec
        var%stream%icount(2) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%realbuff, &
          var%stream%istart(1:2),var%stream%icount(1:2))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%realbuff, &
          var%stream%istart(1:1),var%stream%icount(1:1))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar1d_real

    subroutine stream_writevar1d_integer(var)
      implicit none
      type(ncvariable1d_integer) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = var%stream%irec
        var%stream%icount(2) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%intbuff, &
          var%stream%istart(1:2),var%stream%icount(1:2))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%intbuff, &
          var%stream%istart(1:1),var%stream%icount(1:1))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar1d_integer

    subroutine stream_writevar2d_real(var)
      implicit none
      type(ncvariable2d_real) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      var%stream%buffer%realbuff(1:size(var%stream%buffer%rbuf2d)) = &
        reshape(var%stream%buffer%rbuf2d,(/size(var%stream%buffer%rbuf2d)/))
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = var%stream%irec
        var%stream%icount(3) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%realbuff, &
          var%stream%istart(1:3),var%stream%icount(1:3))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%realbuff, &
          var%stream%istart(1:2),var%stream%icount(1:2))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar2d_real

    subroutine stream_writevar2d_integer(var)
      implicit none
      type(ncvariable2d_integer) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      var%stream%buffer%intbuff(1:size(var%stream%buffer%ibuf2d)) = &
        reshape(var%stream%buffer%ibuf2d,(/size(var%stream%buffer%ibuf2d)/))
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = var%stream%irec
        var%stream%icount(3) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%intbuff, &
          var%stream%istart(1:3),var%stream%icount(1:3))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%intbuff, &
          var%stream%istart(1:2),var%stream%icount(1:2))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar2d_integer

    subroutine stream_writevar3d_real(var)
      implicit none
      type(ncvariable3d_real) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      var%stream%buffer%realbuff(1:size(var%stream%buffer%rbuf3d)) = &
        reshape(var%stream%buffer%rbuf3d,(/size(var%stream%buffer%rbuf3d)/))
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = 1
        var%stream%icount(3) = var%nval(3)
        var%stream%istart(4) = var%stream%irec
        var%stream%icount(4) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%realbuff, &
          var%stream%istart(1:4),var%stream%icount(1:4))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = 1
        var%stream%icount(3) = var%nval(3)
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%realbuff, &
          var%stream%istart(1:3),var%stream%icount(1:3))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar3d_real

    subroutine stream_writevar3d_integer(var)
      implicit none
      type(ncvariable3d_integer) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      var%stream%buffer%intbuff(1:size(var%stream%buffer%ibuf3d)) = &
        reshape(var%stream%buffer%ibuf3d,(/size(var%stream%buffer%ibuf3d)/))
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = 1
        var%stream%icount(3) = var%nval(3)
        var%stream%istart(4) = var%stream%irec
        var%stream%icount(4) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%intbuff, &
          var%stream%istart(1:4),var%stream%icount(1:4))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = 1
        var%stream%icount(3) = var%nval(3)
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%intbuff, &
          var%stream%istart(1:3),var%stream%icount(1:3))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar3d_integer

    subroutine stream_writevar4d_real(var)
      implicit none
      type(ncvariable4d_real) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      var%stream%buffer%realbuff(1:size(var%stream%buffer%rbuf4d)) = &
        reshape(var%stream%buffer%rbuf4d,(/size(var%stream%buffer%rbuf4d)/))
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = 1
        var%stream%icount(3) = var%nval(3)
        var%stream%istart(4) = 1
        var%stream%icount(4) = var%nval(4)
        var%stream%istart(5) = var%stream%irec
        var%stream%icount(5) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%realbuff, &
          var%stream%istart(1:5),var%stream%icount(1:5))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = 1
        var%stream%icount(3) = var%nval(3)
        var%stream%istart(4) = 1
        var%stream%icount(4) = var%nval(4)
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%realbuff, &
          var%stream%istart(1:4),var%stream%icount(1:4))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar4d_real

    subroutine stream_writevar4d_integer(var)
      implicit none
      type(ncvariable4d_integer) , intent(inout) :: var
      if ( .not. associated(var%stream) ) then
        call die('nc_stream', &
          'Unassociated stream for variable '//var%vname,1)
      end if
      if ( .not. var%stream%l_enabled ) then
        call die('nc_stream', &
          'File stream '//trim(var%stream%filename)//' not enabled', 1)
      end if
      var%stream%buffer%intbuff(1:size(var%stream%buffer%ibuf4d)) = &
        reshape(var%stream%buffer%ibuf4d,(/size(var%stream%buffer%ibuf4d)/))
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = 1
        var%stream%icount(3) = var%nval(3)
        var%stream%istart(4) = 1
        var%stream%icount(4) = var%nval(4)
        var%stream%istart(5) = var%stream%irec
        var%stream%icount(5) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%intbuff, &
          var%stream%istart(1:5),var%stream%icount(1:5))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      else
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = 1
        var%stream%icount(3) = var%nval(3)
        var%stream%istart(4) = 1
        var%stream%icount(4) = var%nval(4)
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%intbuff, &
          var%stream%istart(1:4),var%stream%icount(1:4))
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream', &
            'Cannot write variable '//trim(var%vname)// &
            ' in file '//trim(var%stream%filename), 1)
        end if
      end if
      call stream_sync(var%stream)
    end subroutine stream_writevar4d_integer

    subroutine cdumlogical(cdum,yesno)
      implicit none
      character(len=*) , intent(out) :: cdum
      logical , intent(in) :: yesno
      if (yesno) then
        write (cdum,'(a)') 'Yes'
      else
        write (cdum,'(a)') 'No'
      end if
    end subroutine cdumlogical

    subroutine add_common_global_params(ncp)
      implicit none
      type(ncstream_p) , intent(inout) :: ncp
      character(maxstring) :: history
      real(rk8) , dimension(2) :: trlat
      integer(ik4) , dimension(8) :: tvals

      call add_globalatt(ncp, &
        ncattribute_string('title','ICTP Regional Climatic model V4'))
      call add_globalatt(ncp,ncattribute_string('institution','ICTP'))
      call add_globalatt(ncp, &
        ncattribute_string('source','RegCM Model output file'))
      call add_globalatt(ncp,ncattribute_string('Conventions','CF-1.4'))
      call add_globalatt(ncp,ncattribute_string('references', &
        'http://gforge.ictp.it/gf/project/regcm'))
      call add_globalatt(ncp,ncattribute_string('model_revision',SVN_REV))
      call date_and_time(values=tvals)
      write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)') &
        tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,          &
        tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,                &
        ' : Created by RegCM '//trim(ncp%xs%progname)//' program'
      call add_globalatt(ncp,ncattribute_string('history',history))
      call add_globalatt(ncp,ncattribute_string('experiment',domname))
      call add_globalatt(ncp,ncattribute_string('projection',iproj))
      if ( ncp%xs%l_subgrid ) then
        call add_globalatt(ncp,ncattribute_real8('grid_size_in_meters', &
          (ds*1000.0)/dble(nsg)))
        call add_globalatt(ncp,ncattribute_string('model_subgrid','Yes'))
      else
        call add_globalatt(ncp,ncattribute_real8('grid_size_in_meters', &
          (ds*1000.0)))
      end if
      call add_globalatt(ncp, &
        ncattribute_real8('latitude_of_projection_origin',clat))
      call add_globalatt(ncp, &
        ncattribute_real8('longitude_of_projection_origin',clon))
      if ( iproj == 'ROTMER' ) then
        call add_globalatt(ncp, &
          ncattribute_real8('grid_north_pole_latitude',plat))
        call add_globalatt(ncp, &
          ncattribute_real8('grid_north_pole_longitude',plon))
      else if ( iproj == 'LAMCON' ) then
        trlat(1) = truelatl
        trlat(2) = truelath
        call add_globalatt(ncp, &
          ncattribute_real8_array('standard_parallel',trlat,2))
      end if
      call add_globalatt(ncp,ncattribute_real8('grid_factor',xcone))
    end subroutine add_common_global_params

end module mod_ncstream

#ifdef TESTNCSTREAM
subroutine myabort
  call abort
end subroutine myabort

program test
  use mod_ncstream
  use mod_dynparam

  type(ncstream_p) :: ncout
  type(ncvariable0d_integer) :: var0dint
  type(ncvariable1d_real) :: var1dreal
  type(ncvariable2d_real) :: var2dreal
  type(ncvariable3d_real) :: var3dreal

  real(rk8) , dimension(18) :: sigma
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
  iproj = 'LAMCON'
  xcone = 0.71
  clat = 43.0
  clon = 15.0
  truelatl = 30.0
  truelath = 60.0

  var0dint%vname = 'time'
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

  call stream_setup(ncout,'testfile.nc','prog',l_bound=.true.)
  call stream_addvar(ncout,var0dint)
  call stream_addvar(ncout,var1dreal)
  call stream_addvar(ncout,var2dreal)
  call stream_addvar(ncout,var3dreal)

  call stream_enable(ncout)

  ncout%xs%buffer%realbuff(1:size(sigma)) = real(sigma(:))
  call stream_writevar(var1dreal)

  ncout%xs%buffer%rbuf2d(1:jx,1:iy) = 1.0
  ncout%xs%buffer%rbuf2d(jx/2,iy/2) = 2.0
  call stream_writevar(var2dreal)

  var0dint%ival(1) = 12
  call stream_writevar(var0dint)

  ncout%xs%buffer%rbuf3d(1:jx,1:iy,:) = 1.0
  do k = 1 , kz
    ncout%xs%buffer%rbuf3d(jx/4,iy/4,k) = k
  end do
  call stream_writevar(var3dreal)
  call stream_addrec(ncout)

  var0dint%ival(1) = 13
  call stream_writevar(var0dint)
  ncout%xs%buffer%rbuf3d(1:jx,1:iy,:) = 1.0
  do k = 1 , kz
    ncout%xs%buffer%rbuf3d(jx/4,iy/4,k) = k*1.5
  end do
  call stream_writevar(var3dreal)

  call stream_dispose(ncout)

end program test

#endif

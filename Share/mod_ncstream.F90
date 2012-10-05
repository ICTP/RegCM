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

  use netcdf
  use mod_realkinds
  use mod_stdio
  use mod_constants
  use mod_memutil
  use mod_dynparam
  use mod_message

  private

  integer(ik4) , parameter :: maxstreams = 32
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
    integer(ik4) , dimension(:) , pointer :: ibuf1d
    integer(ik4) , dimension(:,:) , pointer :: ibuf2d
    integer(ik4) , dimension(:,:,:) , pointer :: ibuf3d
    integer(ik4) , dimension(:,:,:,:) , pointer :: ibuf4d
    real(rk4) , dimension(:) , pointer :: rbuf1d
    real(rk4) , dimension(:,:) , pointer :: rbuf2d
    real(rk4) , dimension(:,:,:) , pointer :: rbuf3d
    real(rk4) , dimension(:,:,:,:) , pointer :: rbuf4d
  end type internal_iobuffer

  type ncstream
    character(len=maxpath) :: filename
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
    ! DO NOT TOUCH THIS.
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
    integer(ik4) , dimension(4) :: istart , icount
  end type ncstream
!
  type ncvariable_standard
    integer(ik4) :: id
    integer(ik4) :: nctype
    character(len=maxname) :: vname
    character(len=maxunit) :: vunit
    character(len=maxstring) :: long_name
    character(len=maxstring) :: standard_name
    character(len=maxstring) :: cell_method = "time: point"
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
!    module procedure stream_writevar3d_real
!    module procedure stream_writevar3d_integer
!    module procedure stream_writevar4d_real
!    module procedure stream_writevar4d_integer
  end interface stream_writevar

  public :: ncstream
  public :: ncvariable0d_real , ncvariable0d_integer
  public :: ncvariable1d_real , ncvariable1d_integer
  public :: ncvariable2d_real , ncvariable2d_integer
  public :: ncvariable3d_real , ncvariable3d_integer
  public :: ncvariable4d_real , ncvariable4d_integer

  public :: stream_setup , stream_enable , stream_dispose
  public :: stream_addvar , stream_addrec
  public :: stream_writevar

  contains

    subroutine stream_setup(stream,fname,l_bound,l_subgrid,l_full_sigma)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      character(len=*) , intent(in) :: fname
      logical , intent(in) , optional :: l_bound
      logical , intent(in) , optional :: l_subgrid
      logical , intent(in) , optional :: l_full_sigma
      type(ncstream) , pointer :: xs

      if ( associated(stream) ) call stream_dispose(stream)
      allocate(xs)
      stream => xs
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
      nullify(stream%buffer%ibuf1d)
      nullify(stream%buffer%ibuf2d)
      nullify(stream%buffer%ibuf3d)
      nullify(stream%buffer%ibuf4d)
      nullify(stream%buffer%rbuf1d)
      nullify(stream%buffer%rbuf2d)
      nullify(stream%buffer%rbuf3d)
      nullify(stream%buffer%rbuf4d)
      stream%irec = 1
    end subroutine stream_setup

    subroutine stream_dispose(stream)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      if ( .not. associated(stream) ) then
        nullify(stream)
        return
      end if
      if ( stream%id > 0 ) then
        i_nc_status = nf90_close(stream%id)
        if ( i_nc_status /= nf90_noerr ) then
          write (stderr,*) nf90_strerror(i_nc_status)
          call die('nc_stream','Cannot close file '//trim(stream%filename),1)
        end if
      end if
      deallocate(stream)
      nullify(stream)
    end subroutine stream_dispose

    subroutine stream_enable(stream)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      if ( .not. associated(stream) ) then
        call die('nc_stream','Void stream to be enabled !!!',1)
      end if
      if ( stream%l_enabled ) return
      i_nc_status = nf90_enddef(stream%id)
      if ( i_nc_status /= nf90_noerr ) then
        write (stderr,*) nf90_strerror(i_nc_status)
        call die('nc_stream','Cannot enable file '//trim(stream%filename),1)
      end if
      !
      ! Allocate buffer space shared by all vars
      !
      if ( stream%buffer%lhas1dint ) then
        call getmem1d(stream%buffer%ibuf1d,1,stream%buffer%max1d_int(1), &
          'nc_stream:stream%buffer%ibuf1d')
      end if
      if ( stream%buffer%lhas1dreal ) then
        call getmem1d(stream%buffer%rbuf1d,1,stream%buffer%max1d_real(1), &
          'nc_stream:stream%buffer%rbuf1d')
      end if
      if ( stream%buffer%lhas2dint ) then
        call getmem2d(stream%buffer%ibuf2d,1,stream%buffer%max2d_int(1), &
          1,stream%buffer%max2d_int(2),'nc_stream:stream%buffer%ibuf2d')
      end if
      if ( stream%buffer%lhas2dreal ) then
        call getmem2d(stream%buffer%rbuf2d,1,stream%buffer%max2d_real(1), &
          1,stream%buffer%max2d_real(2),'nc_stream:stream%buffer%rbuf2d')
      end if
      if ( stream%buffer%lhas3dint ) then
        call getmem3d(stream%buffer%ibuf3d,1,stream%buffer%max3d_int(1), &
          1,stream%buffer%max3d_int(2),1,stream%buffer%max3d_int(3), &
          'nc_stream:stream%buffer%ibuf3d')
      end if
      if ( stream%buffer%lhas3dreal ) then
        call getmem3d(stream%buffer%rbuf3d,1,stream%buffer%max3d_real(1), &
          1,stream%buffer%max3d_real(2),1,stream%buffer%max3d_real(3), &
          'nc_stream:stream%buffer%rbuf3d')
      end if
      if ( stream%buffer%lhas4dint ) then
        call getmem4d(stream%buffer%ibuf4d,1,stream%buffer%max4d_int(1), &
          1,stream%buffer%max4d_int(2),1,stream%buffer%max4d_int(3), &
          1,stream%buffer%max4d_int(4),'nc_stream:stream%buffer%ibuf4d')
      end if
      if ( stream%buffer%lhas4dreal ) then
        call getmem4d(stream%buffer%rbuf4d,1,stream%buffer%max4d_real(1), &
          1,stream%buffer%max4d_real(2),1,stream%buffer%max4d_real(3), &
          1,stream%buffer%max4d_real(4),'nc_stream:stream%buffer%rbuf4d')
      end if
      stream%l_enabled = .true.
#ifdef DEBUG
      write(stderr,*) 'Enabled netCDF stream ',trim(stream%filename)
      if ( associated(stream%buffer%ibuf1d) ) &
        write(stderr,*) 'Total buffer int  1d size :',size(stream%buffer%ibuf1d)
      if ( associated(stream%buffer%rbuf1d) ) &
        write(stderr,*) 'Total buffer real 1d size :',size(stream%buffer%rbuf1d)
      if ( associated(stream%buffer%ibuf2d) ) &
        write(stderr,*) 'Total buffer int  2d size :',size(stream%buffer%ibuf2d)
      if ( associated(stream%buffer%rbuf2d) ) &
        write(stderr,*) 'Total buffer real 2d size :',size(stream%buffer%rbuf2d)
      if ( associated(stream%buffer%ibuf3d) ) &
        write(stderr,*) 'Total buffer int  3d size :',size(stream%buffer%ibuf3d)
      if ( associated(stream%buffer%rbuf3d) ) &
        write(stderr,*) 'Total buffer real 3d size :',size(stream%buffer%rbuf3d)
      if ( associated(stream%buffer%ibuf4d) ) &
        write(stderr,*) 'Total buffer int  4d size :',size(stream%buffer%ibuf4d)
      if ( associated(stream%buffer%rbuf4d) ) &
        write(stderr,*) 'Total buffer real 4d size :',size(stream%buffer%rbuf4d)
#endif
    end subroutine stream_enable

    subroutine stream_sync(stream)
      implicit none
      type(ncstream) , pointer , intent(in) :: stream
      if ( .not. associated(stream) ) return
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

    subroutine stream_addrec(stream)
      type(ncstream) , pointer , intent(in) :: stream
      if ( .not. associated(stream) ) return
      if ( .not. stream%l_enabled ) return
      stream%irec = stream%irec+1
    end subroutine stream_addrec

    subroutine add_dimension(stream,dname)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      character(len=*) , intent(in) :: dname
      character(len=maxname) :: the_name
      integer(ik4) :: pdim , num

      if ( stream%id < 0 ) then
        write (stderr,*) nf90_strerror(i_nc_status)
        call die('nc_stream','File not opened yet',1)
      end if
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

    subroutine add_variable(stream,var,ndims)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_standard) , intent(inout) :: var
      integer(ik4) , intent(in) :: ndims
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

    subroutine stream_addvar0d_norec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_0d) , intent(inout) :: var
      var%lgridded = .false.
      var%laddmethod = .false.
      call add_variable(stream,var,0)
    end subroutine stream_addvar0d_norec

    subroutine stream_addvar1d_norec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_1d) , intent(inout) :: var
      call dimlist(stream,var%axis)
      var%lgridded = .false.
      var%laddmethod = .false.
      call add_variable(stream,var,1)
    end subroutine stream_addvar1d_norec

    subroutine stream_addvar2d_norec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_2d) , intent(inout) :: var
      call dimlist(stream,var%axis)
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      var%laddmethod = .false.
      call add_variable(stream,var,2)
    end subroutine stream_addvar2d_norec

    subroutine stream_addvar3d_norec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_3d) , intent(inout) :: var
      call dimlist(stream,var%axis)
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      var%laddmethod = .false.
      call add_variable(stream,var,3)
    end subroutine stream_addvar3d_norec

    subroutine stream_addvar4d_norec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_4d) , intent(inout) :: var
      call dimlist(stream,var%axis)
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      var%laddmethod = .false.
      call add_variable(stream,var,4)
    end subroutine stream_addvar4d_norec

    subroutine stream_addvar0d_rec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_0d) , intent(inout) :: var
      call dimlist(stream,'t')
      var%lgridded = .false.
      call add_variable(stream,var,1)
    end subroutine stream_addvar0d_rec

    subroutine stream_addvar1d_rec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_1d) , intent(inout) :: var
      call dimlist(stream,var%axis//'t')
      var%lgridded = .false.
      call add_variable(stream,var,2)
    end subroutine stream_addvar1d_rec

    subroutine stream_addvar2d_rec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_2d) , intent(inout) :: var
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      call dimlist(stream,var%axis//'t')
      call add_variable(stream,var,3)
    end subroutine stream_addvar2d_rec

    subroutine stream_addvar3d_rec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_3d) , intent(inout) :: var
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      call dimlist(stream,var%axis//'t')
      call add_variable(stream,var,4)
    end subroutine stream_addvar3d_rec

    subroutine stream_addvar4d_rec(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      class(ncvariable_4d) , intent(inout) :: var
      if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 .and. &
           var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
        var%lgridded = .true.
      end if
      call dimlist(stream,var%axis//'t')
      call add_variable(stream,var,5)
    end subroutine stream_addvar4d_rec

    subroutine stream_addvar0d_real(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable0d_real) , intent(inout) :: var
      var%nctype = nf90_float
      if ( var%lrecords ) then
        call stream_addvar0d_rec(stream,var)
      else
        call stream_addvar0d_norec(stream,var)
      end if
    end subroutine stream_addvar0d_real

    subroutine stream_addvar0d_integer(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable0d_integer) , intent(inout) :: var
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar0d_rec(stream,var)
      else
        call stream_addvar0d_norec(stream,var)
      end if
    end subroutine stream_addvar0d_integer

    subroutine stream_addvar1d_real(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable1d_real) , intent(inout) :: var
      var%nctype = nf90_float
      if ( var%lrecords ) then
        call stream_addvar1d_rec(stream,var)
      else
        call stream_addvar1d_norec(stream,var)
      end if
      var%nval(1) = len_dim(1)
      stream%buffer%lhas1dreal = .true.
      stream%buffer%max1d_real(1) = max(stream%buffer%max1d_real(1),len_dim(1))
    end subroutine stream_addvar1d_real

    subroutine stream_addvar1d_integer(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable1d_integer) , intent(inout) :: var
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar1d_rec(stream,var)
      else
        call stream_addvar1d_norec(stream,var)
      end if
      var%nval(1) = len_dim(1)
      stream%buffer%lhas1dint = .true.
      stream%buffer%max1d_int(1) = max(stream%buffer%max1d_int(1),len_dim(1))
    end subroutine stream_addvar1d_integer

    subroutine stream_addvar2d_real(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable2d_real) , intent(inout) :: var
      var%nctype = nf90_real
      if ( var%lrecords ) then
        call stream_addvar2d_rec(stream,var)
      else
        call stream_addvar2d_norec(stream,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      stream%buffer%lhas2dreal = .true.
      stream%buffer%max2d_real(1) = max(stream%buffer%max2d_real(1),len_dim(1))
      stream%buffer%max2d_real(2) = max(stream%buffer%max2d_real(2),len_dim(2))
    end subroutine stream_addvar2d_real

    subroutine stream_addvar2d_integer(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable2d_integer) , intent(inout) :: var
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar2d_rec(stream,var)
      else
        call stream_addvar2d_norec(stream,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      stream%buffer%lhas2dint = .true.
      stream%buffer%max2d_int(1) = max(stream%buffer%max2d_int(1),len_dim(1))
      stream%buffer%max2d_int(2) = max(stream%buffer%max2d_int(2),len_dim(2))
    end subroutine stream_addvar2d_integer

    subroutine stream_addvar3d_real(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable3d_real) , intent(inout) :: var
      var%nctype = nf90_float
      if ( var%lrecords ) then
        call stream_addvar3d_rec(stream,var)
      else
        call stream_addvar3d_norec(stream,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%nval(3) = len_dim(3)
      stream%buffer%lhas3dreal = .true.
      stream%buffer%max3d_real(1) = max(stream%buffer%max3d_real(1),len_dim(1))
      stream%buffer%max3d_real(2) = max(stream%buffer%max3d_real(2),len_dim(2))
      stream%buffer%max3d_real(3) = max(stream%buffer%max3d_real(3),len_dim(3))
    end subroutine stream_addvar3d_real

    subroutine stream_addvar3d_integer(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable3d_integer) , intent(inout) :: var
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar3d_rec(stream,var)
      else
        call stream_addvar3d_norec(stream,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%nval(3) = len_dim(3)
      stream%buffer%lhas3dint = .true.
      stream%buffer%max3d_int(1) = max(stream%buffer%max3d_int(1),len_dim(1))
      stream%buffer%max3d_int(2) = max(stream%buffer%max3d_int(2),len_dim(2))
      stream%buffer%max3d_int(3) = max(stream%buffer%max3d_int(3),len_dim(3))
    end subroutine stream_addvar3d_integer

    subroutine stream_addvar4d_real(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable4d_real) , intent(inout) :: var
      var%nctype = nf90_float
      if ( var%lrecords ) then
        call stream_addvar4d_rec(stream,var)
      else
        call stream_addvar4d_norec(stream,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%nval(3) = len_dim(3)
      var%nval(4) = len_dim(4)
      stream%buffer%lhas4dreal = .true.
      stream%buffer%max4d_real(1) = max(stream%buffer%max4d_real(1),len_dim(1))
      stream%buffer%max4d_real(2) = max(stream%buffer%max4d_real(2),len_dim(2))
      stream%buffer%max4d_real(3) = max(stream%buffer%max4d_real(3),len_dim(3))
      stream%buffer%max4d_real(4) = max(stream%buffer%max4d_real(4),len_dim(4))
    end subroutine stream_addvar4d_real

    subroutine stream_addvar4d_integer(stream,var)
      implicit none
      type(ncstream) , pointer , intent(inout) :: stream
      type(ncvariable4d_integer) , intent(inout) :: var
      var%nctype = nf90_int
      if ( var%lrecords ) then
        call stream_addvar4d_rec(stream,var)
      else
        call stream_addvar4d_norec(stream,var)
      end if
      var%nval(1) = len_dim(1)
      var%nval(2) = len_dim(2)
      var%nval(3) = len_dim(3)
      var%nval(4) = len_dim(4)
      stream%buffer%lhas4dint = .true.
      stream%buffer%max4d_int(1) = max(stream%buffer%max4d_int(1),len_dim(1))
      stream%buffer%max4d_int(2) = max(stream%buffer%max4d_int(2),len_dim(2))
      stream%buffer%max4d_int(3) = max(stream%buffer%max4d_int(3),len_dim(3))
      stream%buffer%max4d_int(4) = max(stream%buffer%max4d_int(4),len_dim(4))
    end subroutine stream_addvar4d_integer

    subroutine dimlist(stream,code)
      implicit none
      character(len=*) , intent(in) :: code
      type(ncstream) , pointer , intent(inout) :: stream

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
              call add_dimension(stream,'time_bound')
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
          var%stream%buffer%rbuf1d, &
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
          var%stream%buffer%rbuf1d, &
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
          var%stream%buffer%ibuf1d, &
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
          var%stream%buffer%ibuf1d, &
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
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = var%stream%irec
        var%stream%icount(3) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%rbuf2d, &
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
          var%stream%buffer%rbuf2d, &
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
      if ( var%lrecords ) then
        var%stream%istart(1) = 1
        var%stream%icount(1) = var%nval(1)
        var%stream%istart(2) = 1
        var%stream%icount(2) = var%nval(2)
        var%stream%istart(3) = var%stream%irec
        var%stream%icount(3) = 1
        i_nc_status = nf90_put_var(var%stream%id,var%id, &
          var%stream%buffer%ibuf2d, &
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
          var%stream%buffer%ibuf2d, &
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

end module mod_ncstream

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
  use mod_dynparam
  use mod_message
  use mod_date
  use mod_ncstream_types
  use netcdf

  private

  integer(ik4) :: ncstat
  logical , parameter :: nocopy = .false.
  type(rcm_time_and_date) , save :: reference_date
  integer(ik4) , dimension(ncmaxdims) :: id_dim
  integer(ik4) , dimension(ncmaxdims) :: len_dim

  interface outstream_addrec
    module procedure outstream_addrec_date
    module procedure outstream_addrec_value
  end interface

  public :: ncoutstream_params
  public :: nc_output_stream
  public :: ncvariable0d_char
  public :: ncvariable_standard
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

  public :: nc_input_stream
  public :: ncinstream_params
  public :: instream_setup , instream_dispose
  public :: instream_findrec
  public :: instream_readvar

  contains

    subroutine instream_setup(ncin,params)
      implicit none
      type(nc_input_stream) , intent(inout) :: ncin
      type(ncinstream_params) , intent(in) :: params
      type(ncinstream) , pointer :: stream
      integer(ik4) :: iomode = nf90_nowrite
      integer(ik4) :: dimtime , i
      type(rcm_time_and_date) :: tt

      if ( associated(ncin%ncp%xs) ) call instream_dispose(ncin)
      allocate(ncin%ncp%xs)
      stream => ncin%ncp%xs
      stream%filename = params%fname
#ifdef NETCDF4_HDF5
      if ( params%mpi_comm /= -1 ) then
        iomode = ior(nf90_nowrite,nf90_share)
        ncstat = nf90_open(stream%filename,iomode, &
          stream%id,comm=params%mpi_comm,info=params%mpi_info)
        stream%l_parallel = .true.
      else
        ncstat = nf90_open(stream%filename,iomode,stream%id)
      end if
#else
      if ( params%mpi_comm /= -1 ) then
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        call die('nc_stream','Parallel netcdf with Pnetcdf crash',1)
      else
        ncstat = nf90_open(stream%filename,iomode,stream%id)
      end if
#endif
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) nf90_strerror(ncstat)
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        call die('nc_stream','Cannot open file '//trim(stream%filename),1)
      end if
      if ( stream%l_parallel ) then
        stream%jparbound(1) = params%global_jstart
        stream%jparbound(2) = params%global_jend
        stream%iparbound(1) = params%global_istart
        stream%iparbound(2) = params%global_iend
        stream%global_nj = params%global_jend-params%global_jstart+1
        stream%global_ni = params%global_iend-params%global_istart+1
        stream%parsize = stream%global_ni*stream%global_nj
#ifdef DEBUG
        write(stdout,*) 'Parallel I/O enabled.'
        write(stdout,*) 'Processor ', myid, ' window: ', &
          stream%global_nj,' x ' , stream%global_ni, ' at ', &
          stream%jparbound(1) , ',', stream%iparbound(1)
#endif
      end if
      ncstat = nf90_inq_dimid(stream%id,'time',dimtime)
      if ( ncstat == nf90_noerr ) then
        ncstat = nf90_inquire_dimension(stream%id,dimtime,len=stream%nrec)
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) nf90_strerror(ncstat)
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream','Error reading time dimension in '// &
            trim(stream%filename),1)
        end if
        if ( stream%nrec > 0 ) then
          ncstat = nf90_inq_varid(stream%id,'time',stream%timeid)
          if ( ncstat /= nf90_noerr ) then
            write(stderr,*) nf90_strerror(ncstat)
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            call die('nc_stream','Error reading time variable in '// &
              trim(stream%filename),1)
          end if
          ncstat = nf90_get_att(stream%id,stream%timeid,'units',stream%tunit)
          if ( ncstat /= nf90_noerr ) then
            write(stderr,*) nf90_strerror(ncstat)
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            write(stderr,*) 'Assuming hours since 1949-12-01 00:00:00 UTC'
            write(stderr,*) 'for file ',trim(stream%filename)
            stream%tunit = 'hours since 1949-12-01 00:00:00 UTC'
          end if
          ncstat = nf90_get_att(stream%id,stream%timeid,'calendar',stream%tcal)
          if ( ncstat /= nf90_noerr ) then
            write(stderr,*) nf90_strerror(ncstat)
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            write(stderr,*) 'Assuming gregorian calendar'
            write(stderr,*) 'for file ',trim(stream%filename)
            stream%tcal = 'gregorian'
          end if
          stream%refdate = timeval2date(0.0D0,stream%tunit,stream%tcal)
          stream%istart(1) = 1
          stream%icount(1) = stream%nrec
          stream%istride(1) = stream%nrec-1
          ncstat = nf90_get_var(stream%id,stream%timeid,stream%xtime, &
            stream%istart(1:1),stream%icount(1:1),stream%istride(1:1))
          ! Transform to hours since refdate
          tt = timeval2date(stream%xtime(1),stream%tunit,stream%tcal)
          stream%xtime(1) = hourdiff(tt,stream%refdate)
          tt = timeval2date(stream%xtime(2),stream%tunit,stream%tcal)
          stream%xtime(2) = hourdiff(tt,stream%refdate)
          ! Slight inprecision if unit is in "months since", as
          ! we here assume a regular increment in time computed as hours
          ! from the reference day.
          stream%deltat = (stream%xtime(2)-stream%xtime(1))/stream%nrec
        end if
      end if
      ncstat = nf90_inquire(stream%id,nDimensions=stream%ndims)
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) nf90_strerror(ncstat)
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        write(stderr,*) 'Cannot get dimensional infos for file '// &
          trim(stream%filename)
        call die('nc_stream','Cannot setup file '//trim(stream%filename),1)
      end if
      allocate(stream%len_dims(stream%ndims))
      do i = 1 , stream%ndims
        ncstat = nf90_inquire_dimension(stream%id,i,len=stream%len_dims(i))
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) nf90_strerror(ncstat)
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          write(stderr,*) 'Cannot get dimensional infos for file '// &
            trim(stream%filename)
          call die('nc_stream','Cannot setup file '//trim(stream%filename),1)
        end if
      end do
      allocate(ncin%ibp%xb)
      allocate(ncin%ibp%xb%realbuff(8192))
      allocate(ncin%ibp%xb%intbuff(8192))
    end subroutine instream_setup

    subroutine outstream_setup(ncout,params)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      type(ncoutstream_params) , intent(in) :: params
      integer(ik4) :: iomode = nf90_clobber
      type(ncoutstream) , pointer :: stream
      type(rcm_time_and_date) :: tt

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
      tt = params%zero_date
      reference_date      = 1949120100
      call setcal(reference_date,ical)
      call setcal(tt,reference_date)
      stream%zero_time     = hourdiff(tt,reference_date)
      stream%l_bound       = params%l_bound
      stream%l_band        = params%l_band
      stream%l_sync        = params%l_sync
      stream%l_subgrid     = params%l_subgrid
      stream%l_full_sigma  = params%l_full_sigma
      stream%l_hasspectral = params%l_specint
      stream%id_dims(:)   = -1
      stream%len_dims(:)  = 0
      call add_common_global_params(ncout)
    end subroutine outstream_setup

    subroutine instream_dispose(ncin)
      implicit none
      type(nc_input_stream) , intent(inout) :: ncin
      type(ncinstream) , pointer :: stream

      if ( .not. associated(ncin%ncp%xs) ) return
      stream => ncin%ncp%xs
      if ( stream%id > 0 ) then
        ncstat = nf90_close(stream%id)
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          write(stderr,*) nf90_strerror(ncstat)
          call die('nc_stream','Cannot close file '//trim(stream%filename),1)
        end if
      end if
      if ( allocated(stream%len_dims) ) deallocate(stream%len_dims)
      call deallocate_ibuffer(ncin%ibp%xb)
#ifdef DEBUG
      write(stdout,*) 'Closed input stream ',trim(stream%filename)
#endif
      deallocate(ncin%ncp%xs)
      deallocate(ncin%ibp%xb)
    end subroutine instream_dispose

    subroutine outstream_dispose(ncout)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      type(ncoutstream) , pointer :: stream

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
      call deallocate_obuffer(ncout%obp%xb)
#ifdef DEBUG
      write(stdout,*) 'Closed output stream ',trim(stream%filename)
#endif
      deallocate(ncout%ncp%xs)
      deallocate(ncout%obp%xb)
      deallocate(ncout%svp%xv)
    end subroutine outstream_dispose

    subroutine deallocate_obuffer(xbf)
      implicit none
      type(internal_obuffer) , pointer :: xbf
      if ( .not. associated(xbf) ) return
      if ( allocated(xbf%intbuff) )  deallocate(xbf%intbuff)
      if ( allocated(xbf%realbuff) ) deallocate(xbf%realbuff)
    end subroutine deallocate_obuffer

    subroutine deallocate_ibuffer(xbf)
      implicit none
      type(internal_ibuffer) , pointer :: xbf
      if ( .not. associated(xbf) ) return
      if ( allocated(xbf%intbuff) )  deallocate(xbf%intbuff)
      if ( allocated(xbf%realbuff) ) deallocate(xbf%realbuff)
    end subroutine deallocate_ibuffer

    subroutine outstream_enable(ncout,sigma)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      real(rk8) , dimension(:) , pointer , intent(in) :: sigma

      type(ncoutstream) , pointer :: stream
      type(internal_obuffer) , pointer :: buffer
      type(basic_variables) , pointer :: stvar
      integer(ik4) :: maxnum_int , maxnum_real , i
      real(rk8) :: xds
      type(ncattribute_string) :: attc
      type(ncattribute_real8) :: attr
      type(ncattribute_real8_array) :: attra

      if ( .not. associated(ncout%ncp%xs) ) return
      stream => ncout%ncp%xs
      buffer => ncout%obp%xb
      stvar  => ncout%svp%xv
      if ( stream%l_enabled ) return

      if ( stream%l_hasrec ) then
        stvar%time_var%vname = 'time'
        stvar%time_var%vunit = 'hours since 1949-12-01 00:00:00 UTC'
        stvar%time_var%long_name = 'time'
        stvar%time_var%standard_name = 'time'
        stvar%time_var%lrecords = .true.
        call outstream_addvar(ncout,stvar%time_var)
        attc%aname = 'calendar'
        attc%theval = calendar
        call add_attribute(stream,attc,stvar%time_var%id,stvar%time_var%vname)
        if ( stream%l_hastbound ) then
          attc%aname = 'bounds'
          attc%theval = 'time_bnds'
          call add_attribute(stream,attc,stvar%time_var%id,stvar%time_var%vname)
          stvar%tbound_var%vname = 'time_bnds'
          stvar%tbound_var%vunit = 'hours since 1949-12-01 00:00:00 UTC'
          stvar%tbound_var%axis = 'b'
          stvar%tbound_var%long_name = ''
          stvar%tbound_var%standard_name = ''
          stvar%tbound_var%lrecords = .true.
          call outstream_addvar(ncout,stvar%tbound_var)
          attc%aname = 'calendar'
          attc%theval = calendar
          call add_attribute(stream,attc,stvar%tbound_var%id, &
            stvar%tbound_var%vname)
        end if
      end if
      if ( stream%l_hasgrid ) then
        stvar%map_var%vname = 'rcm_map'
        stvar%map_var%vunit = ''
        stvar%map_var%long_name = ''
        stvar%map_var%standard_name = ''
        call outstream_addvar(ncout,stvar%map_var)
        select case (iproj)
          case('LAMCON')
            attc%aname = 'grid_mapping_name'
            attc%theval = 'lambert_conformal_conic'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attra%aname = 'standard_parallel'
            attra%theval(1) = truelatl
            attra%theval(2) = truelath
            attra%numval = 2
            call add_attribute(stream,attra,stvar%map_var%id, &
              stvar%map_var%vname)
            attr%aname = 'longitude_of_central_meridian'
            attr%theval = clon
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'latitude_of_projection_origin'
            attr%theval = clat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
          case('POLSTR')
            attc%aname = 'grid_mapping_name'
            attc%theval = 'stereographic'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'latitude_of_projection_origin'
            attr%theval = clat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'longitude_of_projection_origin'
            attr%theval = clon
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
          case('NORMER')
            attc%aname = 'grid_mapping_name'
            attc%theval = 'mercator'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'standard_parallel'
            attr%theval = clat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'latitude_of_projection_origin'
            attr%theval = clat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'longitude_of_projection_origin'
            attr%theval = clon
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
          case('ROTMER')
            attc%aname = 'grid_mapping_name'
            attc%theval = 'rotated_mercator'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'latitude_of_projection_origin'
            attr%theval = plat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'longitude_of_projection_origin'
            attr%theval = plon
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
        end select
        attc%aname = '_CoordinateTransformType'
        attc%theval = 'Projection'
        call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
        attc%aname = '_CoordinateAxisTypes'
        attc%theval = 'GeoX GeoY'
        call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
      end if
      if ( stream%l_has2mlev ) then
        stvar%lev2m_var%vname = 'm2'
        stvar%lev2m_var%vunit = 'm'
        stvar%lev2m_var%axis = '2'
        stvar%lev2m_var%long_name = 'Height level'
        stvar%lev2m_var%standard_name = 'height'
        stvar%lev2m_var%lrecords = .false.
        call outstream_addvar(ncout,stvar%lev2m_var)
      end if
      if ( stream%l_has10mlev ) then
        stvar%lev10m_var%vname = 'm10'
        stvar%lev10m_var%vunit = 'm'
        stvar%lev10m_var%axis = 'w'
        stvar%lev10m_var%long_name = 'Height level'
        stvar%lev10m_var%standard_name = 'height'
        stvar%lev10m_var%lrecords = .false.
        call outstream_addvar(ncout,stvar%lev10m_var)
      end if
      if ( stream%l_hassoillev ) then
        stvar%levsoil_var%vname = 'soil_layer'
        stvar%levsoil_var%vunit = 'm'
        stvar%levsoil_var%axis = 's'
        stvar%levsoil_var%long_name = 'Soil layer levels'
        stvar%levsoil_var%standard_name = 'root_depth'
        stvar%levsoil_var%lrecords = .false.
        call outstream_addvar(ncout,stvar%levsoil_var)
      end if
      if ( stream%l_hasspectral ) then
        stvar%spectral_var%vname = 'wavelen'
        stvar%spectral_var%vunit = 'm'
        stvar%spectral_var%axis = 'SB'
        stvar%spectral_var%long_name = 'Spectral intervals'
        stvar%spectral_var%standard_name = 'radiation_wavelength'
        stvar%spectral_var%lrecords = .false.
        call outstream_addvar(ncout,stvar%spectral_var)
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
      buffer%realbuff(1) = &
        -real(((dble(stream%len_dims(jx_dim))-d_one)/d_two) * xds * d_1000)
      do i = 2 , stream%len_dims(jx_dim)
        buffer%realbuff(i) = real(dble(buffer%realbuff(i-1))+xds*d_1000)
      end do
      call outstream_writevar(ncout,stvar%jx_var,nocopy)
      buffer%realbuff(1) = &
        -real(((dble(stream%len_dims(iy_dim))-d_one)/d_two) * xds * d_1000)
      do i = 2 , stream%len_dims(iy_dim)
        buffer%realbuff(i) = real(dble(buffer%realbuff(i-1))+xds*d_1000)
      end do
      call outstream_writevar(ncout,stvar%iy_var,nocopy)
      buffer%realbuff(1:size(sigma)) = real(sigma)
      call outstream_writevar(ncout,stvar%sigma_var,nocopy)
      stvar%ptop_var%rval(1) = real(ptop*10.0D0)
      call outstream_writevar(ncout,stvar%ptop_var)
      if ( stream%l_has2mlev ) then
        buffer%realbuff(1) = 2.0
        call outstream_writevar(ncout,stvar%lev2m_var,nocopy)
      end if
      if ( stream%l_has10mlev ) then
        buffer%realbuff(1) = 10.0
        call outstream_writevar(ncout,stvar%lev10m_var,nocopy)
      end if
      if ( stream%l_hassoillev ) then
        buffer%realbuff(1) = 0.10
        buffer%realbuff(2) = 1.00
        call outstream_writevar(ncout,stvar%levsoil_var,nocopy)
      end if
      if ( stream%l_hasspectral ) then
        buffer%realbuff(1) = 0.00000025
        buffer%realbuff(2) = 0.00000069
        buffer%realbuff(3) = 0.00000119
        buffer%realbuff(4) = 0.00000238
        buffer%realbuff(5) = 0.00000069
        buffer%realbuff(6) = 0.00000119
        buffer%realbuff(7) = 0.00000238
        buffer%realbuff(8) = 0.00000400
        call outstream_writevar(ncout,stvar%spectral_var,nocopy)
      end if
    end subroutine outstream_enable

    subroutine outstream_sync(stream)
      implicit none
      type(ncoutstream) , pointer , intent(in) :: stream
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

    subroutine instream_findrec(ncin,dtime,record)
      implicit none
      type(nc_input_stream) , intent(inout) :: ncin
      type(rcm_time_and_date) , intent(in) :: dtime
      real(rk8) , intent(out) :: record
      type(ncinstream) , pointer :: stream
      real(rk8) :: search
      record = -1.0D0
      if ( .not. associated(ncin%ncp%xs) ) return
      stream => ncin%ncp%xs
      if ( stream%nrec > 0 ) then
        search = hourdiff(dtime,stream%refdate)
        if ( search < stream%xtime(1) .or. search > stream%xtime(2) ) then
          return
        end if
        record = 1.0D0+(search-stream%xtime(1))/stream%deltat
      end if
    end subroutine instream_findrec

    subroutine outstream_addrec_date(ncout,dtime)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      type(rcm_time_and_date) , intent(in) :: dtime
      real(rk8) :: val
      val = hourdiff(dtime,reference_date)
      call outstream_addrec_value(ncout,val)
    end subroutine outstream_addrec_date

    subroutine outstream_addrec_value(ncout,val)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      real(rk8) , intent(in) :: val
      type(ncoutstream) , pointer :: stream
      type(basic_variables) , pointer :: stvar
      type(internal_obuffer) , pointer :: buffer
      if ( .not. associated(ncout%ncp%xs) ) return
      stream => ncout%ncp%xs
      stvar  => ncout%svp%xv
      buffer => ncout%obp%xb
      if ( .not. stream%l_enabled ) return
      stream%irec = stream%irec+1
      stvar%time_var%rval = real(val)
      call outstream_writevar(ncout,stvar%time_var)
      if ( stream%l_hastbound ) then
        buffer%realbuff(1) = real(stream%zero_time)
        buffer%realbuff(2) = real(val)
        call outstream_writevar(ncout,stvar%tbound_var,nocopy)
        stream%zero_time = real(val)
      end if
    end subroutine outstream_addrec_value

    subroutine add_dimension(stream,dname)
      implicit none
      type(ncoutstream) , pointer , intent(inout) :: stream
      character(len=*) , intent(in) :: dname
      character(len=16) :: the_name , in_name
      integer(ik4) :: pdim = -1, num
      if ( stream%l_enabled ) return
      if ( stream%id < 0 ) return
      in_name = dname
      select case (in_name)
        case ('JX','jx')
          if ( stream%l_bound .or. stream%l_band ) then
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
          stream%l_has2mlev = .true.
        case ('10M','10m','lev10m','10mlev','lev_10m')
          num = 1
          the_name = 'm10'
          pdim = h10m_level_dim
          stream%l_has10mlev = .true.
        case ('NSOIL','nsoil','n_soil_layer','nlay','layers')
          num = 2
          the_name = 'soil_layer'
          pdim = soil_layer_dim
          stream%l_hassoillev = .true.
        case ('DPTH','dpth','NDEPTH','ndepth','depth','DEPTH')
          num = ndpmax
          the_name = 'depth'
          pdim = water_depth_dim
        case ('MONTHS','months','month','mon','mpy')
          num = 12
          the_name = 'month'
          pdim = months_dim
        case ('SPECTRAL','spectral','spec')
          num = 4
          the_name = 'spec_intv'
          pdim = spectral_dim
          stream%l_hasspectral = .true.
        case ('SPECTRAL_BOUND','spectral_bound','specbound')
          num = 2
          the_name = 'spec_bound'
          pdim = spectral_b_dim
          stream%l_hasspectral = .true.
        case default
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', &
            'Cannot add dimension '//trim(the_name)//' to file '// &
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
      class(ncattribute_standard) , intent(in) :: att
      if ( .not. associated(ncout%ncp%xs) ) return
      call add_attribute(ncout%ncp%xs,att)
    end subroutine outstream_addatt

    subroutine outstream_addvaratt(ncout,var,att)
      implicit none
      type(nc_output_stream) , intent(in) :: ncout
      class(ncvariable_standard) , intent(in) :: var
      class(ncattribute_standard) , intent(in) :: att
      if ( .not. associated(ncout%ncp%xs) ) return
      call add_attribute(ncout%ncp%xs,att,var%id,var%vname)
    end subroutine outstream_addvaratt

    subroutine add_attribute(stream,att,iloc,vname)
      implicit none
      type(ncoutstream) , pointer , intent(in) :: stream
      class(ncattribute_standard) , intent(in) :: att
      integer(ik4) , optional :: iloc
      character(len=4) :: cdum
      character(len=*) , optional :: vname
      character(len=32) :: the_name
      integer(ik4) :: iv
      if ( stream%l_enabled ) return
      if ( stream%id < 0 ) return
      if ( present(iloc) ) then
        iv = iloc
      else
        iv = nf90_global
      end if
      the_name = vname
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
            'Cannot add attribute '//trim(att%aname)//'to variable '// &
            trim(the_name)//' in file '//trim(stream%filename), 1)
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
      type(ncoutstream) , pointer , intent(in) :: stream
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
        if ( var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
          if ( stream%l_bound .and. &
              (var%vname == 'u' .or. var%vname == 'v') ) then
            call add_attribute(stream, &
              ncattribute_string('coordinates',coords_dot),var%id,var%vname)
          else
            call add_attribute(stream, &
              ncattribute_string('coordinates',coords_cross),var%id,var%vname)
          end if
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
        if ( var%nctype == nf90_real ) then
          call add_attribute(stream, &
            ncattribute_real4('_FillValue',var%rmissval),var%id,var%vname)
        else
          call add_attribute(stream, &
            ncattribute_real4('_FillValue',var%imissval),var%id,var%vname)
        end if

      end if
    end subroutine add_varatts

    subroutine add_variable(ncout,var,ndims)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      class(ncvariable_standard) , intent(inout) :: var
      integer(ik4) , intent(in) :: ndims
      type(ncoutstream) , pointer :: stream
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
#if ! defined(NETCDF4_HDF5) && defined (NETCDF4_COMPRESS)
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
      type(ncoutstream) , pointer :: stream
      type(internal_obuffer) , pointer :: buffer
      integer(ik4) :: nd
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
          var%lgridded = .false.
          if ( var%lrecords ) then
            nd = 1
            stream%l_hasrec = .true.
            call dimlist(stream,'t')
          else
            nd = 0
            var%laddmethod = .false.
          end if
          call add_variable(ncout,var,nd)
        class is (ncvariable_1d)
          var%lgridded = .false.
          if ( var%lrecords ) then
            nd = 2
            stream%l_hasrec = .true.
            call dimlist(stream,var%axis//'t')
          else
            nd = 1
            var%laddmethod = .false.
            call dimlist(stream,var%axis)
          end if
          call add_variable(ncout,var,nd)
          var%nval(1) = len_dim(1)
          var%totsize = product(var%nval)
        class is (ncvariable_2d)
          if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 ) then
            var%lgridded = .true.
          end if
          if ( var%lrecords ) then
            nd = 3
            stream%l_hasrec = .true.
            call dimlist(stream,var%axis//'t')
          else
            nd = 2
            var%laddmethod = .false.
            call dimlist(stream,var%axis)
          end if
          call add_variable(ncout,var,nd)
          var%nval(1) = len_dim(1)
          var%nval(2) = len_dim(2)
          var%totsize = product(var%nval)
        class is (ncvariable_3d)
          if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 ) then
            var%lgridded = .true.
          end if
          if ( var%lrecords ) then
            nd = 4
            stream%l_hasrec = .true.
            call dimlist(stream,var%axis//'t')
          else
            nd = 3
            var%laddmethod = .false.
            call dimlist(stream,var%axis)
          end if
          call add_variable(ncout,var,nd)
          var%nval(1) = len_dim(1)
          var%nval(2) = len_dim(2)
          var%nval(3) = len_dim(3)
          var%totsize = product(var%nval)
        class is (ncvariable_4d)
          if ( scan(var%axis,'x') > 0 .and. scan(var%axis,'y') > 0 ) then
            var%lgridded = .true.
          end if
          if ( var%lrecords ) then
            nd = 5
            stream%l_hasrec = .true.
            call dimlist(stream,var%axis//'t')
          else
            nd = 4
            var%laddmethod = .false.
            call dimlist(stream,var%axis)
          end if
          call add_variable(ncout,var,nd)
          var%nval(1) = len_dim(1)
          var%nval(2) = len_dim(2)
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
      type(ncoutstream) , intent(inout) , pointer :: stream
      character(len=*) , intent(in) :: code
      character(len=ncmaxdims) :: safecode
      integer(ik4) :: ic

      do ic = 1 , ncmaxdims
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
              call add_dimension(stream,'10m')
            end if
            id_dim(ic) = stream%id_dims(h10m_level_dim)
            len_dim(ic) = stream%len_dims(h10m_level_dim)
          case ('s')
            if ( stream%id_dims(soil_layer_dim) < 0 ) then
              call add_dimension(stream,'nsoil')
            end if
            id_dim(ic) = stream%id_dims(soil_layer_dim)
            len_dim(ic) = stream%len_dims(soil_layer_dim)
          case ('d')
            if ( stream%id_dims(water_depth_dim) < 0 ) then
              call add_dimension(stream,'depth')
            end if
            id_dim(ic) = stream%id_dims(water_depth_dim)
            len_dim(ic) = stream%len_dims(water_depth_dim)
          case ('M')
            if ( stream%id_dims(months_dim) < 0 ) then
              call add_dimension(stream,'months')
            end if
            id_dim(ic) = stream%id_dims(months_dim)
            len_dim(ic) = stream%len_dims(months_dim)
          case ('S')
            if ( stream%id_dims(spectral_dim) < 0 ) then
              call add_dimension(stream,'spectral')
            end if
            id_dim(ic) = stream%id_dims(spectral_dim)
            len_dim(ic) = stream%len_dims(spectral_dim)
          case ('B')
            if ( stream%id_dims(spectral_b_dim) < 0 ) then
              call add_dimension(stream,'spectral_bound')
            end if
            id_dim(ic) = stream%id_dims(spectral_b_dim)
            len_dim(ic) = stream%len_dims(spectral_b_dim)
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

    subroutine outstream_writevar(ncout,var,lcopy,is)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      class(ncvariable_standard) , intent(inout) :: var
      logical , intent(in) , optional :: lcopy
      integer(ik4) , intent(in) , optional :: is
      type(ncoutstream) , pointer :: stream
      type(internal_obuffer) , pointer :: buffer
      integer(ik4) :: nd , totsize
      logical :: docopy
      if ( .not. associated(ncout%ncp%xs) ) return
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
          stream%istart(1) = 1
          stream%icount(1) = var%nval(1)
          nd = 1
          if ( var%lrecords ) then
            stream%istart(2) = stream%irec
            stream%icount(2) = 1
            nd = 2
          end if
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
        class is (ncvariable2d_real)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
            totsize = var%totsize
          end if
          if ( docopy ) then
            if ( .not. associated(var%rval) .and. &
                 .not. associated(var%rval_slice) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            if ( var%j1 > 0 .and. var%j2 > 0 ) then
              if ( (var%j2-var%j1+1) /= stream%icount(1) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes JX different from file'
                write(stderr,*) 'file     : ', stream%icount(1)
                write(stderr,*) 'internal : ', var%j2-var%j1+1
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%j1 = 1
              var%j2 = stream%icount(1)
            end if
            if ( var%i1 > 0 .and. var%i2 > 0 ) then
              if ( (var%i2-var%i1+1) /= stream%icount(2) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes IY different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%i1 = 1
              var%i2 = stream%icount(2)
            end if
            if ( var%is_slice ) then
              if ( .not. present(is) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Requesting slice without giving index...'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
              buffer%realbuff(1:totsize) =                                   &
                real(reshape(var%rval_slice(var%j1:var%j2,var%i1:var%i2,is), &
                (/totsize/)))
            else
              buffer%realbuff(1:totsize) =                          &
                real(reshape(var%rval(var%j1:var%j2,var%i1:var%i2), &
                (/totsize/)))
            end if
          end if
          nd = 2
          if ( var%lrecords ) then
            stream%istart(3) = stream%irec
            stream%icount(3) = 1
            nd = 3
          end if
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
        class is (ncvariable3d_real)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
            totsize = var%totsize
          end if
          stream%istart(3) = 1
          stream%icount(3) = var%nval(3)
          if ( docopy ) then
            if ( .not. associated(var%rval) .and. &
                 .not. associated(var%rval_level) .and. &
                 .not. associated(var%rval_slice) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            if ( var%j1 > 0 .and. var%j2 > 0 ) then
              if ( (var%j2-var%j1+1) /= stream%icount(1) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes JX different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%j1 = 1
              var%j2 = stream%icount(1)
            end if
            if ( var%i1 > 0 .and. var%i2 > 0 ) then
              if ( (var%i2-var%i1+1) /= stream%icount(2) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes IY different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%i1 = 1
              var%i2 = stream%icount(2)
            end if
            if ( var%k1 > 0 .and. var%k2 > 0 ) then
              if ( (var%k2-var%k1+1) /= stream%icount(3) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes KZ different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%k1 = 1
              var%k2 = stream%icount(3)
            end if
            if ( var%is_slice ) then
              if ( .not. present(is) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Requesting slice without giving index...'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
              buffer%realbuff(1:totsize) = &
                real(reshape(var%rval_slice(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2,is), (/totsize/)))
            else
              buffer%realbuff(1:totsize) = &
                real(reshape(var%rval(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2),(/totsize/)))
            end if
          end if
          nd = 3
          if ( var%lrecords ) then
            stream%istart(4) = stream%irec
            stream%icount(4) = 1
            nd = 4
          end if
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
        class is (ncvariable4d_real)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
            totsize = var%totsize
          end if
          stream%istart(3) = 1
          stream%icount(3) = var%nval(3)
          stream%istart(4) = 1
          stream%icount(4) = var%nval(4)
          if ( docopy ) then
            if ( .not. associated(var%rval) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            if ( var%j1 > 0 .and. var%j2 > 0 ) then
              if ( (var%j2-var%j1+1) /= stream%icount(1) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes JX different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%j1 = 1
              var%j2 = stream%icount(1)
            end if
            if ( var%i1 > 0 .and. var%i2 > 0 ) then
              if ( (var%i2-var%i1+1) /= stream%icount(2) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes IY different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%i1 = 1
              var%i2 = stream%icount(2)
            end if
            if ( var%k1 > 0 .and. var%k2 > 0 ) then
              if ( (var%k2-var%k1+1) /= stream%icount(3) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes KZ different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%k1 = 1
              var%k2 = stream%icount(3)
            end if
            if ( var%n1 > 0 .and. var%n2 > 0 ) then
              if ( (var%n2-var%n1+1) /= stream%icount(4) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes 4 different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%n1 = 1
              var%n2 = stream%icount(4)
            end if
            if ( .not. associated(var%rval) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%realbuff(1:totsize) = &
              real(reshape(var%rval(var%j1:var%j2,var%i1:var%i2, &
                var%k1:var%k2,var%n1:var%n2),(/totsize/)))
          end if
          nd = 4
          if ( var%lrecords ) then
            stream%istart(5) = stream%irec
            stream%icount(5) = 1
            nd = 5
          end if
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
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
          stream%istart(1) = 1
          stream%icount(1) = var%nval(1)
          nd = 1
          if ( var%lrecords ) then
            stream%istart(2) = stream%irec
            stream%icount(2) = 1
            nd = 2
          end if
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
        class is (ncvariable2d_integer)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
            totsize = var%totsize
          end if
          if ( docopy ) then
            if ( .not. associated(var%ival) .and. &
                 .not. associated(var%ival_slice) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            if ( var%j1 > 0 .and. var%j2 > 0 ) then
              if ( (var%j2-var%j1+1) /= stream%icount(1) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes JX different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%j1 = 1
              var%j2 = stream%icount(1)
            end if
            if ( var%i1 > 0 .and. var%i2 > 0 ) then
              if ( (var%i2-var%i1+1) /= stream%icount(2) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes IY different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%i1 = 1
              var%i2 = stream%icount(2)
            end if
            if ( var%is_slice ) then
              if ( .not. present(is) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Requesting slice without giving index...'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
              buffer%intbuff(1:totsize) =                               &
                reshape(var%ival_slice(var%j1:var%j2,var%i1:var%i2,is), &
                (/totsize/))
            else
              buffer%intbuff(1:totsize) =                               &
                reshape(var%ival(var%j1:var%j2,var%i1:var%i2),(/totsize/))
            end if
          end if
          nd = 2
          if ( var%lrecords ) then
            stream%istart(3) = stream%irec
            stream%icount(3) = 1
            nd = 3
          end if
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
        class is (ncvariable3d_integer)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
            totsize = var%totsize
          end if
          stream%istart(3) = 1
          stream%icount(3) = var%nval(3)
          if ( docopy ) then
            if ( .not. associated(var%ival) .and. &
                 .not. associated(var%ival_slice) .and. &
                 .not. associated(var%ival_level) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            if ( var%j1 > 0 .and. var%j2 > 0 ) then
              if ( (var%j2-var%j1+1) /= stream%icount(1) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes JX different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%j1 = 1
              var%j2 = stream%icount(1)
            end if
            if ( var%i1 > 0 .and. var%i2 > 0 ) then
              if ( (var%i2-var%i1+1) /= stream%icount(2) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes IY different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%i1 = 1
              var%i2 = stream%icount(2)
            end if
            if ( var%k1 > 0 .and. var%k2 > 0 ) then
              if ( (var%k2-var%k1+1) /= stream%icount(3) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes KZ different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%k1 = 1
              var%k2 = stream%icount(3)
            end if
            if ( var%is_slice ) then
              if ( .not. present(is) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Requesting slice without giving index...'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
              buffer%intbuff(1:totsize) =                           &
                reshape(var%ival_slice(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2,is), (/totsize/))
            else
              buffer%realbuff(1:totsize) =                    &
                reshape(var%ival(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2),(/totsize/))
            end if
          end if
          nd = 3
          if ( var%lrecords ) then
            stream%istart(4) = stream%irec
            stream%icount(4) = 1
            nd = 4
          end if
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
        class is (ncvariable4d_integer)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize
          else
            stream%istart(1) = 1
            stream%icount(1) = var%nval(1)
            stream%istart(2) = 1
            stream%icount(2) = var%nval(2)
            totsize = var%totsize
          end if
          stream%istart(3) = 1
          stream%icount(3) = var%nval(3)
          stream%istart(4) = 1
          stream%icount(4) = var%nval(4)
          if ( docopy ) then
            if ( .not. associated(var%ival) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            if ( var%j1 > 0 .and. var%j2 > 0 ) then
              if ( (var%j2-var%j1+1) /= stream%icount(1) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes JX different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%j1 = 1
              var%j2 = stream%icount(1)
            end if
            if ( var%i1 > 0 .and. var%i2 > 0 ) then
              if ( (var%i2-var%i1+1) /= stream%icount(2) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes IY different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%i1 = 1
              var%i2 = stream%icount(2)
            end if
            if ( var%k1 > 0 .and. var%k2 > 0 ) then
              if ( (var%k2-var%k1+1) /= stream%icount(3) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes KZ different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%k1 = 1
              var%k2 = stream%icount(3)
            end if
            if ( var%n1 > 0 .and. var%n2 > 0 ) then
              if ( (var%n2-var%n1+1) /= stream%icount(4) ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                write(stderr,*) 'Internal indexes 4 different from file'
                call die('nc_stream','Cannot write variable '// &
                  trim(var%vname)//' in file '//trim(stream%filename), 1)
              end if
            else
              var%n1 = 1
              var%n2 = stream%icount(4)
            end if
            if ( .not. associated(var%ival) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%intbuff(1:totsize) = &
              reshape(var%ival(var%j1:var%j2,var%i1:var%i2, &
                var%k1:var%k2,var%n1:var%n2),(/totsize/))
          end if
          nd = 4
          if ( var%lrecords ) then
            stream%istart(5) = stream%irec
            stream%icount(5) = 1
            nd = 5
          end if
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
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
      type(ncoutstream) , pointer :: stream
      type(basic_variables) , pointer :: stvar
      character(256) :: history
      integer(ik4) , dimension(8) :: tvals
      type(ncattribute_string) :: attc
      type(ncattribute_real8) :: attr
      type(ncattribute_real8_array) :: attra

      stream => ncout%ncp%xs
      stvar  => ncout%svp%xv
      attc%aname = 'title'
      attc%theval = 'ICTP Regional Climatic model V4'
      call add_attribute(stream,attc)
      attc%aname = 'institution'
      attc%theval = 'ICTP'
      call add_attribute(stream,attc)
      attc%aname = 'source'
      attc%theval = 'RegCM Model output file'
      call add_attribute(stream,attc)
      attc%aname = 'Conventions'
      attc%theval = 'CF-1.4'
      call add_attribute(stream,attc)
      attc%aname = 'references'
      attc%theval = 'http://gforge.ictp.it/gf/project/regcm'
      call add_attribute(stream,attc)
      attc%aname = 'model_revision'
      attc%theval = SVN_REV
      call add_attribute(stream,attc)
      call date_and_time(values=tvals)
      write(history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)') &
        tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,          &
        tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,                &
        ' : Created by RegCM '//trim(stream%progname)//' program'
      attc%aname = 'history'
      attc%theval = history
      call add_attribute(stream,attc)
      attc%aname = 'experiment'
      attc%theval = domname
      call add_attribute(stream,attc)
      attc%aname = 'projection'
      attc%theval = iproj
      call add_attribute(stream,attc)
      if ( stream%l_subgrid ) then
        attr%aname = 'grid_size_in_meters'
        attr%theval = (ds*d_1000)/dble(nsg)
        call add_attribute(stream,attr)
        attc%aname = 'model_subgrid'
        attc%theval = 'Yes'
        call add_attribute(stream,attc)
      else
        attr%aname = 'grid_size_in_meters'
        attr%theval = ds*d_1000
        call add_attribute(stream,attr)
      end if
      attr%aname = 'latitude_of_projection_origin'
      attr%theval = clat
      call add_attribute(stream,attr)
      attr%aname = 'longitude_of_projection_origin'
      attr%theval = clon
      call add_attribute(stream,attr)
      if ( iproj == 'ROTMER' ) then
        attr%aname = 'grid_north_pole_latitude'
        attr%theval = plat
        call add_attribute(stream,attr)
        attr%aname = 'grid_north_pole_longitude'
        attr%theval = plon
        call add_attribute(stream,attr)
      else if ( iproj == 'LAMCON' ) then
        attra%aname = 'standard_parallel'
        attra%theval(1) = truelatl
        attra%theval(2) = truelath
        attra%numval = 2
        call add_attribute(stream,attra)
      end if
      attr%aname = 'grid_factor'
      attr%theval = xcone
      call add_attribute(stream,attr)
      stvar%jx_var%vname = 'jx'
      stvar%jx_var%vunit = 'm'
      stvar%jx_var%long_name = 'x-coordinate in Cartesian system'
      stvar%jx_var%standard_name = 'projection_x_coordinate'
      stvar%jx_var%axis = 'x'
      stvar%iy_var%vname = 'iy'
      stvar%iy_var%vunit = 'm'
      stvar%iy_var%long_name = 'y-coordinate in Cartesian system'
      stvar%iy_var%standard_name = 'projection_y_coordinate'
      stvar%iy_var%axis = 'y'
      stvar%sigma_var%vname = 'sigma'
      stvar%sigma_var%vunit = '1'
      if ( stream%l_full_sigma ) then
        stvar%sigma_var%long_name = "Sigma at full model layers"
      else
        stvar%sigma_var%long_name = "Sigma at half model layers"
      end if
      stvar%sigma_var%standard_name = 'atmosphere_sigma_coordinate'
      stvar%sigma_var%axis = 'z'
      stvar%ptop_var%vname='ptop'
      stvar%ptop_var%vunit='hPa'
      stvar%ptop_var%long_name = "Pressure at model top"
      stvar%ptop_var%standard_name = 'air_pressure'
      call outstream_addvar(ncout,stvar%jx_var)
      call outstream_addvar(ncout,stvar%iy_var)
      call outstream_addvar(ncout,stvar%sigma_var)
      call outstream_addvar(ncout,stvar%ptop_var)
      attc%aname = 'axis'
      attc%theval = 'X'
      call add_attribute(stream,attc,stvar%jx_var%id,stvar%jx_var%vname)
      attc%theval = 'Y'
      call add_attribute(stream,attc,stvar%iy_var%id,stvar%iy_var%vname)
      attc%theval = 'Z'
      call add_attribute(stream,attc,stvar%sigma_var%id,stvar%sigma_var%vname)
      attc%aname = '_CoordinateAxisType'
      attc%theval = 'GeoX'
      call add_attribute(stream,attc,stvar%jx_var%id,stvar%jx_var%vname)
      attc%aname = '_CoordinateAxisType'
      attc%theval = 'GeoY'
      call add_attribute(stream,attc,stvar%iy_var%id,stvar%iy_var%vname)
      attc%aname = 'positive'
      attc%theval = 'down'
      call add_attribute(stream,attc,stvar%sigma_var%id,stvar%sigma_var%vname)
      attc%aname = 'formula_terms'
      attc%theval = 'sigma: sigma ps: ps ptop: ptop'
      call add_attribute(stream,attc,stvar%sigma_var%id,stvar%sigma_var%vname)
      attc%aname = '_CoordinateAxisType'
      attc%theval = 'GeoZ'
      call add_attribute(stream,attc,stvar%sigma_var%id,stvar%sigma_var%vname)
    end subroutine add_common_global_params

    subroutine instream_readvar(ncin,var,irec,window)
      implicit none
      type(nc_input_stream) , intent(inout) :: ncin
      class(ncvariable_standard) , intent(inout) :: var
      integer(ik4) , intent(in) , optional :: irec
      integer(ik4) , dimension(:) , intent(in) , optional :: window
      type(ncinstream) , pointer :: stream
      type(internal_ibuffer) , pointer :: buffer
      integer(ik4) :: ndims

      if ( .not. associated(ncin%ncp%xs) ) return
      stream => ncin%ncp%xs
      buffer => ncin%ibp%xb
      if ( stream%id < 0 ) return
      if ( var%id < 0 ) then
        ncstat = nf90_inq_varid(stream%id,var%vname,var%id)
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) nf90_strerror(ncstat)
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream','Cannot find variable '//trim(var%vname)// &
            ' in '//trim(stream%filename),1)
        end if
        ncstat = nf90_inquire_variable(stream%id,var%id,ndims=var%ndims, &
          dimids=var%idims)
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) nf90_strerror(ncstat)
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream','Cannot inquire variable '//trim(var%vname)// &
            ' in '//trim(stream%filename),1)
        end if
      end if
      select type(var)
        class is (ncvariable0d_real)
          ncstat = nf90_get_var(stream%id,var%id,var%rval)
        class is (ncvariable1d_real)
          ncstat = nf90_get_var(stream%id,var%id,var%rval)
        class is (ncvariable2d_real)
          ndims = 2
          if ( present(window) ) then
            if ( window(1) < 0 .or. window(2) < 0 .or. &
                 window(1) > window(2) .or. &
                 window(3) < 0 .or. window(4) < 0 .or. &
                 window(3) > window(4) .or. &
                 window(1) > stream%len_dims(var%idims(1)) .or. &
                 window(2) > stream%len_dims(var%idims(1)) .or. &
                 window(3) > stream%len_dims(var%idims(2)) .or. &
                 window(4) > stream%len_dims(var%idims(2)) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes out of range'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            if ( stream%l_parallel .and. var%lgridded ) then
              if ( stream%jparbound(1) > window(2) .or. &
                   stream%jparbound(2) < window(1) .or. &
                   stream%iparbound(1) > window(4) .or. &
                   stream%iparbound(2) < window(3) ) then
                ! We are not requested to read (out of our window)
                return
              end if
              stream%istart(1) = max(stream%jparbound(1),window(1))
              stream%istart(2) = max(stream%iparbound(1),window(3))
              if ( stream%jparbound(2) > window(2) ) then
                var%nval(1) = window(2)-stream%istart(1)+1
              else
                var%nval(1) = stream%jparbound(2)-stream%istart(1)+1
              end if
              if ( stream%iparbound(2) > window(4) ) then
                var%nval(2) = window(4)-stream%istart(2)+1
              else
                var%nval(2) = stream%iparbound(2)-stream%istart(2)+1
              end if
            else
              stream%istart(1) = window(1)
              stream%istart(2) = window(3)
              var%nval(1) = window(2) - window(1) + 1
              var%nval(2) = window(4) - window(3) + 1
              if ( var%j1 < 0 .and. var%j2 < 0 ) then
                var%j1 = window(1)
                var%j2 = window(2)
              end if
              if ( var%i1 < 0 .and. var%i2 < 0 ) then
                var%i1 = window(3)
                var%i2 = window(4)
              end if
            end if
          else
            if ( stream%l_parallel .and. var%lgridded ) then
              stream%istart(1) = stream%jparbound(1)
              stream%istart(2) = stream%iparbound(1)
              var%nval(1) = stream%global_nj
              var%nval(2) = stream%global_ni
            else
              stream%istart(1) = 1
              stream%istart(2) = 1
              var%nval(1) = stream%len_dims(var%idims(1))
              var%nval(2) = stream%len_dims(var%idims(2))
            end if
          end if
          var%totsize = product(var%nval)
          stream%icount(1) = var%nval(1)
          stream%icount(2) = var%nval(2)
          if ( present(irec) ) then
            if ( var%ndims /= 3 ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested variable is not time dependent'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            stream%istart(3) = irec
            stream%icount(3) = 1
            ndims = 3
          end if
          if ( size(buffer%realbuff) < var%totsize ) then
            deallocate(buffer%realbuff)
            allocate(buffer%realbuff(var%totsize))
          end if
          ncstat = nf90_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
          if ( var%j1 > 0 .and. var%j2 > 0 .and. &
               var%i1 > 0 .and. var%i2 > 0 ) then
            if ( (var%j2-var%j1+1) /= var%nval(1) .or. &
                 (var%i2-var%i1+1) /= var%nval(2) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes different from file'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            var%rval(var%j1:var%j2,var%i1:var%i2) = &
              reshape(buffer%realbuff(1:var%totsize),stream%icount(1:2))
          else
            var%rval(1:var%nval(1),1:var%nval(2)) = &
              reshape(buffer%realbuff(1:var%totsize),stream%icount(1:2))
          end if
        class is (ncvariable3d_real)
          ndims = 3
          if ( present(window) ) then
            if ( window(1) < 0 .or. window(2) < 0 .or. &
                 window(1) > window(2) .or. &
                 window(3) < 0 .or. window(4) < 0 .or. &
                 window(3) > window(4) .or. &
                 window(5) < 0 .or. window(6) < 0 .or. &
                 window(5) > window(6) .or. &
                 window(1) > stream%len_dims(var%idims(1)) .or. &
                 window(2) > stream%len_dims(var%idims(1)) .or. &
                 window(3) > stream%len_dims(var%idims(2)) .or. &
                 window(4) > stream%len_dims(var%idims(2)) .or. &
                 window(5) > stream%len_dims(var%idims(3)) .or. &
                 window(6) > stream%len_dims(var%idims(3)) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes out of range'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            if ( stream%l_parallel .and. var%lgridded ) then
              if ( stream%jparbound(1) > window(2) .or. &
                   stream%jparbound(2) < window(1) .or. &
                   stream%iparbound(1) > window(4) .or. &
                   stream%iparbound(2) < window(3) ) then
                ! We are not requested to read (out of our window)
                return
              end if
              stream%istart(1) = max(stream%jparbound(1),window(1))
              stream%istart(2) = max(stream%iparbound(1),window(3))
              if ( stream%jparbound(2) > window(2) ) then
                var%nval(1) = window(2)-stream%istart(1)+1
              else
                var%nval(1) = stream%jparbound(2)-stream%istart(1)+1
              end if
              if ( stream%iparbound(2) > window(4) ) then
                var%nval(2) = window(4)-stream%istart(2)+1
              else
                var%nval(2) = stream%iparbound(2)-stream%istart(2)+1
              end if
            else
              stream%istart(1) = window(1)
              stream%istart(2) = window(3)
              var%nval(1) = window(2) - window(1) + 1
              var%nval(2) = window(4) - window(3) + 1
              if ( var%j1 < 0 .and. var%j2 < 0 ) then
                var%j1 = window(1)
                var%j2 = window(2)
              end if
              if ( var%i1 < 0 .and. var%i2 < 0 ) then
                var%i1 = window(3)
                var%i2 = window(4)
              end if
            end if
            stream%istart(3) = window(5)
            var%nval(3) = window(6) - window(5) + 1
            if ( var%k1 < 0 .and. var%k2 < 0 ) then
              var%k1 = window(5)
              var%k2 = window(6)
            end if
          else
            if ( stream%l_parallel .and. var%lgridded ) then
              stream%istart(1) = stream%jparbound(1)
              stream%istart(2) = stream%iparbound(1)
              var%nval(1) = stream%global_nj
              var%nval(2) = stream%global_ni
            else
              stream%istart(1) = 1
              stream%istart(2) = 1
              var%nval(1) = stream%len_dims(var%idims(1))
              var%nval(2) = stream%len_dims(var%idims(2))
            end if
            stream%istart(3) = 1
            var%nval(3) = stream%len_dims(var%idims(3))
          end if
          var%totsize = product(var%nval)
          stream%icount(1) = var%nval(1)
          stream%icount(2) = var%nval(2)
          stream%icount(3) = var%nval(3)
          if ( present(irec) ) then
            if ( var%ndims /= 4 ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested variable is not time dependent'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            stream%istart(4) = irec
            stream%icount(4) = 1
            ndims = 4
          end if
          if ( size(buffer%realbuff) < var%totsize ) then
            deallocate(buffer%realbuff)
            allocate(buffer%realbuff(var%totsize))
          end if
          ncstat = nf90_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
          if ( var%j1 > 0 .and. var%j2 > 0 .and. &
               var%i1 > 0 .and. var%i2 > 0 .and. &
               var%k1 > 0 .and. var%k2 > 0 ) then
            if ( (var%j2-var%j1+1) /= var%nval(1) .or. &
                 (var%i2-var%i1+1) /= var%nval(2) .or. &
                 (var%k2-var%k1+1) /= var%nval(3) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes different from file'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            var%rval(var%j1:var%j2,var%i1:var%i2,var%k1:var%k2) = &
              reshape(buffer%realbuff(1:var%totsize),stream%icount(1:3))
          else
            var%rval(1:var%nval(1),1:var%nval(2),1:var%nval(3)) = &
              reshape(buffer%realbuff(1:var%totsize),stream%icount(1:3))
          end if
        class is (ncvariable4d_real)
          ndims = 4
          if ( present(window) ) then
            if ( window(1) < 0 .or. window(2) < 0 .or. &
                 window(1) > window(2) .or. &
                 window(3) < 0 .or. window(4) < 0 .or. &
                 window(3) > window(4) .or. &
                 window(5) < 0 .or. window(6) < 0 .or. &
                 window(5) > window(6) .or. &
                 window(7) < 0 .or. window(8) < 0 .or. &
                 window(7) > window(8) .or. &
                 window(1) > stream%len_dims(var%idims(1)) .or. &
                 window(2) > stream%len_dims(var%idims(1)) .or. &
                 window(3) > stream%len_dims(var%idims(2)) .or. &
                 window(4) > stream%len_dims(var%idims(2)) .or. &
                 window(5) > stream%len_dims(var%idims(3)) .or. &
                 window(6) > stream%len_dims(var%idims(3)) .or. &
                 window(7) > stream%len_dims(var%idims(4)) .or. &
                 window(8) > stream%len_dims(var%idims(4)) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes out of range'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            if ( stream%l_parallel .and. var%lgridded ) then
              if ( stream%jparbound(1) > window(2) .or. &
                   stream%jparbound(2) < window(1) .or. &
                   stream%iparbound(1) > window(4) .or. &
                   stream%iparbound(2) < window(3) ) then
                ! We are not requested to read (out of our window)
                return
              end if
              stream%istart(1) = max(stream%jparbound(1),window(1))
              stream%istart(2) = max(stream%iparbound(1),window(3))
              if ( stream%jparbound(2) > window(2) ) then
                var%nval(1) = window(2)-stream%istart(1)+1
              else
                var%nval(1) = stream%jparbound(2)-stream%istart(1)+1
              end if
              if ( stream%iparbound(2) > window(4) ) then
                var%nval(2) = window(4)-stream%istart(2)+1
              else
                var%nval(2) = stream%iparbound(2)-stream%istart(2)+1
              end if
            else
              stream%istart(1) = window(1)
              stream%istart(2) = window(3)
              var%nval(1) = window(2) - window(1) + 1
              var%nval(2) = window(4) - window(3) + 1
              if ( var%j1 < 0 .and. var%j2 < 0 ) then
                var%j1 = window(1)
                var%j2 = window(2)
              end if
              if ( var%i1 < 0 .and. var%i2 < 0 ) then
                var%i1 = window(3)
                var%i2 = window(4)
              end if
            end if
            stream%istart(3) = window(5)
            var%nval(3) = window(6) - window(5) + 1
            if ( var%k1 < 0 .and. var%k2 < 0 ) then
              var%k1 = window(5)
              var%k2 = window(6)
            end if
            stream%istart(4) = window(7)
            var%nval(4) = window(8) - window(7) + 1
            if ( var%n1 < 0 .and. var%n2 < 0 ) then
              var%n1 = window(7)
              var%n2 = window(8)
            end if
          else
            if ( stream%l_parallel .and. var%lgridded ) then
              stream%istart(1) = stream%jparbound(1)
              stream%istart(2) = stream%iparbound(1)
              var%nval(1) = stream%global_nj
              var%nval(2) = stream%global_ni
            else
              stream%istart(1) = 1
              stream%istart(2) = 1
              var%nval(1) = stream%len_dims(var%idims(1))
              var%nval(2) = stream%len_dims(var%idims(2))
            end if
            stream%istart(3) = 1
            var%nval(3) = stream%len_dims(var%idims(3))
            stream%istart(4) = 1
            var%nval(4) = stream%len_dims(var%idims(4))
          end if
          var%totsize = product(var%nval)
          stream%icount(1) = var%nval(1)
          stream%icount(2) = var%nval(2)
          stream%icount(3) = var%nval(3)
          stream%icount(4) = var%nval(4)
          if ( present(irec) ) then
            if ( var%ndims /= 5 ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested variable is not time dependent'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            stream%istart(5) = irec
            stream%icount(5) = 1
            ndims = 5
          end if
          if ( size(buffer%realbuff) < var%totsize ) then
            deallocate(buffer%realbuff)
            allocate(buffer%realbuff(var%totsize))
          end if
          ncstat = nf90_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
          if ( var%j1 > 0 .and. var%j2 > 0 .and. &
               var%i1 > 0 .and. var%i2 > 0 .and. &
               var%k1 > 0 .and. var%k2 > 0 .and. &
               var%n1 > 0 .and. var%n2 > 0 ) then
            if ( (var%j2-var%j1+1) /= var%nval(1) .or. &
                 (var%i2-var%i1+1) /= var%nval(2) .or. &
                 (var%k2-var%k1+1) /= var%nval(3) .or. &
                 (var%n2-var%n1+1) /= var%nval(4) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes different from file'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            var%rval(var%j1:var%j2,var%i1:var%i2, &
                     var%k1:var%k2,var%n1:var%n2) = &
              reshape(buffer%realbuff(1:var%totsize),stream%icount(1:4))
          else
            var%rval(1:var%nval(1),1:var%nval(2), &
                     1:var%nval(3),1:var%nval(4)) = &
              reshape(buffer%realbuff(1:var%totsize),stream%icount(1:4))
          end if
        class is (ncvariable0d_integer)
          ncstat = nf90_get_var(stream%id,var%id,var%ival)
        class is (ncvariable1d_integer)
          ncstat = nf90_get_var(stream%id,var%id,var%ival)
        class is (ncvariable2d_integer)
          ndims = 2
          if ( present(window) ) then
            if ( window(1) < 0 .or. window(2) < 0 .or. &
                 window(1) > window(2) .or. &
                 window(3) < 0 .or. window(4) < 0 .or. &
                 window(3) > window(4) .or. &
                 window(1) > stream%len_dims(var%idims(1)) .or. &
                 window(2) > stream%len_dims(var%idims(1)) .or. &
                 window(3) > stream%len_dims(var%idims(2)) .or. &
                 window(4) > stream%len_dims(var%idims(2)) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes out of range'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            if ( stream%l_parallel .and. var%lgridded ) then
              if ( stream%jparbound(1) > window(2) .or. &
                   stream%jparbound(2) < window(1) .or. &
                   stream%iparbound(1) > window(4) .or. &
                   stream%iparbound(2) < window(3) ) then
                ! We are not requested to read (out of our window)
                return
              end if
              stream%istart(1) = max(stream%jparbound(1),window(1))
              stream%istart(2) = max(stream%iparbound(1),window(3))
              if ( stream%jparbound(2) > window(2) ) then
                var%nval(1) = window(2)-stream%istart(1)+1
              else
                var%nval(1) = stream%jparbound(2)-stream%istart(1)+1
              end if
              if ( stream%iparbound(2) > window(4) ) then
                var%nval(2) = window(4)-stream%istart(2)+1
              else
                var%nval(2) = stream%iparbound(2)-stream%istart(2)+1
              end if
            else
              stream%istart(1) = window(1)
              stream%istart(2) = window(3)
              var%nval(1) = window(2) - window(1) + 1
              var%nval(2) = window(4) - window(3) + 1
              if ( var%j1 < 0 .and. var%j2 < 0 ) then
                var%j1 = window(1)
                var%j2 = window(2)
              end if
              if ( var%i1 < 0 .and. var%i2 < 0 ) then
                var%i1 = window(3)
                var%i2 = window(4)
              end if
            end if
          else
            if ( stream%l_parallel .and. var%lgridded ) then
              stream%istart(1) = stream%jparbound(1)
              stream%istart(2) = stream%iparbound(1)
              var%nval(1) = stream%global_nj
              var%nval(2) = stream%global_ni
            else
              stream%istart(1) = 1
              stream%istart(2) = 1
              var%nval(1) = stream%len_dims(var%idims(1))
              var%nval(2) = stream%len_dims(var%idims(2))
            end if
          end if
          var%totsize = product(var%nval)
          stream%icount(1) = var%nval(1)
          stream%icount(2) = var%nval(2)
          if ( present(irec) ) then
            if ( var%ndims /= 3 ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested variable is not time dependent'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            stream%istart(3) = irec
            stream%icount(3) = 1
            ndims = 3
          end if
          if ( size(buffer%intbuff) < var%totsize ) then
            deallocate(buffer%intbuff)
            allocate(buffer%intbuff(var%totsize))
          end if
          ncstat = nf90_get_var(stream%id,var%id,buffer%intbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
          if ( var%j1 > 0 .and. var%j2 > 0 .and. &
               var%i1 > 0 .and. var%i2 > 0 ) then
            if ( (var%j2-var%j1+1) /= var%nval(1) .or. &
                 (var%i2-var%i1+1) /= var%nval(2) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes different from file'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            var%ival(var%j1:var%j2,var%i1:var%i2) = &
              reshape(buffer%intbuff(1:var%totsize),stream%icount(1:2))
          else
            var%ival(1:var%nval(1),1:var%nval(2)) = &
              reshape(buffer%intbuff(1:var%totsize),stream%icount(1:2))
          end if
        class is (ncvariable3d_integer)
          ndims = 3
          if ( present(window) ) then
            if ( window(1) < 0 .or. window(2) < 0 .or. &
                 window(1) > window(2) .or. &
                 window(3) < 0 .or. window(4) < 0 .or. &
                 window(3) > window(4) .or. &
                 window(5) < 0 .or. window(6) < 0 .or. &
                 window(5) > window(6) .or. &
                 window(1) > stream%len_dims(var%idims(1)) .or. &
                 window(2) > stream%len_dims(var%idims(1)) .or. &
                 window(3) > stream%len_dims(var%idims(2)) .or. &
                 window(4) > stream%len_dims(var%idims(2)) .or. &
                 window(5) > stream%len_dims(var%idims(3)) .or. &
                 window(6) > stream%len_dims(var%idims(3)) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes out of range'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            if ( stream%l_parallel .and. var%lgridded ) then
              if ( stream%jparbound(1) > window(2) .or. &
                   stream%jparbound(2) < window(1) .or. &
                   stream%iparbound(1) > window(4) .or. &
                   stream%iparbound(2) < window(3) ) then
                ! We are not requested to read (out of our window)
                return
              end if
              stream%istart(1) = max(stream%jparbound(1),window(1))
              stream%istart(2) = max(stream%iparbound(1),window(3))
              if ( stream%jparbound(2) > window(2) ) then
                var%nval(1) = window(2)-stream%istart(1)+1
              else
                var%nval(1) = stream%jparbound(2)-stream%istart(1)+1
              end if
              if ( stream%iparbound(2) > window(4) ) then
                var%nval(2) = window(4)-stream%istart(2)+1
              else
                var%nval(2) = stream%iparbound(2)-stream%istart(2)+1
              end if
            else
              stream%istart(1) = window(1)
              stream%istart(2) = window(3)
              var%nval(1) = window(2) - window(1) + 1
              var%nval(2) = window(4) - window(3) + 1
              if ( var%j1 < 0 .and. var%j2 < 0 ) then
                var%j1 = window(1)
                var%j2 = window(2)
              end if
              if ( var%i1 < 0 .and. var%i2 < 0 ) then
                var%i1 = window(3)
                var%i2 = window(4)
              end if
            end if
            stream%istart(3) = window(5)
            var%nval(3) = window(6) - window(5) + 1
            if ( var%k1 < 0 .and. var%k2 < 0 ) then
              var%k1 = window(5)
              var%k2 = window(6)
            end if
          else
            if ( stream%l_parallel .and. var%lgridded ) then
              stream%istart(1) = stream%jparbound(1)
              stream%istart(2) = stream%iparbound(1)
              var%nval(1) = stream%global_nj
              var%nval(2) = stream%global_ni
            else
              stream%istart(1) = 1
              stream%istart(2) = 1
              var%nval(1) = stream%len_dims(var%idims(1))
              var%nval(2) = stream%len_dims(var%idims(2))
            end if
            stream%istart(3) = 1
            var%nval(3) = stream%len_dims(var%idims(3))
          end if
          var%totsize = product(var%nval)
          stream%icount(1) = var%nval(1)
          stream%icount(2) = var%nval(2)
          stream%icount(3) = var%nval(3)
          if ( present(irec) ) then
            if ( var%ndims /= 4 ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested variable is not time dependent'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            stream%istart(4) = irec
            stream%icount(4) = 1
            ndims = 4
          end if
          if ( size(buffer%intbuff) < var%totsize ) then
            deallocate(buffer%intbuff)
            allocate(buffer%intbuff(var%totsize))
          end if
          ncstat = nf90_get_var(stream%id,var%id,buffer%intbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
          if ( var%j1 > 0 .and. var%j2 > 0 .and. &
               var%i1 > 0 .and. var%i2 > 0 .and. &
               var%k1 > 0 .and. var%k2 > 0 ) then
            if ( (var%j2-var%j1+1) /= var%nval(1) .or. &
                 (var%i2-var%i1+1) /= var%nval(2) .or. &
                 (var%k2-var%k1+1) /= var%nval(3) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes different from file'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            var%ival(var%j1:var%j2,var%i1:var%i2,var%k1:var%k2) = &
              reshape(buffer%intbuff(1:var%totsize),stream%icount(1:3))
          else
            var%ival(1:var%nval(1),1:var%nval(2),1:var%nval(3)) = &
              reshape(buffer%intbuff(1:var%totsize),stream%icount(1:3))
          end if
        class is (ncvariable4d_integer)
          ndims = 4
          if ( present(window) ) then
            if ( window(1) < 0 .or. window(2) < 0 .or. &
                 window(1) > window(2) .or. &
                 window(3) < 0 .or. window(4) < 0 .or. &
                 window(3) > window(4) .or. &
                 window(5) < 0 .or. window(6) < 0 .or. &
                 window(5) > window(6) .or. &
                 window(7) < 0 .or. window(8) < 0 .or. &
                 window(7) > window(8) .or. &
                 window(1) > stream%len_dims(var%idims(1)) .or. &
                 window(2) > stream%len_dims(var%idims(1)) .or. &
                 window(3) > stream%len_dims(var%idims(2)) .or. &
                 window(4) > stream%len_dims(var%idims(2)) .or. &
                 window(5) > stream%len_dims(var%idims(3)) .or. &
                 window(6) > stream%len_dims(var%idims(3)) .or. &
                 window(7) > stream%len_dims(var%idims(4)) .or. &
                 window(8) > stream%len_dims(var%idims(4)) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes out of range'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            if ( stream%l_parallel .and. var%lgridded ) then
              if ( stream%jparbound(1) > window(2) .or. &
                   stream%jparbound(2) < window(1) .or. &
                   stream%iparbound(1) > window(4) .or. &
                   stream%iparbound(2) < window(3) ) then
                ! We are not requested to read (out of our window)
                return
              end if
              stream%istart(1) = max(stream%jparbound(1),window(1))
              stream%istart(2) = max(stream%iparbound(1),window(3))
              if ( stream%jparbound(2) > window(2) ) then
                var%nval(1) = window(2)-stream%istart(1)+1
              else
                var%nval(1) = stream%jparbound(2)-stream%istart(1)+1
              end if
              if ( stream%iparbound(2) > window(4) ) then
                var%nval(2) = window(4)-stream%istart(2)+1
              else
                var%nval(2) = stream%iparbound(2)-stream%istart(2)+1
              end if
            else
              stream%istart(1) = window(1)
              stream%istart(2) = window(3)
              var%nval(1) = window(2) - window(1) + 1
              var%nval(2) = window(4) - window(3) + 1
              if ( var%j1 < 0 .and. var%j2 < 0 ) then
                var%j1 = window(1)
                var%j2 = window(2)
              end if
              if ( var%i1 < 0 .and. var%i2 < 0 ) then
                var%i1 = window(3)
                var%i2 = window(4)
              end if
            end if
            stream%istart(3) = window(5)
            var%nval(3) = window(6) - window(5) + 1
            if ( var%k1 < 0 .and. var%k2 < 0 ) then
              var%k1 = window(5)
              var%k2 = window(6)
            end if
            stream%istart(4) = window(7)
            var%nval(4) = window(8) - window(7) + 1
            if ( var%n1 < 0 .and. var%n2 < 0 ) then
              var%n1 = window(7)
              var%n2 = window(8)
            end if
          else
            if ( stream%l_parallel .and. var%lgridded ) then
              stream%istart(1) = stream%jparbound(1)
              stream%istart(2) = stream%iparbound(1)
              var%nval(1) = stream%global_nj
              var%nval(2) = stream%global_ni
            else
              stream%istart(1) = 1
              stream%istart(2) = 1
              var%nval(1) = stream%len_dims(var%idims(1))
              var%nval(2) = stream%len_dims(var%idims(2))
            end if
            stream%istart(3) = 1
            var%nval(3) = stream%len_dims(var%idims(3))
            stream%istart(4) = 1
            var%nval(4) = stream%len_dims(var%idims(4))
          end if
          var%totsize = product(var%nval)
          stream%icount(1) = var%nval(1)
          stream%icount(2) = var%nval(2)
          stream%icount(3) = var%nval(3)
          stream%icount(4) = var%nval(4)
          if ( present(irec) ) then
            if ( var%ndims /= 5 ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested variable is not time dependent'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            stream%istart(5) = irec
            stream%icount(5) = 1
            ndims = 5
          end if
          if ( size(buffer%intbuff) < var%totsize ) then
            deallocate(buffer%intbuff)
            allocate(buffer%intbuff(var%totsize))
          end if
          ncstat = nf90_get_var(stream%id,var%id,buffer%intbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
          if ( var%j1 > 0 .and. var%j2 > 0 .and. &
               var%i1 > 0 .and. var%i2 > 0 .and. &
               var%k1 > 0 .and. var%k2 > 0 .and. &
               var%n1 > 0 .and. var%n2 > 0 ) then
            if ( (var%j2-var%j1+1) /= var%nval(1) .or. &
                 (var%i2-var%i1+1) /= var%nval(2) .or. &
                 (var%k2-var%k1+1) /= var%nval(3) .or. &
                 (var%n2-var%n1+1) /= var%nval(4) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Requested indexes different from file'
              call die('nc_stream','Cannot read variable '// &
                trim(var%vname)//' in file '//trim(stream%filename), 1)
            end if
            var%ival(var%j1:var%j2,var%i1:var%i2, &
                     var%k1:var%k2,var%n1:var%n2) = &
              reshape(buffer%intbuff(1:var%totsize),stream%icount(1:4))
          else
            var%ival(1:var%nval(1),1:var%nval(2), &
                     1:var%nval(3),1:var%nval(4)) = &
              reshape(buffer%intbuff(1:var%totsize),stream%icount(1:4))
          end if
        class default
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', &
            'Cannot read variable '//trim(var%vname)// &
            ' in file '//trim(stream%filename), 1)
      end select
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) nf90_strerror(ncstat)
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        call die('nc_stream','Cannot read variable '//trim(var%vname)// &
          ' in '//trim(stream%filename),1)
      end if
    end subroutine instream_readvar

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
  type(nc_input_stream) :: ncin
  type(ncoutstream_params) :: opar
  type(ncinstream_params) :: ipar

  type(rcm_time_and_date) :: idate
  real(rk8) :: xrec

  type(ncvariable0d_integer) :: var0dint
  type(ncvariable1d_real) :: var1dreal
  type(ncvariable2d_real) :: var2dreal
  type(ncvariable3d_real) :: var3dreal

  type(ncvariable2d_real) :: var2read

  real(rk8) , target , dimension(18) :: sigma
  real(rk8) , target , dimension(16,12) :: d2dvar
  real(rk8) , target , dimension(16,12,18) :: d3dvar
  integer(ik4) :: i , k

  data sigma /0.025000000372529, 0.0750000011175871, 0.129999998956919, &
              0.195000000298023, 0.270000003278255, 0.349999994039536, &
              0.429999992251396, 0.510000005364418, 0.590000003576279, &
              0.669999986886978, 0.744999974966049, 0.809999972581863, &
              0.864999979734421, 0.909999996423721, 0.944999992847443, &
              0.969999998807907, 0.985000014305115, 0.995000004768372 /
  jx = 16
  iy = 12
  kz = 18
  ds = 50.0D0
  domname = 'domname'
  calendar = 'gregorian'
  ical = 1
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
  opar%zero_date = 1979022200
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
  d2dvar(1:jx,1:iy) = 1.0D0
  d2dvar(jx/2,iy/2) = 2.0D0
  var2dreal%rval => d2dvar

  ! Write some static variables
  call outstream_writevar(ncout,var1dreal)
  call outstream_writevar(ncout,var2dreal)

  var0dint%ival(1) = 12
  d3dvar(1:jx,1:iy,:) = 1.0D0
  do k = 1 , kz
    d3dvar(jx/4,iy/4,k) = dble(k)
  end do
  var3dreal%rval => d3dvar

  ! Write variables in the current record step
  idate = 1979022206
  call outstream_addrec(ncout,idate)
  call outstream_writevar(ncout,var0dint)
  call outstream_writevar(ncout,var3dreal)

  var0dint%ival(1) = 13
  d3dvar(1:jx,1:iy,:) = 1.0D0
  do k = 1 , kz
    d3dvar(jx/4,iy/4,k) = dble(k)*1.5D0
  end do

  ! Add a new record
  idate = 1979022212
  call outstream_addrec(ncout,idate)
  call outstream_writevar(ncout,var0dint)
  call outstream_writevar(ncout,var3dreal)

  ! Finally, close the file and cleanup all
  call outstream_dispose(ncout)

  ipar%fname = opar%fname
  call instream_setup(ncin,ipar)

  idate = 1979022207

  call instream_findrec(ncin,idate,xrec)

  print *, 'Record ',trim(tochar(idate)),' is at ', xrec

  d2dvar(:,:) = -1.0D0
  var2read%vname = 'real2dvar'
  var2read%lgridded = .true.
  var2read%rval => d2dvar
  var2read%j1 = 1
  var2read%j2 = 9
  var2read%i1 = 1
  var2read%i2 = 7

  call instream_readvar(ncin,var2read,window=(/4,12,2,8/))

  do i = 1 , iy
    print '(16i2)', int(d2dvar(:,i))
  end do

  call instream_dispose(ncin)

end program test
#endif

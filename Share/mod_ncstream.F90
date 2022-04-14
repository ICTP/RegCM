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
  use mod_zita
#ifdef PNETCDF
  use mpi , only : mpi_comm_self , mpi_info_null
  use pnetcdf
#else
  use netcdf
#endif

  implicit none

  private

  integer(ik4) :: ncstat
  logical , parameter :: nocopy = .false.
  type(rcm_time_and_date) , save :: reference_date
  integer(ik4) , dimension(ncmaxdims) :: id_dim
  integer(ik4) , dimension(ncmaxdims) :: len_dim

  integer(ik4) , parameter :: cordexvtype = nf90_float

  interface outstream_addrec
    module procedure outstream_addrec_date
    module procedure outstream_addrec_value
  end interface

  public :: ncoutstream_params
  public :: nc_output_stream
  public :: ncvariable0d_char
  public :: ncvariable_standard
  public :: ncvariable0d_real
  public :: ncvariable0d_double
  public :: ncvariable0d_mixed
  public :: ncvariable0d_integer
  public :: ncvariable1d_real
  public :: ncvariable1d_double
  public :: ncvariable1d_mixed
  public :: ncvariable1d_integer
  public :: ncvariable2d_real
  public :: ncvariable2d_double
  public :: ncvariable2d_mixed
  public :: ncvariable2d_integer
  public :: ncvariable3d_real
  public :: ncvariable3d_double
  public :: ncvariable3d_mixed
  public :: ncvariable3d_integer
  public :: ncvariable4d_real
  public :: ncvariable4d_double
  public :: ncvariable4d_mixed
  public :: ncvariable4d_integer

  public :: ncattribute_string
  public :: ncattribute_logical
  public :: ncattribute_integer
  public :: ncattribute_real4
  public :: ncattribute_real8
  public :: ncattribute_real4_array
  public :: ncattribute_real8_array

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
      integer(ik4) :: imode = nf90_nowrite
      integer(ik4) :: dimtime , i
      type(rcm_time_and_date) :: tt

      if ( associated(ncin%ncp%xs) ) call instream_dispose(ncin)
      allocate(ncin%ncp%xs)
      stream => ncin%ncp%xs
      stream%filename = params%fname
#ifdef NETCDF4_HDF5
      if ( params%mpi_comm /= -1 ) then
        imode = ior(nf90_nowrite,nf90_share)
        ncstat = nf90_open(stream%filename,imode, &
          stream%id,comm=params%mpi_comm,info=params%mpi_info)
        stream%l_parallel = .true.
      else
        ncstat = nf90_open(stream%filename,imode,stream%id)
      end if
#else
      if ( params%mpi_comm /= -1 ) then
#ifdef PNETCDF
        imode = nf90_nowrite
        ncstat = nf90mpi_open(params%mpi_comm, stream%filename, &
                              imode, params%mpi_info, stream%id)
        stream%l_parallel = .true.
#else
#ifdef PNETCDF_IN_NETCDF
        imode = ior(nf90_nowrite,ior(nf90_share,nf90_pnetcdf))
        ncstat = nf90_open(stream%filename,imode, &
          stream%id,comm=params%mpi_comm,info=params%mpi_info)
        stream%l_parallel = .true.
#else
        ncstat = nf90_open(stream%filename,imode,stream%id)
#endif
#endif
      else
#ifdef PNETCDF
        imode = nf90_nowrite
        ncstat = nf90mpi_open(mpi_comm_self,stream%filename, &
                              imode, mpi_info_null, stream%id)
#else
        ncstat = nf90_open(stream%filename,imode,stream%id)
#endif
      end if
#endif
      if ( ncstat /= nf90_noerr ) then
        call printerror
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
        if ( myid == 0 ) then
          write(stdout,*) 'Parallel I/O enabled.'
        end if
#endif
      end if
#ifdef PNETCDF
      ncstat = nf90mpi_inq_dimid(stream%id,'time',dimtime)
#else
      ncstat = nf90_inq_dimid(stream%id,'time',dimtime)
#endif
      if ( ncstat == nf90_noerr ) then
#ifdef PNETCDF
        ncstat = nf90mpi_inquire_dimension(stream%id,dimtime,len=stream%nrec)
#else
        ncstat = nf90_inquire_dimension(stream%id,dimtime,len=stream%nrec)
#endif
        if ( ncstat /= nf90_noerr ) then
          call printerror
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream','Error reading time dimension in '// &
            trim(stream%filename),1)
        end if
        if ( stream%nrec > 0 ) then
#ifdef PNETCDF
          ncstat = nf90mpi_inq_varid(stream%id,'time',stream%timeid)
#else
          ncstat = nf90_inq_varid(stream%id,'time',stream%timeid)
#endif
          if ( ncstat /= nf90_noerr ) then
            call printerror
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            call die('nc_stream','Error reading time variable in '// &
              trim(stream%filename),1)
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_get_att(stream%id,stream%timeid,'units',stream%tunit)
#else
          ncstat = nf90_get_att(stream%id,stream%timeid,'units',stream%tunit)
#endif
          if ( ncstat /= nf90_noerr ) then
            call printerror
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            write(stderr,*) 'Assuming hours since 1949-12-01 00:00:00 UTC'
            write(stderr,*) 'for file ',trim(stream%filename)
            stream%tunit = 'hours since 1949-12-01 00:00:00 UTC'
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_get_att(stream%id,stream%timeid, &
                                   'calendar',stream%tcal)
#else
          ncstat = nf90_get_att(stream%id,stream%timeid,'calendar',stream%tcal)
#endif
          if ( ncstat /= nf90_noerr ) then
            call printerror
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            write(stderr,*) 'Assuming gregorian calendar'
            write(stderr,*) 'for file ',trim(stream%filename)
            stream%tcal = 'gregorian'
          end if
          stream%refdate = timeval2date(0.0_rkx,stream%tunit,stream%tcal)
          stream%istart(1) = 1
          stream%icount(1) = stream%nrec
          stream%istride(1) = stream%nrec-1
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,stream%timeid,stream%xtime, &
            stream%istart(1:1),stream%icount(1:1),stream%istride(1:1))
#else
          ncstat = nf90_get_var(stream%id,stream%timeid,stream%xtime, &
            stream%istart(1:1),stream%icount(1:1),stream%istride(1:1))
#endif
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
#ifdef PNETCDF
      ncstat = nf90mpi_inquire(stream%id,nDimensions=stream%ndims)
#else
      ncstat = nf90_inquire(stream%id,nDimensions=stream%ndims)
#endif
      if ( ncstat /= nf90_noerr ) then
        call printerror
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        write(stderr,*) 'Cannot get dimensional infos for file '// &
          trim(stream%filename)
        call die('nc_stream','Cannot setup file '//trim(stream%filename),1)
      end if
      allocate(stream%len_dims(stream%ndims))
      do i = 1 , stream%ndims
#ifdef PNETCDF
        ncstat = nf90mpi_inquire_dimension(stream%id,i,len=stream%len_dims(i))
#else
        ncstat = nf90_inquire_dimension(stream%id,i,len=stream%len_dims(i))
#endif
        if ( ncstat /= nf90_noerr ) then
          call printerror
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
      type(ncoutstream) , pointer :: stream
      type(rcm_time_and_date) :: tt
#if defined(NETCDF4_HDF5) || defined(PNETCDF_IN_NETCDF) || defined(PNETCDF)
      integer(ik4) :: imode
#endif

      if ( associated(ncout%ncp%xs) ) call outstream_dispose(ncout)
      ! Allocate all space
      allocate(ncout%ncp%xs)
      allocate(ncout%obp%xb)
      allocate(ncout%svp%xv)
      stream => ncout%ncp%xs
      stream%filename = params%fname
      if ( params%l_keep ) then
        stream%l_keep = params%l_keep
#ifdef NETCDF4_HDF5
        if ( params%mpi_comm /= -1 ) then
          imode = ior(params%mpi_iotype,nf90_write)
          ncstat = nf90_open(stream%filename,imode, &
             stream%id,comm=params%mpi_comm,info=params%mpi_info)
          stream%l_parallel = .true.
        else
          ncstat = nf90_open(stream%filename,nf90_write,stream%id)
        end if
#else
        if ( params%mpi_comm /= -1 ) then
#ifdef PNETCDF
          imode = nf90_write
          ncstat = nf90mpi_open(params%mpi_comm,stream%filename,imode, &
                                params%mpi_info,stream%id)
          stream%l_parallel = .true.
#else
#ifdef PNETCDF_IN_NETCDF
          imode = ior(params%mpi_iotype,nf90_write)
          ncstat = nf90_open_par(stream%filename,imode, &
                          params%mpi_comm,params%mpi_info,stream%id)
          stream%l_parallel = .true.
#else
          ncstat = nf90_open(stream%filename,nf90_write,stream%id)
#endif
#endif
        else
#ifdef PNETCDF
          ncstat = nf90mpi_open(mpi_comm_self,stream%filename, &
                                nf90_write,mpi_info_null,stream%id)
#else
          ncstat = nf90_open(stream%filename,nf90_write,stream%id)
#endif
        end if
#endif
      else
#ifdef NETCDF4_HDF5
        if ( params%mpi_comm /= -1 ) then
          imode = ior(params%mpi_iotype,iomode)
          if ( params%l_sync ) imode = ior(imode,nf90_share)
          ncstat = nf90_create(stream%filename,imode, &
                    comm=params%mpi_comm,info=params%mpi_info,ncid=stream%id)
          stream%l_parallel = .true.
        else
          ncstat = nf90_create(stream%filename,iomode,stream%id)
        end if
#else
        if ( params%mpi_comm /= -1 ) then
#ifdef PNETCDF
          imode = ior(nf90_clobber, nf90_64bit_offset)
          if ( params%l_sync ) imode = ior(imode,nf90_share)
          ncstat = nf90mpi_create(params%mpi_comm,stream%filename, &
                                  imode,params%mpi_info,stream%id)
          stream%l_parallel = .true.
#else
#ifdef PNETCDF_IN_NETCDF
          imode = ior(params%mpi_iotype,iomode)
          if ( params%l_sync ) imode = ior(imode,nf90_share)
          ncstat = nf90_create_par(stream%filename,imode, &
                    params%mpi_comm,params%mpi_info,stream%id)
          stream%l_parallel = .true.
#else
          ncstat = nf90_create(stream%filename,iomode,stream%id)
#endif
#endif
        else
#ifdef PNETCDF
          imode = ior(nf90_clobber, nf90_64bit_offset)
          if ( params%l_sync ) imode = ior(imode,nf90_share)
          ncstat = nf90mpi_create(mpi_comm_self,stream%filename, &
                                  imode,mpi_info_null,stream%id)
#else
          ncstat = nf90_create(stream%filename,iomode,stream%id)
#endif
        end if
#endif
      end if
      if ( ncstat /= nf90_noerr ) then
        call printerror
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        call die('nc_stream', &
                 'Cannot create or open file '//trim(stream%filename),1)
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
        if ( myid == 0 ) then
          write(stdout,*) 'Parallel I/O enabled.'
        end if
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
      stream%l_crm         = params%l_crm
      stream%l_sync        = params%l_sync
      stream%l_subgrid     = params%l_subgrid
      stream%l_full_sigma  = params%l_full_sigma
      stream%l_plev        = params%l_plev
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
      if ( stream%id >= 0 ) then
#ifdef PNETCDF
        ncstat = nf90mpi_close(stream%id)
#else
        ncstat = nf90_close(stream%id)
#endif
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call printerror
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
      if ( stream%id >= 0 ) then
#ifdef PNETCDF
        ncstat = nf90mpi_close(stream%id)
#else
        ncstat = nf90_close(stream%id)
#endif
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call printerror
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
      if ( allocated(xbf%doublebuff) ) deallocate(xbf%doublebuff)
    end subroutine deallocate_obuffer

    subroutine deallocate_ibuffer(xbf)
      implicit none
      type(internal_ibuffer) , pointer :: xbf
      if ( .not. associated(xbf) ) return
      if ( allocated(xbf%intbuff) )  deallocate(xbf%intbuff)
      if ( allocated(xbf%realbuff) ) deallocate(xbf%realbuff)
      if ( allocated(xbf%doublebuff) ) deallocate(xbf%doublebuff)
    end subroutine deallocate_ibuffer

    subroutine outstream_enable(ncout,sigma)
      implicit none
      type(nc_output_stream) , intent(inout) :: ncout
      real(rkx) , dimension(:) , pointer , intent(in) :: sigma
      real(rkx) , dimension(size(sigma)) :: zita

      type(ncoutstream) , pointer :: stream
      type(internal_obuffer) , pointer :: buffer
      type(basic_variables) , pointer :: stvar
      integer(ik4) :: maxnum_int , maxnum_real , maxnum_double , i
#ifdef CLM45
      integer(ik4) :: nl
#endif
      character(len=16) , dimension(8) :: tempstr
      real(rkx) :: xds , x0
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
          stvar%tbound_var%vunit = ''
          stvar%tbound_var%axis = 'b'
          stvar%tbound_var%long_name = ''
          stvar%tbound_var%standard_name = ''
          stvar%tbound_var%lrecords = .true.
          call outstream_addvar(ncout,stvar%tbound_var)
        end if
      end if
      if ( stream%l_hasgrid ) then
        x0 = -ds*1000.0_rkx/2.0_rkx ! Cross grid
        stvar%map_var%vname = 'crs'
        stvar%map_var%vunit = ''
        stvar%map_var%long_name = ''
        stvar%map_var%standard_name = ''
        call outstream_addvar(ncout,stvar%map_var)
        select case (iproj)
          case('LAMCON')
            attc%aname = 'proj4_params'
            write(tempstr(1),'(f7.2)') truelatl
            write(tempstr(2),'(f7.2)') truelath
            write(tempstr(3),'(f7.2)') clat
            write(tempstr(4),'(f7.2)') clon
            write(tempstr(5),'(f10.0)') x0
            write(tempstr(6),'(f10.0)') x0
            write(tempstr(7),'(f9.0)') earthrad
            write(tempstr(8),'(f9.0)') earthrad
            attc%theval = '+proj=lcc +lat_1='//trim(adjustl(tempstr(1)))// &
              ' +lat_2='//trim(adjustl(tempstr(2)))// &
              ' +lat_0='//trim(adjustl(tempstr(3)))// &
              ' +lon_0='//trim(adjustl(tempstr(4)))// &
              ' +x_0='//trim(adjustl(tempstr(5)))// &
              ' +y_0='//trim(adjustl(tempstr(6)))// &
              ' +ellps=sphere +a='//trim(adjustl(tempstr(7)))// &
              ' +b='//trim(adjustl(tempstr(8)))//' +units=m +no_defs'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
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
            attr%aname = 'false_easting'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'false_northing'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attc%aname = 'crs_wkt'
            attc%theval = 'PROJCS["unnamed",'//NEW_LINE('A')//                 &
                '    GEOGCS["Normal Sphere (r='//                              &
                trim(adjustl(tempstr(7)))//')",'//NEW_LINE('A')//              &
                '        DATUM["unknown",'//NEW_LINE('A')//                    &
                '            SPHEROID["sphere",'//                             &
                trim(adjustl(tempstr(7)))//',0]],'//NEW_LINE('A')//            &
                '        PRIMEM["Greenwich",0],'//NEW_LINE('A')//              &
                '        UNIT["degree",0.0174532925199433]],'//NEW_LINE('A')// &
                '    PROJECTION["Lambert_Conformal_Conic_2SP"],'//             &
                NEW_LINE('A')//                                                &
                '    PARAMETER["standard_parallel_1",'//                       &
                trim(adjustl(tempstr(1)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["standard_parallel_2",'//                       &
                trim(adjustl(tempstr(2)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["latitude_of_origin",'//                        &
                trim(adjustl(tempstr(3)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["central_meridian",'//                          &
                trim(adjustl(tempstr(4)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["false_easting",'//                             &
                trim(adjustl(tempstr(5)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["false_northing",'//                            &
                trim(adjustl(tempstr(6)))//'],'//NEW_LINE('A')//               &
                '    UNIT["Meter",1]]'//NEW_LINE('A')
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
          case('POLSTR')
            attc%aname = 'proj4_params'
            write(tempstr(1),'(f7.2)') clat
            write(tempstr(2),'(f7.2)') clon
            write(tempstr(3),'(f10.0)') x0
            write(tempstr(4),'(f10.0)') x0
            write(tempstr(5),'(f9.0)') earthrad
            write(tempstr(6),'(f9.0)') earthrad
            attc%theval = '+proj=stere +lat_0='//trim(adjustl(tempstr(1)))// &
              ' +lon_0='//trim(adjustl(tempstr(2)))// &
              ' +x_0='//trim(adjustl(tempstr(3)))// &
              ' +y_0='//trim(adjustl(tempstr(4)))// &
              ' +ellps=sphere +a='//trim(adjustl(tempstr(5)))// &
              ' +b='//trim(adjustl(tempstr(6)))//' +units=m +no_defs'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attc%aname = 'grid_mapping_name'
            attc%theval = 'stereographic'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'latitude_of_projection_origin'
            attr%theval = clat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'longitude_of_projection_origin'
            attr%theval = clon
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'scale_factor_at_projection_origin'
            attr%theval = 1.0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'false_easting'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'false_northing'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attc%aname = 'crs_wkt'
            attc%theval = 'PROJCS["unnamed",'//NEW_LINE('A')//                 &
                '    GEOGCS["Normal Sphere (r='//                              &
                trim(adjustl(tempstr(5)))//')",'//NEW_LINE('A')//              &
                '        DATUM["unknown",'//NEW_LINE('A')//                    &
                '            SPHEROID["sphere",'//                             &
                trim(adjustl(tempstr(5)))//',0]],'//NEW_LINE('A')//            &
                '        PRIMEM["Greenwich",0],'//NEW_LINE('A')//              &
                '        UNIT["degree",0.0174532925199433]],'//NEW_LINE('A')// &
                '    PROJECTION["Stereographic"],'//                           &
                NEW_LINE('A')//                                                &
                '    PARAMETER["latitude_of_origin",'//                        &
                trim(adjustl(tempstr(1)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["central_meridian",'//                          &
                trim(adjustl(tempstr(2)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["scale_factor",1.],'//NEW_LINE('A')//           &
                '    PARAMETER["false_easting",'//                             &
                trim(adjustl(tempstr(3)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["false_northing",'//                            &
                trim(adjustl(tempstr(4)))//'],'//NEW_LINE('A')//               &
                '    UNIT["Meter",1]]'//NEW_LINE('A')
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
          case('ROTLLR')
            attc%aname = 'proj4_params'
            write(tempstr(1),'(f7.2)') clat
            write(tempstr(2),'(f7.2)') clon
            write(tempstr(3),'(f10.0)') x0
            write(tempstr(4),'(f10.0)') x0
            write(tempstr(5),'(f9.0)') earthrad
            write(tempstr(6),'(f7.2)') plat
            write(tempstr(7),'(f7.2)') plon
            attc%theval = '+proj=ob_tran +o_proj=longlat'// &
              ' +o_lat_p='//trim(adjustl(tempstr(6)))// &
              ' +o_lon_p='//trim(adjustl(tempstr(7)))// &
              ' +R='//trim(adjustl(tempstr(5)))// &
              ' +lon_0=180.0 +to_meter=0.01745329'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attc%aname = 'grid_mapping_name'
            attc%theval = 'rotated_latitude_longitude'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'latitude_of_projection_origin'
            attr%theval = clat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'longitude_of_projection_origin'
            attr%theval = clon
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'grid_north_pole_latitude'
            attr%theval = plat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'grid_north_pole_longitude'
            attr%theval = plon
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'false_easting'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'false_northing'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attc%aname = 'crs_wkt'
            attc%theval = 'PROJCS["unnamed",'//NEW_LINE('A')//                 &
                '    GEOGCS["Normal Sphere (r='//                              &
                trim(adjustl(tempstr(5)))//')",'//NEW_LINE('A')//              &
                '        DATUM["unknown",'//NEW_LINE('A')//                    &
                '            SPHEROID["sphere",'//                             &
                trim(adjustl(tempstr(5)))//',0]],'//NEW_LINE('A')//            &
                '        PRIMEM["Greenwich",0],'//NEW_LINE('A')//              &
                '        UNIT["degree",0.0174532925199433]],'//NEW_LINE('A')// &
                '    PROJECTION["Rotated_Latitude_Longitude"],'//              &
                NEW_LINE('A')//                                                &
                '    PARAMETER["grid_north_pole_latitude",'//                  &
                trim(adjustl(tempstr(6)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["grid_north_pole_longitude",'//                 &
                trim(adjustl(tempstr(7)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["latitude_of_origin",'//                        &
                trim(adjustl(tempstr(1)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["central_meridian",'//                          &
                trim(adjustl(tempstr(2)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["false_easting",'//                             &
                trim(adjustl(tempstr(3)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["false_northing",'//                            &
                trim(adjustl(tempstr(4)))//'],'//NEW_LINE('A')//               &
                '    UNIT["Meter",1]]'//NEW_LINE('A')
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
          case('NORMER')
            attc%aname = 'proj4_params'
            write(tempstr(1),'(f7.2)') clat
            write(tempstr(2),'(f7.2)') clon
            write(tempstr(3),'(f10.0)') x0
            write(tempstr(4),'(f10.0)') x0
            write(tempstr(5),'(f9.0)') earthrad
            write(tempstr(6),'(f9.0)') earthrad
            attc%theval = '+proj=merc +lat_ts='//trim(adjustl(tempstr(1)))// &
              ' +lon_0='//trim(adjustl(tempstr(2)))// &
              ' +x_0='//trim(adjustl(tempstr(3)))// &
              ' +y_0='//trim(adjustl(tempstr(4)))// &
              ' +ellps=sphere +a='//trim(adjustl(tempstr(5)))// &
              ' +b='//trim(adjustl(tempstr(6)))//' +units=m +no_defs'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attc%aname = 'grid_mapping_name'
            attc%theval = 'mercator'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'standard_parallel'
            attr%theval = clat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'longitude_of_projection_origin'
            attr%theval = clon
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'false_easting'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'false_northing'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attc%aname = 'crs_wkt'
            attc%theval = 'PROJCS["unnamed",'//NEW_LINE('A')//                 &
                '    GEOGCS["Normal Sphere (r='//                              &
                trim(adjustl(tempstr(5)))//')",'//NEW_LINE('A')//              &
                '        DATUM["unknown",'//NEW_LINE('A')//                    &
                '            SPHEROID["sphere",'//                             &
                trim(adjustl(tempstr(5)))//',0]],'//NEW_LINE('A')//            &
                '        PRIMEM["Greenwich",0],'//NEW_LINE('A')//              &
                '        UNIT["degree",0.0174532925199433]],'//NEW_LINE('A')// &
                '    PROJECTION["Mercator_2SP"],'//                            &
                NEW_LINE('A')//                                                &
                '    PARAMETER["standard_parallel_1",'//                       &
                trim(adjustl(tempstr(1)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["central_meridian",'//                          &
                trim(adjustl(tempstr(2)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["false_easting",'//                             &
                trim(adjustl(tempstr(3)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["false_northing",'//                            &
                trim(adjustl(tempstr(4)))//'],'//NEW_LINE('A')//               &
                '    UNIT["Meter",1]]'//NEW_LINE('A')
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
          case('ROTMER')
            attc%aname = 'proj4_params'
            write(tempstr(1),'(f7.2)') plat
            write(tempstr(2),'(f7.2)') plon
            write(tempstr(3),'(f10.0)') x0
            write(tempstr(4),'(f10.0)') x0
            write(tempstr(5),'(f9.0)') earthrad
            write(tempstr(6),'(f9.0)') earthrad
            attc%theval = '+proj=omerc +lat_0='//trim(adjustl(tempstr(1)))// &
              ' +alpha=89.999999 +lonc='//trim(adjustl(tempstr(2)))// &
              ' +x_0='//trim(adjustl(tempstr(3)))// &
              ' +y_0='//trim(adjustl(tempstr(4)))// &
              ' +ellps=sphere +a='//trim(adjustl(tempstr(5)))// &
              ' +b='//trim(adjustl(tempstr(6)))//' +units=m +no_defs'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attc%aname = 'grid_mapping_name'
            attc%theval = 'oblique_mercator'
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'latitude_of_projection_origin'
            attr%theval = plat
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'longitude_of_projection_origin'
            attr%theval = plon
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'scale_factor_at_projection_origin'
            attr%theval = 1.0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'azimuth_of_central_line'
            attr%theval = 89.999999
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'false_easting'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attr%aname = 'false_northing'
            attr%theval = x0
            call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
            attc%aname = 'crs_wkt'
            attc%theval = 'PROJCS["unnamed",'//NEW_LINE('A')//                 &
                '    GEOGCS["Normal Sphere (r='//                              &
                trim(adjustl(tempstr(5)))//')",'//NEW_LINE('A')//              &
                '        DATUM["unknown",'//NEW_LINE('A')//                    &
                '            SPHEROID["sphere",'//                             &
                trim(adjustl(tempstr(5)))//',0]],'//NEW_LINE('A')//            &
                '        PRIMEM["Greenwich",0],'//NEW_LINE('A')//              &
                '        UNIT["degree",0.0174532925199433]],'//NEW_LINE('A')// &
                '    PROJECTION["Hotine_Oblique_Mercator_Azimuth_Center"],'//  &
                NEW_LINE('A')//                                                &
                '    PARAMETER["latitude_of_center",'//                        &
                trim(adjustl(tempstr(1)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["longitude_of_center",'//                       &
                trim(adjustl(tempstr(2)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["azimuth",89.999999],'//NEW_LINE('A')//         &
                '    PARAMETER["rectified_grid_angle",89.999999],'//           &
                NEW_LINE('A')//                                                &
                '    PARAMETER["scale_factor",1.],'//NEW_LINE('A')//           &
                '    PARAMETER["false_easting",'//                             &
                trim(adjustl(tempstr(3)))//'],'//NEW_LINE('A')//               &
                '    PARAMETER["false_northing",'//                            &
                trim(adjustl(tempstr(4)))//'],'//NEW_LINE('A')//               &
                '    UNIT["Meter",1]]'//NEW_LINE('A')
            call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
        end select
        attr%aname = 'semi_major_axis'
        attr%theval = earthrad
        call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
        attr%aname = 'inverse_flattening'
        attr%theval = 0.0
        call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
        attr%aname = 'false_easting'
        attr%theval = x0
        call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
        attr%aname = 'false_northing'
        attr%theval = x0
        call add_attribute(stream,attr,stvar%map_var%id,stvar%map_var%vname)
        attc%aname = 'CoordinateTransformType'
        attc%theval = 'Projection'
        call add_attribute(stream,attc,stvar%map_var%id,stvar%map_var%vname)
        attc%aname = 'CoordinateAxisTypes'
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
      if ( stream%l_has100mlev ) then
        stvar%lev100m_var%vname = 'm100'
        stvar%lev100m_var%vunit = 'm'
        stvar%lev100m_var%axis = 'w'
        stvar%lev100m_var%long_name = 'Height level'
        stvar%lev100m_var%standard_name = 'height'
        stvar%lev100m_var%lrecords = .false.
        call outstream_addvar(ncout,stvar%lev100m_var)
      end if
      if ( stream%l_hassoillev ) then
        stvar%levsoil_var%vname = 'soil_layer'
        stvar%levsoil_var%vunit = 'm'
        stvar%levsoil_var%axis = 's'
        stvar%levsoil_var%long_name = 'Soil layer level'
        stvar%levsoil_var%standard_name = 'depth'
        stvar%levsoil_var%lrecords = .false.
        call outstream_addvar(ncout,stvar%levsoil_var)
        attc%aname = 'bounds'
        attc%theval = 'soil_bounds'
        call add_attribute(stream,attc,stvar%levsoil_var%id, &
                           stvar%levsoil_var%vname)
        stvar%levsoilbound_var%vname = 'soil_bounds'
        stvar%levsoilbound_var%vunit = 'm'
        stvar%levsoilbound_var%axis = 'bs'
        stvar%levsoilbound_var%long_name = 'Soil layer level bounds'
        stvar%levsoilbound_var%standard_name = 'depth'
        stvar%levsoilbound_var%lrecords = .false.
        call outstream_addvar(ncout,stvar%levsoilbound_var)
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
      if ( .not. stream%l_keep ) then
#ifdef PNETCDF
        ncstat = nf90mpi_enddef(stream%id)
#else
        ncstat = nf90_enddef(stream%id)
#endif
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call printerror
          call die('nc_stream','Cannot enable file '//trim(stream%filename),1)
        end if
      end if
      !
      ! Allocate buffer space shared by all vars
      !
      maxnum_int  = buffer%max1d_int(1)
      maxnum_real = buffer%max1d_real(1)
      maxnum_double = buffer%max1d_double(1)

      if ( buffer%lhas2dint ) then
        maxnum_int = max(product(buffer%max2d_int),maxnum_int)
      end if
      if ( buffer%lhas2dreal ) then
        maxnum_real = max(product(buffer%max2d_real),maxnum_real)
      end if
      if ( buffer%lhas2ddouble ) then
        maxnum_double = max(product(buffer%max2d_double),maxnum_double)
      end if
      if ( buffer%lhas3dint ) then
        maxnum_int = max(product(buffer%max3d_int),maxnum_int)
      end if
      if ( buffer%lhas3dreal ) then
        maxnum_real = max(product(buffer%max3d_real),maxnum_real)
      end if
      if ( buffer%lhas3ddouble ) then
        maxnum_double = max(product(buffer%max3d_double),maxnum_double)
      end if
      if ( buffer%lhas4dint ) then
        maxnum_int = max(product(buffer%max4d_int),maxnum_int)
      end if
      if ( buffer%lhas4dreal ) then
        maxnum_real = max(product(buffer%max4d_real),maxnum_real)
      end if
      if ( buffer%lhas4ddouble ) then
        maxnum_double = max(product(buffer%max4d_double),maxnum_double)
      end if
      if ( maxnum_int > 0 ) allocate(buffer%intbuff(maxnum_int))
      if ( maxnum_real > 0 ) allocate(buffer%realbuff(maxnum_real))
      if ( maxnum_double > 0 ) allocate(buffer%doublebuff(maxnum_double))
      stream%l_enabled = .true.
#ifdef DEBUG
      if ( myid == 0 ) then
        write(stdout,*) 'Enabled netCDF output stream ',trim(stream%filename)
        if ( allocated(buffer%intbuff) ) &
          write(stdout,*) 'Total buffer integer size :', &
            size(buffer%intbuff)*4
        if ( allocated(buffer%realbuff) ) &
          write(stdout,*) 'Total buffer float size   :', &
            size(buffer%realbuff)*4
          if ( allocated(buffer%doublebuff) ) &
            write(stdout,*) 'Total buffer double size   :', &
            size(buffer%doublebuff)*4
        end if
#endif
      ! Put "basic" information in the file
      if ( stream%l_subgrid ) then
        xds = ds/real(nsg,rkx)
      else
        xds = ds
      end if
      if ( iproj /= 'ROTLLR' ) then
        buffer%doublebuff(1) = &
          -real(((real(stream%len_dims(jx_dim),rkx)-d_one)/d_two) * &
                    xds * d_1000,rk8)
        do i = 2 , stream%len_dims(jx_dim)
          buffer%doublebuff(i) = &
            real(real(buffer%doublebuff(i-1),rkx)+xds*d_1000,rk8)
        end do
        call outstream_writevar(ncout,stvar%jx_var,nocopy)
        buffer%doublebuff(1) = &
          -real(((real(stream%len_dims(iy_dim),rkx)-d_one)/d_two) * &
                    xds * d_1000,rk8)
        do i = 2 , stream%len_dims(iy_dim)
          buffer%doublebuff(i) = &
            real(real(buffer%doublebuff(i-1),rkx)+xds*d_1000,rk8)
        end do
        call outstream_writevar(ncout,stvar%iy_var,nocopy)
      else
        xds = xds/erkm*raddeg
        buffer%doublebuff(1) = &
          -real(((real(stream%len_dims(jx_dim),rkx)-d_one)/d_two) * &
                    xds,rk8)
        do i = 2 , stream%len_dims(jx_dim)
          buffer%doublebuff(i) = &
            real(real(buffer%doublebuff(i-1),rkx)+xds,rk8)
        end do
        where ( buffer%doublebuff > 180.0 )
          buffer%doublebuff = 360.0_rk8 - buffer%doublebuff
        end where
        where ( buffer%doublebuff < -180.0 )
          buffer%doublebuff = 360.0_rk8 + buffer%doublebuff
        end where
        call outstream_writevar(ncout,stvar%jx_var,nocopy)
        buffer%doublebuff(1) = &
          -real(((real(stream%len_dims(iy_dim),rkx)-d_one)/d_two) * &
                    xds,rk8)
        do i = 2 , stream%len_dims(iy_dim)
          buffer%doublebuff(i) = &
            real(real(buffer%doublebuff(i-1),rkx)+xds,rk8)
        end do
        call outstream_writevar(ncout,stvar%iy_var,nocopy)
      end if
      buffer%doublebuff(1:size(sigma)) = real(sigma,rk8)
      call outstream_writevar(ncout,stvar%sigma_var,nocopy)
      if ( .not. stream%l_plev ) then
        if ( idynamic < 3 ) then
          stvar%ptop_var%rval(1) = real(ptop*10.0_rkx,rk8)
          call outstream_writevar(ncout,stvar%ptop_var)
        else
          zita = zitasigma(sigma)
          buffer%doublebuff(1:size(sigma)) = md_ak(zita)
          call outstream_writevar(ncout,stvar%ak_var,nocopy)
          buffer%doublebuff(1:size(sigma)) = md_bk(zita)
          call outstream_writevar(ncout,stvar%bk_var,nocopy)
        end if
      end if
      if ( stream%l_has2mlev ) then
        buffer%doublebuff(1) = 2.0_rk8
        call outstream_writevar(ncout,stvar%lev2m_var,nocopy)
      end if
      if ( stream%l_has10mlev ) then
        buffer%doublebuff(1) = 10.0_rk8
        call outstream_writevar(ncout,stvar%lev10m_var,nocopy)
      end if
      if ( stream%l_has100mlev ) then
        buffer%doublebuff(1) = 100.0_rk8
        call outstream_writevar(ncout,stvar%lev100m_var,nocopy)
      end if
      if ( stream%l_hassoillev ) then
#ifdef CLM45
        do nl = 1 , num_soil_layers
          buffer%doublebuff(nl) = real(scalez*(exp(0.5_rkx * &
            (real(nl,rkx)-0.5_rkx))-1._rkx),rk8)
        end do
#else
        ! Here is not precise, as the depth of levels is function of the
        ! landuse class.
        buffer%doublebuff(1) = 0.05_rk8
        buffer%doublebuff(2) = 0.95_rk8
        buffer%doublebuff(3) = 3.00_rk8
#endif
        call outstream_writevar(ncout,stvar%levsoil_var,nocopy)
#ifdef CLM45
        buffer%doublebuff(1) = 0.0_rk8
        do nl = 2 , num_soil_layers
          buffer%doublebuff(2*nl-1) = 0.5_rk8 * &
            (real(scalez*(exp(0.5_rkx * &
                         (real(nl-1,rkx)-0.5_rkx))-1.0_rkx),rk8) + &
             real(scalez*(exp(0.5_rkx * &
                         (real(nl,rkx)-0.5_rkx))-1.0_rkx),rk8))
        end do
        do nl = 1 , num_soil_layers
          buffer%doublebuff(2*nl) = buffer%doublebuff(2*nl+1)
        end do
        buffer%doublebuff(2*num_soil_layers) = 0.5_rk8 * &
          (real(scalez*(exp(0.5_rkx * &
                 (real(num_soil_layers,rkx)-0.5_rkx))-1.0_rkx),rk8) + &
           real(scalez*(exp(0.5_rkx * &
                 (real(num_soil_layers+1,rkx)-0.5_rkx))-1.0_rkx),rk8))
#else
        buffer%doublebuff(1) = 0.0_rk8
        buffer%doublebuff(3) = 0.1_rk8
        buffer%doublebuff(5) = 2.0_rk8
        buffer%doublebuff(2) = 0.1_rk8
        buffer%doublebuff(4) = 2.0_rk8
        buffer%doublebuff(6) = 4.0_rk8
#endif
        call outstream_writevar(ncout,stvar%levsoilbound_var,nocopy)
      end if
      if ( stream%l_hasspectral ) then
        buffer%doublebuff(1) = 0.00000025_rk8
        buffer%doublebuff(2) = 0.00000069_rk8
        buffer%doublebuff(3) = 0.00000119_rk8
        buffer%doublebuff(4) = 0.00000238_rk8
        buffer%doublebuff(5) = 0.00000069_rk8
        buffer%doublebuff(6) = 0.00000119_rk8
        buffer%doublebuff(7) = 0.00000238_rk8
        buffer%doublebuff(8) = 0.00000400_rk8
        call outstream_writevar(ncout,stvar%spectral_var,nocopy)
      end if
    end subroutine outstream_enable

    subroutine outstream_sync(stream)
      implicit none
      type(ncoutstream) , pointer , intent(in) :: stream
      if ( .not. stream%l_enabled ) return
      if ( stream%l_sync ) then
#ifdef PNETCDF
        ncstat = nf90mpi_sync(stream%id)
#else
        ncstat = nf90_sync(stream%id)
#endif
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call printerror
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
      record = -1.0_rkx
      if ( .not. associated(ncin%ncp%xs) ) return
      stream => ncin%ncp%xs
      if ( stream%nrec > 0 ) then
        search = hourdiff(dtime,stream%refdate)
        if ( search < stream%xtime(1) .or. search > stream%xtime(2) ) then
          return
        end if
        record = 1.0_rk8+(search-stream%xtime(1))/stream%deltat
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
      stvar%time_var%rval = real(val,rk8)
      call outstream_writevar(ncout,stvar%time_var)
      if ( stream%l_hastbound ) then
        buffer%doublebuff(1) = real(stream%zero_time,rk8)
        buffer%doublebuff(2) = real(val,rk8)
        call outstream_writevar(ncout,stvar%tbound_var,nocopy)
        stream%zero_time = real(val,rk8)
      end if
    end subroutine outstream_addrec_value

    subroutine add_dimension(stream,dname)
      implicit none
      type(ncoutstream) , pointer , intent(inout) :: stream
      character(len=*) , intent(in) :: dname
      character(len=16) :: the_name , in_name
#ifdef PNETCDF
      integer(ik4) :: pdim = -1
      integer(kind=mpi_offset_kind) :: num
#else
      integer(ik4) :: pdim = -1, num
#endif
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
          if ( stream%l_bound .or. stream%l_crm ) then
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
          if ( stream%l_plev ) then
            num = npgwlev
          else
            if ( stream%l_full_sigma ) then
              ! this is the number of FULL sigma levels
              num = kz + 1
            else
              ! this is the number of HALF sigma levels
              num = kz
            end if
          end if
          the_name = 'kz'
          pdim = kz_dim
        case ('TIME','time','Time','T')
          num = nf90_unlimited
          the_name = 'time'
          pdim = time_dim
        case ('TBOUND','time_bounds','Time_bounds','TimeBounds','TB','bnds')
          num = 2
          the_name = 'bnds'
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
        case ('100M','100m','lev100m','100mlev','lev_100m')
          num = 1
          the_name = 'm100'
          pdim = h100m_level_dim
          stream%l_has100mlev = .true.
        case ('NSOIL','nsoil','n_soil_layer','nlay','layers')
          num = num_soil_layers
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
          call die('nc_stream', 'Cannot add dimension to file '// &
            trim(stream%filename)//': Undefined in add_dimension', 1)
      end select
      if ( stream%l_keep ) then
#ifdef PNETCDF
        ncstat = nf90mpi_inq_dimid(stream%id,the_name,stream%id_dims(pdim))
#else
        ncstat = nf90_inq_dimid(stream%id,the_name,stream%id_dims(pdim))
#endif
      else
#ifdef PNETCDF
        ncstat = nf90mpi_def_dim(stream%id,the_name,num,stream%id_dims(pdim))
#else
        ncstat = nf90_def_dim(stream%id,the_name,num,stream%id_dims(pdim))
#endif
      end if
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        call printerror
        call die('nc_stream', &
          'Cannot add or find dimension '//trim(the_name)//' to file '// &
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
          if ( stream%l_keep ) then
#ifdef PNETCDF
            ncstat = nf90mpi_inquire_attribute(stream%id,iv,att%aname)
#else
            ncstat = nf90_inquire_attribute(stream%id,iv,att%aname)
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_att(stream%id,iv,att%aname,att%theval)
#else
            ncstat = nf90_put_att(stream%id,iv,att%aname,att%theval)
#endif
          end if
        class is (ncattribute_logical)
          if ( stream%l_keep ) then
#ifdef PNETCDF
            ncstat = nf90mpi_inquire_attribute(stream%id,iv,att%aname)
#else
            ncstat = nf90_inquire_attribute(stream%id,iv,att%aname)
#endif
          else
            call cdumlogical(cdum,att%theval)
#ifdef PNETCDF
            ncstat = nf90mpi_put_att(stream%id,iv,att%aname,cdum)
#else
            ncstat = nf90_put_att(stream%id,iv,att%aname,cdum)
#endif
          end if
        class is (ncattribute_integer)
          if ( stream%l_keep ) then
#ifdef PNETCDF
            ncstat = nf90mpi_inquire_attribute(stream%id,iv,att%aname)
#else
            ncstat = nf90_inquire_attribute(stream%id,iv,att%aname)
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_att(stream%id,iv,att%aname,att%theval)
#else
            ncstat = nf90_put_att(stream%id,iv,att%aname,att%theval)
#endif
          end if
        class is (ncattribute_real4)
          if ( stream%l_keep ) then
#ifdef PNETCDF
            ncstat = nf90mpi_inquire_attribute(stream%id,iv,att%aname)
#else
            ncstat = nf90_inquire_attribute(stream%id,iv,att%aname)
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_att(stream%id,iv,att%aname,att%theval)
#else
            ncstat = nf90_put_att(stream%id,iv,att%aname,att%theval)
#endif
          end if
        class is (ncattribute_real8)
          if ( stream%l_keep ) then
#ifdef PNETCDF
            ncstat = nf90mpi_inquire_attribute(stream%id,iv,att%aname)
#else
            ncstat = nf90_inquire_attribute(stream%id,iv,att%aname)
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_att(stream%id,iv,att%aname,att%theval)
#else
            ncstat = nf90_put_att(stream%id,iv,att%aname,att%theval)
#endif
          end if
        class is (ncattribute_real4_array)
          if ( stream%l_keep ) then
#ifdef PNETCDF
            ncstat = nf90mpi_inquire_attribute(stream%id,iv,att%aname)
#else
            ncstat = nf90_inquire_attribute(stream%id,iv,att%aname)
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_att(stream%id,iv, &
                      att%aname,att%theval(1:att%numval))
#else
            ncstat = nf90_put_att(stream%id,iv, &
                      att%aname,att%theval(1:att%numval))
#endif
          end if
        class is (ncattribute_real8_array)
          if ( stream%l_keep ) then
#ifdef PNETCDF
            ncstat = nf90mpi_inquire_attribute(stream%id,iv,att%aname)
#else
            ncstat = nf90_inquire_attribute(stream%id,iv,att%aname)
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_att(stream%id,iv, &
                       att%aname,att%theval(1:att%numval))
#else
            ncstat = nf90_put_att(stream%id,iv, &
                       att%aname,att%theval(1:att%numval))
#endif
          end if
        class default
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', 'Cannot add attribute of unknow type',1)
      end select
      if ( ncstat /= nf90_noerr ) then
        call printerror
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
      type(ncoutstream) , pointer , intent(inout) :: stream
      class(ncvariable_standard) , intent(in) :: var
      character(len=*) , parameter :: coords_cross = 'xlat xlon'
      character(len=*) , parameter :: coords_depth = 'soil_layer xlat xlon'
      character(len=*) , parameter :: coords_dot   = 'dlat dlon'
      character(len=*) , parameter :: coords_udot  = 'ulat ulon'
      character(len=*) , parameter :: coords_vdot  = 'vlat vlon'
      if ( len_trim(var%long_name) > 0 ) &
        call add_attribute(stream, &
          ncattribute_string('long_name',var%long_name),var%id,var%vname)
      if ( len_trim(var%standard_name) > 0 ) &
        call add_attribute(stream, &
         ncattribute_string('standard_name',var%standard_name),var%id,var%vname)
      if ( len_trim(var%vunit) > 0 ) &
        call add_attribute(stream, &
          ncattribute_string('units',var%vunit),var%id,var%vname)
      if ( len_trim(var%notes) > 0 ) &
        call add_attribute(stream, &
          ncattribute_string('notes',var%notes),var%id,var%vname)
      if ( var%lgridded ) then
        if ( var%vname(2:5) /= 'lat' .and. var%vname(2:5) /= 'lon' ) then
          if ( stream%l_bound ) then
            if ( idynamic < 3 ) then
              if (var%vname == 'u' .or. var%vname == 'v') then
                call add_attribute(stream, &
                        ncattribute_string('coordinates', &
                        coords_dot),var%id,var%vname)
              else
                call add_attribute(stream, &
                        ncattribute_string('coordinates', &
                        coords_cross),var%id,var%vname)
              end if
            else
              if ( var%vname == 'u' ) then
                call add_attribute(stream, &
                        ncattribute_string('coordinates', &
                        coords_udot),var%id,var%vname)
              else if ( var%vname == 'v') then
                call add_attribute(stream, &
                        ncattribute_string('coordinates', &
                        coords_vdot),var%id,var%vname)
              else
                call add_attribute(stream, &
                        ncattribute_string('coordinates', &
                        coords_cross),var%id,var%vname)
              end if
            end if
          else
            if ( var%vname == 'mrsos' ) then
              call add_attribute(stream, &
                ncattribute_string('coordinates',coords_depth),var%id,var%vname)
            else
              call add_attribute(stream, &
                ncattribute_string('coordinates',coords_cross),var%id,var%vname)
            end if
          end if
          call add_attribute(stream, &
            ncattribute_string('grid_mapping','crs'),var%id,var%vname)
          stream%l_hasgrid = .true.
        end if
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
        else if ( var%nctype == nf90_double ) then
          call add_attribute(stream, &
            ncattribute_real8('_FillValue',var%rmissval),var%id,var%vname)
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
        if ( stream%l_keep ) then
#ifdef PNETCDF
          ncstat = nf90mpi_inq_varid(stream%id,var%vname,var%id)
#else
          ncstat = nf90_inq_varid(stream%id,var%vname,var%id)
#endif
        else
#ifdef PNETCDF
          ncstat = nf90mpi_def_var(stream%id,var%vname,var%nctype,var%id)
#else
          ncstat = nf90_def_var(stream%id,var%vname,var%nctype,var%id)
#endif
        end if
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call printerror
          call die('nc_stream', &
            'Cannot add variable '//trim(var%vname)//' to file '// &
            trim(stream%filename), 1)
        end if
      else
        if ( stream%l_keep ) then
#ifdef PNETCDF
          ncstat = nf90mpi_inq_varid(stream%id,var%vname,var%id)
#else
          ncstat = nf90_inq_varid(stream%id,var%vname,var%id)
#endif
        else
#ifdef PNETCDF
          ncstat = nf90mpi_def_var(stream%id,var%vname,var%nctype, &
                                id_dim(1:ndims),var%id)
#else
          ncstat = nf90_def_var(stream%id,var%vname,var%nctype, &
                                id_dim(1:ndims),var%id)
#endif
        end if
        if ( ncstat /= nf90_noerr ) then
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call printerror
          call die('nc_stream', &
            'Cannot define variable '//trim(var%vname)// &
            ' in file '//trim(stream%filename), 1)
        end if
#if defined(NETCDF4_HDF5)
#if defined (NETCDF4_COMPRESS)
        if ( ndims > 3 ) then
          if ( stream%l_keep ) then
            ncstat = nf90_inq_varid(stream%id,var%vname,var%id)
          else
            if ( .not. stream%l_parallel ) then
              ncstat = nf90_def_var_deflate(stream%id,var%id,1,1,deflate_level)
              if ( ncstat /= nf90_noerr ) then
                write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
                call printerror
                call die('nc_stream', &
                  'Cannot set compression on variable '//trim(var%vname)// &
                  ' in file '//trim(stream%filename), 1)
              end if
            end if
          end if
        end if
#endif
        ! This forces collective I/O on time dependent variables.
        if ( stream%l_parallel .and. var%lrecords ) then
          ncstat = nf90_var_par_access(stream%id,var%id,nf90_collective)
          if ( ncstat /= nf90_noerr ) then
            write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
            call printerror
            call die('nc_stream', &
              'Cannot set correct mode (collective) for variable I/O'// &
              'on var'//trim(var%vname)//' in file '// &
              trim(stream%filename), 1)
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
        class is (ncvariable0d_mixed)
          var%nctype = cordexvtype
        class is (ncvariable1d_mixed)
          var%nctype = cordexvtype
        class is (ncvariable2d_mixed)
          var%nctype = cordexvtype
        class is (ncvariable3d_mixed)
          var%nctype = cordexvtype
        class is (ncvariable4d_mixed)
          var%nctype = cordexvtype
        class is (ncvariable0d_real)
          var%nctype = cordexvtype
        class is (ncvariable1d_real)
          var%nctype = cordexvtype
        class is (ncvariable2d_real)
          var%nctype = cordexvtype
        class is (ncvariable3d_real)
          var%nctype = cordexvtype
        class is (ncvariable4d_real)
          var%nctype = cordexvtype
        class is (ncvariable0d_double)
          var%nctype = nf90_double
        class is (ncvariable1d_double)
          var%nctype = nf90_double
        class is (ncvariable2d_double)
          var%nctype = nf90_double
        class is (ncvariable3d_double)
          var%nctype = nf90_double
        class is (ncvariable4d_double)
          var%nctype = nf90_double
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
          call die('nc_stream', 'Cannot add variable of unknown type',1)
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
          call die('nc_stream', 'Cannot add variable of unknown type',1)
      end select
      select type(var)
        class is (ncvariable1d_double)
          buffer%lhas1ddouble = .true.
          buffer%max1d_double(1) = max(buffer%max1d_double(1),len_dim(1))
        class is (ncvariable2d_double)
          buffer%lhas2ddouble = .true.
          if ( stream%l_parallel .and. var%lgridded ) then
            buffer%max2d_double(1) = &
               max(buffer%max2d_double(1),stream%global_nj)
            buffer%max2d_double(2) = &
               max(buffer%max2d_double(2),stream%global_ni)
          else
            buffer%max2d_double(1) = max(buffer%max2d_double(1),len_dim(1))
            buffer%max2d_double(2) = max(buffer%max2d_double(2),len_dim(2))
          end if
        class is (ncvariable3d_double)
          buffer%lhas3ddouble = .true.
          if ( stream%l_parallel .and. var%lgridded ) then
            buffer%max3d_double(1) = &
               max(buffer%max3d_double(1),stream%global_nj)
            buffer%max3d_double(2) = &
               max(buffer%max3d_double(2),stream%global_ni)
          else
            buffer%max3d_double(1) = max(buffer%max3d_double(1),len_dim(1))
            buffer%max3d_double(2) = max(buffer%max3d_double(2),len_dim(2))
          end if
          buffer%max3d_double(3) = max(buffer%max3d_double(3),len_dim(3))
        class is (ncvariable4d_double)
          buffer%lhas4ddouble = .true.
          if ( stream%l_parallel .and. var%lgridded ) then
            buffer%max4d_double(1) = &
               max(buffer%max4d_double(1),stream%global_nj)
            buffer%max4d_double(2) = &
               max(buffer%max4d_double(2),stream%global_ni)
          else
            buffer%max4d_double(1) = max(buffer%max4d_double(1),len_dim(1))
            buffer%max4d_double(2) = max(buffer%max4d_double(2),len_dim(2))
          end if
          buffer%max4d_double(3) = max(buffer%max4d_double(3),len_dim(3))
          buffer%max4d_double(4) = max(buffer%max4d_double(4),len_dim(4))
        class is (ncvariable1d_real)
          buffer%lhas1dreal = .true.
          buffer%max1d_real(1) = max(buffer%max1d_real(1),len_dim(1))
        class is (ncvariable1d_mixed)
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
        class is (ncvariable2d_mixed)
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
        class is (ncvariable3d_mixed)
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
        class is (ncvariable4d_mixed)
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
          case ('W')
            if ( stream%id_dims(h100m_level_dim) < 0 ) then
              call add_dimension(stream,'100m')
            end if
            id_dim(ic) = stream%id_dims(h100m_level_dim)
            len_dim(ic) = stream%len_dims(h100m_level_dim)
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
        class is (ncvariable0d_double)
          if ( var%lrecords ) then
            stream%istart(1) = stream%irec
            stream%icount(1) = 1
            nd = 1
#ifdef PNETCDF
            ncstat = nf90mpi_put_var_all(stream%id,var%id,var%rval, &
              stream%istart(1:nd),stream%icount(1:nd))
#else
            ncstat = nf90_put_var(stream%id,var%id,var%rval, &
              stream%istart(1:nd),stream%icount(1:nd))
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_var_all(stream%id,var%id,var%rval(1))
#else
            ncstat = nf90_put_var(stream%id,var%id,var%rval(1))
#endif
          end if
        class is (ncvariable1d_double)
          if ( docopy ) then
            if ( .not. associated(var%rval) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%doublebuff(1:var%nval(1)) = var%rval(1:var%nval(1))
          end if
          stream%istart(1) = 1
          stream%icount(1) = var%nval(1)
          nd = 1
          if ( var%lrecords ) then
            stream%istart(2) = stream%irec
            stream%icount(2) = 1
            nd = 2
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%doublebuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%doublebuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable2d_double)
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
              buffer%doublebuff(1:totsize) =  &
                reshape(var%rval_slice(var%j1:var%j2,var%i1:var%i2,is), &
                [totsize])
            else
              buffer%doublebuff(1:totsize) =  &
                reshape(var%rval(var%j1:var%j2,var%i1:var%i2),[totsize])
            end if
          end if
          nd = 2
          if ( var%lrecords ) then
            stream%istart(3) = stream%irec
            stream%icount(3) = 1
            nd = 3
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%doublebuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%doublebuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable3d_double)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize*var%nval(3)
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
              buffer%doublebuff(1:totsize) = &
                reshape(var%rval_slice(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2,is), [totsize])
            else
              buffer%doublebuff(1:totsize) = &
                reshape(var%rval(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2),[totsize])
            end if
          end if
          nd = 3
          if ( var%lrecords ) then
            stream%istart(4) = stream%irec
            stream%icount(4) = 1
            nd = 4
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%doublebuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%doublebuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable4d_double)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize*var%nval(3)*var%nval(4)
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
            buffer%doublebuff(1:totsize) = &
              reshape(var%rval(var%j1:var%j2,var%i1:var%i2, &
                var%k1:var%k2,var%n1:var%n2),[totsize])
          end if
          nd = 4
          if ( var%lrecords ) then
            stream%istart(5) = stream%irec
            stream%icount(5) = 1
            nd = 5
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%doublebuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%doublebuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable0d_real)
          if ( var%lrecords ) then
            stream%istart(1) = stream%irec
            stream%icount(1) = 1
#ifdef PNETCDF
            ncstat = nf90mpi_put_var_all(stream%id,var%id,var%rval, &
              stream%istart(1:1),stream%icount(1:1))
#else
            ncstat = nf90_put_var(stream%id,var%id,var%rval, &
              stream%istart(1:1),stream%icount(1:1))
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_var_all(stream%id,var%id,var%rval(1))
#else
            ncstat = nf90_put_var(stream%id,var%id,var%rval(1))
#endif
          end if
        class is (ncvariable0d_mixed)
          if ( var%lrecords ) then
            stream%istart(1) = stream%irec
            stream%icount(1) = 1
#ifdef PNETCDF
            ncstat = nf90mpi_put_var_all(stream%id,var%id,var%rval, &
              stream%istart(1:1),stream%icount(1:1))
#else
            ncstat = nf90_put_var(stream%id,var%id,var%rval, &
              stream%istart(1:1),stream%icount(1:1))
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_var_all(stream%id,var%id,var%rval(1))
#else
            ncstat = nf90_put_var(stream%id,var%id,var%rval(1))
#endif
          end if
        class is (ncvariable1d_real)
          if ( docopy ) then
            if ( .not. associated(var%rval) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%realbuff(1:var%nval(1)) = var%rval(1:var%nval(1))
          end if
          stream%istart(1) = 1
          stream%icount(1) = var%nval(1)
          nd = 1
          if ( var%lrecords ) then
            stream%istart(2) = stream%irec
            stream%icount(2) = 1
            nd = 2
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable1d_mixed)
          if ( docopy ) then
            if ( .not. associated(var%rval) ) then
              write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
              write(stderr,*) 'Unassociated pointer to variable'
              call die('nc_stream','Cannot write variable '//trim(var%vname)// &
                ' in file '//trim(stream%filename), 1)
            end if
            buffer%realbuff(1:var%nval(1)) = real(var%rval(1:var%nval(1)),rk4)
          end if
          stream%istart(1) = 1
          stream%icount(1) = var%nval(1)
          nd = 1
          if ( var%lrecords ) then
            stream%istart(2) = stream%irec
            stream%icount(2) = 1
            nd = 2
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
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
                reshape(var%rval_slice(var%j1:var%j2,var%i1:var%i2,is), &
                [totsize])
            else
              buffer%realbuff(1:totsize) =                          &
                reshape(var%rval(var%j1:var%j2,var%i1:var%i2), &
                [totsize])
            end if
          end if
          nd = 2
          if ( var%lrecords ) then
            stream%istart(3) = stream%irec
            stream%icount(3) = 1
            nd = 3
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable2d_mixed)
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
                [totsize]),rk4)
            else
              buffer%realbuff(1:totsize) =                          &
                real(reshape(var%rval(var%j1:var%j2,var%i1:var%i2), &
                [totsize]),rk4)
            end if
          end if
          nd = 2
          if ( var%lrecords ) then
            stream%istart(3) = stream%irec
            stream%icount(3) = 1
            nd = 3
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable3d_real)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize*var%nval(3)
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
                reshape(var%rval_slice(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2,is), [totsize])
            else
              buffer%realbuff(1:totsize) = &
                reshape(var%rval(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2),[totsize])
            end if
          end if
          nd = 3
          if ( var%lrecords ) then
            stream%istart(4) = stream%irec
            stream%icount(4) = 1
            nd = 4
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable3d_mixed)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize*var%nval(3)
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
                  var%k1:var%k2,is), [totsize]),rk4)
            else
              buffer%realbuff(1:totsize) = &
                real(reshape(var%rval(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2),[totsize]),rk4)
            end if
          end if
          nd = 3
          if ( var%lrecords ) then
            stream%istart(4) = stream%irec
            stream%icount(4) = 1
            nd = 4
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable4d_real)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize*var%nval(3)*var%nval(4)
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
              reshape(var%rval(var%j1:var%j2,var%i1:var%i2, &
                var%k1:var%k2,var%n1:var%n2),[totsize])
          end if
          nd = 4
          if ( var%lrecords ) then
            stream%istart(5) = stream%irec
            stream%icount(5) = 1
            nd = 5
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable4d_mixed)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize*var%nval(3)*var%nval(4)
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
                var%k1:var%k2,var%n1:var%n2),[totsize]),rk4)
          end if
          nd = 4
          if ( var%lrecords ) then
            stream%istart(5) = stream%irec
            stream%icount(5) = 1
            nd = 5
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%realbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable0d_integer)
          if ( var%lrecords ) then
            stream%istart(1) = stream%irec
            stream%icount(1) = 1
#ifdef PNETCDF
            ncstat = nf90mpi_put_var_all(stream%id,var%id,var%ival, &
              stream%istart(1:1),stream%icount(1:1))
#else
            ncstat = nf90_put_var(stream%id,var%id,var%ival, &
              stream%istart(1:1),stream%icount(1:1))
#endif
          else
#ifdef PNETCDF
            ncstat = nf90mpi_put_var_all(stream%id,var%id,var%ival(1))
#else
            ncstat = nf90_put_var(stream%id,var%id,var%ival(1))
#endif
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
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
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
                [totsize])
            else
              buffer%intbuff(1:totsize) =                               &
                reshape(var%ival(var%j1:var%j2,var%i1:var%i2),[totsize])
            end if
          end if
          nd = 2
          if ( var%lrecords ) then
            stream%istart(3) = stream%irec
            stream%icount(3) = 1
            nd = 3
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable3d_integer)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize*var%nval(3)
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
                  var%k1:var%k2,is), [totsize])
            else
              buffer%realbuff(1:totsize) =                    &
                reshape(var%ival(var%j1:var%j2,var%i1:var%i2, &
                  var%k1:var%k2),[totsize])
            end if
          end if
          nd = 3
          if ( var%lrecords ) then
            stream%istart(4) = stream%irec
            stream%icount(4) = 1
            nd = 4
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class is (ncvariable4d_integer)
          if ( stream%l_parallel .and. var%lgridded ) then
            stream%istart(1) = stream%jparbound(1)
            stream%istart(2) = stream%iparbound(1)
            stream%icount(1) = stream%global_nj
            stream%icount(2) = stream%global_ni
            totsize = stream%parsize*var%nval(3)*var%nval(4)
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
                var%k1:var%k2,var%n1:var%n2),[totsize])
          end if
          nd = 4
          if ( var%lrecords ) then
            stream%istart(5) = stream%irec
            stream%icount(5) = 1
            nd = 5
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_put_var_all(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
#else
          ncstat = nf90_put_var(stream%id,var%id, &
            buffer%intbuff,stream%istart(1:nd),stream%icount(1:nd))
#endif
        class default
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream', 'Cannot write variable of unknown type',1)
      end select
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        call printerror
        write(stderr,*) 'In file ',trim(stream%filename)
        write(stderr,*) 'Cannot write variable '//trim(var%vname)
        write(stderr,*) 'Bounds :',stream%istart(1:nd),stream%icount(1:nd)
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
      attc%aname = 'executable_date'
      attc%theval = __DATE__
      call add_attribute(stream,attc)
      attc%aname = 'Conventions'
      attc%theval = 'CF-1.7'
      call add_attribute(stream,attc)
      attc%aname = 'references'
      attc%theval = 'https://github.com/ICTP/RegCM'
      call add_attribute(stream,attc)
      attc%aname = 'model_revision'
      attc%theval = GIT_VER
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
        attr%theval = (ds*d_1000)/real(nsg,rkx)
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
      attra%aname = 'index_of_projection_origin'
      attra%theval(1) = cntrj
      attra%theval(2) = cntri
      attra%numval = 2
      call add_attribute(stream,attra)
      if ( iproj == 'ROTMER' .or. iproj == 'ROTLLR' ) then
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
      if ( iproj /= 'ROTLLR' ) then
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
      else
        stvar%jx_var%vname = 'rlon'
        stvar%jx_var%vunit = 'degrees'
        stvar%jx_var%long_name = 'longitude in rotated pole grid'
        stvar%jx_var%standard_name = 'grid_longitude'
        stvar%jx_var%axis = 'x'
        stvar%iy_var%vname = 'rlat'
        stvar%iy_var%vunit = 'degrees'
        stvar%iy_var%long_name = 'latitude in rotated pole grid'
        stvar%iy_var%standard_name = 'grid_latitude'
        stvar%iy_var%axis = 'y'
      end if
      if ( stream%l_plev ) then
        stvar%sigma_var%vname = 'plev'
        stvar%sigma_var%vunit = 'Pa'
        stvar%sigma_var%long_name = "Pressure"
        stvar%sigma_var%standard_name = 'atmosphere_pressure_coordinate'
        stvar%sigma_var%axis = 'z'
      else
        stvar%sigma_var%vname = 'kz'
        stvar%sigma_var%vunit = '1'
        if ( stream%l_full_sigma ) then
          stvar%sigma_var%long_name = "sigma at layer interfaces"
        else
          stvar%sigma_var%long_name = "sigma at layer midpoints"
        end if
        if ( idynamic < 3 ) then
          stvar%sigma_var%standard_name = 'atmosphere_sigma_coordinate'
        else
          stvar%sigma_var%standard_name = 'atmosphere_hybrid_height_coordinate'
        end if
        stvar%sigma_var%axis = 'z'
        if ( idynamic < 3 ) then
          stvar%ptop_var%vname = 'ptop'
          stvar%ptop_var%vunit = 'hPa'
          stvar%ptop_var%long_name = "Pressure at model top"
          stvar%ptop_var%standard_name = 'air_pressure'
          call outstream_addvar(ncout,stvar%ptop_var)
        else
          stvar%ak_var%vname = 'a'
          stvar%ak_var%vunit = 'm'
          stvar%ak_var%axis = 'z'
          stvar%ak_var%standard_name = "atmosphere_hybrid_height_coordinate"
          stvar%ak_var%long_name = "vertical coordinate formula term a(k)"
          call outstream_addvar(ncout,stvar%ak_var)
          stvar%bk_var%vname = 'b'
          stvar%bk_var%vunit = '1'
          stvar%bk_var%axis = 'z'
          stvar%bk_var%standard_name = "atmosphere_hybrid_height_coordinate"
          stvar%bk_var%long_name = "vertical coordinate formula term b(k)"
          call outstream_addvar(ncout,stvar%bk_var)
        end if
      end if
      call outstream_addvar(ncout,stvar%jx_var)
      call outstream_addvar(ncout,stvar%iy_var)
      call outstream_addvar(ncout,stvar%sigma_var)
      attc%aname = 'axis'
      attc%theval = 'X'
      call add_attribute(stream,attc,stvar%jx_var%id,stvar%jx_var%vname)
      attc%theval = 'Y'
      call add_attribute(stream,attc,stvar%iy_var%id,stvar%iy_var%vname)
      attc%theval = 'Z'
      call add_attribute(stream,attc,stvar%sigma_var%id,stvar%sigma_var%vname)
      attc%aname = 'CoordinateAxisType'
      attc%theval = 'GeoX'
      call add_attribute(stream,attc,stvar%jx_var%id,stvar%jx_var%vname)
      attc%aname = 'CoordinateAxisType'
      attc%theval = 'GeoY'
      call add_attribute(stream,attc,stvar%iy_var%id,stvar%iy_var%vname)
      attc%aname = 'positive'
      attc%theval = 'down'
      call add_attribute(stream,attc,stvar%sigma_var%id,stvar%sigma_var%vname)
      if ( .not. stream%l_plev ) then
        if ( idynamic == 3 ) then
          attc%aname = 'formula'
          attc%theval = 'z(k,j,i) = a(k) + topo(j,i) * b(k)'
        else if ( idynamic == 2 ) then
          attc%aname = 'formula'
          attc%theval = 'p(n,k,j,i) = ptop + kz(k)*(p0(j,i)-ptop)+ppa(n,k,j,i)'
        else
          attc%aname = 'formula'
          attc%theval = 'p(n,k,j,i) = ptop + kz(k)*(ps(n,j,i)-ptop)'
        end if
        call add_attribute(stream,attc,stvar%sigma_var%id,stvar%sigma_var%vname)
      end if
      attc%aname = 'CoordinateAxisType'
      attc%theval = 'GeoZ'
      call add_attribute(stream,attc,stvar%sigma_var%id,stvar%sigma_var%vname)
      if ( idynamic == 3 ) then
        attr%aname = 'zita_factor_a0'
        attr%theval = mo_a0
        call add_attribute(stream,attr)
        attr%aname = 'zita_height_top'
        attr%theval = mo_ztop
        call add_attribute(stream,attr)
        attr%aname = 'zita_atmosphere_h'
        attr%theval = mo_h
        call add_attribute(stream,attr)
      end if
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
#ifdef PNETCDF
        ncstat = nf90mpi_inq_varid(stream%id,var%vname,var%id)
#else
        ncstat = nf90_inq_varid(stream%id,var%vname,var%id)
#endif
        if ( ncstat /= nf90_noerr ) then
          call printerror
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream','Cannot find variable '//trim(var%vname)// &
            ' in '//trim(stream%filename),1)
        end if
#ifdef PNETCDF
        ncstat = nf90mpi_inquire_variable(stream%id,var%id,ndims=var%ndims, &
          dimids=var%idims)
#else
        ncstat = nf90_inquire_variable(stream%id,var%id,ndims=var%ndims, &
          dimids=var%idims)
#endif
        if ( ncstat /= nf90_noerr ) then
          call printerror
          write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
          call die('nc_stream','Cannot inquire variable '//trim(var%vname)// &
            ' in '//trim(stream%filename),1)
        end if
      end if
      select type(var)
        class is (ncvariable0d_double)
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,var%rval)
#else
          ncstat = nf90_get_var(stream%id,var%id,var%rval)
#endif
        class is (ncvariable1d_double)
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,var%rval)
#else
          ncstat = nf90_get_var(stream%id,var%id,var%rval)
#endif
        class is (ncvariable2d_double)
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
          if ( size(buffer%doublebuff) < var%totsize ) then
            deallocate(buffer%doublebuff)
            allocate(buffer%doublebuff(var%totsize))
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%doublebuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%doublebuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
              reshape(buffer%doublebuff(1:var%totsize),stream%icount(1:2))
          else
            var%rval(1:var%nval(1),1:var%nval(2)) = &
              reshape(buffer%doublebuff(1:var%totsize),stream%icount(1:2))
          end if
        class is (ncvariable3d_double)
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
          if ( size(buffer%doublebuff) < var%totsize ) then
            deallocate(buffer%doublebuff)
            allocate(buffer%doublebuff(var%totsize))
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%doublebuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%doublebuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
              reshape(buffer%doublebuff(1:var%totsize),stream%icount(1:3))
          else
            var%rval(1:var%nval(1),1:var%nval(2),1:var%nval(3)) = &
              reshape(buffer%doublebuff(1:var%totsize),stream%icount(1:3))
          end if
        class is (ncvariable4d_double)
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
          if ( size(buffer%doublebuff) < var%totsize ) then
            deallocate(buffer%doublebuff)
            allocate(buffer%doublebuff(var%totsize))
          end if
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%doublebuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%doublebuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
                 reshape(buffer%doublebuff(1:var%totsize),stream%icount(1:4))
          else
            var%rval(1:var%nval(1),1:var%nval(2), &
                     1:var%nval(3),1:var%nval(4)) = &
                 reshape(buffer%doublebuff(1:var%totsize),stream%icount(1:4))
          end if
        class is (ncvariable0d_real)
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,var%rval)
#else
          ncstat = nf90_get_var(stream%id,var%id,var%rval)
#endif
        class is (ncvariable1d_real)
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,var%rval)
#else
          ncstat = nf90_get_var(stream%id,var%id,var%rval)
#endif
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
        class is (ncvariable0d_mixed)
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,var%rval)
#else
          ncstat = nf90_get_var(stream%id,var%id,var%rval)
#endif
        class is (ncvariable1d_mixed)
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,var%rval)
#else
          ncstat = nf90_get_var(stream%id,var%id,var%rval)
#endif
        class is (ncvariable2d_mixed)
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
        class is (ncvariable3d_mixed)
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
        class is (ncvariable4d_mixed)
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%realbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,var%ival)
#else
          ncstat = nf90_get_var(stream%id,var%id,var%ival)
#endif
        class is (ncvariable1d_integer)
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,var%ival)
#else
          ncstat = nf90_get_var(stream%id,var%id,var%ival)
#endif
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%intbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%intbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%intbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%intbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
#ifdef PNETCDF
          ncstat = nf90mpi_get_var(stream%id,var%id,buffer%intbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#else
          ncstat = nf90_get_var(stream%id,var%id,buffer%intbuff, &
            stream%istart(1:ndims),stream%icount(1:ndims))
#endif
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
          call die('nc_stream', 'Cannot read variable of unknown type',1)
      end select
      if ( ncstat /= nf90_noerr ) then
        call printerror
        write(stderr,*) 'In File ',__FILE__,' at line: ',__LINE__
        call die('nc_stream','Cannot read variable '//trim(var%vname)// &
          ' in '//trim(stream%filename),1)
      end if
    end subroutine instream_readvar

    subroutine printerror
      implicit none
#ifdef PNETCDF
      write(stderr, *) nf90mpi_strerror(ncstat)
#else
      write(stderr,*) nf90_strerror(ncstat)
#endif
    end subroutine printerror

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
  real(rkx) :: xrec

  type(ncvariable0d_integer) :: var0dint
  type(ncvariable1d_real) :: var1dreal
  type(ncvariable2d_real) :: var2dreal
  type(ncvariable3d_real) :: var3dreal

  type(ncvariable2d_real) :: var2read

  real(rkx) , target , dimension(18) :: sigma
  real(rkx) , target , dimension(16,12) :: d2dvar
  real(rkx) , target , dimension(16,12,18) :: d3dvar
  integer(ik4) :: i , k

  data sigma /0.025_rkx, 0.075_rkx, 0.130_rkx, &
              0.195_rkx, 0.270_rkx, 0.359_rkx, &
              0.430_rkx, 0.510_rkx, 0.590_rkx, &
              0.670_rkx, 0.750_rkx, 0.810_rkx, &
              0.865_rkx, 0.910_rkx, 0.950_rkx, &
              0.970_rkx, 0.985_rkx, 0.995_rkx /
  jx = 16
  iy = 12
  kz = 18
  ds = 50.0_rkx
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
    ncattribute_real8('minimum_h2o_pct_for_water',50.0_rkx))

  ! Enable the output stream for write
  call outstream_enable(ncout,sigma)

  var1dreal%rval => sigma
  d2dvar(1:jx,1:iy) = 1.0_rkx
  d2dvar(jx/2,iy/2) = 2.0_rkx
  var2dreal%rval => d2dvar

  ! Write some static variables
  call outstream_writevar(ncout,var1dreal)
  call outstream_writevar(ncout,var2dreal)

  var0dint%ival(1) = 12
  d3dvar(1:jx,1:iy,:) = 1.0_rkx
  do k = 1 , kz
    d3dvar(jx/4,iy/4,k) = real(k,rkx)
  end do
  var3dreal%rval => d3dvar

  ! Write variables in the current record step
  idate = 1979022206
  call outstream_addrec(ncout,idate)
  call outstream_writevar(ncout,var0dint)
  call outstream_writevar(ncout,var3dreal)

  var0dint%ival(1) = 13
  d3dvar(1:jx,1:iy,:) = 1.0_rkx
  do k = 1 , kz
    d3dvar(jx/4,iy/4,k) = real(k,rkx)*1.5_rkx
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

  d2dvar(:,:) = -1.0_rkx
  var2read%vname = 'real2dvar'
  var2read%lgridded = .true.
  var2read%rval => d2dvar
  var2read%j1 = 1
  var2read%j2 = 9
  var2read%i1 = 1
  var2read%i2 = 7

  call instream_readvar(ncin,var2read,window=[4,12,2,8])

  do i = 1 , iy
    print '(16i2)', int(d2dvar(:,i))
  end do

  call instream_dispose(ncin)

end program test
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

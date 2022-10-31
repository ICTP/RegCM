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

module mod_sst_cmip6

  use mod_intkinds
  use mod_realkinds
  use mod_cmip6_helper
  use mod_message
  use mod_dynparam
  use mod_memutil
  use mod_sst_grid
  use mod_kdinterp
  use mod_date
  use mod_stdio
  use mod_cmip6_cesm
  use mod_cmip6_canesm
  use mod_cmip6_cnrm
  use mod_cmip6_ecea
  use mod_cmip6_gfdl
  use mod_cmip6_hadmm
  use mod_cmip6_miroc6
  use mod_cmip6_miresl
  use mod_cmip6_mpihr
  use mod_cmip6_normm
  use mod_cmip6_cmcc
  use netcdf

  implicit none

  private

  public :: cmip6_sst

  type(cmip6_2d_var) :: sst

  abstract interface
    subroutine read_cmip6_sst(id,var,lat,lon)
      import
      type(rcm_time_and_date) , intent(in) :: id
      type(cmip6_2d_var) , intent(inout) :: var
      real(rkx) , dimension(:,:) , pointer , intent(in) :: lat , lon
    end subroutine read_cmip6_sst
  end interface

  contains

    subroutine cmip6_sst
      implicit none
      type(rcm_time_and_date) :: idate , idatef , idateo
      type(rcm_time_interval) :: tdif , step
      procedure(read_cmip6_sst) , pointer :: read_func
      integer :: nsteps , n

      idateo = globidate1
      idatef = globidate2
      tdif = idatef-idateo

      select case (cmip6_model)
        case ('MPI-ESM1-2-HR')
          if ( calendar /= 'gregorian' ) then
            write(stderr,*) 'MPI-ESM1-2-HR requires gregorian calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_mpihr
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'HadGEM3-GC31-MM' )
          if ( calendar /= '360_day' ) then
            write(stderr,*) 'HadGEM3-GC31-MM requires 360_day calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_hadmm
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'GFDL-ESM4' )
          if ( calendar /= 'noleap' ) then
            write(stderr,*) 'GFDL-ESM4 requires noleap calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_gfdl
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'NorESM2-MM' )
          if ( calendar /= 'noleap' ) then
            write(stderr,*) 'NorESM2-MM requires noleap calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_normm
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'CNRM-ESM2-1' )
          if ( calendar /= 'gregorian' ) then
            write(stderr,*) 'CNRM-ESM2-1 requires gregorian calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_cnrm
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'EC-Earth3-Veg' )
          if ( calendar /= 'gregorian' ) then
            write(stderr,*) 'EC-Earth3-Veg requires gregorian calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_ecea
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'CESM2' )
          if ( calendar /= 'noleap' ) then
            write(stderr,*) 'CESM2 requires noleap calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_cesm
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'CMCC-ESM2' )
          if ( calendar /= 'noleap' ) then
            write(stderr,*) 'CMCC-ESM2 requires noleap calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_cmcc
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'CanESM5' )
          if ( calendar /= 'noleap' ) then
            write(stderr,*) 'CanESM5 requires noleap calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_canesm
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'MIROC6' )
          if ( calendar /= 'gregorian' ) then
            write(stderr,*) 'MIROC6 requires gregorian calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_miroc6
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case ( 'MIROC-ES2L' )
          if ( calendar /= 'gregorian' ) then
            write(stderr,*) 'MIROC-ES2L requires gregorian calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_miresl
          sst%vname = 'tos'
          idateo = prevmon(globidate1)
          idatef = monfirst(globidate2)
          if (idatef < globidate2) then
            idatef = nextmon(idatef)
          end if
          step = rcm_time_interval(1_ik8,umnt)
          nsteps = imondiff(idatef,idateo) + 1
        case default
          call die('sst','Unknown CMIP6 model: '//trim(cmip6_model),1)
      end select

      write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
      write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
      write (stdout,*) 'NSTEPS = ', nsteps

      call open_sstfile(idateo)

      allocate(sst%hint(1))

      idate = idateo
      do n = 1 , nsteps
        call read_func(idate,sst,xlat,xlon)
        call h_interpolate_cont(sst%hint(1),sst%var,sstmm)
        call writerec(idate)
        write (stdout,*) 'WRITEN OUT SST DATA : ' , tochar(idate)
        idate = idate + step
      end do

      call h_interpolator_destroy(sst%hint(1))
      deallocate(sst%hint)
    end subroutine cmip6_sst

end module mod_sst_cmip6

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

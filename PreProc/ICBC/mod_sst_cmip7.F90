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

module mod_sst_cmip7

  use mod_intkinds
  use mod_realkinds
  use mod_cmip7_helper
  use mod_message
  use mod_dynparam
  use mod_memutil
  use mod_sst_grid
  use mod_kdinterp
  use mod_date
  use mod_stdio
  use mod_cmip7_ecea
  use netcdf

  implicit none

  private

  public :: cmip7_sst

  type(cmip7_2d_var) :: sst

  abstract interface
    subroutine read_cmip7_sst(id,var,lat,lon)
      import
      implicit none
      type(rcm_time_and_date), intent(in) :: id
      type(cmip7_2d_var), intent(inout) :: var
      real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: lat, lon
    end subroutine read_cmip7_sst
  end interface

  contains

    subroutine cmip7_sst
      implicit none
      type(rcm_time_and_date) :: idate, idatef, idateo
      type(rcm_time_interval) :: tdif, step
      procedure(read_cmip7_sst), pointer :: read_func => null( )
      integer :: nsteps, n

      idateo = globidate1
      idatef = globidate2
      tdif = idatef-idateo

      dattyp = ssttyp

      select case (cmip7_model)
        case ( 'EC-Earth3-ESM-1-1' )
          if ( ical /= gregorian ) then
            write(stderr,*) 'EC-Earth3-ESM-1-1 requires gregorian calendar.'
            call die('sst','Calendar mismatch',1)
          end if
          read_func => read_sst_ecea
          sst%vname = 'tos'
          step = 86400
          nsteps = int(tohours(tdif))/24 + 1
        case default
          call die('sst','Unknown CMIP7 model: '//trim(cmip7_model),1)
      end select

      write (stdout,*) 'GLOBIDATE1 : ', tochar(globidate1)
      write (stdout,*) 'GLOBIDATE2 : ', tochar(globidate2)
      write (stdout,*) 'NSTEPS = ', nsteps

      call open_sstfile(idateo)

      allocate(sst%hint(1))

      idate = idateo
      do n = 1, nsteps
        call read_func(idate,sst,xlat,xlon)
        call h_interpolate_cont(sst%hint(1),sst%var,sstmm)
        call writerec(idate)
        write (stdout,*) 'WRITEN OUT SST DATA : ', tochar(idate)
        idate = idate + step
      end do

      call h_interpolator_destroy(sst%hint(1))
      deallocate(sst%hint)
    end subroutine cmip7_sst

end module mod_sst_cmip7

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

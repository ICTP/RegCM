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

module mod_sst_cmip6

  use mod_intkinds, only : ik8
  use mod_realkinds, only : rkx
  use mod_cmip6_helper, only : cmip6_2d_var
  use mod_message, only : die
  use mod_dynparam, only : cmip6_model, dattyp, globidate1, globidate2, ical, &
    pmip4_model, ssttyp
  use mod_sst_grid, only : open_sstfile, sstmm, writerec, xlat, xlon
  use mod_kdinterp, only : h_interpolate_cont, h_interpolator_destroy
  use mod_date, only : gregorian, imondiff, monfirst, nextmon, noleap, &
    prevmon, rcm_time_and_date, rcm_time_interval, tochar, tohours, umnt, &
    y360, assignment(=), operator(+), operator(-), operator(/=), operator(<), &
    operator(<=), operator(==), operator(>), operator(>=)
  use mod_stdio, only : stderr, stdout
  use mod_cmip6_cesm, only : read_sst_cesm
  use mod_cmip6_canesm, only : read_sst_canesm
  use mod_cmip6_cmcc, only : read_sst_cmcc
  use mod_cmip6_cnrm, only : read_sst_cnrm
  use mod_cmip6_ecea, only : read_sst_ecea
  use mod_cmip6_gfdl, only : read_sst_gfdl
  use mod_cmip6_hadmm, only : read_sst_hadmm
  use mod_cmip6_ipsllr, only : read_sst_ipsllr
  use mod_cmip6_miroc6, only : read_sst_miroc6
  use mod_cmip6_miresl, only : read_sst_miresl
  use mod_cmip6_mpihr, only : read_sst_mpihr
  use mod_cmip6_mpilr, only : read_sst_mpilr
  use mod_cmip6_normm, only : read_sst_normm
  use netcdf

  implicit none

  private

  public :: cmip6_sst

  type(cmip6_2d_var) :: sst

  abstract interface
    subroutine read_cmip6_sst(id,var,lat,lon)
      import
      implicit none
      type(rcm_time_and_date), intent(in) :: id
      type(cmip6_2d_var), intent(inout) :: var
      real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: lat, lon
    end subroutine read_cmip6_sst
  end interface

  contains

    subroutine cmip6_sst
      implicit none
      type(rcm_time_and_date) :: idate, idatef, idateo
      type(rcm_time_interval) :: tdif, step
      procedure(read_cmip6_sst), pointer :: read_func => null( )
      integer :: nsteps, n

      idateo = globidate1
      idatef = globidate2
      tdif = idatef-idateo

      dattyp = ssttyp

      if ( ssttyp == 'CMIP6' ) then
        select case (cmip6_model)
          case ('MPI-ESM1-2-HR')
            if ( ical /= gregorian ) then
              write(stderr,*) 'MPI-ESM1-2-HR requires gregorian calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_mpihr
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'HadGEM3-GC31-MM' )
            if ( ical /= y360 ) then
              write(stderr,*) 'HadGEM3-GC31-MM requires 360_day calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_hadmm
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'GFDL-ESM4' )
            if ( ical /= noleap ) then
              write(stderr,*) 'GFDL-ESM4 requires noleap calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_gfdl
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'NorESM2-MM' )
            if ( ical /= noleap ) then
              write(stderr,*) 'NorESM2-MM requires noleap calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_normm
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'CNRM-ESM2-1' )
            if ( ical /= gregorian ) then
              write(stderr,*) 'CNRM-ESM2-1 requires gregorian calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_cnrm
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'EC-Earth3-Veg' )
            if ( ical /= gregorian ) then
              write(stderr,*) 'EC-Earth3-Veg requires gregorian calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_ecea
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'CESM2' )
            if ( ical /= noleap ) then
              write(stderr,*) 'CESM2 requires noleap calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_cesm
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'CMCC-ESM2' )
            if ( ical /= noleap ) then
              write(stderr,*) 'CMCC-ESM2 requires noleap calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_cmcc
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'CanESM5' )
            if ( ical /= noleap ) then
              write(stderr,*) 'CanESM5 requires noleap calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_canesm
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'MIROC6' )
            if ( ical /= gregorian ) then
              write(stderr,*) 'MIROC6 requires gregorian calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_miroc6
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ( 'MIROC-ES2L' )
            if ( ical /= gregorian ) then
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
          case ('MPI-ESM1-2-LR')
            if ( ical /= gregorian ) then
              write(stderr,*) 'MPI-ESM1-2-LR requires gregorian calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_mpilr
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case default
            call die('sst','Unknown CMIP6 model: '//trim(cmip6_model),1)
        end select
      else if ( ssttyp == 'PMIP4' ) then
        select case (pmip4_model)
          case ('MPI-ESM1-2-LR')
            if ( ical /= gregorian ) then
              write(stderr,*) 'MPI-ESM1-2-LR requires gregorian calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_mpilr
            sst%vname = 'tos'
            step = 86400
            nsteps = int(tohours(tdif))/24 + 1
          case ('IPSL-CM6A-LR')
            if ( ical /= gregorian ) then
              write(stderr,*) 'IPSL-CM6A-LR requires gregorian calendar.'
              call die('sst','Calendar mismatch',1)
            end if
            read_func => read_sst_ipsllr
            sst%vname = 'tos'
            idateo = prevmon(globidate1)
            idatef = monfirst(globidate2)
            if (idatef < globidate2) then
              idatef = nextmon(idatef)
            end if
            step = rcm_time_interval(1_ik8,umnt)
            nsteps = imondiff(idatef,idateo) + 1
          case default
            call die('sst','Unknown PMIP4 model: '//trim(pmip4_model),1)
        end select
      else
        call die('sst','No PMIP4/CMIP6 model', 1)
      end if

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
    end subroutine cmip6_sst

end module mod_sst_cmip6

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

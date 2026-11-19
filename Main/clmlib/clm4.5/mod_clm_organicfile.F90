module mod_clm_organicfile
  use mod_realkinds, only : rk8
  use mod_stdio, only : stdout
  use mod_dynparam, only : myid
  use mod_mppparam, only : italk
  use mod_clm_nchelper, only : clm_closefile, clm_filetype, clm_openfile, &
    clm_readvar
  use mod_clm_domain, only : ldomain
  use mod_clm_varctl, only : fsurdat
  use mod_clm_decomp, only : gcomm_gridcell

  implicit none

  private

  save

  public :: organicrd  ! Read organic matter dataset

  contains
  !
  ! Read the organic matter dataset.
  !
  subroutine organicrd(organic)
    implicit none
    real(rk8), pointer, contiguous :: organic(:,:)  ! organic matter density (kg/m3)
    type(clm_filetype) :: ncid                  ! netcdf id
    character(len=32)  :: subname = 'organicrd' ! subroutine name

    ! Initialize data to zero - no organic matter dataset

    organic(:,:)   = 0._rk8

    ! Read data if file was specified in namelist

    if (fsurdat /= ' ') then
      if (myid == italk) then
        write(stdout,*) 'Attempting to read organic matter data .....'
        write(stdout,*) ' ',trim(subname),' : ',trim(fsurdat)
      end if

      call clm_openfile (fsurdat,ncid)

      call clm_readvar(ncid,'ORGANIC',organic,gcomm_gridcell)

      call clm_closefile (ncid)

      if ( myid == italk )then
        write(stdout,*) 'Successfully read organic matter data'
        write(stdout,*)
      end if
    end if
  end subroutine organicrd

end module mod_clm_organicfile
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

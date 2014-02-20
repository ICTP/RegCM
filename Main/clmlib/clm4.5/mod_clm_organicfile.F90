module mod_clm_organicfile
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_mppparam
  use mod_mpmessage
  use mod_clm_nchelper
  use mod_clm_domain , only : ldomain
  use mod_clm_type , only : grlnd
  use mod_clm_varctl , only : fsurdat , single_column

  implicit none

  private

  public :: organicrd  ! Read organic matter dataset

  contains
  !
  ! Read the organic matter dataset.
  !
  subroutine organicrd(organic)
    implicit none
    real(rk8), pointer :: organic(:,:)  ! organic matter density (kg/m3)
    character(len=256) :: locfn                 ! local file name
    type(clm_filetype) :: ncid                  ! netcdf id
    integer(ik4)       :: ni,nj,ns              ! dimension sizes  
    logical            :: isgrid2d              ! true => file is 2d
    logical            :: readvar               ! true => variable is on dataset
    character(len=32)  :: subname = 'organicrd' ! subroutine name

    ! Initialize data to zero - no organic matter dataset

    organic(:,:)   = 0.D0
       
    ! Read data if file was specified in namelist
       
    if (fsurdat /= ' ') then
      if (myid == italk) then
        write(stdout,*) 'Attempting to read organic matter data .....'
        write(stdout,*) subname,trim(fsurdat)
      end if

      call clm_openfile (fsurdat,ncid)

      call clm_check_dims(ncid, ni, nj)
      if ( nj == 1 ) isgrid2d = .false.
      ns = ni*nj
      if ( ldomain%ns /= ns .or. &
           ldomain%ni /= ni .or. &
           ldomain%nj /= nj ) then
        write(stderr,*) &
                trim(subname), 'ldomain and input file do not match dims '
        write(stderr,*) trim(subname), 'ldomain%ni,ni,= ',ldomain%ni,ni
        write(stderr,*) trim(subname), 'ldomain%nj,nj,= ',ldomain%nj,nj
        write(stderr,*) trim(subname), 'ldomain%ns,ns,= ',ldomain%ns,ns
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
       
      call clm_readvar(ncid,'ORGANIC',organic)

      if ( myid == italk )then
        write(stdout,*) 'Successfully read organic matter data'
        write(stdout,*)
      end if
    end if
  end subroutine organicrd

end module mod_clm_organicfile

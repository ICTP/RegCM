module mod_clm_restfile
  !
  ! Reads from or writes to/ the CLM restart file.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mpmessage
  use mod_dynparam
  use mod_mppparam
  use mod_clm_surfrd , only : crop_prog
  use clm_time_manager , only : timemgr_restart_io , get_nstep
  use mod_clm_subgridrest , only : SubgridRest
  use mod_clm_biogeophysrest , only : BiogeophysRest
  use mod_clm_accumul , only : accumulRest
  use mod_clm_slakerest , only : SLakeRest
  use mod_clm_decomp ,  only : get_proc_bounds , get_proc_global
  use mod_clm_varpar , only : nlevsno , nlevlak , nlevgrnd , nlevurb
  use mod_clm_varpar , only : numrad , nlevcan
  use mod_clm_varctl , only : single_column , nsrest , nsrStartup

  implicit none

  private

  public :: restFile_read
  public :: restFile_write
  public :: restFile_open
  public :: restFile_close
  public :: restFile_getfile
  public :: restFile_filename   ! Sets restart filename

  private :: restFile_read_pfile
  ! Writes restart pointer file
  private :: restFile_write_pfile
  ! Close restart file and write restart pointer file
  private :: restFile_closeRestart
  private :: restFile_dimset
  private :: restFile_dimcheck
  private :: restFile_enddef

  contains
  !
  ! Read/write CLM restart file.
  !
  subroutine restFile_write(file,nlend,noptr,rdate)
#if (defined CN)
    use CNRestMod        , only : CNRest
    use CropRestMod      , only : CropRest
#endif
#if (defined LCH4)
    use ch4RestMod       , only : ch4Rest
#endif
    use histFileMod      , only : hist_restart_ncd
    implicit none
    character(len=*) , intent(in) :: file  ! output netcdf restart file
    logical , intent(in) :: nlend          ! if at the end of the simulation
    character(len=*) , intent(in) :: rdate ! restart file time stamp for name
    ! if should NOT write to the restart pointer file
    logical , intent(in) , optional :: noptr
    type(clm_filetype) :: ncid ! netcdf id
    integer(ik4) :: i       ! index
    logical :: ptrfile ! write out the restart pointer file

    if ( present(noptr) ) then
      ptrfile = .not. noptr
    else
      ptrfile = .true.
    end if

    ! --------------------------------------------
    ! Open restart file
    ! --------------------------------------------

    call restFile_open(flag='write',file=file,ncid=ncid)

    ! --------------------------------------------
    ! Define dimensions and variables
    ! --------------------------------------------

    call restFile_dimset(ncid)

    ! Define restart file variables

    call timemgr_restart_io(ncid,flag='define')

    call SubgridRest( ncid, flag='define' )

    call BiogeophysRest( ncid, flag='define' )
#if (defined CN)
    call CNRest( ncid, flag='define' )
    if ( crop_prog ) call CropRest( ncid, flag='define' )
#endif

    call accumulRest( ncid, flag='define' )
    call SLakeRest( ncid, flag='define' )
#if (defined LCH4)
    call ch4Rest ( ncid, flag='define' )
#endif

    call hist_restart_ncd ( ncid, flag='define', rdate=rdate )

    call restFile_enddef( ncid )

    ! --------------------------------------------
    ! Write restart file variables
    ! --------------------------------------------

    call timemgr_restart_io( ncid, flag='write' )

    call SubgridRest( ncid, flag='write' )

    call BiogeophysRest( ncid, flag='write' )

#if (defined CN)
    call CNRest( ncid, flag='write' )
    if ( crop_prog ) call CropRest( ncid, flag='write' )
#endif

    call SLakeRest( ncid, flag='write' )
#if (defined LCH4)
    call ch4Rest ( ncid, flag='write' )
#endif
    call accumulRest( ncid, flag='write' )

    call hist_restart_ncd (ncid, flag='write' )

    ! --------------------------------------------
    ! Close restart file and write restart pointer file
    ! --------------------------------------------

    call restFile_close( ncid )
    call restFile_closeRestart( file, nlend )

    ! Write restart pointer file

    if ( ptrfile ) call restFile_write_pfile( file )

    ! Write out diagnostic info

    if (myid == italk) then
       write(stdout,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
       write(stdout,'(72a1)') ("-",i=1,60)
    end if

  end subroutine restFile_write

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read
!
! !INTERFACE:
  subroutine restFile_read( file )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
    use BiogeophysRestMod, only : BiogeophysRest
#if (defined CN)
    use CNRestMod        , only : CNRest
    use CropRestMod      , only : CropRest
#endif
    use SLakeRestMod  , only : SLakeRest
#if (defined LCH4)
    use ch4RestMod       , only : ch4Rest
#endif
    use accumulMod       , only : accumulRest
    use histFileMod      , only : hist_restart_ncd
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file  ! output netcdf restart file
!
! !CALLED FROM:
! subroutine initialize2
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    type(file_desc_t) :: ncid ! netcdf id
    integer(ik4) :: i              ! index
!-----------------------------------------------------------------------

    ! Open file

    call restFile_open( flag='read', file=file, ncid=ncid )

    ! Read file

    call restFile_dimcheck( ncid )

    call SLakeRest( ncid, flag='read' )

    call BiogeophysRest( ncid, flag='read' )

#if (defined CN)
    call CNRest( ncid, flag='read' )
    if ( crop_prog ) call CropRest( ncid, flag='read' )
#endif

#if (defined LCH4)
    call ch4Rest( ncid, flag='read' )
#endif

    call accumulRest( ncid, flag='read' )

    call hist_restart_ncd (ncid, flag='read')

    ! Close file

    call restFile_close( ncid )

    ! Write out diagnostic info

    if (myid == italk) then
       write(stdout,'(72a1)') ("-",i=1,60)
       write(stdout,*) 'Successfully read restart data for restart run'
       write(stdout,*)
    end if

  end subroutine restFile_read

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_getfile
!
! !INTERFACE:
  subroutine restFile_getfile( file, path )
!
! !DESCRIPTION:
! Determine and obtain netcdf restart file
!
! !USES:
    use clm_varctl, only : caseid, finidat, nrevsn, nsrest, brnch_retain_casename, &
                           nsrContinue, nsrBranch, nsrStartup
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file
!
! !CALLED FROM:
! subroutine initialize2
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer(ik4) :: status                      ! return status
    integer(ik4) :: length                      ! temporary
    character(len=256) :: ftest,ctest      ! temporaries
!-----------------------------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file

    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, file, 0 )
    end if

    ! Branch run:
    ! Restart file pathname is obtained from namelist "nrevsn"
    ! Check case name consistency (case name must be different for branch run,
    ! unless namelist specification states otherwise)

    if (nsrest==nsrBranch) then
       length = len_trim(nrevsn)
       if (nrevsn(length-2:length) == '.nc') then
          path = trim(nrevsn)
       else
          path = trim(nrevsn) // '.nc'
       end if
       call getfil( path, file, 0 )

       ! tcraig, adding xx. and .clm2 makes this more robust
       ctest = 'xx.'//trim(caseid)//'.clm2'
       ftest = 'xx.'//trim(file)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
          write(stderr,*) 'Must change case name on branch run if ',&
               'brnch_retain_casename namelist is not set'
          write(stderr,*) 'previous case filename= ',trim(file),&
               ' current case = ',trim(caseid), ' ctest = ',trim(ctest), &
               ' ftest = ',trim(ftest)
          call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! Initial run:
    ! Restart file pathname is obtained from namelist "finidat"

    if (nsrest==nsrStartup) then
       call getfil( finidat, file, 0 )
    end if

  end subroutine restFile_getfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read_pfile
!
! !INTERFACE:
  subroutine restFile_read_pfile( pnamer )
!
! !DESCRIPTION:
! Setup restart file and perform necessary consistency checks
!
! !USES:
    use fileutils , only : opnfil, getavu, relavu
    use clm_varctl, only : rpntfil, rpntdir, inst_suffix
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: pnamer ! full path of restart file
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer(ik4) :: i                  ! indices
    integer(ik4) :: nio                ! restart unit
    integer(ik4) :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
!-----------------------------------------------------------------------

    ! Obtain the restart file from the restart pointer file
    ! For restart runs, the restart pointer file contains the full pathname
    ! of the restart file. For branch runs, the namelist variable
    ! [nrevsn] contains the full pathname of the restart file.
    ! New history files are always created for branch runs.

    if (myid == italk) then
       write(stdout,*) 'Reading restart pointer file....'
    endif

    nio = getavu()
    locfn = trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
    call opnfil (locfn, nio, 'f')
    read (nio,'(a256)') pnamer
    call relavu (nio)

    if (myid == italk) then
       write(stdout,*) 'Reading restart data.....'
       write(stdout,'(72a1)') ("-",i=1,60)
    end if

  end subroutine restFile_read_pfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_closeRestart
!
! !INTERFACE:
  subroutine restFile_closeRestart( file, nlend )
!
! !DESCRIPTION:
! Close restart file and write restart pointer file if
! in write mode, otherwise just close restart file if in read mode
!
! !USES:
    use clm_time_manager, only : is_last_step
!
! !ARGUMENTS:
    implicit none
    character(len=*) , intent(in) :: file  ! local output filename
    logical,           intent(in) :: nlend
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer(ik4) :: i                   !index
!-----------------------------------------------------------------------

   if (myid == italk) then
      write(stdout,*) 'Successfully wrote local restart file ',trim(file)
      write(stdout,'(72a1)') ("-",i=1,60)
      write(stdout,*)
   end if

 end subroutine restFile_closeRestart

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_write_pfile
!
! !INTERFACE:
  subroutine restFile_write_pfile( fnamer )
!
! !DESCRIPTION:
! Open restart pointer file. Write names of current netcdf restart file.
!
! !USES:
    use clm_varctl, only : rpntdir, rpntfil, inst_suffix
    use fileutils , only : relavu
    use fileutils , only : getavu, opnfil
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fnamer
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer(ik4) :: m                    ! index
    integer(ik4) :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
!-----------------------------------------------------------------------

    if (myid == italk) then
       nio = getavu()
       filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
       call opnfil( filename, nio, 'f' )

       write(nio,'(a)') fnamer
       call relavu( nio )
       write(stdout,*)'Successfully wrote local restart pointer file'
    end if

  end subroutine restFile_write_pfile

!-----------------------------------------------------------------------
  subroutine restFile_open( flag, file, ncid )

    use clm_time_manager, only : get_nstep

    implicit none
    character(len=*),  intent(in) :: flag ! flag to specify read or write
    character(len=*),  intent(in) :: file ! filename
    type(file_desc_t), intent(out):: ncid ! netcdf id

    integer(ik4) :: omode                              ! netCDF dummy variable
    character(len= 32) :: subname='restFile_open' ! subroutine name

    if (flag == 'write') then

       ! Create new netCDF file (in define mode) and set fill mode
       ! to "no fill" to optimize performance

       if (myid == italk) then
          write(stdout,*)
          write(stdout,*)'restFile_open: writing restart dataset at ',&
               trim(file), ' at nstep = ',get_nstep()
          write(stdout,*)
       end if
       call ncd_pio_createfile(ncid, trim(file))

    else if (flag == 'read') then

       ! Open netcdf restart file

       if (myid == italk) then
          write(stdout,*) 'Reading restart dataset'
       end if
       call ncd_pio_openfile (ncid, trim(file), 0)

    end if

  end subroutine restFile_open

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_filename
!
! !INTERFACE:
  character(len=256) function restFile_filename( rdate )
!
! !DESCRIPTION:
!
! !USES:
    use clm_varctl, only : caseid, inst_suffix
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: rdate   ! input date for restart file name
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

    restFile_filename = "./"//trim(caseid)//".clm2"//trim(inst_suffix)//&
                        ".r."//trim(rdate)//".nc"
    if (myid == italk) then
       write(stdout,*)'writing restart file ',trim(restFile_filename),' for model date = ',rdate
    end if

  end function restFile_filename

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_dimset
!
! !INTERFACE:
  subroutine restFile_dimset( ncid )
!
! !DESCRIPTION:
! Read/Write initial data from/to netCDF instantaneous initial data file
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kindD0
    use clm_time_manager, only : get_nstep, get_curr_date
    use spmdMod     , only : mpicom, MPI_LOGICAL
    use clm_varctl  , only : caseid, ctitle, version, username, hostname, fsurdat, &
                             conventions, source
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer(ik4) :: yr                  ! current year (0 -> ...)
    integer(ik4) :: mon                 ! current month (1 -> 12)
    integer(ik4) :: day                 ! current day (1 -> 31)
    integer(ik4) :: mcsec               ! seconds of current date
    integer(ik4) :: mcdate              ! current date
    integer(ik4) :: dimid               ! netCDF dimension id
    integer(ik4) :: numg                ! total number of gridcells across all processors
    integer(ik4) :: numl                ! total number of landunits across all processors
    integer(ik4) :: numc                ! total number of columns across all processors
    integer(ik4) :: nump                ! total number of pfts across all processors
    integer(ik4) :: ier                 ! error status
    integer(ik4) :: strlen_dimid        ! string dimension id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: str
    character(len= 32) :: subname='restFile_dimset' ! subroutine name
!------------------------------------------------------------------------

    call get_proc_global(numg, numl, numc, nump)

    ! Define dimensions

    call ncd_defdim(ncid, 'gridcell', numg           , dimid)
    call ncd_defdim(ncid, 'landunit', numl           , dimid)
    call ncd_defdim(ncid, 'column'  , numc           , dimid)
    call ncd_defdim(ncid, 'pft'     , nump           , dimid)

    call ncd_defdim(ncid, 'levgrnd' , nlevgrnd       , dimid)
    call ncd_defdim(ncid, 'levurb'  , nlevurb        , dimid)
    call ncd_defdim(ncid, 'levlak'  , nlevlak        , dimid)
    call ncd_defdim(ncid, 'levsno'  , nlevsno        , dimid)
    call ncd_defdim(ncid, 'levsno1'  , nlevsno+1     , dimid)
    call ncd_defdim(ncid, 'levtot'  , nlevsno+nlevgrnd, dimid)
    call ncd_defdim(ncid, 'numrad'  , numrad         , dimid)
    call ncd_defdim(ncid, 'levcan'  , nlevcan        , dimid)
    call ncd_defdim(ncid, 'string_length', 64        , dimid)

    ! Define global attributes

    call ncd_putatt(ncid, NCD_GLOBAL, 'Conventions', trim(conventions))
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call ncd_putatt(ncid, NCD_GLOBAL, 'history' , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'username', trim(username))
    call ncd_putatt(ncid, NCD_GLOBAL, 'host'    , trim(hostname))
    call ncd_putatt(ncid, NCD_GLOBAL, 'version' , trim(version))
    call ncd_putatt(ncid, NCD_GLOBAL, 'source'  , trim(source))
    str = '$Id: restFileMod.F90 41292 2012-10-26 13:51:45Z erik $'
    call ncd_putatt(ncid, NCD_GLOBAL, 'revision_id'    , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(caseid))
    call ncd_putatt(ncid, NCD_GLOBAL, 'surface_dataset', trim(fsurdat))
    call ncd_putatt(ncid, NCD_GLOBAL, 'title', &
          'CLM Restart information, required to continue a simulation' )
  end subroutine restFile_dimset
  !
  ! Check dimensions of restart file
  !
  subroutine restFile_dimcheck( ncid )
    implicit none
    type(clm_filetype), intent(inout) :: ncid
    integer(ik4) :: numg  ! total number of gridcells across all processors
    integer(ik4) :: numl  ! total number of landunits across all processors
    integer(ik4) :: numc  ! total number of columns across all processors
    integer(ik4) :: nump  ! total number of pfts across all processors

    ! Get relevant sizes

    if ( .not. single_column .or. nsrest /= nsrStartup ) then
      call get_proc_global(numg, numl, numc, nump)
      if ( .not. clm_check_dimlen(ncid, 'gridcell', numg) ) &
        call fatal(__FILE__,__LINE__,'NUM GRIDCELL DIFFER !')
      if ( .not. clm_check_dimlen(ncid, 'landunit', numl) ) &
        call fatal(__FILE__,__LINE__,'NUM LANDUNIT DIFFER !')
      if ( .not. clm_check_dimlen(ncid, 'column', numc) ) &
        call fatal(__FILE__,__LINE__,'NUM COLUMN DIFFER !')
      if ( .not. clm_check_dimlen(ncid, 'pft', nump) ) &
        call fatal(__FILE__,__LINE__,'NUM PFTS DIFFER !')
    end if
    if ( .not. clm_check_dimlen(ncid, 'levsno', nlevsno) ) &
      call fatal(__FILE__,__LINE__,'NUM NLEVSNO DIFFER !')
    if ( .not. clm_check_dimlen(ncid, 'levgrnd', nlevgrnd) ) &
      call fatal(__FILE__,__LINE__,'NUM LEVGRND DIFFER !')
    if ( .not. clm_check_dimlen(ncid, 'levurb', nlevurb) ) &
      call fatal(__FILE__,__LINE__,'NUM LEVURB DIFFER !')
    if ( .not. clm_check_dimlen(ncid, 'levlak', nlevlak) ) &
      call fatal(__FILE__,__LINE__,'NUM LEVLAK DIFFER !')
  end subroutine restFile_dimcheck
  !
  ! Read a CLM restart file.
  !
  subroutine restFile_enddef( ncid )
    implicit none
    type(clm_filetype), intent(inout) :: ncid
    call clm_enddef(ncid)
  end subroutine restFile_enddef
  !
  ! Read a CLM restart file.
  !
  subroutine restFile_close( ncid )
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    call clm_closefile(ncid)
  end subroutine restFile_close

end module mod_clm_restfile

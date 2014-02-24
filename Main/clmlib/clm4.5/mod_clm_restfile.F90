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
  use mod_runparams
  use mod_clm_nchelper
  use mod_clm_surfrd , only : crop_prog
  use mod_clm_subgridrest , only : SubgridRest
  use mod_clm_biogeophysrest , only : BiogeophysRest
  use mod_clm_accumul , only : accumulRest
  use mod_clm_slakerest , only : SLakeRest
  use mod_clm_decomp ,  only : get_proc_bounds , get_proc_global
  use mod_clm_varpar , only : nlevsno , nlevlak , nlevgrnd , nlevurb
  use mod_clm_varpar , only : numrad , nlevcan
  use mod_clm_varctl , only : caseid , ctitle , version , username ,    &
          hostname , finidat , fsurdat , single_column , nsrest ,       &
          nrevsn , nsrStartup , nsrContinue , nsrBranch , inst_suffix , &
          brnch_retain_casename , conventions , source
  use mod_clm_varctl , only : rpntfil, rpntdir, inst_suffix
#if (defined CN)
  use mod_clm_cnrest , only : CNRest
  use mod_clm_croprest , only : CropRest
#endif
#if (defined LCH4)
  use mod_clm_ch4rest , only : ch4Rest
#endif
  use mod_clm_histfile , only : hist_restart_ncd
  use mod_clm_time_manager

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
  subroutine restFile_write(rfile,nlend,noptr,rdate)
    implicit none
    character(len=*) , intent(in) :: rfile  ! output netcdf restart file
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

    call restFile_open(flag='write',rfile=rfile,ncid=ncid)

    ! --------------------------------------------
    ! Define dimensions and variables
    ! --------------------------------------------

    call restFile_dimset(ncid)

    ! Define restart file variables

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
    call restFile_closeRestart( rfile, nlend )

    ! Write restart pointer file

    if ( ptrfile ) call restFile_write_pfile( rfile )

    ! Write out diagnostic info

    if (myid == italk) then
      write(stdout,*) &
              'Successfully wrote out restart data at nstep = ',ktau
      write(stdout,'(72a1)') ("-",i=1,60)
    end if
  end subroutine restFile_write
  !
  ! Read a CLM restart file.
  !
  subroutine restFile_read( rfile )
    implicit none
    character(len=*), intent(in) :: rfile  ! output netcdf restart file
    type(clm_filetype) :: ncid ! netcdf id
    integer(ik4) :: i          ! index

    ! Open file

    call restFile_open( flag='read', rfile=rfile, ncid=ncid )

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
  !
  ! Determine and obtain netcdf restart file
  !
  subroutine restFile_getfile( rfile, path )
    implicit none
    character(len=*), intent(out) :: rfile  ! name of netcdf restart file
    ! full pathname of netcdf restart file
    character(len=*), intent(out) :: path
    integer(ik4) :: status                      ! return status
    integer(ik4) :: length                      ! temporary
    character(len=256) :: ftest,ctest      ! temporaries

    ! Continue run:
    ! Restart file pathname is read restart pointer file

    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
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

      ! tcraig, adding xx. and .clm2 makes this more robust
      ctest = 'xx.'//trim(caseid)//'.clm2'
      ftest = 'xx.'//trim(rfile)
      status = index(trim(ftest),trim(ctest))
      if (status /= 0 .and. .not.(brnch_retain_casename)) then
        write(stderr,*) 'Must change case name on branch run if ',&
               'brnch_retain_casename namelist is not set'
        write(stderr,*) 'previous case filename= ',trim(rfile),&
               ' current case = ',trim(caseid), ' ctest = ',trim(ctest), &
               ' ftest = ',trim(ftest)
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end if

  end subroutine restFile_getfile
  !
  ! Setup restart file and perform necessary consistency checks
  !
  subroutine restFile_read_pfile( pnamer )
    implicit none
    character(len=*), intent(out) :: pnamer ! full path of restart file
    integer(ik4) :: i      ! indices
    integer(ik4) :: nio    ! restart unit
    integer(ik4) :: ier    ! substring check status
    character(len=256) :: locfn ! Restart pointer file name
    ! Obtain the restart file from the restart pointer file
    ! For restart runs, the restart pointer file contains the full pathname
    ! of the restart file. For branch runs, the namelist variable
    ! [nrevsn] contains the full pathname of the restart file.
    ! New history files are always created for branch runs.
    if (myid == italk) then
       write(stdout,*) 'Reading restart pointer file....'
    endif
    nio = file_getunit( )
    locfn = trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
    open(unit=nio,file=locfn,form='formatted', &
         status='old',action='read',iostat=ier)
    if ( ier /= 0 ) then
      write(stderr,*) 'Cannot open pointer file ',trim(locfn)
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    read (nio,'(a256)') pnamer
    call file_freeunit(nio)
    if (myid == italk) then
      write(stdout,*) 'Reading restart data.....'
      write(stdout,'(72a1)') ("-",i=1,60)
    end if
  end subroutine restFile_read_pfile
  !
  ! Close restart file and write restart pointer file if
  ! in write mode, otherwise just close restart file if in read mode
  !
  subroutine restFile_closeRestart( rfile, nlend )
    implicit none
    character(len=*) , intent(in) :: rfile  ! local output filename
    logical,           intent(in) :: nlend
    integer(ik4) :: i                   !index
    if (myid == italk) then
      write(stdout,*) 'Successfully wrote local restart file ',trim(rfile)
      write(stdout,'(72a1)') ("-",i=1,60)
      write(stdout,*)
    end if
  end subroutine restFile_closeRestart
  !
  ! Open restart pointer file. Write names of current netcdf restart file.
  !
  subroutine restFile_write_pfile( fnamer )
    implicit none
    character(len=*), intent(in) :: fnamer
    integer(ik4) :: m , ier              ! index
    integer(ik4) :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
    if (myid == italk) then
      nio = file_getunit( )
      filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
      open(unit=nio,file=filename,form='formatted', &
           status='replace',action='write',iostat=ier)
      if ( ier /= 0 ) then
        write(stderr,*) 'Cannot open pointer file ',trim(filename)
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
      write(nio,'(a)') fnamer
      call file_freeunit(nio)
      write(stdout,*)'Successfully wrote local restart pointer file'
    end if
  end subroutine restFile_write_pfile

  subroutine restFile_open( flag, rfile, ncid )
    implicit none
    character(len=*),  intent(in) :: flag ! flag to specify read or write
    character(len=*),  intent(in) :: rfile ! filename
    type(clm_filetype), intent(out):: ncid ! netcdf id
    integer(ik4) :: omode                              ! netCDF dummy variable
    if (flag == 'write') then
      ! Create new netCDF file (in define mode) and set fill mode
      ! to "no fill" to optimize performance
      if (myid == italk) then
        write(stdout,*)
        write(stdout,*)'Writing restart dataset at ', trim(rfile), &
                ' at ktau = ', ktau
        write(stdout,*)
      end if
      call clm_createfile(trim(rfile),ncid)
    else if (flag == 'read') then
      ! Open netcdf restart file
      if (myid == italk) then
        write(stdout,*) 'Reading restart dataset at ', trim(rfile)
      end if
      call clm_openfile(trim(rfile),ncid)
    end if
  end subroutine restFile_open

  character(len=256) function restFile_filename( rdate )
    implicit none
    character(len=*) , intent(in) :: rdate ! input date for restart file name
    restFile_filename = "./"//trim(caseid)//".clm2"//trim(inst_suffix)//&
                        ".r."//trim(rdate)//".nc"
    if (myid == italk) then
      write(stdout,*) &
       'writing restart file ',trim(restFile_filename), &
       ' for model date = ',rdate
    end if
  end function restFile_filename
  !
  ! Read/Write initial data from/to netCDF instantaneous initial data file
  !
  subroutine restFile_dimset( ncid )
    implicit none
    type(clm_filetype), intent(inout) :: ncid
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

    call get_proc_global(numg, numl, numc, nump)

    ! Define dimensions

    call clm_adddim(ncid,'gridcell',numg)
    call clm_adddim(ncid,'landunit',numl)
    call clm_adddim(ncid,'column',numc)
    call clm_adddim(ncid,'pft',nump)

    call clm_adddim(ncid,'levgrnd',nlevgrnd)
    call clm_adddim(ncid,'levurb',nlevurb)
    call clm_adddim(ncid,'levlak',nlevlak)
    call clm_adddim(ncid,'levsno',nlevsno)
    call clm_adddim(ncid,'levsno1',nlevsno+1)
    call clm_adddim(ncid,'levtot',nlevsno+nlevgrnd)
    call clm_adddim(ncid,'numrad',numrad)
    call clm_adddim(ncid,'levcan',nlevcan)
    call clm_adddim(ncid,'string_length', 64)

    ! Define global attributes

    call clm_addatt(ncid,'Conventions',trim(conventions))
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call clm_addatt(ncid, 'history',trim(str))
    call clm_addatt(ncid, 'username',trim(username))
    call clm_addatt(ncid, 'host',trim(hostname))
    call clm_addatt(ncid, 'version',trim(version))
    call clm_addatt(ncid, 'source',trim(source))
    call clm_addatt(ncid, 'case_title',trim(ctitle))
    call clm_addatt(ncid, 'case_id', trim(caseid))
    call clm_addatt(ncid, 'surface_dataset',trim(fsurdat))
    call clm_addatt(ncid, 'title', &
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

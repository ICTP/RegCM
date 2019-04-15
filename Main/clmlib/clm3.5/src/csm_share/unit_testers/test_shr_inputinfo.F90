program test_shr_inputinfo
!
! Simple unit test of shared input module:
!
 use shr_kind_mod,  only: shr_kind_cl, shr_kind_cs
 use shr_inputInfo_mod, only: &
              shr_inputInfo_initRead,          shr_inputInfo_initSetDefaults,  &
              shr_inputInfo_initGetData,       shr_inputInfo_initPutData,      &
              shr_inputInfo_initIsRestart,     shr_inputInfo_initRestWrite,    &
              shr_inputInfo_initRestRead,      shr_inputInfo_initType,         &
              shr_inputInfo_initPrint,         shr_inputInfo_initRPointerRead, &
              shr_inputInfo_initRPointerWrite, shr_inputInfo_initIsBranch,     &
              shr_inputInfo_initIsContinue,    shr_inputInfo_initIsModelToRun, &
              shr_inputInfo_initRestArchive,   shr_inputInfo_initIsStartup
  use shr_ncio_mod, only: shr_ncio_setDebug
  use shr_sys_mod,  only: shr_sys_abort
  use shr_mpi_mod,  only: shr_mpi_init, shr_mpi_finalize, shr_mpi_commrank

  implicit none
  type(shr_inputInfo_initType) :: CCSMInit
  character(len=SHR_KIND_CL) :: NLFilename
  character(len=SHR_KIND_CL) :: restart_file
  character(len=SHR_KIND_CL) :: restart_file2
  integer, parameter :: ymd  = 20051225, tod = 1200
  integer, parameter :: ymd2 =      101, tod2=85700
  integer :: perpetual_ymd
  logical :: perpetual
  integer :: mpicom = 0
  logical :: mastertask = .true.
  integer :: mpirank
  character(len=3), parameter :: models(4) = (/"atm", "ice", "lnd", "ocn"/)
  integer :: itemp, ivalue, mss_irt, mss_irt2
  character(len=SHR_KIND_CL) :: string1, string2, start_type, case_desc, case_name, mss_wpass, &
                                archive_dir, restart_pfile
  character(len=SHR_KIND_CS) :: logFilePostFix, outPathRoot
  character(len=SHR_KIND_CL) :: case_desc2, case_name2, mss_wpass2, &
                                archive_dir2, restart_pfile2, mode, start_type2
  character(len=SHR_KIND_CS) :: logFilePostFix2, outPathRoot2
  logical :: atm_adiabatic, atm_ideal_phys, aqua_planet, exists, brnch_retain_casename
#include <mpif.h>

  call shr_mpi_init( )
  mpicom = MPI_COMM_WORLD
  call shr_mpi_commrank( mpicom, mpirank )
  mastertask = (mpirank == 0)
  call shr_ncio_setDebug( 5 )
  if ( mastertask )print *, 'Set to default values'
  NLFilename = "namelist"
  if ( mastertask )print *, 'Read in namelist from: ', trim(NLFilename)
  call shr_inputInfo_initSetDefaults( CCSMInit )
  if ( mastertask )print *, 'Write out default values of CCSMInit: '
  if ( mastertask )call shr_inputInfo_initPrint( CCSMInit )
  if ( mastertask )print *, 'Read namelist: '
  mss_wpass = "test_password"
  start_type="startup"
  case_name="csmrun"
  logFilePostFix = ".cdate.log"
  outpathRoot    = "/ccsm/run/case/logfiles/"
  call shr_inputInfo_initPutData( CCSMInit, start_type=start_type,  &
                                  case_name=case_name,              &
                                  mss_wpass=mss_wpass, logFilePostFix=logFilePostFix, &
                                  outpathRoot=outpathRoot )
  call shr_inputInfo_initGetData( CCSMInit, start_type=start_type2, &
                                  case_name=case_name2, &
                                  mss_wpass=mss_wpass2, logFilePostFix=logFilePostFix2, &
                                  outpathRoot=outpathRoot2 )
  if ( trim(logFilePostFix2) /= trim(logFilePostFix) )then
     call shr_sys_abort( 'Error setting and then getting logFilePostFix' )
  end if
  if ( trim(outPathRoot2) /= trim(outPathRoot) )then
     call shr_sys_abort( 'Error setting and then getting outPathRoot' )
  end if
  if ( trim(mss_wpass2) /= trim(mss_wpass) )then
     call shr_sys_abort( 'Error setting and then getting mss_wpass' )
  end if
  if ( trim(start_type2) /= trim(start_type) )then
     call shr_sys_abort( 'Error setting and then getting start_type' )
  end if
  if ( trim(case_name2) /= trim(case_name) )then
     call shr_sys_abort( 'Error setting and then getting case_name' )
  end if
  call shr_inputInfo_initRead( NLFilename, LogPrint=mastertask, mpicom=mpicom, &
                             MasterTask=mastertask, CCSMInitOut=CCSMInit )
  call shr_inputInfo_initGetData( CCSMInit, case_desc=case_desc, mss_irt=mss_irt, &
                          mss_wpass=mss_wpass, archive_dir=archive_dir,   &
                          restart_pfile=restart_pfile )
  if ( .not. shr_inputInfo_initIsStartup( CCSMInit ) ) then
     call shr_sys_abort( 'This should be a startup case' )
  end if
  if ( shr_inputInfo_initIsRestart( CCSMInit ) ) then
     call shr_sys_abort( 'This should NOT be a restart case' )
  end if
  if ( shr_inputInfo_initIsContinue( CCSMInit ) ) then
     call shr_sys_abort( 'This should NOT be a continue case' )
  end if
  if ( shr_inputInfo_initIsBranch( CCSMInit ) ) then
     call shr_sys_abort( 'This should NOT be a branch case' )
  end if
  call shr_inputInfo_initGetData( CCSMInit, perpetual_run=perpetual, perpetual_ymd=perpetual_ymd, &
                      restart_file=restart_file )
  if ( len_trim(restart_file) /= 0 )then
     call shr_sys_abort( 'Restart file is set on initialization and should NOT be' )
  end if
  if ( mastertask )print *, 'perpetual, perpetual_ymd: ', perpetual, perpetual_ymd
  if ( mastertask )print *, 'Write restart file: '
  call shr_inputInfo_initRPointerWrite( ymd, tod, MPICom, MasterTask, &
                                        CCSMInit=CCSMInit, restart_file=restart_file )
  if ( len_trim(restart_file) == 0 )then
     call shr_sys_abort( 'After WriteRPointer, restart_file should be set and is not' )
  end if
  call shr_inputInfo_initRestWrite( restart_file, mpicom=mpicom, MasterTask=mastertask, &
                                    CCSMInit=CCSMInit )
  if ( mastertask )print *, 'Write restart file again (with different filetype) to make sure we can: '
  call shr_inputInfo_initRPointerWrite( ymd, tod, mpicom, MasterTask, &
                                        FileType="d.r",               &
                                        CCSMInit=CCSMInit, restart_file=restart_file )
  call shr_inputInfo_initRestWrite( restart_file, mpicom=mpicom, MasterTask=mastertask, &
                                   CCSMInit=CCSMInit )
  if ( mastertask )print *, 'Read restart file for continue case: '
  call shr_inputInfo_initSetDefaults( CCSMInit )
  call shr_inputInfo_initPutData( CCSMInit, start_type="continue", case_name="csmrun" )
  call shr_inputInfo_initRead( NLFilename, LogPrint=mastertask, mpicom=mpicom, &
                             MasterTask=mastertask, CCSMInitOut=CCSMInit )
  call shr_inputInfo_initRPointerRead( mpicom=mpicom, MasterTask=mastertask, &
                                   CCSMInit=CCSMInit, restart_file=restart_file )
  call shr_inputInfo_initRestRead( restart_file, mpicom=mpicom, MasterTask=mastertask, &
                                  CCSMInitOut=CCSMinit )
  call shr_inputInfo_initGetData( CCSMInit, case_desc=case_desc2, mss_irt=mss_irt2, &
                          mss_wpass=mss_wpass2, archive_dir=archive_dir2,   &
                          restart_pfile=restart_pfile2 )
  if ( case_desc /= case_desc2 )then
     call shr_sys_abort( 'case_desc Values from namelist not same as after write and read of restart' )
  end if
  if ( mss_irt /= mss_irt)then
     call shr_sys_abort( 'mss_irt Values from namelist not same as after write and read of restart' )
  end if
  if ( mss_wpass /= mss_wpass2 )then
     call shr_sys_abort( 'mss_wpass Values from namelist not same as after write and read of restart' )
  end if
  if ( archive_dir /= archive_dir2 )then
     call shr_sys_abort( 'archive_dir Values from namelist not same as after write and read of restart' )
  end if
  if ( restart_pfile /= restart_pfile2 )then
     call shr_sys_abort( 'restart_pfile Values from namelist not same as after write and read of restart' )
  end if
  if ( shr_inputInfo_initIsStartup( CCSMInit ) ) then
     call shr_sys_abort( 'This should NOT be a startup case' )
  end if
  if ( .not. shr_inputInfo_initIsRestart( CCSMInit ) ) then
     call shr_sys_abort( 'This should be a restart case' )
  end if
  if ( .not. shr_inputInfo_initIsContinue( CCSMInit ) ) then
     call shr_sys_abort( 'This should be a continue case' )
  end if
  if ( shr_inputInfo_initIsBranch( CCSMInit ) ) then
     call shr_sys_abort( 'This should NOT be a branch case' )
  end if
  if ( mastertask )call shr_inputInfo_initPrint( CCSMInit )
  if ( mastertask )print *, 'Read restart file for branch case: '
  call shr_inputInfo_initSetDefaults( CCSMInit )
  itemp = 29
  string1 = "case description changed to this"
  call shr_inputInfo_initPutData( CCSMInit, start_type="branch", case_name="csmrun_branch", &
                              mss_irt=itemp, restart_file=restart_file,                &
                              restart_file_override="mss_irt:case_desc" )
  call shr_inputInfo_initRead( NLFilename, LogPrint=mastertask, mpicom=mpicom, &
                             MasterTask=mastertask, CCSMInitOut=CCSMInit )
  call shr_inputInfo_initPutData( CCSMInit, mss_irt=itemp, case_desc=string1, &
                                  restart_file=restart_file )
  call shr_inputInfo_initRPointerRead( mpicom, mastertask, CCSMInit, restart_file )
  if ( mastertask )print *, 'restart_file: ', trim(restart_file)
  call shr_inputInfo_initRestRead( restart_file, mpicom=mpicom, MasterTask=mastertask, &
                                  CCSMInitOut=CCSMInit )
  call shr_inputInfo_initPutData( CCSMInit, atm_adiabatic=.false., &
                              atm_ideal_phys= .false., aqua_planet=.true., &
                              archive_dir="mss:./" )
  if ( shr_inputInfo_initIsStartup( CCSMInit ) ) then
     call shr_sys_abort( 'This should NOT be a startup case' )
  end if
  if ( .not. shr_inputInfo_initIsRestart( CCSMInit ) ) then
     call shr_sys_abort( 'This should be a restart case' )
  end if
  if ( shr_inputInfo_initIsContinue( CCSMInit ) ) then
     call shr_sys_abort( 'This should NOT be a continue case' )
  end if
  if ( .not. shr_inputInfo_initIsBranch( CCSMInit ) ) then
     call shr_sys_abort( 'This should be a branch case' )
  end if
  if ( mastertask )call shr_inputInfo_initPrint( CCSMInit )
  call shr_inputInfo_initGetData( CCSMInit, perpetual_run=perpetual, perpetual_ymd=perpetual_ymd, &
         mss_irt=ivalue, case_desc=string2, start_type=start_type )
  call shr_inputInfo_initGetData( CCSMInit, atm_adiabatic=atm_adiabatic,   &
                                    atm_ideal_phys=atm_ideal_phys, &
                                    aqua_planet=aqua_planet,       &
                                    brnch_retain_casename=brnch_retain_casename )
  if ( atm_adiabatic )then
     call shr_sys_abort( 'atm_adiabatic was set and should not have been' )
  end if
  if ( atm_ideal_phys )then
     call shr_sys_abort( 'atm_ideal_phys was set and should not have been' )
  end if
  if ( .not. aqua_planet )then
     call shr_sys_abort( 'aqua_planet was NOT set and should have been' )
  end if
  if ( trim(start_type) /= "branch" )then
     call shr_sys_abort( 'start_type was not set correctly' )
  end if
  if ( ivalue /= itemp )then
     call shr_sys_abort( 'mss_irt was not set correctly' )
  end if
  if ( string1 /= string2 )then
     call shr_sys_abort( 'case_desc was not set correctly' )
  end if
  if ( .not. shr_inputInfo_initIsModelToRun( CCSMInit, "atm") )then
     call shr_sys_abort( 'In Aqua-planet mode atm should run' )
  end if
  if ( shr_inputInfo_initIsModelToRun( CCSMInit, "ice") )then
     call shr_sys_abort( 'In Aqua-planet mode ice should NOT run' )
  end if
  if ( shr_inputInfo_initIsModelToRun( CCSMInit, "lnd") )then
     call shr_sys_abort( 'In Aqua-planet mode lnd should NOT run' )
  end if
  if ( .not. shr_inputInfo_initIsModelToRun( CCSMInit, "ocn") )then
     call shr_sys_abort( 'In Aqua-planet mode ocn should run' )
  end if
  if ( .not. perpetual .or. perpetual_ymd /= 321 )then
     call shr_sys_abort( 'Aqua-planet mode should be in perpetual mode' )
  end if
  if ( mastertask )print *, 'perpetual, perpetual_ymd: ', perpetual, perpetual_ymd
  if ( mastertask )print *, 'Write restart file: '
  call shr_inputInfo_initRPointerWrite( ymd2, tod2, MPICom=mpicom, &
                                    MasterTask=mastertask, CCSMInit=CCSMInit, &
                                    restart_file=restart_file2 )
  call shr_inputInfo_initRestWrite( restart_file2, mpicom=mpicom, MasterTask=mastertask, &
                                   CCSMInit=CCSMInit )
  if ( mastertask )print *, 'Restart file: ', trim(restart_file2)
  if ( len_trim(restart_file2) == 0 )then
     call shr_sys_abort( 'After WriteRPointer, restart_file should be set and is not' )
  end if
  if ( trim(restart_file) == trim(restart_file2) )then
     call shr_sys_abort( 'restart file for should be different now and is not' )
  end if
  restart_file = trim(restart_file2)
  if ( MasterTask ) write(6,*) 'archive restart_file = ', restart_file
  call shr_inputInfo_initRestArchive( MasterTask, CCSMInit, async=.false., remove=.true., &
                                      restart_file=restart_file )
  if ( len_trim(restart_file) /= 0 )then
     call shr_sys_abort( 'Local restart_file name should be set to null and is not' )
  end if
  inquire(file=trim(restart_file2), exist=exists )
  if ( exists )then
     call shr_sys_abort( 'local restart file should be removed and is not' )
  end if
  ! atm_adiabatic
  mode = "atm_adiabatic"
  call shr_inputInfo_initPutData( CCSMInit, atm_adiabatic=.true., &
                              atm_ideal_phys= .false., aqua_planet=.false. )
  if ( mastertask )call shr_inputInfo_initPrint( CCSMInit )
  call shr_inputInfo_initGetData( CCSMInit, perpetual_run=perpetual, perpetual_ymd=perpetual_ymd )
  call shr_inputInfo_initGetData( CCSMInit, atm_adiabatic=atm_adiabatic,   &
                                    atm_ideal_phys=atm_ideal_phys, &
                                    aqua_planet=aqua_planet,       &
                                    brnch_retain_casename=brnch_retain_casename )
  if ( .not. atm_adiabatic )then
     call shr_sys_abort( 'atm_adiabatic was NOT set and should have been' )
  end if
  if ( atm_ideal_phys )then
     call shr_sys_abort( 'atm_ideal_phys was set and should not have been' )
  end if
  if (  aqua_planet )then
     call shr_sys_abort( 'aqua_planet was set and should NOT have been' )
  end if
  if ( .not. shr_inputInfo_initIsModelToRun( CCSMInit, "atm") )then
     call shr_sys_abort( 'In '//mode//' mode atm should run' )
  end if
  if ( shr_inputInfo_initIsModelToRun( CCSMInit, "ice") )then
     call shr_sys_abort( 'In '//mode//' mode ice should NOT run' )
  end if
  if ( shr_inputInfo_initIsModelToRun( CCSMInit, "lnd") )then
     call shr_sys_abort( 'In '//mode//' mode lnd should NOT run' )
  end if
  if ( shr_inputInfo_initIsModelToRun( CCSMInit, "ocn") )then
     call shr_sys_abort( 'In '//mode//' mode ocn should NOT run' )
  end if
  if ( perpetual .or. perpetual_ymd /= 0 )then
     call shr_sys_abort( mode//'mode should NOT be in perpetual mode' )
  end if
  ! atm_ideal_phys
  mode = "atm_ideal_phys"
  call shr_inputInfo_initPutData( CCSMInit, atm_adiabatic=.false.,              &
                             atm_ideal_phys= .true., aqua_planet=.false.,  &
                             brnch_retain_casename=.true. )
  if ( mastertask )call shr_inputInfo_initPrint( CCSMInit )
  call shr_inputInfo_initGetData( CCSMInit, perpetual_run=perpetual, perpetual_ymd=perpetual_ymd )
  call shr_inputInfo_initGetData( CCSMInit, atm_adiabatic=atm_adiabatic,   &
                                    atm_ideal_phys=atm_ideal_phys, &
                                    aqua_planet=aqua_planet,       &
                                    brnch_retain_casename=brnch_retain_casename )
  if ( atm_adiabatic )then
     call shr_sys_abort( 'atm_adiabatic was set and should NOT have been' )
  end if
  if ( .not. atm_ideal_phys )then
     call shr_sys_abort( 'atm_ideal_phys was NOT set and should have been' )
  end if
  if (  aqua_planet )then
     call shr_sys_abort( 'aqua_planet was set and should NOT have been' )
  end if
  if ( .not. shr_inputInfo_initIsModelToRun( CCSMInit, "atm") )then
     call shr_sys_abort( 'In '//mode//' mode atm should run' )
  end if
  if ( shr_inputInfo_initIsModelToRun( CCSMInit, "ice") )then
     call shr_sys_abort( 'In '//mode//' mode ice should NOT run' )
  end if
  if ( shr_inputInfo_initIsModelToRun( CCSMInit, "lnd") )then
     call shr_sys_abort( 'In '//mode//' mode lnd should NOT run' )
  end if
  if ( shr_inputInfo_initIsModelToRun( CCSMInit, "ocn") )then
     call shr_sys_abort( 'In '//mode//' mode ocn should NOT run' )
  end if
  if ( perpetual .or. perpetual_ymd /= 0 )then
     call shr_sys_abort( mode//'mode should NOT be in perpetual mode' )
  end if
  if (  .not. brnch_retain_casename )then
     call shr_sys_abort( 'brnch_retain_casename was NOT set and should have been' )
  end if
  ! default
  mode = "default"
  call shr_inputInfo_initPutData( CCSMInit, atm_adiabatic=.false., &
                              atm_ideal_phys= .false., aqua_planet=.false., &
                              brnch_retain_casename=.false. )
  if ( mastertask )call shr_inputInfo_initPrint( CCSMInit )
  call shr_inputInfo_initGetData( CCSMInit, perpetual_run=perpetual, perpetual_ymd=perpetual_ymd )
  call shr_inputInfo_initGetData( CCSMInit, atm_adiabatic=atm_adiabatic,   &
                                    atm_ideal_phys=atm_ideal_phys, &
                                    aqua_planet=aqua_planet,       &
                                    brnch_retain_casename=brnch_retain_casename )
  if ( atm_adiabatic )then
     call shr_sys_abort( 'atm_adiabatic was set and should NOT have been' )
  end if
  if ( atm_ideal_phys )then
     call shr_sys_abort( 'atm_ideal_phys was set and should NOT have been' )
  end if
  if (  aqua_planet )then
     call shr_sys_abort( 'aqua_planet was set and should NOT have been' )
  end if
  if (  brnch_retain_casename )then
     call shr_sys_abort( 'brnch_retain_casename was set and should NOT have been' )
  end if
  if ( .not. shr_inputInfo_initIsModelToRun( CCSMInit, "atm") )then
     call shr_sys_abort( 'In '//mode//' mode atm should run' )
  end if
  if ( .not. shr_inputInfo_initIsModelToRun( CCSMInit, "ice") )then
     call shr_sys_abort( 'In '//mode//' mode ice should run' )
  end if
  if ( .not. shr_inputInfo_initIsModelToRun( CCSMInit, "lnd") )then
     call shr_sys_abort( 'In '//mode//' mode lnd should run' )
  end if
  if ( .not. shr_inputInfo_initIsModelToRun( CCSMInit, "ocn") )then
     call shr_sys_abort( 'In '//mode//' mode ocn should run' )
  end if
  if ( perpetual .or. perpetual_ymd /= 0 )then
     call shr_sys_abort( mode//'mode should NOT be in perpetual mode' )
  end if

  if ( mastertask )print *, 'Testing passed!!!'
  call shr_mpi_finalize( )
end program test_shr_inputinfo


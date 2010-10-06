!!>
!!c*******************************************************************
!!c
!!c   MODULE: mod_service
!!c   PACKAGE VERSION:  regcm4  coming from DLPROTEIN-2.1  package.. 
!!c   ACTION: basic service module: 
!!c           Timing/Debugging/Error routines 
!!c           defines sp and dp parameters 
!!c*******************************************************************
!! 
!!                                  
!!
!! *How to use this module for TIMING*:
!! 1) First the module need to be initialized with : 
!!        call activate_debug
!!
!! 2) The timer is activated when in each subroutine is inserted 
!!        use mod_service 
!!    and having declared the two following local variables: 
!!        character (len=50) :: sub='name_of_your_subroutine' 
!!        integer :: index=0
!!    and as a first instruction
!!        call time_begin(subroutine_name,index)
!!
!! 3) If the subroutine is not a communication subroutine
!!    as last line insert the following call:
!!        call time_end(subroutine_name,index) 
!! 4) If the subroutine is a communication subroutine:
!!    compute the length of the messages exchanged and as
!!    the last instruction insert the following call:
!!        call time_end(subroutine_name,index,length)
!!
!! 5) To print out times during the run:
!!    - insert "call time_print" to print out at desired location: 
!!      it is possible to pass a string to this call for labelling
!!    - insert "call time_reset" to reset and start a new timing section.
!!
!! *How to use this module for DEBUGGING:*
!!  
!! Debugging facilities are enabled by the call
!!      call start_debug
!! and are disabled by the call
!!      call stop_debug
!!                       
!!
!!
!!* Module main variables:*
!!  - _TYPE timing_info_: data structure to store all informations
!!                        on the calling soubroutine
!!  -_type info_serial_ : informations for serial subroutines
!!  -_type info_comm_   : informations  for communication routines
!!  -_integer number_of_entry_ : counter ( unique for the previous two types)
!!  -_real time_et,time_bt_    : array where to store timing info
!!
!!<
MODULE mod_service

  use mod_dynparam , only : debug_level

!!! definition of single and double precision 

    integer, parameter, public :: sp = kind( 1.0 )
    integer, parameter, public :: dp = kind( 1.0d0 )

!! all constant/numbers used in the code should be defined as
!!  DD.DD_dp for double precision numbers: 
!!  DD.DD_sp for single precision numbers: 
!! 
!! examples: 
!!      3.745566_dp : double precision numbers 
!!      4.2_sp      : single precision numbers  
!! 

  TYPE timing_info
     INTEGER :: n_of_time
     CHARACTER (len=50) :: name_of_section
     REAL (kind=8) :: total_time
     INTEGER :: total_size
  END TYPE timing_info

  TYPE (timing_info) , DIMENSION (100) :: info_serial,info_comm
  REAL (kind=8)  :: time_et(100),time_bt(100) 
  INTEGER :: n_of_ENTRY=0
  INTEGER :: nrite=6

  CHARACTER (len=120) :: errmsg   !! a string where to compose an error message

  !! some global variable for debugging purposes 
  !! set by prepare_debug, start_debug and stop_debug subroutines

  INTEGER :: ndebug=28 !! unit for debugging files..
  INTEGER :: Nlevel=0  !! level of depth in printing calling tree.. 
  LOGICAL :: ldebug=.FALSE. !! if true debug is enabled ( set by activate_debug)
!  INTEGER :: debug_level=0  !! Level of information to be printed out 

!!! 

  INTEGER, PRIVATE :: node=0,mxnode=1
  INTEGER, PRIVATE :: called
  INTEGER, PRIVATE, ALLOCATABLE :: a_tmp(:)
  LOGICAL,PRIVATE :: called_mpi=.TRUE.

  !! interface for write_info subroutine

  INTERFACE WRITE_info
     MODULE PROCEDURE printa_i 
     MODULE PROCEDURE printa_r
  END INTERFACE

  !! interface for check_memory subroutine

  INTERFACE check_memory
     MODULE PROCEDURE check_memory_r 
     MODULE PROCEDURE check_memory_i
  END INTERFACE


CONTAINS

  !!>
  !!   ROUTINE : activate_DEBUG
  !!   ACTION : Open Debugging files and activates debugging
  !!            Set variables needed in the timing routines 
  !!            It is executed only if the code is compiled with DEBUG enabled
  !!     
  !!< 
  SUBROUTINE activate_debug(level)

    IMPLICIT NONE
    INTEGER, optional :: LEVEL
    CHARACTER(len=3) ::  np='   '
    CHARACTER(len=9) ::  string
    CHARACTER (len=100) :: record
    CHARACTER (len=50) :: sub='activate_debug'
    LOGICAL :: safe=.TRUE.
    INTEGER :: ierr1,idum
    INTEGER,EXTERNAL :: intstr
#ifdef MPP1

    INCLUDE 'mpif.h'  

    ! check if MPI is on.
    called_mpi=.FALSE.
    CALL MPI_initialized(called_mpi,ierr1)
    IF (.NOT.called_mpi)  CALL MPI_init(ierr1)
    CALL MPI_initialized(called_mpi,ierr1)

    ! Number of processes 
    CALL MPI_COMM_RANK( MPI_COMM_WORLD,node,ierr1)
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD,mxnode,ierr1)

    !! allocate and initialize this vector needed in timing routines..
    ALLOCATE(a_tmp(0:mxnode-1))
    a_tmp=0 


#endif 

    !! the following procedure accounts up to 999 processors  
    IF (mxnode.LT.10) THEN
       WRITE(np(1:1),'(i1)') node
    ELSEIF(mxnode.LT.100)THEN
       IF (node.LT.10) THEN
          np ='0'
          WRITE(np(2:2),'(i1)') node
       ELSE
          WRITE(np(1:2),'(i2)') node
       ENDIF
    ELSE
       IF (node.LT.10) THEN
          np ='00'
          WRITE(np(3:3),'(i1)') node
       ELSEIF(node.LT.100)THEN
          np ='0'
          WRITE(np(2:3),'(i2)') node
       ELSE
          WRITE(np,'(i3)') node
       ENDIF
    END IF

    WRITE(string,'(A6,A3)') "DEBUG_",np
    OPEN(ndebug+node,FILE=string,status='unknown')
     idum=ndebug+node
     CALL write_info(sub, "DEBUGGING FILE CORRECTLY OPEN:unit is",idum)
     if (present(level)) debug_level=level
     CALL write_info(sub,"debug_level not specified in regcm.in file:default is ",debug_level)
     ldebug=.true.

  END SUBROUTINE activate_debug



  !!>
  !!   ROUTINE : START_DEBUG
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : start debug facilites on specific portion of code
  !!< 
  SUBROUTINE start_debug(level,sub,line)

    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(in) :: LEVEL
    CHARACTER*(*), OPTIONAL, INTENT(in) :: sub
    INTEGER, OPTIONAL, INTENT(in) :: line
    CHARACTER(len=80) :: string='   '
    CHARACTER(len=50) :: substr='not specified'
    CHARACTER(len=8)  :: sline =' no spec'

#ifdef DEBUG
    string=' '
    substr='not specified'
    sline =' no spec'
    IF (.NOT.CALLED_MPI) THEN
       CALL error_prot(sub='start_debug',err_code=3001,  &
            message='error : process unkown to debugging facility')
    END IF
    IF (PRESENT(SUB)) substr=sub
    IF (PRESENT(LINE)) WRITE(SLINE,'(1X,I6)') LINE
    WRITE(string,'(A29,A20,A6,A8)') &
         'Debugging started in routine:',substr(1:20),', line',SLINE  
    WRITE(ndebug+NODE,'(A80)') string
    CALL flusha(ndebug+node)
    IF (LDEBUG) THEN 
       WRITE(ndebug+NODE,*) 'Debug already active, debug level=',debug_level
    ELSE 
       LDEBUG=.TRUE.
       IF (PRESENT(level)) debug_level=level
    END IF

#endif
  END SUBROUTINE START_DEBUG



  !!>
  !!   ROUTINE : STOP_DEBUG
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : stop debug facilities
  !!<   
  SUBROUTINE STOP_debug(level,sub,line)

    IMPLICIT NONE
    INTEGER, OPTIONAL :: LEVEL
    CHARACTER*(*), OPTIONAL, INTENT(in) :: sub
    INTEGER, OPTIONAL, INTENT(in) :: line
    CHARACTER(len=80) :: string=' '
    CHARACTER(len=50) :: substr=' not specified'
    CHARACTER(len=8)  :: sline =' no spec'

#ifdef DEBUG
    string=' '
    substr=' not specified'
    sline =' no spec'
    IF (PRESENT(level)) debug_level=0 
    IF (PRESENT(level)) debug_level=level
    IF (PRESENT(SUB)) substr=sub
    IF (PRESENT(LINE)) WRITE(SLINE,'(1X,I6)') LINE


    WRITE(string,'(A24,A20,A6,A8)') &
         'Ending debug in routine:',substr(1:20),", line",SLINE  
    LDEBUG=.FALSE.
    debug_level=0
    WRITE(ndebug+NODE,'(A80)') string
    CALL flusha(ndebug+node)
#endif
  END SUBROUTINE STOP_DEBUG


  !!>
  !!   ROUTINE : TIME_BEGIN
  !!   ACTION : start recording time for the calling subroutine
  !!<
  SUBROUTINE time_begin(name,index)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER, INTENT(INOUT) :: index
    CHARACTER (len=50) :: name
    CHARACTER (len=50) :: stringa='                                      '

#ifdef DEBUG
    IF (index==0) THEN 
       n_of_ENTRY=n_of_ENTRY+1
       index=n_of_ENTRY
    END IF

    time_bt(index)=timer()

    ! debugging stuff: if debug_level is greater than 3 print out name of the  
    ! routine ( just 4 levels depth..)

    IF (ldebug.AND.(debug_level.GT.3)) THEN
       nlevel=nlevel+1
       IF (nlevel.LE.4) THEN 
          WRITE(ndebug+node,*) stringa(1:nlevel*2),name(1:20),"(in)",nlevel
          CALL flusha(ndebug+node)
       ENDIF
    END IF
#endif 
  END SUBROUTINE time_begin

  !!> 
  !!   ROUTINE : TIME_END
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : end recording time for the calling subroutine
  !!<
  SUBROUTINE time_end(name_of_section,index,isize)
    IMPLICIT NONE
    CHARACTER (len=50) ::  name_of_section
    INTEGER, INTENT(IN) :: index
    INTEGER, OPTIONAL :: isize
    REAL (Kind=8)  :: time_CALL
    CHARACTER (len=50) :: name
    CHARACTER (len=50)  :: stringa=' '

#ifdef DEBUG    

    ! check for index: should not be less than zero
    IF (index<0) THEN 
       CALL error_prot(sub='time_end',err_code=3502,  &
            message='index less then 0: bad thing !! ')
    END IF

    ! debugging stuff 
    IF (ldebug.AND.(debug_level.GT.3)) THEN
       IF (nlevel.LE.4) THEN 
          WRITE(ndebug+node,*) &
               stringa(1:nlevel*2),name_of_section(1:20),"(out)",nlevel
          CALL flusha(ndebug+node)
       END IF
       nlevel=nlevel-1
    END IF

    !
    time_et(index)=timer()
    time_CALL=time_et(index)-time_bt(index)

    IF (PRESENT(isize)) THEN 
       info_comm(index)%total_size=info_comm(index)%total_size+isize
       info_comm(index)%name_of_section=name_of_section
       info_comm(index)%n_of_time=info_comm(index)%n_of_time+1
       info_comm(index)%total_time=info_comm(index)%total_time+time_CALL
    ELSE
       info_serial(index)%name_of_section=name_of_section
       info_serial(index)%n_of_time=info_serial(index)%n_of_time+1
       info_serial(index)%total_time=info_serial(index)%total_time+time_CALL
    END IF
#endif 

  END SUBROUTINE time_end


  !!>
  !!   ROUTINE : TIME_PRINT
  !!   ACTION : print out times collected among processors
  !!            on unit iunit (meant to be the OUTPUT file)
  !!<

  SUBROUTINE time_print(iunit,name_of_section)

    IMPLICIT NONE
    ! arguments:
    CHARACTER (len=*),OPTIONAL :: name_of_section
    INTEGER :: iunit
    ! local variables:  
    INTEGER :: ENTRY
    INTEGER :: imin,imax
    INTEGER :: i,test,len
    REAL (kind=8)  :: avg,xmin,xmax
    REAL (kind=8)  , ALLOCATABLE :: array_tmp(:)
    INTEGER , ALLOCATABLE :: array_entries(:)
    LOGICAL :: called=.TRUE.,FLAG_ENTRY=.FALSE.,L_TIMES_on_PE=.FALSE.
    LOGICAL :: L_ENTRY=.TRUE.
    INTEGER :: nwords,ierr,nword_send,nword_receive,ierr1
    REAL    :: total_comm_time=0.0
    REAL    :: avg_value
    CHARACTER (len=50) :: name
    CHARACTER (len=50) :: sub='time_print'

#ifdef DEBUG

#ifdef MPP1
    INCLUDE 'mpif.h'  
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr1)
#endif  

    L_ENTRY=.TRUE.
    IF (node==0) THEN
       ! WRITE(iunit,"(80('-'))")
       IF(PRESENT(name_of_section))  THEN
          len=len_strim(name_of_section)

          WRITE(iunit,"(/,1x,10('!'),' Specific TIMING for section: ',a30,1x,10('!'))") &    
               name_of_section(1:len)
       ELSE
          WRITE(iunit,"(/,1x,10('!'),' Specific TIMING up to checkpoint ',30('!'))")
       END IF
       WRITE(iunit,"(5('!'),19x,a30,/)") '!Serial subroutines : '
       WRITE(iunit,105) 'section','times','avg-time','max(PE)','min(PE)'
    END IF
105 FORMAT(1x,'!',A7,9x,A10,1x,A9,1x,A15,1x,A10)   

    ALLOCATE(array_tmp(0:mxnode-1))
    ALLOCATE(array_entries(0:mxnode-1))

#ifdef MPP1
    ! check if the calling tree is the same on all nodes 
    CALL gather_i(array_entries,n_of_ENTRY)
    test=array_entries(0)  
    DO i=1,mxnode-1
       IF (array_entries(i)/=test) THEN 
          WRITE(ndebug+node,*) ' Warning:Different trees on different pe:',  &
               n_of_ENTRY
          CALL error_prot(sub='time_end',err_code=-3503,  &
               message='different trees on different pe! ')
          L_ENTRY=.FALSE.
       END IF
    END DO
#endif  


    ! if the calling tree is the same print out gathered data 
    ! on iunit (OUTPUT file)
    IF (L_ENTRY) THEN
       DO ENTRY=1,n_of_ENTRY

#ifdef MPP1 
          L_times_on_pe=.FALSE.
          CALL gather_i(a_tmp,info_serial(ENTRY)%n_of_time)
          ! check the number of time
          test=a_tmp(0)
          avg_value= SUM(a_tmp)/mxnode 
          DO i=0,mxnode-1
             IF (a_tmp(i)/=test) THEN 
                L_times_on_pe=.TRUE.
             END IF
             CALL MPI_barrier(MPI_COMM_WORLD,ierr1)
          END DO
#endif
          ! set to zero times less then 0.1 microseconds
          IF (info_serial(ENTRY)%total_time<=0.0000001) & 
               info_serial(ENTRY)%total_time=0.0000001

          IF (info_serial(ENTRY)%n_of_time>=1) THEN
             CALL gather(array_tmp,info_serial(ENTRY)%total_time)
             ! compute avg, max and min 
             CALL av_max_MIN(array_tmp,avg,xmax,imax,xmin,imin)
             IF (node==0) THEN
                NAME=info_serial(ENTRY)%name_of_section
                IF (L_times_on_pe) &
                     NAME = "*" // info_serial(ENTRY)%name_of_section 
                WRITE(iunit,100) NAME(1:15), & 
                     info_serial(ENTRY)%n_of_time, &
                     avg,xmax,imax-1,xmin,imin-1
             END IF

          END IF
       END DO
    ELSE 
       IF (node==0) THEN
          WRITE (iunit,*) 'No global times for serial routines: different calling tree'
       ENDIF
    ENDIF

    !save anyway the data for each pe on debug unit:
    if (ldebug) then  
    CALL write_info(sub,'Specific local time up to checkpoint for node',node) 
    DO ENTRY=1,n_of_ENTRY
       NAME=info_serial(ENTRY)%name_of_section
       IF (L_times_on_pe) NAME = "*" // info_serial(ENTRY)%name_of_section
       IF (info_serial(ENTRY)%n_of_time>=1) THEN 
          WRITE(ndebug+NODE,102) NAME, & 
               info_serial(ENTRY)%n_of_time, &
               info_serial(ENTRY)%total_time
          call flusha(ndebug+NODE)
       ENDIF
    END DO
    endif 

#ifdef MPP1

    IF (node==0) THEN
       WRITE(iunit,"(1x,'!',19x,a30,/)") 'Communication subroutines :   '
       WRITE(iunit,*) &
            '!! section       times avg-time  max(PE)  min(PE) data  ratio(Mb/sec)'
    END IF
    !!   
    DO ENTRY=1,n_of_ENTRY
       ! set to zero times less then 0.1 microseconds
       IF (info_comm(ENTRY)%total_time<=0.0000001) & 
            info_comm(ENTRY)%total_time=0.0000001

       IF (info_comm(ENTRY)%n_of_time>=1) THEN
          CALL gather(array_tmp,info_comm(ENTRY)%total_time)
          CALL av_max_MIN(array_tmp,avg,xmax,imax,xmin,imin)
          IF (node==0) THEN
             WRITE(iunit,101) info_comm(ENTRY)%name_of_section, & 
                  info_comm(ENTRY)%n_of_time, &
                  avg,xmax,imax-1,xmin,imin-1, & 
                  info_comm(ENTRY)%total_size, &
                  info_comm(ENTRY)%total_size/(info_comm(ENTRY)%total_time*1e6)
             total_comm_time=total_comm_time+avg
          END IF
       END IF
    END DO
    IF(node==0) THEN  
       WRITE(iunit,"(1x,'!',a40,3x,f9.4)") &
            'total communication time:',total_comm_time
       total_comm_time=0.0

    END IF
#endif
    IF(node==0) THEN  
       WRITE(iunit,"(/,1x,10('!'),' End of Specific TIMING ', 47('!'),/)")
       call flusha(iunit)
    ENDIF
100 FORMAT(1x,'!',A15,1x,I10,1x,F9.4,1x,F9.4,'(',i3,')',1x,f9.4,'(',i3,')')   

101 FORMAT(1x,'!',A15,1x,I10,1x,F9.4,1x,F9.4,'(',i3,')',1x,f9.4,'(',i3,')',I12,1X,F13.3)   
102 FORMAT(1x,'!',A15,1x,I10,1x,F9.4,1x,F9.4)   

    DEALLOCATE(array_tmp) 
    DEALLOCATE(array_entries) 

#endif
  END SUBROUTINE time_print


  !!>
  !!   ROUTINE : TIME_RESET
  !!   ACTION :reset all the data structures 
  !!<
  SUBROUTINE time_reset
    IMPLICIT NONE 
    INTEGER :: ENTRY

#ifdef DEBUG

    DO ENTRY=1,100
       info_serial(ENTRY)%n_of_time=0
       info_serial(ENTRY)%total_time=0.0
       info_serial(ENTRY)%total_size=0
       info_comm(ENTRY)%n_of_time=0
       info_comm(ENTRY)%total_time=0.0
       info_comm(ENTRY)%total_size=0

    END DO
#endif
  END SUBROUTINE time_reset


  !!>
  !!   ROUTINE : GATHER
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : gather and fill an array of data from different processors 
  !!<

  SUBROUTINE gather(f_collect,f_sub)
    IMPLICIT NONE 
    REAL (kind=8), DIMENSION(:)  :: f_collect 
    REAL (kind=8) :: f_sub

#ifdef DEBUG
#ifndef MPP1 
    f_collect(1)=f_sub
#else
    INCLUDE 'mpif.h'  
    LOGICAL :: called
    INTEGER :: nwords,ierr,nword_send,nword_receive
    Nword_send=1
    Nword_receive=Nword_send
    CALL MPI_Allgather(f_sub,Nword_send,MPI_DOUBLE_PRECISION,f_collect, &
         Nword_receive,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    IF (ierr /=0) THEN
       CALL error_prot(sub='gather_in_timing_mod',err_code=ierr,  &
            message=' error in MPI_allgather!! ')

    END IF
#endif 
#endif 
  END SUBROUTINE gather



#ifdef MPP1 
  !!>
  !!   ROUTINE : GATHER_I
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : another gathering routine
  !!<
  SUBROUTINE gather_i(f_collect,f_sub)
    IMPLICIT NONE 
    INCLUDE 'mpif.h'  
    ! assumed shaped array... 
    INTEGER , DIMENSION(:)  :: f_collect 
    INTEGER  :: f_sub
    LOGICAL :: called
    INTEGER :: nwords,ierr,nword_send,nword_receive

#ifdef DEBUG
    Nword_send=1
    Nword_receive=Nword_send
    CALL MPI_Allgather(f_sub,Nword_send,MPI_INTEGER,f_collect, &
         Nword_receive,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    IF (ierr /=0) THEN
       CALL error_prot(sub='gather_in_timing_mod',err_code=ierr,  &
            message=' error in MPI_allgather!! ')
    END IF
#endif
  END SUBROUTINE gather_i
#endif

  !!>
  !!
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : compute average, maximum and minimum of array and indices
  !!<
  SUBROUTINE av_max_MIN(array,avg,xmax,index_max,xmin,index_min) 
    INTEGER :: index_min,index_max,n_elements
    REAL (kind=8)    :: xmax,xmin,avg
    REAL (kind=8), DIMENSION(:) :: array
    INTEGER :: i

#ifdef DEBUG
    xmax=0.0
    xmin=1.e6
    avg=0.0
    index_max=0
    index_min=0
    n_elements=SIZE(array)
    DO i=1,n_elements 
       avg=avg+array(i)
       IF (array(i).GT.xmax) THEN 
          xmax=array(i)
          index_max=i
       END IF
       IF (array(i).LT.xmin) THEN 
          xmin=array(i)
          index_min=i
       END IF
    END DO

    avg=avg/float(n_elements)

#endif

  END SUBROUTINE av_max_min


  !!>
  !!   ROUTINE : ERROR_PROT
  !!   PACKAGE VERSION : regcm4 (inerithed by dlprotein2.1) 
  !!   ACTION : general purpose error/debug routine. Arguments are
  !!            sub: name of calling subroutine
  !!            ierr > 0: warning. Returns control to calling subroutine.
  !!            ierr < 0: fatal error. Stops.
  !!            message (optional): error message, as explicit as possible
  !!            line (optional): source line (cpp '__LINE__' macro)
  !!<
  SUBROUTINE error_prot(sub,err_code,message,line) 


    IMPLICIT NONE
#ifdef MPP1 
    INCLUDE 'mpif.h'    
#endif
    CHARACTER*(*), INTENT(in) :: sub
    INTEGER, INTENT(in) :: err_code
    CHARACTER*(*), OPTIONAL, INTENT(in) :: message
    INTEGER, OPTIONAL, INTENT(in) :: line
    ! local
    CHARACTER(len=11) :: error_TYPE
    INTEGER :: ierr 

    error_TYPE = 'Error'


    IF (err_code < 0) error_TYPE = 'Warning'

    IF (node==0) THEN

       IF (err_code > 0) THEN

          WRITE (nrite,'(/1x,78(''*''))')
          WRITE (nrite,'(3x,a8,'' no.'',i5,'' from '',a)') &
               error_TYPE, err_code, sub

          ! optional strings
          IF (PRESENT(line)) WRITE(nrite,'(5x,''line:    '',i6)') line
          IF (PRESENT(message)) WRITE(nrite,'(/1x,''>>> '',a)') message

          WRITE (nrite,'(/1x,78(''*''))')

          WRITE (nrite,'(1x,a)') 'Stopping the code...'
       ELSE

          WRITE (nrite,'(3x,a8,'' no.'',i5,'' from '',a,'' : '',a)') &
               error_TYPE, err_code, sub, message

       ENDIF
       CALL flusha(nrite)

    ENDIF

    !! close all file and shutdown 
    IF (err_code > 0) THEN

       CLOSE (nrite)
#ifndef  MPP1 
       STOP 
#else
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       CALL MPI_finalize(ierr) 
#endif 
    ENDIF


  END SUBROUTINE error_prot



  !!>
  !!   ROUTINE : PRINTA_I
  !!   ACTION : general purpose error printing routine for integers
  !!<
  SUBROUTINE printa_i(sub,variable,value,line)

    IMPLICIT NONE
    CHARACTER*(*), INTENT(in) :: sub
    CHARACTER*(*), INTENT(in) :: variable
    INTEGER, INTENT (in) :: value
    INTEGER, OPTIONAL, INTENT(in) :: line

#ifdef DEBUG
    IF (PRESENT(line)) THEN 
       WRITE(ndebug+node,'(A20,'':at line'',I6,A,i10)') &
            sub(1:20),line,variable,value
    ELSE
       WRITE(ndebug+node,'(A20,'':'',A,i10)') sub(1:20),variable,value
    END IF
    CALL flusha(ndebug+node)

#endif
  END SUBROUTINE printa_i

  !!>
  !!   ROUTINE : PRINTA_R
  !!   ACTION : general purpose error printing routine for reals
  !!<

  SUBROUTINE printa_r(sub,variable,value,line)


    IMPLICIT NONE
    CHARACTER*(*), INTENT(in) :: sub
    CHARACTER*(*), INTENT(in) :: variable
    REAL(kind=8), INTENT (in) :: value
    INTEGER, OPTIONAL, INTENT(in) :: line

#ifdef DEBUG
    IF (PRESENT(line)) THEN 
       WRITE(ndebug+node,'(A20,'':at line'',I6,A,F20.10)') &
            sub(1:20),line,variable,value
    ELSE
       WRITE(ndebug+node,'(A20,'':'',A,f20.10)') sub(1:20),variable,value
    END IF
    CALL flusha(ndebug+node)

#endif
  END SUBROUTINE printa_r

  !!>
  !!   ROUTINE : E_ALLOCA
  !!   ACTION :signal error on allocation
  !!<                                                                     
  SUBROUTINE e_alloca(sub,variable,line)

    IMPLICIT NONE
    CHARACTER*(*), INTENT(in) :: sub
    CHARACTER*(*), INTENT(in) :: variable
    CHARACTER*10  :: sub_e='e_alloca'
    CHARACTER*80  :: string=' ' 
    INTEGER, OPTIONAL, INTENT(in) :: line

    IF (PRESENT(line)) THEN 
       WRITE(ndebug+node,'(A20,'':at line'',I6,A24,A)') & 
            sub(1:20),line,'error allocating ', variable
       WRITE(string,'(A10,'':at line'',I6,A24,A20)') & 
            sub(1:20),line,'error allocating ', variable(1:20)
    ELSE
       WRITE(ndebug+node,'(A20,'':'',A24,A20)') &
            sub(1:20),'error allocating: ',variable(1:20)
       WRITE(string,'(A20,'':'',A24,A20)') &
            sub(1:20),'error allocating: ',variable(1:20)
    END IF
    CALL flusha(ndebug+node)
    !! if an error occours stops the program... 
    CALL error_prot(sub_e,3000,string)
  END SUBROUTINE e_alloca


  !!>
  !! ROUTINE : CHECK_MEMORY_R
  !! ACTION : check the address of real variable
  !!<
  SUBROUTINE check_memory_r(variable,name_of_variable,sub,line)
    IMPLICIT NONE
    CHARACTER*(*), OPTIONAL, INTENT(in) :: name_of_variable
    CHARACTER*(*), OPTIONAL, INTENT(in) :: sub
    REAL(kind=8),            INTENT(in) :: variable
    INTEGER,       OPTIONAL, INTENT(in) :: line
    CHARACTER(len=10) :: varia
    CHARACTER(len=20) :: subro
    CHARACTER(len=7) :: sline
    CHARACTER(len=36) :: string

    INTEGER i_addr
    string='    '
    subro= '    '
    varia= '    '
    string='    '
    sline ='    '

    !!i_addr=loc(variable)

    WRITE(string,"(A20,I15)") 'CM_mes>address is =' ,i_addr
    IF (PRESENT(name_of_variable)) varia=name_of_variable 
    IF (PRESENT(sub)) subro=sub
    IF (PRESENT(line))  WRITE(SLINE,'(1X,I6)') LINE
    WRITE(ndebug+node,"(A36,'for ',A10, ' in sub ',A20,'at line=',A7)") &
         string,varia,subro,sline

    CALL flusha(ndebug+node) 

  END SUBROUTINE check_memory_r

  !!>
  !!   ROUTINE : CHECK_MEMORY_I
  !!   ACTION :check the address of real variable
  !!           It may use the loc function (not architecture-universal)
  !!<
  SUBROUTINE check_memory_i(variable,name_of_variable,sub,line)
    IMPLICIT NONE
    CHARACTER*(*), OPTIONAL, INTENT(in) :: name_of_variable
    CHARACTER*(*), OPTIONAL, INTENT(in) :: sub
    INTEGER,            INTENT(in) :: variable
    INTEGER,       OPTIONAL, INTENT(in) :: line
    CHARACTER(len=10) :: varia
    CHARACTER(len=20) :: subro
    CHARACTER(len=7) :: sline
    CHARACTER(len=36) :: string

    INTEGER i_addr
    string='    '
    subro= '    '
    varia= '    '
    string='    '
    sline ='    '

    !! i_addr=loc(variable)

    WRITE(string,"(A20,I15)") 'CM_mes>address is =' ,i_addr
    IF (PRESENT(name_of_variable)) varia=name_of_variable 
    IF (PRESENT(sub)) subro=sub
    IF (PRESENT(line))  WRITE(SLINE,'(1X,I6)') LINE
    WRITE(ndebug+node,"(A36,'for ',A10, ' in sub ',A20,'at line=',A7)") &
         string,varia,subro,sline

    CALL flusha(ndebug+node) 

  END SUBROUTINE check_memory_i


  !!>
  !!   ROUTINE : LEN_STRIM
  !!   ACTION : trimming of string
  !!< 
  FUNCTION len_strim (string) RESULT (len_trim_RESULT)
    IMPLICIT NONE 
    CHARACTER  (len=*), INTENT(IN) :: string
    INTEGER :: len_trim_RESULT,k

    len_trim_RESULT =0 
    DO k= LEN(string),1,-1
       IF (string(k:k) /=' ') THEN
          len_trim_RESULT =k
          EXIT
       END IF
    END DO

  END FUNCTION len_strim
 
  double precision function timer()
  integer(8) C,R,M
    CALL SYSTEM_CLOCK (COUNT=C, COUNT_RATE=R, COUNT_MAX=M)
    timer= dble(C)/dble(R)
    return
  end function

  subroutine flusha( lunit, errno )
  integer, intent(in)  :: lunit
  integer, intent(out),optional :: ERRNO
  if ( PRESENT(ERRNO) ) errno=0

! If whe have a FLUSH, use it 
!  On IBM, flush is flush_ 
#ifdef IBM 
	call FLUSH_(lunit)  
#else
	call flush(lunit)
#endif        
 end subroutine flusha


END MODULE  mod_service

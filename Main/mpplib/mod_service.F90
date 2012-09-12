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

module mod_service

#ifdef DEBUG
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam , only : mycomm , myid , nproc , debug_level

  private

  integer(ik4) , parameter , public :: dbgslen = 64
  integer(ik4) , parameter , public :: dbglinelen = 80
  integer(ik4) , parameter , public :: ndebug = 28 !! unit for debugging files

  type timing_info
    character(len=dbgslen) :: name_of_section
    real(rk8) :: total_time
    integer(ik4) :: n_of_time
    integer(ik4) :: total_size
  end type timing_info

  integer(ik4) , parameter :: maxcall = 100
  type(timing_info) , dimension(maxcall) :: info_serial , info_comm
  real(rk8) :: time_et(maxcall) , time_bt(maxcall) 
  integer(ik4) :: n_of_entry = 0

  character(len=120) :: errmsg   !! a string where to compose an error message

  !! some global variable for debugging purposes 
  !! set by prepare_debug, start_debug and stop_debug subroutines

  integer(ik4) :: nlevel = 0  !! level of depth in printing calling tree.. 
  logical :: ldebug = .false. !! if true debug is enabled

  integer(ik4) :: node = 0
  integer(ik4) :: mxnode = 1
  integer(ik4) , allocatable , dimension(:) :: a_tmp

  !! interface for write_info subroutine

  interface write_info
    module procedure printa_i 
    module procedure printa_r
  end interface

  !! interface for check_memory subroutine

  interface check_memory
    module procedure check_memory_r 
    module procedure check_memory_i
  end interface

  public :: activate_debug , start_debug , stop_debug
  public :: time_begin , time_end , time_reset , time_print

  contains

    subroutine activate_debug(level)
      implicit none
      include 'mpif.h'  
      integer(ik4) , optional :: level
      character(len=3) :: np = '   '
      character(len=9) :: string
      character(len=dbgslen) :: sub = 'activate_debug'
      integer(ik4) :: ierr1 , idum

      ! Number of processes
      node = myid
      mxnode = nproc

      ! allocate and initialize this vector needed in timing routines..
      allocate(a_tmp(0:mxnode-1))
      a_tmp = 0 

      ! the following procedure accounts up to 999 processors  
      if ( mxnode < 10 ) then
        write(np(1:1),'(i1)') node
      else if ( mxnode < 100 ) then
        if (node.LT.10) THEN
          np ='0'
          write(np(2:2),'(i1)') node
        else
          write(np(1:2),'(i2)') node
        end if
      else
        if ( node < 10 ) then
          np ='00'
          write(np(3:3),'(i1)') node
        else if ( node < 100 ) THEN
          np ='0'
          write(np(2:3),'(i2)') node
        else
          write(np,'(i3)') node
        end if
      end if
      write(string,'(a6,a3)') "debug_",np
      open(ndebug+node,file=string,status='unknown')
      idum = ndebug+node
      call write_info(sub,'DEBUGGING FILE CORRECTLY OPEN: unit is ', idum)
      if ( present(level) ) debug_level = level
      call write_info(sub, &
         'debug_level not specified in regcm.in file : default is ',debug_level)
      ldebug = .true.
    end subroutine activate_debug

    subroutine start_debug(level,sub,line)
      implicit none
      integer(ik4) , optional , intent(in) :: level
      character(len=*) , optional , intent(in) :: sub
      integer(ik4) , optional , intent(in) :: line
      character(len=dbglinelen) :: string = '   '
      character(len=dbgslen) :: substr = 'not specified'
      character(len=8) :: sline = ' no spec'
      string = ' '
      substr = 'not specified'
      sline = ' no spec'
      if ( present(sub) ) substr = sub
      if ( present(line) ) write(sline,'(1x,i6)') line
      write(string,'(A29,A20,A6,A8)') &
           'Debugging started in routine:',substr(1:20),', line',sline  
      write(ndebug+node,'(a80)') string
      call flusha(ndebug+node)
      if ( ldebug ) then 
        write(ndebug+node,*) 'Debug already active, debug level = ', debug_level
      else 
        ldebug = .true.
        if ( present(level) ) debug_level = level
      end if
    end subroutine start_debug

    subroutine stop_debug(level,sub,line)
      implicit none
      integer(ik4) , optional :: level
      character(len=*) , optional , intent(in) :: sub
      integer(ik4) , optional , intent(in) :: line
      character(len=dbglinelen) :: string = ' '
      character(len=dbgslen) :: substr = ' not specified'
      character(len=8) :: sline = ' no spec'
      string = ' '
      substr = ' not specified'
      sline = ' no spec'
      if ( present(level) ) debug_level = level
      if ( present(sub) ) substr = sub
      if ( present(line) ) write(SLINE,'(1X,I6)') LINE
      write(string,'(a24,a20,a6,a8)') &
           'Ending debug in routine:', substr(1:20), ", line", sline  
      ldebug = .false.
      debug_level = 0
      write(ndebug+node,'(a80)') string
      call flusha(ndebug+node)
    end subroutine stop_debug

    subroutine time_begin(name,indx)
      implicit none
      integer(ik4) , intent(inout) :: indx
      character(len=dbgslen) :: name
      character(len=dbgslen) :: stringa
      if ( indx == 0 ) then 
        n_of_entry = n_of_entry + 1
        indx = n_of_entry
      end if
      time_bt(indx) = timer()
      ! debugging stuff: if debug_level is greater than 3 print out name of the
      ! routine ( just 4 levels depth..)
      if ( ldebug .and. (debug_level > 3) ) then
        nlevel = nlevel+1
        if ( nlevel <= 4 ) then 
          write(ndebug+node,*) stringa(1:nlevel*2),name(1:20),"(in)",nlevel
          call flusha(ndebug+node)
        end if
      end if
    end subroutine time_begin

  !!> 
  !!   ROUTINE : TIME_END
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : end recording time for the calling subroutine
  !!<
  subroutine time_end(name_of_section,indx,isize)
    implicit none
    character (len=dbgslen) ::  name_of_section
    integer(ik4), intent(IN) :: indx
    integer(ik4), optional :: isize
    real (Kind=8)  :: time_call
    character (len=dbgslen)  :: stringa=' '

    ! check for indx: should not be less than zero
    IF (indx<0) THEN 
       call error_prot(sub='time_end',err_code=3502,  &
            message='indx less then 0: bad thing !! ')
    END IF

    ! debugging stuff 
    IF (ldebug.AND.(debug_level.GT.3)) THEN
       IF (nlevel.LE.4) THEN 
          write(ndebug+node,*) &
               stringa(1:nlevel*2),name_of_section(1:20),"(out)",nlevel
          call flusha(ndebug+node)
       END IF
       nlevel=nlevel-1
    END IF

    !
    time_et(indx)=timer()
    time_call=time_et(indx)-time_bt(indx)

    IF (present(isize)) THEN 
       info_comm(indx)%total_size=info_comm(indx)%total_size+isize
       info_comm(indx)%name_of_section=name_of_section
       info_comm(indx)%n_of_time=info_comm(indx)%n_of_time+1
       info_comm(indx)%total_time=info_comm(indx)%total_time+time_call
    ELSE
       info_serial(indx)%name_of_section=name_of_section
       info_serial(indx)%n_of_time=info_serial(indx)%n_of_time+1
       info_serial(indx)%total_time=info_serial(indx)%total_time+time_call
    END IF

  end subroutine time_end


  !!>
  !!   ROUTINE : TIME_PRINT
  !!   ACTION : print out times collected among processors
  !!            on unit iunit (meant to be the OUTPUT file)
  !!<

  subroutine time_print(iunit,name_of_section)

    implicit none
    INCLUDE 'mpif.h'  
    ! arguments:
    character (len=*),optional :: name_of_section
    integer(ik4) :: iunit
    ! local variables:  
    integer(ik4) :: ENTRY
    integer(ik4) :: imin,imax
    integer(ik4) :: i,test,len
    real (rk8)  :: avg,xmin,xmax
    real (rk8)  , ALLOCATABLE :: array_tmp(:)
    integer(ik4) , ALLOCATABLE :: array_entries(:)
    LOGICAL :: L_TIMES_on_PE=.FALSE.
    LOGICAL :: L_ENTRY=.TRUE.
    integer(ik4) :: ierr1
    real(dp):: total_comm_time=0.0D0
    real    :: avg_value
    character (len=128) :: name
    character (len=dbgslen) :: sub='time_print'

    call MPI_BARRIER(mycomm,ierr1)

    L_ENTRY=.TRUE.
    IF (node==0) THEN
       ! write(iunit,"(80('-'))")
       IF(present(name_of_section))  THEN
          len=len_strim(name_of_section)

          write(iunit,"(/,1x,10('!'),' Specific TIMING for section: ',a30,1x,10('!'))") &    
               name_of_section(1:len)
       ELSE
          write(iunit,"(/,1x,10('!'),' Specific TIMING up to checkpoint ',30('!'))")
       END IF
       write(iunit,"(5('!'),19x,a30,/)") '!Serial subroutines : '
       write(iunit,105) 'section','times','avg-time','max(PE)','min(PE)'
    END IF
105 FORMAT(1x,'!',A7,9x,A10,1x,A9,1x,A15,1x,A10)   

    ALLOCATE(array_tmp(0:mxnode-1))
    ALLOCATE(array_entries(0:mxnode-1))

    ! check if the calling tree is the same on all nodes 
    call gather_i(array_entries,n_of_ENTRY)
    test=array_entries(0)  
    DO i=1,mxnode-1
       IF (array_entries(i)/=test) THEN 
          write(ndebug+node,*) ' Warning:Different trees on different pe:',  &
               n_of_ENTRY
          call error_prot(sub='time_end',err_code=-3503,  &
               message='different trees on different pe! ')
          L_ENTRY=.FALSE.
       END IF
    END DO

    ! if the calling tree is the same print out gathered data 
    ! on iunit (OUTPUT file)
    IF (L_ENTRY) THEN
       DO ENTRY=1,n_of_ENTRY

          L_times_on_pe=.FALSE.
          call gather_i(a_tmp,info_serial(ENTRY)%n_of_time)
          ! check the number of time
          test=a_tmp(0)
          avg_value= SUM(a_tmp)/mxnode 
          DO i=0,mxnode-1
             IF (a_tmp(i)/=test) THEN 
                L_times_on_pe=.TRUE.
             END IF
             call MPI_barrier(mycomm,ierr1)
          END DO
          ! set to zero times less then 0.1 microseconds
          IF (info_serial(ENTRY)%total_time<=0.0000001) & 
               info_serial(ENTRY)%total_time=0.0000001

          IF (info_serial(ENTRY)%n_of_time>=1) THEN
             call gather(array_tmp,info_serial(ENTRY)%total_time)
             ! compute avg, max and min 
             call av_max_MIN(array_tmp,avg,xmax,imax,xmin,imin)
             IF (node==0) THEN
                NAME=info_serial(ENTRY)%name_of_section
                IF (L_times_on_pe) &
                     NAME = "*" // info_serial(ENTRY)%name_of_section 
                write(iunit,100) NAME(1:15), & 
                     info_serial(ENTRY)%n_of_time, &
                     avg,xmax,imax-1,xmin,imin-1
             END IF

          END IF
       END DO
    ELSE 
       IF (node==0) THEN
          write (iunit,*) 'No global times for serial routines: different calling tree'
       ENDIF
    ENDIF

    !save anyway the data for each pe on debug unit:
    if (ldebug) then  
    call write_info(sub,'Specific local time up to checkpoint for node',node) 
    DO ENTRY=1,n_of_ENTRY
       NAME=info_serial(ENTRY)%name_of_section
       IF (L_times_on_pe) NAME = "*" // info_serial(ENTRY)%name_of_section
       IF (info_serial(ENTRY)%n_of_time>=1) THEN 
          write(ndebug+NODE,102) NAME, & 
               info_serial(ENTRY)%n_of_time, &
               info_serial(ENTRY)%total_time
          call flusha(ndebug+NODE)
       ENDIF
    END DO
    endif 

    IF (node==0) THEN
       write(iunit,"(1x,'!',19x,a30,/)") 'Communication subroutines :   '
       write(iunit,*) &
            '!! section       times avg-time  max(PE)  min(PE) data  ratio(Mb/sec)'
    END IF
    !!   
    DO ENTRY=1,n_of_ENTRY
       ! set to zero times less then 0.1 microseconds
       IF (info_comm(ENTRY)%total_time<=0.0000001) & 
            info_comm(ENTRY)%total_time=0.0000001

       IF (info_comm(ENTRY)%n_of_time>=1) THEN
          call gather(array_tmp,info_comm(ENTRY)%total_time)
          call av_max_MIN(array_tmp,avg,xmax,imax,xmin,imin)
          IF (node==0) THEN
             write(iunit,101) info_comm(ENTRY)%name_of_section, & 
                  info_comm(ENTRY)%n_of_time, &
                  avg,xmax,imax-1,xmin,imin-1, & 
                  info_comm(ENTRY)%total_size, &
                  info_comm(ENTRY)%total_size/(info_comm(ENTRY)%total_time*1e6)
             total_comm_time=total_comm_time+avg
          END IF
       END IF
    END DO
    IF(node==0) THEN  
       write(iunit,"(1x,'!',a40,3x,f9.4)") &
            'total communication time:',total_comm_time
       total_comm_time=0.0

    END IF
    IF(node==0) THEN  
       write(iunit,"(/,1x,10('!'),' End of Specific TIMING ', 47('!'),/)")
       call flusha(iunit)
    ENDIF
100 FORMAT(1x,'!',A15,1x,I10,1x,F9.4,1x,F9.4,'(',i3,')',1x,f9.4,'(',i3,')')   

101 FORMAT(1x,'!',A15,1x,I10,1x,F9.4,1x,F9.4,'(',i3,')',1x,f9.4,'(',i3,')',I12,1X,F13.3)   
102 FORMAT(1x,'!',A15,1x,I10,1x,F9.4,1x,F9.4)   

    DEALLOCATE(array_tmp) 
    DEALLOCATE(array_entries) 

  end subroutine time_print


  !!>
  !!   ROUTINE : TIME_RESET
  !!   ACTION :reset all the data structures 
  !!<
  subroutine time_reset
    implicit none 
    integer(ik4) :: ENTRY

    DO ENTRY=1,100
       info_serial(ENTRY)%n_of_time=0
       info_serial(ENTRY)%total_time=0.0
       info_serial(ENTRY)%total_size=0
       info_comm(ENTRY)%n_of_time=0
       info_comm(ENTRY)%total_time=0.0
       info_comm(ENTRY)%total_size=0

    END DO
  end subroutine time_reset


  !!>
  !!   ROUTINE : GATHER
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : gather and fill an array of data from different processors 
  !!<

  subroutine gather(f_collect,f_sub)
    implicit none 
    INCLUDE 'mpif.h'  
    real (rk8), DIMENSION(:)  :: f_collect 
    real (rk8) :: f_sub

    integer(ik4) :: ierr,nword_send,nword_receive
    Nword_send=1
    Nword_receive=Nword_send
    call MPI_Allgather(f_sub,Nword_send,MPI_DOUBLE_PRECISION,f_collect, &
         Nword_receive,MPI_DOUBLE_PRECISION,mycomm,ierr)
    IF (ierr /=0) THEN
       call error_prot(sub='gather_in_timing_mod',err_code=ierr,  &
            message=' error in MPI_allgather!! ')

    END IF
  end subroutine gather

  !!>
  !!   ROUTINE : GATHER_I
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : another gathering routine
  !!<
  subroutine gather_i(f_collect,f_sub)
    implicit none 
    INCLUDE 'mpif.h'  
    ! assumed shaped array... 
    integer(ik4) , DIMENSION(:)  :: f_collect 
    integer(ik4)  :: f_sub
    integer(ik4) :: ierr,nword_send,nword_receive

    Nword_send=1
    Nword_receive=Nword_send
    call MPI_Allgather(f_sub,Nword_send,MPI_integer,f_collect, &
         Nword_receive,MPI_integer,mycomm,ierr)
    IF (ierr /=0) THEN
       call error_prot(sub='gather_in_timing_mod',err_code=ierr,  &
            message=' error in MPI_allgather!! ')
    END IF
  end subroutine gather_i

  !!>
  !!
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : compute average, maximum and minimum of array and indices
  !!<
  subroutine av_max_MIN(array,avg,xmax,indx_max,xmin,indx_min) 
    integer(ik4) :: indx_min,indx_max,n_elements
    real (rk8)    :: xmax,xmin,avg
    real (rk8), DIMENSION(:) :: array
    integer(ik4) :: i

    xmax=0.0
    xmin=1.e6
    avg=0.0
    indx_max=0
    indx_min=0
    n_elements=SIZE(array)
    DO i=1,n_elements 
       avg=avg+array(i)
       IF (array(i).GT.xmax) THEN 
          xmax=array(i)
          indx_max=i
       END IF
       IF (array(i).LT.xmin) THEN 
          xmin=array(i)
          indx_min=i
       END IF
    END DO

    avg=avg/dble(n_elements)

  end subroutine av_max_min


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
  subroutine error_prot(sub,err_code,message,line) 

    implicit none
    INCLUDE 'mpif.h'  
    character*(*), intent(in) :: sub
    integer(ik4), intent(in) :: err_code
    character*(*), optional, intent(in) :: message
    integer(ik4), optional, intent(in) :: line
    ! local
    character(len=11) :: error_TYPE
    integer(ik4) :: ierr 

    error_TYPE = 'Error'


    IF (err_code < 0) error_TYPE = 'Warning'

    IF (node==0) THEN

       IF (err_code > 0) THEN

          write (stderr,'(/1x,78(''*''))')
          write (stderr,'(3x,a8,'' no.'',i5,'' from '',a)') &
               error_TYPE, err_code, sub

          ! optional strings
          IF (present(line)) write(stderr,'(5x,''line:    '',i6)') line
          IF (present(message)) write(stderr,'(/1x,''>>> '',a)') message

          write (stderr,'(/1x,78(''*''))')

          write (stderr,'(1x,a)') 'Stopping the code...'
       ELSE

          write (stderr,'(3x,a8,'' no.'',i5,'' from '',a,'' : '',a)') &
               error_TYPE, err_code, sub, message

       ENDIF
       call flusha(stderr)

    ENDIF

    !! close all file and shutdown 
    IF (err_code > 0) THEN

       CLOSE (stderr)
       call MPI_BARRIER(mycomm,ierr) 
       call MPI_finalize(ierr) 
    ENDIF

  end subroutine error_prot



  !!>
  !!   ROUTINE : PRINTA_I
  !!   ACTION : general purpose error printing routine for integers
  !!<
  subroutine printa_i(sub,variable,value,line)

    implicit none
    character*(*), intent(in) :: sub
    character*(*), intent(in) :: variable
    integer(ik4), intent (in) :: value
    integer(ik4), optional, intent(in) :: line

    IF (present(line)) THEN 
       write(ndebug+node,'(A20,'':at line'',I6,A,i10)') &
            sub(1:20),line,variable,value
    ELSE
       write(ndebug+node,'(A20,'':'',A,i10)') sub(1:20),variable,value
    END IF
    call flusha(ndebug+node)

  end subroutine printa_i

  !!>
  !!   ROUTINE : PRINTA_R
  !!   ACTION : general purpose error printing routine for reals
  !!<

  subroutine printa_r(sub,variable,value,line)

    implicit none
    character*(*), intent(in) :: sub
    character*(*), intent(in) :: variable
    real(rk8), intent (in) :: value
    integer(ik4), optional, intent(in) :: line

    IF (present(line)) THEN 
       write(ndebug+node,'(A20,'':at line'',I6,A,F20.10)') &
            sub(1:20),line,variable,value
    ELSE
       write(ndebug+node,'(A20,'':'',A,f20.10)') sub(1:20),variable,value
    END IF
    call flusha(ndebug+node)

  end subroutine printa_r

  !!>
  !!   ROUTINE : E_ALLOCA
  !!   ACTION :signal error on allocation
  !!<                                                                     
  subroutine e_alloca(sub,variable,line)

    implicit none
    character*(*), intent(in) :: sub
    character*(*), intent(in) :: variable
    character(len=10)  :: sub_e='e_alloca'
    character(len=dbglinelen)  :: string=' ' 
    integer(ik4), optional, intent(in) :: line

    IF (present(line)) THEN 
       write(ndebug+node,'(A20,'':at line'',I6,A24,A)') & 
            sub(1:20),line,'error allocating ', variable
       write(string,'(A10,'':at line'',I6,A24,A20)') & 
            sub(1:20),line,'error allocating ', variable(1:20)
    ELSE
       write(ndebug+node,'(A20,'':'',A24,A20)') &
            sub(1:20),'error allocating: ',variable(1:20)
       write(string,'(A20,'':'',A24,A20)') &
            sub(1:20),'error allocating: ',variable(1:20)
    END IF
    call flusha(ndebug+node)
    !! if an error occours stops the program... 
    call error_prot(sub_e,3000,string)
  end subroutine e_alloca


  !!>
  !! ROUTINE : CHECK_MEMORY_R
  !! ACTION : check the address of real variable
  !!<
  subroutine check_memory_r(variable,name_of_variable,sub,line)
    implicit none
    character*(*), optional, intent(in) :: name_of_variable
    character*(*), optional, intent(in) :: sub
    real(rk8),            intent(in) :: variable
    integer(ik4),       optional, intent(in) :: line
    character(len=10) :: varia
    character(len=20) :: subro
    character(len=7) :: sline
    character(len=36) :: string

    integer(ik8) i_addr
    string='    '
    subro= '    '
    varia= '    '
    string='    '
    sline ='    '

    i_addr=loc(variable)

    write(string,"(A20,I15)") 'CM_mes>address is =' ,i_addr
    IF (present(name_of_variable)) varia=name_of_variable 
    IF (present(sub)) subro=sub
    IF (present(line))  write(SLINE,'(1X,I6)') LINE
    write(ndebug+node,"(A36,'for ',A10, ' in sub ',A20,'at line=',A7)") &
         string,varia,subro,sline

    call flusha(ndebug+node) 

  end subroutine check_memory_r

  !!>
  !!   ROUTINE : CHECK_MEMORY_I
  !!   ACTION :check the address of real variable
  !!           It may use the loc function (not architecture-universal)
  !!<
  subroutine check_memory_i(variable,name_of_variable,sub,line)
    implicit none
    character*(*), optional, intent(in) :: name_of_variable
    character*(*), optional, intent(in) :: sub
    integer(ik4),            intent(in) :: variable
    integer(ik4),       optional, intent(in) :: line
    character(len=10) :: varia
    character(len=20) :: subro
    character(len=7) :: sline
    character(len=36) :: string

    integer(ik8) i_addr
    string='    '
    subro= '    '
    varia= '    '
    string='    '
    sline ='    '

    i_addr=loc(variable)

    write(string,"(A20,I15)") 'CM_mes>address is =' ,i_addr
    IF (present(name_of_variable)) varia=name_of_variable 
    IF (present(sub)) subro=sub
    IF (present(line))  write(SLINE,'(1X,I6)') LINE
    write(ndebug+node,"(A36,'for ',A10, ' in sub ',A20,'at line=',A7)") &
         string,varia,subro,sline

    call flusha(ndebug+node) 

  end subroutine check_memory_i


  !!>
  !!   ROUTINE : LEN_STRIM
  !!   ACTION : trimming of string
  !!< 
  FUNCTION len_strim (string) RESULT (len_trim_RESULT)
    implicit none 
    character  (len=*), intent(IN) :: string
    integer(ik4) :: len_trim_RESULT,k

    len_trim_RESULT =0 
    DO k= LEN(string),1,-1
       IF (string(k:k) /=' ') THEN
          len_trim_RESULT =k
          EXIT
       END IF
    END DO

  END FUNCTION len_strim
 
    double precision function timer()
      implicit none
      integer(ik8) :: c , r , m
      call system_clock(count=c,count_rate=r,count_max=m)
      timer= dble(c)/dble(r)
    end function

    subroutine flusha(lunit)
      implicit none
      integer(ik4) , intent(in) :: lunit
      ! If whe have a FLUSH, use it 
      ! On IBM, flush is flush_ 
#ifdef IBM 
      call flush_(lunit)  
#else
      call flush(lunit)
#endif        
    end subroutine flusha

#endif

end module  mod_service

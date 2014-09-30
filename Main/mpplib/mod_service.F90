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
  use mod_mpmessage
  use mod_stdio
  use mod_mppparam
  use mod_dynparam , only : mycomm , myid , nproc , debug_level
  use mpi

  implicit none

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

  integer(ik4) , parameter :: maxnsubs = 100
  type(timing_info) , dimension(maxnsubs) :: info_serial , info_comm
  real(rk8) , dimension(maxnsubs) :: time_et , time_bt
  integer(ik4) :: n_of_nsubs = 0

  character(len=120) :: errmsg   !! a string where to compose an error message

  !! some global variable for debugging purposes
  !! set by prepare_debug, start_debug and stop_debug subroutines

  integer(ik4) :: nlevel = 0  !! level of depth in printing calling tree
  logical :: ldebug = .false. !! if true debug is enabled

  integer(ik4) :: node = 0
  integer(ik4) :: mxnode = 1
  integer(ik4) , allocatable , dimension(:) :: a_tmp
  character(len=dbgslen) :: stringa

  public :: activate_debug , start_debug , stop_debug
  public :: time_begin , time_end , time_reset , time_print

  contains

    subroutine activate_debug(level)
      implicit none
      integer(ik4) , optional :: level
      character(len=3) :: np = '   '
      character(len=9) :: string
      character(len=dbgslen) :: sub = 'activate_debug'
      integer(ik4) :: idum , iretval

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
        if (node.LT.10) then
          np ='0'
          write(np(2:2),'(i1)') node
        else
          write(np(1:2),'(i2)') node
        end if
      else
        if ( node < 10 ) then
          np ='00'
          write(np(3:3),'(i1)') node
        else if ( node < 100 ) then
          np ='0'
          write(np(2:3),'(i2)') node
        else
          write(np,'(i3)') node
        end if
      end if
      write(string,'(a6,a3)') "debug_",np
      open(ndebug+node,file=string,status='replace', &
              action='write',form='formatted',iostat=iretval)
      if ( iretval /= 0 ) then
        call fatal(__FILE__,__LINE__,'Cannot open debug files!')
      end if
      idum = ndebug+node
      write(ndebug+node,'(A20,'':'',A,i10)') &
          sub(1:20), 'DEBUGGING FILE CORRECTLY OPEN: unit is ', idum
      if ( present(level) ) debug_level = level
      write(ndebug+node,'(A20,'':'',A,i10)') &
          sub(1:20), 'Default debug_level is ', debug_level
      ldebug = .true.
      stringa(1:dbgslen) = ' '
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
      if ( indx == 0 ) then
        n_of_nsubs = n_of_nsubs + 1
        indx = n_of_nsubs
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

    subroutine time_end(name_of_section,indx,isize)
      implicit none
      character(len=dbgslen) :: name_of_section
      integer(ik4) , intent(in) :: indx
      integer(ik4) , optional :: isize
      real(rk8) :: time_call

      ! check for indx: should not be less than zero
      if ( indx < 0 ) then
        call fatal(__FILE__,__LINE__,'indx less then 0!!!')
      end if

      ! debugging stuff
      if ( ldebug .and. (debug_level > 3) ) then
        if ( nlevel <= 4 ) then
          write(ndebug+node,*) &
               stringa(1:nlevel*2),name_of_section(1:20),"(out)",nlevel
          call flusha(ndebug+node)
        end if
        nlevel = nlevel-1
      end if
      !
      time_et(indx) = timer()
      time_call = time_et(indx)-time_bt(indx)
      if ( present(isize) ) then
        info_comm(indx)%total_size = info_comm(indx)%total_size+isize
        info_comm(indx)%name_of_section = name_of_section
        info_comm(indx)%n_of_time = info_comm(indx)%n_of_time+1
        info_comm(indx)%total_time = info_comm(indx)%total_time+time_call
      else
        info_serial(indx)%name_of_section = name_of_section
        info_serial(indx)%n_of_time = info_serial(indx)%n_of_time+1
        info_serial(indx)%total_time = info_serial(indx)%total_time+time_call
      end if
    end subroutine time_end

    subroutine time_print(iunit,name_of_section)
      implicit none
      character(len=*) , optional :: name_of_section
      integer(ik4) :: iunit
      integer(ik4) :: nsubs , imin , imax , i , test , ilen , ierr
      real(rk8) :: avg , xmin , xmax
      real(rk8) , allocatable :: array_tmp(:)
      integer(ik4) , allocatable :: array_entries(:)
      logical :: l_times_on_pe = .false.
      logical :: l_nsubs = .true.
      real(rk8) :: total_comm_time = 0.0D0
      real(rk8) :: avg_value = 0.0D0
      character(len=128) :: cname
      character(len=dbgslen) :: sub='time_print'

      call mpi_barrier(mycomm,ierr)

      l_nsubs=.TRUE.
      if ( myid == italk ) then
        if ( present(name_of_section) )  then
          ilen = len_strim(name_of_section)
          write(iunit, &
             "(/,1x,5('!'),' Specific TIMING for section: ',a30,1x,5('!'))") &
             name_of_section(1:ilen)
        else
           write(iunit, &
             "(/,1x,5('!'),' Specific TIMING up to checkpoint ',25('!'))")
        end if
        write(iunit,105) 'section','times','avg-time','max(PE)','min(PE)'
      end if
105   format(1x,'!',a7,9x,a10,1x,a9,1x,a15,1x,a10)

      allocate(array_tmp(0:mxnode-1))
      allocate(array_entries(0:mxnode-1))

      ! check if the calling tree is the same on all nodes
      call allgather_i(array_entries,n_of_nsubs)
      test = array_entries(0)
      do i = 1 , mxnode-1
        if ( array_entries(i) /= test ) then
          write(ndebug+node,*) 'Warning:Different trees on different pe:',  &
               n_of_nsubs
          call fatal(__FILE__,__LINE__,'different trees on different pe!')
        end if
      end do

      ! if the calling tree is the same print out gathered data on OUTPUT file
      if ( l_nsubs ) then
        do nsubs = 1 , n_of_nsubs
          l_times_on_pe = .false.
          call allgather_i(a_tmp,info_serial(nsubs)%n_of_time)
          ! check the number of time
          test = a_tmp(0)
          avg_value = sum(a_tmp)/mxnode
          do i = 0 , mxnode-1
            if ( a_tmp(i) /= test ) then
              l_times_on_pe = .true.
            end if
            call MPI_barrier(mycomm,ierr)
          end do
          ! set to zero times less then 0.1 microseconds
          if ( info_serial(nsubs)%total_time <= 0.0000001D0 ) &
               info_serial(nsubs)%total_time = 0.0000001D0

          if ( info_serial(nsubs)%n_of_time >= 1 ) then
            call allgather_r(array_tmp,info_serial(nsubs)%total_time)
            call av_max_min(array_tmp,avg,xmax,imax,xmin,imin)
            ! compute avg, max and min
            if ( myid == italk ) then
              cname = info_serial(nsubs)%name_of_section
              if ( l_times_on_pe ) then
                cname = "*" // info_serial(nsubs)%name_of_section
              end if
              write(iunit,100) cname(1:15), &
                     info_serial(nsubs)%n_of_time, &
                     avg,xmax,imax-1,xmin,imin-1
            end if
          end if
        end do
      else
        if ( myid == italk ) then
          write (iunit,*) &
             'No global times for serial routines: different calling tree'
        end if
      end if

      !save anyway the data for each pe on debug unit:
      if ( ldebug ) then
        write(ndebug+node,'(A20,'':'',A,i10)') &
          sub(1:20), 'Specific local time up to checkpoint for node', node
        do nsubs = 1 , n_of_nsubs
          cname = info_serial(nsubs)%name_of_section
          if ( l_times_on_pe ) then
            cname = "*" // info_serial(nsubs)%name_of_section
          end if
          if ( info_serial(nsubs)%n_of_time >= 1 ) then
            write(ndebug+node,102) cname, &
               info_serial(nsubs)%n_of_time, &
               info_serial(nsubs)%total_time
            call flusha(ndebug+node)
          end if
        end do
      endif
      if ( myid == italk ) then
        write(iunit,"(1x,'!',19x,a30,/)") 'Communication subroutines :   '
        write(iunit,*) &
         '! section       times avg-time  max(PE)  min(PE) data  ratio(Mb/sec)'
      end if
      do nsubs = 1 , n_of_nsubs
        ! set to zero times less then 0.1 microseconds
        if ( info_comm(nsubs)%total_time<=0.0000001D0 ) &
              info_comm(nsubs)%total_time=0.0000001D0
        if ( info_comm(nsubs)%n_of_time >= 1 ) then
          call allgather_r(array_tmp,info_comm(nsubs)%total_time)
          call av_max_min(array_tmp,avg,xmax,imax,xmin,imin)
          if ( myid == italk ) then
            write(iunit,101) info_comm(nsubs)%name_of_section, &
                info_comm(nsubs)%n_of_time, &
                avg,xmax,imax-1,xmin,imin-1, &
                info_comm(nsubs)%total_size, &
                info_comm(nsubs)%total_size/(info_comm(nsubs)%total_time*1D6)
            total_comm_time = total_comm_time+avg
          end if
        end if
      end do
      if ( myid == italk ) then
        write(iunit,"(1x,'!',a40,3x,f9.4)") &
            'total communication time:', total_comm_time
        total_comm_time = 0.0D0
      end if
      if ( myid == italk ) then
        write(iunit,"(/,1x,10('!'),' End of Specific TIMING ', 37('!'),/)")
        call flusha(iunit)
      ENDIF
      deallocate(array_tmp)
      deallocate(array_entries)
100 format(1x,'!',a15,1x,i10,1x,f9.4,1x,f9.4,'(',i3,')',1x,f9.4,'(',i3,')')
101 format(1x,'!',a15,1x,i10,1x,f9.4,1x,f9.4,'(',i3,')', &
           1x,f9.4,'(',i3,')',i12,1x,f13.3)
102 format(1x,'!',a15,1x,i10,1x,f9.4,1x,f9.4)
    end subroutine time_print

    subroutine time_reset
      implicit none
      integer(ik4) :: nsubs
      do nsubs = 1 , maxnsubs
        info_serial(nsubs)%n_of_time = 0
        info_serial(nsubs)%total_time = 0.0D0
        info_serial(nsubs)%total_size = 0
        info_comm(nsubs)%n_of_time = 0
        info_comm(nsubs)%total_time = 0.0D0
        info_comm(nsubs)%total_size = 0
      end do
    end subroutine time_reset

    subroutine av_max_min(array,avg,xmax,indx_max,xmin,indx_min)
      implicit none
      real(rk8) , dimension(:) , intent(in) :: array
      integer(ik4) , intent(out) :: indx_min , indx_max
      real(rk8) , intent(out) :: xmax , xmin , avg
      integer(ik4) :: i , n_elements
      xmax = 0.0D0
      xmin = 1.D6
      avg = 0.0D0
      indx_max = 0
      indx_min = 0
      n_elements = size(array)
      do i = 1 , n_elements
        avg = avg+array(i)
        if ( array(i) > xmax ) then
          xmax = array(i)
          indx_max = i
        end if
        if ( array(i) < xmin ) then
          xmin = array(i)
          indx_min = i
        end if
      end do
      avg = avg/dble(n_elements)
    end subroutine av_max_min

    integer(ik4) function len_strim (string) result (len_trim_result)
      implicit none
      character (len=*) , intent(in) :: string
      integer(ik4) :: k
      len_trim_result = 0
      do k = len(string) , 1 , -1
        if ( string(k:k) /=' ' ) then
          len_trim_result = k
          exit
        end if
      end do
    end function len_strim

    real(rk8) function timer()
      implicit none
      integer(ik8) :: c , r , m
      call system_clock(count=c,count_rate=r,count_max=m)
      timer = dble(c)/dble(r)
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
#else
    character(len=4) , public :: unised_module
#endif

end module  mod_service

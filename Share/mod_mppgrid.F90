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

module mod_mppgrid

  use mpi
  use mod_realkinds
  use mod_message
  use mod_stdio
  use mod_memutil

  private

  logical :: mpi_debug = .false.
  integer :: mpilogunit = 233
  character(len=256) :: mpilogfile
  integer , parameter :: talkproc = 0
  integer , parameter :: masterproc = 0  ! Do not touch this
  integer , parameter :: mindimsize = 3  ! Do not touch this
  integer , parameter :: minpatchpoints = mindimsize*mindimsize

  logical , parameter , public :: global_grid = .true.
  logical , parameter , public :: proc_grid = .false.
  logical , parameter , public :: exchange_grid = .true.
  logical , parameter , public :: local_grid = .false.
  logical , parameter , public :: dot_grid = .true.
  logical , parameter , public :: cross_grid = .false.

  logical :: is_setup = .false.

  type model_processors
    integer :: total_cpus = 1
    integer , dimension(2) :: cpus_per_dim = (/1,1/)
    logical , dimension(2) :: dim_period = (/.false.,.false./)
  end type model_processors

  ! A domain is the global space of the 3D model.
  type model_domain
    integer :: g_i = mindimsize ! Global NS dimension
    integer :: g_j = mindimsize ! Global EW dimension
    integer :: totalpoints ! Product of the above, calculate once
    integer :: global_communicator = mpi_comm_world
    integer :: global_rank = -1
  end type model_domain

  type masternode
    integer , dimension(:) , pointer :: pgis
    integer , dimension(:) , pointer :: pgjs
    integer , dimension(:) , pointer :: pgie
    integer , dimension(:) , pointer :: pgje
    integer , dimension(:) , pointer :: pgsize
    real(sp) , pointer , dimension(:) :: excbuf2dr
    real(dp) , pointer , dimension(:) :: excbuf2dd
    integer , pointer , dimension(:) :: excbuf2di
    integer(2) , pointer , dimension(:) :: excbuf2ds
    logical , pointer , dimension(:) :: excbuf2dl
  end type masternode

  ! The part this processor is going to compute the tendencies on
  type processor_domain
    integer :: p_i = mindimsize ! Local to processor i points
    integer :: p_j = mindimsize ! Local to processor j points
    integer :: g_i1 = 1
    integer :: g_j1 = 1
    integer :: totalpoints ! Product of the above, calculate once
    integer :: cartesian_communicator = mpi_comm_world
    integer :: cartesian_rank = -1
    integer , dimension(2) :: location
    integer :: top = mpi_proc_null
    integer :: lht = mpi_proc_null
    integer :: rht = mpi_proc_null
    integer :: btm = mpi_proc_null
    integer :: rhb = mpi_proc_null
    integer :: lhb = mpi_proc_null
    integer :: lhs = mpi_proc_null
    integer :: rhs = mpi_proc_null
  end type processor_domain

  type procbounds
    integer :: icompstart
    integer :: icompend
    integer :: icompnp
    integer :: jcompstart
    integer :: jcompend
    integer :: jcompnp
  end type procbounds

  integer :: csize                                                          
  integer :: gsize                                                          
  real(sp) , dimension(1) :: sndcrnpr
  real(sp) , dimension(1) :: rcvcrnpr
  real(sp) , pointer , dimension(:) :: sndbuf1dr
  real(sp) , pointer , dimension(:) :: rcvbuf1dr
  real(sp) , pointer , dimension(:) :: excbuf2dr
  real(dp) , dimension(1) :: sndcrnpd
  real(dp) , dimension(1) :: rcvcrnpd
  real(dp) , pointer , dimension(:) :: sndbuf1dd
  real(dp) , pointer , dimension(:) :: rcvbuf1dd
  real(dp) , pointer , dimension(:) :: excbuf2dd
  integer , dimension(1) :: sndcrnpi
  integer , dimension(1) :: rcvcrnpi
  integer , pointer , dimension(:) :: sndbuf1di
  integer , pointer , dimension(:) :: rcvbuf1di
  integer , pointer , dimension(:) :: excbuf2di
  integer(2) , dimension(1) :: sndcrnps
  integer(2) , dimension(1) :: rcvcrnps
  integer(2) , pointer , dimension(:) :: sndbuf1ds
  integer(2) , pointer , dimension(:) :: rcvbuf1ds
  integer(2) , pointer , dimension(:) :: excbuf2ds
  logical , dimension(1) :: sndcrnpl
  logical , dimension(1) :: rcvcrnpl
  logical , pointer , dimension(:) :: sndbuf1dl
  logical , pointer , dimension(:) :: rcvbuf1dl
  logical , pointer , dimension(:) :: excbuf2dl

  type (model_processors) :: xproc
  type (model_domain) :: gspace
  type (processor_domain) :: pspace
  type (masternode) :: mnode
  type (procbounds) , protected :: pbnds

  interface getgrid
    module procedure getgrid2d_d
    module procedure getgrid2d_r
    module procedure getgrid2d_i
    module procedure getgrid2d_s
    module procedure getgrid2d_l
    module procedure getgrid3d_d
    module procedure getgrid3d_r
    module procedure getgrid3d_i
    module procedure getgrid3d_s
    module procedure getgrid3d_l
    module procedure getgrid4d_d
    module procedure getgrid4d_r
    module procedure getgrid4d_i
    module procedure getgrid4d_s
    module procedure getgrid4d_l
  end interface getgrid

  interface exchange_internal
    module procedure exchange_internal2d_d
    module procedure exchange_internal2d_r
    module procedure exchange_internal2d_i
    module procedure exchange_internal2d_s
    module procedure exchange_internal2d_l
    module procedure exchange_internal3d_d
    module procedure exchange_internal3d_r
    module procedure exchange_internal3d_i
    module procedure exchange_internal3d_s
    module procedure exchange_internal3d_l
    module procedure exchange_internal4d_d
    module procedure exchange_internal4d_r
    module procedure exchange_internal4d_i
    module procedure exchange_internal4d_s
    module procedure exchange_internal4d_l
  end interface exchange_internal

  interface global_to_proc
    module procedure global_to_proc2d_d
    module procedure global_to_proc2d_r
    module procedure global_to_proc2d_i
    module procedure global_to_proc2d_s
    module procedure global_to_proc2d_l
    module procedure global_to_proc3d_d
    module procedure global_to_proc3d_r
    module procedure global_to_proc3d_i
    module procedure global_to_proc3d_s
    module procedure global_to_proc3d_l
    module procedure global_to_proc4d_d
    module procedure global_to_proc4d_r
    module procedure global_to_proc4d_i
    module procedure global_to_proc4d_s
    module procedure global_to_proc4d_l
  end interface global_to_proc

  interface proc_to_global
    module procedure proc_to_global2d_d
    module procedure proc_to_global2d_r
    module procedure proc_to_global2d_i
    module procedure proc_to_global2d_s
    module procedure proc_to_global2d_l
    module procedure proc_to_global3d_d
    module procedure proc_to_global3d_r
    module procedure proc_to_global3d_i
    module procedure proc_to_global3d_s
    module procedure proc_to_global3d_l
    module procedure proc_to_global4d_d
    module procedure proc_to_global4d_r
    module procedure proc_to_global4d_i
    module procedure proc_to_global4d_s
    module procedure proc_to_global4d_l
  end interface proc_to_global

  interface master_to_nodes
    module procedure master_to_nodes_d
    module procedure master_to_nodes_r
    module procedure master_to_nodes_i
    module procedure master_to_nodes_s
    module procedure master_to_nodes_l
    module procedure master_to_nodes1d_d
    module procedure master_to_nodes1d_r
    module procedure master_to_nodes1d_i
    module procedure master_to_nodes1d_s
    module procedure master_to_nodes1d_l
    module procedure master_to_nodes2d_d
    module procedure master_to_nodes2d_r
    module procedure master_to_nodes2d_i
    module procedure master_to_nodes2d_s
    module procedure master_to_nodes2d_l
    module procedure master_to_nodes3d_d
    module procedure master_to_nodes3d_r
    module procedure master_to_nodes3d_i
    module procedure master_to_nodes3d_s
    module procedure master_to_nodes3d_l
    module procedure master_to_nodes4d_d
    module procedure master_to_nodes4d_r
    module procedure master_to_nodes4d_i
    module procedure master_to_nodes4d_s
    module procedure master_to_nodes4d_l
  end interface master_to_nodes

  interface nodes_to_master
    module procedure nodes_to_master2d_d
    module procedure nodes_to_master2d_r
    module procedure nodes_to_master2d_i
    module procedure nodes_to_master2d_s
    module procedure nodes_to_master2d_l
    module procedure nodes_to_master3d_d
    module procedure nodes_to_master3d_r
    module procedure nodes_to_master3d_i
    module procedure nodes_to_master3d_s
    module procedure nodes_to_master3d_l
    module procedure nodes_to_master4d_d
    module procedure nodes_to_master4d_r
    module procedure nodes_to_master4d_i
    module procedure nodes_to_master4d_s
    module procedure nodes_to_master4d_l
  end interface nodes_to_master

  interface global_sum
    module procedure global_sum_d
    module procedure global_sum_r
    module procedure global_sum_i
  end interface global_sum

  public :: am_i_master , cantalk , toggle_mpi_debug
  public :: setup_domain , delete_domain
  public :: getgrid
  public :: exchange_internal
  public :: global_to_proc , proc_to_global
  public :: master_to_nodes , nodes_to_master
  public :: global_sum

  contains

    subroutine toggle_mpi_debug(logpath)
      implicit none
      character(len=*) , intent(in) , optional :: logpath
      integer :: ierr
      if ( .not. mpi_debug ) then
        if ( present(logpath) ) then
          write(mpilogfile,'(a,a,a,i0.4,a)') trim(logpath), '/', 'mpilog_', &
                    gspace%global_rank+1, '.txt'
        else
          write(mpilogfile,'(a,i0.4,a)') './mpilog_', &
                    gspace%global_rank+1, '.txt'
        end if
        open(mpilogunit,file=mpilogfile,status='replace',action='write', &
             iostat=ierr)
        if ( ierr /= 0 ) mpilogunit = stderr
        write(mpilogunit,'(a,i0)') 'Total CPUS      ', xproc%total_cpus
        write(mpilogunit,'(a,i0)') 'Global Rank     ', gspace%global_rank
        write(mpilogunit,'(a,i0)') 'Global Start I  ', pspace%g_i1
        write(mpilogunit,'(a,i0)') 'Global Start J  ', pspace%g_j1
        write(mpilogunit,'(a,i0)') 'Mine number I   ', pspace%p_i
        write(mpilogunit,'(a,i0)') 'Mine number J   ', pspace%p_j
        write(mpilogunit,'(a,i0)') 'BottomLeft CPU  ', pspace%lhb
        write(mpilogunit,'(a,i0)') 'Left CPU        ', pspace%lhs
        write(mpilogunit,'(a,i0)') 'TopLeft CPU     ', pspace%lht
        write(mpilogunit,'(a,i0)') 'Top CPU         ', pspace%top
        write(mpilogunit,'(a,i0)') 'TopRight CPU    ', pspace%rht
        write(mpilogunit,'(a,i0)') 'Right CPU       ', pspace%rhs
        write(mpilogunit,'(a,i0)') 'BottomRight CPU ', pspace%rhb
        write(mpilogunit,'(a,i0)') 'Bottom CPU      ', pspace%btm
        call flush(mpilogunit)
      else
        if ( mpilogunit /= stderr ) close(mpilogunit)
      end if
      mpi_debug = .not. mpi_debug
    end subroutine toggle_mpi_debug

    logical function cantalk( )
      implicit none
      cantalk = (gspace%global_rank == talkproc)
    end function cantalk

    logical function am_i_master( )
      implicit none
      am_i_master = (gspace%global_rank == masterproc)
    end function am_i_master

    subroutine setup_domain(ni,nj,iband,comm)
      implicit none
      integer , intent(in) :: ni , nj , iband
      integer , intent(in) , optional :: comm
      integer :: ierr , jbi , jbj , imaxcpus , max_pi , max_pj , max_p
      integer :: inode , maxec , maxcpu
      integer , dimension(2) :: search_coord
      real :: dimfac
      if ( present(comm) ) then
        gspace%global_communicator = comm
      end if
      call mpi_comm_size(gspace%global_communicator,xproc%total_cpus,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_comm_rank(gspace%global_communicator,gspace%global_rank,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( ni < mindimsize .or. nj < mindimsize ) then
        if ( cantalk( ) ) then
          write(stderr,*) 'NI :', ni
          write(stderr,*) 'NJ :', nj
          write(stderr,*) 'A Dimension cannot be less than ', mindimsize
        end if
        call fatal(__FILE__,__LINE__,'Dimension error')
      end if
      gspace%g_i = ni
      gspace%g_j = nj
      gspace%totalpoints = gspace%g_i*gspace%g_j
      if ( gspace%totalpoints/xproc%total_cpus < minpatchpoints ) then
        if ( cantalk() ) then
          write(stderr,*) 'Have too much workers and too little a job...'
          write(stderr,*) 'Total Points : ', gspace%totalpoints
          write(stderr,*) 'Total CPUS   : ', xproc%total_cpus
          imaxcpus = gspace%totalpoints/minpatchpoints
          if ( imaxcpus < 1 ) imaxcpus = 1
          write(stderr,*) 'Suggested number of CPUS : ', imaxcpus
        end if
        call fatal(__FILE__,__LINE__,'Too much processors')
      end if
      ! Here we divide the work among all cpus
      if ( xproc%total_cpus > 1 ) then
        ! Some trivial cases
        if ( xproc%total_cpus == 2 ) then
          xproc%cpus_per_dim(1) = 2
          xproc%cpus_per_dim(2) = 1
        else if ( xproc%total_cpus == 3 ) then
          if ( gspace%g_i > gspace%g_j ) then
            xproc%cpus_per_dim(1) = 1
            xproc%cpus_per_dim(2) = 3
          else
            xproc%cpus_per_dim(1) = 3
            xproc%cpus_per_dim(2) = 1
          end if
        else
          dimfac = real(gspace%g_j)/real(gspace%g_i+real(gspace%g_j))
          xproc%cpus_per_dim(1) = int(real(xproc%total_cpus)**dimfac)
          xproc%cpus_per_dim(2) = xproc%total_cpus / xproc%cpus_per_dim(1)
        end if
      end if
      imaxcpus = xproc%cpus_per_dim(1)*xproc%cpus_per_dim(2)
      if ( mod(xproc%total_cpus,imaxcpus) > 0 ) then
        if ( cantalk() ) then
          write(stderr,*) 'Work does not evenly divide.'
          write(stderr,*) 'Suggested number of CPUS : ', imaxcpus
        end if
        call fatal(__FILE__,__LINE__,'Scheme not working')
      end if
      if ( iband == 1 ) then
        xproc%dim_period(1) = .true.
      end if
      call mpi_cart_create(gspace%global_communicator,2,               &
                           xproc%cpus_per_dim,xproc%dim_period,.true., &
                           pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_comm_rank(pspace%cartesian_communicator, &
                         pspace%cartesian_rank,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_cart_coords(pspace%cartesian_communicator, &
                           pspace%cartesian_rank,2,pspace%location,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      jbj = gspace%g_j / xproc%cpus_per_dim(1)
      jbi = gspace%g_i / xproc%cpus_per_dim(2)
      pspace%g_j1 = pspace%location(1) * jbj + 1
      pspace%g_i1 = pspace%location(2) * jbi + 1
      max_pj = xproc%cpus_per_dim(1) - 1
      max_pi = xproc%cpus_per_dim(2) - 1
      if ( pspace%location(1) == max_pj ) then
        pspace%p_j = jbj + mod(gspace%g_j,jbj)
      else
        pspace%p_j = jbj
      end if
      if ( pspace%location(2) == max_pi ) then
        pspace%p_i = jbi + mod(gspace%g_i,jbi)
      else
        pspace%p_i = jbi
      end if
      pspace%totalpoints = pspace%p_i*pspace%p_j
      pbnds%icompnp = pspace%p_i
      pbnds%jcompnp = pspace%p_j
      if ( xproc%total_cpus > 1 ) then
        ! Search neigh processors ranks
        if ( pspace%location(2) < max_pi ) then
          search_coord(1) = pspace%location(1)
          search_coord(2) = pspace%location(2) + 1
          call mpi_cart_rank(pspace%cartesian_communicator,search_coord, &
                             pspace%top,ierr)
          pbnds%icompend = pspace%p_i
        else
          pbnds%icompend = pspace%p_i-1
        end if
        if ( pspace%location(2) < max_pi .and. &
            (pspace%location(1) > 0 .or. iband == 1) ) then
          search_coord(1) = pspace%location(1) - 1
          search_coord(2) = pspace%location(2) + 1
          call mpi_cart_rank(pspace%cartesian_communicator,search_coord, &
                             pspace%lht,ierr)
        end if
        if ( pspace%location(2) < max_pi .and. &
            (pspace%location(1) < max_pj .or. iband == 1) ) then
          search_coord(1) = pspace%location(1) + 1
          search_coord(2) = pspace%location(2) + 1
          call mpi_cart_rank(pspace%cartesian_communicator,search_coord, &
                             pspace%rht,ierr)
        end if
        if ( pspace%location(1) < max_pj .or. iband == 1 ) then
          search_coord(1) = pspace%location(1) + 1
          search_coord(2) = pspace%location(2)
          call mpi_cart_rank(pspace%cartesian_communicator,search_coord, &
                             pspace%rhs,ierr)
          pbnds%jcompend = pspace%p_j
        else
          pbnds%jcompend = pspace%p_j-1
        end if
        if ( pspace%location(1) > 0 .or. iband == 1 ) then
          search_coord(1) = pspace%location(1) - 1
          search_coord(2) = pspace%location(2)
          call mpi_cart_rank(pspace%cartesian_communicator,search_coord, &
                             pspace%lhs,ierr)
          pbnds%jcompstart = 1
        else
          pbnds%jcompstart = 2
        end if
        if ( pspace%location(2) > 0 ) then
          search_coord(1) = pspace%location(1)
          search_coord(2) = pspace%location(2) - 1
          call mpi_cart_rank(pspace%cartesian_communicator,search_coord, &
                             pspace%btm,ierr)
          pbnds%icompstart = 1
        else
          pbnds%icompstart = 2
        end if
        if ( pspace%location(2) > 0 .and. &
            (pspace%location(1) < max_pj .or. iband == 1) ) then
          search_coord(1) = pspace%location(1) + 1
          search_coord(2) = pspace%location(2) - 1
          call mpi_cart_rank(pspace%cartesian_communicator,search_coord, &
                           pspace%rhb,ierr)
        end if
        if ( pspace%location(2) > 0 .and. &
            (pspace%location(1) > 0 .or. iband == 1) ) then
          search_coord(1) = pspace%location(1) - 1
          search_coord(2) = pspace%location(2) - 1
          call mpi_cart_rank(pspace%cartesian_communicator,search_coord, &
                             pspace%lhb,ierr)
        end if
      else
        if ( iband == 1 ) then
          pspace%rhs = pspace%cartesian_rank
          pspace%lhs = pspace%cartesian_rank
          pbnds%jcompstart = 1
          pbnds%jcompend = pspace%p_j
        else
          pbnds%jcompstart = 2
          pbnds%jcompend = pspace%p_j-1
        end if
        pbnds%icompstart = 2
        pbnds%icompend = pspace%p_i-1
      end if
      if ( cantalk( ) .and. .false. ) then
        write(stdout, *) '----------------------------------------------'
        write(stdout,'(a,i0,a,i0)') 'Domain decomp               : ', &
                        xproc%cpus_per_dim(1),' x ',xproc%cpus_per_dim(2)
        write(stdout,'(a,i0,a,i0)') 'Processor load points       : ', &
                        jbi ,' x ',jbj
        write(stdout,'(a,i0)') 'Corner North procs overhead : ', &
                        mod(gspace%g_i,jbi)
        write(stdout,'(a,i0)') 'Corner East procs overhead  : ', &
                        mod(gspace%g_j,jbj)
      end if
      max_p = max(pspace%p_i,pspace%p_j)+1
      csize = max_p
      gsize = gspace%g_i*gspace%g_j
      call getmem1d(sndbuf1dr,0,max_p,__FILE__)
      call getmem1d(rcvbuf1dr,0,max_p,__FILE__)
      call getmem1d(sndbuf1dd,0,max_p,__FILE__)
      call getmem1d(rcvbuf1dd,0,max_p,__FILE__)
      call getmem1d(sndbuf1di,0,max_p,__FILE__)
      call getmem1d(rcvbuf1di,0,max_p,__FILE__)
      call getmem1d(sndbuf1ds,0,max_p,__FILE__)
      call getmem1d(rcvbuf1ds,0,max_p,__FILE__)
      call getmem1d(sndbuf1dl,0,max_p,__FILE__)
      call getmem1d(rcvbuf1dl,0,max_p,__FILE__)
      if ( am_i_master( ) ) then
        maxcpu = xproc%total_cpus-1
        call getmem1d(mnode%pgis,0,maxcpu,__FILE__)
        call getmem1d(mnode%pgjs,0,maxcpu,__FILE__)
        call getmem1d(mnode%pgie,0,maxcpu,__FILE__)
        call getmem1d(mnode%pgje,0,maxcpu,__FILE__)
        call getmem1d(mnode%pgsize,0,maxcpu,__FILE__)
        do inode = 0 , maxcpu
          if ( inode /= masterproc ) then
            call mpi_recv(mnode%pgis(inode),1,mpi_integer,inode,0, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            call mpi_recv(mnode%pgjs(inode),1,mpi_integer,inode,0, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            call mpi_recv(mnode%pgie(inode),1,mpi_integer,inode,0, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            call mpi_recv(mnode%pgje(inode),1,mpi_integer,inode,0, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            call mpi_recv(mnode%pgsize(inode),1,mpi_integer,inode,0, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            mnode%pgie(inode) = mnode%pgie(inode)+mnode%pgis(inode)-1
            mnode%pgje(inode) = mnode%pgje(inode)+mnode%pgjs(inode)-1
          else
            mnode%pgis(inode) = pspace%g_i1
            mnode%pgjs(inode) = pspace%g_j1
            mnode%pgie(inode) = pspace%g_i1+pspace%p_i-1
            mnode%pgje(inode) = pspace%g_j1+pspace%p_j-1
            mnode%pgsize(inode) = pspace%totalpoints
          end if
          if ( cantalk( ) .and. .false. ) then
            write(stdout, *) '###################################'
            write(stdout,'(a,i0)') 'CPU number: ', inode
            write(stdout,'(a,i0)') 'Start   i : ', mnode%pgis(inode)
            write(stdout,'(a,i0)') 'End     i : ', mnode%pgie(inode)
            write(stdout,'(a,i0)') 'Start   j : ', mnode%pgjs(inode)
            write(stdout,'(a,i0)') 'End     j : ', mnode%pgje(inode)
            write(stdout, *) '###################################'
          end if
        end do
        maxec = maxval(mnode%pgsize)
        call getmem1d(mnode%excbuf2dr,1,maxec,__FILE__)
        call getmem1d(mnode%excbuf2dd,1,maxec,__FILE__)
        call getmem1d(mnode%excbuf2di,1,maxec,__FILE__)
        call getmem1d(mnode%excbuf2ds,1,maxec,__FILE__)
        call getmem1d(mnode%excbuf2dl,1,maxec,__FILE__)
      else
        call mpi_send(pspace%g_i1,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        call mpi_send(pspace%g_j1,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        call mpi_send(pspace%p_i,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        call mpi_send(pspace%p_j,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        call mpi_send(pspace%totalpoints,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        if ( cantalk( ) .and. .false. ) then
          write(stdout, *) '###################################'
          write(stdout,'(a,i0)') 'CPU number: ', gspace%global_rank
          write(stdout,'(a,i0)') 'Start   i : ', pspace%g_i1
          write(stdout,'(a,i0)') 'End     i : ', pspace%g_i1+pspace%p_i-1
          write(stdout,'(a,i0)') 'Start   j : ', pspace%g_j1
          write(stdout,'(a,i0)') 'End     j : ', pspace%g_j1+pspace%p_j-1
          write(stdout, *) '###################################'
        end if
      end if
      call getmem1d(excbuf2dr,1,pspace%totalpoints,__FILE__)
      call getmem1d(excbuf2dd,1,pspace%totalpoints,__FILE__)
      call getmem1d(excbuf2di,1,pspace%totalpoints,__FILE__)
      call getmem1d(excbuf2ds,1,pspace%totalpoints,__FILE__)
      call getmem1d(excbuf2dl,1,pspace%totalpoints,__FILE__)
      if ( cantalk( ) .and. .false. ) then
        write(stdout, *) '----------------------------------------------'
      end if
      is_setup = .true.
    end subroutine setup_domain

    subroutine delete_domain
      implicit none
      if ( .not. is_setup ) return
      call relmem1d(sndbuf1dr)
      call relmem1d(rcvbuf1dr)
      call relmem1d(sndbuf1dd)
      call relmem1d(rcvbuf1dd)
      call relmem1d(sndbuf1di)
      call relmem1d(rcvbuf1di)
      call relmem1d(sndbuf1ds)
      call relmem1d(rcvbuf1ds)
      call relmem1d(sndbuf1dl)
      call relmem1d(rcvbuf1dl)
      call relmem1d(excbuf2dr)
      call relmem1d(excbuf2dd)
      call relmem1d(excbuf2di)
      call relmem1d(excbuf2ds)
      call relmem1d(excbuf2dl)
      if ( am_i_master( ) ) then
        call relmem1d(mnode%pgis)
        call relmem1d(mnode%pgie)
        call relmem1d(mnode%pgjs)
        call relmem1d(mnode%pgje)
        call relmem1d(mnode%pgsize)
        call relmem1d(mnode%excbuf2dr)
        call relmem1d(mnode%excbuf2dd)
        call relmem1d(mnode%excbuf2di)
        call relmem1d(mnode%excbuf2ds)
        call relmem1d(mnode%excbuf2dl)
      end if
      is_setup = .false.
    end subroutine delete_domain

    subroutine getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      implicit none
      logical , intent(in) :: isglobal , isstagger , isexchange
      integer , intent(out) :: i1 , i2 , j1 , j2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling getgrid before domain_setup')
      end if
      if ( isglobal ) then
        i1 = 1
        i2 = gspace%g_i
        j1 = 1
        j2 = gspace%g_j
      else
        i1 = 1
        i2 = pspace%p_i
        j1 = 1
        j2 = pspace%p_j
        if ( isexchange ) then
          i1 = i1 - 1
          j1 = j1 - 1
          i2 = i2 + 1
          j2 = j2 + 1
        end if
      end if
      if ( isstagger ) then
        i2 = i2 + 1
        j2 = j2 + 1
      end if
    end subroutine getextrema

    subroutine getgrid2d_d(g,lglobal,lexchange,lstagger)
      implicit none
      real(dp) , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_d

    subroutine getgrid2d_r(g,lglobal,lexchange,lstagger)
      implicit none
      real(sp) , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_r

    subroutine getgrid2d_i(g,lglobal,lexchange,lstagger)
      implicit none
      integer , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_i

    subroutine getgrid2d_s(g,lglobal,lexchange,lstagger)
      implicit none
      integer(2) , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_s

    subroutine getgrid2d_l(g,lglobal,lexchange,lstagger)
      implicit none
      logical , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_l

    subroutine getgrid3d_d(g,k1,k2,lglobal,lexchange,lstagger)
      implicit none
      real(dp) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_d

    subroutine getgrid3d_r(g,k1,k2,lglobal,lexchange,lstagger)
      implicit none
      real(sp) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_r

    subroutine getgrid3d_i(g,k1,k2,lglobal,lexchange,lstagger)
      implicit none
      integer , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_i

    subroutine getgrid3d_s(g,k1,k2,lglobal,lexchange,lstagger)
      implicit none
      integer(2) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_s

    subroutine getgrid3d_l(g,k1,k2,lglobal,lexchange,lstagger)
      implicit none
      logical , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_l

    subroutine getgrid4d_d(g,k1,k2,t1,t2,lglobal,lexchange,lstagger)
      implicit none
      real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_d

    subroutine getgrid4d_r(g,k1,k2,t1,t2,lglobal,lexchange,lstagger)
      implicit none
      real(sp) , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_r

    subroutine getgrid4d_i(g,k1,k2,t1,t2,lglobal,lexchange,lstagger)
      implicit none
      integer , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_i

    subroutine getgrid4d_s(g,k1,k2,t1,t2,lglobal,lexchange,lstagger)
      implicit none
      integer(2) , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_s

    subroutine getgrid4d_l(g,k1,k2,t1,t2,lglobal,lexchange,lstagger)
      implicit none
      logical , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lexchange
      logical , intent(in) , optional :: lstagger
      logical :: isglobal = .false.
      logical :: isstagger = .false.
      logical :: isexchange = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lexchange) ) isexchange = lexchange
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(isglobal,isexchange,isstagger,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_l

    subroutine exchange_internal2d_d(l)
      implicit none
      real(dp) , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      njj = ubound(l,1) - 1
      nii = ubound(l,2) - 1
      if ( pspace%btm /= mpi_proc_null ) then
        sndbuf1dd(1:njj) = l(1:njj,1)
      end if
      call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%btm,1, &
                        rcvbuf1dd,csize,mpi_real8,pspace%top,1, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%top /= mpi_proc_null ) then
        l(1:njj,nii+1) = rcvbuf1dd(1:njj)
      end if
      if ( pspace%rhs /= mpi_proc_null ) then
        sndbuf1dd(1:nii) = l(njj,1:nii)
      end if
      call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%rhs,2, &
                        rcvbuf1dd,csize,mpi_real8,pspace%lhs,2, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhs /= mpi_proc_null ) then
        l(0,1:nii) = rcvbuf1dd(1:nii)
      end if
      if ( pspace%lhs /= mpi_proc_null ) then
        sndbuf1dd(1:nii) = l(1,1:nii)
      end if
      call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%lhs,3, &
                        rcvbuf1dd,csize,mpi_real8,pspace%rhs,3, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhs /= mpi_proc_null ) then
        l(njj+1,1:nii) = rcvbuf1dd(1:nii)
      end if
      if ( pspace%top /= mpi_proc_null ) then
        sndbuf1dd(1:njj) = l(1:njj,nii)
      end if
      call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%top,4, &
                        rcvbuf1dd,csize,mpi_real8,pspace%btm,4, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%btm /= mpi_proc_null ) then
        l(1:njj,0) = rcvbuf1dd(1:njj)
      end if
      if ( pspace%rht /= mpi_proc_null ) then
        sndcrnpd(1) = l(njj,nii)
      end if
      call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%rht,1, &
                        rcvcrnpd,1,mpi_real8,pspace%lhb,1, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        l(0,0) = rcvcrnpd(1)
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        sndcrnpd(1) = l(1,1)
      end if
      call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%lhb,2, &
                        rcvcrnpd,1,mpi_real8,pspace%rht,2, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        l(njj+1,nii+1) = rcvcrnpd(1)
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        sndcrnpd(1) = l(njj,1)
      end if
      call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%rhb,3, &
                        rcvcrnpd,1,mpi_real8,pspace%lht,3, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        l(0,nii+1) = rcvcrnpd(1)
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        sndcrnpd(1) = l(1,nii)
      end if
      call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%lht,4, &
                        rcvcrnpd,1,mpi_real8,pspace%rhb,4, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        l(njj+1,0) = rcvcrnpd(1)
      end if
    end subroutine exchange_internal2d_d

    subroutine exchange_internal2d_r(l)
      implicit none
      real(sp) , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      if ( pspace%btm /= mpi_proc_null ) then
        sndbuf1dr(1:njj) = l(1:njj,1)
      end if
      call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%btm,1, &
                        rcvbuf1dr,csize,mpi_real4,pspace%top,1, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%top /= mpi_proc_null ) then
        l(1:njj,nii+1) = rcvbuf1dr(1:njj)
      end if
      if ( pspace%rhs /= mpi_proc_null ) then
        sndbuf1dr(1:nii) = l(njj,1:nii)
      end if
      call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%rhs,2, &
                        rcvbuf1dr,csize,mpi_real4,pspace%lhs,2, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhs /= mpi_proc_null ) then
        l(0,1:nii) = rcvbuf1dr(1:nii)
      end if
      if ( pspace%lhs /= mpi_proc_null ) then
        sndbuf1dr(1:nii) = l(1,1:nii)
      end if
      call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%lhs,3, &
                        rcvbuf1dr,csize,mpi_real4,pspace%rhs,3, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhs /= mpi_proc_null ) then
        l(njj+1,1:nii) = rcvbuf1dr(1:nii)
      end if
      if ( pspace%top /= mpi_proc_null ) then
        sndbuf1dr(1:njj) = l(1:njj,nii)
      end if
      call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%top,4, &
                        rcvbuf1dr,csize,mpi_real4,pspace%btm,4, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%btm /= mpi_proc_null ) then
        l(1:njj,0) = rcvbuf1dr(1:njj)
      end if
      if ( pspace%rht /= mpi_proc_null ) then
        sndcrnpr(1) = l(njj,nii)
      end if
      call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%rht,1, &
                        rcvcrnpr,1,mpi_real4,pspace%lhb,1, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        l(0,0) = rcvcrnpr(1)
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        sndcrnpr(1) = l(1,1)
      end if
      call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%lhb,2, &
                        rcvcrnpr,1,mpi_real4,pspace%rht,2, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        l(njj+1,nii+1) = rcvcrnpr(1)
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        sndcrnpr(1) = l(njj,1)
      end if
      call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%rhb,3, &
                        rcvcrnpr,1,mpi_real4,pspace%lht,3, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        l(0,nii+1) = rcvcrnpr(1)
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        sndcrnpr(1) = l(1,nii)
      end if
      call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%lht,4, &
                        rcvcrnpr,1,mpi_real4,pspace%rhb,4, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        l(njj+1,0) = rcvcrnpr(1)
      end if
    end subroutine exchange_internal2d_r

    subroutine exchange_internal2d_i(l)
      implicit none
      integer , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      if ( pspace%btm /= mpi_proc_null ) then
        sndbuf1di(1:njj) = l(1:njj,1)
      end if
      call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%btm,1,&
                        rcvbuf1di,csize,mpi_integer,pspace%top,1,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%top /= mpi_proc_null ) then
        l(1:njj,nii+1) = rcvbuf1di(1:njj)
      end if
      if ( pspace%rhs /= mpi_proc_null ) then
        sndbuf1di(1:nii) = l(njj,1:nii)
      end if
      call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%rhs,2,&
                        rcvbuf1di,csize,mpi_integer,pspace%lhs,2,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhs /= mpi_proc_null ) then
        l(0,1:nii) = rcvbuf1di(1:nii)
      end if
      if ( pspace%lhs /= mpi_proc_null ) then
        sndbuf1di(1:nii) = l(1,1:nii)
      end if
      call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%lhs,3,&
                        rcvbuf1di,csize,mpi_integer,pspace%rhs,3,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhs /= mpi_proc_null ) then
        l(njj+1,1:nii) = rcvbuf1di(1:nii)
      end if
      if ( pspace%top /= mpi_proc_null ) then
        sndbuf1di(1:njj) = l(1:njj,nii)
      end if
      call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%top,4,&
                        rcvbuf1di,csize,mpi_integer,pspace%btm,4,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%btm /= mpi_proc_null ) then
        l(1:njj,0) = rcvbuf1di(1:njj)
      end if
      if ( pspace%rht /= mpi_proc_null ) then
        sndcrnpi(1) = l(njj,nii)
      end if
      call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%rht,1, &
                        rcvcrnpi,1,mpi_integer,pspace%lhb,1, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        l(0,0) = rcvcrnpi(1)
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        sndcrnpi(1) = l(1,1)
      end if
      call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%lhb,2, &
                        rcvcrnpi,1,mpi_integer,pspace%rht,2, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        l(njj+1,nii+1) = rcvcrnpi(1)
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        sndcrnpi(1) = l(njj,1)
      end if
      call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%rhb,3, &
                        rcvcrnpi,1,mpi_integer,pspace%lht,3, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        l(0,nii+1) = rcvcrnpi(1)
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        sndcrnpi(1) = l(1,nii)
      end if
      call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%lht,4, &
                        rcvcrnpi,1,mpi_integer,pspace%rhb,4, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        l(njj+1,0) = rcvcrnpi(1)
      end if
    end subroutine exchange_internal2d_i

    subroutine exchange_internal2d_s(l)
      implicit none
      integer(2) , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      if ( pspace%btm /= mpi_proc_null ) then
        sndbuf1ds(1:njj) = l(1:njj,1)
      end if
      call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%btm,1,&
                        rcvbuf1ds,csize,mpi_short,pspace%top,1,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%top /= mpi_proc_null ) then
        l(1:njj,nii+1) = rcvbuf1ds(1:njj)
      end if
      if ( pspace%rhs /= mpi_proc_null ) then
        sndbuf1ds(1:nii) = l(njj,1:nii)
      end if
      call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%rhs,2,&
                        rcvbuf1ds,csize,mpi_short,pspace%lhs,2,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhs /= mpi_proc_null ) then
        l(0,1:nii) = rcvbuf1ds(1:nii)
      end if
      if ( pspace%lhs /= mpi_proc_null ) then
        sndbuf1ds(1:nii) = l(1,1:nii)
      end if
      call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%lhs,3,&
                        rcvbuf1ds,csize,mpi_short,pspace%rhs,3,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhs /= mpi_proc_null ) then
        l(njj+1,1:nii) = rcvbuf1ds(1:nii)
      end if
      if ( pspace%top /= mpi_proc_null ) then
        sndbuf1ds(1:njj) = l(1:njj,nii)
      end if
      call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%top,4,&
                        rcvbuf1ds,csize,mpi_short,pspace%btm,4,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%btm /= mpi_proc_null ) then
        l(1:njj,0) = rcvbuf1ds(1:njj)
      end if
      if ( pspace%rht /= mpi_proc_null ) then
        sndcrnps(1) = l(njj,nii)
      end if
      call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%rht,1, &
                        rcvcrnps,1,mpi_short,pspace%lhb,1, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        l(0,0) = rcvcrnps(1)
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        sndcrnps(1) = l(1,1)
      end if
      call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%lhb,2, &
                        rcvcrnps,1,mpi_short,pspace%rht,2, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        l(njj+1,nii+1) = rcvcrnps(1)
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        sndcrnps(1) = l(njj,1)
      end if
      call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%rhb,3, &
                        rcvcrnps,1,mpi_short,pspace%lht,3, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        l(0,nii+1) = rcvcrnps(1)
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        sndcrnps(1) = l(1,nii)
      end if
      call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%lht,4, &
                        rcvcrnps,1,mpi_short,pspace%rhb,4, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        l(njj+1,0) = rcvcrnps(1)
      end if
    end subroutine exchange_internal2d_s

    subroutine exchange_internal2d_l(l)
      implicit none
      logical , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      if ( pspace%btm /= mpi_proc_null ) then
        sndbuf1dl(1:njj) = l(1:njj,1)
      end if
      call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%btm,1,&
                        rcvbuf1dl,csize,mpi_logical,pspace%top,1,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%top /= mpi_proc_null ) then
        l(1:njj,nii+1) = rcvbuf1dl(1:njj)
      end if
      if ( pspace%rhs /= mpi_proc_null ) then
        sndbuf1dl(1:nii) = l(njj,1:nii)
      end if
      call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%rhs,2,&
                        rcvbuf1dl,csize,mpi_logical,pspace%lhs,2,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhs /= mpi_proc_null ) then
        l(0,1:nii) = rcvbuf1dl(1:nii)
      end if
      if ( pspace%lhs /= mpi_proc_null ) then
        sndbuf1dl(1:nii) = l(1,1:nii)
      end if
      call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%lhs,3,&
                        rcvbuf1dl,csize,mpi_logical,pspace%rhs,3,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhs /= mpi_proc_null ) then
        l(njj+1,1:nii) = rcvbuf1dl(1:nii)
      end if
      if ( pspace%top /= mpi_proc_null ) then
        sndbuf1dl(1:njj) = l(1:njj,nii)
      end if
      call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%top,4,&
                        rcvbuf1dl,csize,mpi_logical,pspace%btm,4,&
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%btm /= mpi_proc_null ) then
        l(1:njj,0) = rcvbuf1dl(1:njj)
      end if
      if ( pspace%rht /= mpi_proc_null ) then
        sndcrnpl(1) = l(njj,nii)
      end if
      call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%rht,1, &
                        rcvcrnpl,1,mpi_logical,pspace%lhb,1, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        l(0,0) = rcvcrnpl(1)
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        sndcrnpl(1) = l(1,1)
      end if
      call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%lhb,2, &
                        rcvcrnpl,1,mpi_logical,pspace%rht,2, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        l(njj+1,nii+1) = rcvcrnpl(1)
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        sndcrnpl(1) = l(njj,1)
      end if
      call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%rhb,3, &
                        rcvcrnpl,1,mpi_logical,pspace%lht,3, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        l(0,nii+1) = rcvcrnpl(1)
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        sndcrnpl(1) = l(1,nii)
      end if
      call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%lht,4, &
                        rcvcrnpl,1,mpi_logical,pspace%rhb,4, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        l(njj+1,0) = rcvcrnpl(1)
      end if
    end subroutine exchange_internal2d_l

    subroutine exchange_internal3d_d(l)
      implicit none
      real(dp) , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do k = nk1 , nk2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1dd(1:njj) = l(1:njj,1,k)
        end if
        call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%btm,1, &
                          rcvbuf1dd,csize,mpi_real8,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+1,k) = rcvbuf1dd(1:njj)
        end if
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1dd(1:nii) = l(njj,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%rhs,2, &
                          rcvbuf1dd,csize,mpi_real8,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(0,1:nii,k) = rcvbuf1dd(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1dd(1:nii) = l(1,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%lhs,3, &
                          rcvbuf1dd,csize,mpi_real8,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+1,1:nii,k) = rcvbuf1dd(1:nii)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1dd(1:njj) = l(1:njj,nii,k)
        end if
        call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%top,4, &
                          rcvbuf1dd,csize,mpi_real8,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,0,k) = rcvbuf1dd(1:njj)
        end if
        if ( pspace%rht /= mpi_proc_null ) then
          sndcrnpd(1) = l(njj,nii,k)
        end if
        call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%rht,1, &
                          rcvcrnpd,1,mpi_real8,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          l(0,0,k) = rcvcrnpd(1)
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          sndcrnpd(1) = l(1,1,k)
        end if
        call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%lhb,2, &
                          rcvcrnpd,1,mpi_real8,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          l(njj+1,nii+1,k) = rcvcrnpd(1)
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          sndcrnpd(1) = l(njj,1,k)
        end if
        call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%rhb,3, &
                          rcvcrnpd,1,mpi_real8,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          l(0,nii+1,k) = rcvcrnpd(1)
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          sndcrnpd(1) = l(1,nii,k)
        end if
        call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%lht,4, &
                          rcvcrnpd,1,mpi_real8,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          l(njj+1,0,k) = rcvcrnpd(1)
        end if
      end do
    end subroutine exchange_internal3d_d

    subroutine exchange_internal3d_r(l)
      implicit none
      real(sp) , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do k = nk1 , nk2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1dr(1:njj) = l(1:njj,1,k)
        end if
        call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%btm,1, &
                          rcvbuf1dr,csize,mpi_real4,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+1,k) = rcvbuf1dr(1:njj)
        end if
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1dr(1:nii) = l(njj,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%rhs,2, &
                          rcvbuf1dr,csize,mpi_real4,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(0,1:nii,k) = rcvbuf1dr(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1dr(1:nii) = l(1,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%lhs,3, &
                          rcvbuf1dr,csize,mpi_real4,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+1,1:nii,k) = rcvbuf1dr(1:nii)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1dr(1:njj) = l(1:njj,nii,k)
        end if
        call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%top,4, &
                          rcvbuf1dr,csize,mpi_real4,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,0,k) = rcvbuf1dr(1:njj)
        end if
        if ( pspace%rht /= mpi_proc_null ) then
          sndcrnpr(1) = l(njj,nii,k)
        end if
        call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%rht,1, &
                          rcvcrnpr,1,mpi_real4,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          l(0,0,k) = rcvcrnpr(1)
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          sndcrnpr(1) = l(1,1,k)
        end if
        call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%lhb,2, &
                          rcvcrnpr,1,mpi_real4,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          l(njj+1,nii+1,k) = rcvcrnpr(1)
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          sndcrnpr(1) = l(njj,1,k)
        end if
        call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%rhb,3, &
                          rcvcrnpr,1,mpi_real4,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          l(0,nii+1,k) = rcvcrnpr(1)
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          sndcrnpr(1) = l(1,nii,k)
        end if
        call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%lht,4, &
                          rcvcrnpr,1,mpi_real4,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          l(njj+1,0,k) = rcvcrnpr(1)
        end if
      end do
    end subroutine exchange_internal3d_r

    subroutine exchange_internal3d_i(l)
      implicit none
      integer , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do k = nk1 , nk2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1di(1:njj) = l(1:njj,1,k)
        end if
        call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%btm,1, &
                          rcvbuf1di,csize,mpi_integer,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+1,k) = rcvbuf1di(1:njj)
        end if
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1di(1:nii) = l(njj,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%rhs,2, &
                          rcvbuf1di,csize,mpi_integer,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(0,1:nii,k) = rcvbuf1di(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1di(1:nii) = l(1,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%lhs,3, &
                          rcvbuf1di,csize,mpi_integer,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+1,1:nii,k) = rcvbuf1di(1:nii)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1di(1:njj) = l(1:njj,nii,k)
        end if
        call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%top,4, &
                          rcvbuf1di,csize,mpi_integer,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,0,k) = rcvbuf1di(1:njj)
        end if
        if ( pspace%rht /= mpi_proc_null ) then
          sndcrnpi(1) = l(njj,nii,k)
        end if
        call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%rht,1, &
                          rcvcrnpi,1,mpi_integer,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          l(0,0,k) = rcvcrnpi(1)
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          sndcrnpi(1) = l(1,1,k)
        end if
        call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%lhb,2, &
                          rcvcrnpi,1,mpi_integer,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          l(njj+1,nii+1,k) = rcvcrnpi(1)
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          sndcrnpi(1) = l(njj,1,k)
        end if
        call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%rhb,3, &
                          rcvcrnpi,1,mpi_integer,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          l(0,nii+1,k) = rcvcrnpi(1)
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          sndcrnpi(1) = l(1,nii,k)
        end if
        call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%lht,4, &
                          rcvcrnpi,1,mpi_integer,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          l(njj+1,0,k) = rcvcrnpi(1)
        end if
      end do
    end subroutine exchange_internal3d_i

    subroutine exchange_internal3d_s(l)
      implicit none
      integer(2) , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do k = nk1 , nk2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1ds(1:njj) = l(1:njj,1,k)
        end if
        call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%btm,1, &
                          rcvbuf1ds,csize,mpi_short,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+1,k) = rcvbuf1ds(1:njj)
        end if
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1ds(1:nii) = l(njj,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%rhs,2, &
                          rcvbuf1ds,csize,mpi_short,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(0,1:nii,k) = rcvbuf1ds(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1ds(1:nii) = l(1,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%lhs,3, &
                          rcvbuf1ds,csize,mpi_short,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+1,1:nii,k) = rcvbuf1ds(1:nii)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1ds(1:njj) = l(1:njj,nii,k)
        end if
        call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%top,4, &
                          rcvbuf1ds,csize,mpi_short,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,0,k) = rcvbuf1ds(1:njj)
        end if
        if ( pspace%rht /= mpi_proc_null ) then
          sndcrnps(1) = l(njj,nii,k)
        end if
        call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%rht,1, &
                          rcvcrnps,1,mpi_short,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          l(0,0,k) = rcvcrnps(1)
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          sndcrnps(1) = l(1,1,k)
        end if
        call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%lhb,2, &
                          rcvcrnps,1,mpi_short,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          l(njj+1,nii+1,k) = rcvcrnps(1)
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          sndcrnps(1) = l(njj,1,k)
        end if
        call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%rhb,3, &
                          rcvcrnps,1,mpi_short,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          l(0,nii+1,k) = rcvcrnps(1)
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          sndcrnps(1) = l(1,nii,k)
        end if
        call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%lht,4, &
                          rcvcrnps,1,mpi_short,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          l(njj+1,0,k) = rcvcrnps(1)
        end if
      end do
    end subroutine exchange_internal3d_s

    subroutine exchange_internal3d_l(l)
      implicit none
      logical , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do k = nk1 , nk2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1dl(1:njj) = l(1:njj,1,k)
        end if
        call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%btm,1, &
                          rcvbuf1dl,csize,mpi_logical,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+1,k) = rcvbuf1dl(1:njj)
        end if
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1dl(1:nii) = l(njj,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%rhs,2, &
                          rcvbuf1dl,csize,mpi_logical,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(0,1:nii,k) = rcvbuf1dl(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1dl(1:nii) = l(1,1:nii,k)
        end if
        call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%lhs,3, &
                          rcvbuf1dl,csize,mpi_logical,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+1,1:nii,k) = rcvbuf1dl(1:nii)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1dl(1:njj) = l(1:njj,nii,k)
        end if
        call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%top,4, &
                          rcvbuf1dl,csize,mpi_logical,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,0,k) = rcvbuf1dl(1:njj)
        end if
        if ( pspace%rht /= mpi_proc_null ) then
          sndcrnpl(1) = l(njj,nii,k)
        end if
        call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%rht,1, &
                          rcvcrnpl,1,mpi_logical,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          l(0,0,k) = rcvcrnpl(1)
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          sndcrnpl(1) = l(1,1,k)
        end if
        call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%lhb,2, &
                          rcvcrnpl,1,mpi_logical,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          l(njj+1,nii+1,k) = rcvcrnpl(1)
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          sndcrnpl(1) = l(njj,1,k)
        end if
        call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%rhb,3, &
                          rcvcrnpl,1,mpi_logical,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          l(0,nii+1,k) = rcvcrnpl(1)
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          sndcrnpl(1) = l(1,nii,k)
        end if
        call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%lht,4, &
                          rcvcrnpl,1,mpi_logical,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          l(njj+1,0,k) = rcvcrnpl(1)
        end if
      end do
    end subroutine exchange_internal3d_l

    subroutine exchange_internal4d_d(l)
      implicit none
      real(dp) , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do t = nt1 , nt2
        do k = nk1 , nk2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1dd(1:njj) = l(1:njj,1,k,t)
          end if
          call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%btm,1, &
                          rcvbuf1dd,csize,mpi_real8,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+1,k,t) = rcvbuf1dd(1:njj)
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1dd(1:nii) = l(njj,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%rhs,2, &
                          rcvbuf1dd,csize,mpi_real8,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(0,1:nii,k,t) = rcvbuf1dd(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1dd(1:nii) = l(1,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%lhs,3, &
                          rcvbuf1dd,csize,mpi_real8,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+1,1:nii,k,t) = rcvbuf1dd(1:nii)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1dd(1:njj) = l(1:njj,nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%top,4, &
                          rcvbuf1dd,csize,mpi_real8,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,0,k,t) = rcvbuf1dd(1:njj)
          end if
          if ( pspace%rht /= mpi_proc_null ) then
            sndcrnpd(1) = l(njj,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%rht,1, &
                          rcvcrnpd,1,mpi_real8,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            l(0,0,k,t) = rcvcrnpd(1)
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            sndcrnpd(1) = l(1,1,k,t)
          end if
          call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%lhb,2, &
                          rcvcrnpd,1,mpi_real8,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            l(njj+1,nii+1,k,t) = rcvcrnpd(1)
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            sndcrnpd(1) = l(njj,1,k,t)
          end if
          call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%rhb,3, &
                          rcvcrnpd,1,mpi_real8,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            l(0,nii+1,k,t) = rcvcrnpd(1)
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            sndcrnpd(1) = l(1,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnpd,1,mpi_real8,pspace%lht,4, &
                          rcvcrnpd,1,mpi_real8,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            l(njj+1,0,k,t) = rcvcrnpd(1)
          end if
        end do
      end do
    end subroutine exchange_internal4d_d

    subroutine exchange_internal4d_r(l)
      implicit none
      real(sp) , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do t = nt1 , nt2
        do k = nk1 , nk2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1dr(1:njj) = l(1:njj,1,k,t)
          end if
          call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%btm,1, &
                          rcvbuf1dr,csize,mpi_real4,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+1,k,t) = rcvbuf1dr(1:njj)
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1dr(1:nii) = l(njj,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%rhs,2, &
                          rcvbuf1dr,csize,mpi_real4,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(0,1:nii,k,t) = rcvbuf1dr(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1dr(1:nii) = l(1,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%lhs,3, &
                          rcvbuf1dr,csize,mpi_real4,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+1,1:nii,k,t) = rcvbuf1dr(1:nii)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1dr(1:njj) = l(1:njj,nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%top,4, &
                          rcvbuf1dr,csize,mpi_real4,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,0,k,t) = rcvbuf1dr(1:njj)
          end if
          if ( pspace%rht /= mpi_proc_null ) then
            sndcrnpr(1) = l(njj,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%rht,1, &
                          rcvcrnpr,1,mpi_real4,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            l(0,0,k,t) = rcvcrnpr(1)
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            sndcrnpr(1) = l(1,1,k,t)
          end if
          call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%lhb,2, &
                          rcvcrnpr,1,mpi_real4,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            l(njj+1,nii+1,k,t) = rcvcrnpr(1)
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            sndcrnpr(1) = l(njj,1,k,t)
          end if
          call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%rhb,3, &
                          rcvcrnpr,1,mpi_real4,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            l(0,nii+1,k,t) = rcvcrnpr(1)
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            sndcrnpr(1) = l(1,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnpr,1,mpi_real4,pspace%lht,4, &
                          rcvcrnpr,1,mpi_real4,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            l(njj+1,0,k,t) = rcvcrnpr(1)
          end if
        end do
      end do
    end subroutine exchange_internal4d_r

    subroutine exchange_internal4d_i(l)
      implicit none
      integer , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do t = nt1 , nt2
        do k = nk1 , nk2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1di(1:njj) = l(1:njj,1,k,t)
          end if
          call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%btm,1, &
                          rcvbuf1di,csize,mpi_integer,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+1,k,t) = rcvbuf1di(1:njj)
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1di(1:nii) = l(njj,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%rhs,2, &
                          rcvbuf1di,csize,mpi_integer,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(0,1:nii,k,t) = rcvbuf1di(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1di(1:nii) = l(1,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%lhs,3, &
                          rcvbuf1di,csize,mpi_integer,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+1,1:nii,k,t) = rcvbuf1di(1:nii)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1di(1:njj) = l(1:njj,nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%top,4, &
                          rcvbuf1di,csize,mpi_integer,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,0,k,t) = rcvbuf1di(1:njj)
          end if
          if ( pspace%rht /= mpi_proc_null ) then
            sndcrnpi(1) = l(njj,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%rht,1, &
                          rcvcrnpi,1,mpi_integer,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            l(0,0,k,t) = rcvcrnpi(1)
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            sndcrnpi(1) = l(1,1,k,t)
          end if
          call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%lhb,2, &
                          rcvcrnpi,1,mpi_integer,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            l(njj+1,nii+1,k,t) = rcvcrnpi(1)
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            sndcrnpi(1) = l(njj,1,k,t)
          end if
          call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%rhb,3, &
                          rcvcrnpi,1,mpi_integer,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            l(0,nii+1,k,t) = rcvcrnpi(1)
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            sndcrnpi(1) = l(1,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnpi,1,mpi_integer,pspace%lht,4, &
                          rcvcrnpi,1,mpi_integer,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            l(njj+1,0,k,t) = rcvcrnpi(1)
          end if
        end do
      end do
    end subroutine exchange_internal4d_i

    subroutine exchange_internal4d_s(l)
      implicit none
      integer(2) , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do t = nt1 , nt2
        do k = nk1 , nk2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1ds(1:njj) = l(1:njj,1,k,t)
          end if
          call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%btm,1, &
                          rcvbuf1ds,csize,mpi_short,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+1,k,t) = rcvbuf1ds(1:njj)
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1ds(1:nii) = l(njj,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%rhs,2, &
                          rcvbuf1ds,csize,mpi_short,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(0,1:nii,k,t) = rcvbuf1ds(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1ds(1:nii) = l(1,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%lhs,3, &
                          rcvbuf1ds,csize,mpi_short,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+1,1:nii,k,t) = rcvbuf1ds(1:nii)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1ds(1:njj) = l(1:njj,nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1ds,csize,mpi_short,pspace%top,4, &
                          rcvbuf1ds,csize,mpi_short,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,0,k,t) = rcvbuf1ds(1:njj)
          end if
          if ( pspace%rht /= mpi_proc_null ) then
            sndcrnps(1) = l(njj,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%rht,1, &
                          rcvcrnps,1,mpi_short,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            l(0,0,k,t) = rcvcrnps(1)
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            sndcrnps(1) = l(1,1,k,t)
          end if
          call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%lhb,2, &
                          rcvcrnps,1,mpi_short,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            l(njj+1,nii+1,k,t) = rcvcrnps(1)
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            sndcrnps(1) = l(njj,1,k,t)
          end if
          call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%rhb,3, &
                          rcvcrnps,1,mpi_short,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            l(0,nii+1,k,t) = rcvcrnps(1)
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            sndcrnps(1) = l(1,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnps,1,mpi_short,pspace%lht,4, &
                          rcvcrnps,1,mpi_short,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            l(njj+1,0,k,t) = rcvcrnps(1)
          end if
        end do
      end do
    end subroutine exchange_internal4d_s

    subroutine exchange_internal4d_l(l)
      implicit none
      logical , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      if ( lbound(l,1) /= 0 .and. lbound(l,2) /= 0 ) then
        return
      end if
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      nii = ubound(l,2) - 1
      njj = ubound(l,1) - 1
      do t = nt1 , nt2
        do k = nk1 , nk2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1dl(1:njj) = l(1:njj,1,k,t)
          end if
          call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%btm,1, &
                          rcvbuf1dl,csize,mpi_logical,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+1,k,t) = rcvbuf1dl(1:njj)
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1dl(1:nii) = l(njj,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%rhs,2, &
                          rcvbuf1dl,csize,mpi_logical,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(0,1:nii,k,t) = rcvbuf1dl(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1dl(1:nii) = l(1,1:nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%lhs,3, &
                          rcvbuf1dl,csize,mpi_logical,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+1,1:nii,k,t) = rcvbuf1dl(1:nii)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1dl(1:njj) = l(1:njj,nii,k,t)
          end if
          call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%top,4, &
                          rcvbuf1dl,csize,mpi_logical,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,0,k,t) = rcvbuf1dl(1:njj)
          end if
          if ( pspace%rht /= mpi_proc_null ) then
            sndcrnpl(1) = l(njj,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%rht,1, &
                          rcvcrnpl,1,mpi_logical,pspace%lhb,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            l(0,0,k,t) = rcvcrnpl(1)
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            sndcrnpl(1) = l(1,1,k,t)
          end if
          call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%lhb,2, &
                          rcvcrnpl,1,mpi_logical,pspace%rht,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            l(njj+1,nii+1,k,t) = rcvcrnpl(1)
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            sndcrnpl(1) = l(njj,1,k,t)
          end if
          call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%rhb,3, &
                          rcvcrnpl,1,mpi_logical,pspace%lht,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            l(0,nii+1,k,t) = rcvcrnpl(1)
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            sndcrnpl(1) = l(1,nii,k,t)
          end if
          call mpi_sendrecv(sndcrnpl,1,mpi_logical,pspace%lht,4, &
                          rcvcrnpl,1,mpi_logical,pspace%rhb,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            l(njj+1,0,k,t) = rcvcrnpl(1)
          end if
        end do
      end do
    end subroutine exchange_internal4d_l

    subroutine global_to_proc2d_d(g,l)
      implicit none
      real(dp) , pointer , dimension(:,:) , intent(in) :: g
      real(dp) , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_d

    subroutine proc_to_global2d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:) , intent(in) :: l
      real(dp) , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_d

    subroutine global_to_proc2d_r(g,l)
      implicit none
      real(sp) , pointer , dimension(:,:) , intent(in) :: g
      real(sp) , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_r

    subroutine proc_to_global2d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:) , intent(in) :: l
      real(sp) , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_r

    subroutine global_to_proc2d_i(g,l)
      implicit none
      integer , pointer , dimension(:,:) , intent(in) :: g
      integer , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_i

    subroutine proc_to_global2d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:) , intent(in) :: l
      integer , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_i

    subroutine global_to_proc2d_s(g,l)
      implicit none
      integer(2) , pointer , dimension(:,:) , intent(in) :: g
      integer(2) , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_s

    subroutine proc_to_global2d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:) , intent(in) :: l
      integer(2) , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_s

    subroutine global_to_proc2d_l(g,l)
      implicit none
      logical , pointer , dimension(:,:) , intent(in) :: g
      logical , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_l

    subroutine proc_to_global2d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:) , intent(in) :: l
      logical , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_l

    subroutine global_to_proc3d_d(g,l)
      implicit none
      real(dp) , pointer , dimension(:,:,:) , intent(in) :: g
      real(dp) , pointer , dimension(:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        l(1:njj,1:nii,k) = g(lbgjj:ubgjj,lbgii:ubgii,k)
      end do
    end subroutine global_to_proc3d_d

    subroutine proc_to_global3d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:,:) , intent(in) :: l
      real(dp) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        g(lbgjj:ubgjj,lbgii:ubgii,k) = l(1:njj,1:nii,k)
      end do
    end subroutine proc_to_global3d_d

    subroutine global_to_proc3d_r(g,l)
      implicit none
      real(sp) , pointer , dimension(:,:,:) , intent(in) :: g
      real(sp) , pointer , dimension(:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        l(1:njj,1:nii,k) = g(lbgjj:ubgjj,lbgii:ubgii,k)
      end do
    end subroutine global_to_proc3d_r

    subroutine proc_to_global3d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:,:) , intent(in) :: l
      real(sp) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        g(lbgjj:ubgjj,lbgii:ubgii,k) = l(1:njj,1:nii,k)
      end do
    end subroutine proc_to_global3d_r

    subroutine global_to_proc3d_i(g,l)
      implicit none
      integer , pointer , dimension(:,:,:) , intent(in) :: g
      integer , pointer , dimension(:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        l(1:njj,1:nii,k) = g(lbgjj:ubgjj,lbgii:ubgii,k)
      end do
    end subroutine global_to_proc3d_i

    subroutine proc_to_global3d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:,:) , intent(in) :: l
      integer , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        g(lbgjj:ubgjj,lbgii:ubgii,k) = l(1:njj,1:nii,k)
      end do
    end subroutine proc_to_global3d_i

    subroutine global_to_proc3d_s(g,l)
      implicit none
      integer(2) , pointer , dimension(:,:,:) , intent(in) :: g
      integer(2) , pointer , dimension(:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        l(1:njj,1:nii,k) = g(lbgjj:ubgjj,lbgii:ubgii,k)
      end do
    end subroutine global_to_proc3d_s

    subroutine proc_to_global3d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:,:) , intent(in) :: l
      integer(2) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        g(lbgjj:ubgjj,lbgii:ubgii,k) = l(1:njj,1:nii,k)
      end do
    end subroutine proc_to_global3d_s

    subroutine global_to_proc3d_l(g,l)
      implicit none
      logical , pointer , dimension(:,:,:) , intent(in) :: g
      logical , pointer , dimension(:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        l(1:njj,1:nii,k) = g(lbgjj:ubgjj,lbgii:ubgii,k)
      end do
    end subroutine global_to_proc3d_l

    subroutine proc_to_global3d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:,:) , intent(in) :: l
      logical , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do k = k1 , k2
        g(lbgjj:ubgjj,lbgii:ubgii,k) = l(1:njj,1:nii,k)
      end do
    end subroutine proc_to_global3d_l

    subroutine global_to_proc4d_d(g,l)
      implicit none
      real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: g
      real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          l(1:njj,1:nii,k,t) = g(lbgjj:ubgjj,lbgii:ubgii,k,t)
        end do
      end do
    end subroutine global_to_proc4d_d

    subroutine proc_to_global4d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: l
      real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          g(lbgjj:ubgjj,lbgii:ubgii,k,t) = l(1:njj,1:nii,k,t)
        end do
      end do
    end subroutine proc_to_global4d_d

    subroutine global_to_proc4d_r(g,l)
      implicit none
      real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: g
      real(sp) , pointer , dimension(:,:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          l(1:njj,1:nii,k,t) = g(lbgjj:ubgjj,lbgii:ubgii,k,t)
        end do
      end do
    end subroutine global_to_proc4d_r

    subroutine proc_to_global4d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: l
      real(sp) , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          g(lbgjj:ubgjj,lbgii:ubgii,k,t) = l(1:njj,1:nii,k,t)
        end do
      end do
    end subroutine proc_to_global4d_r

    subroutine global_to_proc4d_i(g,l)
      implicit none
      integer , pointer , dimension(:,:,:,:) , intent(in) :: g
      integer , pointer , dimension(:,:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          l(1:njj,1:nii,k,t) = g(lbgjj:ubgjj,lbgii:ubgii,k,t)
        end do
      end do
    end subroutine global_to_proc4d_i

    subroutine proc_to_global4d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:,:,:) , intent(in) :: l
      integer , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          g(lbgjj:ubgjj,lbgii:ubgii,k,t) = l(1:njj,1:nii,k,t)
        end do
      end do
    end subroutine proc_to_global4d_i

    subroutine global_to_proc4d_s(g,l)
      implicit none
      integer(2) , pointer , dimension(:,:,:,:) , intent(in) :: g
      integer(2) , pointer , dimension(:,:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          l(1:njj,1:nii,k,t) = g(lbgjj:ubgjj,lbgii:ubgii,k,t)
        end do
      end do
    end subroutine global_to_proc4d_s

    subroutine proc_to_global4d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:,:,:) , intent(in) :: l
      integer(2) , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          g(lbgjj:ubgjj,lbgii:ubgii,k,t) = l(1:njj,1:nii,k,t)
        end do
      end do
    end subroutine proc_to_global4d_s

    subroutine global_to_proc4d_l(g,l)
      implicit none
      logical , pointer , dimension(:,:,:,:) , intent(in) :: g
      logical , pointer , dimension(:,:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          l(1:njj,1:nii,k,t) = g(lbgjj:ubgjj,lbgii:ubgii,k,t)
        end do
      end do
    end subroutine global_to_proc4d_l

    subroutine proc_to_global4d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:,:,:) , intent(in) :: l
      logical , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer :: k , k1 , k2 , t , t1 , t2 , nii , njj
      integer :: lbgii , lbgjj , ubgii , ubgjj
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( lbound(l,1) == 0 ) then
        nii = ubound(l,2)-1
        njj = ubound(l,1)-1
      else
        nii = ubound(l,2)
        njj = ubound(l,1)
      end if
      lbgii = pspace%g_i1
      lbgjj = pspace%g_j1
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      do t = t1 , t2
        do k = k1 , k2
          g(lbgjj:ubgjj,lbgii:ubgii,k,t) = l(1:njj,1:nii,k,t)
        end do
      end do
    end subroutine proc_to_global4d_l

    subroutine master_to_nodes_d(v)
      implicit none
      real(dp) , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,1,mpi_real8,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes_d

    subroutine master_to_nodes_r(v)
      implicit none
      real(sp) , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,1,mpi_real4,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes_r

    subroutine master_to_nodes_i(v)
      implicit none
      integer , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,1,mpi_integer,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes_i

    subroutine master_to_nodes_s(v)
      implicit none
      integer(2) , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,1,mpi_short,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes_s

    subroutine master_to_nodes_l(v)
      implicit none
      logical , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,1,mpi_logical,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes_l

    subroutine master_to_nodes1d_d(v)
      implicit none
      real(dp) , dimension(:) , pointer , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,size(v),mpi_real8,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes1d_d

    subroutine master_to_nodes1d_r(v)
      implicit none
      real(sp) , dimension(:) , pointer , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,size(v),mpi_real4,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes1d_r

    subroutine master_to_nodes1d_i(v)
      implicit none
      integer , dimension(:) , pointer , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,size(v),mpi_integer,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes1d_i

    subroutine master_to_nodes1d_s(v)
      implicit none
      integer(2) , dimension(:) , pointer , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,size(v),mpi_short,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes1d_s

    subroutine master_to_nodes1d_l(v)
      implicit none
      logical , dimension(:) , pointer , intent(inout) :: v
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      call mpi_bcast(v,size(v),mpi_logical,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
    end subroutine master_to_nodes1d_l

    subroutine master_to_nodes2d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:) , intent(inout) :: l
      real(dp) , pointer , dimension(:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                mnode%excbuf2dd(icount) = g(j,i)
                icount = icount + 1
              end do
            end do
            call mpi_send(mnode%excbuf2dd,mnode%pgsize(inode),mpi_real8, &
                          inode,0,pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          else
            call global_to_proc2d_d(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        call mpi_recv(excbuf2dd,pspace%totalpoints,mpi_real8,masterproc,0, &
                      pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            l(j,i) = excbuf2dd(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine master_to_nodes2d_d

    subroutine master_to_nodes2d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:) , intent(inout) :: l
      real(sp) , pointer , dimension(:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                mnode%excbuf2dr(icount) = g(j,i)
                icount = icount + 1
              end do
            end do
            call mpi_send(mnode%excbuf2dr,mnode%pgsize(inode),mpi_real4, &
                          inode,0,pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          else
            call global_to_proc2d_r(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        call mpi_recv(excbuf2dr,pspace%totalpoints,mpi_real4,masterproc,0, &
                      pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            l(j,i) = excbuf2dr(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine master_to_nodes2d_r

    subroutine master_to_nodes2d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:) , intent(inout) :: l
      integer , pointer , dimension(:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                mnode%excbuf2di(icount) = g(j,i)
                icount = icount + 1
              end do
            end do
            call mpi_send(mnode%excbuf2di,mnode%pgsize(inode),mpi_integer, &
                          inode,0,pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          else
            call global_to_proc2d_i(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        call mpi_recv(excbuf2di,pspace%totalpoints,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            l(j,i) = excbuf2di(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine master_to_nodes2d_i

    subroutine master_to_nodes2d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:) , intent(inout) :: l
      integer(2) , pointer , dimension(:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                mnode%excbuf2ds(icount) = g(j,i)
                icount = icount + 1
              end do
            end do
            call mpi_send(mnode%excbuf2ds,mnode%pgsize(inode),mpi_short, &
                          inode,0,pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          else
            call global_to_proc2d_s(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        call mpi_recv(excbuf2ds,pspace%totalpoints,mpi_short,masterproc,0, &
                      pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            l(j,i) = excbuf2ds(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine master_to_nodes2d_s

    subroutine master_to_nodes2d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:) , intent(inout) :: l
      logical , pointer , dimension(:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                mnode%excbuf2dl(icount) = g(j,i)
                icount = icount + 1
              end do
            end do
            call mpi_send(mnode%excbuf2dl,mnode%pgsize(inode),mpi_logical, &
                          inode,0,pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          else
            call global_to_proc2d_l(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        call mpi_recv(excbuf2dl,pspace%totalpoints,mpi_logical,masterproc,0, &
                      pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            l(j,i) = excbuf2dl(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine master_to_nodes2d_l

    subroutine master_to_nodes3d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:,:) , intent(inout) :: l
      real(dp) , pointer , dimension(:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  mnode%excbuf2dd(icount) = g(j,i,k)
                  icount = icount + 1
                end do
              end do
              call mpi_send(mnode%excbuf2dd,mnode%pgsize(inode),mpi_real8, &
                            inode,0,pspace%cartesian_communicator,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            end do
          else
            call global_to_proc3d_d(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          call mpi_recv(excbuf2dd,pspace%totalpoints,mpi_real8,masterproc,0, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              l(j,i,k) = excbuf2dd(icount)
              icount = icount + 1
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes3d_d

    subroutine master_to_nodes3d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:,:) , intent(inout) :: l
      real(sp) , pointer , dimension(:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  mnode%excbuf2dr(icount) = g(j,i,k)
                  icount = icount + 1
                end do
              end do
              call mpi_send(mnode%excbuf2dr,mnode%pgsize(inode),mpi_real4, &
                            inode,0,pspace%cartesian_communicator,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            end do
          else
            call global_to_proc3d_r(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          call mpi_recv(excbuf2dr,pspace%totalpoints,mpi_real4,masterproc,0, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              l(j,i,k) = excbuf2dr(icount)
              icount = icount + 1
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes3d_r

    subroutine master_to_nodes3d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:,:) , intent(inout) :: l
      integer , pointer , dimension(:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  mnode%excbuf2di(icount) = g(j,i,k)
                  icount = icount + 1
                end do
              end do
              call mpi_send(mnode%excbuf2di,mnode%pgsize(inode),mpi_integer, &
                            inode,0,pspace%cartesian_communicator,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            end do
          else
            call global_to_proc3d_i(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          call mpi_recv(excbuf2di,pspace%totalpoints,mpi_integer,masterproc,0, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              l(j,i,k) = excbuf2di(icount)
              icount = icount + 1
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes3d_i

    subroutine master_to_nodes3d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:,:) , intent(inout) :: l
      integer(2) , pointer , dimension(:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  mnode%excbuf2ds(icount) = g(j,i,k)
                  icount = icount + 1
                end do
              end do
              call mpi_send(mnode%excbuf2ds,mnode%pgsize(inode),mpi_short, &
                            inode,0,pspace%cartesian_communicator,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            end do
          else
            call global_to_proc3d_s(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          call mpi_recv(excbuf2ds,pspace%totalpoints,mpi_short,masterproc,0, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              l(j,i,k) = excbuf2ds(icount)
              icount = icount + 1
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes3d_s

    subroutine master_to_nodes3d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:,:) , intent(inout) :: l
      logical , pointer , dimension(:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  mnode%excbuf2dl(icount) = g(j,i,k)
                  icount = icount + 1
                end do
              end do
              call mpi_send(mnode%excbuf2dl,mnode%pgsize(inode),mpi_logical, &
                            inode,0,pspace%cartesian_communicator,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            end do
          else
            call global_to_proc3d_l(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          call mpi_recv(excbuf2dl,pspace%totalpoints,mpi_logical,masterproc,0, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              l(j,i,k) = excbuf2dl(icount)
              icount = icount + 1
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes3d_l

    subroutine master_to_nodes4d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: l
      real(dp) , pointer , dimension(:,:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    mnode%excbuf2dd(icount) = g(j,i,k,t)
                    icount = icount + 1
                  end do
                end do
                call mpi_send(mnode%excbuf2dd,mnode%pgsize(inode),mpi_real8, &
                              inode,0,pspace%cartesian_communicator,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
              end do
            end do
          else
            call global_to_proc4d_d(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            call mpi_recv(excbuf2dd,pspace%totalpoints,mpi_real8,masterproc,0, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                l(j,i,k,t) = excbuf2dd(icount)
                icount = icount + 1
              end do
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes4d_d

    subroutine master_to_nodes4d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:,:,:) , intent(inout) :: l
      real(sp) , pointer , dimension(:,:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    mnode%excbuf2dr(icount) = g(j,i,k,t)
                    icount = icount + 1
                  end do
                end do
                call mpi_send(mnode%excbuf2dr,mnode%pgsize(inode),mpi_real4, &
                              inode,0,pspace%cartesian_communicator,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
              end do
            end do
          else
            call global_to_proc4d_r(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            call mpi_recv(excbuf2dr,pspace%totalpoints,mpi_real4,masterproc,0, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                l(j,i,k,t) = excbuf2dr(icount)
                icount = icount + 1
              end do
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes4d_r

    subroutine master_to_nodes4d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:,:,:) , intent(inout) :: l
      integer , pointer , dimension(:,:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    mnode%excbuf2di(icount) = g(j,i,k,t)
                    icount = icount + 1
                  end do
                end do
                call mpi_send(mnode%excbuf2di,mnode%pgsize(inode),mpi_integer, &
                              inode,0,pspace%cartesian_communicator,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
              end do
            end do
          else
            call global_to_proc4d_i(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            call mpi_recv(excbuf2di,pspace%totalpoints,mpi_integer, &
                          masterproc,0,pspace%cartesian_communicator, &
                          mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                l(j,i,k,t) = excbuf2di(icount)
                icount = icount + 1
              end do
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes4d_i

    subroutine master_to_nodes4d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:,:,:) , intent(inout) :: l
      integer(2) , pointer , dimension(:,:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    mnode%excbuf2ds(icount) = g(j,i,k,t)
                    icount = icount + 1
                  end do
                end do
                call mpi_send(mnode%excbuf2ds,mnode%pgsize(inode),mpi_short, &
                              inode,0,pspace%cartesian_communicator,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
              end do
            end do
          else
            call global_to_proc4d_s(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            call mpi_recv(excbuf2ds,pspace%totalpoints,mpi_short, &
                          masterproc,0,pspace%cartesian_communicator, &
                          mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                l(j,i,k,t) = excbuf2ds(icount)
                icount = icount + 1
              end do
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes4d_s

    subroutine master_to_nodes4d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:,:,:) , intent(inout) :: l
      logical , pointer , dimension(:,:,:,:) , intent(in) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    mnode%excbuf2dl(icount) = g(j,i,k,t)
                    icount = icount + 1
                  end do
                end do
                call mpi_send(mnode%excbuf2dl,mnode%pgsize(inode),mpi_logical, &
                              inode,0,pspace%cartesian_communicator,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
              end do
            end do
          else
            call global_to_proc4d_l(g,l)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            call mpi_recv(excbuf2dl,pspace%totalpoints,mpi_logical, &
                          masterproc,0,pspace%cartesian_communicator, &
                          mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                l(j,i,k,t) = excbuf2dl(icount)
                icount = icount + 1
              end do
            end do
          end do
        end do
      end if
    end subroutine master_to_nodes4d_l

    subroutine nodes_to_master2d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:) , intent(in) :: l
      real(dp) , pointer , dimension(:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using nodes_to_master as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in nodes_to_master from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            call mpi_recv(mnode%excbuf2dd,mnode%pgsize(inode),mpi_real8, &
                          inode,0,pspace%cartesian_communicator, &
                          mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                g(j,i) = mnode%excbuf2dd(icount)
                icount = icount + 1
              end do
            end do
          else
            call proc_to_global2d_d(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using nodes_to_master as non master')
        end if
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            excbuf2dd(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
        call mpi_send(excbuf2dd,pspace%totalpoints,mpi_real8,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine nodes_to_master2d_d

    subroutine nodes_to_master2d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:) , intent(in) :: l
      real(sp) , pointer , dimension(:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using nodes_to_master as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in nodes_to_master from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            call mpi_recv(mnode%excbuf2dr,mnode%pgsize(inode),mpi_real4, &
                          inode,0,pspace%cartesian_communicator, &
                          mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                g(j,i) = mnode%excbuf2dr(icount)
                icount = icount + 1
              end do
            end do
          else
            call proc_to_global2d_r(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using nodes_to_master as non master')
        end if
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            excbuf2dr(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
        call mpi_send(excbuf2dr,pspace%totalpoints,mpi_real4,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine nodes_to_master2d_r

    subroutine nodes_to_master2d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:) , intent(in) :: l
      integer , pointer , dimension(:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using nodes_to_master as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in nodes_to_master from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            call mpi_recv(mnode%excbuf2di,mnode%pgsize(inode),mpi_integer, &
                          inode,0,pspace%cartesian_communicator, &
                          mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                g(j,i) = mnode%excbuf2di(icount)
                icount = icount + 1
              end do
            end do
          else
            call proc_to_global2d_i(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using nodes_to_master as non master')
        end if
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            excbuf2di(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
        call mpi_send(excbuf2di,pspace%totalpoints,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine nodes_to_master2d_i

    subroutine nodes_to_master2d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:) , intent(in) :: l
      integer(2) , pointer , dimension(:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using nodes_to_master as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in nodes_to_master from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            call mpi_recv(mnode%excbuf2ds,mnode%pgsize(inode),mpi_short, &
                          inode,0,pspace%cartesian_communicator, &
                          mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                g(j,i) = mnode%excbuf2ds(icount)
                icount = icount + 1
              end do
            end do
          else
            call proc_to_global2d_s(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using nodes_to_master as non master')
        end if
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            excbuf2ds(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
        call mpi_send(excbuf2ds,pspace%totalpoints,mpi_short,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine nodes_to_master2d_s

    subroutine nodes_to_master2d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:) , intent(in) :: l
      logical , pointer , dimension(:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using nodes_to_master as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in nodes_to_master from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            call mpi_recv(mnode%excbuf2dl,mnode%pgsize(inode),mpi_logical, &
                          inode,0,pspace%cartesian_communicator, &
                          mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            icount = 1
            do i = mnode%pgis(inode) , mnode%pgie(inode)
              do j = mnode%pgjs(inode) , mnode%pgje(inode)
                g(j,i) = mnode%excbuf2dl(icount)
                icount = icount + 1
              end do
            end do
          else
            call proc_to_global2d_l(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using nodes_to_master as non master')
        end if
        icount = 1
        do i = 1 , pspace%p_i
          do j = 1 , pspace%p_j
            excbuf2dl(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
        call mpi_send(excbuf2dl,pspace%totalpoints,mpi_logical,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine nodes_to_master2d_l

    subroutine nodes_to_master3d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:,:) , intent(in) :: l
      real(dp) , pointer , dimension(:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              call mpi_recv(mnode%excbuf2dd,mnode%pgsize(inode),mpi_real8, &
                            inode,0,pspace%cartesian_communicator, &
                            mpi_status_ignore,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  g(j,i,k) = mnode%excbuf2dd(icount)
                  icount = icount + 1
                end do
              end do
            end do
          else
            call proc_to_global3d_d(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              excbuf2dd(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
          call mpi_send(excbuf2dd,pspace%totalpoints,mpi_real8,masterproc,0, &
                        pspace%cartesian_communicator,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        end do
      end if
    end subroutine nodes_to_master3d_d

    subroutine nodes_to_master3d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:,:) , intent(in) :: l
      real(sp) , pointer , dimension(:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              call mpi_recv(mnode%excbuf2dr,mnode%pgsize(inode),mpi_real4, &
                            inode,0,pspace%cartesian_communicator, &
                            mpi_status_ignore,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  g(j,i,k) = mnode%excbuf2dr(icount)
                  icount = icount + 1
                end do
              end do
            end do
          else
            call proc_to_global3d_r(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              excbuf2dr(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
          call mpi_send(excbuf2dr,pspace%totalpoints,mpi_real4,masterproc,0, &
                        pspace%cartesian_communicator,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        end do
      end if
    end subroutine nodes_to_master3d_r

    subroutine nodes_to_master3d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:,:) , intent(in) :: l
      integer , pointer , dimension(:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              call mpi_recv(mnode%excbuf2di,mnode%pgsize(inode),mpi_integer, &
                            inode,0,pspace%cartesian_communicator, &
                            mpi_status_ignore,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  g(j,i,k) = mnode%excbuf2di(icount)
                  icount = icount + 1
                end do
              end do
            end do
          else
            call proc_to_global3d_i(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              excbuf2di(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
          call mpi_send(excbuf2di,pspace%totalpoints,mpi_integer,masterproc,0, &
                        pspace%cartesian_communicator,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        end do
      end if
    end subroutine nodes_to_master3d_i

    subroutine nodes_to_master3d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:,:) , intent(in) :: l
      integer(2) , pointer , dimension(:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              call mpi_recv(mnode%excbuf2ds,mnode%pgsize(inode),mpi_short, &
                            inode,0,pspace%cartesian_communicator, &
                            mpi_status_ignore,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  g(j,i,k) = mnode%excbuf2ds(icount)
                  icount = icount + 1
                end do
              end do
            end do
          else
            call proc_to_global3d_s(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              excbuf2ds(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
          call mpi_send(excbuf2ds,pspace%totalpoints,mpi_short,masterproc,0, &
                        pspace%cartesian_communicator,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        end do
      end if
    end subroutine nodes_to_master3d_s

    subroutine nodes_to_master3d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:,:) , intent(in) :: l
      logical , pointer , dimension(:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do k = gk1 , gk2
              call mpi_recv(mnode%excbuf2dl,mnode%pgsize(inode),mpi_logical, &
                            inode,0,pspace%cartesian_communicator, &
                            mpi_status_ignore,ierr)
              if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
              icount = 1
              do i = mnode%pgis(inode) , mnode%pgie(inode)
                do j = mnode%pgjs(inode) , mnode%pgje(inode)
                  g(j,i,k) = mnode%excbuf2dl(icount)
                  icount = icount + 1
                end do
              end do
            end do
          else
            call proc_to_global3d_l(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do k = k1 , k2
          icount = 1
          do i = 1 , pspace%p_i
            do j = 1 , pspace%p_j
              excbuf2dl(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
          call mpi_send(excbuf2dl,pspace%totalpoints,mpi_logical,masterproc,0, &
                        pspace%cartesian_communicator,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        end do
      end if
    end subroutine nodes_to_master3d_l

    subroutine nodes_to_master4d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: l
      real(dp) , pointer , dimension(:,:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                call mpi_recv(mnode%excbuf2dd,mnode%pgsize(inode),mpi_real8, &
                              inode,0,pspace%cartesian_communicator, &
                              mpi_status_ignore,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    g(j,i,k,t) = mnode%excbuf2dd(icount)
                    icount = icount + 1
                  end do
                end do
              end do
            end do
          else
            call proc_to_global4d_d(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                excbuf2dd(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
            call mpi_send(excbuf2dd,pspace%totalpoints,mpi_real8,masterproc,0, &
                          pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          end do
        end do
      end if
    end subroutine nodes_to_master4d_d

    subroutine nodes_to_master4d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: l
      real(sp) , pointer , dimension(:,:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                call mpi_recv(mnode%excbuf2dr,mnode%pgsize(inode),mpi_real4, &
                              inode,0,pspace%cartesian_communicator, &
                              mpi_status_ignore,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    g(j,i,k,t) = mnode%excbuf2dr(icount)
                    icount = icount + 1
                  end do
                end do
              end do
            end do
          else
            call proc_to_global4d_r(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                excbuf2dr(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
            call mpi_send(excbuf2dr,pspace%totalpoints,mpi_real4,masterproc,0, &
                          pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          end do
        end do
      end if
    end subroutine nodes_to_master4d_r

    subroutine nodes_to_master4d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:,:,:) , intent(in) :: l
      integer , pointer , dimension(:,:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                call mpi_recv(mnode%excbuf2di,mnode%pgsize(inode),mpi_integer, &
                              inode,0,pspace%cartesian_communicator, &
                              mpi_status_ignore,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    g(j,i,k,t) = mnode%excbuf2di(icount)
                    icount = icount + 1
                  end do
                end do
              end do
            end do
          else
            call proc_to_global4d_i(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                excbuf2di(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
            call mpi_send(excbuf2di,pspace%totalpoints,mpi_integer, &
                          masterproc,0,pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          end do
        end do
      end if
    end subroutine nodes_to_master4d_i

    subroutine nodes_to_master4d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:,:,:) , intent(in) :: l
      integer(2) , pointer , dimension(:,:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                call mpi_recv(mnode%excbuf2ds,mnode%pgsize(inode),mpi_short, &
                              inode,0,pspace%cartesian_communicator, &
                              mpi_status_ignore,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    g(j,i,k,t) = mnode%excbuf2ds(icount)
                    icount = icount + 1
                  end do
                end do
              end do
            end do
          else
            call proc_to_global4d_s(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                excbuf2ds(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
            call mpi_send(excbuf2ds,pspace%totalpoints,mpi_short, &
                          masterproc,0,pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          end do
        end do
      end if
    end subroutine nodes_to_master4d_s

    subroutine nodes_to_master4d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:,:,:) , intent(in) :: l
      logical , pointer , dimension(:,:,:,:) , intent(inout) , optional :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: t , t1 , t2 , gt1 , gt2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( present(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using master_to_nodes as master')
        end if
        gt1 = lbound(g,4)
        gt2 = ubound(g,4)
        gk1 = lbound(g,3)
        gk2 = ubound(g,3)
        if ( gk1 /= k1 .or. gk2 /= k2 .or. gt1 /= t1 .or. gt2 /= t2 ) then
          call fatal(__FILE__,__LINE__, &
          'Not equal vertical bounds in master_to_nodes')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in master_to_nodes from master')
        end if
        do inode = 0 , xproc%total_cpus-1
          if ( inode /= masterproc ) then
            do t = gt1 , gt2
              do k = gk1 , gk2
                call mpi_recv(mnode%excbuf2dl,mnode%pgsize(inode),mpi_logical, &
                              inode,0,pspace%cartesian_communicator, &
                              mpi_status_ignore,ierr)
                if ( ierr /= mpi_success ) then
                  call mpi_fatal(__FILE__,__LINE__,ierr)
                end if
                icount = 1
                do i = mnode%pgis(inode) , mnode%pgie(inode)
                  do j = mnode%pgjs(inode) , mnode%pgje(inode)
                    g(j,i,k,t) = mnode%excbuf2dl(icount)
                    icount = icount + 1
                  end do
                end do
              end do
            end do
          else
            call proc_to_global4d_l(l,g)
          end if
        end do
      else
        if ( am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Master node using master_to_nodes as non master')
        end if
        do t = t1 , t2
          do k = k1 , k2
            icount = 1
            do i = 1 , pspace%p_i
              do j = 1 , pspace%p_j
                excbuf2dl(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
            call mpi_send(excbuf2dl,pspace%totalpoints,mpi_logical, &
                          masterproc,0,pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          end do
        end do
      end if
    end subroutine nodes_to_master4d_l

    subroutine global_sum_d(s)
      implicit none
      real(dp) , intent(inout) :: s
      real(dp) :: retval
      integer :: ierr
      call mpi_allreduce(s,retval,1,mpi_real8,mpi_sum, &
                         pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      s = retval
    end subroutine global_sum_d

    subroutine global_sum_r(s)
      implicit none
      real(sp) , intent(inout) :: s
      real(sp) :: retval
      integer :: ierr
      call mpi_allreduce(s,retval,1,mpi_real4,mpi_sum, &
                         pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      s = retval
    end subroutine global_sum_r

    subroutine global_sum_i(s)
      implicit none
      integer , intent(inout) :: s
      integer :: retval
      integer :: ierr
      call mpi_allreduce(s,retval,1,mpi_integer,mpi_sum, &
                         pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      s = retval
    end subroutine global_sum_i

    subroutine mpi_fatal(f,l,iecode)
      implicit none
      integer , intent(in) :: iecode
      integer , intent(in) :: l
      character(len=*) , intent(in) :: f
      integer :: emlen , ierr
      character(len=256) :: mpiemsg
      emlen = 256
      if ( cantalk( ) ) then
        call mpi_error_string(iecode,mpiemsg,emlen,ierr)
        write (stderr, *) '----------------------------------------------'
        write (stderr, *) '             FATAL CONDITION MET'
        write (stderr, *) '----------------------------------------------'
        write (stderr, *) 'File : ',f
        write (stderr, *) 'Line : ',l
        write (stderr, *) 'Cond : ',mpiemsg(1:emlen)
        write (stderr, *) '----------------------------------------------'
      end if
      call mpi_abort(gspace%global_communicator,233,ierr)
    end subroutine mpi_fatal

end module mod_mppgrid

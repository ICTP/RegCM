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
  use mod_mpmessage
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

  logical , parameter , public :: global_grid    = .true.
  logical , parameter , public :: processor_grid = .false.

  integer , parameter , public :: ghost_none     = 0
  integer , parameter , public :: ghost_1_point  = 1
  integer , parameter , public :: ghost_2_points = 2

  logical , parameter , public :: dot_grid   = .true.
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
    real(dp) , pointer , dimension(:) :: excbuf2dd
    real(sp) , pointer , dimension(:) :: excbuf2dr
    integer , pointer , dimension(:) :: excbuf2di
    integer(2) , pointer , dimension(:) :: excbuf2ds
    logical , pointer , dimension(:) :: excbuf2dl
  end type masternode

  ! The part this processor is going to compute the tendencies on
  type processor_domain
    integer :: p_i = mindimsize ! Local to processor i points
    integer :: p_j = mindimsize ! Local to processor j points
    integer :: g_is = 1
    integer :: g_js = 1
    integer :: g_ie = 1
    integer :: g_je = 1
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

  type global_boundary4d_d
    real(dp) , pointer , dimension(:,:,:,:) :: north
    real(dp) , pointer , dimension(:,:,:,:) :: south
    real(dp) , pointer , dimension(:,:,:,:) :: east
    real(dp) , pointer , dimension(:,:,:,:) :: west
  end type global_boundary4d_d

  type global_boundary3d_d
    real(dp) , pointer , dimension(:,:,:) :: north
    real(dp) , pointer , dimension(:,:,:) :: south
    real(dp) , pointer , dimension(:,:,:) :: east
    real(dp) , pointer , dimension(:,:,:) :: west
  end type global_boundary3d_d

  type global_boundary2d_d
    real(dp) , pointer , dimension(:,:) :: north
    real(dp) , pointer , dimension(:,:) :: south
    real(dp) , pointer , dimension(:,:) :: east
    real(dp) , pointer , dimension(:,:) :: west
  end type global_boundary2d_d

  type global_boundary4d_r
    real(sp) , pointer , dimension(:,:,:,:) :: north
    real(sp) , pointer , dimension(:,:,:,:) :: south
    real(sp) , pointer , dimension(:,:,:,:) :: east
    real(sp) , pointer , dimension(:,:,:,:) :: west
  end type global_boundary4d_r

  type global_boundary3d_r
    real(sp) , pointer , dimension(:,:,:) :: north
    real(sp) , pointer , dimension(:,:,:) :: south
    real(sp) , pointer , dimension(:,:,:) :: east
    real(sp) , pointer , dimension(:,:,:) :: west
  end type global_boundary3d_r

  type global_boundary2d_r
    real(sp) , pointer , dimension(:,:) :: north
    real(sp) , pointer , dimension(:,:) :: south
    real(sp) , pointer , dimension(:,:) :: east
    real(sp) , pointer , dimension(:,:) :: west
  end type global_boundary2d_r

  type global_boundary4d_i
    integer , pointer , dimension(:,:,:,:) :: north
    integer , pointer , dimension(:,:,:,:) :: south
    integer , pointer , dimension(:,:,:,:) :: east
    integer , pointer , dimension(:,:,:,:) :: west
  end type global_boundary4d_i

  type global_boundary3d_i
    integer , pointer , dimension(:,:,:) :: north
    integer , pointer , dimension(:,:,:) :: south
    integer , pointer , dimension(:,:,:) :: east
    integer , pointer , dimension(:,:,:) :: west
  end type global_boundary3d_i

  type global_boundary2d_i
    integer , pointer , dimension(:,:) :: north
    integer , pointer , dimension(:,:) :: south
    integer , pointer , dimension(:,:) :: east
    integer , pointer , dimension(:,:) :: west
  end type global_boundary2d_i

  type global_boundary4d_s
    integer(2) , pointer , dimension(:,:,:,:) :: north
    integer(2) , pointer , dimension(:,:,:,:) :: south
    integer(2) , pointer , dimension(:,:,:,:) :: east
    integer(2) , pointer , dimension(:,:,:,:) :: west
  end type global_boundary4d_s

  type global_boundary3d_s
    integer(2) , pointer , dimension(:,:,:) :: north
    integer(2) , pointer , dimension(:,:,:) :: south
    integer(2) , pointer , dimension(:,:,:) :: east
    integer(2) , pointer , dimension(:,:,:) :: west
  end type global_boundary3d_s

  type global_boundary2d_s
    integer(2) , pointer , dimension(:,:) :: north
    integer(2) , pointer , dimension(:,:) :: south
    integer(2) , pointer , dimension(:,:) :: east
    integer(2) , pointer , dimension(:,:) :: west
  end type global_boundary2d_s

  type global_boundary4d_l
    logical , pointer , dimension(:,:,:,:) :: north
    logical , pointer , dimension(:,:,:,:) :: south
    logical , pointer , dimension(:,:,:,:) :: east
    logical , pointer , dimension(:,:,:,:) :: west
  end type global_boundary4d_l

  type global_boundary3d_l
    logical , pointer , dimension(:,:,:) :: north
    logical , pointer , dimension(:,:,:) :: south
    logical , pointer , dimension(:,:,:) :: east
    logical , pointer , dimension(:,:,:) :: west
  end type global_boundary3d_l

  type global_boundary2d_l
    logical , pointer , dimension(:,:) :: north
    logical , pointer , dimension(:,:) :: south
    logical , pointer , dimension(:,:) :: east
    logical , pointer , dimension(:,:) :: west
  end type global_boundary2d_l

  integer :: csize = 0
  integer :: gsize = 0
  integer :: bdysize = 0
  integer :: btmbdy1 = 0
  integer :: btmbdy2 = 0
  integer :: topbdy1 = 0
  integer :: topbdy2 = 0
  integer :: lhsbdy1 = 0
  integer :: lhsbdy2 = 0
  integer :: rhsbdy1 = 0
  integer :: rhsbdy2 = 0
  integer :: gbtmbdy1 = 0
  integer :: gbtmbdy2 = 0
  integer :: gtopbdy1 = 0
  integer :: gtopbdy2 = 0
  integer :: glhsbdy1 = 0
  integer :: glhsbdy2 = 0
  integer :: grhsbdy1 = 0
  integer :: grhsbdy2 = 0
  logical :: jperiodic = .false.
  real(dp) , pointer , dimension(:) :: sndbuf1dd
  real(dp) , pointer , dimension(:) :: rcvbuf1dd
  real(dp) , pointer , dimension(:) :: excbuf2dd
  real(sp) , pointer , dimension(:) :: sndbuf1dr
  real(sp) , pointer , dimension(:) :: rcvbuf1dr
  real(sp) , pointer , dimension(:) :: excbuf2dr
  integer , pointer , dimension(:) :: sndbuf1di
  integer , pointer , dimension(:) :: rcvbuf1di
  integer , pointer , dimension(:) :: excbuf2di
  integer(2) , pointer , dimension(:) :: sndbuf1ds
  integer(2) , pointer , dimension(:) :: rcvbuf1ds
  integer(2) , pointer , dimension(:) :: excbuf2ds
  logical , pointer , dimension(:) :: sndbuf1dl
  logical , pointer , dimension(:) :: rcvbuf1dl
  logical , pointer , dimension(:) :: excbuf2dl

  type (model_processors) :: xproc
  type (model_domain) :: gspace
  type (processor_domain) :: pspace
  type (masternode) :: mnode

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

  interface getbdy
    module procedure getbdy2d_d
    module procedure getbdy2d_r
    module procedure getbdy2d_i
    module procedure getbdy2d_s
    module procedure getbdy2d_l
    module procedure getbdy3d_d
    module procedure getbdy3d_r
    module procedure getbdy3d_i
    module procedure getbdy3d_s
    module procedure getbdy3d_l
    module procedure getbdy4d_d
    module procedure getbdy4d_r
    module procedure getbdy4d_i
    module procedure getbdy4d_s
    module procedure getbdy4d_l
  end interface getbdy

  interface global_to_globbdy
    module procedure global_to_globbdy2d_d
    module procedure global_to_globbdy2d_r
    module procedure global_to_globbdy2d_i
    module procedure global_to_globbdy2d_s
    module procedure global_to_globbdy2d_l
    module procedure global_to_globbdy3d_d
    module procedure global_to_globbdy3d_r
    module procedure global_to_globbdy3d_i
    module procedure global_to_globbdy3d_s
    module procedure global_to_globbdy3d_l
    module procedure global_to_globbdy4d_d
    module procedure global_to_globbdy4d_r
    module procedure global_to_globbdy4d_i
    module procedure global_to_globbdy4d_s
    module procedure global_to_globbdy4d_l
  end interface global_to_globbdy

  interface set_external
    module procedure set_external2d_d
    module procedure set_external2d_r
    module procedure set_external2d_i
    module procedure set_external2d_s
    module procedure set_external2d_l
    module procedure set_external3d_d
    module procedure set_external3d_r
    module procedure set_external3d_i
    module procedure set_external3d_s
    module procedure set_external3d_l
    module procedure set_external4d_d
    module procedure set_external4d_r
    module procedure set_external4d_i
    module procedure set_external4d_s
    module procedure set_external4d_l
  end interface set_external

  type (procbounds) , protected :: pbnds

  public :: am_i_master , cantalk , toggle_mpi_debug
  public :: setup_domain , delete_domain
  public :: getgrid , getbdy
  public :: exchange_internal
  public :: global_to_proc , proc_to_global , global_to_globbdy
  public :: master_to_nodes , nodes_to_master
  public :: set_external
  public :: global_sum
  public :: global_boundary4d_d , global_boundary3d_d , global_boundary2d_d
  public :: global_boundary4d_r , global_boundary3d_r , global_boundary2d_r
  public :: global_boundary4d_i , global_boundary3d_i , global_boundary2d_i
  public :: global_boundary4d_s , global_boundary3d_s , global_boundary2d_s
  public :: global_boundary4d_l , global_boundary3d_l , global_boundary2d_l

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
        write(mpilogunit,'(a,i0)') 'Global Start I  ', pspace%g_is
        write(mpilogunit,'(a,i0)') 'Global Start J  ', pspace%g_js
        write(mpilogunit,'(a,i0)') 'Mine number I   ', pspace%p_i
        write(mpilogunit,'(a,i0)') 'Mine number J   ', pspace%p_j
        write(mpilogunit,'(a,i0)') 'Totalpoints     ', pspace%totalpoints
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

    subroutine setup_domain(ni,nj,nbdy,iband,comm)
      implicit none
      integer , intent(in) :: ni , nj , nbdy , iband
      integer , intent(in) , optional :: comm
      integer :: ierr , jbi , jbj , imaxcpus , max_pi , max_pj , max_p
      integer :: inode , maxec , maxgbl , maxcpu
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
          write(stderr,*) 'Suggested minimum number of CPUS : ', imaxcpus
        end if
        call fatal(__FILE__,__LINE__,'Scheme not working')
      end if
      if ( iband == 1 ) then
        xproc%dim_period(1) = .true.
        jperiodic = .true.
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
      pspace%g_js = pspace%location(1) * jbj + 1
      pspace%g_is = pspace%location(2) * jbi + 1
      pspace%g_je = pspace%g_js + jbj - 1
      pspace%g_ie = pspace%g_is + jbi - 1
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
          else
            mnode%pgis(inode) = pspace%g_is
            mnode%pgjs(inode) = pspace%g_js
            mnode%pgie(inode) = pspace%g_ie
            mnode%pgje(inode) = pspace%g_je
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
        maxgbl = max(ni,nj)
        call getmem1d(mnode%excbuf2dd,1,maxec,__FILE__)
        call getmem1d(mnode%excbuf2dr,1,maxec,__FILE__)
        call getmem1d(mnode%excbuf2di,1,maxec,__FILE__)
        call getmem1d(mnode%excbuf2ds,1,maxec,__FILE__)
        call getmem1d(mnode%excbuf2dl,1,maxec,__FILE__)
      else
        call mpi_send(pspace%g_is,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        call mpi_send(pspace%g_js,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        call mpi_send(pspace%g_ie,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        call mpi_send(pspace%g_je,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        call mpi_send(pspace%totalpoints,1,mpi_integer,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        if ( cantalk( ) .and. .false. ) then
          write(stdout, *) '###################################'
          write(stdout,'(a,i0)') 'CPU number: ', gspace%global_rank
          write(stdout,'(a,i0)') 'Start   i : ', pspace%g_is
          write(stdout,'(a,i0)') 'End     i : ', pspace%g_ie
          write(stdout,'(a,i0)') 'Start   j : ', pspace%g_js
          write(stdout,'(a,i0)') 'End     j : ', pspace%g_je
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
      bdysize = nbdy
      if ( pspace%btm == mpi_proc_null ) then
        btmbdy1 = 1
        btmbdy2 = bdysize
      end if
      if ( pspace%top == mpi_proc_null ) then
        topbdy1 = pspace%p_i-bdysize+1
        topbdy2 = pspace%p_i
      end if
      if ( .not. jperiodic ) then
        if ( pspace%rhs == mpi_proc_null ) then
          rhsbdy1 = pspace%p_j-bdysize+1
          rhsbdy2 = pspace%p_j
        end if
        if ( pspace%lhs == mpi_proc_null ) then
          lhsbdy1 = 1
          lhsbdy2 = bdysize
        end if
      end if
      gbtmbdy1 = 1
      gbtmbdy2 = bdysize
      gtopbdy1 = gspace%g_i-bdysize+1
      gtopbdy2 = gspace%g_i
      if ( .not. jperiodic ) then
        grhsbdy1 = gspace%g_j-bdysize+1
        grhsbdy2 = gspace%g_j
        glhsbdy1 = 1
        glhsbdy2 = bdysize
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

    subroutine getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      implicit none
      logical , intent(in) :: isglobal , isstagger
      integer , intent(in) :: ghostp
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
        i1 = 1 - ghostp
        i2 = pspace%p_i + ghostp
        j1 = 1 - ghostp
        j2 = pspace%p_j + ghostp
      end if
      if ( isstagger ) then
        i2 = i2 + 1
        j2 = j2 + 1
      end if
    end subroutine getextrema

    subroutine getgrid2d_d(g,lglobal,lstagger,nghost)
      implicit none
      real(dp) , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_d

    subroutine getgrid2d_r(g,lglobal,lstagger,nghost)
      implicit none
      real(sp) , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_r

    subroutine getgrid2d_i(g,lglobal,lstagger,nghost)
      implicit none
      integer , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_i

    subroutine getgrid2d_s(g,lglobal,lstagger,nghost)
      implicit none
      integer(2) , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_s

    subroutine getgrid2d_l(g,lglobal,lstagger,nghost)
      implicit none
      logical , pointer , dimension(:,:) , intent(inout) :: g
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem2d(g,j1,j2,i1,i2,__FILE__)
    end subroutine getgrid2d_l

    subroutine getgrid3d_d(g,k1,k2,lglobal,lstagger,nghost)
      implicit none
      real(dp) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_d

    subroutine getgrid3d_r(g,k1,k2,lglobal,lstagger,nghost)
      implicit none
      real(sp) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_r

    subroutine getgrid3d_i(g,k1,k2,lglobal,lstagger,nghost)
      implicit none
      integer , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_i

    subroutine getgrid3d_s(g,k1,k2,lglobal,lstagger,nghost)
      implicit none
      integer(2) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_s

    subroutine getgrid3d_l(g,k1,k2,lglobal,lstagger,nghost)
      implicit none
      logical , pointer , dimension(:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem3d(g,j1,j2,i1,i2,k1,k2,__FILE__)
    end subroutine getgrid3d_l

    subroutine getgrid4d_d(g,k1,k2,t1,t2,lglobal,lstagger,nghost)
      implicit none
      real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_d

    subroutine getgrid4d_r(g,k1,k2,t1,t2,lglobal,lstagger,nghost)
      implicit none
      real(sp) , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_r

    subroutine getgrid4d_i(g,k1,k2,t1,t2,lglobal,lstagger,nghost)
      implicit none
      integer , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_i

    subroutine getgrid4d_s(g,k1,k2,t1,t2,lglobal,lstagger,nghost)
      implicit none
      integer(2) , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_s

    subroutine getgrid4d_l(g,k1,k2,t1,t2,lglobal,lstagger,nghost)
      implicit none
      logical , pointer , dimension(:,:,:,:) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lglobal
      logical , intent(in) , optional :: lstagger
      integer , intent(in) , optional :: nghost
      logical :: isglobal = processor_grid
      logical :: isstagger = cross_grid
      integer :: ghostp = ghost_none
      integer :: i1 , i2 , j1 , j2
      if ( present(lglobal) )   isglobal   = lglobal
      if ( present(lstagger) )  isstagger  = lstagger
      if ( present(nghost) )    ghostp = nghost
      call getextrema(isglobal,isstagger,ghostp,i1,i2,j1,j2)
      call getmem4d(g,j1,j2,i1,i2,k1,k2,t1,t2,__FILE__)
    end subroutine getgrid4d_l

    subroutine exchange_internal2d_d(l)
      implicit none
      real(dp) , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj , icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do i = 1 , esize2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1dd(1:njj) = l(1:njj,i)
        end if
        call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%btm,1, &
                          rcvbuf1dd,csize,mpi_real8,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+i) = rcvbuf1dd(1:njj)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1dd(1:njj) = l(1:njj,nii-i+1)
        end if
        call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%top,4, &
                          rcvbuf1dd,csize,mpi_real8,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,-i+1) = rcvbuf1dd(1:njj)
        end if
      end do
      do j = 1 , esize1
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1dd(1:nii) = l(njj-j+1,1:nii)
        end if
        call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%rhs,2, &
                          rcvbuf1dd,csize,mpi_real8,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(-j+1,1:nii) = rcvbuf1dd(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1dd(1:nii) = l(j,1:nii)
        end if
        call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%lhs,3, &
                          rcvbuf1dd,csize,mpi_real8,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+j,1:nii) = rcvbuf1dd(1:nii)
        end if
      end do
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dd(icount) = l(njj-j+1,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%rht,5, &
                        rcvbuf1dd,bsize,mpi_real8,pspace%lhb,5, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,-i+1) = rcvbuf1dd(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dd(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%lhb,6, &
                        rcvbuf1dd,bsize,mpi_real8,pspace%rht,6, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,nii+i) = rcvbuf1dd(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dd(icount) = l(njj-j+1,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%rhb,7, &
                        rcvbuf1dd,bsize,mpi_real8,pspace%lht,7, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,nii+i) = rcvbuf1dd(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dd(icount) = l(j,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%lht,8, &
                        rcvbuf1dd,bsize,mpi_real8,pspace%rhb,8, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,-i+1) = rcvbuf1dd(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine exchange_internal2d_d

    subroutine exchange_internal2d_r(l)
      implicit none
      real(sp) , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj , icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do i = 1 , esize2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1dr(1:njj) = l(1:njj,i)
        end if
        call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%btm,1, &
                          rcvbuf1dr,csize,mpi_real4,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+i) = rcvbuf1dr(1:njj)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1dr(1:njj) = l(1:njj,nii-i+1)
        end if
        call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%top,4, &
                          rcvbuf1dr,csize,mpi_real4,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,-i+1) = rcvbuf1dr(1:njj)
        end if
      end do
      do j = 1 , esize1
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1dr(1:nii) = l(njj-j+1,1:nii)
        end if
        call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%rhs,2, &
                          rcvbuf1dr,csize,mpi_real4,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(-j+1,1:nii) = rcvbuf1dr(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1dr(1:nii) = l(j,1:nii)
        end if
        call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%lhs,3, &
                          rcvbuf1dr,csize,mpi_real4,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+j,1:nii) = rcvbuf1dr(1:nii)
        end if
      end do
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dr(icount) = l(njj-j+1,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%rht,5, &
                        rcvbuf1dr,bsize,mpi_real4,pspace%lhb,5, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,-i+1) = rcvbuf1dr(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dr(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%lhb,6, &
                        rcvbuf1dr,bsize,mpi_real4,pspace%rht,6, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,nii+i) = rcvbuf1dr(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dr(icount) = l(njj-j+1,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%rhb,7, &
                        rcvbuf1dr,bsize,mpi_real4,pspace%lht,7, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,nii+i) = rcvbuf1dr(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dr(icount) = l(j,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%lht,8, &
                        rcvbuf1dr,bsize,mpi_real4,pspace%rhb,8, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,-i+1) = rcvbuf1dr(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine exchange_internal2d_r

    subroutine exchange_internal2d_i(l)
      implicit none
      integer , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj , icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do i = 1 , esize2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1di(1:njj) = l(1:njj,i)
        end if
        call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%btm,1, &
                          rcvbuf1di,csize,mpi_integer,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+i) = rcvbuf1di(1:njj)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1di(1:njj) = l(1:njj,nii-i+1)
        end if
        call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%top,4, &
                          rcvbuf1di,csize,mpi_integer,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,-i+1) = rcvbuf1di(1:njj)
        end if
      end do
      do j = 1 , esize1
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1di(1:nii) = l(njj-j+1,1:nii)
        end if
        call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%rhs,2, &
                          rcvbuf1di,csize,mpi_integer,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(-j+1,1:nii) = rcvbuf1di(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1di(1:nii) = l(j,1:nii)
        end if
        call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%lhs,3, &
                          rcvbuf1di,csize,mpi_integer,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+j,1:nii) = rcvbuf1di(1:nii)
        end if
      end do
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1di(icount) = l(njj-j+1,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%rht,5, &
                        rcvbuf1di,bsize,mpi_integer,pspace%lhb,5, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,-i+1) = rcvbuf1di(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1di(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%lhb,6, &
                        rcvbuf1di,bsize,mpi_integer,pspace%rht,6, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,nii+i) = rcvbuf1di(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1di(icount) = l(njj-j+1,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%rhb,7, &
                        rcvbuf1di,bsize,mpi_integer,pspace%lht,7, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,nii+i) = rcvbuf1di(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1di(icount) = l(j,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%lht,8, &
                        rcvbuf1di,bsize,mpi_integer,pspace%rhb,8, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,-i+1) = rcvbuf1di(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine exchange_internal2d_i

    subroutine exchange_internal2d_s(l)
      implicit none
      integer(2) , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj , icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do i = 1 , esize2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1ds(1:njj) = l(1:njj,i)
        end if
        call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%btm,1, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+i) = rcvbuf1ds(1:njj)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1ds(1:njj) = l(1:njj,nii-i+1)
        end if
        call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%top,4, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,-i+1) = rcvbuf1ds(1:njj)
        end if
      end do
      do j = 1 , esize1
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1ds(1:nii) = l(njj-j+1,1:nii)
        end if
        call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%rhs,2, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(-j+1,1:nii) = rcvbuf1ds(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1ds(1:nii) = l(j,1:nii)
        end if
        call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%lhs,3, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+j,1:nii) = rcvbuf1ds(1:nii)
        end if
      end do
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1ds(icount) = l(njj-j+1,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%rht,5, &
                        rcvbuf1ds,bsize,mpi_integer2,pspace%lhb,5, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,-i+1) = rcvbuf1ds(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1ds(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%lhb,6, &
                        rcvbuf1ds,bsize,mpi_integer2,pspace%rht,6, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,nii+i) = rcvbuf1ds(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1ds(icount) = l(njj-j+1,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%rhb,7, &
                        rcvbuf1ds,bsize,mpi_integer2,pspace%lht,7, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,nii+i) = rcvbuf1ds(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1ds(icount) = l(j,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%lht,8, &
                        rcvbuf1ds,bsize,mpi_integer2,pspace%rhb,8, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,-i+1) = rcvbuf1ds(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine exchange_internal2d_s

    subroutine exchange_internal2d_l(l)
      implicit none
      logical , dimension(:,:) , pointer , intent(inout) :: l
      integer :: ierr
      integer :: nii , njj , icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do i = 1 , esize2
        if ( pspace%btm /= mpi_proc_null ) then
          sndbuf1dl(1:njj) = l(1:njj,i)
        end if
        call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%btm,1, &
                          rcvbuf1dl,csize,mpi_logical,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%top /= mpi_proc_null ) then
          l(1:njj,nii+i) = rcvbuf1dl(1:njj)
        end if
        if ( pspace%top /= mpi_proc_null ) then
          sndbuf1dl(1:njj) = l(1:njj,nii-i+1)
        end if
        call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%top,4, &
                          rcvbuf1dl,csize,mpi_logical,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%btm /= mpi_proc_null ) then
          l(1:njj,-i+1) = rcvbuf1dl(1:njj)
        end if
      end do
      do j = 1 , esize1
        if ( pspace%rhs /= mpi_proc_null ) then
          sndbuf1dl(1:nii) = l(njj-j+1,1:nii)
        end if
        call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%rhs,2, &
                          rcvbuf1dl,csize,mpi_logical,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhs /= mpi_proc_null ) then
          l(-j+1,1:nii) = rcvbuf1dl(1:nii)
        end if
        if ( pspace%lhs /= mpi_proc_null ) then
          sndbuf1dl(1:nii) = l(j,1:nii)
        end if
        call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%lhs,3, &
                          rcvbuf1dl,csize,mpi_logical,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhs /= mpi_proc_null ) then
          l(njj+j,1:nii) = rcvbuf1dl(1:nii)
        end if
      end do
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dl(icount) = l(njj-j+1,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%rht,5, &
                        rcvbuf1dl,bsize,mpi_logical,pspace%lhb,5, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,-i+1) = rcvbuf1dl(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dl(icount) = l(j,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%lhb,6, &
                        rcvbuf1dl,bsize,mpi_logical,pspace%rht,6, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,nii+i) = rcvbuf1dl(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dl(icount) = l(njj-j+1,i)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%rhb,7, &
                        rcvbuf1dl,bsize,mpi_logical,pspace%lht,7, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(-j+1,nii+i) = rcvbuf1dl(icount)
            icount = icount + 1
          end do
        end do
      end if
      if ( pspace%lht /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            sndbuf1dl(icount) = l(j,nii-i+1)
            icount = icount + 1
          end do
        end do
      end if
      call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%lht,8, &
                        rcvbuf1dl,bsize,mpi_logical,pspace%rhb,8, &
                        pspace%cartesian_communicator,mpi_status_ignore,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( pspace%rhb /= mpi_proc_null ) then
        icount = 0
        do i = 1 , esize2
          do j = 1 , esize1
            l(njj+j,-i+1) = rcvbuf1dl(icount)
            icount = icount + 1
          end do
        end do
      end if
    end subroutine exchange_internal2d_l

    subroutine exchange_internal3d_d(l)
      implicit none
      real(dp) , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do k = nk1 , nk2
        do i = 1 , esize2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1dd(1:njj) = l(1:njj,i,k)
          end if
          call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%btm,1, &
                          rcvbuf1dd,csize,mpi_real8,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+i,k) = rcvbuf1dd(1:njj)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1dd(1:njj) = l(1:njj,nii-i+1,k)
          end if
          call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%top,4, &
                          rcvbuf1dd,csize,mpi_real8,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,-i+1,k) = rcvbuf1dd(1:njj)
          end if
        end do
        do j = 1 , esize1
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1dd(1:nii) = l(njj-j+1,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%rhs,2, &
                          rcvbuf1dd,csize,mpi_real8,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(-j+1,1:nii,k) = rcvbuf1dd(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1dd(1:nii) = l(j,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%lhs,3, &
                          rcvbuf1dd,csize,mpi_real8,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+j,1:nii,k) = rcvbuf1dd(1:nii)
          end if
        end do
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dd(icount) = l(njj-j+1,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%rht,5, &
                          rcvbuf1dd,bsize,mpi_real8,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,-i+1,k) = rcvbuf1dd(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dd(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%lhb,6, &
                          rcvbuf1dd,bsize,mpi_real8,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,nii+i,k) = rcvbuf1dd(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dd(icount) = l(njj-j+1,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%rhb,7, &
                          rcvbuf1dd,bsize,mpi_real8,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,nii+i,k) = rcvbuf1dd(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dd(icount) = l(j,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%lht,8, &
                          rcvbuf1dd,bsize,mpi_real8,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,-i+1,k) = rcvbuf1dd(icount)
              icount = icount + 1
            end do
          end do
        end if
      end do
    end subroutine exchange_internal3d_d

    subroutine exchange_internal3d_r(l)
      implicit none
      real(sp) , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do k = nk1 , nk2
        do i = 1 , esize2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1dr(1:njj) = l(1:njj,i,k)
          end if
          call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%btm,1, &
                          rcvbuf1dr,csize,mpi_real4,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+i,k) = rcvbuf1dr(1:njj)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1dr(1:njj) = l(1:njj,nii-i+1,k)
          end if
          call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%top,4, &
                          rcvbuf1dr,csize,mpi_real4,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,-i+1,k) = rcvbuf1dr(1:njj)
          end if
        end do
        do j = 1 , esize1
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1dr(1:nii) = l(njj-j+1,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%rhs,2, &
                          rcvbuf1dr,csize,mpi_real4,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(-j+1,1:nii,k) = rcvbuf1dr(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1dr(1:nii) = l(j,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%lhs,3, &
                          rcvbuf1dr,csize,mpi_real4,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+j,1:nii,k) = rcvbuf1dr(1:nii)
          end if
        end do
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dr(icount) = l(njj-j+1,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%rht,5, &
                          rcvbuf1dr,bsize,mpi_real4,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,-i+1,k) = rcvbuf1dr(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dr(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%lhb,6, &
                          rcvbuf1dr,bsize,mpi_real4,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,nii+i,k) = rcvbuf1dr(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dr(icount) = l(njj-j+1,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%rhb,7, &
                          rcvbuf1dr,bsize,mpi_real4,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,nii+i,k) = rcvbuf1dr(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dr(icount) = l(j,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%lht,8, &
                          rcvbuf1dr,bsize,mpi_real4,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,-i+1,k) = rcvbuf1dr(icount)
              icount = icount + 1
            end do
          end do
        end if
      end do
    end subroutine exchange_internal3d_r

    subroutine exchange_internal3d_i(l)
      implicit none
      integer , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do k = nk1 , nk2
        do i = 1 , esize2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1di(1:njj) = l(1:njj,i,k)
          end if
          call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%btm,1, &
                          rcvbuf1di,csize,mpi_integer,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+i,k) = rcvbuf1di(1:njj)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1di(1:njj) = l(1:njj,nii-i+1,k)
          end if
          call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%top,4, &
                          rcvbuf1di,csize,mpi_integer,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,-i+1,k) = rcvbuf1di(1:njj)
          end if
        end do
        do j = 1 , esize1
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1di(1:nii) = l(njj-j+1,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%rhs,2, &
                          rcvbuf1di,csize,mpi_integer,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(-j+1,1:nii,k) = rcvbuf1di(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1di(1:nii) = l(j,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%lhs,3, &
                          rcvbuf1di,csize,mpi_integer,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+j,1:nii,k) = rcvbuf1di(1:nii)
          end if
        end do
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1di(icount) = l(njj-j+1,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%rht,5, &
                          rcvbuf1di,bsize,mpi_integer,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,-i+1,k) = rcvbuf1di(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1di(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%lhb,6, &
                          rcvbuf1di,bsize,mpi_integer,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,nii+i,k) = rcvbuf1di(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1di(icount) = l(njj-j+1,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%rhb,7, &
                          rcvbuf1di,bsize,mpi_integer,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,nii+i,k) = rcvbuf1di(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1di(icount) = l(j,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%lht,8, &
                          rcvbuf1di,bsize,mpi_integer,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,-i+1,k) = rcvbuf1di(icount)
              icount = icount + 1
            end do
          end do
        end if
      end do
    end subroutine exchange_internal3d_i

    subroutine exchange_internal3d_s(l)
      implicit none
      integer(2) , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do k = nk1 , nk2
        do i = 1 , esize2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1ds(1:njj) = l(1:njj,i,k)
          end if
          call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%btm,1, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+i,k) = rcvbuf1ds(1:njj)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1ds(1:njj) = l(1:njj,nii-i+1,k)
          end if
          call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%top,4, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,-i+1,k) = rcvbuf1ds(1:njj)
          end if
        end do
        do j = 1 , esize1
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1ds(1:nii) = l(njj-j+1,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%rhs,2, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(-j+1,1:nii,k) = rcvbuf1ds(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1ds(1:nii) = l(j,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%lhs,3, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+j,1:nii,k) = rcvbuf1ds(1:nii)
          end if
        end do
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1ds(icount) = l(njj-j+1,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%rht,5, &
                          rcvbuf1ds,bsize,mpi_integer2,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,-i+1,k) = rcvbuf1ds(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1ds(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%lhb,6, &
                          rcvbuf1ds,bsize,mpi_integer2,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,nii+i,k) = rcvbuf1ds(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1ds(icount) = l(njj-j+1,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%rhb,7, &
                          rcvbuf1ds,bsize,mpi_integer2,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,nii+i,k) = rcvbuf1ds(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1ds(icount) = l(j,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%lht,8, &
                          rcvbuf1ds,bsize,mpi_integer2,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,-i+1,k) = rcvbuf1ds(icount)
              icount = icount + 1
            end do
          end do
        end if
      end do
    end subroutine exchange_internal3d_s

    subroutine exchange_internal3d_l(l)
      implicit none
      logical , dimension(:,:,:) , pointer , intent(inout) :: l
      integer :: k , nii , njj , nk1 , nk2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do k = nk1 , nk2
        do i = 1 , esize2
          if ( pspace%btm /= mpi_proc_null ) then
            sndbuf1dl(1:njj) = l(1:njj,i,k)
          end if
          call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%btm,1, &
                          rcvbuf1dl,csize,mpi_logical,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%top /= mpi_proc_null ) then
            l(1:njj,nii+i,k) = rcvbuf1dl(1:njj)
          end if
          if ( pspace%top /= mpi_proc_null ) then
            sndbuf1dl(1:njj) = l(1:njj,nii-i+1,k)
          end if
          call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%top,4, &
                          rcvbuf1dl,csize,mpi_logical,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%btm /= mpi_proc_null ) then
            l(1:njj,-i+1,k) = rcvbuf1dl(1:njj)
          end if
        end do
        do j = 1 , esize1
          if ( pspace%rhs /= mpi_proc_null ) then
            sndbuf1dl(1:nii) = l(njj-j+1,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%rhs,2, &
                          rcvbuf1dl,csize,mpi_logical,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhs /= mpi_proc_null ) then
            l(-j+1,1:nii,k) = rcvbuf1dl(1:nii)
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            sndbuf1dl(1:nii) = l(j,1:nii,k)
          end if
          call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%lhs,3, &
                          rcvbuf1dl,csize,mpi_logical,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhs /= mpi_proc_null ) then
            l(njj+j,1:nii,k) = rcvbuf1dl(1:nii)
          end if
        end do
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dl(icount) = l(njj-j+1,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%rht,5, &
                          rcvbuf1dl,bsize,mpi_logical,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,-i+1,k) = rcvbuf1dl(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dl(icount) = l(j,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%lhb,6, &
                          sndbuf1dl,bsize,mpi_logical,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,nii+i,k) = rcvbuf1dl(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dl(icount) = l(njj-j+1,i,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%rhb,7, &
                          sndbuf1dl,bsize,mpi_logical,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(-j+1,nii+i,k) = rcvbuf1dl(icount)
              icount = icount + 1
            end do
          end do
        end if
        if ( pspace%lht /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              sndbuf1dl(icount) = l(j,nii-i+1,k)
              icount = icount + 1
            end do
          end do
        end if
        call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%lht,8, &
                          sndbuf1dl,bsize,mpi_logical,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        if ( pspace%rhb /= mpi_proc_null ) then
          icount = 0
          do i = 1 , esize2
            do j = 1 , esize1
              l(njj+j,-i+1,k) = rcvbuf1dl(icount)
              icount = icount + 1
            end do
          end do
        end if
      end do
    end subroutine exchange_internal3d_l

    subroutine exchange_internal4d_d(l)
      implicit none
      real(dp) , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do t = nt1 , nt2
        do k = nk1 , nk2
          do i = 1 , esize2
            if ( pspace%btm /= mpi_proc_null ) then
              sndbuf1dd(1:njj) = l(1:njj,i,k,t)
            end if
            call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%btm,1, &
                          rcvbuf1dd,csize,mpi_real8,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%top /= mpi_proc_null ) then
              l(1:njj,nii+i,k,t) = rcvbuf1dd(1:njj)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              sndbuf1dd(1:njj) = l(1:njj,nii-i+1,k,t)
            end if
            call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%top,4, &
                          rcvbuf1dd,csize,mpi_real8,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%btm /= mpi_proc_null ) then
              l(1:njj,-i+1,k,t) = rcvbuf1dd(1:njj)
            end if
          end do
          do j = 1 , esize1
            if ( pspace%rhs /= mpi_proc_null ) then
              sndbuf1dd(1:nii) = l(njj-j+1,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%rhs,2, &
                          rcvbuf1dd,csize,mpi_real8,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%lhs /= mpi_proc_null ) then
              l(-j+1,1:nii,k,t) = rcvbuf1dd(1:nii)
            end if
            if ( pspace%lhs /= mpi_proc_null ) then
              sndbuf1dd(1:nii) = l(j,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1dd,csize,mpi_real8,pspace%lhs,3, &
                          rcvbuf1dd,csize,mpi_real8,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%rhs /= mpi_proc_null ) then
              l(njj+j,1:nii,k,t) = rcvbuf1dd(1:nii)
            end if
          end do
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dd(icount) = l(njj-j+1,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%rht,5, &
                          rcvbuf1dd,bsize,mpi_real8,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,-i+1,k,t) = rcvbuf1dd(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dd(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%lhb,6, &
                          rcvbuf1dd,bsize,mpi_real8,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,nii+i,k,t) = rcvbuf1dd(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dd(icount) = l(njj-j+1,i,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%rhb,7, &
                          rcvbuf1dd,bsize,mpi_real8,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,nii+i,k,t) = rcvbuf1dd(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dd(icount) = l(j,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dd,bsize,mpi_real8,pspace%lht,8, &
                          rcvbuf1dd,bsize,mpi_real8,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,-i+1,k,t) = rcvbuf1dd(icount)
                icount = icount + 1
              end do
            end do
          end if
        end do
      end do
    end subroutine exchange_internal4d_d

    subroutine exchange_internal4d_r(l)
      implicit none
      real(sp) , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do t = nt1 , nt2
        do k = nk1 , nk2
          do i = 1 , esize2
            if ( pspace%btm /= mpi_proc_null ) then
              sndbuf1dr(1:njj) = l(1:njj,i,k,t)
            end if
            call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%btm,1, &
                          rcvbuf1dr,csize,mpi_real4,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%top /= mpi_proc_null ) then
              l(1:njj,nii+i,k,t) = rcvbuf1dr(1:njj)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              sndbuf1dr(1:njj) = l(1:njj,nii-i+1,k,t)
            end if
            call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%top,4, &
                          rcvbuf1dr,csize,mpi_real4,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%btm /= mpi_proc_null ) then
              l(1:njj,-i+1,k,t) = rcvbuf1dr(1:njj)
            end if
          end do
          do j = 1 , esize1
            if ( pspace%rhs /= mpi_proc_null ) then
              sndbuf1dr(1:nii) = l(njj-j+1,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%rhs,2, &
                          rcvbuf1dr,csize,mpi_real4,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%lhs /= mpi_proc_null ) then
              l(-j+1,1:nii,k,t) = rcvbuf1dr(1:nii)
            end if
            if ( pspace%lhs /= mpi_proc_null ) then
              sndbuf1dr(1:nii) = l(j,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1dr,csize,mpi_real4,pspace%lhs,3, &
                          rcvbuf1dr,csize,mpi_real4,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%rhs /= mpi_proc_null ) then
              l(njj+j,1:nii,k,t) = rcvbuf1dr(1:nii)
            end if
          end do
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dr(icount) = l(njj-j+1,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%rht,5, &
                          rcvbuf1dr,bsize,mpi_real4,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,-i+1,k,t) = rcvbuf1dr(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dr(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%lhb,6, &
                          rcvbuf1dr,bsize,mpi_real4,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,nii+i,k,t) = rcvbuf1dr(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dr(icount) = l(njj-j+1,i,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%rhb,7, &
                          rcvbuf1dr,bsize,mpi_real4,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,nii+i,k,t) = rcvbuf1dr(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dr(icount) = l(j,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dr,bsize,mpi_real4,pspace%lht,8, &
                          rcvbuf1dr,bsize,mpi_real4,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,-i+1,k,t) = rcvbuf1dr(icount)
                icount = icount + 1
              end do
            end do
          end if
        end do
      end do
    end subroutine exchange_internal4d_r

    subroutine exchange_internal4d_i(l)
      implicit none
      integer , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do t = nt1 , nt2
        do k = nk1 , nk2
          do i = 1 , esize2
            if ( pspace%btm /= mpi_proc_null ) then
              sndbuf1di(1:njj) = l(1:njj,i,k,t)
            end if
            call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%btm,1, &
                          rcvbuf1di,csize,mpi_integer,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%top /= mpi_proc_null ) then
              l(1:njj,nii+i,k,t) = rcvbuf1di(1:njj)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              sndbuf1di(1:njj) = l(1:njj,nii-i+1,k,t)
            end if
            call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%top,4, &
                          rcvbuf1di,csize,mpi_integer,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%btm /= mpi_proc_null ) then
              l(1:njj,-i+1,k,t) = rcvbuf1di(1:njj)
            end if
          end do
          do j = 1 , esize1
            if ( pspace%rhs /= mpi_proc_null ) then
              sndbuf1di(1:nii) = l(njj-j+1,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%rhs,2, &
                          rcvbuf1di,csize,mpi_integer,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%lhs /= mpi_proc_null ) then
              l(-j+1,1:nii,k,t) = rcvbuf1di(1:nii)
            end if
            if ( pspace%lhs /= mpi_proc_null ) then
              sndbuf1di(1:nii) = l(j,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1di,csize,mpi_integer,pspace%lhs,3, &
                          rcvbuf1di,csize,mpi_integer,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%rhs /= mpi_proc_null ) then
              l(njj+j,1:nii,k,t) = rcvbuf1di(1:nii)
            end if
          end do
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1di(icount) = l(njj-j+1,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%rht,5, &
                          rcvbuf1di,bsize,mpi_integer,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,-i+1,k,t) = rcvbuf1di(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1di(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%lhb,6, &
                          rcvbuf1di,bsize,mpi_integer,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,nii+i,k,t) = rcvbuf1di(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1di(icount) = l(njj-j+1,i,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%rhb,7, &
                          rcvbuf1di,bsize,mpi_integer,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,nii+i,k,t) = rcvbuf1di(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1di(icount) = l(j,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1di,bsize,mpi_integer,pspace%lht,8, &
                          rcvbuf1di,bsize,mpi_integer,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,-i+1,k,t) = rcvbuf1di(icount)
                icount = icount + 1
              end do
            end do
          end if
        end do
      end do
    end subroutine exchange_internal4d_i

    subroutine exchange_internal4d_s(l)
      implicit none
      integer(2) , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do t = nt1 , nt2
        do k = nk1 , nk2
          do i = 1 , esize2
            if ( pspace%btm /= mpi_proc_null ) then
              sndbuf1ds(1:njj) = l(1:njj,i,k,t)
            end if
            call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%btm,1, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%top /= mpi_proc_null ) then
              l(1:njj,nii+i,k,t) = rcvbuf1ds(1:njj)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              sndbuf1ds(1:njj) = l(1:njj,nii-i+1,k,t)
            end if
            call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%top,4, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%btm /= mpi_proc_null ) then
              l(1:njj,-i+1,k,t) = rcvbuf1ds(1:njj)
            end if
          end do
          do j = 1 , esize1
            if ( pspace%rhs /= mpi_proc_null ) then
              sndbuf1ds(1:nii) = l(njj-j+1,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%rhs,2, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%lhs /= mpi_proc_null ) then
              l(-j+1,1:nii,k,t) = rcvbuf1ds(1:nii)
            end if
            if ( pspace%lhs /= mpi_proc_null ) then
              sndbuf1ds(1:nii) = l(j,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1ds,csize,mpi_integer2,pspace%lhs,3, &
                          rcvbuf1ds,csize,mpi_integer2,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%rhs /= mpi_proc_null ) then
              l(njj+j,1:nii,k,t) = rcvbuf1ds(1:nii)
            end if
          end do
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1ds(icount) = l(njj-j+1,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%rht,5, &
                          rcvbuf1ds,bsize,mpi_integer2,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,-i+1,k,t) = rcvbuf1ds(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1ds(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%lhb,6, &
                          rcvbuf1ds,bsize,mpi_integer2,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,nii+i,k,t) = rcvbuf1ds(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1ds(icount) = l(njj-j+1,i,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%rhb,7, &
                          rcvbuf1ds,bsize,mpi_integer2,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,nii+i,k,t) = rcvbuf1ds(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1ds(icount) = l(j,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1ds,bsize,mpi_integer2,pspace%lht,8, &
                          rcvbuf1ds,bsize,mpi_integer2,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,-i+1,k,t) = rcvbuf1ds(icount)
                icount = icount + 1
              end do
            end do
          end if
        end do
      end do
    end subroutine exchange_internal4d_s

    subroutine exchange_internal4d_l(l)
      implicit none
      logical , dimension(:,:,:,:) , pointer , intent(inout) :: l
      integer :: k , t , nii , njj , nk1 , nk2 , nt1 , nt2 , ierr
      integer :: icount , i , j , esize1 , esize2 , bsize
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling exchange_internal before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      if ( esize1 == 0 .and. esize2 == 0 ) return
      nk1 = lbound(l,3)
      nk2 = ubound(l,3)
      nt1 = lbound(l,4)
      nt2 = ubound(l,4)
      bsize = esize1*esize2
      njj = ubound(l,1) - esize1
      nii = ubound(l,2) - esize2
      do t = nt1 , nt2
        do k = nk1 , nk2
          do i = 1 , esize2
            if ( pspace%btm /= mpi_proc_null ) then
              sndbuf1dl(1:njj) = l(1:njj,i,k,t)
            end if
            call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%btm,1, &
                          rcvbuf1dl,csize,mpi_logical,pspace%top,1, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%top /= mpi_proc_null ) then
              l(1:njj,nii+i,k,t) = rcvbuf1dl(1:njj)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              sndbuf1dl(1:njj) = l(1:njj,nii-i+1,k,t)
            end if
            call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%top,4, &
                          rcvbuf1dl,csize,mpi_logical,pspace%btm,4, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%btm /= mpi_proc_null ) then
              l(1:njj,-i+1,k,t) = rcvbuf1dl(1:njj)
            end if
          end do
          do j = 1 , esize1
            if ( pspace%rhs /= mpi_proc_null ) then
              sndbuf1dl(1:nii) = l(njj-j+1,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%rhs,2, &
                          rcvbuf1dl,csize,mpi_logical,pspace%lhs,2, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%lhs /= mpi_proc_null ) then
              l(-j+1,1:nii,k,t) = rcvbuf1dl(1:nii)
            end if
            if ( pspace%lhs /= mpi_proc_null ) then
              sndbuf1dl(1:nii) = l(j,1:nii,k,t)
            end if
            call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%lhs,3, &
                          rcvbuf1dl,csize,mpi_logical,pspace%rhs,3, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
            if ( pspace%rhs /= mpi_proc_null ) then
              l(njj+j,1:nii,k,t) = rcvbuf1dl(1:nii)
            end if
          end do
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dl(icount) = l(njj-j+1,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dl,csize,mpi_logical,pspace%rht,5, &
                          rcvbuf1dl,csize,mpi_logical,pspace%lhb,5, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,-i+1,k,t) = rcvbuf1dl(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dl(icount) = l(j,i,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%lhb,6, &
                          rcvbuf1dl,bsize,mpi_logical,pspace%rht,6, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,nii+i,k,t) = rcvbuf1dl(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dl(icount) = l(njj-j+1,-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%rhb,7, &
                          rcvbuf1dl,bsize,mpi_logical,pspace%lht,7, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(-j+1,nii+i,k,t) = rcvbuf1dl(icount)
                icount = icount + 1
              end do
            end do
          end if
          if ( pspace%lht /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                sndbuf1dl(icount) = l(j,nii-i+1,k,t)
                icount = icount + 1
              end do
            end do
          end if
          call mpi_sendrecv(sndbuf1dl,bsize,mpi_logical,pspace%lht,8, &
                          rcvbuf1dl,bsize,mpi_logical,pspace%rhb,8, &
                          pspace%cartesian_communicator,mpi_status_ignore,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          if ( pspace%rhb /= mpi_proc_null ) then
            icount = 0
            do i = 1 , esize2
              do j = 1 , esize1
                l(njj+j,-i+1,k,t) = rcvbuf1dl(icount)
                icount = icount + 1
              end do
            end do
          end if
        end do
      end do
    end subroutine exchange_internal4d_l

    subroutine global_to_proc2d_d(g,l)
      implicit none
      real(dp) , pointer , dimension(:,:) , intent(in) :: g
      real(dp) , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_d

    subroutine proc_to_global2d_d(l,g)
      implicit none
      real(dp) , pointer , dimension(:,:) , intent(in) :: l
      real(dp) , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_d

    subroutine global_to_proc2d_r(g,l)
      implicit none
      real(sp) , pointer , dimension(:,:) , intent(in) :: g
      real(sp) , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_r

    subroutine proc_to_global2d_r(l,g)
      implicit none
      real(sp) , pointer , dimension(:,:) , intent(in) :: l
      real(sp) , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_r

    subroutine global_to_proc2d_i(g,l)
      implicit none
      integer , pointer , dimension(:,:) , intent(in) :: g
      integer , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_i

    subroutine proc_to_global2d_i(l,g)
      implicit none
      integer , pointer , dimension(:,:) , intent(in) :: l
      integer , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_i

    subroutine global_to_proc2d_s(g,l)
      implicit none
      integer(2) , pointer , dimension(:,:) , intent(in) :: g
      integer(2) , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_s

    subroutine proc_to_global2d_s(l,g)
      implicit none
      integer(2) , pointer , dimension(:,:) , intent(in) :: l
      integer(2) , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_s

    subroutine global_to_proc2d_l(g,l)
      implicit none
      logical , pointer , dimension(:,:) , intent(in) :: g
      logical , pointer , dimension(:,:) , intent(inout) :: l
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      l(1:njj,1:nii) = g(lbgjj:ubgjj,lbgii:ubgii)
    end subroutine global_to_proc2d_l

    subroutine proc_to_global2d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:) , intent(in) :: l
      logical , pointer , dimension(:,:) , intent(inout) :: g
      integer :: nii , njj , lbgii , lbgjj , ubgii , ubgjj , esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
      ubgii = lbgii+nii-1
      ubgjj = lbgjj+njj-1
      g(lbgjj:ubgjj,lbgii:ubgii) = l(1:njj,1:nii)
    end subroutine proc_to_global2d_l

    subroutine global_to_proc3d_d(g,l)
      implicit none
      real(dp) , pointer , dimension(:,:,:) , intent(in) :: g
      real(dp) , pointer , dimension(:,:,:) , intent(inout) :: l
      integer :: k , k1 , k2 , nii , njj , lbgii , lbgjj , ubgii , ubgjj
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_proc before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      integer :: esize1 , esize2
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling proc_to_global before domain_setup')
      end if
      t1 = lbound(l,4)
      t2 = ubound(l,4)
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      esize1 = 1-lbound(l,1)
      esize2 = 1-lbound(l,2)
      njj = ubound(l,1)-esize1
      nii = ubound(l,2)-esize2
      lbgii = pspace%g_is
      lbgjj = pspace%g_js
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
      call mpi_bcast(v,1,mpi_integer2,masterproc, &
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
      call mpi_bcast(v,size(v),mpi_integer2,masterproc, &
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
      real(dp) , pointer , dimension(:,:) , intent(in) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( associated(g) ) then
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
      real(sp) , pointer , dimension(:,:) , intent(in) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( associated(g) ) then
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
      integer , pointer , dimension(:,:) , intent(in) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( associated(g) ) then
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
      integer(2) , pointer , dimension(:,:) , intent(in) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( associated(g) ) then
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
            call mpi_send(mnode%excbuf2ds,mnode%pgsize(inode),mpi_integer2, &
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
        call mpi_recv(excbuf2ds,pspace%totalpoints,mpi_integer2,masterproc,0, &
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
      logical , pointer , dimension(:,:) , intent(in) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      if ( associated(g) ) then
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
      real(dp) , pointer , dimension(:,:,:) , intent(in) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
      real(sp) , pointer , dimension(:,:,:) , intent(in) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
      integer , pointer , dimension(:,:,:) , intent(in) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
      integer(2) , pointer , dimension(:,:,:) , intent(in) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
              call mpi_send(mnode%excbuf2ds,mnode%pgsize(inode),mpi_integer2, &
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
          call mpi_recv(excbuf2ds,pspace%totalpoints,mpi_integer2,masterproc, &
                        0,pspace%cartesian_communicator,mpi_status_ignore,ierr)
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
      logical , pointer , dimension(:,:,:) , intent(in) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
      real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: g
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
      if ( associated(g) ) then
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
      real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: g
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
      if ( associated(g) ) then
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
      integer , pointer , dimension(:,:,:,:) , intent(in) :: g
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
      if ( associated(g) ) then
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
      integer(2) , pointer , dimension(:,:,:,:) , intent(in) :: g
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
      if ( associated(g) ) then
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
                call mpi_send(mnode%excbuf2ds,mnode%pgsize(inode), &
                              mpi_integer2,inode,0, &
                              pspace%cartesian_communicator,ierr)
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
            call mpi_recv(excbuf2ds,pspace%totalpoints,mpi_integer2, &
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
      logical , pointer , dimension(:,:,:,:) , intent(in) :: g
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
      if ( associated(g) ) then
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
      real(dp) , pointer , dimension(:,:) , intent(inout) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( associated(g) ) then
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
      real(sp) , pointer , dimension(:,:) , intent(inout) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( associated(g) ) then
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
      integer , pointer , dimension(:,:) , intent(inout) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( associated(g) ) then
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
      integer(2) , pointer , dimension(:,:) , intent(inout) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( associated(g) ) then
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
            call mpi_recv(mnode%excbuf2ds,mnode%pgsize(inode),mpi_integer2, &
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
        call mpi_send(excbuf2ds,pspace%totalpoints,mpi_integer2,masterproc,0, &
                      pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine nodes_to_master2d_s

    subroutine nodes_to_master2d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:) , intent(in) :: l
      logical , pointer , dimension(:,:) , intent(inout) :: g
      integer :: inode , icount , i , j
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling nodes_to_master before domain_setup')
      end if
      if ( associated(g) ) then
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
      real(dp) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
      real(sp) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
      integer , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
      integer(2) , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
              call mpi_recv(mnode%excbuf2ds,mnode%pgsize(inode),mpi_integer2, &
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
          call mpi_send(excbuf2ds,pspace%totalpoints,mpi_integer2,masterproc, &
                        0,pspace%cartesian_communicator,ierr)
          if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        end do
      end if
    end subroutine nodes_to_master3d_s

    subroutine nodes_to_master3d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:,:) , intent(in) :: l
      logical , pointer , dimension(:,:,:) , intent(inout) :: g
      integer :: inode , icount , i , j , k , k1 , k2 , gk1 , gk2
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling master_to_nodes before domain_setup')
      end if
      k1 = lbound(l,3)
      k2 = ubound(l,3)
      if ( associated(g) ) then
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
      real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: g
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
      if ( associated(g) ) then
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
      real(sp) , pointer , dimension(:,:,:,:) , intent(inout) :: g
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
      if ( associated(g) ) then
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
      integer , pointer , dimension(:,:,:,:) , intent(inout) :: g
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
      if ( associated(g) ) then
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
      integer(2) , pointer , dimension(:,:,:,:) , intent(inout) :: g
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
      if ( associated(g) ) then
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
                call mpi_recv(mnode%excbuf2ds,mnode%pgsize(inode), &
                              mpi_integer2,inode,0, &
                              pspace%cartesian_communicator, &
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
            call mpi_send(excbuf2ds,pspace%totalpoints,mpi_integer2, &
                          masterproc,0,pspace%cartesian_communicator,ierr)
            if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
          end do
        end do
      end if
    end subroutine nodes_to_master4d_s

    subroutine nodes_to_master4d_l(l,g)
      implicit none
      logical , pointer , dimension(:,:,:,:) , intent(in) :: l
      logical , pointer , dimension(:,:,:,:) , intent(inout) :: g
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
      if ( associated(g) ) then
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

    subroutine getbdy2d_d(g,lstagger)
      implicit none
      type(global_boundary2d_d) , intent(inout) :: g
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem2d(g%north,j1,j2,1,bdysize,__FILE__)
      call getmem2d(g%south,j1,j2,1,bdysize,__FILE__)
      if ( .not. jperiodic ) then
        call getmem2d(g%east,1,bdysize,i1,i2,__FILE__)
        call getmem2d(g%west,1,bdysize,i1,i2,__FILE__)
      end if
    end subroutine getbdy2d_d

    subroutine getbdy3d_d(g,k1,k2,lstagger)
      implicit none
      type(global_boundary3d_d) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem3d(g%north,j1,j2,1,bdysize,k1,k2,__FILE__)
      call getmem3d(g%south,j1,j2,1,bdysize,k1,k2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem3d(g%east,1,bdysize,i1,i2,k1,k2,__FILE__)
        call getmem3d(g%west,1,bdysize,i1,i2,k1,k2,__FILE__)
      end if
    end subroutine getbdy3d_d

    subroutine getbdy4d_d(g,k1,k2,t1,t2,lstagger)
      implicit none
      type(global_boundary4d_d) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem4d(g%north,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      call getmem4d(g%south,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem4d(g%east,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
        call getmem4d(g%west,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
      end if
    end subroutine getbdy4d_d

    subroutine getbdy2d_r(g,lstagger)
      implicit none
      type(global_boundary2d_r) , intent(inout) :: g
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem2d(g%north,j1,j2,1,bdysize,__FILE__)
      call getmem2d(g%south,j1,j2,1,bdysize,__FILE__)
      if ( .not. jperiodic ) then
        call getmem2d(g%east,1,bdysize,i1,i2,__FILE__)
        call getmem2d(g%west,1,bdysize,i1,i2,__FILE__)
      end if
    end subroutine getbdy2d_r

    subroutine getbdy3d_r(g,k1,k2,lstagger)
      implicit none
      type(global_boundary3d_r) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem3d(g%north,j1,j2,1,bdysize,k1,k2,__FILE__)
      call getmem3d(g%south,j1,j2,1,bdysize,k1,k2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem3d(g%east,1,bdysize,i1,i2,k1,k2,__FILE__)
        call getmem3d(g%west,1,bdysize,i1,i2,k1,k2,__FILE__)
      end if
    end subroutine getbdy3d_r

    subroutine getbdy4d_r(g,k1,k2,t1,t2,lstagger)
      implicit none
      type(global_boundary4d_r) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem4d(g%north,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      call getmem4d(g%south,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem4d(g%east,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
        call getmem4d(g%west,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
      end if
    end subroutine getbdy4d_r

    subroutine getbdy2d_i(g,lstagger)
      implicit none
      type(global_boundary2d_i) , intent(inout) :: g
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem2d(g%north,j1,j2,1,bdysize,__FILE__)
      call getmem2d(g%south,j1,j2,1,bdysize,__FILE__)
      if ( .not. jperiodic ) then
        call getmem2d(g%east,1,bdysize,i1,i2,__FILE__)
        call getmem2d(g%west,1,bdysize,i1,i2,__FILE__)
      end if
    end subroutine getbdy2d_i

    subroutine getbdy3d_i(g,k1,k2,lstagger)
      implicit none
      type(global_boundary3d_i) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem3d(g%north,j1,j2,1,bdysize,k1,k2,__FILE__)
      call getmem3d(g%south,j1,j2,1,bdysize,k1,k2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem3d(g%east,1,bdysize,i1,i2,k1,k2,__FILE__)
        call getmem3d(g%west,1,bdysize,i1,i2,k1,k2,__FILE__)
      end if
    end subroutine getbdy3d_i

    subroutine getbdy4d_i(g,k1,k2,t1,t2,lstagger)
      implicit none
      type(global_boundary4d_i) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem4d(g%north,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      call getmem4d(g%south,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem4d(g%east,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
        call getmem4d(g%west,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
      end if
    end subroutine getbdy4d_i

    subroutine getbdy2d_s(g,lstagger)
      implicit none
      type(global_boundary2d_s) , intent(inout) :: g
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem2d(g%north,j1,j2,1,bdysize,__FILE__)
      call getmem2d(g%south,j1,j2,1,bdysize,__FILE__)
      if ( .not. jperiodic ) then
        call getmem2d(g%east,1,bdysize,i1,i2,__FILE__)
        call getmem2d(g%west,1,bdysize,i1,i2,__FILE__)
      end if
    end subroutine getbdy2d_s

    subroutine getbdy3d_s(g,k1,k2,lstagger)
      implicit none
      type(global_boundary3d_s) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem3d(g%north,j1,j2,1,bdysize,k1,k2,__FILE__)
      call getmem3d(g%south,j1,j2,1,bdysize,k1,k2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem3d(g%east,1,bdysize,i1,i2,k1,k2,__FILE__)
        call getmem3d(g%west,1,bdysize,i1,i2,k1,k2,__FILE__)
      end if
    end subroutine getbdy3d_s

    subroutine getbdy4d_s(g,k1,k2,t1,t2,lstagger)
      implicit none
      type(global_boundary4d_s) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem4d(g%north,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      call getmem4d(g%south,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem4d(g%east,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
        call getmem4d(g%west,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
      end if
    end subroutine getbdy4d_s

    subroutine getbdy2d_l(g,lstagger)
      implicit none
      type(global_boundary2d_l) , intent(inout) :: g
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem2d(g%north,j1,j2,1,bdysize,__FILE__)
      call getmem2d(g%south,j1,j2,1,bdysize,__FILE__)
      if ( .not. jperiodic ) then
        call getmem2d(g%east,1,bdysize,i1,i2,__FILE__)
        call getmem2d(g%west,1,bdysize,i1,i2,__FILE__)
      end if
    end subroutine getbdy2d_l

    subroutine getbdy3d_l(g,k1,k2,lstagger)
      implicit none
      type(global_boundary3d_l) , intent(inout) :: g
      integer , intent(in) :: k1 , k2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem3d(g%north,j1,j2,1,bdysize,k1,k2,__FILE__)
      call getmem3d(g%south,j1,j2,1,bdysize,k1,k2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem3d(g%east,1,bdysize,i1,i2,k1,k2,__FILE__)
        call getmem3d(g%west,1,bdysize,i1,i2,k1,k2,__FILE__)
      end if
    end subroutine getbdy3d_l

    subroutine getbdy4d_l(g,k1,k2,t1,t2,lstagger)
      implicit none
      type(global_boundary4d_l) , intent(inout) :: g
      integer , intent(in) :: k1 , k2 , t1 , t2
      logical , intent(in) , optional :: lstagger
      logical :: isstagger = .false.
      integer :: i1 , i2 , j1 , j2
      if ( present(lstagger) )  isstagger  = lstagger
      call getextrema(global_grid,isstagger,ghost_none,i1,i2,j1,j2)
      call getmem4d(g%north,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      call getmem4d(g%south,j1,j2,1,bdysize,k1,k2,t1,t2,__FILE__)
      if ( .not. jperiodic ) then
        call getmem4d(g%east,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
        call getmem4d(g%west,1,bdysize,i1,i2,k1,k2,t1,t2,__FILE__)
      end if
    end subroutine getbdy4d_l

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

    subroutine set_external2d_d(l,bdy)
      implicit none
      real(dp) , dimension(:,:) , pointer , intent(inout) :: l
      type(global_boundary2d_d) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2) = bdy%south(pspace%g_js:pspace%g_je,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2) = bdy%north(pspace%g_js:pspace%g_je,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i) = bdy%west(:,pspace%g_is:pspace%g_ie)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i) = bdy%east(:,pspace%g_is:pspace%g_ie)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2) = bdy%south(pspace%g_js-1,:)
            else
              l(0,btmbdy1:btmbdy2) = bdy%south(gspace%g_j,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(1,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2) = bdy%north(pspace%g_js-1,:)
            else
              l(0,topbdy1:topbdy2) = bdy%north(gspace%g_j,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(1,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
        end if
      end if
    end subroutine set_external2d_d

    subroutine set_external2d_r(l,bdy)
      implicit none
      real(sp) , dimension(:,:) , pointer , intent(inout) :: l
      type(global_boundary2d_r) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2) = bdy%south(pspace%g_js:pspace%g_je,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2) = bdy%north(pspace%g_js:pspace%g_je,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i) = bdy%west(:,pspace%g_is:pspace%g_ie)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i) = bdy%east(:,pspace%g_is:pspace%g_ie)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2) = bdy%south(pspace%g_js-1,:)
            else
              l(0,btmbdy1:btmbdy2) = bdy%south(gspace%g_j,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(1,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2) = bdy%north(pspace%g_js-1,:)
            else
              l(0,topbdy1:topbdy2) = bdy%north(gspace%g_j,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(1,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
        end if
      end if
    end subroutine set_external2d_r

    subroutine set_external2d_i(l,bdy)
      implicit none
      integer , dimension(:,:) , pointer , intent(inout) :: l
      type(global_boundary2d_i) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2) = bdy%south(pspace%g_js:pspace%g_je,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2) = bdy%north(pspace%g_js:pspace%g_je,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i) = bdy%west(:,pspace%g_is:pspace%g_ie)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i) = bdy%east(:,pspace%g_is:pspace%g_ie)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2) = bdy%south(pspace%g_js-1,:)
            else
              l(0,btmbdy1:btmbdy2) = bdy%south(gspace%g_j,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(1,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2) = bdy%north(pspace%g_js-1,:)
            else
              l(0,topbdy1:topbdy2) = bdy%north(gspace%g_j,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(1,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
        end if
      end if
    end subroutine set_external2d_i

    subroutine set_external2d_s(l,bdy)
      implicit none
      integer(2) , dimension(:,:) , pointer , intent(inout) :: l
      type(global_boundary2d_s) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2) = bdy%south(pspace%g_js:pspace%g_je,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2) = bdy%north(pspace%g_js:pspace%g_je,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i) = bdy%west(:,pspace%g_is:pspace%g_ie)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i) = bdy%east(:,pspace%g_is:pspace%g_ie)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2) = bdy%south(pspace%g_js-1,:)
            else
              l(0,btmbdy1:btmbdy2) = bdy%south(gspace%g_j,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(1,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2) = bdy%north(pspace%g_js-1,:)
            else
              l(0,topbdy1:topbdy2) = bdy%north(gspace%g_j,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(1,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
        end if
      end if
    end subroutine set_external2d_s

    subroutine set_external2d_l(l,bdy)
      implicit none
      logical , dimension(:,:) , pointer , intent(inout) :: l
      type(global_boundary2d_l) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2) = bdy%south(pspace%g_js:pspace%g_je,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2) = bdy%north(pspace%g_js:pspace%g_je,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i) = bdy%west(:,pspace%g_is:pspace%g_ie)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i) = bdy%east(:,pspace%g_is:pspace%g_ie)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2) = bdy%south(pspace%g_js-1,:)
            else
              l(0,btmbdy1:btmbdy2) = bdy%south(gspace%g_j,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2) = bdy%south(1,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2) = bdy%north(pspace%g_js-1,:)
            else
              l(0,topbdy1:topbdy2) = bdy%north(gspace%g_j,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(pspace%g_je+1,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2) = bdy%north(1,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0) = bdy%west(:,pspace%g_is-1)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1) = bdy%east(:,pspace%g_ie+1)
            end if
          end if
        end if
      end if
    end subroutine set_external2d_l

    subroutine set_external3d_d(l,bdy)
      implicit none
      real(dp) , dimension(:,:,:) , pointer , intent(inout) :: l
      type(global_boundary3d_d) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_js-1,:,:)
            else
              l(0,btmbdy1:btmbdy2,:) = bdy%south(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(1,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:) = bdy%north(pspace%g_js-1,:,:)
            else
              l(0,topbdy1:topbdy2,:) = bdy%north(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(1,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
        end if
      end if
    end subroutine set_external3d_d

    subroutine set_external3d_r(l,bdy)
      implicit none
      real(sp) , dimension(:,:,:) , pointer , intent(inout) :: l
      type(global_boundary3d_r) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_js-1,:,:)
            else
              l(0,btmbdy1:btmbdy2,:) = bdy%south(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(1,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:) = bdy%north(pspace%g_js-1,:,:)
            else
              l(0,topbdy1:topbdy2,:) = bdy%north(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(1,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
        end if
      end if
    end subroutine set_external3d_r

    subroutine set_external3d_i(l,bdy)
      implicit none
      integer , dimension(:,:,:) , pointer , intent(inout) :: l
      type(global_boundary3d_i) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_js-1,:,:)
            else
              l(0,btmbdy1:btmbdy2,:) = bdy%south(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(1,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:) = bdy%north(pspace%g_js-1,:,:)
            else
              l(0,topbdy1:topbdy2,:) = bdy%north(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(1,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
        end if
      end if
    end subroutine set_external3d_i

    subroutine set_external3d_s(l,bdy)
      implicit none
      integer(2) , dimension(:,:,:) , pointer , intent(inout) :: l
      type(global_boundary3d_s) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_js-1,:,:)
            else
              l(0,btmbdy1:btmbdy2,:) = bdy%south(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(1,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:) = bdy%north(pspace%g_js-1,:,:)
            else
              l(0,topbdy1:topbdy2,:) = bdy%north(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(1,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
        end if
      end if
    end subroutine set_external3d_s

    subroutine set_external3d_l(l,bdy)
      implicit none
      logical , dimension(:,:,:) , pointer , intent(inout) :: l
      type(global_boundary3d_l) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_js-1,:,:)
            else
              l(0,btmbdy1:btmbdy2,:) = bdy%south(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:) = bdy%south(1,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:) = bdy%north(pspace%g_js-1,:,:)
            else
              l(0,topbdy1:topbdy2,:) = bdy%north(gspace%g_j,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(pspace%g_je+1,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:) = bdy%north(1,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:) = bdy%west(:,pspace%g_is-1,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:) = bdy%east(:,pspace%g_ie+1,:)
            end if
          end if
        end if
      end if
    end subroutine set_external3d_l

    subroutine set_external4d_d(l,bdy)
      implicit none
      real(dp) , dimension(:,:,:,:) , pointer , intent(inout) :: l
      type(global_boundary4d_d) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(pspace%g_js-1,:,:,:)
            else
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = &
                          bdy%south(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = bdy%south(1,:,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:,:) = bdy%north(pspace%g_js-1,:,:,:)
            else
              l(0,topbdy1:topbdy2,:,:) = bdy%north(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = &
                          bdy%north(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = bdy%north(1,:,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
        end if
      end if
    end subroutine set_external4d_d

    subroutine set_external4d_r(l,bdy)
      implicit none
      real(sp) , dimension(:,:,:,:) , pointer , intent(inout) :: l
      type(global_boundary4d_r) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(pspace%g_js-1,:,:,:)
            else
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = &
                          bdy%south(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = bdy%south(1,:,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:,:) = bdy%north(pspace%g_js-1,:,:,:)
            else
              l(0,topbdy1:topbdy2,:,:) = bdy%north(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = &
                          bdy%north(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = bdy%north(1,:,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
        end if
      end if
    end subroutine set_external4d_r

    subroutine set_external4d_i(l,bdy)
      implicit none
      integer , dimension(:,:,:,:) , pointer , intent(inout) :: l
      type(global_boundary4d_i) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(pspace%g_js-1,:,:,:)
            else
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = &
                          bdy%south(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = bdy%south(1,:,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:,:) = bdy%north(pspace%g_js-1,:,:,:)
            else
              l(0,topbdy1:topbdy2,:,:) = bdy%north(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = &
                          bdy%north(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = bdy%north(1,:,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
        end if
      end if
    end subroutine set_external4d_i

    subroutine set_external4d_s(l,bdy)
      implicit none
      integer(2) , dimension(:,:,:,:) , pointer , intent(inout) :: l
      type(global_boundary4d_s) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(pspace%g_js-1,:,:,:)
            else
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = &
                          bdy%south(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = bdy%south(1,:,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:,:) = bdy%north(pspace%g_js-1,:,:,:)
            else
              l(0,topbdy1:topbdy2,:,:) = bdy%north(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = &
                          bdy%north(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = bdy%north(1,:,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
        end if
      end if
    end subroutine set_external4d_s

    subroutine set_external4d_l(l,bdy)
      implicit none
      logical , dimension(:,:,:,:) , pointer , intent(inout) :: l
      type(global_boundary4d_l) , intent(in) :: bdy
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__,'Calling set_externa before domain_setup')
      end if
      ! internal points
      if ( pspace%btm == mpi_proc_null ) then
        l(1:pspace%p_j,btmbdy1:btmbdy2,:,:) = &
                      bdy%south(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( pspace%top == mpi_proc_null ) then
        l(1:pspace%p_j,topbdy1:topbdy2,:,:) = &
                      bdy%north(pspace%g_js:pspace%g_je,:,:,:)
      end if
      if ( .not. jperiodic ) then
        if ( pspace%lhs == mpi_proc_null ) then
          l(lhsbdy1:lhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%west(:,pspace%g_is:pspace%g_ie,:,:)
        end if
        if ( pspace%rhs == mpi_proc_null ) then
          l(rhsbdy1:rhsbdy2,1:pspace%p_i,:,:) = &
                      bdy%east(:,pspace%g_is:pspace%g_ie,:,:)
        end if
      end if
      ! exchange grid
      if ( lbound(l,1) == 0 ) then
        if ( pspace%btm == mpi_proc_null ) then
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(pspace%g_js-1,:,:,:)
            else
              l(0,btmbdy1:btmbdy2,:,:) = bdy%south(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = &
                          bdy%south(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,btmbdy1:btmbdy2,:,:) = bdy%south(1,:,:,:)
            end if
          end if
        end if
        if ( pspace%top == mpi_proc_null ) then
          if ( pspace%lhs /= mpi_proc_null ) then
            if ( pspace%g_js-1 > 1 ) then
              l(0,topbdy1:topbdy2,:,:) = bdy%north(pspace%g_js-1,:,:,:)
            else
              l(0,topbdy1:topbdy2,:,:) = bdy%north(gspace%g_j,:,:,:)
            end if
          end if
          if ( pspace%rhs /= mpi_proc_null ) then
            if ( pspace%g_je+1 <= gspace%g_j ) then
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = &
                          bdy%north(pspace%g_je+1,:,:,:)
            else
              l(pspace%p_j+1,topbdy1:topbdy2,:,:) = bdy%north(1,:,:,:)
            end if
          end if
        end if
        if ( .not. jperiodic ) then
          if ( pspace%lhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(lhsbdy1:lhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
          if ( pspace%rhs == mpi_proc_null ) then
            if ( pspace%btm /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,0,:,:) = bdy%west(:,pspace%g_is-1,:,:)
            end if
            if ( pspace%top /= mpi_proc_null ) then
              l(rhsbdy1:rhsbdy2,pspace%p_i+1,:,:) = &
                          bdy%east(:,pspace%g_ie+1,:,:)
            end if
          end if
        end if
      end if
    end subroutine set_external4d_l

    subroutine global_to_globbdy2d_d(bdy,g)
      implicit none
      type(global_boundary2d_d) , intent(out) :: bdy
      real(dp) , dimension(:,:) , pointer , intent(in) :: g
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:) = g(:,gtopbdy1:gtopbdy2)
        bdy%south(:,:) = g(:,gbtmbdy1:gbtmbdy2)
        if ( .not. jperiodic ) then
          bdy%west(:,:) = g(glhsbdy1:glhsbdy2,:)
          bdy%east(:,:) = g(grhsbdy1:grhsbdy2,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*bdysize,mpi_real8,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*bdysize,mpi_real8,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*bdysize,mpi_real8,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*bdysize,mpi_real8,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy2d_d

    subroutine global_to_globbdy2d_r(bdy,g)
      implicit none
      type(global_boundary2d_r) , intent(out) :: bdy
      real(sp) , dimension(:,:) , pointer , intent(in) :: g
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:) = g(:,gtopbdy1:gtopbdy2)
        bdy%south(:,:) = g(:,gbtmbdy1:gbtmbdy2)
        if ( .not. jperiodic ) then
          bdy%west(:,:) = g(glhsbdy1:glhsbdy2,:)
          bdy%east(:,:) = g(grhsbdy1:grhsbdy2,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*bdysize,mpi_real4,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*bdysize,mpi_real4,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*bdysize,mpi_real4,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*bdysize,mpi_real4,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy2d_r

    subroutine global_to_globbdy2d_i(bdy,g)
      implicit none
      type(global_boundary2d_i) , intent(out) :: bdy
      integer , dimension(:,:) , pointer , intent(in) :: g
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:) = g(:,gtopbdy1:gtopbdy2)
        bdy%south(:,:) = g(:,gbtmbdy1:gbtmbdy2)
        if ( .not. jperiodic ) then
          bdy%west(:,:) = g(glhsbdy1:glhsbdy2,:)
          bdy%east(:,:) = g(grhsbdy1:grhsbdy2,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*bdysize,mpi_integer,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*bdysize,mpi_integer,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic) then
        call mpi_bcast(bdy%east,gspace%g_j*bdysize,mpi_integer,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*bdysize,mpi_integer,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy2d_i

    subroutine global_to_globbdy2d_s(bdy,g)
      implicit none
      type(global_boundary2d_s) , intent(out) :: bdy
      integer(2) , dimension(:,:) , pointer , intent(in) :: g
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:) = g(:,gtopbdy1:gtopbdy2)
        bdy%south(:,:) = g(:,gbtmbdy1:gbtmbdy2)
        if ( .not. jperiodic ) then
          bdy%west(:,:) = g(glhsbdy1:glhsbdy2,:)
          bdy%east(:,:) = g(grhsbdy1:grhsbdy2,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*bdysize,mpi_integer2,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*bdysize,mpi_integer2,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*bdysize,mpi_integer2,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*bdysize,mpi_integer2,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy2d_s

    subroutine global_to_globbdy2d_l(bdy,g)
      implicit none
      type(global_boundary2d_l) , intent(out) :: bdy
      logical , dimension(:,:) , pointer , intent(in) :: g
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g)) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:) = g(:,gtopbdy1:gtopbdy2)
        bdy%south(:,:) = g(:,gbtmbdy1:gbtmbdy2)
        if ( .not. jperiodic ) then
          bdy%west(:,:) = g(glhsbdy1:glhsbdy2,:)
          bdy%east(:,:) = g(grhsbdy1:grhsbdy2,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*bdysize,mpi_logical,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*bdysize,mpi_logical,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*bdysize,mpi_logical,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*bdysize,mpi_logical,masterproc,  &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy2d_l

    subroutine global_to_globbdy3d_d(bdy,g)
      implicit none
      type(global_boundary3d_d) , intent(out) :: bdy
      real(dp) , dimension(:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      nk = size(g,3)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:) = g(:,gtopbdy1:gtopbdy2,:)
        bdy%south(:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:) = g(glhsbdy1:glhsbdy2,:,:)
          bdy%east(:,:,:) = g(grhsbdy1:grhsbdy2,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*bdysize,mpi_real8,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*bdysize,mpi_real8,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*bdysize,mpi_real8,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*bdysize,mpi_real8,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy3d_d

    subroutine global_to_globbdy3d_r(bdy,g)
      implicit none
      type(global_boundary3d_r) , intent(out) :: bdy
      real(sp) , dimension(:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      nk = size(g,3)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:) = g(:,gtopbdy1:gtopbdy2,:)
        bdy%south(:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:) = g(glhsbdy1:glhsbdy2,:,:)
          bdy%east(:,:,:) = g(grhsbdy1:grhsbdy2,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*bdysize,mpi_real4,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*bdysize,mpi_real4,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*bdysize,mpi_real4,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*bdysize,mpi_real4,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy3d_r

    subroutine global_to_globbdy3d_i(bdy,g)
      implicit none
      type(global_boundary3d_i) , intent(out) :: bdy
      integer , dimension(:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      nk = size(g,3)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:) = g(:,gtopbdy1:gtopbdy2,:)
        bdy%south(:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:) = g(glhsbdy1:glhsbdy2,:,:)
          bdy%east(:,:,:) = g(grhsbdy1:grhsbdy2,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*bdysize,mpi_integer,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*bdysize,mpi_integer,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*bdysize,mpi_integer,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*bdysize,mpi_integer,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy3d_i

    subroutine global_to_globbdy3d_s(bdy,g)
      implicit none
      type(global_boundary3d_s) , intent(out) :: bdy
      integer(2) , dimension(:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      nk = size(g,3)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:) = g(:,gtopbdy1:gtopbdy2,:)
        bdy%south(:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:) = g(glhsbdy1:glhsbdy2,:,:)
          bdy%east(:,:,:) = g(grhsbdy1:grhsbdy2,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*bdysize,mpi_integer2,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*bdysize,mpi_integer2,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*bdysize,mpi_integer2,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*bdysize,mpi_integer2,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy3d_s

    subroutine global_to_globbdy3d_l(bdy,g)
      implicit none
      type(global_boundary3d_l) , intent(out) :: bdy
      logical , dimension(:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      nk = size(g,3)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:) = g(:,gtopbdy1:gtopbdy2,:)
        bdy%south(:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:) = g(glhsbdy1:glhsbdy2,:,:)
          bdy%east(:,:,:) = g(grhsbdy1:grhsbdy2,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*bdysize,mpi_logical,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*bdysize,mpi_logical,masterproc, &
                     pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*bdysize,mpi_logical,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*bdysize,mpi_logical,masterproc, &
                       pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy3d_l

    subroutine global_to_globbdy4d_d(bdy,g)
      implicit none
      type(global_boundary4d_d) , intent(out) :: bdy
      real(dp) , dimension(:,:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk , t1 , t2 , nt
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      t1 = lbound(g,4)
      t2 = ubound(g,4)
      nk = size(g,3)
      nt = size(g,4)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( t1 /= lbound(bdy%north,4) .or. t2 /= ubound(bdy%north,4) ) then
        call fatal(__FILE__,__LINE__, &
          'Fourth dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:,:) = g(:,gtopbdy1:gtopbdy2,:,:)
        bdy%south(:,:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:,:) = g(glhsbdy1:glhsbdy2,:,:,:)
          bdy%east(:,:,:,:) = g(grhsbdy1:grhsbdy2,:,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*nt*bdysize,mpi_real8, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*nt*bdysize,mpi_real8, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*nt*bdysize,mpi_real8, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*nt*bdysize,mpi_real8, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy4d_d

    subroutine global_to_globbdy4d_r(bdy,g)
      implicit none
      type(global_boundary4d_r) , intent(out) :: bdy
      real(sp) , dimension(:,:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk , t1 , t2 , nt
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      t1 = lbound(g,4)
      t2 = ubound(g,4)
      nk = size(g,3)
      nt = size(g,4)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( t1 /= lbound(bdy%north,4) .or. t2 /= ubound(bdy%north,4) ) then
        call fatal(__FILE__,__LINE__, &
          'Fourth dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:,:) = g(:,gtopbdy1:gtopbdy2,:,:)
        bdy%south(:,:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:,:) = g(glhsbdy1:glhsbdy2,:,:,:)
          bdy%east(:,:,:,:) = g(grhsbdy1:grhsbdy2,:,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*nt*bdysize,mpi_real4, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*nt*bdysize,mpi_real4, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*nt*bdysize,mpi_real4, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*nt*bdysize,mpi_real4, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy4d_r

    subroutine global_to_globbdy4d_i(bdy,g)
      implicit none
      type(global_boundary4d_i) , intent(out) :: bdy
      integer , dimension(:,:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk , t1 , t2 , nt
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      t1 = lbound(g,4)
      t2 = ubound(g,4)
      nk = size(g,3)
      nt = size(g,4)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( t1 /= lbound(bdy%north,4) .or. t2 /= ubound(bdy%north,4) ) then
        call fatal(__FILE__,__LINE__, &
          'Fourth dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:,:) = g(:,gtopbdy1:gtopbdy2,:,:)
        bdy%south(:,:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:,:) = g(glhsbdy1:glhsbdy2,:,:,:)
          bdy%east(:,:,:,:) = g(grhsbdy1:grhsbdy2,:,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*nt*bdysize,mpi_integer, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*nt*bdysize,mpi_integer, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*nt*bdysize,mpi_integer, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*nt*bdysize,mpi_integer, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy4d_i

    subroutine global_to_globbdy4d_s(bdy,g)
      implicit none
      type(global_boundary4d_s) , intent(out) :: bdy
      integer(2) , dimension(:,:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk , t1 , t2 , nt
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      t1 = lbound(g,4)
      t2 = ubound(g,4)
      nk = size(g,3)
      nt = size(g,4)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( t1 /= lbound(bdy%north,4) .or. t2 /= ubound(bdy%north,4) ) then
        call fatal(__FILE__,__LINE__, &
          'Fourth dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:,:) = g(:,gtopbdy1:gtopbdy2,:,:)
        bdy%south(:,:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:,:) = g(glhsbdy1:glhsbdy2,:,:,:)
          bdy%east(:,:,:,:) = g(grhsbdy1:grhsbdy2,:,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*nt*bdysize,mpi_integer2, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*nt*bdysize,mpi_integer2, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*nt*bdysize,mpi_integer2, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*nt*bdysize,mpi_integer2, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy4d_s

    subroutine global_to_globbdy4d_l(bdy,g)
      implicit none
      type(global_boundary4d_l) , intent(out) :: bdy
      logical , dimension(:,:,:,:) , pointer , intent(in) :: g
      integer :: k1 , k2 , nk , t1 , t2 , nt
      integer :: ierr
      if ( .not. is_setup ) then
        call fatal(__FILE__,__LINE__, &
          'Calling global_to_globbdy before domain_setup')
      end if
      k1 = lbound(g,3)
      k2 = ubound(g,3)
      t1 = lbound(g,4)
      t2 = ubound(g,4)
      nk = size(g,3)
      nt = size(g,4)
      if ( k1 /= lbound(bdy%north,3) .or. k2 /= ubound(bdy%north,3) ) then
        call fatal(__FILE__,__LINE__, &
          'Vertical dimension mismatch in global_to_globbdy')
      end if
      if ( t1 /= lbound(bdy%north,4) .or. t2 /= ubound(bdy%north,4) ) then
        call fatal(__FILE__,__LINE__, &
          'Fourth dimension mismatch in global_to_globbdy')
      end if
      if ( associated(g) ) then
        if ( .not. am_i_master( ) ) then
          call fatal(__FILE__,__LINE__, &
          'Non master node using global_to_globbdy as master')
        end if
        if ( product(shape(g(:,:,1,1))) /= gsize ) then
          call fatal(__FILE__,__LINE__, &
          'Non global data in global_to_globbdy from master')
        end if
        bdy%north(:,:,:,:) = g(:,gtopbdy1:gtopbdy2,:,:)
        bdy%south(:,:,:,:) = g(:,gbtmbdy1:gbtmbdy2,:,:)
        if ( .not. jperiodic ) then
          bdy%west(:,:,:,:) = g(glhsbdy1:glhsbdy2,:,:,:)
          bdy%east(:,:,:,:) = g(grhsbdy1:grhsbdy2,:,:,:)
        end if
      end if
      call mpi_bcast(bdy%north,gspace%g_i*nk*nt*bdysize,mpi_logical, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      call mpi_bcast(bdy%south,gspace%g_i*nk*nt*bdysize,mpi_logical, &
                     masterproc,pspace%cartesian_communicator,ierr)
      if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      if ( .not. jperiodic ) then
        call mpi_bcast(bdy%east,gspace%g_j*nk*nt*bdysize,mpi_logical, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
        call mpi_bcast(bdy%west,gspace%g_j*nk*nt*bdysize,mpi_logical, &
                       masterproc,pspace%cartesian_communicator,ierr)
        if ( ierr /= mpi_success ) call mpi_fatal(__FILE__,__LINE__,ierr)
      end if
    end subroutine global_to_globbdy4d_l

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

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Regression driver for Main/mpplib/mod_mppparam.F90.
!
! Calls every public generic of the module, for every kind and rank it
! provides, and writes each resulting array (halos included) to a per
! processor file. Running it with two builds of the library and comparing
! the dumps record by record (compare_mpp_regression.py) shows any change
! in behaviour; run_mpp_regression.sh does both.
!
! Usage: mpirun -np N mpp_regression <outdir> <case>
!   case 0: limited area, nsg = 1    case 1: limited area, nsg = 2
!   case 2: i_band,       nsg = 1    case 3: i_crm,        nsg = 2
!

module mpp_regression_dump
  use mod_intkinds
  use mod_realkinds
  implicit none
  integer :: iu = -1
  interface dump
    module procedure d_r8_1, d_r8_2, d_r8_3, d_r8_4
    module procedure d_r4_1, d_r4_2, d_r4_3, d_r4_4
    module procedure d_i4_1, d_i4_2, d_i4_3, d_i4_4
    module procedure d_l_1, d_l_2, d_l_3
    module procedure d_r8_0, d_r4_0, d_i4_0, d_l_0
  end interface dump
  contains
  subroutine hdr(label,n,es)
    character(len=*), intent(in) :: label
    integer, intent(in) :: n, es
    character(len=48) :: lab
    lab = label
    write(iu) lab, n, es
  end subroutine hdr
  subroutine d_r8_0(l,a)
    character(len=*), intent(in) :: l
    real(rk8), intent(in) :: a
    call hdr(l,1,storage_size(a)/8) ; write(iu) a
  end subroutine d_r8_0
  subroutine d_r4_0(l,a)
    character(len=*), intent(in) :: l
    real(rk4), intent(in) :: a
    call hdr(l,1,storage_size(a)/8) ; write(iu) a
  end subroutine d_r4_0
  subroutine d_i4_0(l,a)
    character(len=*), intent(in) :: l
    integer(ik4), intent(in) :: a
    call hdr(l,1,storage_size(a)/8) ; write(iu) a
  end subroutine d_i4_0
  subroutine d_l_0(l,a)
    character(len=*), intent(in) :: l
    logical, intent(in) :: a
    call hdr(l,1,4) ; write(iu) merge(1,0,a)
  end subroutine d_l_0
  subroutine d_r8_1(l,a)
    character(len=*), intent(in) :: l
    real(rk8), dimension(:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_r8_1
  subroutine d_r8_2(l,a)
    character(len=*), intent(in) :: l
    real(rk8), dimension(:,:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_r8_2
  subroutine d_r8_3(l,a)
    character(len=*), intent(in) :: l
    real(rk8), dimension(:,:,:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_r8_3
  subroutine d_r8_4(l,a)
    character(len=*), intent(in) :: l
    real(rk8), dimension(:,:,:,:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_r8_4
  subroutine d_r4_1(l,a)
    character(len=*), intent(in) :: l
    real(rk4), dimension(:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_r4_1
  subroutine d_r4_2(l,a)
    character(len=*), intent(in) :: l
    real(rk4), dimension(:,:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_r4_2
  subroutine d_r4_3(l,a)
    character(len=*), intent(in) :: l
    real(rk4), dimension(:,:,:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_r4_3
  subroutine d_r4_4(l,a)
    character(len=*), intent(in) :: l
    real(rk4), dimension(:,:,:,:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_r4_4
  subroutine d_i4_1(l,a)
    character(len=*), intent(in) :: l
    integer(ik4), dimension(:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_i4_1
  subroutine d_i4_2(l,a)
    character(len=*), intent(in) :: l
    integer(ik4), dimension(:,:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_i4_2
  subroutine d_i4_3(l,a)
    character(len=*), intent(in) :: l
    integer(ik4), dimension(:,:,:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_i4_3
  subroutine d_i4_4(l,a)
    character(len=*), intent(in) :: l
    integer(ik4), dimension(:,:,:,:), intent(in) :: a
    call hdr(l,size(a),storage_size(a)/8) ; write(iu) a
  end subroutine d_i4_4
  subroutine d_l_1(l,a)
    character(len=*), intent(in) :: l
    logical, dimension(:), intent(in) :: a
    call hdr(l,size(a),4) ; write(iu) merge(1,0,a)
  end subroutine d_l_1
  subroutine d_l_2(l,a)
    character(len=*), intent(in) :: l
    logical, dimension(:,:), intent(in) :: a
    call hdr(l,size(a),4) ; write(iu) merge(1,0,a)
  end subroutine d_l_2
  subroutine d_l_3(l,a)
    character(len=*), intent(in) :: l
    logical, dimension(:,:,:), intent(in) :: a
    call hdr(l,size(a),4) ; write(iu) merge(1,0,a)
  end subroutine d_l_3
end module mpp_regression_dump

program mpp_regression
  use mpi_f08
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_mppparam
  use mod_regcm_types
  use mod_date
  use mpp_regression_dump
  implicit none
  integer(ik4) :: ierr, icase, nex, j, i, k, n, kk, nn, iv
  character(len=256) :: outdir, arg
  character(len=16) :: tag
  integer(ik4), parameter :: nk = 3, nl = 2

  real(rk8), pointer, contiguous :: r8a2(:,:), r8a3(:,:,:), r8a4(:,:,:,:)
  real(rk4), pointer, contiguous :: r4a2(:,:), r4a3(:,:,:), r4a4(:,:,:,:)
  real(rk8), pointer, contiguous :: r8g2(:,:), r8g3(:,:,:), r8g4(:,:,:,:)
  real(rk4), pointer, contiguous :: r4g2(:,:), r4g3(:,:,:), r4g4(:,:,:,:)
  integer(ik4), pointer, contiguous :: i4a2(:,:), i4a3(:,:,:), i4a4(:,:,:,:)
  integer(ik4), pointer, contiguous :: i4g2(:,:), i4g3(:,:,:), i4g4(:,:,:,:)
  logical, pointer, contiguous :: la2(:,:), la3(:,:,:), la4(:,:,:,:)
  logical, pointer, contiguous :: lg2(:,:), lg3(:,:,:), lg4(:,:,:,:)
  logical, pointer, contiguous :: lmask(:,:,:)
  integer(ik4), pointer, contiguous :: imask(:,:,:)
  real(rk8), pointer, contiguous :: r8v1(:), r8v2(:,:)
  real(rk4), pointer, contiguous :: r4v1(:), r4v2(:,:)
  integer(ik4), pointer, contiguous :: i4v1(:)
  logical, pointer, contiguous :: lv1(:)
  real(rk8), pointer, contiguous :: gm(:,:), sgm(:,:,:)
  real(rk8) :: s8, t8
  real(rk4) :: s4, t4
  integer(ik4) :: si, ti
  logical :: sl, tl
  real(rk8), dimension(3) :: a8
  real(rk4), dimension(3) :: a4
  integer(ik4), dimension(3) :: ai
  logical, dimension(3) :: al
  character(len=8), dimension(2) :: ac
  real(rk8), dimension(2,3) :: m8
  real(rk4), dimension(2,3) :: m4
  real(rk8), allocatable :: g8(:)
  integer(ik4), allocatable :: gi(:)
  real(rk8), allocatable :: v8(:), w8(:)
  real(rk4), allocatable :: v4(:), w4(:)
  type(masked_comm) :: cl
  type(rcm_time_and_date) :: td
  integer(ik4) :: npt

  call mpi_init(ierr)
  call mpi_comm_dup(mpi_comm_world,mycomm,ierr)
  call mpi_comm_rank(mycomm,myid,ierr)
  call mpi_comm_size(mycomm,nproc,ierr)
  call get_command_argument(1,outdir)
  call get_command_argument(2,arg)
  read(arg,*) icase

  jx = 24 ; iy = 18 ; kz = 4
  nsg = 1 ; i_band = 0 ; i_crm = 0
  select case (icase)
    case (1) ; nsg = 2
    case (2) ; i_band = 1
    case (3) ; i_crm = 1 ; i_band = 1 ; nsg = 2
  end select
  njxcpus = 0 ; niycpus = 0
  call derived_dims
  call memory_init
  call set_nproc
  call setup_indexes

  write(tag,'(i0)') myid
  open(newunit=iu,file=trim(outdir)//'/res_'//trim(tag)//'.bin', &
       access='stream',form='unformatted',status='replace')

  call dump('jde',[jde1,jde2,ide1,ide2,jce1,jce2,ice1,ice2, &
                   jci1,jci2,ici1,ici2])

  !======================= halo exchanges =======================
  do nex = 1, 2
    ! ---- real8
    allocate(r8a2(jde1-2:jde2+2,ide1-2:ide2+2))
    allocate(r8a3(jde1-2:jde2+2,ide1-2:ide2+2,nk))
    allocate(r8a4(jde1-2:jde2+2,ide1-2:ide2+2,nk,nl))
    allocate(r4a2(jde1-2:jde2+2,ide1-2:ide2+2))
    allocate(r4a3(jde1-2:jde2+2,ide1-2:ide2+2,nk))
    allocate(r4a4(jde1-2:jde2+2,ide1-2:ide2+2,nk,nl))
    do kk = 1, 6
      call fillx
      select case (kk)
      case (1)
        call exchange(r8a2,nex,jde1,jde2,ide1,ide2)
        call exchange(r8a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange(r8a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
        call exchange(r4a2,nex,jde1,jde2,ide1,ide2)
        call exchange(r4a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange(r4a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
      case (2)
        call exchange_lrbt(r8a2,nex,jde1,jde2,ide1,ide2)
        call exchange_lrbt(r8a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_lrbt(r8a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
        call exchange_lrbt(r4a2,nex,jde1,jde2,ide1,ide2)
        call exchange_lrbt(r4a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_lrbt(r4a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
      case (3)
        call exchange_lr(r8a2,nex,jde1,jde2,ide1,ide2)
        call exchange_lr(r8a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_lr(r8a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
        call exchange_lr(r4a2,nex,jde1,jde2,ide1,ide2)
        call exchange_lr(r4a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_lr(r4a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
      case (4)
        call exchange_bt(r8a2,nex,jde1,jde2,ide1,ide2)
        call exchange_bt(r8a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_bt(r8a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
        call exchange_bt(r4a2,nex,jde1,jde2,ide1,ide2)
        call exchange_bt(r4a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_bt(r4a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
      case (5)
        call exchange_lb(r8a2,nex,jde1,jde2,ide1,ide2)
        call exchange_lb(r8a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_lb(r8a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
        call exchange_lb(r4a2,nex,jde1,jde2,ide1,ide2)
        call exchange_lb(r4a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_lb(r4a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
      case (6)
        call exchange_rt(r8a2,nex,jde1,jde2,ide1,ide2)
        call exchange_rt(r8a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_rt(r8a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
        call exchange_rt(r4a2,nex,jde1,jde2,ide1,ide2)
        call exchange_rt(r4a3,nex,jde1,jde2,ide1,ide2,1,nk)
        call exchange_rt(r4a4,nex,jde1,jde2,ide1,ide2,1,nk,1,nl)
      end select
      call dumpx('exch')
      ! Same on the cross grid sub-range
      call fillx
      select case (kk)
      case (1)
        call exchange(r8a3,nex,jce1,jce2,ice1,ice2,1,nk)
        call exchange(r4a3,nex,jci1,jci2,ici1,ici2,2,nk)
      case (2)
        call exchange_lrbt(r8a3,nex,jce1,jce2,ice1,ice2,1,nk)
        call exchange_lrbt(r4a3,nex,jci1,jci2,ici1,ici2,2,nk)
      case (3)
        call exchange_lr(r8a3,nex,jce1,jce2,ice1,ice2,1,nk)
        call exchange_lr(r4a3,nex,jci1,jci2,ici1,ici2,2,nk)
      case (4)
        call exchange_bt(r8a3,nex,jce1,jce2,ice1,ice2,1,nk)
        call exchange_bt(r4a3,nex,jci1,jci2,ici1,ici2,2,nk)
      case (5)
        call exchange_lb(r8a3,nex,jce1,jce2,ice1,ice2,1,nk)
        call exchange_lb(r4a3,nex,jci1,jci2,ici1,ici2,2,nk)
      case (6)
        call exchange_rt(r8a3,nex,jce1,jce2,ice1,ice2,1,nk)
        call exchange_rt(r4a3,nex,jci1,jci2,ici1,ici2,2,nk)
      end select
      call dump('exch_cross_r8',r8a3)
      call dump('exch_cross_r4',r4a3)
    end do
    deallocate(r8a2,r8a3,r8a4,r4a2,r4a3,r4a4)
  end do

  ! boundary exchanges
  allocate(r8a2(jde1-1:jde2+1,1:kz), r4a2(jde1-1:jde2+1,1:kz))
  r8a2 = -1.0_rk8 ; r4a2 = -1.0_rk4
  do k = 1, kz
    do j = jde1, jde2
      r8a2(j,k) = gv(j,k,k,1) ; r4a2(j,k) = real(gv(j,k,k,1),rk4)
    end do
  end do
  call exchange_bdy_lr(r8a2,1,kz)
  call exchange_bdy_lr(r4a2,2,kz)
  call dump('bdy_lr_r8',r8a2) ; call dump('bdy_lr_r4',r4a2)
  deallocate(r8a2,r4a2)
  allocate(r8a2(ide1-1:ide2+1,1:kz), r4a2(ide1-1:ide2+1,1:kz))
  r8a2 = -1.0_rk8 ; r4a2 = -1.0_rk4
  do k = 1, kz
    do i = ide1, ide2
      r8a2(i,k) = gv(k,i,k,2) ; r4a2(i,k) = real(gv(k,i,k,2),rk4)
    end do
  end do
  call exchange_bdy_bt(r8a2,1,kz)
  call exchange_bdy_bt(r4a2,1,kz-1)
  call dump('bdy_bt_r8',r8a2) ; call dump('bdy_bt_r4',r4a2)
  deallocate(r8a2,r4a2)

  !======================= distribute / collect =======================
  allocate(r8g2(jx,iy),r8g3(jx,iy,nk),r8g4(jx,iy,nk,nl))
  allocate(r4g2(jx,iy),r4g3(jx,iy,nk),r4g4(jx,iy,nk,nl))
  allocate(i4g2(jx,iy),i4g3(jx,iy,nk),i4g4(jx,iy,nk,nl))
  allocate(lg2(jx,iy),lg3(jx,iy,nk),lg4(jx,iy,nk,nl))
  do n = 1, nl
    do k = 1, nk
      do i = 1, iy
        do j = 1, jx
          r8g4(j,i,k,n) = gv(j,i,k,n)
          r4g4(j,i,k,n) = real(gv(j,i,k,n),rk4)
          i4g4(j,i,k,n) = nint(gv(j,i,k,n))
          lg4(j,i,k,n) = mod(j+2*i+k+n,3) == 0
        end do
      end do
    end do
  end do
  r8g3 = r8g4(:,:,:,2) ; r4g3 = r4g4(:,:,:,2)
  i4g3 = i4g4(:,:,:,2) ; lg3 = lg4(:,:,:,2)
  r8g2 = r8g4(:,:,3,1) ; r4g2 = r4g4(:,:,3,1)
  i4g2 = i4g4(:,:,3,1) ; lg2 = lg4(:,:,3,1)

  allocate(r8a2(jde1:jde2,ide1:ide2),r8a3(jde1:jde2,ide1:ide2,nk))
  allocate(r8a4(jde1:jde2,ide1:ide2,nk,nl))
  allocate(r4a2(jde1:jde2,ide1:ide2),r4a3(jde1:jde2,ide1:ide2,nk))
  allocate(r4a4(jde1:jde2,ide1:ide2,nk,nl))
  allocate(i4a2(jde1:jde2,ide1:ide2),i4a3(jde1:jde2,ide1:ide2,nk))
  allocate(i4a4(jde1:jde2,ide1:ide2,nk,nl))
  allocate(la2(jde1:jde2,ide1:ide2),la3(jde1:jde2,ide1:ide2,nk))
  allocate(la4(jde1:jde2,ide1:ide2,nk,nl))
  call clear_local
  call grid_distribute(r8g2,r8a2,jde1,jde2,ide1,ide2)
  call grid_distribute(r8g3,r8a3,jde1,jde2,ide1,ide2,1,nk)
  call grid_distribute(r8g4,r8a4,jde1,jde2,ide1,ide2,1,nk,1,nl)
  call grid_distribute(r4g2,r4a2,jce1,jce2,ice1,ice2)
  call grid_distribute(r4g3,r4a3,jde1,jde2,ide1,ide2,2,nk)
  call grid_distribute(r4g4,r4a4,jde1,jde2,ide1,ide2,1,nk,2,nl)
  call grid_distribute(i4g2,i4a2,jde1,jde2,ide1,ide2)
  call grid_distribute(i4g3,i4a3,jci1,jci2,ici1,ici2,1,nk)
  call grid_distribute(i4g4,i4a4,jde1,jde2,ide1,ide2,1,nk,1,nl)
  call grid_distribute(lg2,la2,jde1,jde2,ide1,ide2)
  call grid_distribute(lg3,la3,jde1,jde2,ide1,ide2,1,nk)
  call grid_distribute(lg4,la4,jde1,jde2,ide1,ide2,1,nk,1,nl)
  call dump_local('dist')

  ! collect back into cleared globals
  call clear_global
  call grid_collect(r8a2,r8g2,jde1,jde2,ide1,ide2)
  call grid_collect(r8a3,r8g3,jde1,jde2,ide1,ide2,1,nk)
  call grid_collect(r8a4,r8g4,jde1,jde2,ide1,ide2,1,nk,1,nl)
  call grid_collect(r4a2,r4g2,jce1,jce2,ice1,ice2)
  call grid_collect(r4a3,r4g3,jde1,jde2,ide1,ide2,2,nk)
  call grid_collect(r4a4,r4g4,jde1,jde2,ide1,ide2,1,nk,2,nl)
  call grid_collect(i4a2,i4g2,jde1,jde2,ide1,ide2)
  call grid_collect(i4a3,i4g3,jci1,jci2,ici1,ici2,1,nk)
  call grid_collect(i4a4,i4g4,jde1,jde2,ide1,ide2,1,nk,1,nl)
  call grid_collect(la2,lg2,jde1,jde2,ide1,ide2)
  call dump_global('coll')
  call clear_global
  call grid_collect(r8a2,r8g3,jde1,jde2,ide1,ide2,2)
  call grid_collect(r8a2,r8g3,jde1,jde2,ide1,ide2)
  call grid_collect(r8a3,r8g2,jde1,jde2,ide1,ide2,3)
  call grid_collect(r8a4,r8g2,jde1,jde2,ide1,ide2,2,2)
  call grid_collect(r4a2,r4g3,jde1,jde2,ide1,ide2,3)
  call grid_collect(r4a3,r4g2,jde1,jde2,ide1,ide2,2)
  call grid_collect(r4a4,r4g2,jde1,jde2,ide1,ide2,1,2)
  call dump('coll_2d3d_r8',r8g3) ; call dump('coll_3d2d_r8',r8g2)
  call dump('coll_2d3d_r4',r4g3) ; call dump('coll_4d2d_r4',r4g2)
  call clear_global
  call grid_collect(r8a3,r8g2,jde1,jde2,ide1,ide2,1)
  call grid_collect(r4a4,r4g2,jde1,jde2,ide1,ide2,3,1)
  call dump('coll_3d2d_r8b',r8g2) ; call dump('coll_4d2d_r4b',r4g2)

  ! grid_fill
  call clear_global
  call grid_fill(r8a2,r8g2)
  call grid_fill(r4a2,r4g2)
  call dump('fill1_r8',r8g2) ; call dump('fill1_r4',r4g2)
  call clear_global
  call grid_fill(r8a2,r8g2,jce1,jce2,ice1,ice2)
  call grid_fill(r4a2,r4g2,jci1,jci2,ici1,ici2)
  call dump('fill2_r8',r8g2) ; call dump('fill2_r4',r4g2)
  deallocate(r8a2,r8a3,r8a4,r4a2,r4a3,r4a4,i4a2,i4a3,i4a4,la2,la3,la4)
  deallocate(r8g2,r8g3,r8g4,r4g2,r4g3,r4g4,i4g2,i4g3,i4g4,lg2,lg3,lg4)

  !======================= subgrid distribute / collect =======================
  allocate(r8g3(nnsg,jx,iy),r8g4(nnsg,jx,iy,nk))
  allocate(r4g3(nnsg,jx,iy),r4g4(nnsg,jx,iy,nk))
  allocate(i4g3(nnsg,jx,iy),i4g4(nnsg,jx,iy,nk))
  allocate(lg3(nnsg,jx,iy))
  do k = 1, nk
    do i = 1, iy
      do j = 1, jx
        do n = 1, nnsg
          r8g4(n,j,i,k) = gv(j,i,k,n)
          r4g4(n,j,i,k) = real(gv(j,i,k,n),rk4)
          i4g4(n,j,i,k) = nint(gv(j,i,k,n))
        end do
      end do
    end do
  end do
  r8g3 = r8g4(:,:,:,2) ; r4g3 = r4g4(:,:,:,2) ; i4g3 = i4g4(:,:,:,3)
  do i = 1, iy
    do j = 1, jx
      do n = 1, nnsg
        lg3(n,j,i) = mod(n+j+i,2) == 0
      end do
    end do
  end do
  allocate(r8a3(nnsg,jde1:jde2,ide1:ide2),r8a4(nnsg,jde1:jde2,ide1:ide2,nk))
  allocate(r4a3(nnsg,jde1:jde2,ide1:ide2),r4a4(nnsg,jde1:jde2,ide1:ide2,nk))
  allocate(i4a3(nnsg,jde1:jde2,ide1:ide2),i4a4(nnsg,jde1:jde2,ide1:ide2,nk))
  allocate(la3(nnsg,jde1:jde2,ide1:ide2))
  allocate(lmask(nnsg,jde1:jde2,ide1:ide2))
  do i = ide1, ide2
    do j = jde1, jde2
      do n = 1, nnsg
        lmask(n,j,i) = mod(n+j*i,3) /= 0
      end do
    end do
  end do
  do kk = 1, 2
    r8a3 = -1.0_rk8 ; r8a4 = -1.0_rk8 ; r4a3 = -1.0_rk4 ; r4a4 = -1.0_rk4
    i4a3 = -1 ; i4a4 = -1 ; la3 = .false.
    if ( kk == 1 ) then
      call subgrid_distribute(r8g3,r8a3,jde1,jde2,ide1,ide2)
      call subgrid_distribute(r8g4,r8a4,jci1,jci2,ici1,ici2,1,nk)
      call subgrid_distribute(r4g3,r4a3,jci1,jci2,ici1,ici2)
      call subgrid_distribute(r4g4,r4a4,jde1,jde2,ide1,ide2,2,nk)
      call subgrid_distribute(i4g3,i4a3,jde1,jde2,ide1,ide2)
      call subgrid_distribute(i4g4,i4a4,jci1,jci2,ici1,ici2,1,nk)
      call subgrid_distribute(lg3,la3,jde1,jde2,ide1,ide2)
    else
      call subgrid_distribute(r8g3,r8a3,jde1,jde2,ide1,ide2,lmask)
      call subgrid_distribute(r8g4,r8a4,jci1,jci2,ici1,ici2,1,nk,lmask)
      call subgrid_distribute(r4g3,r4a3,jci1,jci2,ici1,ici2,lmask)
      call subgrid_distribute(r4g4,r4a4,jde1,jde2,ide1,ide2,2,nk,lmask)
      call subgrid_distribute(i4g3,i4a3,jde1,jde2,ide1,ide2,lmask)
      call subgrid_distribute(i4g4,i4a4,jci1,jci2,ici1,ici2,1,nk,lmask)
      call subgrid_distribute(lg3,la3,jde1,jde2,ide1,ide2,lmask)
    end if
    call dump('sdist_r8_3',r8a3) ; call dump('sdist_r8_4',r8a4)
    call dump('sdist_r4_3',r4a3) ; call dump('sdist_r4_4',r4a4)
    call dump('sdist_i4_3',i4a3) ; call dump('sdist_i4_4',i4a4)
    call dump('sdist_l_3',la3)
  end do
  r8g3 = 0.0_rk8 ; r8g4 = 0.0_rk8 ; r4g3 = 0.0_rk4 ; r4g4 = 0.0_rk4
  i4g3 = 0 ; i4g4 = 0 ; lg3 = .false.
  call subgrid_collect(r8a3,r8g3,jde1,jde2,ide1,ide2)
  call subgrid_collect(r8a4,r8g4,jci1,jci2,ici1,ici2,1,nk)
  call subgrid_collect(r4a3,r4g3,jci1,jci2,ici1,ici2)
  call subgrid_collect(r4a4,r4g4,jde1,jde2,ide1,ide2,2,nk)
  call subgrid_collect(i4a3,i4g3,jde1,jde2,ide1,ide2)
  call subgrid_collect(i4a4,i4g4,jci1,jci2,ici1,ici2,1,nk)
  call subgrid_collect(la3,lg3,jde1,jde2,ide1,ide2)
  call dump('scoll_r8_3',r8g3) ; call dump('scoll_r8_4',r8g4)
  call dump('scoll_r4_3',r4g3) ; call dump('scoll_r4_4',r4g4)
  call dump('scoll_i4_3',i4g3) ; call dump('scoll_i4_4',i4g4)
  call dump('scoll_l_3',lg3)
  deallocate(r8a3,r8a4,r4a3,r4a4,i4a3,i4a4,la3,lmask)
  deallocate(r8g3,r8g4,r4g3,r4g4,i4g3,i4g4,lg3)

  !======================= reorder =======================
  allocate(r8a3(nnsg,jci1:jci2,ici1:ici2),r8a4(nnsg,jci1:jci2,ici1:ici2,nk))
  allocate(r4a3(nnsg,jci1:jci2,ici1:ici2),r4a4(nnsg,jci1:jci2,ici1:ici2,nk))
  allocate(imask(nnsg,jci1:jci2,ici1:ici2), la3(nnsg,jci1:jci2,ici1:ici2))
  do k = 1, nk
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          r8a4(n,j,i,k) = gv(j,i,k,n)
          r4a4(n,j,i,k) = real(gv(j,i,k,n),rk4)
          imask(n,j,i) = mod(n+j+2*i,3)
          la3(n,j,i) = mod(n*j+i,2) == 0
        end do
      end do
    end do
  end do
  r8a3 = r8a4(:,:,:,1) ; r4a3 = r4a4(:,:,:,2)
  allocate(r8g2((jci1-1)*nsg+1:jci2*nsg,(ici1-1)*nsg+1:ici2*nsg))
  allocate(r8g3((jci1-1)*nsg+1:jci2*nsg,(ici1-1)*nsg+1:ici2*nsg,nk))
  allocate(r4g2((jci1-1)*nsg+1:jci2*nsg,(ici1-1)*nsg+1:ici2*nsg))
  allocate(r4g3((jci1-1)*nsg+1:jci2*nsg,(ici1-1)*nsg+1:ici2*nsg,nk))
  allocate(lg2((jci1-1)*nsg+1:jci2*nsg,(ici1-1)*nsg+1:ici2*nsg))
  do kk = 1, 2
    r8g2 = 1.0_rk8 ; r8g3 = 1.0_rk8 ; r4g2 = 1.0_rk4 ; r4g3 = 1.0_rk4
    if ( kk == 1 ) then
      call reorder_subgrid(r8a3,r8g2)
      call reorder_subgrid(r8a3,r8g3,2)
      call reorder_subgrid(r8a3,r8g3)
      call reorder_subgrid(r4a3,r4g2)
      call reorder_subgrid(r4a3,r4g3,3)
    else
      call reorder_subgrid(r8a3,r8g2,mask=imask)
      call reorder_subgrid(r8a3,r8g3,2,imask)
      call reorder_subgrid(r4a3,r4g2,mask=imask)
      call reorder_subgrid(r4a3,r4g3,3,imask)
    end if
    call dump('reo_2d_r8',r8g2) ; call dump('reo_2d3d_r8',r8g3)
    call dump('reo_2d_r4',r4g2) ; call dump('reo_2d3d_r4',r4g3)
    r8g2 = 1.0_rk8 ; r4g2 = 1.0_rk4
    if ( kk == 1 ) then
      call reorder_subgrid(r8a4,r8g2,3)
      call reorder_subgrid(r4a4,r4g2,2)
      call reorder_subgrid(r8a4,r8g3)
      call reorder_subgrid(r4a4,r4g3)
    else
      call reorder_subgrid(r8a4,r8g2,3,imask)
      call reorder_subgrid(r4a4,r4g2,2,imask)
      call reorder_subgrid(r8a4,r8g3,imask)
      call reorder_subgrid(r4a4,r4g3,imask)
    end if
    call dump('reo_3d_r8',r8g2) ; call dump('reo_4d_r8',r8g3)
    call dump('reo_3d_r4',r4g2) ; call dump('reo_4d_r4',r4g3)
    r8g2 = 1.0_rk8 ; r8g3 = 1.0_rk8 ; r4g2 = 1.0_rk4 ; r4g3 = 1.0_rk4
    if ( kk == 1 ) then
      call reorder_add_subgrid(r8a3,r8g2)
      call reorder_add_subgrid(r8a3,r8g3,2)
      call reorder_add_subgrid(r8a4,r8g2,3)
      call reorder_add_subgrid(r4a3,r4g2)
      call reorder_add_subgrid(r4a3,r4g3)
      call reorder_add_subgrid(r4a4,r4g2,1)
    else
      call reorder_add_subgrid(r8a3,r8g2,mask=imask)
      call reorder_add_subgrid(r8a3,r8g3,2,imask)
      call reorder_add_subgrid(r8a4,r8g2,3,imask)
      call reorder_add_subgrid(r4a3,r4g2,mask=imask)
      call reorder_add_subgrid(r4a3,r4g3,1,imask)
      call reorder_add_subgrid(r4a4,r4g2,1,imask)
    end if
    call dump('reoadd_r8_2',r8g2) ; call dump('reoadd_r8_3',r8g3)
    call dump('reoadd_r4_2',r4g2) ; call dump('reoadd_r4_3',r4g3)
  end do
  lg2 = .false.
  call reorder_subgrid(la3,lg2)
  call dump('reo_l',lg2)
  ! input_reorder: from a global-subgrid indexed 2D field
  r8a3 = -1.0_rk8 ; r4a3 = -1.0_rk4
  allocate(gm(jxsg,iysg))
  do i = 1, iysg
    do j = 1, jxsg
      gm(j,i) = gv(j,i,1,1)
    end do
  end do
  call input_reorder(gm,r8a3,jci1,jci2,ici1,ici2)
  call input_reorder(real(gm,rk4),r4a3,jci1,jci2,ici1,ici2)
  call dump('inreo_r8',r8a3) ; call dump('inreo_r4',r4a3)
  deallocate(gm)
  deallocate(r8a3,r8a4,r4a3,r4a4,imask,la3,r8g2,r8g3,r4g2,r4g3,lg2)

  !======================= cartesian <-> linear =======================
  allocate(gm(jde1:jde2,ide1:ide2),sgm(nnsg,jde1:jde2,ide1:ide2))
  do i = ide1, ide2
    do j = jde1, jde2
      gm(j,i) = merge(1.0_rk8,0.0_rk8,mod(j+i,3) /= 0)
      do n = 1, nnsg
        sgm(n,j,i) = merge(1.0_rk8,0.0_rk8,mod(j+i+n,4) /= 0)
      end do
    end do
  end do
  ! Without subgrid the model has sgmask == gmask (clset relies on it)
  if ( nsg == 1 ) sgm(1,:,:) = gm
  do kk = 1, 2
    call cl_setup(cl,gm,sgm,kk==2)
    call dump('cl_npg',cl%linear_npoint_g)
    call dump('cl_npsg',cl%linear_npoint_sg)
    call dump('cl_sgid',cl%sgmask_id)
    call dump('cl_gid',cl%gmask_id)
    npt = cl%linear_npoint_sg(myid+1)
    allocate(r8v1(npt),r4v1(npt),i4v1(npt),lv1(npt))
    allocate(r8v2(npt,nk),r4v2(npt,nk))
    allocate(r8a2(jci1:jci2,ici1:ici2),r4a2(jci1:jci2,ici1:ici2))
    allocate(i4a2(jci1:jci2,ici1:ici2),la2(jci1:jci2,ici1:ici2))
    allocate(r8a3(nnsg,jci1:jci2,ici1:ici2),r4a3(nnsg,jci1:jci2,ici1:ici2))
    allocate(i4a3(nnsg,jci1:jci2,ici1:ici2),la3(nnsg,jci1:jci2,ici1:ici2))
    allocate(r8a4(nnsg,jci1:jci2,ici1:ici2,nk))
    allocate(r4a4(nnsg,jci1:jci2,ici1:ici2,nk))
    do k = 1, nk
      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1, nnsg
            r8a4(n,j,i,k) = gv(j,i,k,n)
            r4a4(n,j,i,k) = real(gv(j,i,k,n),rk4)
          end do
        end do
      end do
    end do
    r8a3 = r8a4(:,:,:,1) ; r4a3 = r4a4(:,:,:,2)
    i4a3 = nint(r8a3) ; la3 = mod(i4a3,3) == 1
    r8a2 = r8a3(1,:,:) ; r4a2 = r4a3(1,:,:)
    i4a2 = i4a3(1,:,:) ; la2 = la3(1,:,:)
    ! cartesian -> linear
    call c2l_gs(cl,r8a2,r8v1) ; call dump('c2lgs_r8',r8v1)
    call c2l_gs(cl,r4a2,r4v1) ; call dump('c2lgs_r4',r4v1)
    call c2l_gs(cl,i4a2,i4v1) ; call dump('c2lgs_i4',i4v1)
    call c2l_gs(cl,la2,lv1) ; call dump('c2lgs_l',lv1)
    call c2l_ss(cl,r8a3,r8v1) ; call dump('c2lss_r8',r8v1)
    call c2l_ss(cl,r4a3,r4v1) ; call dump('c2lss_r4',r4v1)
    call c2l_ss(cl,i4a3,i4v1) ; call dump('c2lss_i4',i4v1)
    call c2l_ss(cl,la3,lv1) ; call dump('c2lss_l',lv1)
    call c2l_ss(cl,r8a4,r8v2) ; call dump('c2lss4_r8',r8v2)
    call c2l_ss(cl,r4a4,r4v2) ; call dump('c2lss4_r4',r4v2)
    ! linear -> cartesian
    r8a3 = -9.0_rk8 ; r4a3 = -9.0_rk4 ; i4a3 = -9 ; la3 = .false.
    r8a4 = -9.0_rk8 ; r4a4 = -9.0_rk4
    call l2c_ss(cl,r8v1,r8a3) ; call dump('l2css_r8',r8a3)
    call l2c_ss(cl,r4v1,r4a3) ; call dump('l2css_r4',r4a3)
    call l2c_ss(cl,i4v1,i4a3) ; call dump('l2css_i4',i4a3)
    call l2c_ss(cl,lv1,la3) ; call dump('l2css_l',la3)
    call l2c_ss(cl,r8v2,r8a4) ; call dump('l2css4_r8',r8a4)
    call l2c_ss(cl,r4v2,r4a4) ; call dump('l2css4_r4',r4a4)
    ! global variants
    do k = 1, nk
      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1, nnsg
            r8a4(n,j,i,k) = gv(j,i,k,n) + 0.5_rk8
            r4a4(n,j,i,k) = real(gv(j,i,k,n),rk4) + 0.25_rk4
          end do
        end do
      end do
    end do
    r8a3 = r8a4(:,:,:,2) ; r4a3 = r4a4(:,:,:,1)
    i4a3 = nint(r8a3) ; la3 = mod(i4a3,4) == 1
    r8a2 = r8a3(1,:,:) ; r4a2 = r4a3(1,:,:)
    i4a2 = i4a3(1,:,:) ; la2 = la3(1,:,:)
    r8v1 = 0.0_rk8 ; r4v1 = 0.0_rk4 ; i4v1 = 0 ; lv1 = .false.
    call glb_c2l_gs(cl,r8a2,r8v1) ; call dump('gc2lgs_r8',r8v1)
    call glb_c2l_gs(cl,r4a2,r4v1) ; call dump('gc2lgs_r4',r4v1)
    call glb_c2l_gs(cl,i4a2,i4v1) ; call dump('gc2lgs_i4',i4v1)
    call glb_c2l_gs(cl,la2,lv1) ; call dump('gc2lgs_l',lv1)
    call glb_c2l_gs(cl,r4a2,r8v1) ; call dump('gc2lgs_r48',r8v1)
    call glb_c2l_ss(cl,r8a3,r8v1) ; call dump('gc2lss_r8',r8v1)
    call glb_c2l_ss(cl,r4a3,r4v1) ; call dump('gc2lss_r4',r4v1)
    call glb_c2l_ss(cl,i4a3,i4v1) ; call dump('gc2lss_i4',i4v1)
    call glb_c2l_ss(cl,la3,lv1) ; call dump('gc2lss_l',lv1)
    call glb_c2l_ss(cl,r8a4,r8v2) ; call dump('gc2lss4_r8',r8v2)
    call glb_c2l_ss(cl,r4a4,r4v2) ; call dump('gc2lss4_r4',r4v2)
    call glb_c2l_ss(cl,r4a3,r8v1) ; call dump('gc2lss_r48',r8v1)
    call glb_c2l_ss(cl,r4a4,r8v2) ; call dump('gc2lss4_r48',r8v2)
    r8a3 = -9.0_rk8 ; r4a3 = -9.0_rk4 ; i4a3 = -9 ; la3 = .false.
    r8a4 = -9.0_rk8 ; r4a4 = -9.0_rk4
    call glb_l2c_ss(cl,r8v1,r8a3) ; call dump('gl2css_r8',r8a3)
    call glb_l2c_ss(cl,r4v1,r4a3) ; call dump('gl2css_r4',r4a3)
    call glb_l2c_ss(cl,i4v1,i4a3) ; call dump('gl2css_i4',i4a3)
    call glb_l2c_ss(cl,lv1,la3) ; call dump('gl2css_l',la3)
    call glb_l2c_ss(cl,r8v2,r8a4) ; call dump('gl2css4_r8',r8a4)
    call glb_l2c_ss(cl,r4v2,r4a4) ; call dump('gl2css4_r4',r4a4)
    r4a3 = -9.0_rk4 ; r4a4 = -9.0_rk4
    call glb_l2c_ss(cl,r8v1,r4a3) ; call dump('gl2css_r84',r4a3)
    call glb_l2c_ss(cl,r8v2,r4a4) ; call dump('gl2css4_r84',r4a4)
    deallocate(r8v1,r4v1,i4v1,lv1,r8v2,r4v2,r8a2,r4a2,i4a2,la2)
    deallocate(r8a3,r4a3,i4a3,la3,r8a4,r4a4)
    call cl_dispose(cl)
  end do
  deallocate(gm,sgm)

  !======================= collectives =======================
  s8 = 0.0_rk8 ; s4 = 0.0_rk4 ; si = 0 ; sl = .false.
  a8 = 0.0_rk8 ; a4 = 0.0_rk4 ; ai = 0 ; al = .false.
  ac = '' ; m8 = 0.0_rk8 ; m4 = 0.0_rk4
  if ( myid == 0 ) then
    s8 = 1.25_rk8 ; s4 = 2.5_rk4 ; si = 7 ; sl = .true.
    a8 = [1.0_rk8,2.0_rk8,3.0_rk8] ; a4 = [4.0_rk4,5.0_rk4,6.0_rk4]
    ai = [7,8,9] ; al = [.true.,.false.,.true.]
    ac = ['alpha   ','beta    ']
    m8 = reshape([1,2,3,4,5,6],[2,3]) ; m4 = reshape([6,5,4,3,2,1],[2,3])
    td = 1950010106
  end if
  call bcast(s8) ; call bcast(s4) ; call bcast(si) ; call bcast(sl)
  call bcast(a8) ; call bcast(a4) ; call bcast(ai) ; call bcast(al)
  call bcast(ac,8) ; call bcast(m8) ; call bcast(m4) ; call bcast(td)
  call dump('bc_s8',s8) ; call dump('bc_s4',s4) ; call dump('bc_si',si)
  call dump('bc_sl',sl) ; call dump('bc_a8',a8) ; call dump('bc_a4',a4)
  call dump('bc_ai',ai) ; call dump('bc_al',al)
  call dump('bc_m8',m8) ; call dump('bc_m4',m4)
  call dump('bc_ac',ichar(ac(2)(1:1)))
  call dump('bc_td',int(toint10(td),ik4))
  s8 = 1.0_rk8 + myid ; s4 = 2.0_rk4 + myid ; si = 3 + myid
  call sumall(s8,t8) ; call dump('sum_r8',t8)
  call sumall(s4,t4) ; call dump('sum_r4',t4)
  call sumall(si,ti) ; call dump('sum_i4',ti)
  call maxall(s8,t8) ; call dump('max_r8',t8)
  call maxall(s4,t4) ; call dump('max_r4',t4)
  call maxall(si,ti) ; call dump('max_i4',ti)
  call minall(s8,t8) ; call dump('min_r8',t8)
  call minall(s4,t4) ; call dump('min_r4',t4)
  call minall(si,ti) ; call dump('min_i4',ti)
  call meanall(s8,t8) ; call dump('mean_r8',t8)
  call meanall(s4,t4) ; call dump('mean_r4',t4)
  allocate(v8(3),w8(3),v4(3),w4(3))
  v8 = [1.0_rk8,2.0_rk8,3.0_rk8]*(myid+1) ; v4 = real(v8,rk4)
  call meanall(v8,w8,3) ; call dump('mean1d_r8',w8)
  call meanall(v4,w4,3) ; call dump('mean1d_r4',w4)
  call trueforall(mod(myid,2)==0,tl) ; call dump('trueforall',tl)
  call trueforall(.true.,tl) ; call dump('trueforall2',tl)
  allocate(g8(nproc),gi(nproc))
  g8 = 0.0_rk8 ; gi = 0
  call gather_r(g8,s8) ; call dump('gather_r',g8)
  call gather_i(gi,si) ; call dump('gather_i',gi)
  call allgather_r(g8,s8+1.0_rk8) ; call dump('allgather_r',g8)
  call allgather_i(gi,si+1) ; call dump('allgather_i',gi)
  allocate(r8a3(2,3,kz),r8g3(2,3,kz),r4a3(2,3,kz),r4g3(2,3,kz))
  r8a3 = real(myid+1,rk8) ; r4a3 = real(myid+2,rk4)
  call column_reduce(r8a3,r8g3) ; call dump('colred_r8',r8g3)
  call column_reduce(r4a3,r4g3) ; call dump('colred_r4',r4g3)
  call row_reduce(r8a3,r8g3) ; call dump('rowred_r8',r8g3)
  call row_reduce(r4a3,r4g3) ; call dump('rowred_r4',r4g3)

  close(iu)
  call mpi_barrier(mycomm,ierr)
  if ( myid == 0 ) write(*,*) 'mpp_regression case ', icase, ' done on ', nproc
  call mpi_finalize(ierr)

  contains

  pure real(rk8) function gv(j,i,k,n)
    integer(ik4), intent(in) :: j, i, k, n
    gv = real(j,rk8) + 100.0_rk8*i + 10000.0_rk8*k + 100000.0_rk8*n
  end function gv

  subroutine fillx
    r8a2 = -1.0_rk8 ; r8a3 = -1.0_rk8 ; r8a4 = -1.0_rk8
    r4a2 = -1.0_rk4 ; r4a3 = -1.0_rk4 ; r4a4 = -1.0_rk4
    do n = 1, nl
      do k = 1, nk
        do i = ide1, ide2
          do j = jde1, jde2
            r8a4(j,i,k,n) = gv(j,i,k,n)
            r4a4(j,i,k,n) = real(gv(j,i,k,n),rk4)
          end do
        end do
      end do
    end do
    do k = 1, nk
      do i = ide1, ide2
        do j = jde1, jde2
          r8a3(j,i,k) = gv(j,i,k,3)
          r4a3(j,i,k) = real(gv(j,i,k,3),rk4)
        end do
      end do
    end do
    do i = ide1, ide2
      do j = jde1, jde2
        r8a2(j,i) = gv(j,i,5,0)
        r4a2(j,i) = real(gv(j,i,5,0),rk4)
      end do
    end do
  end subroutine fillx

  subroutine dumpx(l)
    character(len=*), intent(in) :: l
    call dump(l//'_r8_2',r8a2) ; call dump(l//'_r8_3',r8a3)
    call dump(l//'_r8_4',r8a4) ; call dump(l//'_r4_2',r4a2)
    call dump(l//'_r4_3',r4a3) ; call dump(l//'_r4_4',r4a4)
  end subroutine dumpx

  subroutine clear_local
    r8a2 = -1.0_rk8 ; r8a3 = -1.0_rk8 ; r8a4 = -1.0_rk8
    r4a2 = -1.0_rk4 ; r4a3 = -1.0_rk4 ; r4a4 = -1.0_rk4
    i4a2 = -1 ; i4a3 = -1 ; i4a4 = -1
    la2 = .false. ; la3 = .false. ; la4 = .false.
  end subroutine clear_local

  subroutine clear_global
    r8g2 = 0.0_rk8 ; r8g3 = 0.0_rk8 ; r8g4 = 0.0_rk8
    r4g2 = 0.0_rk4 ; r4g3 = 0.0_rk4 ; r4g4 = 0.0_rk4
    i4g2 = 0 ; i4g3 = 0 ; i4g4 = 0
    lg2 = .false. ; lg3 = .false. ; lg4 = .false.
  end subroutine clear_global

  subroutine dump_local(l)
    character(len=*), intent(in) :: l
    call dump(l//'_r8_2',r8a2) ; call dump(l//'_r8_3',r8a3)
    call dump(l//'_r8_4',r8a4) ; call dump(l//'_r4_2',r4a2)
    call dump(l//'_r4_3',r4a3) ; call dump(l//'_r4_4',r4a4)
    call dump(l//'_i4_2',i4a2) ; call dump(l//'_i4_3',i4a3)
    call dump(l//'_i4_4',i4a4) ; call dump(l//'_l_2',la2)
    call dump(l//'_l_3',la3)
    call dump(l//'_l_4',reshape(la4,[size(la4)]))
  end subroutine dump_local

  subroutine dump_global(l)
    character(len=*), intent(in) :: l
    call dump(l//'_r8_2',r8g2) ; call dump(l//'_r8_3',r8g3)
    call dump(l//'_r8_4',r8g4) ; call dump(l//'_r4_2',r4g2)
    call dump(l//'_r4_3',r4g3) ; call dump(l//'_r4_4',r4g4)
    call dump(l//'_i4_2',i4g2) ; call dump(l//'_i4_3',i4g3)
    call dump(l//'_i4_4',i4g4) ; call dump(l//'_l_2',lg2)
  end subroutine dump_global

  subroutine derived_dims
    iym1 = iy - 1 ; iym2 = iy - 2 ; iym3 = iy - 3
    jxm1 = jx - 1 ; jxm2 = jx - 2 ; jxm3 = jx - 3
    kzm1 = kz - 1 ; kzm2 = kz - 2 ; kzp1 = kz + 1
    kzp2 = kz + 2 ; kzp3 = kz + 3 ; kzp4 = kz + 4
    iysg = iy * nsg ; jxsg = jx * nsg
    nnsg = nsg*nsg
    jdot1 = 1 ; jdot2 = jx ; jcross1 = 1
    if ( i_band == 1 ) then
      jcross2 = jx ; jout1 = 1 ; jout2 = jx
      joutsg1 = 1 ; joutsg2 = jx*nsg
    else
      jcross2 = jxm1 ; jout1 = 2 ; jout2 = jxm2
      joutsg1 = nsg+1 ; joutsg2 = jxm2*nsg
    end if
    idot1 = 1 ; idot2 = iy ; icross1 = 1
    if ( i_crm == 1 ) then
      icross2 = iy ; iout1 = 1 ; iout2 = iy
      ioutsg1 = 1 ; ioutsg2 = iy*nsg
    else
      icross2 = iym1 ; iout1 = 2 ; iout2 = iym2
      ioutsg1 = nsg+1 ; ioutsg2 = iym2*nsg
    end if
  end subroutine derived_dims

  ! Copy of mod_atm_interface:setup_model_indexes (index part only)
  subroutine setup_indexes
    jde1 = global_dot_jstart ; jdi1 = global_dot_jstart
    jde2 = global_dot_jend   ; jdi2 = global_dot_jend
    ide1 = global_dot_istart ; idi1 = global_dot_istart
    ide2 = global_dot_iend   ; idi2 = global_dot_iend
    if ( ma%has_bdyleft ) jdi1 = jde1 + 1
    if ( ma%has_bdyright ) jdi2 = jde2 - 1
    if ( ma%has_bdybottom ) idi1 = ide1 + 1
    if ( ma%has_bdytop ) idi2 = ide2 - 1
    jce1 = global_cross_jstart ; jci1 = global_cross_jstart
    jce2 = global_cross_jend   ; jci2 = global_cross_jend
    ice1 = global_cross_istart ; ici1 = global_cross_istart
    ice2 = global_cross_iend   ; ici2 = global_cross_iend
    if ( ma%has_bdyleft ) jci1 = jce1 + 1
    if ( ma%has_bdyright ) jci2 = jce2 - 1
    if ( ma%has_bdybottom ) ici1 = ice1 + 1
    if ( ma%has_bdytop ) ici2 = ice2 - 1
  end subroutine setup_indexes

end program mpp_regression

! The library calls myabort (Main/abort.F90) on fatal errors.
subroutine myabort
  use mpi_f08
  implicit none
  integer :: ierr
  call mpi_abort(mpi_comm_world,1,ierr)
end subroutine myabort

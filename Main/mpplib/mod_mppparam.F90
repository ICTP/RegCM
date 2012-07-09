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
!
module mod_mppparam

#ifndef IBM
  use mpi
#endif
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_memutil
  use mod_date
  use mod_stdio
  use netcdf

  private

#ifdef IBM
  include 'mpif.h'
#endif

  public :: set_nproc , broadcast_params , date_bcast

  integer :: cartesian_communicator

  type model_area
    logical :: bandflag
    logical :: has_bdyleft , has_bdyright , has_bdytop , has_bdybottom
    integer , dimension(2) :: location
    integer :: left , right , top , bottom
    integer :: topleft , topright , bottomleft , bottomright
    integer :: ibt1 , ibt2 , ibb1 , ibb2
    integer :: jbl1 , jbl2 , jbr1 , jbr2
  end type model_area

  type grid_nc_var2d
    character(len=64) :: varname
    integer :: irec = -1
    integer :: ncid = -1
    integer :: varid = 1
    integer :: nx = 0
    integer :: ny = 0
    integer :: mynx1 = 0
    integer :: mynx2 = 0
    integer :: myny1 = 0
    integer :: myny2 = 0
    real(dp) , pointer , dimension(:,:) :: val => null()
    real(dp) , pointer , dimension(:,:) :: iobuf => null()
  end type grid_nc_var2d

  type grid_nc_var3d
    character(len=64) :: varname
    integer :: irec = -1
    integer :: ncid = -1
    integer :: varid = 1
    integer :: nx = 0
    integer :: ny = 0
    integer :: mynx1 = 0
    integer :: mynx2 = 0
    integer :: myny1 = 0
    integer :: myny2 = 0
    integer :: nz = 0
    real(dp) , pointer , dimension(:,:,:) :: val => null()
    real(dp) , pointer , dimension(:,:,:) :: iobuf => null()
  end type grid_nc_var3d

  type grid_nc_var4d
    character(len=64) :: varname
    integer :: irec = -1
    integer :: ncid = -1
    integer :: varid = 1
    integer :: nx = 0
    integer :: ny = 0
    integer :: mynx1 = 0
    integer :: mynx2 = 0
    integer :: myny1 = 0
    integer :: myny2 = 0
    integer :: nz = 0
    integer :: nl = 0
    real(dp) , pointer , dimension(:,:,:,:) :: val => null()
    real(dp) , pointer , dimension(:,:,:,:) :: iobuf => null()
  end type grid_nc_var4d

  public :: grid_nc_var2d , grid_nc_var3d , grid_nc_var4d

  interface grid_nc_create
    module procedure grid_nc_create_var2d , &
                     grid_nc_create_var3d , &
                     grid_nc_create_var4d
  end interface grid_nc_create

  interface grid_nc_write
    module procedure grid_nc_write_var2d , &
                     grid_nc_write_var3d , &
                     grid_nc_write_var4d
  end interface grid_nc_write

  interface grid_nc_destroy
    module procedure grid_nc_destroy_var2d , &
                     grid_nc_destroy_var3d , &
                     grid_nc_destroy_var4d
  end interface grid_nc_destroy

  public :: grid_nc_create , grid_nc_write , grid_nc_destroy

  interface grid_distribute
    module procedure real8_2d_distribute ,   &
                     real8_3d_distribute ,   &
                     real8_4d_distribute ,   &
                     real4_2d_distribute ,   &
                     real4_3d_distribute ,   &
                     real4_4d_distribute ,   &
                     integer_2d_distribute , &
                     integer_3d_distribute , &
                     integer_4d_distribute
  end interface grid_distribute

  interface subgrid_distribute
    module procedure real8_2d_sub_distribute ,   &
                     real8_3d_sub_distribute ,   &
                     real4_2d_sub_distribute ,   &
                     real4_3d_sub_distribute ,   &
                     integer_2d_sub_distribute , &
                     integer_3d_sub_distribute
  end interface subgrid_distribute

  interface grid_collect
    module procedure real8_2d_collect ,   &
                     real8_3d_collect ,   &
                     real8_4d_collect ,   &
                     real4_2d_collect ,   &
                     real4_3d_collect ,   &
                     real4_4d_collect ,   &
                     integer_2d_collect , &
                     integer_3d_collect , &
                     integer_4d_collect
  end interface grid_collect

  interface subgrid_collect
    module procedure real8_2d_sub_collect ,   &
                     real8_3d_sub_collect ,   &
                     real4_2d_sub_collect ,   &
                     real4_3d_sub_collect ,   &
                     integer_2d_sub_collect , &
                     integer_3d_sub_collect
  end interface subgrid_collect

  interface exchange
    module procedure real8_2d_exchange , &
                     real8_3d_exchange , &
                     real8_4d_exchange
  end interface exchange

  interface exchange_lb
    module procedure real8_2d_exchange_left_bottom , &
                     real8_3d_exchange_left_bottom , &
                     real8_4d_exchange_left_bottom
  end interface exchange_lb

  interface exchange_rt
    module procedure real8_2d_exchange_right_top , &
                     real8_3d_exchange_right_top , &
                     real8_4d_exchange_right_top
  end interface exchange_rt

  interface exchange_lr
    module procedure real8_bdy_exchange_left_right
  end interface exchange_lr

  interface exchange_tb
    module procedure real8_bdy_exchange_top_bottom
  end interface exchange_tb

  interface grid_fill
    module procedure real8_2d_grid_fill
  end interface grid_fill

#ifdef DEBUG
  type(grid_nc_var4d) , public :: qqxp
#endif

  public :: model_area
  type(model_area) , public :: ma
!
  real(dp) , pointer , dimension(:) :: r8vector1
  real(dp) , pointer , dimension(:) :: r8vector2
  real(sp) , pointer , dimension(:) :: r4vector1
  real(sp) , pointer , dimension(:) :: r4vector2
  integer , pointer , dimension(:) :: i4vector1
  integer , pointer , dimension(:) :: i4vector2
  integer , dimension(4) :: window
  integer :: mpierr
!
  public :: exchange , exchange_lb , exchange_rt , exchange_lr , exchange_tb
  public :: grid_distribute , grid_collect , grid_fill
  public :: subgrid_distribute , subgrid_collect
  public :: uvcross2dot , psc2psd
!
  contains
!
  subroutine set_nproc(ncpu)
    implicit none
    integer , intent(in) :: ncpu
    integer , dimension(2) :: cpus_per_dim
    logical , dimension(2) :: dim_period
    integer , dimension(2) :: isearch
    integer :: imaxcpus
    real(dp) :: dimfac
    data dim_period /.false.,.false./

    nproc = ncpu

    if ( i_band == 1 ) dim_period(1) = .true.

    if ( nproc == 1 ) then
      cpus_per_dim(1) = 1
      cpus_per_dim(2) = 1
      jxp =  jx
      iyp =  iy
      cartesian_communicator = mycomm
      ma%location(1) = 0
      ma%location(2) = 0
      ma%top         = mpi_proc_null
      ma%bottom      = mpi_proc_null
      ma%left        = mpi_proc_null
      ma%right       = mpi_proc_null
      ma%bottomleft  = mpi_proc_null
      ma%bottomright = mpi_proc_null
      ma%topleft     = mpi_proc_null
      ma%topright    = mpi_proc_null
    else
      if ( nproc < 4 ) then
        cpus_per_dim(1) = ncpu
        cpus_per_dim(2) = 1
      else if ( nproc >= 4 ) then
        if ( mod(nproc,2) /= 0 ) then
          write(stderr,*) 'Work does not evenly divide.'
          write(stderr,*) 'Required number of CPUS must be even.'
          call fatal(__FILE__,__LINE__,'CPU/WORK mismatch')
        end if
        if ( iy > jx ) then
          dimfac = dble(iy)/dble(iy+jx)
          cpus_per_dim(2) = (int(dble(ncpu)*dimfac)/2)*2
          cpus_per_dim(1) = ncpu / cpus_per_dim(2)
        else
          dimfac = dble(jx)/dble(iy+jx)
          cpus_per_dim(1) = (int(dble(ncpu)*dimfac)/2)*2
          cpus_per_dim(2) = ncpu / cpus_per_dim(1)
        end if
        imaxcpus = cpus_per_dim(1)*cpus_per_dim(2)
        if ( mod(ncpu,imaxcpus) /= 0 ) then
          write(stderr,*) 'Work does not evenly divide.'
          write(stderr,*) 'Suggested minimum number of CPUS : ', imaxcpus
          call fatal(__FILE__,__LINE__,'CPU/WORK mismatch')
        end if
      end if
      call mpi_cart_create(mycomm,2,cpus_per_dim,dim_period,.true., &
                           cartesian_communicator,mpierr)
      call mpi_comm_rank(cartesian_communicator,myid,mpierr)
      call mpi_cart_coords(cartesian_communicator,myid,2,ma%location,mpierr)

      ! Set coordinates in the grid for the other processors
      isearch(1) = ma%location(1)
      isearch(2) = ma%location(2)+1
      if ( isearch(2) < cpus_per_dim(2) ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%top,mpierr)
      else
        ma%top = mpi_proc_null
      end if
      isearch(1) = ma%location(1)
      isearch(2) = ma%location(2)-1
      if ( isearch(2) >= 0 ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%bottom,mpierr)
      else
        ma%bottom = mpi_proc_null
      end if
      isearch(1) = ma%location(1)-1
      isearch(2) = ma%location(2)
      if ( isearch(1) >= 0 ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%left,mpierr)
      else
        ma%left = mpi_proc_null
      end if
      isearch(1) = ma%location(1)+1
      isearch(2) = ma%location(2)
      if ( isearch(1) < cpus_per_dim(1) .or. i_band == 1 ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%right,mpierr)
      else
        ma%right = mpi_proc_null
      end if
      isearch(1) = ma%location(1)+1
      isearch(2) = ma%location(2)+1
      if ( ( isearch(1) < cpus_per_dim(1) .or. i_band == 1 ) .and. &
             isearch(2) < cpus_per_dim(2) ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%topright,mpierr)
      else
        ma%topright = mpi_proc_null
      end if
      isearch(1) = ma%location(1)-1
      isearch(2) = ma%location(2)+1
      if ( ( isearch(1) >= 0 .or. i_band == 1 ) .and. &
             isearch(2) < cpus_per_dim(2) ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%topleft,mpierr)
      else
        ma%topleft = mpi_proc_null
      end if
      isearch(1) = ma%location(1)+1
      isearch(2) = ma%location(2)-1
      if ( ( isearch(1) < cpus_per_dim(1) .or. i_band == 1 ) .and. &
             isearch(2) >= 0 ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%bottomright,mpierr)
      else
        ma%bottomright = mpi_proc_null
      end if
      isearch(1) = ma%location(1)-1
      isearch(2) = ma%location(2)-1
      if ( ( isearch(1) >= 0 .or. i_band == 1 ) .and. isearch(2) >= 0 ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%bottomleft,mpierr)
      else
        ma%bottomleft = mpi_proc_null
      end if
    end if
  
    ma%has_bdytop    = (ma%top    == mpi_proc_null)
    ma%has_bdybottom = (ma%bottom == mpi_proc_null)
    ma%has_bdyright  = (ma%right  == mpi_proc_null)
    ma%has_bdyleft   = (ma%left   == mpi_proc_null)

    jxp =  jx/cpus_per_dim(1)
    iyp =  iy/cpus_per_dim(2)

    global_jstart = ma%location(1)*jxp+1
    global_istart = ma%location(2)*iyp+1

    if ( ma%has_bdytop ) then
      if ( mod(iy,iyp) /= 0 ) then
        global_iend = global_istart+iyp-1
        if ( global_iend < iy ) then
          iyp = iy-global_istart+1
        else
          iyp = mod(iy,iyp)
        end if
      end if
    end if
    if ( ma%has_bdyright ) then
      if ( mod(jx,jxp) /= 0 ) then
        global_jend = global_jstart+jxp-1
        if ( global_jend < jx ) then
          jxp = jx-global_jstart+1
        else
          jxp = mod(jx,jxp)
        end if
      end if
    end if

    ! Check the results
    if ( jxp < 3 .or. iyp < 3 ) then
      write(stderr,*) 'Cannot have one processor with less than 9 points.'
      write(stderr,*) 'Processor ',myid,' has ',jxp*iyp,' (',jxp,'x',iyp,')'
      call fatal(__FILE__,__LINE__,'Too much processors')
    end if

    global_jend = global_jstart+jxp-1
    global_iend = global_istart+iyp-1
    jxpsg  = jxp * nsg
    iypsg  = iyp * nsg
  end subroutine set_nproc

  subroutine broadcast_params
    implicit none

    call mpi_barrier(mycomm,mpierr)

    call mpi_bcast(iy,1,mpi_integer,0,mycomm,mpierr)
    call mpi_bcast(jx,1,mpi_integer,0,mycomm,mpierr)
    call mpi_bcast(kz,1,mpi_integer,0,mycomm,mpierr)
    call mpi_bcast(nsg,1,mpi_integer,0,mycomm,mpierr)
    call mpi_bcast(nveg,1,mpi_integer,0,mycomm,mpierr)

    call mpi_bcast(iproj,6,mpi_character,0,mycomm,mpierr)
    call mpi_bcast(ds,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(ptop,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(clat,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(clon,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(plat,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(plon,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(truelatl,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(truelath,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(i_band,1,mpi_integer,0,mycomm,mpierr)

    call mpi_bcast(domname,64,mpi_character,0,mycomm,mpierr)

    call mpi_bcast(ibyte,1,mpi_integer,0,mycomm,mpierr)

    call mpi_bcast(debug_level,1,mpi_integer,0,mycomm,mpierr)
    call mpi_bcast(dbgfrq,1,mpi_integer,0,mycomm,mpierr)

    call mpi_bcast(nspgx,1,mpi_integer,0,mycomm,mpierr)
    call mpi_bcast(nspgd,1,mpi_integer,0,mycomm,mpierr)
    call mpi_bcast(high_nudge,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(medium_nudge,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(low_nudge,1,mpi_real8,0,mycomm,mpierr)

    call mpi_bcast(calendar,12,mpi_character,0,mycomm,mpierr)
    call mpi_bcast(ical,1,mpi_integer,0,mycomm,mpierr)
    call mpi_bcast(dayspy,1,mpi_real8,0,mycomm,mpierr)
    call mpi_bcast(dpd,1,mpi_real8,0,mycomm,mpierr)

    call mpi_bcast(nsplit,1,mpi_integer,0,mycomm,mpierr)

    call mpi_bcast(ibdyfrq,1,mpi_integer,0,mycomm,mpierr)

    ! Setup all convenience dimensions

    if ( myid /= 0) then
      iym1 = iy - 1
      iym2 = iy - 2
      iym3 = iy - 3
      jxm1 = jx - 1
      jxm2 = jx - 2
      jxm3 = jx - 3
      kzm1 = kz - 1
      kzm2 = kz - 2
      kzp1 = kz + 1
      kzp2 = kz + 2
      kzp3 = kz + 3
      kzp4 = kz + 4
      iysg = iy * nsg
      jxsg = jx * nsg
      iym1sg = (iy-1) * nsg
      jxm1sg = (jx-1) * nsg
      iym2sg = (iy-2) * nsg
      jxm2sg = (jx-2) * nsg
      iym3sg = (iy-3) * nsg
      jxm3sg = (jx-3) * nsg
      nnsg = nsg*nsg
      jdot1 = 1
      jdot2 = jx
      jcross1 = 1
      if ( i_band == 1 ) then
        jcross2 = jx
        jout1 = 1
        jout2 = jx
      else
        jcross2 = jxm1
        jout1 = 2
        jout2 = jxm2
      end if
      idot1 = 1
      idot2 = iy
      icross1 = 1
      icross2 = iym1
      iout1 = 2
      iout2 = iym2
      njcross = jcross2-jcross1+1
      nicross = icross2-icross1+1
      njdot = jdot2-jdot1+1
      nidot = idot2-idot1+1
      njout = jout2-jout1+1
      niout = iout2-iout1+1
    end if

    call mpi_barrier(mycomm,mpierr)

  end subroutine broadcast_params

  subroutine date_bcast(x,from,comm,mpierr)
    type (rcm_time_and_date) , intent(inout) :: x
    integer , intent(in) :: from , comm
    integer , intent(out) :: mpierr
    integer :: lerr
    mpierr = 0
    call mpi_bcast(x%calendar,1,mpi_integer,from,comm,lerr)
    mpierr = mpierr+lerr
    call mpi_bcast(x%days_from_reference,1,mpi_integer,from,comm,lerr)
    mpierr = mpierr+lerr
    call mpi_bcast(x%second_of_day,1,mpi_integer,from,comm,lerr)
    mpierr = mpierr+lerr
  end subroutine date_bcast
!
  subroutine real8_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = mg(global_jstart+j-1,global_istart+i-1)
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_distribute')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            r8vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
        call mpi_send(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_distribute')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine real8_2d_distribute
!
  subroutine real8_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = mg(global_jstart+j-1,global_istart+i-1,k)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_distribute')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              r8vector1(ib) = mg(j,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call mpi_send(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_distribute')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_3d_distribute
!
  subroutine real8_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , isize , jsize , ksize , nsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = mg(global_jstart+j-1,global_istart+i-1,k,n)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_4d_distribute')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_4d_distribute')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                r8vector1(ib) = mg(j,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call mpi_send(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_4d_distribute')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_4d_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_4d_distribute
!
  subroutine real4_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(sp) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = mg(global_jstart+j-1,global_istart+i-1)
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_distribute')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            r4vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
        call mpi_send(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_distribute')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine real4_2d_distribute
!
  subroutine real4_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = mg(global_jstart+j-1,global_istart+i-1,k)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_distribute')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              r4vector1(ib) = mg(j,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call mpi_send(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_distribute')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real4_3d_distribute
!
  subroutine real4_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , isize , jsize , ksize , nsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = mg(global_jstart+j-1,global_istart+i-1,k,n)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_4d_distribute')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_4d_distribute')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                r4vector1(ib) = mg(j,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call mpi_send(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_4d_distribute')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_4d_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real4_4d_distribute
!
  subroutine integer_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    integer , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = mg(global_jstart+j-1,global_istart+i-1)
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_distribute')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            i4vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
        call mpi_send(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_distribute')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = i4vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine integer_2d_distribute
!
  subroutine integer_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = mg(global_jstart+j-1,global_istart+i-1,k)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_distribute')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              i4vector1(ib) = mg(j,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call mpi_send(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_distribute')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine integer_3d_distribute
!
  subroutine integer_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , isize , jsize , ksize , nsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = mg(global_jstart+j-1,global_istart+i-1,k,n)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_4d_distribute')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_4d_distribute')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                i4vector1(ib) = mg(j,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call mpi_send(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_4d_distribute')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_4d_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = i4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine integer_4d_distribute
!
  subroutine real8_2d_sub_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = mg(n,global_jstart+j-1,global_istart+i-1)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        lsize = isize*jsize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_sub_distribute')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_sub_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              r8vector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
        call mpi_send(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      lsize = isize*jsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_sub_distribute')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_sub_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_2d_sub_distribute
!
  subroutine real8_3d_sub_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = mg(n,global_jstart+j-1,global_istart+i-1,k)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_sub_distribute')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_sub_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                r8vector1(ib) = mg(n,j,i,k)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call mpi_send(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_sub_distribute')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_sub_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_3d_sub_distribute
!
  subroutine real4_2d_sub_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = mg(n,global_jstart+j-1,global_istart+i-1)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        lsize = isize*jsize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_sub_distribute')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_sub_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              r4vector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
        call mpi_send(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      lsize = isize*jsize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_sub_distribute')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_sub_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real4_2d_sub_distribute
!
  subroutine real4_3d_sub_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = mg(n,global_jstart+j-1,global_istart+i-1,k)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_sub_distribute')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_sub_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                r4vector1(ib) = mg(n,j,i,k)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call mpi_send(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_sub_distribute')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_sub_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real4_3d_sub_distribute
!
  subroutine integer_2d_sub_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = mg(n,global_jstart+j-1,global_istart+i-1)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        lsize = isize*jsize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_sub_distribute')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_sub_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              i4vector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
        call mpi_send(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      lsize = isize*jsize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_sub_distribute')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_sub_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine integer_2d_sub_distribute
!
  subroutine integer_3d_sub_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = mg(n,global_jstart+j-1,global_istart+i-1,k)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_sub_distribute')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_sub_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                i4vector1(ib) = mg(n,j,i,k)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call mpi_send(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_sub_distribute')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_sub_distribute')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = i4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine integer_3d_sub_distribute
!
  subroutine real8_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          mg(global_jstart+j-1,global_istart+i-1) = ml(j,i)
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_collect')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_collect')
        end if
        call mpi_recv(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_collect')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          r8vector2(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call mpi_send(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_2d_collect
!
  subroutine real8_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            mg(global_jstart+j-1,global_istart+i-1,k) = ml(j,i,k)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_collect')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_collect')
        end if
        call mpi_recv(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              mg(j,i,k) = r8vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_collect')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            r8vector2(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_3d_collect
!
  subroutine real8_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , isize , jsize , ksize , nsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              mg(global_jstart+j-1,global_istart+i-1,k,n) = ml(j,i,k,n)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_4d_collect')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_4d_collect')
        end if
        call mpi_recv(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                mg(j,i,k,n) = r8vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_4d_collect')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_4d_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              r8vector2(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_4d_collect
!
  subroutine real4_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(sp) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          mg(global_jstart+j-1,global_istart+i-1) = ml(j,i)
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_collect')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_collect')
        end if
        call mpi_recv(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = r4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_collect')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          r4vector2(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call mpi_send(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real4_2d_collect
!
  subroutine real4_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            mg(global_jstart+j-1,global_istart+i-1,k) = ml(j,i,k)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_collect')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_collect')
        end if
        call mpi_recv(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              mg(j,i,k) = r4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_collect')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            r4vector2(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real4_3d_collect
!
  subroutine real4_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , isize , jsize , ksize , nsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              mg(global_jstart+j-1,global_istart+i-1,k,n) = ml(j,i,k,n)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_4d_collect')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_4d_collect')
        end if
        call mpi_recv(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                mg(j,i,k,n) = r4vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_4d_collect')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_4d_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              r4vector2(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real4_4d_collect
!
  subroutine integer_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    integer , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          mg(global_jstart+j-1,global_istart+i-1) = ml(j,i)
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_collect')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_collect')
        end if
        call mpi_recv(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = i4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_collect')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          i4vector2(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call mpi_send(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine integer_2d_collect
!
  subroutine integer_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            mg(global_jstart+j-1,global_istart+i-1,k) = ml(j,i,k)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_collect')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_collect')
        end if
        call mpi_recv(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              mg(j,i,k) = i4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_collect')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            i4vector2(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine integer_3d_collect
!
  subroutine integer_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , isize , jsize , ksize , nsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              mg(global_jstart+j-1,global_istart+i-1,k,n) = ml(j,i,k,n)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_4d_collect')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_4d_collect')
        end if
        call mpi_recv(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                mg(j,i,k,n) = i4vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_4d_collect')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_4d_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              i4vector2(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine integer_4d_collect
!
  subroutine real8_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            mg(n,global_jstart+j-1,global_istart+i-1) = ml(n,j,i)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        lsize = isize*jsize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_sub_collect')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_sub_collect')
        end if
        call mpi_recv(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              mg(n,j,i) = r8vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      lsize = isize*jsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_sub_collect')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_sub_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            r8vector2(ib) = ml(n,j,i)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_2d_sub_collect
!
  subroutine real8_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              mg(n,global_jstart+j-1,global_istart+i-1,k) = ml(n,j,i,k)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(r8vector1) ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_sub_collect')
        else if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_sub_collect')
        end if
        call mpi_recv(r8vector1,lsize,mpi_real8,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                mg(n,j,i,k) = r8vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_sub_collect')
      else if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_sub_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              r8vector2(ib) = ml(n,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector2,lsize,mpi_real8,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_3d_sub_collect
!
  subroutine real4_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            mg(n,global_jstart+j-1,global_istart+i-1) = ml(n,j,i)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        lsize = isize*jsize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_sub_collect')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_sub_collect')
        end if
        call mpi_recv(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              mg(n,j,i) = r4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      lsize = isize*jsize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_sub_collect')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_sub_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            r4vector2(ib) = ml(n,j,i)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real4_2d_sub_collect
!
  subroutine real4_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              mg(n,global_jstart+j-1,global_istart+i-1,k) = ml(n,j,i,k)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(r4vector1) ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_sub_collect')
        else if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_sub_collect')
        end if
        call mpi_recv(r4vector1,lsize,mpi_real4,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                mg(n,j,i,k) = r4vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(r4vector2) ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_sub_collect')
      else if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_sub_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              r4vector2(ib) = ml(n,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(r4vector2,lsize,mpi_real4,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real4_3d_sub_collect
!
  subroutine integer_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            mg(n,global_jstart+j-1,global_istart+i-1) = ml(n,j,i)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        lsize = isize*jsize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_sub_collect')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_sub_collect')
        end if
        call mpi_recv(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              mg(n,j,i) = i4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      lsize = isize*jsize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_sub_collect')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_sub_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            i4vector2(ib) = ml(n,j,i)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine integer_2d_sub_collect
!
  subroutine integer_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == 0 ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              mg(n,global_jstart+j-1,global_istart+i-1,k) = ml(n,j,i,k)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call mpi_recv(window,4,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        isize = (window(2)-window(1)+1)*nsg
        jsize = (window(4)-window(3)+1)*nsg
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( .not. associated(i4vector1) ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_sub_collect')
        else if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_sub_collect')
        end if
        call mpi_recv(i4vector1,lsize,mpi_integer,icpu,0, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                mg(n,j,i,k) = i4vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = (i2-i1+1)*nsg
      jsize = (j2-j1+1)*nsg
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( .not. associated(i4vector2) ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_sub_collect')
      else if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_sub_collect')
      end if
      window(1) = global_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_jstart+j1-1
      window(4) = window(3)+jsize-1
      call mpi_send(window,4,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              i4vector2(ib) = ml(n,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(i4vector2,lsize,mpi_integer,0,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine integer_3d_sub_collect
!
  subroutine real8_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: isize , jsize , ssize , j , i , ib , jb
    isize = i2-i1+1
    jsize = j2-j1+1
    if ( ma%right /= mpi_proc_null) then
      ssize = nex*isize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml((jxp-j)+1,i)
          ib = ib + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml(jxp+j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%left /= mpi_proc_null ) then
      ssize = nex*isize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml((1-j),i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(jb) = ml(j,(iyp-i)+1)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,iyp+i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,1-i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(jb) = ml(j,i)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%topleft /= mpi_proc_null ) then
      ssize = nex*nex
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(jb) = ml(j,iyp-i+1)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%topleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%topleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(1-j,iyp+i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(jb) = ml(jxp-j+1,iyp-i+1)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(jxp+j,iyp+i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(jb) = ml(j,i)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(1-j,1-i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
    if ( ma%bottomright /= mpi_proc_null ) then
      ssize = nex*nex
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(jb) = ml(jxp-j+1,i)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottomright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottomright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(jxp+j,1-i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
  end subroutine real8_2d_exchange
!
  subroutine real8_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: isize , jsize , ksize , ssize , j , i , k , ib , jb
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    if ( ma%right /= mpi_proc_null) then
      ssize = nex*isize*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml((jxp-j)+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            ml(jxp+j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%left /= mpi_proc_null ) then
      ssize = nex*isize*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            ml((1-j),i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(jb) = ml(j,(iyp-i)+1,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            ml(j,iyp+i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            ml(j,1-i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(jb) = ml(j,i,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%topleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(jb) = ml(j,iyp-i+1,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%topleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%topleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(1-j,iyp+i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(jb) = ml(jxp-j+1,iyp-i+1,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(jxp+j,iyp+i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(jb) = ml(j,i,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(1-j,1-i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
    end if
    if ( ma%bottomright /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(jb) = ml(jxp-j+1,i,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottomright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottomright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(jxp+j,1-i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
    end if
  end subroutine real8_3d_exchange
!
  subroutine real8_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: isize , jsize , ksize , nsize , ssize , j , i , k , n , ib , jb
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    if ( ma%right /= mpi_proc_null) then
      ssize = nex*isize*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml((jxp-j)+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(jxp+j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%left /= mpi_proc_null ) then
      ssize = nex*isize*ksize*nsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml((1-j),i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              r8vector1(jb) = ml(j,(iyp-i)+1,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              ml(j,iyp+i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              ml(j,1-i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              r8vector1(jb) = ml(j,i,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%topleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(jb) = ml(j,iyp-i+1,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%topleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%topleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(1-j,iyp+i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(jb) = ml(jxp-j+1,iyp-i+1,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(jxp+j,iyp+i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(jb) = ml(j,i,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(1-j,1-i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%bottomright /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(jb) = ml(jxp-j+1,i,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottomright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottomright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(jxp+j,1-i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_4d_exchange
!
  subroutine real8_2d_exchange_left_bottom(ml,nex,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: isize , jsize , ssize , j , i , ib , jb
    isize = i2-i1+1
    jsize = j2-j1+1
    if ( ma%right /= mpi_proc_null) then
      ssize = nex*isize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml((jxp-j)+1,i)
          ib = ib + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%left /= mpi_proc_null ) then
      ssize = nex*isize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_bottom')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml((1-j),i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(jb) = ml(j,(iyp-i)+1)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_bottom')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,1-i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(jb) = ml(jxp-j+1,iyp-i+1)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(1-j,1-i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
  end subroutine real8_2d_exchange_left_bottom
!
  subroutine real8_3d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: isize , jsize , ksize , ssize , j , i , k , ib , jb
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    if ( ma%right /= mpi_proc_null) then
      ssize = nex*isize*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml((jxp-j)+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%left /= mpi_proc_null ) then
      ssize = nex*isize*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_bottom')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            ml((1-j),i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(jb) = ml(j,(iyp-i)+1,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_bottom')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            ml(j,1-i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(jb) = ml(jxp-j+1,iyp-i+1,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(1-j,1-i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
    end if
  end subroutine real8_3d_exchange_left_bottom
!
  subroutine real8_4d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: isize , jsize , ksize , nsize , ssize , j , i , k , n , ib , jb
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    if ( ma%right /= mpi_proc_null) then
      ssize = nex*isize*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml((jxp-j)+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%left /= mpi_proc_null ) then
      ssize = nex*isize*ksize*nsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_bottom')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml((1-j),i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              r8vector1(jb) = ml(j,(iyp-i)+1,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_bottom')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              ml(j,1-i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(jb) = ml(jxp-j+1,iyp-i+1,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(1-j,1-i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_4d_exchange_left_bottom
!
  subroutine real8_2d_exchange_right_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: isize , jsize , ssize , j , i , ib , jb
    isize = i2-i1+1
    jsize = j2-j1+1
    if ( ma%right /= mpi_proc_null) then
      ssize = nex*isize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right_top')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right_top')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml(jxp+j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%left /= mpi_proc_null ) then
      ssize = nex*isize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right_top')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right_top')
      end if
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right_top')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right_top')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,iyp+i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right_top')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right_top')
      end if
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(jb) = ml(j,i)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(jxp+j,iyp+i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      jb = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(jb) = ml(j,i)
          jb = jb + 1
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_2d_exchange_right_top
!
  subroutine real8_3d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: isize , jsize , ksize , ssize , j , i , k , ib , jb
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    if ( ma%right /= mpi_proc_null) then
      ssize = nex*isize*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right_top')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right_top')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            ml(jxp+j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%left /= mpi_proc_null ) then
      ssize = nex*isize*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right_top')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right_top')
      end if
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right_top')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right_top')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            ml(j,iyp+i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right_top')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right_top')
      end if
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(jb) = ml(j,i,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(jxp+j,iyp+i,k) = r8vector2(jb)
            jb = jb + 1
          end do
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(jb) = ml(j,i,k)
            jb = jb + 1
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_3d_exchange_right_top
!
  subroutine real8_4d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: isize , jsize , ksize , nsize , ssize , j , i , k , n , ib , jb
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    if ( ma%right /= mpi_proc_null) then
      ssize = nex*isize*ksize*nsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right_top')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right_top')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(jxp+j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%left /= mpi_proc_null ) then
      ssize = nex*isize*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right_top')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right_top')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right_top')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right_top')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              ml(j,iyp+i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right_top')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right_top')
      end if
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              r8vector1(jb) = ml(j,i,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      call mpi_recv(r8vector2,ssize,mpi_real8,ma%topright,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(jxp+j,iyp+i,k,n) = r8vector2(jb)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      else if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(jb) = ml(j,i,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
      call mpi_send(r8vector1,ssize,mpi_real8,ma%bottomleft,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_4d_exchange_right_top
!
  subroutine real8_bdy_exchange_left_right(ml,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: k1 , k2
    integer :: ksize , k , ib
    ksize = k2-k1+1
    if ( ma%right /= mpi_proc_null) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_left_right')
      else if ( size(r8vector1) < ksize ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_left_right')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_left_right')
      else if ( size(r8vector2) < ksize ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_left_right')
      end if
      ib = 1
      do k = k1 , k2
        r8vector1(ib) = ml(jxp,k)
        ib = ib + 1
      end do
      call mpi_send(r8vector1,ksize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ksize,mpi_real8,ma%right,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        ml(jxp+1,k) = r8vector2(ib)
        ib = ib + 1
      end do
    end if
    if ( ma%left /= mpi_proc_null ) then
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_left_right')
      else if ( size(r8vector2) < ksize ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_left_right')
      end if
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_left_right')
      else if ( size(r8vector1) < ksize ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_left_right')
      end if
      call mpi_recv(r8vector2,ksize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        ml(0,k) = r8vector2(ib)
        ib = ib + 1
      end do
      ib = 1
      do k = k1 , k2
        r8vector1(ib) = ml(1,k)
        ib = ib + 1
      end do
      call mpi_send(r8vector1,ksize,mpi_real8,ma%left,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_bdy_exchange_left_right
!
  subroutine real8_bdy_exchange_top_bottom(ml,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: k1 , k2
    integer :: ksize , k , ib
    ksize = k2-k1+1
    if ( ma%top /= mpi_proc_null) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_top_bottom')
      else if ( size(r8vector1) < ksize ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_top_bottom')
      end if
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_top_bottom')
      else if ( size(r8vector2) < ksize ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_top_bottom')
      end if
      ib = 1
      do k = k1 , k2
        r8vector1(ib) = ml(iyp,k)
        ib = ib + 1
      end do
      call mpi_send(r8vector1,ksize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      call mpi_recv(r8vector2,ksize,mpi_real8,ma%top,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        ml(iyp+1,k) = r8vector2(ib)
        ib = ib + 1
      end do
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      if ( .not. associated(r8vector2) ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_top_bottom')
      else if ( size(r8vector2) < ksize ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_top_bottom')
      end if
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_top_bottom')
      else if ( size(r8vector1) < ksize ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_top_bottom')
      end if
      call mpi_recv(r8vector2,ksize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
      ib = 1
      do k = k1 , k2
        ml(0,k) = r8vector2(ib)
        ib = ib + 1
      end do
      ib = 1
      do k = k1 , k2
        r8vector1(ib) = ml(1,k)
        ib = ib + 1
      end do
      call mpi_send(r8vector1,ksize,mpi_real8,ma%bottom,0, &
                    cartesian_communicator,mpi_status_ignore,mpierr)
    end if
  end subroutine real8_bdy_exchange_top_bottom
!
  subroutine real8_2d_grid_fill(a,b)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(in) :: a
    real(dp) , pointer , dimension(:,:) , intent(out) :: b
    call grid_collect(a,b,1,jxp,1,iyp)
    call mpi_bcast(b,iy*jx,mpi_real8,0,cartesian_communicator,mpierr)
  end subroutine real8_2d_grid_fill
!
! Takes u and v on the cross grid (the same grid as t, qv, qc, etc.)
! and interpolates the u and v to the dot grid.
! This routine sheilds the user of the function from the need to worry
! about the details of the domain decomposition.
!
! Written by Travis A. O'Brien 01/04/11.
!
  subroutine uvcross2dot(ux,vx,ud,vd)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ux , vx
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ud , vd
    integer :: i , j

    ! TODO:  It might make sense to encapsulate the following code
    ! in to a standard routine, since this boundary sending code is
    ! ubiquitous throughout the RegCM code and it is domain
    ! decomposition-dependent.

    ! Send the right-edge of the u/v tendencies to the left
    ! edge of the next process's u/v tendencies (so that
    ! invar%u(i,k,0) holds invar%u(i,k,jxp) of the parallel
    ! chunk next door)

    call exchange_lb(ux,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange_lb(vx,1,jde1,jde2,ide1,ide2,1,kz)

    !
    !     x     x     x     x     x     x
    !
    !        o     o     o     o     o
    !         (i-1,j-1)     (i,j-1)
    !     x     x     x-----x     x     x
    !                 |(i,j)|
    !        o     o  |  o  |  o     o
    !                 |     |
    !     x     x     x-----x     x     x
    !           (i-1,j)     (i,j)
    !
    !        o     o     o     o     o
    !
    !     x     x     x     x     x     x
    !

    ! Perform the bilinear interpolation necessary
    ! to put the u and v variables on the dot grid.

    do i = idi1 , idi2
      do j = jdi1 , jdi2
        ud(j,i,:) =  ud(j,i,:) +             &
          d_rfour*(ux(j,i,:) + ux(j-1,i,:) +   &
                   ux(j,i-1,:) + ux(j-1,i-1,:))
        vd(j,i,:) =  vd(j,i,:) +             &
          d_rfour*(vx(j,i,:) + vx(j-1,i,:) +   &
                   vx(j,i-1,:) + vx(j-1,i-1,:))
      end do
    end do
  end subroutine uvcross2dot
!
  subroutine psc2psd(pc,pd)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(in)  :: pc
    real(dp) , pointer , dimension(:,:) , intent(out) :: pd
    integer :: i , j
    !
    ! Internal points
    !
    do i = idi1 , idi2
      do j = jdi1 , jdi2
        pd(j,i) = (pc(j,i)+pc(j,i-1)+pc(j-1,i)+pc(j-1,i-1))*d_rfour
      end do
    end do
    if ( .not. ma%bandflag ) then
      !
      ! Corner points
      !
      if ( ma%has_bdyleft .and. ma%has_bdybottom ) then
        pd(jde1,ide1) = pc(jce1,ice1)
      end if
      if ( ma%has_bdyleft .and. ma%has_bdytop ) then
        pd(jde1,ide2) = pc(jce1,ice2)
      end if
      if ( ma%has_bdyright .and. ma%has_bdybottom ) then
        pd(jde2,ide1) = pc(jce2,ice1)
      end if
      if ( ma%has_bdyright .and. ma%has_bdytop ) then
        pd(jde2,ide2) = pc(jce2,ice2)
      end if
      !
      ! Boundaries
      !
      if ( ma%has_bdyleft ) then
        do i = idi1 , idi2
          pd(jde1,i) = (pc(jce1,i)+pc(jce1,i-1))*d_half
        end do
      end if
      if ( ma%has_bdyright ) then
        do i = idi1 , idi2
          pd(jde2,i) = (pc(jce2,i)+pc(jce2,i-1))*d_half
        end do
      end if
      if ( ma%has_bdybottom ) then
        do j = jdi1 , jdi2
          pd(j,ide1)  = (pc(j,ice1)+pc(j-1,ice1))*d_half
        end do
      end if
      if ( ma%has_bdytop ) then
        do j = jdi1 , jdi2
          pd(j,ide2) = (pc(j,ice2)+pc(j-1,ice2))*d_half
        end do
      end if
    else
      !
      ! Band (no east and west)
      !
      if ( ma%has_bdybottom ) then
        do j = jce1 , jce2
          pd(j,ide1)  = (pc(j,ice1)+pc(j-1,ice1))*d_half
        end do
      end if
      if ( ma%has_bdytop ) then
        do j = jce1 , jce2
          pd(j,ide2) = (pc(j,ice2)+pc(j-1,ice2))*d_half
        end do
      end if
    end if
  end subroutine psc2psd

  subroutine grid_nc_create_var2d(varname,ldot,val,xvar)
    implicit none
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(dp) , pointer , dimension(:,:) , intent(in) :: val
    type (grid_nc_var2d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(3) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
      xvar%mynx1 = jde1
      xvar%mynx2 = jde2
      xvar%myny1 = ide1
      xvar%myny2 = ide2
    else
      xvar%nx = njcross
      xvar%ny = nicross
      xvar%mynx1 = jce1
      xvar%mynx2 = jce2
      xvar%myny1 = ice1
      xvar%myny2 = ice2
    end if
    xvar%irec = 1
    if ( myid /= 0 ) return
    istat = nf90_create(trim(xvar%varname)//'.nc',nf90_clobber,xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'JX', xvar%nx, idims(1))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'IY', xvar%ny, idims(2))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KTAU', nf90_unlimited, idims(3))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_var(xvar%ncid,xvar%varname,nf90_double,idims,xvar%varid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_enddef(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    call getmem2d(xvar%iobuf,1,xvar%nx,1,xvar%ny,'var2d:iobuf')
  end subroutine grid_nc_create_var2d

  subroutine grid_nc_write_var2d(xvar)
    implicit none
    type (grid_nc_var2d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(3) :: istart , icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    call grid_collect(xvar%val,xvar%iobuf, &
                      xvar%mynx1,xvar%mynx2,xvar%myny1,xvar%myny2)
    if ( myid == 0 ) then
      istart(3) = xvar%irec
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = xvar%ny
      icount(1) = xvar%nx
      istat = nf90_put_var(xvar%ncid,xvar%varid,xvar%iobuf,istart,icount)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
      istat = nf90_sync(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var2d

  subroutine grid_nc_destroy_var2d(xvar)
    implicit none
    type (grid_nc_var2d) , intent(inout) :: xvar
    integer :: istat
    if ( myid == 0 ) then
      istat = nf90_close(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem2d(xvar%iobuf)
  end subroutine grid_nc_destroy_var2d

  subroutine grid_nc_create_var3d(varname,ldot,val,xvar)
    implicit none
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: val
    type (grid_nc_var3d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(4) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
      xvar%mynx1 = jde1
      xvar%mynx2 = jde2
      xvar%myny1 = ide1
      xvar%myny2 = ide2
    else
      xvar%nx = njcross
      xvar%ny = nicross
      xvar%mynx1 = jce1
      xvar%mynx2 = jce2
      xvar%myny1 = ice1
      xvar%myny2 = ice2
    end if
    xvar%nz = size(xvar%val,3)
    xvar%irec = 1
    if ( myid /= 0 ) return
    istat = nf90_create(trim(xvar%varname)//'.nc',nf90_clobber,xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'JX', xvar%nx, idims(1))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'IY', xvar%ny, idims(2))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KZ', xvar%nz, idims(3))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KTAU', nf90_unlimited, idims(4))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_var(xvar%ncid,xvar%varname,nf90_double,idims,xvar%varid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_enddef(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    call getmem3d(xvar%iobuf,1,xvar%nx,1,xvar%ny,1,xvar%nz,'var3d:iobuf')
  end subroutine grid_nc_create_var3d

  subroutine grid_nc_write_var3d(xvar)
    implicit none
    type (grid_nc_var3d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(4) :: istart , icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    call grid_collect(xvar%val,xvar%iobuf, &
                      xvar%mynx1,xvar%mynx2,xvar%myny1,xvar%myny2,1,xvar%nz)
    if ( myid == 0 ) then
      istart(4) = xvar%irec
      istart(3) = 1
      istart(2) = 1
      istart(1) = 1
      icount(4) = 1
      icount(3) = xvar%nz
      icount(2) = xvar%ny
      icount(1) = xvar%nx
      istat = nf90_put_var(xvar%ncid,xvar%varid,xvar%iobuf,istart,icount)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
      istat = nf90_sync(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var3d

  subroutine grid_nc_destroy_var3d(xvar)
    implicit none
    type (grid_nc_var3d) , intent(inout) :: xvar
    integer :: istat
    if ( myid == 0 ) then
      istat = nf90_close(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem3d(xvar%iobuf)
  end subroutine grid_nc_destroy_var3d

  subroutine grid_nc_create_var4d(varname,ldot,val,xvar)
    implicit none
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: val
    type (grid_nc_var4d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(5) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
      xvar%mynx1 = jde1
      xvar%mynx2 = jde2
      xvar%myny1 = ide1
      xvar%myny2 = ide2
    else
      xvar%nx = njcross
      xvar%ny = nicross
      xvar%mynx1 = jce1
      xvar%mynx2 = jce2
      xvar%myny1 = ice1
      xvar%myny2 = ice2
    end if
    xvar%nz = size(xvar%val,3)
    xvar%nl = size(xvar%val,4)
    xvar%irec = 1
    if ( myid /= 0 ) return
    istat = nf90_create(trim(xvar%varname)//'.nc',nf90_clobber,xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'JX', xvar%nx, idims(1))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'IY', xvar%ny, idims(2))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KZ', xvar%nz, idims(3))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'LL', xvar%nl, idims(4))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KTAU', nf90_unlimited, idims(5))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_var(xvar%ncid,xvar%varname,nf90_double,idims,xvar%varid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_enddef(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    call getmem4d(xvar%iobuf,1,xvar%nx,1,xvar%ny,1,xvar%nz, &
                  1,xvar%nl,'var3d:iobuf')
  end subroutine grid_nc_create_var4d

  subroutine grid_nc_write_var4d(xvar)
    implicit none
    type (grid_nc_var4d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(5) :: istart , icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    call grid_collect(xvar%val,xvar%iobuf, &
                      xvar%mynx1,xvar%mynx2,xvar%myny1,xvar%myny2, &
                      1,xvar%nz,1,xvar%nl)
    if ( myid == 0 ) then
      istart(5) = xvar%irec
      istart(4) = 1
      istart(3) = 1
      istart(2) = 1
      istart(1) = 1
      icount(5) = 1
      icount(4) = xvar%nl
      icount(3) = xvar%nz
      icount(2) = xvar%ny
      icount(1) = xvar%nx
      istat = nf90_put_var(xvar%ncid,xvar%varid,xvar%iobuf,istart,icount)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
      istat = nf90_sync(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var4d

  subroutine grid_nc_destroy_var4d(xvar)
    implicit none
    type (grid_nc_var4d) , intent(inout) :: xvar
    integer :: istat
    if ( myid == 0 ) then
      istat = nf90_close(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem4d(xvar%iobuf)
  end subroutine grid_nc_destroy_var4d

end module mod_mppparam

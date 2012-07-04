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

  integer , public , parameter :: nqx = 2
  integer , public , parameter :: iqv = 1
  integer , public , parameter :: iqc = 2

  public :: set_nproc , broadcast_params , date_bcast

  type model_area
    logical :: bandflag
    logical :: hasleft , hasright , hastop , hasbottom
    integer :: left , right , top , bottom
    integer :: ibt1 , ibt2 , ibb1 , ibb2
    integer :: jbl1 , jbl2 , jbr1 , jbr2
  end type model_area

  type deco1d_nc_var2d
    character(len=64) :: varname
    integer :: irec = -1
    integer :: ncid = -1
    integer :: varid = 1
    integer :: nx = 0
    integer :: ny = 0
    real(dp) , pointer , dimension(:,:) :: val => null()
    real(dp) , pointer , dimension(:,:) :: iobuf => null()
  end type deco1d_nc_var2d

  type deco1d_nc_var3d
    character(len=64) :: varname
    integer :: irec = -1
    integer :: ncid = -1
    integer :: varid = 1
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    real(dp) , pointer , dimension(:,:,:) :: val => null()
    real(dp) , pointer , dimension(:,:,:) :: iobuf => null()
  end type deco1d_nc_var3d

  type deco1d_nc_var4d
    character(len=64) :: varname
    integer :: irec = -1
    integer :: ncid = -1
    integer :: varid = 1
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    integer :: nl = 0
    real(dp) , pointer , dimension(:,:,:,:) :: val => null()
    real(dp) , pointer , dimension(:,:,:,:) :: iobuf => null()
  end type deco1d_nc_var4d

  public :: deco1d_nc_var2d , deco1d_nc_var3d , deco1d_nc_var4d

  interface deco1d_nc_create
    module procedure deco1d_nc_create_var2d , &
                     deco1d_nc_create_var3d , &
                     deco1d_nc_create_var4d
  end interface deco1d_nc_create

  interface deco1d_nc_write
    module procedure deco1d_nc_write_var2d , &
                     deco1d_nc_write_var3d , &
                     deco1d_nc_write_var4d
  end interface deco1d_nc_write

  interface deco1d_nc_destroy
    module procedure deco1d_nc_destroy_var2d , &
                     deco1d_nc_destroy_var3d , &
                     deco1d_nc_destroy_var4d
  end interface deco1d_nc_destroy

  public :: deco1d_nc_create , deco1d_nc_write , deco1d_nc_destroy

  interface deco1_scatter
    module procedure deco1_1d_real8_scatter ,   &
                     deco1_2d_real8_scatter ,   &
                     deco1_3d_real8_scatter ,   &
                     deco1_4d_real8_scatter ,   &
                     deco1_2d_real4_scatter ,   &
                     deco1_3d_real4_scatter ,   &
                     deco1_4d_real4_scatter ,   &
                     deco1_2d_integer_scatter , &
                     deco1_3d_integer_scatter , &
                     deco1_4d_integer_scatter
  end interface deco1_scatter

  interface deco1_gather
    module procedure deco1_1d_real8_gather ,   &
                     deco1_2d_real8_gather ,   &
                     deco1_3d_real8_gather ,   &
                     deco1_4d_real8_gather ,   &
                     deco1_2d_real4_gather ,   &
                     deco1_3d_real4_gather ,   &
                     deco1_4d_real4_gather ,   &
                     deco1_2d_integer_gather , &
                     deco1_3d_integer_gather , &
                     deco1_4d_integer_gather
  end interface deco1_gather

  interface subgrid_deco1_scatter
    module procedure subgrid_deco1_2d_real8_scatter ,   &
                     subgrid_deco1_3d_real8_scatter ,   &
                     subgrid_deco1_2d_real4_scatter ,   &
                     subgrid_deco1_3d_real4_scatter ,   &
                     subgrid_deco1_2d_integer_scatter , &
                     subgrid_deco1_3d_integer_scatter
  end interface subgrid_deco1_scatter

  interface subgrid_deco1_gather
    module procedure subgrid_deco1_2d_real8_gather ,   &
                     subgrid_deco1_3d_real8_gather ,   &
                     subgrid_deco1_2d_real4_gather ,   &
                     subgrid_deco1_3d_real4_gather ,   &
                     subgrid_deco1_2d_integer_gather , &
                     subgrid_deco1_3d_integer_gather
  end interface subgrid_deco1_gather

  interface exchange
    module procedure real8_2d_exchange ,  &
                     real8_3d_exchange,   &
                     real8_4d_exchange
  end interface exchange

  interface exchange_left
    module procedure real8_2d_exchange_left ,  &
                     real8_3d_exchange_left,   &
                     real8_4d_exchange_left
  end interface exchange_left

  interface exchange_right
    module procedure real8_2d_exchange_right , &
                     real8_3d_exchange_right,  &
                     real8_4d_exchange_right
  end interface exchange_right

  interface exchange_top
    module procedure real8_2d_exchange_top ,  &
                     real8_3d_exchange_top,   &
                     real8_4d_exchange_top
  end interface exchange_top

  interface exchange_bottom
    module procedure real8_2d_exchange_bottom , &
                     real8_3d_exchange_bottom,  &
                     real8_4d_exchange_bottom
  end interface exchange_bottom

  interface deco1_allgat
    module procedure deco1d_2d_real8_allgat
  end interface deco1_allgat

#ifdef DEBUG
  type(deco1d_nc_var4d) , public :: qqxp
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
  integer :: mpierr
!
  public :: deco1_scatter , deco1_gather
  public :: subgrid_deco1_scatter , subgrid_deco1_gather
  public :: exchange
  public :: exchange_left , exchange_right , exchange_top , exchange_bottom
  public :: deco1_allgat
  public :: uvcross2dot , psc2psd
!
  contains
!
  subroutine set_nproc(ncpu)
    implicit none
    integer , intent(in) :: ncpu
    nproc = ncpu 
    jxp =  jx/nproc
    iyp =  iy
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
  subroutine deco1_1d_real8_scatter(mg,ml)
    implicit none
    real(dp) , pointer , dimension(:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:) , intent(out) :: ml ! model local
    call mpi_scatter(mg,jxp,mpi_real8,ml,jxp,mpi_real8,0,mycomm,mpierr)
  end subroutine deco1_1d_real8_scatter
!
  subroutine deco1_1d_real8_gather(ml,mg)
    implicit none
    real(dp) , pointer , dimension(:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:) , intent(out) :: mg ! model global
    call mpi_gather(ml,jxp,mpi_real8,mg,jxp,mpi_real8,0,mycomm,mpierr)
  end subroutine deco1_1d_real8_gather
!
  subroutine deco1_2d_real8_scatter(mg,ml,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_2d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_2d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_2d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_2d_real8_scatter')
      end if
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          r8vector1(ib) = mg(j,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        ml(j,i) = r8vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine deco1_2d_real8_scatter
!
  subroutine deco1_2d_real8_gather(ml,mg,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_2d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_2d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        r8vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_2d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_2d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          mg(j,i) = r8vector1(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine deco1_2d_real8_gather
!
  subroutine deco1_3d_real8_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_3d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_3d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_3d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_3d_real8_scatter')
      end if
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            r8vector1(ib) = mg(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          ml(j,i,k) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine deco1_3d_real8_scatter
!
  subroutine deco1_3d_real8_gather(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_3d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_3d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          r8vector2(ib) = ml(j,i,k)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_3d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_3d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            mg(j,i,k) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine deco1_3d_real8_gather
!
  subroutine deco1_4d_real8_scatter(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_4d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_4d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_4d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_4d_real8_scatter')
      end if
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              r8vector1(ib) = mg(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            ml(j,i,k,n) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine deco1_4d_real8_scatter
!
  subroutine deco1_4d_real8_gather(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_4d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_4d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            r8vector2(ib) = ml(j,i,k,n)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_4d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_4d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              mg(j,i,k,n) = r8vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine deco1_4d_real8_gather
!
  subroutine deco1_2d_real4_scatter(mg,ml,j1,j2,i1,i2)
    implicit none
    real(sp) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_2d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_2d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_2d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_2d_real4_scatter')
      end if
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          r4vector1(ib) = mg(j,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        ml(j,i) = r4vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine deco1_2d_real4_scatter
!
  subroutine deco1_2d_real4_gather(ml,mg,j1,j2,i1,i2)
    implicit none
    real(sp) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_2d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_2d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        r4vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_2d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_2d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          mg(j,i) = r4vector1(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine deco1_2d_real4_gather
!
  subroutine deco1_3d_real4_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_3d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_3d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_3d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_3d_real4_scatter')
      end if
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            r4vector1(ib) = mg(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          ml(j,i,k) = r4vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine deco1_3d_real4_scatter
!
  subroutine deco1_3d_real4_gather(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_3d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_3d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          r4vector2(ib) = ml(j,i,k)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_3d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_3d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            mg(j,i,k) = r4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine deco1_3d_real4_gather
!
  subroutine deco1_4d_real4_scatter(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_4d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_4d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_4d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_4d_real4_scatter')
      end if
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              r4vector1(ib) = mg(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            ml(j,i,k,n) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine deco1_4d_real4_scatter
!
  subroutine deco1_4d_real4_gather(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_4d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_4d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            r4vector2(ib) = ml(j,i,k,n)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_4d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_4d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              mg(j,i,k,n) = r4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine deco1_4d_real4_gather
!
  subroutine deco1_2d_integer_scatter(mg,ml,j1,j2,i1,i2)
    implicit none
    integer , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_2d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_2d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_2d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_2d_integer_scatter')
      end if
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          i4vector1(ib) = mg(j,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        ml(j,i) = i4vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine deco1_2d_integer_scatter
!
  subroutine deco1_2d_integer_gather(ml,mg,j1,j2,i1,i2)
    implicit none
    integer , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_2d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_2d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        i4vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_2d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_2d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          mg(j,i) = i4vector1(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine deco1_2d_integer_gather
!
  subroutine deco1_3d_integer_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_3d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_3d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_3d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_3d_integer_scatter')
      end if
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            i4vector1(ib) = mg(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          ml(j,i,k) = i4vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine deco1_3d_integer_scatter
!
  subroutine deco1_3d_integer_gather(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_3d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_3d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          i4vector2(ib) = ml(j,i,k)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_3d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_3d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            mg(j,i,k) = i4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine deco1_3d_integer_gather
!
  subroutine deco1_4d_integer_scatter(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_4d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_4d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_4d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_4d_integer_scatter')
      end if
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              i4vector1(ib) = mg(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            ml(j,i,k,n) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine deco1_4d_integer_scatter
!
  subroutine deco1_4d_integer_gather(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_4d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_4d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            i4vector2(ib) = ml(j,i,k,n)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_4d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_4d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              mg(j,i,k,n) = i4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine deco1_4d_integer_gather
!
  subroutine subgrid_deco1_2d_real8_scatter(mg,ml,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*nnsg*jx
    lsize = isize*nnsg*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_2d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_2d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_2d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_2d_real8_scatter')
      end if
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            r8vector1(ib) = mg(nn,j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          ml(nn,j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine subgrid_deco1_2d_real8_scatter
!
  subroutine subgrid_deco1_2d_real8_gather(ml,mg,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*nnsg*jx
    lsize = isize*nnsg*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_2d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_2d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          r8vector2(ib) = ml(nn,j,i)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_2d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_2d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            mg(nn,j,i) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_2d_real8_gather
!
  subroutine subgrid_deco1_3d_real8_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , nn , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*nnsg*jx
    lsize = isize*ksize*nnsg*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_3d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_3d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_3d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_3d_real8_scatter')
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              r8vector1(ib) = mg(nn,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            ml(nn,j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine subgrid_deco1_3d_real8_scatter
!
  subroutine subgrid_deco1_3d_real8_gather(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , nn , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_3d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_3d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            r8vector2(ib) = ml(nn,j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_3d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_3d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              mg(nn,j,i,k) = r8vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_3d_real8_gather
!
  subroutine subgrid_deco1_2d_real4_scatter(mg,ml,j1,j2,i1,i2)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx*nnsg
    lsize = isize*jxp*nnsg
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_2d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_2d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_2d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_2d_real4_scatter')
      end if
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            r4vector1(ib) = mg(nn,j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          ml(nn,j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine subgrid_deco1_2d_real4_scatter
!
  subroutine subgrid_deco1_2d_real4_gather(ml,mg,j1,j2,i1,i2)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx*nnsg
    lsize = isize*jxp*nnsg
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_2d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_2d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          r4vector2(ib) = ml(nn,j,i)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_2d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_2d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            mg(nn,j,i) = r4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_2d_real4_gather
!
  subroutine subgrid_deco1_3d_real4_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , nn , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_3d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_3d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_3d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_3d_real4_scatter')
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              r4vector1(ib) = mg(nn,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            ml(nn,j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine subgrid_deco1_3d_real4_scatter
!
  subroutine subgrid_deco1_3d_real4_gather(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , k , nn , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_3d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_3d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            r4vector2(ib) = ml(nn,j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_3d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_3d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              mg(nn,j,i,k) = r4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_3d_real4_gather
!
  subroutine subgrid_deco1_2d_integer_scatter(mg,ml,j1,j2,i1,i2)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx*nnsg
    lsize = isize*jxp*nnsg
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_2d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_2d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_2d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_2d_integer_scatter')
      end if
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            i4vector1(ib) = mg(nn,j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          ml(nn,j,i) = i4vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine subgrid_deco1_2d_integer_scatter
!
  subroutine subgrid_deco1_2d_integer_gather(ml,mg,j1,j2,i1,i2)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx*nnsg
    lsize = isize*jxp*nnsg
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_2d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_2d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          i4vector2(ib) = ml(nn,j,i)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_2d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_2d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            mg(nn,j,i) = i4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_2d_integer_gather
!
  subroutine subgrid_deco1_3d_integer_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , nn , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_3d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_3d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_3d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_3d_integer_scatter')
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              i4vector1(ib) = mg(nn,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,mpierr)
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            ml(nn,j,i,k) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine subgrid_deco1_3d_integer_scatter
!
  subroutine subgrid_deco1_3d_integer_gather(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ib , i , j , nn , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_3d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_3d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      if ( nproc == 1 ) then
        je = j2
      else
        je = jxp
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            i4vector2(ib) = ml(nn,j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_3d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_3d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,mpierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              mg(nn,j,i,k) = i4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_3d_integer_gather
!
  subroutine real8_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2  , i1 , i2
    call real8_2d_exchange_left(ml,nex,i1,i2)
    call real8_2d_exchange_right(ml,nex,i1,i2)
    call real8_2d_exchange_top(ml,nex,j1,j2)
    call real8_2d_exchange_bottom(ml,nex,j1,j2)
  end subroutine real8_2d_exchange
!
  subroutine real8_2d_exchange_right(ml,nex,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: nex , i1 , i2
    integer :: isize , ssize , i , j , ib
    if ( ma%left == mpi_proc_null .and. ma%right == mpi_proc_null) return
    isize = i2-i1+1
    ssize = nex*isize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right')
    end if
    if ( ma%left /= mpi_proc_null ) then
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%left,2, &
                      r8vector2,ssize,mpi_real8,ma%right,2, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%right /= mpi_proc_null ) then
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml(jxp+j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine real8_2d_exchange_right
!
  subroutine real8_2d_exchange_left(ml,nex,i1,i2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: nex , i1 , i2
    integer :: isize , ssize , j , i , ib
    if ( ma%left == mpi_proc_null .and. ma%right == mpi_proc_null) return
    isize = i2-i1+1
    ssize = nex*isize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left')
    end if
    if ( ma%right /= mpi_proc_null ) then
      ib = 1 
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml((jxp-j)+1,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%right,1, &
                      r8vector2,ssize,mpi_real8,ma%left,1, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%left /= mpi_proc_null ) then
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml((1-j),i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine real8_2d_exchange_left
!
  subroutine real8_2d_exchange_bottom(ml,nex,j1,j2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2
    integer :: jsize , ssize , i , j , jb
    if ( ma%top == mpi_proc_null .and. ma%bottom == mpi_proc_null) return
    jsize = j2-j1+1
    ssize = nex*jsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_bottom')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_bottom')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_bottom')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_bottom')
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(jb) = ml(j,i)
          jb = jb + 1
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%bottom,2, &
                      r8vector2,ssize,mpi_real8,ma%top,2, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%top /= mpi_proc_null ) then
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,iyp+i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
  end subroutine real8_2d_exchange_bottom
!
  subroutine real8_2d_exchange_top(ml,nex,j1,j2)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2
    integer :: jsize , ssize , i , j , jb
    if ( ma%top == mpi_proc_null .and. ma%bottom == mpi_proc_null) return
    jsize = j2-j1+1
    ssize = nex*jsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_top')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_top')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_top')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_top')
    end if
    if ( ma%top /= mpi_proc_null ) then
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(jb) = ml(j,iyp-i+1)
          jb = jb + 1
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%top,2, &
                      r8vector2,ssize,mpi_real8,ma%bottom,2, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%bottom /= mpi_proc_null ) then
      jb = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,1-i) = r8vector2(jb)
          jb = jb + 1
        end do
      end do
    end if
  end subroutine real8_2d_exchange_top
!
  subroutine real8_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2 , i1 , i2 , k1 , k2
    call real8_3d_exchange_left(ml,nex,i1,i2,k1,k2)
    call real8_3d_exchange_right(ml,nex,i1,i2,k1,k2)
    call real8_3d_exchange_top(ml,nex,j1,j2,k1,k2)
    call real8_3d_exchange_bottom(ml,nex,j1,j2,k1,k2)
  end subroutine real8_3d_exchange
!
  subroutine real8_3d_exchange_right(ml,nex,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , i1 , i2 , k1 , k2
    integer :: isize , ksize , ssize , hsize , i , j , k , ib
    if ( ma%left == mpi_proc_null .and. ma%right == mpi_proc_null) return
    isize = i2-i1+1
    ksize = k2-k1+1
    hsize = isize*ksize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right')
    end if
    if ( ma%left /= mpi_proc_null ) then
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%left,2, &
                      r8vector2,ssize,mpi_real8,ma%right,2, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%right /= mpi_proc_null ) then
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
  end subroutine real8_3d_exchange_right
!
  subroutine real8_3d_exchange_left(ml,nex,i1,i2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , i1 , i2 , k1 , k2
    integer :: isize , ksize , ssize , hsize , i , j , k , ib
    if ( ma%left == mpi_proc_null .and. ma%right == mpi_proc_null) return
    isize = i2-i1+1
    ksize = k2-k1+1
    hsize = isize*ksize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left')
    end if
    if ( ma%right /= mpi_proc_null ) then
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml((jxp-j)+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%right,1, &
                      r8vector2,ssize,mpi_real8,ma%left,1, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%left /= mpi_proc_null ) then
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
  end subroutine real8_3d_exchange_left
!
  subroutine real8_3d_exchange_bottom(ml,nex,j1,j2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2 , k1 , k2
    integer :: jsize , ksize , ssize , hsize , i , j , k , jb
    if ( ma%top == mpi_proc_null .and. ma%bottom == mpi_proc_null) return
    jsize = j2-j1+1
    ksize = k2-k1+1
    hsize = jsize*ksize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_bottom')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_bottom')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_bottom')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_bottom')
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(jb) = ml(j,i,k)
            jb = jb + 1
          end do
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%bottom,1, &
                      r8vector2,ssize,mpi_real8,ma%top,1, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%top /= mpi_proc_null ) then
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
  end subroutine real8_3d_exchange_bottom
!
  subroutine real8_3d_exchange_top(ml,nex,j1,j2,k1,k2)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2 , k1 , k2
    integer :: jsize , ksize , ssize , hsize , i , j , k , jb
    if ( ma%top == mpi_proc_null .and. ma%bottom == mpi_proc_null) return
    jsize = j2-j1+1
    ksize = k2-k1+1
    hsize = jsize*ksize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_top')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_top')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_top')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_top')
    end if
    if ( ma%top /= mpi_proc_null ) then
      jb = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(jb) = ml(j,iyp-i+1,k)
            jb = jb + 1
          end do
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%top,1, &
                      r8vector2,ssize,mpi_real8,ma%bottom,1, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%bottom /= mpi_proc_null ) then
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
  end subroutine real8_3d_exchange_top
!
  subroutine real8_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    call real8_4d_exchange_left(ml,nex,i1,i2,k1,k2,n1,n2)
    call real8_4d_exchange_right(ml,nex,i1,i2,k1,k2,n1,n2)
    call real8_4d_exchange_top(ml,nex,j1,j2,k1,k2,n1,n2)
    call real8_4d_exchange_bottom(ml,nex,j1,j2,k1,k2,n1,n2)
  end subroutine real8_4d_exchange
!
  subroutine real8_4d_exchange_right(ml,nex,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , i1 , i2 , k1 , k2 , n1 , n2
    integer :: isize , ssize , ksize , nsize , vsize , hsize , ib
    integer :: i , j , k , n
    if ( ma%left == mpi_proc_null .and. ma%right == mpi_proc_null) return
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    vsize = isize*ksize
    hsize = vsize*nsize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right')
    end if
    if ( ma%left /= mpi_proc_null ) then
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
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%left,2, &
                      r8vector2,ssize,mpi_real8,ma%right,2, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%right /= mpi_proc_null ) then
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
  end subroutine real8_4d_exchange_right
!
  subroutine real8_4d_exchange_left(ml,nex,i1,i2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , i1 , i2 , k1 , k2 , n1 , n2
    integer :: isize , ssize , ksize , nsize , vsize , hsize , ib
    integer :: i , j , k , n
    if ( ma%left == mpi_proc_null .and. ma%right == mpi_proc_null) return
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    vsize = isize*ksize
    hsize = vsize*nsize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left')
    end if
    if ( ma%right /= mpi_proc_null ) then
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
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%right,1, &
                      r8vector2,ssize,mpi_real8,ma%left,1, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%left /= mpi_proc_null ) then
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
  end subroutine real8_4d_exchange_left
!
  subroutine real8_4d_exchange_bottom(ml,nex,j1,j2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2 , k1 , k2 , n1 , n2
    integer :: jsize , ssize , ksize , nsize , vsize , hsize , jb
    integer :: i , j , k , n
    if ( ma%top == mpi_proc_null .and. ma%bottom == mpi_proc_null) return
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    vsize = jsize*ksize
    hsize = vsize*nsize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_bottom')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_bottom')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_bottom')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_bottom')
    end if
    if ( ma%bottom /= mpi_proc_null ) then
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
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%bottom,2, &
                      r8vector2,ssize,mpi_real8,ma%top,2, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%top /= mpi_proc_null ) then
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
  end subroutine real8_4d_exchange_bottom
!
  subroutine real8_4d_exchange_top(ml,nex,j1,j2,k1,k2,n1,n2)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer , intent(in) :: nex , j1 , j2 , k1 , k2 , n1 , n2
    integer :: jsize , ssize , ksize , nsize , vsize , hsize , jb
    integer :: i , j , k , n
    if ( ma%top == mpi_proc_null .and. ma%bottom == mpi_proc_null) return
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    vsize = jsize*ksize
    hsize = vsize*nsize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_top')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_top')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_top')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_top')
    end if
    if ( ma%top /= mpi_proc_null ) then
      jb = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              r8vector1(jb) = ml(j,iyp-i+1,k,n)
              jb = jb + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ma%top,2, &
                      r8vector2,ssize,mpi_real8,ma%bottom,2, &
                      mycomm,mpi_status_ignore,mpierr)
    if ( ma%bottom /= mpi_proc_null ) then
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
  end subroutine real8_4d_exchange_top
!
  subroutine deco1d_2d_real8_allgat(a,b)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(in) :: a
    real(dp) , pointer , dimension(:,:) , intent(out) :: b
    integer :: ii , i , j , ssize , gsize
    ssize = jxp*iy
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'deco1d_real8_allgat')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'deco1d_real8_allgat')
    end if
    gsize = jx*iy
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,gsize,'deco1d_real8_allgat')
    else if ( size(r8vector2) < gsize ) then
      call getmem1d(r8vector2,1,gsize,'deco1d_real8_allgat')
    end if
    ii = 1
    do j = 1 , jxp
      do i = 1 , iy
        r8vector1(ii) = a(j,i)
        ii = ii + 1
      end do
    end do
    call mpi_allgather(r8vector1,jxp*iy,mpi_real8, &
                       r8vector2,jxp*iy,mpi_real8,mycomm,mpierr)
    ii = 1
    do j = 1 , jx
      do i = 1 , iy
        b(j,i) = r8vector2(ii)
        ii = ii + 1
      end do
    end do
  end subroutine deco1d_2d_real8_allgat
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

    call exchange_left(ux,1,ide1,ide2,1,kz)
    call exchange_left(vx,1,ide1,ide2,1,kz)

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
      if ( ma%hasleft .and. ma%hasbottom ) then
        pd(jde1,ide1) = pc(jce1,ice1)
      end if
      if ( ma%hasleft .and. ma%hastop ) then
        pd(jde1,ide2) = pc(jce1,ice2)
      end if
      if ( ma%hasright .and. ma%hasbottom ) then
        pd(jde2,ide1) = pc(jce2,ice1)
      end if
      if ( ma%hasright .and. ma%hastop ) then
        pd(jde2,ide2) = pc(jce2,ice2)
      end if
      !
      ! Boundaries
      !
      if ( ma%hasleft ) then
        do i = idi1 , idi2
          pd(jde1,i) = (pc(jce1,i)+pc(jce1,i-1))*d_half
        end do
      end if
      if ( ma%hasright ) then
        do i = idi1 , idi2
          pd(jde2,i) = (pc(jce2,i)+pc(jce2,i-1))*d_half
        end do
      end if
      if ( ma%hasbottom ) then
        do j = jdi1 , jdi2
          pd(j,ide1)  = (pc(j,ice1)+pc(j-1,ice1))*d_half
        end do
      end if
      if ( ma%hastop ) then
        do j = jdi1 , jdi2
          pd(j,ide2) = (pc(j,ice2)+pc(j-1,ice2))*d_half
        end do
      end if
    else
      !
      ! Band (no east and west)
      !
      if ( ma%hasbottom ) then
        do j = jce1 , jce2
          pd(j,ide1)  = (pc(j,ice1)+pc(j-1,ice1))*d_half
        end do
      end if
      if ( ma%hastop ) then
        do j = jce1 , jce2
          pd(j,ide2) = (pc(j,ice2)+pc(j-1,ice2))*d_half
        end do
      end if
    end if
  end subroutine psc2psd

  subroutine deco1d_nc_create_var2d(varname,ldot,val,xvar)
    implicit none
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(dp) , pointer , dimension(:,:) , intent(in) :: val
    type (deco1d_nc_var2d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(3) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
    else
      xvar%nx = njcross
      xvar%ny = nicross
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
  end subroutine deco1d_nc_create_var2d

  subroutine deco1d_nc_write_var2d(xvar)
    implicit none
    type (deco1d_nc_var2d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(3) :: istart , icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    call deco1_gather(xvar%val,xvar%iobuf,1,xvar%nx,1,xvar%ny)
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
  end subroutine deco1d_nc_write_var2d

  subroutine deco1d_nc_destroy_var2d(xvar)
    implicit none
    type (deco1d_nc_var2d) , intent(inout) :: xvar
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
  end subroutine deco1d_nc_destroy_var2d

  subroutine deco1d_nc_create_var3d(varname,ldot,val,xvar)
    implicit none
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: val
    type (deco1d_nc_var3d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(4) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
    else
      xvar%nx = njcross
      xvar%ny = nicross
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
  end subroutine deco1d_nc_create_var3d

  subroutine deco1d_nc_write_var3d(xvar)
    implicit none
    type (deco1d_nc_var3d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(4) :: istart , icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    call deco1_gather(xvar%val,xvar%iobuf,1,xvar%nx,1,xvar%ny,1,xvar%nz)
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
  end subroutine deco1d_nc_write_var3d

  subroutine deco1d_nc_destroy_var3d(xvar)
    implicit none
    type (deco1d_nc_var3d) , intent(inout) :: xvar
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
  end subroutine deco1d_nc_destroy_var3d

  subroutine deco1d_nc_create_var4d(varname,ldot,val,xvar)
    implicit none
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: val
    type (deco1d_nc_var4d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(5) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
    else
      xvar%nx = njcross
      xvar%ny = nicross
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
  end subroutine deco1d_nc_create_var4d

  subroutine deco1d_nc_write_var4d(xvar)
    implicit none
    type (deco1d_nc_var4d) , intent(inout) :: xvar
    integer :: istat
    integer , dimension(5) :: istart , icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    call deco1_gather(xvar%val,xvar%iobuf,1,xvar%nx,1,xvar%ny, &
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
  end subroutine deco1d_nc_write_var4d

  subroutine deco1d_nc_destroy_var4d(xvar)
    implicit none
    type (deco1d_nc_var4d) , intent(inout) :: xvar
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
  end subroutine deco1d_nc_destroy_var4d

end module mod_mppparam

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

! Distribution of global fields to the processors and collection back.
!
! Template included by mod_mppparam.F90 once per kind: the TK_*
! macros (type, MPI datatype, buffers and routine names) are set by
! the mpp_tk_<kind>.inc headers; routines not meaningful
! for every kind are guarded by TK_REAL (real kinds only) or
! TK_LOGICAL (not for logical).
! It must be included AFTER the type/kind is included

  subroutine TK_2D_DO_DISTRIBUTE(mg,ml,j1,j2,i1,i2,tsize)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: mg
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, iboffset, i, j, icpu, w1, w2, w3, w4
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        ml(j,i) = mg(j,i)
      end do
      return
    end if
    !$acc data create(TK_VEC2,TK_VEC1)
    if ( ccid == ccio ) then
      iboffset = 0
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        w1 = window(1)
        w2 = window(2)
        w3 = window(3)
        w4 = window(4)
        do concurrent ( j = w3:w4, i = w1:w2 )
          ib = iboffset + (i-w1)*(w4-w3+1) + (j-w3) + 1
          TK_VEC1(ib) = mg(j,i)
        end do
        iboffset = iboffset + (w4-w3+1)*(w2-w1+1)
      end do
    end if
    !$acc host_data use_device(TK_VEC1,TK_VEC2)
    call mpi_scatterv(TK_VEC1,wincount,windispl,TK_MPI, &
                      TK_VEC2,tsize,TK_MPI,ccio,mycomm,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    do concurrent ( j = j1:j2, i = i1:i2 )
      ib = (i-i1)*(j2-j1+1) + (j-j1) + 1
      ml(j,i) = TK_VEC2(ib)
    end do
    !$acc end data
  end subroutine TK_2D_DO_DISTRIBUTE

  subroutine TK_2D_DISTRIBUTE(mg,ml,j1,j2,i1,i2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: mg
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call TK_2D_DO_DISTRIBUTE(mg,ml,j1,j2,i1,i2,tsize)
  end subroutine TK_2D_DISTRIBUTE

  subroutine TK_3D_DISTRIBUTE(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: mg
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    TK_TYPE, pointer, contiguous, dimension(:,:) :: mg2 => null( )
    TK_TYPE, pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call TK_2D_DO_DISTRIBUTE(mg2,ml2,j1,j2,i1,i2,tsize)
    end do
  end subroutine TK_3D_DISTRIBUTE

  subroutine TK_4D_DISTRIBUTE(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: mg
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2, n1, n2
    TK_TYPE, pointer, contiguous, dimension(:,:) :: mg2 => null( )
    TK_TYPE, pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize, k, n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1, n2
      do k = k1, k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call TK_2D_DO_DISTRIBUTE(mg2,ml2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine TK_4D_DISTRIBUTE

  subroutine TK_2D_DO_SUB_DISTRIBUTE(mg,ml,j1,j2,i1,i2,tsize,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: mg
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    integer(ik4) :: w1, w2, w3, w4, iboffset
    if ( nproc == 1 ) then
      if ( present(mask) ) then
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          if ( mask(n,j,i) ) ml(n,j,i) = mg(n,j,i)
        end do
      else
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          ml(n,j,i) = mg(n,j,i)
        end do
      end if
      return
    end if
    !$acc data create(TK_VEC2,TK_VEC1)
    if ( ccid == ccio ) then
      iboffset = 0
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        w1 = window(1)
        w2 = window(2)
        w3 = window(3)
        w4 = window(4)
        do concurrent ( n = 1:nnsg, j = w3:w4, i = w1:w2 )
          ib = iboffset + (i-w1)*(nnsg)*(w4-w3+1) + (j-w3)*(nnsg) + n
          TK_VEC1(ib) = mg(n,j,i)
        end do
        iboffset = iboffset + (nnsg)*(w4-w3+1)*(w2-w1+1)
      end do
    end if
    !$acc host_data use_device(TK_VEC1,TK_VEC2)
    call mpi_scatterv(TK_VEC1,wincount,windispl,TK_MPI, &
                      TK_VEC2,tsize,TK_MPI,ccio,mycomm,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    if ( present(mask) ) then
      do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
        ib = (i-i1)*(j2-j1+1)*(nnsg) + (j-j1)*(nnsg) + n
        if ( mask(n,j,i) ) ml(n,j,i) = TK_VEC2(ib)
      end do
    else
      do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
        ib = (i-i1)*(j2-j1+1)*(nnsg) + (j-j1)*(nnsg) + n
        ml(n,j,i) = TK_VEC2(ib)
      end do
    end if
    !$acc end data
  end subroutine TK_2D_DO_SUB_DISTRIBUTE

  subroutine TK_2D_SUB_DISTRIBUTE(mg,ml,j1,j2,i1,i2,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: mg
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call TK_2D_DO_SUB_DISTRIBUTE(mg,ml,j1,j2,i1,i2,tsize,mask)
  end subroutine TK_2D_SUB_DISTRIBUTE

#ifndef TK_LOGICAL
  subroutine TK_3D_SUB_DISTRIBUTE(mg,ml,j1,j2,i1,i2,k1,k2,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: mg
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    TK_TYPE, pointer, contiguous, dimension(:,:,:) :: mg2 => null()
    TK_TYPE, pointer, contiguous, dimension(:,:,:) :: ml2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call TK_2D_DO_SUB_DISTRIBUTE(mg2,ml2,j1,j2,i1,i2,tsize,mask)
    end do
  end subroutine TK_3D_SUB_DISTRIBUTE
#endif

  subroutine TK_2D_DO_COLLECT(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, icpu, ibstart, w1, w2, w3, w4
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        mg(j,i) = ml(j,i)
      end do
      return
    end if
    !$acc data create(TK_VEC2,TK_VEC1)
    do concurrent ( j = j1:j2, i = i1:i2 )
      ib = (i-i1)*(j2-j1+1)+j-j1+1
      TK_VEC2(ib) = ml(j,i)
    end do
    !$acc host_data use_device(TK_VEC2,TK_VEC1)
    call mpi_gatherv(TK_VEC2,tsize,TK_MPI, &
                     TK_VEC1,wincount,windispl,TK_MPI, &
                     ccio,mycomm,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      ibstart = 0
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        w1 = window(1)
        w2 = window(2)
        w3 = window(3)
        w4 = window(4)
        do concurrent ( j = w3:w4, i = w1:w2 )
          ib = ibstart + (i-w1)*(w4-w3+1)+j-w3+1
          mg(j,i) = TK_VEC1(ib)
        end do
        ibstart = ibstart + (w4-w3+1)*(w2-w1+1)
      end do
    end if
    !$acc end data
  end subroutine TK_2D_DO_COLLECT

  subroutine TK_2D_COLLECT(ml,mg,j1,j2,i1,i2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call TK_2D_DO_COLLECT(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine TK_2D_COLLECT

#ifndef TK_LOGICAL
  subroutine TK_3D_COLLECT(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    TK_TYPE, pointer, contiguous, dimension(:,:) :: ml2 => null( )
    TK_TYPE, pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1, k2
      if ( ccid == ccio )  call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call TK_2D_DO_COLLECT(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine TK_3D_COLLECT

  subroutine TK_4D_COLLECT(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2, n1, n2
    TK_TYPE, pointer, contiguous, dimension(:,:) :: ml2 => null( )
    TK_TYPE, pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k, n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1, n2
      do k = k1, k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call TK_2D_DO_COLLECT(ml2,mg2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine TK_4D_COLLECT
#endif

#ifdef TK_REAL
  subroutine TK_2D_3D_COLLECT(ml,mg,j1,j2,i1,i2,k)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4), intent(in), optional :: k
    TK_TYPE, pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: kk
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    kk = 1
    if ( present(k) ) kk = k
    if ( ccid == ccio )  call assignpnt(mg,mg2,kk)
    call TK_2D_DO_COLLECT(ml,mg2,j1,j2,i1,i2,tsize)
  end subroutine TK_2D_3D_COLLECT

  subroutine TK_3D_2D_COLLECT(ml,mg,j1,j2,i1,i2,k)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2, k
    TK_TYPE, pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k)
    call TK_2D_DO_COLLECT(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine TK_3D_2D_COLLECT

  subroutine TK_4D_2D_COLLECT(ml,mg,j1,j2,i1,i2,k,n)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2, k, n
    TK_TYPE, pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k,n)
    call TK_2D_DO_COLLECT(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine TK_4D_2D_COLLECT
#endif

  subroutine TK_2D_DO_SUB_COLLECT(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    integer(ik4) :: w1, w2, w3, w4, iboffset
    if ( nproc == 1 ) then
      do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
        mg(n,j,i) = ml(n,j,i)
      end do
      return
    end if
    !$acc data create(TK_VEC2,TK_VEC1)
    do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
      ib = (i-i1)*(nnsg)*(j2-j1+1) + (j-j1)*(nnsg) + n
      TK_VEC2(ib) = ml(n,j,i)
    end do
    !$acc host_data use_device(TK_VEC2,TK_VEC1)
    call mpi_gatherv(TK_VEC2,tsize,TK_MPI, &
                     TK_VEC1,wincount,windispl,TK_MPI, &
                     ccio,mycomm,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      iboffset = 0
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        w1 = window(1)
        w2 = window(2)
        w3 = window(3)
        w4 = window(4)
        do concurrent ( n = 1:nnsg, j = w3:w4, i = w1:w2 )
          ib = iboffset + (i-w1)*(nnsg)*(w4-w3+1) + (j-w3)*(nnsg) + n
          mg(n,j,i) = TK_VEC1(ib)
        end do
        iboffset = iboffset + (nnsg)*(w4-w3+1)*(w2-w1+1)
      end do
    end if
    !$acc end data
  end subroutine TK_2D_DO_SUB_COLLECT

  subroutine TK_2D_SUB_COLLECT(ml,mg,j1,j2,i1,i2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call TK_2D_DO_SUB_COLLECT(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine TK_2D_SUB_COLLECT

  subroutine TK_3D_SUB_COLLECT(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mg
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    TK_TYPE, pointer, contiguous, dimension(:,:,:) :: ml2 => null( )
    TK_TYPE, pointer, contiguous, dimension(:,:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call TK_2D_DO_SUB_COLLECT(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine TK_3D_SUB_COLLECT

#ifdef TK_REAL
  subroutine TK_2D_GRID_FILL_EXTEND1(a,b)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: a
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: b
    call grid_collect(a,b,max(jce1,lbound(a,1)),min(jce2,ubound(a,1)), &
                          max(ice1,lbound(a,2)),min(ice2,ubound(a,2)))
    ! Extend on a 'fake' filled cross grid (+1)
    b(1,:) = b(2,:)
    b(jx,:) = b(jx-2,:)
    b(jx-1,:) = b(jx-2,:)
    b(:,1) = b(:,2)
    b(:,iy) = b(:,iy-2)
    b(:,iy-1) = b(:,iy-2)
    b(1,1) = b(2,2)
    b(jx,1) = b(jx-1,2)
    b(jx-1,1) = b(jx-1,2)
    b(1,iy) = b(2,iy-2)
    b(1,iy-1) = b(2,iy-2)
    b(jx,iy) = b(jx-2,iy-2)
    b(jx-1,iy-1) = b(jx-2,iy-2)
    call mpi_bcast(b,iy*jx,TK_MPI,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine TK_2D_GRID_FILL_EXTEND1

  subroutine TK_2D_GRID_FILL_EXTEND2(a,b,i1,i2,j1,j2)
    implicit none
    integer(ik4), intent(in) :: i1, i2, j1, j2
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: a
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: b
    call grid_collect(a,b,i1,i2,j1,j2)
    call mpi_bcast(b,product(shape(b)),TK_MPI,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine TK_2D_GRID_FILL_EXTEND2
#endif

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

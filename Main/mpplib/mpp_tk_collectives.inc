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

! Broadcasts, reductions and point to point helpers.
!
! Template included by mod_mppparam.F90 once per real kind: the TK_*
! macros (type, MPI datatype, buffers and routine names) are set by
! mpp_tk_real8.inc or mpp_tk_real4.inc.
! It must be included AFTER the type/kind is included

  subroutine BCAST_TK(rval)
    implicit none
    TK_TYPE, intent(inout) :: rval
    call mpi_bcast(rval,1,TK_MPI,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine BCAST_TK

  subroutine BCAST_ARR_TK(rval)
    implicit none
    TK_TYPE, dimension(:), intent(inout) :: rval
    call mpi_bcast(rval,size(rval),TK_MPI,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine BCAST_ARR_TK

  subroutine BCAST_MATR_TK(rval)
    implicit none
    TK_TYPE, dimension(:,:), intent(inout) :: rval
    call mpi_bcast(rval,size(rval,1)*size(rval,2), &
                   TK_MPI,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine BCAST_MATR_TK

  subroutine SUMALL_TK(rlval,rtval)
    implicit none
    TK_TYPE, intent(in) :: rlval
    TK_TYPE, intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,TK_MPI,mpi_sum,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine SUMALL_TK

  subroutine MAXALL_TK(rlval,rtval)
    implicit none
    TK_TYPE, intent(in) :: rlval
    TK_TYPE, intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,TK_MPI,mpi_max,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine MAXALL_TK

  subroutine MEANALL_1D_TK(rlval,rtval,elem)
    implicit none
    integer(ik4), intent(in) :: elem
    TK_TYPE, dimension(:), intent(in) :: rlval
    TK_TYPE, dimension(:), intent(out) :: rtval
    integer(ik4) :: i
    !$acc data copyin(rlval(1:elem)) copyout(rtval(1:elem))
    !$acc host_data use_device(rlval,rtval)
    call mpi_allreduce(rlval,rtval,elem,TK_MPI,mpi_sum,mycomm,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
    do concurrent( i = 1:elem )
      rtval(i) = rtval(i)/real(nproc,TK_KIND)
    end do
    !$acc end data
  end subroutine MEANALL_1D_TK

  subroutine MEANALL_TK(rlval,rtval)
    implicit none
    TK_TYPE, intent(in) :: rlval
    TK_TYPE, intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,TK_MPI,mpi_sum,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
    rtval = rtval/real(nproc,TK_KIND)
  end subroutine MEANALL_TK

  subroutine MINALL_TK(rlval,rtval)
    implicit none
    TK_TYPE, intent(in) :: rlval
    TK_TYPE, intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,TK_MPI,mpi_min,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine MINALL_TK

  subroutine EXCHANGE_ARRAY_TK(rv1,rv2,isize,icpu,tag1,tag2,srq,rrq)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, tag1, tag2
    TK_TYPE, dimension(isize), intent(in) :: rv1
    TK_TYPE, dimension(isize), intent(inout) :: rv2
    type(mpi_request), intent(out) :: srq, rrq
    call mpi_irecv(rv2,isize,TK_MPI,icpu,tag1, &
                   cartesian_communicator,rrq,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_irecv error.')
    end if
#endif
    call mpi_isend(rv1,isize,TK_MPI,icpu,tag2, &
                  cartesian_communicator,srq,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_isend error.')
    end if
#endif
  end subroutine EXCHANGE_ARRAY_TK

  subroutine CYCLIC_EXCHANGE_ARRAY_TK(rv1,rv2,isize,icpu1,icpu2,itag)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:), intent(in) :: rv1
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: rv2
    integer(ik4), intent(in) :: isize, icpu1, icpu2, itag
    call mpi_sendrecv(rv1,isize,TK_MPI,icpu1,itag, &
                      rv2,isize,TK_MPI,icpu2,itag, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_sendrecv error.')
    end if
#endif
  end subroutine CYCLIC_EXCHANGE_ARRAY_TK

  subroutine TK_COLUMN_REDUCE(m,g)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: m
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: g
    integer(ik4) :: nk
    nk = size(m)
    call mpi_allreduce(m,g,nk,TK_MPI,mpi_sum,&
                       cartesian_column_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine TK_COLUMN_REDUCE

  subroutine TK_ROW_REDUCE(m,g)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: m
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: g
    integer(ik4) :: nk
    nk = size(m)
    call mpi_allreduce(m,g,nk,TK_MPI,mpi_sum,&
                       cartesian_row_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine TK_ROW_REDUCE

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

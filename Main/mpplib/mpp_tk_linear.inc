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

! Masked transfers between the cartesian and linear decompositions.
!
! Template included by mod_mppparam.F90 once per kind: the TK_*
! macros (type, MPI datatype, buffers and routine names) are set by
! the mpp_tk_<kind>.inc headers; routines not meaningful
! for every kind are guarded by TK_REAL (real kinds only) or
! TK_LOGICAL (not for logical).
! It must be included AFTER the type/kind is included

#ifdef TK_REAL
  subroutine CL_SETUP_TK(cl,gmask,sgmask,lrev)
    implicit none
    type(masked_comm), intent(inout) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: gmask
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: sgmask
    logical, optional, intent(in) :: lrev
    integer(ik4) :: ncart_tot_g, ncart_tot_sg
    integer(ik4) :: i, j, n, iv
    if ( .not. associated(cl%linear_npoint_g) ) then
      call mpi_comm_dup(mycomm,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_comm_dup error.')
      end if
#endif
      call getmem(cl%linear_npoint_g,1,nproc,'cl:linear_npoint_g')
      call getmem(cl%linear_displ_g,1,nproc,'cl:linear_displ_g')
      call getmem(cl%cartesian_npoint_g,1,nproc,'cl:cartesian_npoint_g')
      call getmem(cl%cartesian_displ_g,1,nproc,'cl:cartesian_displ_g')
      call getmem(cl%linear_npoint_sg,1,nproc,'cl:linear_npoint_sg')
      call getmem(cl%linear_displ_sg,1,nproc,'cl:linear_displ_sg')
      call getmem(cl%cartesian_npoint_sg,1,nproc,'cl:cartesian_npoint_sg')
      call getmem(cl%cartesian_displ_sg,1,nproc,'cl:cartesian_displ_sg')
      call getmem(cl%gmask,jci1,jci2,ici1,ici2,'cl:gmask')
      call getmem(cl%gmask_id,jci1,jci2,ici1,ici2,'cl:gmask_id')
      call getmem(cl%sgmask,1,nnsg,jci1,jci2,ici1,ici2,'cl:sgmask')
      call getmem(cl%sgmask_id,1,nnsg,jci1,jci2,ici1,ici2,'cl:sgmask_id')
      call getmem(cl%global_gmask,jout1,jout2,iout1,iout2,'cl:global_gmask')
      call getmem(cl%global_gmask_id,jout1,jout2, &
                                   iout1,iout2,'cl:global_gmask_id')
      call getmem(cl%global_sgmask,1,nnsg,jout1,jout2, &
                                   iout1,iout2,'cl:global_sgmask')
      call getmem(cl%global_sgmask_id,1,nnsg,jout1,jout2, &
                                   iout1,iout2,'cl:global_sgmask_id')
      call getmem(cl%global_out_sgmask,joutsg1,joutsg2, &
                                   ioutsg1,ioutsg2,'cl:global_out_sgmask')
    end if
    cl%gmask = gmask(jci1:jci2,ici1:ici2) > 0.0
    cl%sgmask = sgmask(1:nnsg,jci1:jci2,ici1:ici2) > 0.0
    if ( present(lrev) ) then
      if ( lrev ) then
        cl%gmask = .not. cl%gmask
        cl%sgmask = .not. cl%sgmask
      end if
    end if
    call grid_collect(cl%gmask,cl%global_gmask,jci1,jci2,ici1,ici2)
    call subgrid_collect(cl%sgmask,cl%global_sgmask,jci1,jci2,ici1,ici2)
    call reorder_glb_subgrid(cl%global_sgmask,cl%global_out_sgmask)
    ncart_tot_g = count(cl%gmask)
    ncart_tot_sg = count(cl%sgmask)
    call clset(ncart_tot_g,ncart_tot_sg,cl)
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          cl%global_gmask_id(j,i) = iv
          iv = iv + 1
        else
          cl%global_gmask_id(j,i) = 0
        end if
      end do
    end do
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            cl%global_sgmask_id(n,j,i) = iv
            iv = iv + 1
          else
            cl%global_sgmask_id(n,j,i) = 0
          end if
        end do
      end do
    end do
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            cl%sgmask_id(n,j,i) = iv
            iv = iv + 1
          else
            cl%sgmask_id(n,j,i) = 0
          end if
        end do
      end do
    end do
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          cl%gmask_id(j,i) = iv
          iv = iv + 1
        else
          cl%gmask_id(j,i) = 0
        end if
      end do
    end do
  end subroutine CL_SETUP_TK
#endif

  subroutine CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call mypack(cl,matrix,vector)
      return
    end if
    !$acc data create(TK_VEC1,TK_VEC2)
    if ( nval > 0 ) then
      call mypack(cl,matrix,TK_VEC1)
    end if
    !$acc host_data use_device(TK_VEC1,TK_VEC2)
    call mpi_gatherv(TK_VEC1,nval,TK_MPI,                               &
                     TK_VEC2,cl%cartesian_npoint_sg,cl%cartesian_displ_sg, &
                     TK_MPI,ccio,cartesian_communicator,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    !$acc host_data use_device(TK_VEC2,vector)
    call mpi_scatterv(TK_VEC2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      TK_MPI,vector,npt,TK_MPI,                   &
                      iocpu,cl%linear_communicator,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    !$acc end data
  end subroutine CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID

  subroutine LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix)
      return
    end if
    !$acc data create(TK_VEC1,TK_VEC2)
    !$acc host_data use_device(vector,TK_VEC2)
    call mpi_gatherv(vector,npt,TK_MPI,                             &
                     TK_VEC2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     TK_MPI,iocpu,cl%linear_communicator,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    !$acc host_data use_device(TK_VEC2,TK_VEC1)
    call mpi_scatterv(TK_VEC2,cl%cartesian_npoint_sg, &
                      cl%cartesian_displ_sg,TK_MPI,  &
                      TK_VEC1,nval,TK_MPI,         &
                      ccio,cartesian_communicator,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    if ( nval > 0 ) then
      call myunpack(cl,TK_VEC1,matrix)
    end if
    !$acc end data
  end subroutine LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID

  subroutine CARTESIAN_TO_LINEAR_TK_GRID_SUBGRID(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: matrix
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    !$acc data create(TK_SUBGRID)
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      TK_SUBGRID(n,j,i) = matrix(j,i)
    end do
    call CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID(cl,TK_SUBGRID,vector)
    !$acc end data
  end subroutine CARTESIAN_TO_LINEAR_TK_GRID_SUBGRID

  subroutine GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: npt
    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector)
      return
    end if
    !$acc data create(TK_GLOBAL_SUBGRID)
    call subgrid_collect(matrix,TK_GLOBAL_SUBGRID,jci1,jci2,ici1,ici2)
    !$acc data create(TK_VEC1)
    call mypack_global(cl,TK_GLOBAL_SUBGRID,TK_VEC1)
    npt  = cl%linear_npoint_sg(myid+1)
    !$acc host_data use_device(TK_VEC1,vector)
    call mpi_scatterv(TK_VEC1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      TK_MPI,vector,npt,TK_MPI,                   &
                      iocpu,cl%linear_communicator,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    !$acc end data
    !$acc end data
  end subroutine GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID

  subroutine LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID(cl,vector,matrix)
    !@acc use nvtx
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    !@acc call nvtxStartRange("l2g_"//TK_NAME//"_ss")
    if ( nproc == 1 ) then
      call myunpack_global(cl,vector,matrix)
      return
    end if
    !$acc data create(TK_GLOBAL_SUBGRID)
    !$acc data create(TK_VEC1)
    !@acc call nvtxStartRange("mpi_gatherv")
    !$acc host_data use_device(vector,TK_VEC1)
    call mpi_gatherv(vector,cl%linear_npoint_sg(myid+1),TK_MPI,  &
                     TK_VEC1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     TK_MPI,iocpu,cl%linear_communicator,mpierr)
    !$acc end host_data
    !@acc call nvtxEndRange
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    !@acc call nvtxStartRange("myunpack_global")
    call myunpack_global(cl,TK_VEC1,TK_GLOBAL_SUBGRID)
    !@acc call nvtxEndRange
    !$acc end data
    !@acc call nvtxStartRange("subgrid_distribute")
    call subgrid_distribute(TK_GLOBAL_SUBGRID,matrix, &
            jci1,jci2,ici1,ici2,cl%sgmask)
    !@acc call nvtxEndRange
    !$acc end data
    !@acc call nvtxEndRange
  end subroutine LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID

  subroutine GLOBAL_TO_LINEAR_TK_GRID_SUBGRID(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: matrix
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    !$acc data create(TK_SUBGRID)
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      TK_SUBGRID(n,j,i) = matrix(j,i)
    end do
    call GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID(cl,TK_SUBGRID,vector)
    !$acc end data
  end subroutine GLOBAL_TO_LINEAR_TK_GRID_SUBGRID

  subroutine MYPACK_TK_GRID(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4) :: i, j, iv
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      iv = cl%gmask_id(j,i)
      if ( iv > 0 ) then
        vector(iv) = matrix(j,i)
      end if
    end do
  end subroutine MYPACK_TK_GRID

  subroutine MYPACK_TK_SUBGRID(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      iv = cl%sgmask_id(n,j,i)
      if ( iv > 0 ) then
        vector(iv) = matrix(n,j,i)
      end if
    end do
  end subroutine MYPACK_TK_SUBGRID

  subroutine MYUNPACK_TK_GRID(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      iv = cl%gmask_id(j,i)
      if ( iv > 0 ) then
        matrix(j,i) = vector(iv)
      end if
    end do
  end subroutine MYUNPACK_TK_GRID

  subroutine MYUNPACK_TK_SUBGRID(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      iv = cl%sgmask_id(n,j,i)
      if ( iv > 0 ) then
        matrix(n,j,i) = vector(iv)
      end if
    end do
  end subroutine MYUNPACK_TK_SUBGRID

  subroutine MYPACK_GLOBAL_TK_GRID(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: matrix
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, iv
    do concurrent ( j = jout1:jout2, i = iout1:iout2 )
      iv = cl%global_gmask_id(j,i)
      if ( iv > 0 ) then
        vector(iv) = matrix(j,i)
      end if
    end do
  end subroutine MYPACK_GLOBAL_TK_GRID

  subroutine MYPACK_GLOBAL_TK_SUBGRID(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    do concurrent ( n = 1:nnsg, j = jout1:jout2, i = iout1:iout2 )
      iv = cl%global_sgmask_id(n,j,i)
      if ( iv > 0 ) then
        vector(iv) = matrix(n,j,i)
      end if
    end do
  end subroutine MYPACK_GLOBAL_TK_SUBGRID

  subroutine MYUNPACK_GLOBAL_TK_GRID(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    do concurrent ( j = jout1:jout2, i = iout1:iout2 )
      iv = cl%global_gmask_id(j,i)
      if ( iv > 0 ) then
        matrix(j,i) = vector(iv)
      end if
    end do
  end subroutine MYUNPACK_GLOBAL_TK_GRID

  subroutine MYUNPACK_GLOBAL_TK_SUBGRID(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    do concurrent ( n = 1:nnsg, j = jout1:jout2, i = iout1:iout2 )
      iv = cl%global_sgmask_id(n,j,i)
      if ( iv > 0 ) then
        matrix(n,j,i) = vector(iv)
      end if
    end do
  end subroutine MYUNPACK_GLOBAL_TK_SUBGRID

#ifdef TK_REAL
  subroutine CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID_4D(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: vector
    integer(ik4) :: nval, npt, nlev, k
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    nlev = size(matrix,4)
    if ( nproc == 1 ) then
      call mypack(cl,matrix,vector,nlev)
      return
    end if
    !$acc data create(TK_VEC1,TK_VEC2)
    do k = 1, nlev
      if ( nval > 0 ) then
        call mypack(cl,matrix,TK_VEC1,k)
      end if
      !$acc host_data use_device(TK_VEC1,TK_VEC2)
      call mpi_gatherv(TK_VEC1,nval,TK_MPI,                               &
                       TK_VEC2,cl%cartesian_npoint_sg,cl%cartesian_displ_sg, &
                       TK_MPI,ccio,cartesian_communicator,mpierr)
      !$acc end host_data
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      !$acc host_data use_device(TK_VEC2,vector)
      call mpi_scatterv(TK_VEC2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                        TK_MPI,vector(:,k),npt,TK_MPI,              &
                        iocpu,cl%linear_communicator,mpierr)
      !$acc end host_data
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
    end do
    !$acc end data
  end subroutine CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID_4D

  subroutine LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID_4D(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4) :: nval, npt, nlev, k
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    nlev = size(vector,2)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix,nlev)
      return
    end if
    !$acc data create(TK_VEC1,TK_VEC2)
    do k = 1, nlev
      !$acc host_data use_device(vector,TK_VEC2)
      call mpi_gatherv(vector(:,k),npt,TK_MPI,                        &
                       TK_VEC2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                       TK_MPI,iocpu,cl%linear_communicator,mpierr)
      !$acc end host_data
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      !$acc host_data use_device(TK_VEC2,TK_VEC1)
      call mpi_scatterv(TK_VEC2,cl%cartesian_npoint_sg, &
                        cl%cartesian_displ_sg,TK_MPI,  &
                        TK_VEC1,nval,TK_MPI,         &
                        ccio,cartesian_communicator,mpierr)
      !$acc end host_data
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
      if ( nval > 0 ) then
        call myunpack(cl,TK_VEC1,matrix,k)
      end if
    end do
    !$acc end data
  end subroutine LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID_4D

  subroutine GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID_4D(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: vector
    integer(ik4) :: npt, k, nlev
    integer(ik4) :: i, j, n
    nlev = size(matrix,4)
    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector,nlev)
      return
    end if
    npt  = cl%linear_npoint_sg(myid+1)
    !$acc data create(TK_SUBGRID,TK_GLOBAL_SUBGRID)
    do k = 1, nlev
      do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
        TK_SUBGRID(n,j,i) = matrix(n,j,i,k)
      end do
      call subgrid_collect(TK_SUBGRID,TK_GLOBAL_SUBGRID,jci1,jci2,ici1,ici2)
      !$acc data create(TK_VEC1)
      call mypack_global(cl,TK_GLOBAL_SUBGRID,TK_VEC1)
      !$acc host_data use_device(TK_VEC1,vector)
      call mpi_scatterv(TK_VEC1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                        TK_MPI,vector(:,k),npt,TK_MPI,              &
                        iocpu,cl%linear_communicator,mpierr)
      !$acc end host_data
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
      !$acc end data
    end do
    !$acc end data
  end subroutine GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID_4D

  subroutine LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID_4D(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4) :: npt, nlev, k
    integer(ik4) :: i, j, n
    nlev = size(vector,2)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix,nlev)
      return
    end if
    npt  = cl%linear_npoint_sg(myid+1)
    !$acc data create(TK_GLOBAL_SUBGRID,TK_SUBGRID)
    do k = 1, nlev
      !$acc data create(TK_VEC1)
      !$acc host_data use_device(vector,TK_VEC1)
      call mpi_gatherv(vector(:,k),npt,TK_MPI,                        &
                       TK_VEC1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                       TK_MPI,iocpu,cl%linear_communicator,mpierr)
      !$acc end host_data
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      call myunpack_global(cl,TK_VEC1,TK_GLOBAL_SUBGRID)
      !$acc end data
      call subgrid_distribute(TK_GLOBAL_SUBGRID,TK_SUBGRID, &
              jci1,jci2,ici1,ici2,cl%sgmask)
      do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
        if ( cl%sgmask(n,j,i) ) matrix(n,j,i,k) = TK_SUBGRID(n,j,i)
      end do
    end do
    !$acc end data
  end subroutine LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID_4D

  subroutine MYPACK_TK_SUBGRID_4D(cl,matrix,vector,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
        iv = cl%sgmask_id(n,j,i)
        if ( iv > 0 ) then
          vector(iv,k) = matrix(n,j,i,k)
        end if
      end do
    end do
  end subroutine MYPACK_TK_SUBGRID_4D

  subroutine MYPACK_TK_SUBGRID_SLICE(cl,matrix,vector,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      iv = cl%sgmask_id(n,j,i)
      if ( iv > 0 ) then
        vector(iv) = matrix(n,j,i,k)
      end if
    end do
  end subroutine MYPACK_TK_SUBGRID_SLICE

  subroutine MYUNPACK_TK_SUBGRID_4D(cl,vector,matrix,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
        iv = cl%sgmask_id(n,j,i)
        if ( iv > 0 ) then
          matrix(n,j,i,k) = vector(iv,k)
        end if
      end do
    end do
  end subroutine MYUNPACK_TK_SUBGRID_4D

  subroutine MYUNPACK_TK_SUBGRID_SLICE(cl,vector,matrix,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      iv = cl%sgmask_id(n,j,i)
      if ( iv > 0 ) then
        matrix(n,j,i,k) = vector(iv)
      end if
    end do
  end subroutine MYUNPACK_TK_SUBGRID_SLICE

  subroutine MYPACK_GLOBAL_TK_SUBGRID_4D(cl,matrix,vector,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      do concurrent ( n = 1:nnsg, j = jout1:jout2, i = iout1:iout2 )
        iv = cl%global_sgmask_id(n,j,i)
        if ( iv > 0 ) then
          vector(iv,k) = matrix(n,j,i,k)
        end if
      end do
    end do
  end subroutine MYPACK_GLOBAL_TK_SUBGRID_4D

  subroutine MYPACK_GLOBAL_TK_SUBGRID_SLICE(cl,matrix,vector,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(inout) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    do concurrent ( n = 1:nnsg, j = jout1:jout2, i = iout1:iout2 )
      iv = cl%global_sgmask_id(n,j,i)
      if ( iv > 0 ) then
        vector(iv) = matrix(n,j,i,k)
      end if
    end do
  end subroutine MYPACK_GLOBAL_TK_SUBGRID_SLICE

  subroutine MYUNPACK_GLOBAL_TK_SUBGRID_4D(cl,vector,matrix,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      do concurrent ( n = 1:nnsg, j = jout1:jout2, i = iout1:iout2 )
        iv = cl%global_sgmask_id(n,j,i)
        if ( iv > 0 ) then
          matrix(n,j,i,k) = vector(iv,k)
        end if
      end do
    end do
  end subroutine MYUNPACK_GLOBAL_TK_SUBGRID_4D

  subroutine MYUNPACK_GLOBAL_TK_SUBGRID_SLICE(cl,vector,matrix,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    TK_TYPE, pointer, contiguous, dimension(:), intent(in) :: vector
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    do concurrent ( n = 1:nnsg, j = jout1:jout2, i = iout1:iout2 )
      iv = cl%global_sgmask_id(n,j,i)
      if ( iv > 0 ) then
        matrix(n,j,i,k) = vector(iv)
      end if
    end do
  end subroutine MYUNPACK_GLOBAL_TK_SUBGRID_SLICE
#endif

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

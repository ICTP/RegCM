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

! Halo exchanges between neighbouring processors.
!
! Template included by mod_mppparam.F90 once per real kind: the TK_*
! macros (type, MPI datatype, buffers and routine names) are set by
! mpp_tk_real8.inc or mpp_tk_real4.inc.
! It must be included AFTER the type/kind is included

  subroutine TK_2D_EXCHANGE(ml,nex,j1,j2,i1,i2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: nx, ny
    integer(ik4) :: lb, rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey
    logical :: do_lt, do_rt, do_bt, do_tp
    !@acc call nvtxStartRange(TK_NAME//"_2d_exchange")
    do_lt = .true.
    do_rt = .true.
    do_bt = .true.
    do_tp = .true.
    if ( ma%left   == mpi_proc_null ) do_lt = .false.
    if ( ma%right  == mpi_proc_null ) do_rt = .false.
    if ( ma%bottom == mpi_proc_null ) do_bt = .false.
    if ( ma%top    == mpi_proc_null ) do_tp = .false.
    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    rb = nex
    if ( ma%left  == mpi_proc_null ) lb = 0
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdatax
    TK_TYPE, dimension(ndy) :: sdatay
    TK_TYPE, dimension(ndx), volatile :: rdatax
    TK_TYPE, dimension(ndy), volatile :: rdatay
    integer(ik4) :: id, iex, i, j

    !$acc data copy(ml) create(sdatax,sdatay,rdatax,rdatay)
    do concurrent ( i = i1:i2, iex = 1:nex )
      id = (iex-1)*tx + 1 + (i-i1)
      sdatax(id      ) = ml(j1+(iex-1),i)
      sdatax(id+sizex) = ml(j2-(iex-1),i)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdatax,rdatax)
    call mpi_neighbor_alltoallv(sdatax, counts, displs, TK_MPI, &
        rdatax, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    do concurrent ( i = i1:i2, iex = 1:nex )
      id = (iex-1)*tx + 1 + (i-i1)
      if ( do_lt ) ml(j1-iex,i) = rdatax(id      )
      if ( do_rt ) ml(j2+iex,i) = rdatax(id+sizex)
    end do

    do concurrent ( j = j1-lb:j2+rb, iex = 1:nex )
      id = (iex-1)*ty + 1 + (j-j1+lb)
      sdatay(id      ) = ml(j,i1+(iex-1))
      sdatay(id+sizey) = ml(j,i2-(iex-1))
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdatay,rdatay)
    call mpi_neighbor_alltoallv(sdatay, counts, displs, TK_MPI, &
        rdatay, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    do concurrent ( j = j1-lb:j2+rb, iex = 1:nex )
      id = (iex-1)*ty + 1 + (j-j1+lb)
      if ( do_bt ) ml(j,i1-iex) = rdatay(id      )
      if ( do_tp ) ml(j,i2+iex) = rdatay(id+sizey)
    end do
    !$acc end data

    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_2D_EXCHANGE

  subroutine TK_3D_EXCHANGE(ml,nex,j1,j2,i1,i2,k1,k2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: nx, ny, nk
    integer(ik4) :: lb, rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey
    logical :: do_lt, do_rt, do_bt, do_tp
    !@acc call nvtxStartRange(TK_NAME//"_3d_exchange")
    do_lt = .true.
    do_rt = .true.
    do_bt = .true.
    do_tp = .true.
    if ( ma%left   == mpi_proc_null ) do_lt = .false.
    if ( ma%right  == mpi_proc_null ) do_rt = .false.
    if ( ma%bottom == mpi_proc_null ) do_bt = .false.
    if ( ma%top    == mpi_proc_null ) do_tp = .false.
    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    rb = nex
    if ( ma%left  == mpi_proc_null ) lb = 0
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdatax
    TK_TYPE, dimension(ndy) :: sdatay
    TK_TYPE, dimension(ndx), volatile :: rdatax
    TK_TYPE, dimension(ndy), volatile :: rdatay
    integer(ik4) :: id, iex, i, j, k

    !$acc data copy(ml) create(sdatax,sdatay,rdatax,rdatay)
    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*tx + 1 + (i-i1)
      sdatax(id      ) = ml(j1+(iex-1),i,k)
      sdatax(id+sizex) = ml(j2-(iex-1),i,k)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdatax,rdatax)
    call mpi_neighbor_alltoallv(sdatax, counts, displs, TK_MPI, &
        rdatax, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif
    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*tx + 1 + (i-i1)
      if ( do_lt ) ml(j1-iex,i,k) = rdatax(id      )
      if ( do_rt ) ml(j2+iex,i,k) = rdatax(id+sizex)
    end do

    do concurrent ( j = j1-lb:j2+rb, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*ty + 1 + (j-j1+lb)
      sdatay(id      ) = ml(j,i1+(iex-1),k)
      sdatay(id+sizey) = ml(j,i2-(iex-1),k)
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdatay,rdatay)
    call mpi_neighbor_alltoallv(sdatay, counts, displs, TK_MPI, &
        rdatay, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif
    do concurrent ( j = j1-lb:j2+rb, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*ty + 1 + (j-j1+lb)
      if ( do_bt ) ml(j,i1-iex,k) = rdatay(id      )
      if ( do_tp ) ml(j,i2+iex,k) = rdatay(id+sizey)
    end do
    !$acc end data
    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_3D_EXCHANGE

  subroutine TK_4D_EXCHANGE(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2, i1, i2, k1, k2, n1, n2
    integer(ik4) :: nx, ny, nk, nn
    integer(ik4) :: lb, rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey
    logical :: do_lt, do_rt, do_bt, do_tp
    !@acc call nvtxStartRange(TK_NAME//"_4d_exchange")
    do_lt = .true.
    do_rt = .true.
    do_bt = .true.
    do_tp = .true.
    if ( ma%left   == mpi_proc_null ) do_lt = .false.
    if ( ma%right  == mpi_proc_null ) do_rt = .false.
    if ( ma%bottom == mpi_proc_null ) do_bt = .false.
    if ( ma%top    == mpi_proc_null ) do_tp = .false.
    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    rb = nex
    if ( ma%left  == mpi_proc_null ) lb = 0
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdatax
    TK_TYPE, dimension(ndy) :: sdatay
    TK_TYPE, dimension(ndx), volatile :: rdatax
    TK_TYPE, dimension(ndy), volatile :: rdatay
    integer(ik4) :: id, iex, i, j, k, n

    !$acc data copy(ml) create(sdatax,sdatay,rdatax,rdatay)
    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2, n = n1:n2 )
      id = (((n-n1)*nk+(k-k1))*nex+iex-1)*tx + 1 + (i-i1)
      sdatax(id      ) = ml(j1+(iex-1),i,k,n)
      sdatax(id+sizex) = ml(j2-(iex-1),i,k,n)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdatax,rdatax)
    call mpi_neighbor_alltoallv(sdatax, counts, displs, TK_MPI, &
        rdatax, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif
    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2, n = n1:n2 )
      id = (((n-n1)*nk+(k-k1))*nex+iex-1)*tx + 1 + (i-i1)
      if ( do_lt ) ml(j1-iex,i,k,n) = rdatax(id      )
      if ( do_rt ) ml(j2+iex,i,k,n) = rdatax(id+sizex)
    end do

    do concurrent ( j = j1-lb:j2+rb, iex = 1:nex, k = k1:k2, n = n1:n2 )
      id = (((n-n1)*nk+(k-k1))*nex+iex-1)*ty + 1 + (j-j1+lb)
      sdatay(id      ) = ml(j,i1+(iex-1),k,n)
      sdatay(id+sizey) = ml(j,i2-(iex-1),k,n)
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdatay,rdatay)
    call mpi_neighbor_alltoallv(sdatay, counts, displs, TK_MPI, &
        rdatay, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif
    do concurrent ( j = j1-lb:j2+rb, iex = 1:nex, k = k1:k2, n = n1:n2 )
      id = (((n-n1)*nk+(k-k1))*nex+iex-1)*ty + 1 + (j-j1+lb)
      if ( do_bt ) ml(j,i1-iex,k,n) = rdatay(id      )
      if ( do_tp ) ml(j,i2+iex,k,n) = rdatay(id+sizey)
    end do
    !$acc end data
    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_4D_EXCHANGE

  subroutine TK_2D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP(ml,nex,j1,j2,i1,i2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: ndx, ndy, nx, ny, tx, ty, sizex, sizey
    !@acc call nvtxStartRange(TK_NAME//"_2d_exchange_left_right_bottom_top")
    nx = j2-j1+1
    ny = i2-i1+1
    tx = ny
    ty = nx
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx+ndy) :: sdata
    TK_TYPE, dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: iex, i, j, id, offset

    offset = 0
    do concurrent ( i = i1:i2, iex = 1:nex )
      id = (iex-1)*tx + 1 + (i-i1) + offset
      sdata(id) = ml(j1+(iex-1),i)
    end do
    offset = offset + sizex
    do concurrent ( i = i1:i2, iex = 1:nex )
      id = (iex-1)*tx + 1 + (i-i1) + offset
      sdata(id) = ml(j2-(iex-1),i)
    end do
    offset = offset + sizex
    do concurrent ( j = j1:j2, iex = 1:nex )
      id = (iex-1)*ty + 1 + (j-j1) + offset
      sdata(id) = ml(j,i1+(iex-1))
    end do
    offset = offset + sizey
    do concurrent ( j = j1:j2, iex = 1:nex )
      id = (iex-1)*ty + 1 + (j-j1) + offset
      sdata(id) = ml(j,i2-(iex-1))
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, TK_MPI, &
        rdata, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    offset = 0
    if ( ma%left /= mpi_proc_null ) then
      do concurrent ( i = i1:i2, iex = 1:nex )
        id = (iex-1)*tx + 1 + (i-i1) + offset
        ml(j1-iex,i) = rdata(id)
      end do
    end if
    offset = offset + sizex
    if ( ma%right /= mpi_proc_null ) then
      do concurrent ( i = i1:i2, iex = 1:nex )
        id = (iex-1)*tx + 1 + (i-i1) + offset
        ml(j2+iex,i) = rdata(id)
      end do
    end if
    offset = offset + sizex
    if ( ma%bottom /= mpi_proc_null ) then
      do concurrent ( j = j1:j2, iex = 1:nex )
        id = (iex-1)*ty + 1 + (j-j1) + offset
        ml(j,i1-iex) = rdata(id)
      end do
    end if
    offset = offset + sizey
    if ( ma%top /= mpi_proc_null ) then
      do concurrent ( j = j1:j2, iex = 1:nex )
        id = (iex-1)*ty + 1 + (j-j1) + offset
        ml(j,i2+iex) = rdata(id)
      end do
    end if

    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_2D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP

  subroutine TK_3D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP(ml,nex,j1,j2,i1,i2,k1,k2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: ndx, ndy, nx, ny, nk, tx, ty, sizex, sizey
    logical :: do_lt, do_rt, do_bt, do_tp
    !@acc call nvtxStartRange(TK_NAME//"_3d_exchange_left_right_bottom_top")
    do_lt = .true.
    do_rt = .true.
    do_bt = .true.
    do_tp = .true.
    if ( ma%left   == mpi_proc_null ) do_lt = .false.
    if ( ma%right  == mpi_proc_null ) do_rt = .false.
    if ( ma%bottom == mpi_proc_null ) do_bt = .false.
    if ( ma%top    == mpi_proc_null ) do_tp = .false.
    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx+ndy) :: sdata
    TK_TYPE, dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: id, iex, i, j, k

    !$acc data copy(ml) create(sdata,rdata)
    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*tx + 1 + (i-i1)
      sdata(id      ) = ml(j1+(iex-1),i,k)
      sdata(id+sizex) = ml(j2-(iex-1),i,k)
    end do
    do concurrent ( j = j1:j2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*ty + 1 + (j-j1) + sizex + sizex
      sdata(id      ) = ml(j,i1+(iex-1),k)
      sdata(id+sizey) = ml(j,i2-(iex-1),k)
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    !$acc host_data use_device(sdata,rdata)
    call mpi_neighbor_alltoallv(sdata, counts, displs, TK_MPI, &
        rdata, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif
    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*tx + 1 + (i-i1)
      if ( do_lt ) ml(j1-iex,i,k) = rdata(id      )
      if ( do_rt ) ml(j2+iex,i,k) = rdata(id+sizex)
    end do
    do concurrent ( j = j1:j2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*ty + 1 + (j-j1) + sizex + sizex
      if ( do_bt ) ml(j,i1-iex,k) = rdata(id      )
      if ( do_tp ) ml(j,i2+iex,k) = rdata(id+sizey)
    end do
    !$acc end data
    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_3D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP

  subroutine TK_4D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP(ml,nex,j1,j2, &
                                                     i1,i2,k1,k2,n1,n2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: ndx, ndy, nx, ny, nk, nn, tx, ty, sizex, sizey
    logical :: do_lt, do_rt, do_bt, do_tp
    !@acc call nvtxStartRange(TK_NAME//"_4d_exchange_left_right_bottom_top")
    do_lt = .true.
    do_rt = .true.
    do_bt = .true.
    do_tp = .true.
    if ( ma%left   == mpi_proc_null ) do_lt = .false.
    if ( ma%right  == mpi_proc_null ) do_rt = .false.
    if ( ma%bottom == mpi_proc_null ) do_bt = .false.
    if ( ma%top    == mpi_proc_null ) do_tp = .false.
    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx+ndy) :: sdata
    TK_TYPE, dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: id, iex, i, j, k, n

    !$acc data copy(ml) create(sdata,rdata)
    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2, n = n1:n2 )
      id = (((n-n1)*nk+(k-k1))*nex+iex-1)*tx + 1 + (i-i1)
      sdata(id      ) = ml(j1+(iex-1),i,k,n)
      sdata(id+sizex) = ml(j2-(iex-1),i,k,n)
    end do
    do concurrent ( j = j1:j2, iex = 1:nex, k = k1:k2, n = n1:n2 )
      id = (((n-n1)*nk+(k-k1))*nex+iex-1)*ty + 1 + (j-j1) + &
           sizex + sizex
      sdata(id      ) = ml(j,i1+(iex-1),k,n)
      sdata(id+sizey) = ml(j,i2-(iex-1),k,n)
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    !$acc host_data use_device(sdata,rdata)
    call mpi_neighbor_alltoallv(sdata, counts, displs, TK_MPI, &
        rdata, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2, n = n1:n2 )
      id = (((n-n1)*nk+(k-k1))*nex+iex-1)*tx + 1 + (i-i1)
      if ( do_lt ) ml(j1-iex,i,k,n) = rdata(id      )
      if ( do_rt ) ml(j2+iex,i,k,n) = rdata(id+sizex)
    end do
    do concurrent ( j = j1:j2, iex = 1:nex, k = k1:k2, n = n1:n2 )
      id = (((n-n1)*nk+(k-k1))*nex+iex-1)*ty + 1 + (j-j1) + &
           sizex + sizex
      if ( do_bt ) ml(j,i1-iex,k,n) = rdata(id      )
      if ( do_tp ) ml(j,i2+iex,k,n) = rdata(id+sizey)
    end do
    !$acc end data

    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_4D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP

  subroutine TK_2D_EXCHANGE_LEFT_RIGHT(ml,nex,j1,j2,i1,i2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: ndx, ny, tx, sizex
    !@acc call nvtxStartRange(TK_NAME//"_2d_exchange_left_right")
    ny = i2-i1+1
    tx = ny
    sizex = nex*tx
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdata
    TK_TYPE, dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, TK_MPI, &
        rdata, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_2D_EXCHANGE_LEFT_RIGHT

  subroutine TK_3D_EXCHANGE_LEFT_RIGHT(ml,nex,j1,j2,i1,i2,k1,k2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: ndx, ny, nk, tx, sizex
    logical :: do_lt, do_rt
    !@acc call nvtxStartRange(TK_NAME//"_3d_exchange_left_right")
    do_lt = .true.
    do_rt = .true.
    if ( ma%left   == mpi_proc_null ) do_lt = .false.
    if ( ma%right  == mpi_proc_null ) do_rt = .false.
    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    sizex = nex*tx*nk
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdata
    TK_TYPE, dimension(ndx), volatile :: rdata
    integer(ik4) :: id, iex, i, k

    !$acc data copy(ml) create(sdata,rdata)
    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*tx + 1 + (i-i1)
      sdata(id      ) = ml(j1+(iex-1),i,k)
      sdata(id+sizex) = ml(j2-(iex-1),i,k)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdata,rdata)
    call mpi_neighbor_alltoallv(sdata, counts, displs, TK_MPI, &
        rdata, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    do concurrent ( i = i1:i2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*tx + 1 + (i-i1)
      if ( do_lt ) ml(j1-iex,i,k) = rdata(id      )
      if ( do_rt ) ml(j2+iex,i,k) = rdata(id+sizex)
    end do
    !$acc end data
    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_3D_EXCHANGE_LEFT_RIGHT

  subroutine TK_4D_EXCHANGE_LEFT_RIGHT(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: ndx, ny, nk, nn, tx, sizex
    !@acc call nvtxStartRange(TK_NAME//"_4d_exchange_left_right")
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    sizex = nex*tx*nk*nn
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdata
    TK_TYPE, dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, TK_MPI, &
        rdata, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_4D_EXCHANGE_LEFT_RIGHT

  subroutine TK_2D_EXCHANGE_BOTTOM_TOP(ml,nex,j1,j2,i1,i2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: ndy, nx, ty, sizey
    !@acc call nvtxStartRange(TK_NAME//"_2d_exchange_bottom_top")
    nx = j2-j1+1
    ty = nx
    sizey = nex*ty
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndy) :: sdata
    TK_TYPE, dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, TK_MPI, &
        rdata, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_2D_EXCHANGE_BOTTOM_TOP

  subroutine TK_3D_EXCHANGE_BOTTOM_TOP(ml,nex,j1,j2,i1,i2,k1,k2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: ndy, nx, nk, ty, sizey
    logical :: do_bt, do_tp
    !@acc call nvtxStartRange(TK_NAME//"_3d_exchange_bottom_top")
    do_bt = .true.
    do_tp = .true.
    if ( ma%bottom == mpi_proc_null ) do_bt = .false.
    if ( ma%top    == mpi_proc_null ) do_tp = .false.
    nx = j2-j1+1
    nk = k2-k1+1
    ty = nx
    sizey = nex*ty*nk
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndy) :: sdata
    TK_TYPE, dimension(ndy), volatile :: rdata
    integer(ik4) :: id, iex, j, k

    !$acc data copy(ml) create(sdata,rdata)
    do concurrent ( j = j1:j2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*ty + 1 + (j-j1)
      sdata(id      ) = ml(j,i1+(iex-1),k)
      sdata(id+sizey) = ml(j,i2-(iex-1),k)
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdata,rdata)
    call mpi_neighbor_alltoallv(sdata, counts, displs, TK_MPI, &
        rdata, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif
    do concurrent ( j = j1:j2, iex = 1:nex, k = k1:k2 )
      id = ((k-k1)*nex+iex-1)*ty + 1 + (j-j1)
      if ( do_bt ) ml(j,i1-iex,k) = rdata(id      )
      if ( do_tp ) ml(j,i2+iex,k) = rdata(id+sizey)
    end do
    !$acc end data
    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_3D_EXCHANGE_BOTTOM_TOP

  subroutine TK_4D_EXCHANGE_BOTTOM_TOP(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    !@acc use nvtx
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: ndy, nx, nk, nn, ty, sizey
    !@acc call nvtxStartRange(TK_NAME//"_4d_exchange_bottom_top")
    nx = j2-j1+1
    nk = k2-k1+1
    nn = n2-n1+1
    ty = nx
    sizey = nex*ty*nk*nn
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndy) :: sdata
    TK_TYPE, dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, TK_MPI, &
        rdata, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange
    !@acc call nvtxEndRange
  end subroutine TK_4D_EXCHANGE_BOTTOM_TOP

  subroutine TK_2D_EXCHANGE_LEFT_BOTTOM(ml,nex,j1,j2,i1,i2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: nx, ny
    integer(ik4) :: lb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdatax
    TK_TYPE, dimension(ndy) :: sdatay
    TK_TYPE, dimension(ndx), volatile :: rdatax
    TK_TYPE, dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex

    !$acc data copy(ml) create(sdatax,sdatay,rdatax,rdatay)

    do concurrent ( iex = 1:nex )
      ib1 = 1 + (iex-1)*(tx)
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do concurrent ( iex = 1:nex )
      ib1 = 1 + (iex-1)*(tx) + sizex
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdatax,rdatax)
    call mpi_neighbor_alltoallv(sdatax, counts, displs, TK_MPI, &
        rdatax, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%left /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex )
        ib1 = 1 + (iex-1)*(tx)
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdatax(ib1:ib2)
      end do
    end if

    do concurrent ( iex = 1:nex )
      ib1 = 1 + (iex-1)*(ty)
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1))
    end do
    do concurrent ( iex = 1:nex )
      ib1 = 1 + (iex-1)*(ty) + sizey
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdatay,rdatay)
    call mpi_neighbor_alltoallv(sdatay, counts, displs, TK_MPI, &
        rdatay, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%bottom /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex )
        ib1 = 1 + (iex-1)*(ty)
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2,i1-iex) = rdatay(ib1:ib2)
      end do
    end if
    !$acc end data
    end block exchange

  end subroutine TK_2D_EXCHANGE_LEFT_BOTTOM

  subroutine TK_3D_EXCHANGE_LEFT_BOTTOM(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: nx, ny, nk
    integer(ik4) :: lb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdatax
    TK_TYPE, dimension(ndy) :: sdatay
    TK_TYPE, dimension(ndx), volatile :: rdatax
    TK_TYPE, dimension(ndy), volatile :: rdatay
    integer(ik4) :: id, ib1, ib2, iex, k

    !$acc data copy(ml) create(sdatax,sdatay,rdatax,rdatay)

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx)
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx) + sizex
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdatax,rdatax)
    call mpi_neighbor_alltoallv(sdatax, counts, displs, TK_MPI, &
        rdatax, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%left /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(tx)
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2,k) = rdatax(ib1:ib2)
      end do
    end if

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty)
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty) + sizey
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k)
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdatay,rdatay)
    call mpi_neighbor_alltoallv(sdatay, counts, displs, TK_MPI, &
        rdatay, counts, displs, TK_MPI, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%bottom /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty)
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2,i1-iex,k) = rdatay(ib1:ib2)
      end do
    end if
    !$acc end data
    end block exchange

  end subroutine TK_3D_EXCHANGE_LEFT_BOTTOM

  subroutine TK_4D_EXCHANGE_LEFT_BOTTOM(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: nx, ny, nk, nn
    integer(ik4) :: lb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdatax
    TK_TYPE, dimension(ndy) :: sdatay
    TK_TYPE, dimension(ndx), volatile :: rdatax
    TK_TYPE, dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, TK_MPI, &
        rdatax, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, TK_MPI, &
        rdatay, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2,i1-iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine TK_4D_EXCHANGE_LEFT_BOTTOM

  subroutine TK_2D_EXCHANGE_RIGHT_TOP(ml,nex,j1,j2,i1,i2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: nx, ny
    integer(ik4) :: rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    rb = nex
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdatax
    TK_TYPE, dimension(ndy) :: sdatay
    TK_TYPE, dimension(ndx), volatile :: rdatax
    TK_TYPE, dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, TK_MPI, &
        rdatax, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizex
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, TK_MPI, &
        rdatay, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizey
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2+rb,i2+iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine TK_2D_EXCHANGE_RIGHT_TOP

  subroutine TK_3D_EXCHANGE_RIGHT_TOP(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: nx, ny, nk
    integer(ik4) :: rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    rb = nex
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdatax
    TK_TYPE, dimension(ndy) :: sdatay
    TK_TYPE, dimension(ndx), volatile :: rdatax
    TK_TYPE, dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, TK_MPI, &
        rdatax, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizex
    if ( ma%right /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, TK_MPI, &
        rdatay, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizey
    if ( ma%top /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2+rb,i2+iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine TK_3D_EXCHANGE_RIGHT_TOP

  subroutine TK_4D_EXCHANGE_RIGHT_TOP(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: nx, ny, nk, nn
    integer(ik4) :: rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    rb = nex
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    TK_TYPE, dimension(ndx) :: sdatax
    TK_TYPE, dimension(ndy) :: sdatay
    TK_TYPE, dimension(ndx), volatile :: rdatax
    TK_TYPE, dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, TK_MPI, &
        rdatax, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizex
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, TK_MPI, &
        rdatay, counts, displs, TK_MPI, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizey
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2+rb,i2+iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine TK_4D_EXCHANGE_RIGHT_TOP

  subroutine TK_BDY_EXCHANGE_LEFT_RIGHT(ml,k1,k2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: k1, k2
    integer(ik4) :: ssize, ksize, k, ib, irc, ipos
    type(mpi_request), dimension(4) :: req
    req(:) = mpi_request_null
    ksize = k2-k1+1
    irc = 0
    ipos = 1
    ssize = 2*ksize
    ! Not allocated by set_nproc when running on one processor
    if ( .not. associated(TK_VEC1) ) then
      call getmem(TK_VEC1,1,ssize,TK_NAME//'_bdy_exchange_left_right')
    else if ( size(TK_VEC1) < ssize ) then
      call getmem(TK_VEC1,1,ssize,TK_NAME//'_bdy_exchange_left_right')
    end if
    if ( .not. associated(TK_VEC2) ) then
      call getmem(TK_VEC2,1,ssize,TK_NAME//'_bdy_exchange_left_right')
    else if ( size(TK_VEC2) < ssize ) then
      call getmem(TK_VEC2,1,ssize,TK_NAME//'_bdy_exchange_left_right')
    end if
    if ( ma%bandflag ) then
      ib = 1
      do k = k1, k2
        TK_VEC1(ib) = ml(jde2,k)
        ib = ib + 1
      end do
      call cyclic_exchange_array(TK_VEC1,TK_VEC2, &
                                 ksize,ma%left,ma%right,__LINE__)
      ib = 1
      do k = k1, k2
        ml(jde1-1,k) = TK_VEC2(ib)
        ib = ib + 1
      end do
      ib = 1
      do k = k1, k2
        TK_VEC1(ib) = ml(jde1,k)
        ib = ib + 1
      end do
      call cyclic_exchange_array(TK_VEC1,TK_VEC2, &
                                 ksize,ma%right,ma%left,__LINE__)
      ib = 1
      do k = k1, k2
        ml(jde2+1,k) = TK_VEC2(ib)
        ib = ib + 1
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        ib = ipos
        do k = k1, k2
          TK_VEC1(ib) = ml(jde2,k)
          ib = ib + 1
        end do
        irc = irc + 1
        call exchange_array(TK_VEC1(ipos:ipos+ksize-1), &
                            TK_VEC2(ipos:ipos+ksize-1),ksize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+2))
        ipos = ipos + ksize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ib = ipos
        do k = k1, k2
          TK_VEC1(ib) = ml(jde1,k)
          ib = ib + 1
        end do
        irc = irc + 1
        call exchange_array(TK_VEC1(ipos:ipos+ksize-1), &
                            TK_VEC2(ipos:ipos+ksize-1),ksize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+2))
        ipos = ipos + ksize
      end if
      if ( irc /= 0 ) then
        call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_waitall error.')
        end if
#endif
        ipos = 1
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          do k = k1, k2
            ml(jde2+1,k) = TK_VEC2(ib)
            ib = ib + 1
          end do
          ipos = ipos + ksize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          do k = k1, k2
            ml(jde1-1,k) = TK_VEC2(ib)
            ib = ib + 1
          end do
          ipos = ipos + ksize
        end if
      end if
    end if
  end subroutine TK_BDY_EXCHANGE_LEFT_RIGHT

  subroutine TK_BDY_EXCHANGE_BOTTOM_TOP(ml,k1,k2)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: k1, k2
    integer(ik4) :: ssize, ksize, k, ib, irc, ipos
    type(mpi_request), dimension(4) :: req
    req(:) = mpi_request_null
    ksize = k2-k1+1
    ssize = 2*ksize
    ! Not allocated by set_nproc when running on one processor
    if ( .not. associated(TK_VEC1) ) then
      call getmem(TK_VEC1,1,ssize,TK_NAME//'_bdy_exchange_bottom_top')
    else if ( size(TK_VEC1) < ssize ) then
      call getmem(TK_VEC1,1,ssize,TK_NAME//'_bdy_exchange_bottom_top')
    end if
    if ( .not. associated(TK_VEC2) ) then
      call getmem(TK_VEC2,1,ssize,TK_NAME//'_bdy_exchange_bottom_top')
    else if ( size(TK_VEC2) < ssize ) then
      call getmem(TK_VEC2,1,ssize,TK_NAME//'_bdy_exchange_bottom_top')
    end if
    irc = 0
    ipos = 1
    if ( ma%top /= mpi_proc_null) then
      ib = ipos
      do k = k1, k2
        TK_VEC1(ib) = ml(ide2,k)
        ib = ib + 1
      end do
      irc = irc + 1
      call exchange_array(TK_VEC1(ipos:ipos+ksize-1), &
                          TK_VEC2(ipos:ipos+ksize-1),ksize, &
                          ma%top,tag_tb,tag_bt, &
                          req(irc),req(irc+2))
      ipos = ipos + ksize
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      ib = ipos
      do k = k1, k2
        TK_VEC1(ib) = ml(ide1,k)
        ib = ib + 1
      end do
      irc = irc + 1
      call exchange_array(TK_VEC1(ipos:ipos+ksize-1), &
                          TK_VEC2(ipos:ipos+ksize-1),ksize, &
                          ma%bottom,tag_bt,tag_tb, &
                          req(irc),req(irc+2))
      ipos = ipos + ksize
    end if
    if ( irc /= 0 ) then
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      ipos = 1
      if ( ma%top /= mpi_proc_null) then
        ib = ipos
        do k = k1, k2
          ml(ide2+1,k) = TK_VEC2(ib)
          ib = ib + 1
        end do
        ipos = ipos + ksize
      end if
      if ( ma%bottom /= mpi_proc_null ) then
        ib = ipos
        do k = k1, k2
          ml(ide1-1,k) = TK_VEC2(ib)
          ib = ib + 1
        end do
        ipos = ipos + ksize
      end if
    end if
  end subroutine TK_BDY_EXCHANGE_BOTTOM_TOP

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

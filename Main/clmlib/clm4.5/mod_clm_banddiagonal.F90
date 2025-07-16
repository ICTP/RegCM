module mod_clm_banddiagonal
  !
  ! Band Diagonal matrix solution
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_stdio
  use lapack_dgbsv

  implicit none

  private

  save

  public :: BandDiagonal

  interface BandDiagonal
#if defined(STDPAR) || defined(OPENACC) || defined(_OPENACC)
    module procedure BandDiagonal_gpu
#else
    module procedure BandDiagonal_cpu
#endif
  end interface

  contains
  !
  ! Tridiagonal matrix solution
  !
  subroutine BandDiagonal_cpu(lbc, ubc, lbj, ubj, jtop, jbot, numf, &
                          filter, nband, b, r, u)
    implicit none
    ! lbinning and ubing column indices
    integer(ik4), intent(in)    :: lbc, ubc
    ! lbinning and ubing level indices
    integer(ik4), intent(in)    :: lbj, ubj
    ! top level for each column
    integer(ik4), intent(in)    :: jtop(lbc:ubc)
    !scs:  add jbot
    ! bottom level for each column
    integer(ik4), intent(in)    :: jbot(lbc:ubc)
    !scs
    integer(ik4), intent(in)    :: numf   ! filter dimension
    integer(ik4), intent(in)    :: nband  ! band width
    integer(ik4), intent(in)    :: filter(ubc-lbc+1)       ! filter
    ! compact band matrix
    real(rk8), intent(in)    :: b(lbc:ubc,nband,lbj:ubj)
    ! "r" rhs of linear system
    real(rk8), intent(in)    :: r(lbc:ubc, lbj:ubj)
    ! solution
    real(rk8), intent(inout) :: u(lbc:ubc, lbj:ubj)

    integer(ik4)  :: j,ci,fc,info,m,n          !indices
    integer(ik4)  :: kl,ku                     !number of sub/super diagonals
    integer(ik4), allocatable :: ipiv(:)       !temporary
    real(rk8),allocatable :: ab(:,:),temp(:,:) !compact storage array
    real(rk8),allocatable :: result(:)

!!$     SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
!!$*
!!$*  -- LAPACK driver routine (version 3.1) --
!!$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!!$*     November 2006
!!$*
!!$*     .. Scalar Arguments ..
!!$      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!!$*     ..
!!$*     .. Array Arguments ..
!!$      INTEGER            IPIV( * )
!!$      REAL               AB( LDAB, * ), B( LDB, * )
!!$*     ..
!!$*
!!$*  Purpose
!!$*  =======
!!$*
!!$*  DGBSV computes the solution to a real system of linear equations
!!$*  A * X = B, where A is a band matrix of order N with KL subdiagonals
!!$*  and KU superdiagonals, and X and B are N-by-NRHS matrices.
!!$*
!!$*  The LU decomposition with partial pivoting and row interchanges is
!!$*  used to factor A as A = L * U, where L is a product of permutation
!!$*  and unit lower triangular matrices with KL subdiagonals, and U is
!!$*  upper triangular with KL+KU superdiagonals.  The factored form of A
!!$*  is then used to solve the system of equations A * X = B.
!!$*
!!$*  Arguments
!!$*  =========
!!$*
!!$*  N       (input) INTEGER
!!$*          The number of linear equations, i.e., the order of the
!!$*          matrix A.  N >= 0.
!!$*
!!$*  KL      (input) INTEGER
!!$*          The number of subdiagonals within the band of A.  KL >= 0.
!!$*
!!$*  KU      (input) INTEGER
!!$*          The number of superdiagonals within the band of A.  KU >= 0.
!!$*
!!$*  NRHS    (input) INTEGER
!!$*          The number of right hand sides, i.e., the number of columns
!!$*          of the matrix B.  NRHS >= 0.
!!$*
!!$*  AB      (input/output) REAL array, dimension (LDAB,N)
!!$*          On entry, the matrix A in band storage, in rows KL+1 to
!!$*          2*KL+KU+1; rows 1 to KL of the array need not be set.
!!$*          The j-th column of A is stored in the j-th column of the
!!$*          array AB as follows:
!!$*          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
!!$*          On exit, details of the factorization: U is stored as an
!!$*          upper triangular band matrix with KL+KU superdiagonals in
!!$*          rows 1 to KL+KU+1, and the multipliers used during the
!!$*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!!$*          See below for further details.
!!$*
!!$*  LDAB    (input) INTEGER
!!$*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!!$*
!!$*  IPIV    (output) INTEGER array, dimension (N)
!!$*          The pivot indices that define the permutation matrix P;
!!$*          row i of the matrix was interchanged with row IPIV(i).
!!$*
!!$*  B       (input/output) REAL array, dimension (LDB,NRHS)
!!$*          On entry, the N-by-NRHS right hand side matrix B.
!!$*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!!$*
!!$*  LDB     (input) INTEGER
!!$*          The leading dimension of the array B.  LDB >= max(1,N).
!!$*
!!$*  INFO    (output) INTEGER
!!$*          = 0:  successful exit
!!$*          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!!$*                has been completed, but the factor U is exactly
!!$*                singular, and the solution has not been computed.
!!$*
!!$*  Further Details
!!$*  ===============
!!$*
!!$*  The band storage scheme is illustrated by the following example, when
!!$*  M = N = 6, KL = 2, KU = 1:
!!$*
!!$*  On entry:                       On exit:
!!$*
!!$*      *    *    *    +    +    +       *    *    *   u14  u25  u36
!!$*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!!$*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!!$*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!!$*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!!$*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!!$*
!!$*  Array elements marked * are not used by the routine; elements marked
!!$*  + need not be set on entry, but are required by the routine to store
!!$*  elements of U because of fill-in resulting from the row interchanges.

    !OPEN(unit=10, file='filter.dat', form='unformatted')        
    !write (10) numf,lbc, ubc,lbj, ubj,filter 
    !close(10)
!
!    OPEN(unit=10, file='jtop.dat', form='unformatted')        
!    write (10) numf,jtop 
!    close(10)
!
!    OPEN(unit=10, file='jbot.dat', form='unformatted')        
!    write (10) numf,jbot 
!    close(10)
!
!    OPEN(unit=10, file='b.dat', form='unformatted')        
!    write (10) b 
!    close(10)
!
!    OPEN(unit=10, file='r.dat', form='unformatted')        
!    write (10) r 
!    close(10)

    !Set up input matrix AB
    !An m-by-n band matrix with kl subdiagonals and ku superdiagonals
    !may be stored compactly in a two-dimensional array with
    !kl+ku+1 rows and n columns
    !AB(KL+KU+1+i-j,j) = A(i,j)

    !print *," In banddiagonal numf",numf
    do fc = 1, numf
      ci = filter(fc)

      kl = (nband-1)/2
      ku = kl
      ! m is the number of rows required for storage space by dgbsv
      m = 2*kl+ku+1
      ! n is the number of levels (snow/soil)
      !scs: replace ubj with jbot
      n = jbot(ci)-jtop(ci)+1
      !print *," In banddiagonal fc:  m,n,ci",fc,m,n,ci

      allocate(ab(m,n))
      ab = 0.0_rk8

      ab(kl+ku-1,3:n) = b(ci,1,jtop(ci):jbot(ci)-2)   ! 2nd superdiagonal
      ab(kl+ku+0,2:n) = b(ci,2,jtop(ci):jbot(ci)-1)   ! 1st superdiagonal
      ab(kl+ku+1,1:n) = b(ci,3,jtop(ci):jbot(ci))     ! diagonal
      ab(kl+ku+2,1:n-1) = b(ci,4,jtop(ci)+1:jbot(ci)) ! 1st subdiagonal
      ab(kl+ku+3,1:n-2) = b(ci,5,jtop(ci)+2:jbot(ci)) ! 2nd subdiagonal

      allocate(temp(m,n))
      temp = ab

      allocate(ipiv(n))
      allocate(result(n))

      ! on input result is rhs, on output result is solution vector
      result(:) = r(ci,jtop(ci):jbot(ci))

      ! DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
      call dgbsv( n, kl, ku, 1, ab, m, ipiv, result, n, info )
      u(ci,jtop(ci):jbot(ci))=result(:)

      if ( info /= 0 ) then
        write(stdout,*)'index: ', ci
        write(stdout,*)'n,kl,ku,m ',n,kl,ku,m
        write(stdout,*)'dgbsv info: ',ci,info
        write(stdout,*) ' '
        write(stdout,*) 'ab matrix'
        do j = 1, n
          ! write(stdout,'(i2,7f18.7)') j,temp(:,j)
          write(stdout,'(i2,5f18.7)') j,temp(3:7,j)
        end do
        write(stdout,*) ' '
        call fatal(__FILE__,__LINE__,'Linear algebra error')
      end if
      deallocate(temp)
      deallocate(ab)
      deallocate(ipiv)
      deallocate(result)
    end do
!    OPEN(unit=10, file='u.dat', form='unformatted')        
!    write (10) u 
!    close(10)
!    stop
  end subroutine BandDiagonal_cpu


#if defined(STDPAR) || defined(OPENACC) || defined(_OPENACC)

subroutine BandDiagonal_gpu(lbc, ubc, lbj, ubj, jtop, jbot, numf, &
                           filter, nband, b, r, u)
    use cusparse
    use cudafor
    use nvtx
    implicit none
    ! lbinning and ubing column indices
    integer(ik4), intent(in)    :: lbc, ubc
    ! lbinning and ubing level indices
    integer(ik4), intent(in)    :: lbj, ubj
    ! top level for each column
    integer(ik4), intent(in)    :: jtop(lbc:ubc)
    !scs:  add jbot
    ! bottom level for each column
    integer(ik4), intent(in)    :: jbot(lbc:ubc)
    !scs
    integer(ik4), intent(in)    :: numf   ! filter dimension
    integer(ik4), intent(in)    :: nband  ! band width
    integer(ik4), intent(in)    :: filter(ubc-lbc+1)       ! filter
    ! compact band matrix
    real(rk8), intent(in)    :: b(lbc:ubc,nband,lbj:ubj)
    ! "r" rhs of linear system
    real(rk8), intent(in)    :: r(lbc:ubc, lbj:ubj)
    ! solution
    real(rk8), intent(inout) :: u(lbc:ubc, lbj:ubj)

    integer(ik4)  :: j,ci,fc,info,m,n,i       !indices
    integer(ik4)  :: kl,ku                     !number of sub/super diagonals
    integer(ik4), allocatable :: cnt(:)        !array for histogram
    integer(ik4), allocatable :: ind(:)        !array for indices
    integer(ik4) :: nSets                      !How many 
    integer :: maxbuffer,maxdim           ! Values to dimension device arrays
    logical, allocatable :: mask(:)            !mask used to compact array
    integer(ik4), dimension(:), allocatable :: indices, sizes,locind

    ! cuSPARSE handles and buffers
    type(cusparseHandle) :: cusparse_h
    integer :: algo = 0  ! QR factorization (only one supported at this time)
    integer :: k,stat, batchSize
    integer(ik8) :: bufferSize
    real(rk8),  allocatable :: d_buffer(:)
    logical, save :: first_call=.true.

    ! Matrix diagonals and rhs on GPU
    real(rk8),allocatable :: d_S(:) ! 2nd subdiagonal
    real(rk8),allocatable :: d_L(:) ! 1st subdiagonal
    real(rk8),allocatable :: d_M(:) ! diagonal
    real(rk8),allocatable :: d_U(:) ! 1st superdiagonal
    real(rk8),allocatable :: d_W(:) ! 2nd superdiagonal
    real(rk8),allocatable :: d_X(:) ! RHS in input, solution in output

    call nvtxStartRange("BandDiagonal_gpu")

    ! Initialize cuSparse ( this can be done outside or just at the first call)
    if (first_call) then
      stat = cusparseCreate(cusparse_h)
      if (stat /= CUSPARSE_STATUS_SUCCESS) then
        write(*,*) 'ERROR: cuSPARSE initialization failed'
        stop
      end if
      first_call=.false.
    end if

    ! Find the set of possible values for N in the pentadiagonal solver
    ! building an histagram of the lenghts.
    ! They will be the positions in the cnt array with non-zero value.
    ! The non-zero value at each position will indicate how many pentadiagonal
    ! system of that size we can batch together
    allocate (cnt(numf))
    allocate (ind(numf))
    !allocate (sizes(numf))
    !allocate (locind(numf))

    !$acc kernels
    cnt=0
    !$acc end kernels

    !$acc parallel loop copyout(cnt)
    do fc = 1, numf
      ci = filter(fc)
      n = jbot(ci)-jtop(ci)+1
      ! store n for each position
      ind(fc)=n
      ! count histogram of values in input array
      !$acc atomic update
      cnt(n)=cnt(n)+1
    end do

    ! Find all the non-zero elements of the histogram.
    ! The positions m  will be the size of a batch cnt(m)  of pentadiagonal system
    ! Standard F90:
     indices = Merge( 0, [ ( i, i = 1, Size( cnt ) ) ], cnt == 0 )
     sizes = Pack( indices, indices /= 0 )   ! will contain the indices of the non zero cnt element
     nSets=size(sizes)

    ! HPC SDK has a special intrinsic, packloc, to find the locations in an array corresponding to particular
    ! value. It will also return the count

    !!$acc host_data use_device(sizes,cnt)
    !sizes = packloc( cnt .gt. 0, count=nSets )   ! will contain the indices of the non zero cnt element
    !!$acc end host_data

    ! Allocate device buffers large enough to hold all the required quantities
    ! for cusparseDgpsvInterleavedBatch ( both batches of diagonals, rhs and auxiliary buffer)

    maxbuffer=0
    maxdim=0
    do i=1,nSets
       n=sizes(i)
       batchSize=cnt(n)
       maxdim=max(maxdim,n*batchSize)
       stat = cusparseDgpsvInterleavedBatch_bufferSizeExt(cusparse_h, algo, n, &
                                                          d_S, d_L, d_M, &
                                                          d_U, d_W, d_X, &
                                                          batchSize, bufferSize)
       maxbuffer=max(maxbuffer,bufferSize)
    end do
     allocate(d_buffer(maxbuffer))
     ! We are allocating these buffers as 1D array
     allocate(d_W(maxdim),d_U(maxdim),d_M(maxdim))
     allocate(d_L(maxdim),d_S(maxdim),d_X(maxdim))

    !$acc data create(d_buffer,d_W,d_U,d_M,d_L,d_S,d_X)

    ! Loop over the sets of pentadiagonal systems batching all the ones with same dimension
    ! starting from the larger one 
    do i=1,nSets
       ! find out current n and batchsize
       n=sizes(i)
       ! Standard F90 on filter
       ! In the gather/scatter: ci=locind(k)
        batchSize=cnt(n)
        mask = ( ind ==n ) 
        locind=pack(filter,mask)

       ! Using packloc ( but we are processin ind instead of filter, so we will need an extra
       ! indirection for ci later on: ci=filter(locind(k))

       !!$acc host_data use_device(locind,ind)
       !locind = packloc( ind .eq. n, count=batchSize )   ! will contain the indices of the 
       !!$acc end host_data

       ! Pack the coefficients for cuSparse, expecting arrays(batchsize,n)
       ! Since they are defined as 1D, we will do the indexing as
       ! arr(k,j)=arr(k+(j-1)*batchSize)

       !$acc parallel loop collapse(2)
       do j=1,n
       do k=1,batchSize
        ci = locind(k)
        !ci = filter(locind(k))
        d_W (k+(j-1)*batchSize) = b(ci,1,jtop(ci)+j-1)   ! 2nd superdiagonal
        d_U (k+(j-1)*batchSize) = b(ci,2,jtop(ci)+j-1)   ! 1st superdiagonal
        d_M (k+(j-1)*batchSize) = b(ci,3,jtop(ci)+j-1)     ! diagonal
        d_L (k+(j-1)*batchSize) = b(ci,4,jtop(ci)+j-1) ! 1st subdiagonal
        d_S (k+(j-1)*batchSize) = b(ci,5,jtop(ci)+j-1) ! 2nd subdiagonal
        d_X (k+(j-1)*batchSize) = r(ci,  jtop(ci)+j-1)
       end do
       end do

       ! Solve the batch
       stat = cusparseDgpsvInterleavedBatch(cusparse_h, algo, n, &
              d_S, d_L, d_M, d_U, d_W, d_X, batchSize, d_buffer)

       ! Copy the solution back to the proper index

       !$acc parallel loop collapse(2)
       do j=1,n
       do k=1,batchSize
        ci = locind(k)
        u(ci,jtop(ci)+j-1) = d_X(k+(j-1)*batchSize)
        end do
       end do

    end do

       !Deallocate buffers
        deallocate(d_W,d_U,d_M)
        deallocate(d_L,d_S,d_X)
        deallocate(d_buffer)
        deallocate (cnt)
        deallocate (ind)

     !$acc end data

     ! We are creating the handle the first time the routine is called
     ! No need to destroy it
     ! stat = cusparseDestroy(cusparse_h)

     call nvtxEndRange()

    end subroutine bandDiagonal_gpu

#endif

end module mod_clm_banddiagonal
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

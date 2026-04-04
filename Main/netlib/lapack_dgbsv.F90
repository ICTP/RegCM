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

module lapack_dgbsv

  use mod_realkinds, only: rk8

  implicit none

  private

  public :: dgbsv

contains

!> \brief <b> DGBSV computes the solution to system of linear equations A * X = B for GB matrices</b> (simple driver)
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       real(rk8)   AB( LDAB, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBSV computes the solution to a real system of linear equations
!> A * X = B, where A is a band matrix of order N with KL subdiagonals
!> and KU superdiagonals, and X and B are N-by-NRHS matrices.
!>
!> The LU decomposition with partial pivoting and row interchanges is
!> used to factor A as A = L * U, where L is a product of permutation
!> and unit lower triangular matrices with KL subdiagonals, and U is
!> upper triangular with KL+KU superdiagonals.  The factored form of A
!> is then used to solve the system of equations A * X = B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of linear equations, i.e., the order of the
!>          matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is real(rk8) array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows KL+1 to
!>          2*KL+KU+1; rows 1 to KL of the array need not be set.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
!>          On exit, details of the factorization: U is stored as an
!>          upper triangular band matrix with KL+KU superdiagonals in
!>          rows 1 to KL+KU+1, and the multipliers used during the
!>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!>          See below for further details.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices that define the permutation matrix P;
!>          row i of the matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is real(rk8) array, dimension (LDB,NRHS)
!>          On entry, the N-by-NRHS right hand side matrix B.
!>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!>                has been completed, but the factor U is exactly
!>                singular, and the solution has not been computed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup doubleGBsolve
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  M = N = 6, KL = 2, KU = 1:
!>
!>  On entry:                       On exit:
!>
!>      *    *    *    +    +    +       *    *    *   u14  u25  u36
!>      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!>
!>  Array elements marked * are not used by the routine; elements marked
!>  + need not be set on entry, but are required by the routine to store
!>  elements of U because of fill-in resulting from the row interchanges.
!> \endverbatim
!>
!  =====================================================================
  subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    implicit none
!
!  -- LAPACK driver routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    integer :: info, kl, ku, ldab, ldb, n, nrhs
!     ..
!     .. Array Arguments ..
    integer :: ipiv(*)
    real (rk8) :: ab(ldab, *), b(ldb, *)
!     ..
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
    intrinsic :: max
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
    info = 0
    if (n<0) then
      info = -1
    else if (kl<0) then
      info = -2
    else if (ku<0) then
      info = -3
    else if (nrhs<0) then
      info = -4
    else if (ldab<2*kl+ku+1) then
      info = -6
    else if (ldb<max(n,1)) then
      info = -9
    end if
    if (info/=0) then
      call xerbla('DGBSV ', -info)
      return
    end if
!
!     Compute the LU factorization of the band matrix A.
!
    call dgbtrf(n, n, kl, ku, ab, ldab, ipiv, info)
    if (info==0) then
!
!        Solve the system A*X = B, overwriting B with X.
!
      call dgbtrs('No transpose', n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, &
        info)
    end if
    return
!
!     End of DGBSV
!
  end subroutine dgbsv

!> \brief \b DGBTF2 computes the LU factorization of a general band matrix using the unblocked version of the algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBTF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDAB, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       real(rk8)   AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBTF2 computes an LU factorization of a real m-by-n band matrix A
!> using partial pivoting with row interchanges.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is real(rk8) array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows KL+1 to
!>          2*KL+KU+1; rows 1 to KL of the array need not be set.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!>
!>          On exit, details of the factorization: U is stored as an
!>          upper triangular band matrix with KL+KU superdiagonals in
!>          rows 1 to KL+KU+1, and the multipliers used during the
!>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!>          See below for further details.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!>               has been completed, but the factor U is exactly
!>               singular, and division by zero will occur if it is used
!>               to solve a system of equations.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date September 2012
!
!> \ingroup doubleGBcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  M = N = 6, KL = 2, KU = 1:
!>
!>  On entry:                       On exit:
!>
!>      *    *    *    +    +    +       *    *    *   u14  u25  u36
!>      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!>
!>  Array elements marked * are not used by the routine; elements marked
!>  + need not be set on entry, but are required by the routine to store
!>  elements of U, because of fill-in resulting from the row
!>  interchanges.
!> \endverbatim
!>
!  =====================================================================
  subroutine dgbtf2(m, n, kl, ku, ab, ldab, ipiv, info)
    implicit none
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
    integer :: info, kl, ku, ldab, m, n
!     ..
!     .. Array Arguments ..
    integer :: ipiv(*)
    real (rk8) :: ab(ldab, *)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    real (rk8) :: one, zero
    parameter (one=1.0_rk8, zero=0.0_rk8)
!     ..
!     .. Local Scalars ..
    integer :: i, j, jp, ju, km, kv
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: max, min
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
    kv = ku + kl
!
!     Test the input parameters.
!
    info = 0
    if (m<0) then
      info = -1
    else if (n<0) then
      info = -2
    else if (kl<0) then
      info = -3
    else if (ku<0) then
      info = -4
    else if (ldab<kl+kv+1) then
      info = -6
    end if
    if (info/=0) then
      call xerbla('DGBTF2', -info)
      return
    end if
!
!     Quick return if possible
!
    if (m==0 .or. n==0) return
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
    do j = ku + 2, min(kv, n)
      do i = kv - j + 2, kl
        ab(i, j) = zero
      end do
    end do
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
    ju = 1
!
    do j = 1, min(m, n)
!
!        Set fill-in elements in column J+KV to zero.
!
      if (j+kv<=n) then
        do i = 1, kl
          ab(i, j+kv) = zero
        end do
      end if
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
      km = min(kl, m-j)
      jp = idamax(km+1, ab(kv+1,j), 1)
      ipiv(j) = jp + j - 1
      if (ab(kv+jp,j)/=zero) then
        ju = max(ju, min(j+ku+jp-1,n))
!
!           Apply interchange to columns J to JU.
!
        if (jp/=1) call dswap(ju-j+1, ab(kv+jp,j), ldab-1, ab(kv+1,j), ldab-1)
!
        if (km>0) then
!
!              Compute multipliers.
!
          call dscal(km, one/ab(kv+1,j), ab(kv+2,j), 1)
!
!              Update trailing submatrix within the band.
!
          if (ju>j) call dger(km, ju-j, -one, ab(kv+2,j), 1, ab(kv,j+1), &
            ldab-1, ab(kv+1,j+1), ldab-1)
        end if
      else
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
        if (info==0) info = j
      end if
    end do
    return
!
!     End of DGBTF2
!
  end subroutine dgbtf2
!> \brief \b DGBTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDAB, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       real(rk8)   AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBTRF computes an LU factorization of a real m-by-n band matrix A
!> using partial pivoting with row interchanges.
!>
!> This is the blocked version of the algorithm, calling Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is real(rk8) array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows KL+1 to
!>          2*KL+KU+1; rows 1 to KL of the array need not be set.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!>
!>          On exit, details of the factorization: U is stored as an
!>          upper triangular band matrix with KL+KU superdiagonals in
!>          rows 1 to KL+KU+1, and the multipliers used during the
!>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!>          See below for further details.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!>               has been completed, but the factor U is exactly
!>               singular, and division by zero will occur if it is used
!>               to solve a system of equations.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup doubleGBcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  M = N = 6, KL = 2, KU = 1:
!>
!>  On entry:                       On exit:
!>
!>      *    *    *    +    +    +       *    *    *   u14  u25  u36
!>      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!>
!>  Array elements marked * are not used by the routine; elements marked
!>  + need not be set on entry, but are required by the routine to store
!>  elements of U because of fill-in resulting from the row interchanges.
!> \endverbatim
!>
!  =====================================================================
  subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
    implicit none
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    integer :: info, kl, ku, ldab, m, n
!     ..
!     .. Array Arguments ..
    integer :: ipiv(*)
    real (rk8) :: ab(ldab, *)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    real (rk8) :: one, zero
    parameter (one=1.0_rk8, zero=0.0_rk8)
    integer :: nbmax, ldwork
    parameter (nbmax=64, ldwork=nbmax+1)
!     ..
!     .. Local Scalars ..
    integer :: i, i2, i3, ii, ip, j, j2, j3, jb, jj, jm, jp, ju, k2, km, kv, &
      nb, nw
    real (rk8) :: temp
!     ..
!     .. Local Arrays ..
    real (rk8) :: work13(ldwork, nbmax), work31(ldwork, nbmax)
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: max, min
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
    kv = ku + kl
!
!     Test the input parameters.
!
    info = 0
    if (m<0) then
      info = -1
    else if (n<0) then
      info = -2
    else if (kl<0) then
      info = -3
    else if (ku<0) then
      info = -4
    else if (ldab<kl+kv+1) then
      info = -6
    end if
    if (info/=0) then
      call xerbla('DGBTRF', -info)
      return
    end if
!
!     Quick return if possible
!
    if (m==0 .or. n==0) return
!
!     Determine the block size for this environment
!
    nb = ilaenv(1, 'DGBTRF', m, n, kl, ku)
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
    nb = min(nb, nbmax)
!
    if (nb<=1 .or. nb>kl) then
!
!        Use unblocked code
!
      call dgbtf2(m, n, kl, ku, ab, ldab, ipiv, info)
    else
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
      do j = 1, nb
        do i = 1, j - 1
          work13(i, j) = zero
        end do
      end do
!
!        Zero the subdiagonal elements of the work array WORK31
!
      do j = 1, nb
        do i = j + 1, nb
          work31(i, j) = zero
        end do
      end do
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
      do j = ku + 2, min(kv, n)
        do i = kv - j + 2, kl
          ab(i, j) = zero
        end do
      end do
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
      ju = 1
!
      do j = 1, min(m, n), nb
        jb = min(nb, min(m,n)-j+1)
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
        i2 = min(kl-jb, m-j-jb+1)
        i3 = min(jb, m-j-kl+1)
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
        do jj = j, j + jb - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
          if (jj+kv<=n) then
            do i = 1, kl
              ab(i, jj+kv) = zero
            end do
          end if
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
          km = min(kl, m-jj)
          jp = idamax(km+1, ab(kv+1,jj), 1)
          ipiv(jj) = jp + jj - j
          if (ab(kv+jp,jj)/=zero) then
            ju = max(ju, min(jj+ku+jp-1,n))
            if (jp/=1) then
!
!                    Apply interchange to columns J to J+JB-1
!
              if (jp+jj-1<j+kl) then
!
                call dswap(jb, ab(kv+1+jj-j,j), ldab-1, ab(kv+jp+jj-j,j), &
                  ldab-1)
              else
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                call dswap(jj-j, ab(kv+1+jj-j,j), ldab-1, &
                  work31(jp+jj-j-kl,1), ldwork)
                call dswap(j+jb-jj, ab(kv+1,jj), ldab-1, ab(kv+jp,jj), ldab-1)
              end if
            end if
!
!                 Compute multipliers
!
            call dscal(km, one/ab(kv+1,jj), ab(kv+2,jj), 1)
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
            jm = min(ju, j+jb-1)
            if (jm>jj) call dger(km, jm-jj, -one, ab(kv+2,jj), 1, ab(kv,jj+1), &
              ldab-1, ab(kv+1,jj+1), ldab-1)
          else
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
            if (info==0) info = jj
          end if
!
!              Copy current column of A31 into the work array WORK31
!
          nw = min(jj-j+1, i3)
          if (nw>0) call dcopy(nw, ab(kv+kl+1-jj+j,jj), 1, work31(1,jj-j+1), &
            1)
        end do
        if (j+jb<=n) then
!
!              Apply the row interchanges to the other blocks.
!
          j2 = min(ju-j+1, kv) - jb
          j3 = max(0, ju-j-kv+1)
!
!              Use DLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
          call dlaswp(j2, ab(kv+1-jb,j+jb), ldab-1, 1, jb, ipiv(j), 1)
!
!              Adjust the pivot indices.
!
          do i = j, j + jb - 1
            ipiv(i) = ipiv(i) + j - 1
          end do
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
          k2 = j - 1 + jb + j2
          do i = 1, j3
            jj = k2 + i
            do ii = j + i - 1, j + jb - 1
              ip = ipiv(ii)
              if (ip/=ii) then
                temp = ab(kv+1+ii-jj, jj)
                ab(kv+1+ii-jj, jj) = ab(kv+1+ip-jj, jj)
                ab(kv+1+ip-jj, jj) = temp
              end if
            end do
          end do
!
!              Update the relevant part of the trailing submatrix
!
          if (j2>0) then
!
!                 Update A12
!
            call dtrsm('Left', 'Lower', 'No transpose', 'Unit', jb, j2, one, &
              ab(kv+1,j), ldab-1, ab(kv+1-jb,j+jb), ldab-1)
!
            if (i2>0) then
!
!                    Update A22
!
              call dgemm('No transpose', 'No transpose', i2, j2, jb, -one, &
                ab(kv+1+jb,j), ldab-1, ab(kv+1-jb,j+jb), ldab-1, one, &
                ab(kv+1,j+jb), ldab-1)
            end if
!
            if (i3>0) then
!
!                    Update A32
!
              call dgemm('No transpose', 'No transpose', i3, j2, jb, -one, &
                work31, ldwork, ab(kv+1-jb,j+jb), ldab-1, one, &
                ab(kv+kl+1-jb,j+jb), ldab-1)
            end if
          end if
!
          if (j3>0) then
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
            do jj = 1, j3
              do ii = jj, jb
                work13(ii, jj) = ab(ii-jj+1, jj+j+kv-1)
              end do
            end do
!
!                 Update A13 in the work array
!
            call dtrsm('Left', 'Lower', 'No transpose', 'Unit', jb, j3, one, &
              ab(kv+1,j), ldab-1, work13, ldwork)
!
            if (i2>0) then
!
!                    Update A23
!
              call dgemm('No transpose', 'No transpose', i2, j3, jb, -one, &
                ab(kv+1+jb,j), ldab-1, work13, ldwork, one, ab(1+jb,j+kv), &
                ldab-1)
            end if
!
            if (i3>0) then
!
!                    Update A33
!
              call dgemm('No transpose', 'No transpose', i3, j3, jb, -one, &
                work31, ldwork, work13, ldwork, one, ab(1+kl,j+kv), ldab-1)
            end if
!
!                 Copy the lower triangle of A13 back into place
!
            do jj = 1, j3
              do ii = jj, jb
                ab(ii-jj+1, jj+j+kv-1) = work13(ii, jj)
              end do
            end do
          end if
        else
!
!              Adjust the pivot indices.
!
          do i = j, j + jb - 1
            ipiv(i) = ipiv(i) + j - 1
          end do
        end if
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
        do jj = j + jb - 1, j, -1
          jp = ipiv(jj) - jj + 1
          if (jp/=1) then
!
!                 Apply interchange to columns J to JJ-1
!
            if (jp+jj-1<j+kl) then
!
!                    The interchange does not affect A31
!
              call dswap(jj-j, ab(kv+1+jj-j,j), ldab-1, ab(kv+jp+jj-j,j), &
                ldab-1)
            else
!
!                    The interchange does affect A31
!
              call dswap(jj-j, ab(kv+1+jj-j,j), ldab-1, work31(jp+jj-j-kl,1), &
                ldwork)
            end if
          end if
!
!              Copy the current column of A31 back into place
!
          nw = min(i3, jj-j+1)
          if (nw>0) call dcopy(nw, work31(1,jj-j+1), 1, ab(kv+kl+1-jj+j,jj), &
            1)
        end do
      end do
    end if
!
    return
!
!     End of DGBTRF
!
  end subroutine dgbtrf
!> \brief \b DGBTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       real(rk8)   AB( LDAB, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBTRS solves a system of linear equations
!>    A * X = B  or  A**T * X = B
!> with a general band matrix A using the LU factorization computed
!> by DGBTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations.
!>          = 'N':  A * X = B  (No transpose)
!>          = 'T':  A**T* X = B  (Transpose)
!>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is real(rk8) array, dimension (LDAB,N)
!>          Details of the LU factorization of the band matrix A, as
!>          computed by DGBTRF.  U is stored as an upper triangular band
!>          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!>          the multipliers used during the factorization are stored in
!>          rows KL+KU+2 to 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices; for 1 <= i <= N, row i of the matrix was
!>          interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is real(rk8) array, dimension (LDB,NRHS)
!>          On entry, the right hand side matrix B.
!>          On exit, the solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup doubleGBcomputational
!
!  =====================================================================
  subroutine dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    implicit none
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    character :: trans
    integer :: info, kl, ku, ldab, ldb, n, nrhs
!     ..
!     .. Array Arguments ..
    integer :: ipiv(*)
    real (rk8) :: ab(ldab, *), b(ldb, *)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    real (rk8) :: one
    parameter (one=1.0_rk8)
!     ..
!     .. Local Scalars ..
    logical :: lnoti, notran
    integer :: i, j, kd, l, lm
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: max, min
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
    info = 0
    notran = lsame(trans, 'N')
    if (.not. notran .and. .not. lsame(trans,'T') .and. &
      .not. lsame(trans,'C')) then
      info = -1
    else if (n<0) then
      info = -2
    else if (kl<0) then
      info = -3
    else if (ku<0) then
      info = -4
    else if (nrhs<0) then
      info = -5
    else if (ldab<(2*kl+ku+1)) then
      info = -7
    else if (ldb<max(1,n)) then
      info = -10
    end if
    if (info/=0) then
      call xerbla('DGBTRS', -info)
      return
    end if
!
!     Quick return if possible
!
    if (n==0 .or. nrhs==0) return
!
    kd = ku + kl + 1
    lnoti = kl > 0
!
    if (notran) then
!
!        Solve  A*X = B.
!
!        Solve L*X = B, overwriting B with X.
!
!        L is represented as a product of permutations and unit lower
!        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!        where each transformation L(i) is a rank-one modification of
!        the identity matrix.
!
      if (lnoti) then
        do j = 1, n - 1
          lm = min(kl, n-j)
          l = ipiv(j)
          if (l/=j) call dswap(nrhs, b(l,1), ldb, b(j,1), ldb)
          call dger(lm, nrhs, -one, ab(kd+1,j), 1, b(j,1), ldb, b(j+1,1), ldb)
        end do
      end if
!
      do i = 1, nrhs
!
!           Solve U*X = B, overwriting B with X.
!
        call dtbsv('Upper', 'No transpose', 'Non-unit', n, kl+ku, ab, ldab, &
          b(1,i), 1)
      end do
!
    else
!
!        Solve A**T*X = B.
!
      do i = 1, nrhs
!
!           Solve U**T*X = B, overwriting B with X.
!
        call dtbsv('Upper', 'Transpose', 'Non-unit', n, kl+ku, ab, ldab, &
          b(1,i), 1)
      end do
!
!        Solve L**T*X = B, overwriting B with X.
!
      if (lnoti) then
        do j = n - 1, 1, -1
          lm = min(kl, n-j)
          call dgemv('Transpose', lm, nrhs, -one, b(j+1,1), ldb, ab(kd+1,j), &
            1, one, b(j,1), ldb)
          l = ipiv(j)
          if (l/=j) call dswap(nrhs, b(l,1), ldb, b(j,1), ldb)
        end do
      end if
    end if
    return
!
!     End of DGBTRS
!
  end subroutine dgbtrs
!> \brief \b DLASWP performs a series of row interchanges on a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASWP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaswp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaswp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaswp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, K1, K2, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       real(rk8)   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASWP performs a series of row interchanges on the matrix A.
!> One row interchange is initiated for each of rows K1 through K2 of A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is real(rk8) array, dimension (LDA,N)
!>          On entry, the matrix of column dimension N to which the row
!>          interchanges will be applied.
!>          On exit, the permuted matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[in] K1
!> \verbatim
!>          K1 is INTEGER
!>          The first element of IPIV for which a row interchange will
!>          be done.
!> \endverbatim
!>
!> \param[in] K2
!> \verbatim
!>          K2 is INTEGER
!>          The last element of IPIV for which a row interchange will
!>          be done.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (K2*abs(INCX))
!>          The vector of pivot indices.  Only the elements in positions
!>          K1 through K2 of IPIV are accessed.
!>          IPIV(K) = L implies rows K and L are to be interchanged.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of IPIV.  If IPIV
!>          is negative, the pivots are applied in reverse order.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date September 2012
!
!> \ingroup doubleOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified by
!>   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!> \endverbatim
!>
!  =====================================================================
  subroutine dlaswp(n, a, lda, k1, k2, ipiv, incx)
    implicit none
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
    integer :: incx, k1, k2, lda, n
!     ..
!     .. Array Arguments ..
    integer :: ipiv(*)
    real (rk8) :: a(lda, *)
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
    integer :: i, i1, i2, inc, ip, ix, ix0, j, k, n32
    real (rk8) :: temp
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
    if (incx>0) then
      ix0 = k1
      i1 = k1
      i2 = k2
      inc = 1
    else if (incx<0) then
      ix0 = 1 + (1-k2)*incx
      i1 = k2
      i2 = k1
      inc = -1
    else
      return
    end if
!
    n32 = (n/32)*32
    if (n32/=0) then
      do j = 1, n32, 32
        ix = ix0
        do i = i1, i2, inc
          ip = ipiv(ix)
          if (ip/=i) then
            do k = j, j + 31
              temp = a(i, k)
              a(i, k) = a(ip, k)
              a(ip, k) = temp
            end do
          end if
          ix = ix + incx
        end do
      end do
    end if
    if (n32/=n) then
      n32 = n32 + 1
      ix = ix0
      do i = i1, i2, inc
        ip = ipiv(ix)
        if (ip/=i) then
          do k = n32, n
            temp = a(i, k)
            a(i, k) = a(ip, k)
            a(ip, k) = temp
          end do
        end if
        ix = ix + incx
      end do
    end if
!
    return
!
!     End of DLASWP
!
  end subroutine dlaswp
!> \brief \b IEEECK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download IEEECK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!       .. Scalar Arguments ..
!       INTEGER            ISPEC
!       REAL               ONE, ZERO
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> IEEECK is called from the ILAENV to verify that Infinity and
!> possibly NaN arithmetic is safe (i.e. will not trap).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies whether to test just for inifinity arithmetic
!>          or whether to test for infinity and NaN arithmetic.
!>          = 0: Verify infinity arithmetic only.
!>          = 1: Verify infinity and NaN arithmetic.
!> \endverbatim
!>
!> \param[in] ZERO
!> \verbatim
!>          ZERO is REAL
!>          Must contain the value 0.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!> \endverbatim
!>
!> \param[in] ONE
!> \verbatim
!>          ONE is REAL
!>          Must contain the value 1.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!>
!>  RETURN VALUE:  INTEGER
!>          = 0:  Arithmetic failed to produce the correct answers
!>          = 1:  Arithmetic produced the correct answers
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
  integer function ieeeck(ispec, zero, one)
    implicit none
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    integer :: ispec
    real :: one, zero
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
    real :: nan1, nan2, nan3, nan4, nan5, nan6, neginf, negzro, newzro, posinf
!     ..
!     .. Executable Statements ..
    ieeeck = 1
!
    posinf = one/zero
    if (posinf<=one) then
      ieeeck = 0
      return
    end if
!
    neginf = -one/zero
    if (neginf>=zero) then
      ieeeck = 0
      return
    end if
!
    negzro = one/(neginf+one)
    if (negzro/=zero) then
      ieeeck = 0
      return
    end if
!
    neginf = one/negzro
    if (neginf>=zero) then
      ieeeck = 0
      return
    end if
!
    newzro = negzro + zero
    if (newzro/=zero) then
      ieeeck = 0
      return
    end if
!
    posinf = one/newzro
    if (posinf<=one) then
      ieeeck = 0
      return
    end if
!
    neginf = neginf*posinf
    if (neginf>=zero) then
      ieeeck = 0
      return
    end if
!
    posinf = posinf*posinf
    if (posinf<=one) then
      ieeeck = 0
      return
    end if
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
    if (ispec==0) return
!
    nan1 = posinf + neginf
!
    nan2 = posinf/neginf
!
    nan3 = posinf/posinf
!
    nan4 = posinf*zero
!
    nan5 = neginf*negzro
!
    nan6 = nan5*zero
!
    if (nan1==nan1) then
      ieeeck = 0
      return
    end if
!
    if (nan2==nan2) then
      ieeeck = 0
      return
    end if
!
    if (nan3==nan3) then
      ieeeck = 0
      return
    end if
!
    if (nan4==nan4) then
      ieeeck = 0
      return
    end if
!
    if (nan5==nan5) then
      ieeeck = 0
      return
    end if
!
    if (nan6==nan6) then
      ieeeck = 0
      return
    end if
!
    return
  end function ieeeck
!> \brief \b ILAENV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILAENV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!       .. Scalar Arguments ..
!       CHARACTER*( * )    NAME, OPTS
!       INTEGER            ISPEC, N1, N2, N3, N4
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAENV is called from the LAPACK routines to choose problem-dependent
!> parameters for the local environment.  See ISPEC for a description of
!> the parameters.
!>
!> ILAENV returns an INTEGER
!> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!>
!> This version provides a set of parameters which should give good,
!> but not optimal, performance on many of the currently available
!> computers.  Users are encouraged to modify this subroutine to set
!> the tuning parameters for their particular machine using the option
!> and problem size information in the arguments.
!>
!> This routine will not function correctly if it is converted to all
!> lower case.  Converting it to all upper case is allowed.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies the parameter to be returned as the value of
!>          ILAENV.
!>          = 1: the optimal blocksize; if this value is 1, an unblocked
!>               algorithm will give the best performance.
!>          = 2: the minimum block size for which the block routine
!>               should be used; if the usable block size is less than
!>               this value, an unblocked routine should be used.
!>          = 3: the crossover point (in a block routine, for N less
!>               than this value, an unblocked routine should be used)
!>          = 4: the number of shifts, used in the nonsymmetric
!>               eigenvalue routines (DEPRECATED)
!>          = 5: the minimum column dimension for blocking to be used;
!>               rectangular blocks must have dimension at least k by m,
!>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!>          = 6: the crossover point for the SVD (when reducing an m by n
!>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!>               this value, a QR factorization is used first to reduce
!>               the matrix to a triangular form.)
!>          = 7: the number of processors
!>          = 8: the crossover point for the multishift QR method
!>               for nonsymmetric eigenvalue problems (DEPRECATED)
!>          = 9: maximum size of the subproblems at the bottom of the
!>               computation tree in the divide-and-conquer algorithm
!>               (used by xGELSD and xGESDD)
!>          =10: ieee NaN arithmetic can be trusted not to trap
!>          =11: infinity arithmetic can be trusted not to trap
!>          12 <= ISPEC <= 16:
!>               xHSEQR or one of its subroutines,
!>               see IPARMQ for detailed explanation
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is CHARACTER*(*)
!>          The name of the calling subroutine, in either upper case or
!>          lower case.
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER*(*)
!>          The character options to the subroutine NAME, concatenated
!>          into a single character string.  For example, UPLO = 'U',
!>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!>          be specified as OPTS = 'UTN'.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!> \endverbatim
!>
!> \param[in] N3
!> \verbatim
!>          N3 is INTEGER
!> \endverbatim
!>
!> \param[in] N4
!> \verbatim
!>          N4 is INTEGER
!>          Problem dimensions for the subroutine NAME; these may not all
!>          be required.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The following conventions have been used when calling ILAENV from the
!>  LAPACK routines:
!>  1)  OPTS is a concatenation of all of the character options to
!>      subroutine NAME, in the same order that they appear in the
!>      argument list for NAME, even if they are not used in determining
!>      the value of the parameter specified by ISPEC.
!>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!>      that they appear in the argument list for NAME.  N1 is used
!>      first, N2 second, and so on, and unused problem dimensions are
!>      passed a value of -1.
!>  3)  The parameter value returned by ILAENV is checked for validity in
!>      the calling subroutine.  For example, ILAENV is used to retrieve
!>      the optimal blocksize for STRTRI as follows:
!>
!>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!>      IF( NB.LE.1 ) NB = MAX( 1, N )
!> \endverbatim
!>
!  =====================================================================
  integer function ilaenv(ispec, name, n1, n2, n3, n4)
    implicit none
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    character (len=*) :: name
    integer :: ispec, n1, n2, n3, n4
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
    integer :: i, ic, iz, nb, nbmin, nx
    logical :: cname, sname
    character :: c1*1, c2*2, c4*2, c3*3, subnam*6
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: char, ichar, int, min, real
!     ..
!     .. Executable Statements ..
!
    select case (ispec)
    case (1, 2, 3)
      go to 100
    case (4)
      go to 140
    case (5)
      go to 150
    case (6, 7)
      go to 160
    case (8)
      go to 170
    case (9)
      go to 180
    case (10)
      go to 190
    case (11)
      go to 200
    case (12, 13, 14, 15, 16)
      go to 210
    end select
!
!     Invalid value for ISPEC
!
    ilaenv = -1
    return
!
100 continue
!
!     Convert NAME to upper case if the first character is lower case.
!
    ilaenv = 1
    subnam = name
    ic = ichar(subnam(1:1))
    iz = ichar('Z')
    if (iz==90 .or. iz==122) then
!
!        ASCII character set
!
      if (ic>=97 .and. ic<=122) then
        subnam(1:1) = char(ic-32)
        do i = 2, 6
          ic = ichar(subnam(i:i))
          if (ic>=97 .and. ic<=122) subnam(i:i) = char(ic-32)
        end do
      end if
!
    else if (iz==233 .or. iz==169) then
!
!        EBCDIC character set
!
      if ((ic>=129 .and. ic<=137) .or. (ic>=145 .and. ic<=153) .or. &
        (ic>=162 .and. ic<=169)) then
        subnam(1:1) = char(ic+64)
        do i = 2, 6
          ic = ichar(subnam(i:i))
          if ((ic>=129 .and. ic<=137) .or. (ic>=145 .and. ic<=153) .or. &
            (ic>=162 .and. ic<=169)) subnam(i:i) = char(ic+64)
        end do
      end if
!
    else if (iz==218 .or. iz==250) then
!
!        Prime machines:  ASCII+128
!
      if (ic>=225 .and. ic<=250) then
        subnam(1:1) = char(ic-32)
        do i = 2, 6
          ic = ichar(subnam(i:i))
          if (ic>=225 .and. ic<=250) subnam(i:i) = char(ic-32)
        end do
      end if
    end if
!
    c1 = subnam(1:1)
    sname = c1 == 'S' .or. c1 == 'D'
    cname = c1 == 'C' .or. c1 == 'Z'
    if (.not. (cname .or. sname)) return
    c2 = subnam(2:3)
    c3 = subnam(4:6)
    c4 = c3(2:3)
!
    select case (ispec)
    case (1)
      go to 110
    case (2)
      go to 120
    case (3)
      go to 130
    end select
!
110 continue
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
    nb = 1
!
    if (c2=='GE') then
      if (c3=='TRF') then
        if (sname) then
          nb = 64
        else
          nb = 64
        end if
      else if (c3=='QRF' .or. c3=='RQF' .or. c3=='LQF' .or. c3=='QLF') then
        if (sname) then
          nb = 32
        else
          nb = 32
        end if
      else if (c3=='HRD') then
        if (sname) then
          nb = 32
        else
          nb = 32
        end if
      else if (c3=='BRD') then
        if (sname) then
          nb = 32
        else
          nb = 32
        end if
      else if (c3=='TRI') then
        if (sname) then
          nb = 64
        else
          nb = 64
        end if
      end if
    else if (c2=='PO') then
      if (c3=='TRF') then
        if (sname) then
          nb = 64
        else
          nb = 64
        end if
      end if
    else if (c2=='SY') then
      if (c3=='TRF') then
        if (sname) then
          nb = 64
        else
          nb = 64
        end if
      else if (sname .and. c3=='TRD') then
        nb = 32
      else if (sname .and. c3=='GST') then
        nb = 64
      end if
    else if (cname .and. c2=='HE') then
      if (c3=='TRF') then
        nb = 64
      else if (c3=='TRD') then
        nb = 32
      else if (c3=='GST') then
        nb = 64
      end if
    else if (sname .and. c2=='OR') then
      if (c3(1:1)=='G') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nb = 32
        end if
      else if (c3(1:1)=='M') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nb = 32
        end if
      end if
    else if (cname .and. c2=='UN') then
      if (c3(1:1)=='G') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nb = 32
        end if
      else if (c3(1:1)=='M') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nb = 32
        end if
      end if
    else if (c2=='GB') then
      if (c3=='TRF') then
        if (sname) then
          if (n4<=64) then
            nb = 1
          else
            nb = 32
          end if
        else
          if (n4<=64) then
            nb = 1
          else
            nb = 32
          end if
        end if
      end if
    else if (c2=='PB') then
      if (c3=='TRF') then
        if (sname) then
          if (n2<=64) then
            nb = 1
          else
            nb = 32
          end if
        else
          if (n2<=64) then
            nb = 1
          else
            nb = 32
          end if
        end if
      end if
    else if (c2=='TR') then
      if (c3=='TRI') then
        if (sname) then
          nb = 64
        else
          nb = 64
        end if
      end if
    else if (c2=='LA') then
      if (c3=='UUM') then
        if (sname) then
          nb = 64
        else
          nb = 64
        end if
      end if
    else if (sname .and. c2=='ST') then
      if (c3=='EBZ') then
        nb = 1
      end if
    end if
    ilaenv = nb
    return
!
120 continue
!
!     ISPEC = 2:  minimum block size
!
    nbmin = 2
    if (c2=='GE') then
      if (c3=='QRF' .or. c3=='RQF' .or. c3=='LQF' .or. c3=='QLF') then
        if (sname) then
          nbmin = 2
        else
          nbmin = 2
        end if
      else if (c3=='HRD') then
        if (sname) then
          nbmin = 2
        else
          nbmin = 2
        end if
      else if (c3=='BRD') then
        if (sname) then
          nbmin = 2
        else
          nbmin = 2
        end if
      else if (c3=='TRI') then
        if (sname) then
          nbmin = 2
        else
          nbmin = 2
        end if
      end if
    else if (c2=='SY') then
      if (c3=='TRF') then
        if (sname) then
          nbmin = 8
        else
          nbmin = 8
        end if
      else if (sname .and. c3=='TRD') then
        nbmin = 2
      end if
    else if (cname .and. c2=='HE') then
      if (c3=='TRD') then
        nbmin = 2
      end if
    else if (sname .and. c2=='OR') then
      if (c3(1:1)=='G') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nbmin = 2
        end if
      else if (c3(1:1)=='M') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nbmin = 2
        end if
      end if
    else if (cname .and. c2=='UN') then
      if (c3(1:1)=='G') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nbmin = 2
        end if
      else if (c3(1:1)=='M') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nbmin = 2
        end if
      end if
    end if
    ilaenv = nbmin
    return
!
130 continue
!
!     ISPEC = 3:  crossover point
!
    nx = 0
    if (c2=='GE') then
      if (c3=='QRF' .or. c3=='RQF' .or. c3=='LQF' .or. c3=='QLF') then
        if (sname) then
          nx = 128
        else
          nx = 128
        end if
      else if (c3=='HRD') then
        if (sname) then
          nx = 128
        else
          nx = 128
        end if
      else if (c3=='BRD') then
        if (sname) then
          nx = 128
        else
          nx = 128
        end if
      end if
    else if (c2=='SY') then
      if (sname .and. c3=='TRD') then
        nx = 32
      end if
    else if (cname .and. c2=='HE') then
      if (c3=='TRD') then
        nx = 32
      end if
    else if (sname .and. c2=='OR') then
      if (c3(1:1)=='G') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nx = 128
        end if
      end if
    else if (cname .and. c2=='UN') then
      if (c3(1:1)=='G') then
        if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. &
          c4=='HR' .or. c4=='TR' .or. c4=='BR') then
          nx = 128
        end if
      end if
    end if
    ilaenv = nx
    return
!
140 continue
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
    ilaenv = 6
    return
!
150 continue
!
!     ISPEC = 5:  minimum column dimension (not used)
!
    ilaenv = 2
    return
!
!  100 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
    ilaenv = int(real(min(n1,n2))*1.6E0)
    return
!
160 continue
!
!     ISPEC = 7:  number of processors (not used)
!
    ilaenv = 1
    return
!
170 continue
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
    ilaenv = 50
    return
!
180 continue
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
    ilaenv = 25
    return
!
190 continue
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
    ilaenv = 1
    if (ilaenv==1) then
      ilaenv = ieeeck(1, 0.0, 1.0)
    end if
    return
!
200 continue
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
    ilaenv = 1
    if (ilaenv==1) then
      ilaenv = ieeeck(0, 0.0, 1.0)
    end if
    return
!
210 continue
!
!     12 <= ISPEC <= 16: xHSEQR or one of its subroutines.
!
    ilaenv = iparmq(ispec, n2, n3)
    return
!
!     End of ILAENV
!
  end function ilaenv
!> \brief \b IPARMQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download IPARMQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, ISPEC, LWORK, N
!       CHARACTER          NAME*( * ), OPTS*( * )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      This program sets problem and machine dependent parameters
!>      useful for xHSEQR and its subroutines. It is called whenever
!>      ILAENV is called with 12 <= ISPEC <= 16
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is integer scalar
!>              ISPEC specifies which tunable parameter IPARMQ should
!>              return.
!>
!>              ISPEC=12: (INMIN)  Matrices of order nmin or less
!>                        are sent directly to xLAHQR, the implicit
!>                        double shift QR algorithm.  NMIN must be
!>                        at least 11.
!>
!>              ISPEC=13: (INWIN)  Size of the deflation window.
!>                        This is best set greater than or equal to
!>                        the number of simultaneous shifts NS.
!>                        Larger matrices benefit from larger deflation
!>                        windows.
!>
!>              ISPEC=14: (INIBL) Determines when to stop nibbling and
!>                        invest in an (expensive) multi-shift QR sweep.
!>                        If the aggressive early deflation subroutine
!>                        finds LD converged eigenvalues from an order
!>                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
!>                        then the next QR sweep is skipped and early
!>                        deflation is applied immediately to the
!>                        remaining active diagonal block.  Setting
!>                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!>                        multi-shift QR sweep whenever early deflation
!>                        finds a converged eigenvalue.  Setting
!>                        IPARMQ(ISPEC=14) greater than or equal to 100
!>                        prevents TTQRE from skipping a multi-shift
!>                        QR sweep.
!>
!>              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!>                        a multi-shift QR iteration.
!>
!>              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!>                        following meanings.
!>                        0:  During the multi-shift QR sweep,
!>                            xLAQR5 does not accumulate reflections and
!>                            does not use matrix-matrix multiply to
!>                            update the far-from-diagonal matrix
!>                            entries.
!>                        1:  During the multi-shift QR sweep,
!>                            xLAQR5 and/or xLAQRaccumulates reflections and uses
!>                            matrix-matrix multiply to update the
!>                            far-from-diagonal matrix entries.
!>                        2:  During the multi-shift QR sweep.
!>                            xLAQR5 accumulates reflections and takes
!>                            advantage of 2-by-2 block structure during
!>                            matrix-matrix multiplies.
!>                        (If xTRMM is slower than xGEMM, then
!>                        IPARMQ(ISPEC=16)=1 may be more efficient than
!>                        IPARMQ(ISPEC=16)=2 despite the greater level of
!>                        arithmetic work implied by the latter choice.)
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is character string
!>               Name of the calling subroutine
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is character string
!>               This is a concatenation of the string arguments to
!>               TTQRE.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer scalar
!>               N is the order of the Hessenberg matrix H.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>               It is assumed that H is already upper triangular
!>               in rows and columns 1:ILO-1 and IHI+1:N.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer scalar
!>               The amount of workspace available.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>       Little is known about how best to choose these parameters.
!>       It is possible to use different values of the parameters
!>       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!>
!>       It is probably best to choose different parameters for
!>       different matrices and different parameters at different
!>       times during the iteration, but this has not been
!>       implemented --- yet.
!>
!>
!>       The best choices of most of the parameters depend
!>       in an ill-understood way on the relative execution
!>       rate of xLAQR3 and xLAQR5 and on the nature of each
!>       particular eigenvalue problem.  Experiment may be the
!>       only practical way to determine which choices are most
!>       effective.
!>
!>       Following is a list of default values supplied by IPARMQ.
!>       These defaults may be adjusted in order to attain better
!>       performance in any particular computational environment.
!>
!>       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!>                        Default: 75. (Must be at least 11.)
!>
!>       IPARMQ(ISPEC=13) Recommended deflation window size.
!>                        This depends on ILO, IHI and NS, the
!>                        number of simultaneous shifts returned
!>                        by IPARMQ(ISPEC=15).  The default for
!>                        (IHI-ILO+1).LE.500 is NS.  The default
!>                        for (IHI-ILO+1).GT.500 is 3*NS/2.
!>
!>       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!>
!>       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!>                        a multi-shift QR iteration.
!>
!>                        If IHI-ILO+1 is ...
!>
!>                        greater than      ...but less    ... the
!>                        or equal to ...      than        default is
!>
!>                                0               30       NS =   2+
!>                               30               60       NS =   4+
!>                               60              150       NS =  10
!>                              150              590       NS =  **
!>                              590             3000       NS =  64
!>                             3000             6000       NS = 128
!>                             6000             infinity   NS = 256
!>
!>                    (+)  By default matrices of this order are
!>                         passed to the implicit double shift routine
!>                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!>                         values of NS are used only in case of a rare
!>                         xLAHQR failure.
!>
!>                    (**) The asterisks (**) indicate an ad-hoc
!>                         function increasing from 10 to 64.
!>
!>       IPARMQ(ISPEC=16) Select structured matrix multiply.
!>                        (See ISPEC=16 above for details.)
!>                        Default: 3.
!> \endverbatim
!>
!  =====================================================================
  integer function iparmq(ispec, ilo, ihi)
    implicit none
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    integer :: ihi, ilo, ispec
!
!  ================================================================
!     .. Parameters ..
    integer :: inmin, inwin, inibl, ishfts, iacc22
    parameter (inmin=12, inwin=13, inibl=14, ishfts=15, iacc22=16)
    integer :: nmin, k22min, kacmin, nibble, knwswp
    parameter (nmin=75, k22min=14, kacmin=14, nibble=14, knwswp=500)
    real :: two
    parameter (two=2.0)
!     ..
!     .. Local Scalars ..
    integer :: nh, ns
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: log, max, mod, nint, real
!     ..
!     .. Executable Statements ..
    if ((ispec==ishfts) .or. (ispec==inwin) .or. (ispec==iacc22)) then
!
!        ==== Set the number simultaneous shifts ====
!
      nh = ihi - ilo + 1
      ns = 2
      if (nh>=30) ns = 4
      if (nh>=60) ns = 10
      if (nh>=150) ns = max(10, nh/nint(log(real(nh))/log(two)))
      if (nh>=590) ns = 64
      if (nh>=3000) ns = 128
      if (nh>=6000) ns = 256
      ns = max(2, ns-mod(ns,2))
    end if
!
    if (ispec==inmin) then
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
      iparmq = nmin
!
    else if (ispec==inibl) then
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
      iparmq = nibble
!
    else if (ispec==ishfts) then
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
      iparmq = ns
!
    else if (ispec==inwin) then
!
!        ==== NW: deflation window size.  ====
!
      if (nh<=knwswp) then
        iparmq = ns
      else
        iparmq = 3*ns/2
      end if
!
    else if (ispec==iacc22) then
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
      iparmq = 0
      if (ns>=kacmin) iparmq = 1
      if (ns>=k22min) iparmq = 2
!
    else
!        ===== invalid value of ispec =====
      iparmq = -1
!
    end if
!
!     ==== End of IPARMQ ====
!
  end function iparmq
!> \brief \b LSAME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      LOGICAL FUNCTION LSAME( CA, CB )
!
!     .. Scalar Arguments ..
!      CHARACTER          CA, CB
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAME returns .TRUE. if CA is the same letter as CB regardless of
!> case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CA
!> \verbatim
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CA and CB specify the single characters to be compared.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
  logical function lsame(ca, cb)
    implicit none
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    character :: ca, cb
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
    intrinsic :: ichar
!     ..
!     .. Local Scalars ..
    integer :: inta, intb, zcode
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
    lsame = ca == cb
    if (lsame) return
!
!     Now test for equivalence if both characters are alphabetic.
!
    zcode = ichar('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
    inta = ichar(ca)
    intb = ichar(cb)
!
    if (zcode==90 .or. zcode==122) then
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
      if (inta>=97 .and. inta<=122) inta = inta - 32
      if (intb>=97 .and. intb<=122) intb = intb - 32
!
    else if (zcode==233 .or. zcode==169) then
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
      if (inta>=129 .and. inta<=137 .or. inta>=145 .and. inta<=153 .or. &
        inta>=162 .and. inta<=169) inta = inta + 64
      if (intb>=129 .and. intb<=137 .or. intb>=145 .and. intb<=153 .or. &
        intb>=162 .and. intb<=169) intb = intb + 64
!
    else if (zcode==218 .or. zcode==250) then
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
      if (inta>=225 .and. inta<=250) inta = inta - 32
      if (intb>=225 .and. intb<=250) intb = intb - 32
    end if
    lsame = inta == intb
!
!     RETURN
!
!     End of LSAME
!
  end function lsame
!> \brief \b XERBLA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download XERBLA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/xerbla.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/xerbla.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE XERBLA( SRNAME, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER*(*)      SRNAME
!       INTEGER            INFO
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> XERBLA  is an error handler for the LAPACK routines.
!> It is called by an LAPACK routine if an input parameter has an
!> invalid value.  A message is printed and execution stops.
!>
!> Installers may consider modifying the STOP statement in order to
!> call system-specific exception-handling facilities.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SRNAME
!> \verbatim
!>          SRNAME is CHARACTER*(*)
!>          The name of the routine which called XERBLA.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>          The position of the invalid parameter in the parameter list
!>          of the calling routine.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
  subroutine xerbla(srname, info)
    implicit none
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
    character (len=*) :: srname
    integer :: info
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
    intrinsic :: len_trim
!     ..
!     .. Executable Statements ..
!
    write (*, fmt=100) srname(1:len_trim(srname)), info
!
    stop
!
100 format (' ** On entry to ', a, ' parameter number ', i2, ' had ', &
      'an illegal value')
!
!     End of XERBLA
!
  end subroutine xerbla

  integer function idamax(n, dx, incx)
    implicit none
!     .. Scalar Arguments ..
    integer :: incx, n
!     ..
!     .. Array Arguments ..
    real (rk8) :: dx(*)
!     ..
!
!  Purpose
!  =======
!
!     IDAMAX finds the index of element having max. absolute value.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  =====================================================================
!
!     .. Local Scalars ..
    real (rk8) :: dmax
    integer :: i, ix
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: abs
!     ..
    idamax = 0
    if (n<1 .or. incx<=0) return
    idamax = 1
    if (n==1) return
    if (incx==1) then
!
!        code for increment equal to 1
!
      dmax = abs(dx(1))
      do i = 2, n
        if (abs(dx(i))>dmax) then
          idamax = i
          dmax = abs(dx(i))
        end if
      end do
    else
!
!        code for increment not equal to 1
!
      ix = 1
      dmax = abs(dx(1))
      ix = ix + incx
      do i = 2, n
        if (abs(dx(ix))>dmax) then
          idamax = i
          dmax = abs(dx(ix))
        end if
        ix = ix + incx
      end do
    end if
    return
  end function idamax

  subroutine dswap(n, dx, incx, dy, incy)
    implicit none
!     .. Scalar Arguments ..
    integer :: incx, incy, n
!     ..
!     .. Array Arguments ..
    real (rk8) :: dx(*), dy(*)
!     ..
!
!  Purpose
!  =======
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  =====================================================================
!
!     .. Local Scalars ..
    real (rk8) :: dtemp
    integer :: i, ix, iy, m, mp1
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: mod
!     ..
    if (n<=0) return
    if (incx==1 .and. incy==1) then
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
      m = mod(n, 3)
      if (m/=0) then
        do i = 1, m
          dtemp = dx(i)
          dx(i) = dy(i)
          dy(i) = dtemp
        end do
        if (n<3) return
      end if
      mp1 = m + 1
      do i = mp1, n, 3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i+1)
        dx(i+1) = dy(i+1)
        dy(i+1) = dtemp
        dtemp = dx(i+2)
        dx(i+2) = dy(i+2)
        dy(i+2) = dtemp
      end do
    else
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if (incx<0) ix = (-n+1)*incx + 1
      if (incy<0) iy = (-n+1)*incy + 1
      do i = 1, n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do
    end if
    return
  end subroutine dswap

  subroutine dscal(n, da, dx, incx)
    implicit none
!     .. Scalar Arguments ..
    real (rk8) :: da
    integer :: incx, n
!     ..
!     .. Array Arguments ..
    real (rk8) :: dx(*)
!     ..
!
!  Purpose
!  =======
!
!     DSCAL scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  =====================================================================
!
!     .. Local Scalars ..
    integer :: i, m, mp1, nincx
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: mod
!     ..
    if (n<=0 .or. incx<=0) return
    if (incx==1) then
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
      m = mod(n, 5)
      if (m/=0) then
        do i = 1, m
          dx(i) = da*dx(i)
        end do
        if (n<5) return
      end if
      mp1 = m + 1
      do i = mp1, n, 5
        dx(i) = da*dx(i)
        dx(i+1) = da*dx(i+1)
        dx(i+2) = da*dx(i+2)
        dx(i+3) = da*dx(i+3)
        dx(i+4) = da*dx(i+4)
      end do
    else
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do i = 1, nincx, incx
        dx(i) = da*dx(i)
      end do
    end if
    return
  end subroutine dscal

  subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
    implicit none
!     .. Scalar Arguments ..
    real (rk8) :: alpha
    integer :: incx, incy, lda, m, n
!     ..
!     .. Array Arguments ..
    real (rk8) :: a(lda, *), x(*), y(*)
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y**T + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Arguments
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real(rk8).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real(rk8) array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - real(rk8) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - real(rk8) array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================
!
!     .. Parameters ..
    real (rk8) :: zero
    parameter (zero=0.0_rk8)
!     ..
!     .. Local Scalars ..
    real (rk8) :: temp
    integer :: i, info, ix, j, jy, kx
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: max
!     ..
!
!     Test the input parameters.
!
    info = 0
    if (m<0) then
      info = 1
    else if (n<0) then
      info = 2
    else if (incx==0) then
      info = 5
    else if (incy==0) then
      info = 7
    else if (lda<max(1,m)) then
      info = 9
    end if
    if (info/=0) then
      call xerbla('DGER  ', info)
      return
    end if
!
!     Quick return if possible.
!
    if ((m==0) .or. (n==0) .or. (alpha==zero)) return
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
    if (incy>0) then
      jy = 1
    else
      jy = 1 - (n-1)*incy
    end if
    if (incx==1) then
      do j = 1, n
        if (y(jy)/=zero) then
          temp = alpha*y(jy)
          do i = 1, m
            a(i, j) = a(i, j) + x(i)*temp
          end do
        end if
        jy = jy + incy
      end do
    else
      if (incx>0) then
        kx = 1
      else
        kx = 1 - (m-1)*incx
      end if
      do j = 1, n
        if (y(jy)/=zero) then
          temp = alpha*y(jy)
          ix = kx
          do i = 1, m
            a(i, j) = a(i, j) + x(ix)*temp
            ix = ix + incx
          end do
        end if
        jy = jy + incy
      end do
    end if
!
    return
!
!     End of DGER  .
!
  end subroutine dger

  subroutine dcopy(n, dx, incx, dy, incy)
    implicit none
!     .. Scalar Arguments ..
    integer :: incx, incy, n
!     ..
!     .. Array Arguments ..
    real (rk8) :: dx(*), dy(*)
!     ..
!
!  Purpose
!  =======
!
!     DCOPY copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  =====================================================================
!
!     .. Local Scalars ..
    integer :: i, ix, iy, m, mp1
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: mod
!     ..
    if (n<=0) return
    if (incx==1 .and. incy==1) then
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
      m = mod(n, 7)
      if (m/=0) then
        do i = 1, m
          dy(i) = dx(i)
        end do
        if (n<7) return
      end if
      mp1 = m + 1
      do i = mp1, n, 7
        dy(i) = dx(i)
        dy(i+1) = dx(i+1)
        dy(i+2) = dx(i+2)
        dy(i+3) = dx(i+3)
        dy(i+4) = dx(i+4)
        dy(i+5) = dx(i+5)
        dy(i+6) = dx(i+6)
      end do
    else
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if (incx<0) ix = (-n+1)*incx + 1
      if (incy<0) iy = (-n+1)*incy + 1
      do i = 1, n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
    end if
    return
  end subroutine dcopy

  subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, &
    ldc)
    implicit none
!     .. Scalar Arguments ..
    real (rk8) :: alpha, beta
    integer :: k, lda, ldb, ldc, m, n
    character :: transa, transb
!     ..
!     .. Array Arguments ..
    real (rk8) :: a(lda, *), b(ldb, *), c(ldc, *)
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X**T,
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A**T.
!
!              TRANSA = 'C' or 'c',  op( A ) = A**T.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B**T.
!
!              TRANSB = 'C' or 'c',  op( B ) = B**T.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - real(rk8).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real(rk8) array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - real(rk8) array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - real(rk8).
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - real(rk8) array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
    intrinsic :: max
!     ..
!     .. Local Scalars ..
    real (rk8) :: temp
    integer :: i, info, j, l, nrowa, nrowb
    logical :: nota, notb
!     ..
!     .. Parameters ..
    real (rk8) :: one, zero
    parameter (one=1.0_rk8, zero=0.0_rk8)
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
    nota = lsame(transa, 'N')
    notb = lsame(transb, 'N')
    if (nota) then
      nrowa = m
    else
      nrowa = k
    end if
    if (notb) then
      nrowb = k
    else
      nrowb = n
    end if
!
!     Test the input parameters.
!
    info = 0
    if ((.not. nota) .and. (.not. lsame(transa,'C')) .and. (.not. lsame(transa &
      ,'T'))) then
      info = 1
    else if ((.not. notb) .and. (.not. lsame(transb, &
        'C')) .and. (.not. lsame(transb,'T'))) then
      info = 2
    else if (m<0) then
      info = 3
    else if (n<0) then
      info = 4
    else if (k<0) then
      info = 5
    else if (lda<max(1,nrowa)) then
      info = 8
    else if (ldb<max(1,nrowb)) then
      info = 10
    else if (ldc<max(1,m)) then
      info = 13
    end if
    if (info/=0) then
      call xerbla('DGEMM ', info)
      return
    end if
!
!     Quick return if possible.
!
    if ((m==0) .or. (n==0) .or. (((alpha==zero) .or. (k==0)) .and. (beta== &
      one))) return
!
!     And if  alpha.eq.zero.
!
    if (alpha==zero) then
      if (beta==zero) then
        do j = 1, n
          do i = 1, m
            c(i, j) = zero
          end do
        end do
      else
        do j = 1, n
          do i = 1, m
            c(i, j) = beta*c(i, j)
          end do
        end do
      end if
      return
    end if
!
!     Start the operations.
!
    if (notb) then
      if (nota) then
!
!           Form  C := alpha*A*B + beta*C.
!
        do j = 1, n
          if (beta==zero) then
            do i = 1, m
              c(i, j) = zero
            end do
          else if (beta/=one) then
            do i = 1, m
              c(i, j) = beta*c(i, j)
            end do
          end if
          do l = 1, k
            if (b(l,j)/=zero) then
              temp = alpha*b(l, j)
              do i = 1, m
                c(i, j) = c(i, j) + temp*a(i, l)
              end do
            end if
          end do
        end do
      else
!
!           Form  C := alpha*A**T*B + beta*C
!
        do j = 1, n
          do i = 1, m
            temp = zero
            do l = 1, k
              temp = temp + a(l, i)*b(l, j)
            end do
            if (beta==zero) then
              c(i, j) = alpha*temp
            else
              c(i, j) = alpha*temp + beta*c(i, j)
            end if
          end do
        end do
      end if
    else
      if (nota) then
!
!           Form  C := alpha*A*B**T + beta*C
!
        do j = 1, n
          if (beta==zero) then
            do i = 1, m
              c(i, j) = zero
            end do
          else if (beta/=one) then
            do i = 1, m
              c(i, j) = beta*c(i, j)
            end do
          end if
          do l = 1, k
            if (b(j,l)/=zero) then
              temp = alpha*b(j, l)
              do i = 1, m
                c(i, j) = c(i, j) + temp*a(i, l)
              end do
            end if
          end do
        end do
      else
!
!           Form  C := alpha*A**T*B**T + beta*C
!
        do j = 1, n
          do i = 1, m
            temp = zero
            do l = 1, k
              temp = temp + a(l, i)*b(j, l)
            end do
            if (beta==zero) then
              c(i, j) = alpha*temp
            else
              c(i, j) = alpha*temp + beta*c(i, j)
            end if
          end do
        end do
      end if
    end if
!
    return
!
!     End of DGEMM .
!
  end subroutine dgemm

  subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    implicit none
!     .. Scalar Arguments ..
    real (rk8) :: alpha
    integer :: lda, ldb, m, n
    character :: diag, side, transa, uplo
!     ..
!     .. Array Arguments ..
    real (rk8) :: a(lda, *), b(ldb, *)
!     ..
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A**T.
!
!  The matrix X is overwritten on B.
!
!  Arguments
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A**T.
!
!              TRANSA = 'C' or 'c'   op( A ) = A**T.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real(rk8).
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - real(rk8) array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - real(rk8) array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
    intrinsic :: max
!     ..
!     .. Local Scalars ..
    real (rk8) :: temp
    integer :: i, info, j, k, nrowa
    logical :: lside, nounit, upper
!     ..
!     .. Parameters ..
    real (rk8) :: one, zero
    parameter (one=1.0_rk8, zero=0.0_rk8)
!     ..
!
!     Test the input parameters.
!
    lside = lsame(side, 'L')
    if (lside) then
      nrowa = m
    else
      nrowa = n
    end if
    nounit = lsame(diag, 'N')
    upper = lsame(uplo, 'U')
!
    info = 0
    if ((.not. lside) .and. (.not. lsame(side,'R'))) then
      info = 1
    else if ((.not. upper) .and. (.not. lsame(uplo,'L'))) then
      info = 2
    else if ((.not. lsame(transa,'N')) .and. (.not. lsame(transa, &
        'T')) .and. (.not. lsame(transa,'C'))) then
      info = 3
    else if ((.not. lsame(diag,'U')) .and. (.not. lsame(diag,'N'))) then
      info = 4
    else if (m<0) then
      info = 5
    else if (n<0) then
      info = 6
    else if (lda<max(1,nrowa)) then
      info = 9
    else if (ldb<max(1,m)) then
      info = 11
    end if
    if (info/=0) then
      call xerbla('DTRSM ', info)
      return
    end if
!
!     Quick return if possible.
!
    if (m==0 .or. n==0) return
!
!     And when  alpha.eq.zero.
!
    if (alpha==zero) then
      do j = 1, n
        do i = 1, m
          b(i, j) = zero
        end do
      end do
      return
    end if
!
!     Start the operations.
!
    if (lside) then
      if (lsame(transa,'N')) then
!
!           Form  B := alpha*inv( A )*B.
!
        if (upper) then
          do j = 1, n
            if (alpha/=one) then
              do i = 1, m
                b(i, j) = alpha*b(i, j)
              end do
            end if
            do k = m, 1, -1
              if (b(k,j)/=zero) then
                if (nounit) b(k, j) = b(k, j)/a(k, k)
                do i = 1, k - 1
                  b(i, j) = b(i, j) - b(k, j)*a(i, k)
                end do
              end if
            end do
          end do
        else
          do j = 1, n
            if (alpha/=one) then
              do i = 1, m
                b(i, j) = alpha*b(i, j)
              end do
            end if
            do k = 1, m
              if (b(k,j)/=zero) then
                if (nounit) b(k, j) = b(k, j)/a(k, k)
                do i = k + 1, m
                  b(i, j) = b(i, j) - b(k, j)*a(i, k)
                end do
              end if
            end do
          end do
        end if
      else
!
!           Form  B := alpha*inv( A**T )*B.
!
        if (upper) then
          do j = 1, n
            do i = 1, m
              temp = alpha*b(i, j)
              do k = 1, i - 1
                temp = temp - a(k, i)*b(k, j)
              end do
              if (nounit) temp = temp/a(i, i)
              b(i, j) = temp
            end do
          end do
        else
          do j = 1, n
            do i = m, 1, -1
              temp = alpha*b(i, j)
              do k = i + 1, m
                temp = temp - a(k, i)*b(k, j)
              end do
              if (nounit) temp = temp/a(i, i)
              b(i, j) = temp
            end do
          end do
        end if
      end if
    else
      if (lsame(transa,'N')) then
!
!           Form  B := alpha*B*inv( A ).
!
        if (upper) then
          do j = 1, n
            if (alpha/=one) then
              do i = 1, m
                b(i, j) = alpha*b(i, j)
              end do
            end if
            do k = 1, j - 1
              if (a(k,j)/=zero) then
                do i = 1, m
                  b(i, j) = b(i, j) - a(k, j)*b(i, k)
                end do
              end if
            end do
            if (nounit) then
              temp = one/a(j, j)
              do i = 1, m
                b(i, j) = temp*b(i, j)
              end do
            end if
          end do
        else
          do j = n, 1, -1
            if (alpha/=one) then
              do i = 1, m
                b(i, j) = alpha*b(i, j)
              end do
            end if
            do k = j + 1, n
              if (a(k,j)/=zero) then
                do i = 1, m
                  b(i, j) = b(i, j) - a(k, j)*b(i, k)
                end do
              end if
            end do
            if (nounit) then
              temp = one/a(j, j)
              do i = 1, m
                b(i, j) = temp*b(i, j)
              end do
            end if
          end do
        end if
      else
!
!           Form  B := alpha*B*inv( A**T ).
!
        if (upper) then
          do k = n, 1, -1
            if (nounit) then
              temp = one/a(k, k)
              do i = 1, m
                b(i, k) = temp*b(i, k)
              end do
            end if
            do j = 1, k - 1
              if (a(j,k)/=zero) then
                temp = a(j, k)
                do i = 1, m
                  b(i, j) = b(i, j) - temp*b(i, k)
                end do
              end if
            end do
            if (alpha/=one) then
              do i = 1, m
                b(i, k) = alpha*b(i, k)
              end do
            end if
          end do
        else
          do k = 1, n
            if (nounit) then
              temp = one/a(k, k)
              do i = 1, m
                b(i, k) = temp*b(i, k)
              end do
            end if
            do j = k + 1, n
              if (a(j,k)/=zero) then
                temp = a(j, k)
                do i = 1, m
                  b(i, j) = b(i, j) - temp*b(i, k)
                end do
              end if
            end do
            if (alpha/=one) then
              do i = 1, m
                b(i, k) = alpha*b(i, k)
              end do
            end if
          end do
        end if
      end if
    end if
!
    return
!
!     End of DTRSM .
!
  end subroutine dtrsm

  subroutine dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    implicit none
!     .. Scalar Arguments ..
    integer :: incx, k, lda, n
    character :: diag, trans, uplo
!     ..
!     .. Array Arguments ..
    real (rk8) :: a(lda, *), x(*)
!     ..
!
!  Purpose
!  =======
!
!  DTBSV  solves one of the systems of equations
!
!     A*x = b,   or   A**T*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )
!  diagonals.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Arguments
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A**T*x = b.
!
!              TRANS = 'C' or 'c'   A**T*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0 .le. K.
!           Unchanged on exit.
!
!  A      - real(rk8) array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = K + 1 - J
!                    DO 10, I = MAX( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = 1 - J
!                    DO 10, I = J, MIN( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - real(rk8) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================
!
!     .. Parameters ..
    real (rk8) :: zero
    parameter (zero=0.0_rk8)
!     ..
!     .. Local Scalars ..
    real (rk8) :: temp
    integer :: i, info, ix, j, jx, kplus1, kx, l
    logical :: nounit
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: max, min
!     ..
!
!     Test the input parameters.
!
    info = 0
    if (.not. lsame(uplo,'U') .and. .not. lsame(uplo,'L')) then
      info = 1
    else if (.not. lsame(trans,'N') .and. .not. lsame(trans,'T') .and. &
        .not. lsame(trans,'C')) then
      info = 2
    else if (.not. lsame(diag,'U') .and. .not. lsame(diag,'N')) then
      info = 3
    else if (n<0) then
      info = 4
    else if (k<0) then
      info = 5
    else if (lda<(k+1)) then
      info = 7
    else if (incx==0) then
      info = 9
    end if
    if (info/=0) then
      call xerbla('DTBSV ', info)
      return
    end if
!
!     Quick return if possible.
!
    if (n==0) return
!
    nounit = lsame(diag, 'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
    if (incx<=0) then
      kx = 1 - (n-1)*incx
    else if (incx/=1) then
      kx = 1
    end if
!
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
    if (lsame(trans,'N')) then
!
!        Form  x := inv( A )*x.
!
      if (lsame(uplo,'U')) then
        kplus1 = k + 1
        if (incx==1) then
          do j = n, 1, -1
            if (x(j)/=zero) then
              l = kplus1 - j
              if (nounit) x(j) = x(j)/a(kplus1, j)
              temp = x(j)
              do i = j - 1, max(1, j-k), -1
                x(i) = x(i) - temp*a(l+i, j)
              end do
            end if
          end do
        else
          kx = kx + (n-1)*incx
          jx = kx
          do j = n, 1, -1
            kx = kx - incx
            if (x(jx)/=zero) then
              ix = kx
              l = kplus1 - j
              if (nounit) x(jx) = x(jx)/a(kplus1, j)
              temp = x(jx)
              do i = j - 1, max(1, j-k), -1
                x(ix) = x(ix) - temp*a(l+i, j)
                ix = ix - incx
              end do
            end if
            jx = jx - incx
          end do
        end if
      else
        if (incx==1) then
          do j = 1, n
            if (x(j)/=zero) then
              l = 1 - j
              if (nounit) x(j) = x(j)/a(1, j)
              temp = x(j)
              do i = j + 1, min(n, j+k)
                x(i) = x(i) - temp*a(l+i, j)
              end do
            end if
          end do
        else
          jx = kx
          do j = 1, n
            kx = kx + incx
            if (x(jx)/=zero) then
              ix = kx
              l = 1 - j
              if (nounit) x(jx) = x(jx)/a(1, j)
              temp = x(jx)
              do i = j + 1, min(n, j+k)
                x(ix) = x(ix) - temp*a(l+i, j)
                ix = ix + incx
              end do
            end if
            jx = jx + incx
          end do
        end if
      end if
    else
!
!        Form  x := inv( A**T)*x.
!
      if (lsame(uplo,'U')) then
        kplus1 = k + 1
        if (incx==1) then
          do j = 1, n
            temp = x(j)
            l = kplus1 - j
            do i = max(1, j-k), j - 1
              temp = temp - a(l+i, j)*x(i)
            end do
            if (nounit) temp = temp/a(kplus1, j)
            x(j) = temp
          end do
        else
          jx = kx
          do j = 1, n
            temp = x(jx)
            ix = kx
            l = kplus1 - j
            do i = max(1, j-k), j - 1
              temp = temp - a(l+i, j)*x(ix)
              ix = ix + incx
            end do
            if (nounit) temp = temp/a(kplus1, j)
            x(jx) = temp
            jx = jx + incx
            if (j>k) kx = kx + incx
          end do
        end if
      else
        if (incx==1) then
          do j = n, 1, -1
            temp = x(j)
            l = 1 - j
            do i = min(n, j+k), j + 1, -1
              temp = temp - a(l+i, j)*x(i)
            end do
            if (nounit) temp = temp/a(1, j)
            x(j) = temp
          end do
        else
          kx = kx + (n-1)*incx
          jx = kx
          do j = n, 1, -1
            temp = x(jx)
            ix = kx
            l = 1 - j
            do i = min(n, j+k), j + 1, -1
              temp = temp - a(l+i, j)*x(ix)
              ix = ix - incx
            end do
            if (nounit) temp = temp/a(1, j)
            x(jx) = temp
            jx = jx - incx
            if ((n-j)>=k) kx = kx - incx
          end do
        end if
      end if
    end if
!
    return
!
!     End of DTBSV .
!
  end subroutine dtbsv

  subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    implicit none
!     .. Scalar Arguments ..
    real (rk8) :: alpha, beta
    integer :: incx, incy, lda, m, n
    character :: trans
!     ..
!     .. Array Arguments ..
    real (rk8) :: a(lda, *), x(*), y(*)
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real(rk8).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real(rk8) array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - real(rk8) array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real(rk8).
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real(rk8) array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 2 Blas routine.
!  The vector and matrix arguments are not referenced when N = 0, or M = 0
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================
!
!     .. Parameters ..
    real (rk8) :: one, zero
    parameter (one=1.0_rk8, zero=0.0_rk8)
!     ..
!     .. Local Scalars ..
    real (rk8) :: temp
    integer :: i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
!     ..
!     .. Intrinsic Functions ..
    intrinsic :: max
!     ..
!
!     Test the input parameters.
!
    info = 0
    if (.not. lsame(trans,'N') .and. .not. lsame(trans,'T') .and. &
      .not. lsame(trans,'C')) then
      info = 1
    else if (m<0) then
      info = 2
    else if (n<0) then
      info = 3
    else if (lda<max(1,m)) then
      info = 6
    else if (incx==0) then
      info = 8
    else if (incy==0) then
      info = 11
    end if
    if (info/=0) then
      call xerbla('DGEMV ', info)
      return
    end if
!
!     Quick return if possible.
!
    if ((m==0) .or. (n==0) .or. ((alpha==zero) .and. (beta==one))) return
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
    if (lsame(trans,'N')) then
      lenx = n
      leny = m
    else
      lenx = m
      leny = n
    end if
    if (incx>0) then
      kx = 1
    else
      kx = 1 - (lenx-1)*incx
    end if
    if (incy>0) then
      ky = 1
    else
      ky = 1 - (leny-1)*incy
    end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
    if (beta/=one) then
      if (incy==1) then
        if (beta==zero) then
          do i = 1, leny
            y(i) = zero
          end do
        else
          do i = 1, leny
            y(i) = beta*y(i)
          end do
        end if
      else
        iy = ky
        if (beta==zero) then
          do i = 1, leny
            y(iy) = zero
            iy = iy + incy
          end do
        else
          do i = 1, leny
            y(iy) = beta*y(iy)
            iy = iy + incy
          end do
        end if
      end if
    end if
    if (alpha==zero) return
    if (lsame(trans,'N')) then
!
!        Form  y := alpha*A*x + y.
!
      jx = kx
      if (incy==1) then
        do j = 1, n
          if (x(jx)/=zero) then
            temp = alpha*x(jx)
            do i = 1, m
              y(i) = y(i) + temp*a(i, j)
            end do
          end if
          jx = jx + incx
        end do
      else
        do j = 1, n
          if (x(jx)/=zero) then
            temp = alpha*x(jx)
            iy = ky
            do i = 1, m
              y(iy) = y(iy) + temp*a(i, j)
              iy = iy + incy
            end do
          end if
          jx = jx + incx
        end do
      end if
    else
!
!        Form  y := alpha*A**T*x + y.
!
      jy = ky
      if (incx==1) then
        do j = 1, n
          temp = zero
          do i = 1, m
            temp = temp + a(i, j)*x(i)
          end do
          y(jy) = y(jy) + alpha*temp
          jy = jy + incy
        end do
      else
        do j = 1, n
          temp = zero
          ix = kx
          do i = 1, m
            temp = temp + a(i, j)*x(ix)
            ix = ix + incx
          end do
          y(jy) = y(jy) + alpha*temp
          jy = jy + incy
        end do
      end if
    end if
!
    return
!
!     End of DGEMV .
!
  end subroutine dgemv

end module lapack_dgbsv

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

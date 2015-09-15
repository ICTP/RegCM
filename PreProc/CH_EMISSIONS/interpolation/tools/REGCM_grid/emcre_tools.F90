
MODULE EMCRE_TOOLS

  INTRINSIC :: ADJUSTL, ALLOCATED, ASSOCIATED, MAXVAL, MINVAL, SIZE     &
             , SUM, TRIM, UBOUND, LBOUND, REAL, SIN, SELECTED_REAL_KIND

  ! CONSTANTS
  INTEGER,  PARAMETER :: dp       = SELECTED_REAL_KIND(12,307)
  REAL(DP), PARAMETER :: r_earth  = 6371000.0_dp! radius of the earth in m
  REAL(DP), PARAMETER :: pi       = 3.14159265358979323846_dp
  REAL(DP), PARAMETER :: N_A      = 6.022045E23_dp ! Avogadro constant [1/mol]

  ! PARAMETER
  INTEGER, PARAMETER  :: MAX_NCLASS  = 100 ! max. number of emission-classes
  INTEGER, PARAMETER  :: MAX_SPEC  = 100 ! max. number of emission-classes
  INTEGER, PARAMETER  :: MAX_HEIGHTS = 100 ! max. number of emission heights
  ! ... INTERNAL
  INTEGER, PARAMETER :: str_short = 30  ! length of short strings
  INTEGER, PARAMETER :: str_long  = 200 ! length of long strings
  INTEGER, PARAMETER :: str_vlong = 10000 ! length of very long strings
  INTEGER, PARAMETER :: iou  = 21       ! I/O unit

CONTAINS

  ! ---------------------------------------------------------------------
  SUBROUTINE spline_interpolation_base(n,x,y,y2)

     ! From numerical recepies
     IMPLICIT NONE
     INTRINSIC SIZE

     ! I/O
     INTEGER,  INTENT(IN)  :: n
     REAL(DP), INTENT(IN)  :: x(:)
     REAL(DP), INTENT(IN)  :: y(:,:,:)
     REAL(DP), INTENT(OUT) :: y2(:,:,:)

     !LOCAL
     REAL(DP)  :: sig
     REAL(DP)  :: u(SIZE(y,1),SIZE(y,2),n)
     REAL(DP)  :: p(SIZE(y,1),SIZE(y,2))
     INTEGER  :: i,k

     y2(:,:,1)=0.0_dp
     u(:,:,1) =0.0_dp

     do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p(:,:)=sig*y2(:,:,i-1)+2
        y2(:,:,i)=(sig-1._dp)/p
        u(:,:,i)=(6._dp*((y(:,:,i+1)-y(:,:,i))/(x(i+1)-x(i))-(y(:,:,i)-y(:,:,i-1))  &
             / (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(:,:,i-1))/p(:,:)
     enddo

      y2(:,:,n) =0._dp

     do k=n-1,1,-1
        y2(:,:,k)=y2(:,:,k)*y2(:,:,k+1)+u(:,:,k)
     enddo

  END SUBROUTINE spline_interpolation_base

  ! ---------------------------------------------------------------------

  SUBROUTINE spline_cubic_val(n,xa,ya,y2a,x,y)

     ! From numerical recepies
     IMPLICIT NONE

     ! I/O
     INTEGER,  INTENT(IN)  :: n
     REAL(DP), INTENT(IN)  :: xa(:)
     REAL(DP), INTENT(IN)  :: ya(:,:,:)
     REAL(DP), INTENT(IN)  :: y2a(:,:,:)
     REAL(DP), INTENT(IN)  :: x
     REAL(DP), INTENT(OUT) :: y(:,:)

     !LOCAL
     INTEGER  :: khi, klo
     REAL(DP) :: a,b,h
     INTEGER :: k = 0

     klo=1
     khi=n

     k=2
     do while (xa(k).lt.x)
        k=k+1
     enddo
     klo=k-1
     khi=k

     h=xa(khi)-xa(klo)

     if (h.eq.0) then
        write(*,*) " BAD XA INPUT SPLINE"
        STOP
     endif

     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h
     y(:,:)=a*ya(:,:,klo)+b*ya(:,:,khi)+ &
       ((a**3-a)*y2a(:,:,klo)+(b**3-b)*y2a(:,:,khi))*(h**2)/6._dp


  END SUBROUTINE spline_cubic_val

  ! ---------------------------------------------------------------------

  SUBROUTINE strcrack(str, ch, el, n)

    IMPLICIT NONE

    INTRINSIC :: INDEX, LEN_TRIM

    ! I/O
    CHARACTER(LEN=*),               INTENT(IN)  :: str
    CHARACTER,                      INTENT(IN)  :: ch
    CHARACTER(LEN=*), DIMENSION(:), POINTER     :: el
    INTEGER,                        INTENT(OUT) :: n

    ! LOCAL
    INTEGER :: idx1, idx2

    ! INIT
    IF (ASSOCIATED(el)) DEALLOCATE(el)
    NULLIFY(el)
    n = 0

    ! EMPTY STRING
    IF ( (TRIM(str) == '') .OR. (TRIM(str) == ch) ) RETURN

    idx1 = 0
    idx2 = 0
    DO
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE

       n = n + 1

    END DO

    ! ALLOCATE SPACE
    ALLOCATE(el(n))

    n = 0
    idx1 = 0
    idx2 = 0
    DO
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE

       n = n + 1

       el(n) = str(idx1:idx2-1)

    END DO

  END SUBROUTINE strcrack
  ! ---------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  FUNCTION is_numeric(string)

    IMPLICIT NONE

    INTRINSIC :: INDEX, LEN_TRIM

    ! I/O
    LOGICAL                      :: is_numeric
    CHARACTER(LEN=*), INTENT(IN) :: string

    ! LOCAL
    INTEGER                     :: n, i

    is_numeric = .true.

    n  = LEN_TRIM(string)

    IF (INDEX(string,'--') /= 0) THEN
       is_numeric = .false.
       RETURN
    END IF

    DO i=1, n
       SELECT CASE(string(i:i))
       CASE ('0','1','2','3','4','5','6','7','8','9')
       CASE ('E','e','.','-','+',',',' ')
       CASE DEFAULT
          is_numeric = .false.
          EXIT
       END SELECT
    END DO

  END FUNCTION is_numeric
  ! ------------------------------------------------------------------------


END MODULE EMCRE_TOOLS
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

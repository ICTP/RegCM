MODULE class_componentIndex
  IMPLICIT NONE

  TYPE ComponentIndex

     INTEGER :: ncomp
     INTEGER, ALLOCATABLE :: ind(:)
     CHARACTER(len=3), ALLOCATABLE :: comp(:)

  END TYPE ComponentIndex

  CONTAINS
  
    SUBROUTINE ComponentIndexConstructor(SELF, ncomp, nlist, listcomp)
      IMPLICIT NONE

      TYPE(ComponentIndex), INTENT(inout) :: SELF
      INTEGER, INTENT(in) :: ncomp, nlist
      CHARACTER(len=3), INTENT(in) :: listcomp(nlist)

      INTEGER :: i, jj

      SELF%ncomp = ncomp
      ALLOCATE(SELF%ind(ncomp),SELF%comp(ncomp))
      DO i = 1,ncomp
         SELF%ind(i) = i
      END DO

      jj = 1
      DO i = 1,nlist
         IF (listcomp(i) == '') CYCLE
         SELF%comp(jj) = listcomp(i)
         jj = jj+1
      END DO

    END SUBROUTINE ComponentIndexConstructor

    ! ---------------------------------

    INTEGER FUNCTION GetIndex(SELF,incomp)
      IMPLICIT NONE

      TYPE(ComponentIndex), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: incomp

      INTEGER :: i

      IF ( ANY(SELF%comp == incomp) ) THEN

         i = 1
         DO WHILE ( (SELF%comp(i) /= incomp) )
            i = i + 1
         END DO
         GetIndex = i
      ELSE IF ( incomp == 'H2O' ) THEN
         GetIndex = SELF%ncomp + 1
      ELSE
         write(6,*) 'getIndex: FAILED, no such component - ', incomp
         STOP
      END IF

      RETURN

    END FUNCTION GetIndex

    ! -------------------------------------

    INTEGER FUNCTION GetNcomp(SELF)

      TYPE(ComponentIndex), INTENT(in) :: SELF

      GetNcomp = SELF%ncomp

      RETURN

    END FUNCTION GetNcomp
      
    ! -------------------------------------

    LOGICAL FUNCTION IsUsed(SELF,incomp)
      IMPLICIT NONE
      
      TYPE(ComponentIndex), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: incomp

      IF ( ANY(SELF%comp == incomp) ) THEN
         IsUsed = .TRUE.
      ELSE
         IsUsed = .FALSE.
      END IF

      RETURN

    END FUNCTION
      




END MODULE class_componentIndex

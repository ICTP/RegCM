module mod_clm_tridiagonal
  !
  ! Tridiagonal matrix solution
  !
  use mod_intkinds
  use mod_realkinds

  implicit none

  private

  save

  public :: Tridiagonal

  interface Tridiagonal
#if defined(STDPAR) || defined(OPENACC) || defined(_OPENACC)
    module procedure Tridiagonal_dispatcher
#else
    module procedure Tridiagonal_cpu
#endif
  end interface

  contains
  !
  ! Tridiagonal matrix solution
  !
  
#if defined(STDPAR) || defined(OPENACC) || defined(_OPENACC)

  subroutine Tridiagonal_dispatcher(lbc, ubc, lbj, ubj, jtop, numf, filter, &
                          a, b, c, r, u)
    implicit none
    integer(ik4), intent(in)    :: lbc, ubc ! lbinning and ubing column indices
    integer(ik4), intent(in)    :: lbj, ubj ! lbinning and ubing level indices
    integer(ik4), intent(in)    :: jtop(lbc:ubc) ! top level for each column
    integer(ik4), intent(in)    :: numf          ! filter dimension
    integer(ik4), intent(in)    :: filter(1:numf)  ! filter
    ! "a" left off diagonal of tridiagonal matrix
    real(rk8), intent(in)    :: a(lbc:ubc, lbj:ubj)
    ! "b" diagonal column for tridiagonal matrix
    real(rk8), intent(in)    :: b(lbc:ubc, lbj:ubj)
    ! "c" right off diagonal tridiagonal matrix
    real(rk8), intent(in)    :: c(lbc:ubc, lbj:ubj)
    ! "r" forcing term of tridiagonal matrix
    real(rk8), intent(in)    :: r(lbc:ubc, lbj:ubj)
    real(rk8), intent(inout) :: u(lbc:ubc, lbj:ubj)    ! solution

    if (numf < 128) then
       ! Dispatch on CPU     
       call Tridiagonal_cpu (lbc, ubc, lbj, ubj, jtop, numf, filter, &
                             a, b, c, r, u)
    else
       ! Dispatch on GPU     
       call Tridiagonal_gpu (lbc, ubc, lbj, ubj, jtop, numf, filter, &
                             a, b, c, r, u)
    endif


  end subroutine Tridiagonal_dispatcher

  subroutine Tridiagonal_gpu (lbc, ubc, lbj, ubj, jtop, numf, filter, &
                          a, b, c, r, u)
    use mod_clm_type
    use mod_clm_varpar, only : nlevurb
    use mod_clm_varcon, only : icol_roof, icol_sunwall, icol_shadewall
    implicit none
    integer(ik4), intent(in)    :: lbc, ubc ! lbinning and ubing column indices
    integer(ik4), intent(in)    :: lbj, ubj ! lbinning and ubing level indices
    integer(ik4), intent(in)    :: jtop(lbc:ubc) ! top level for each column
    integer(ik4), intent(in)    :: numf          ! filter dimension
    integer(ik4), intent(in)    :: filter(1:numf)  ! filter
    ! "a" left off diagonal of tridiagonal matrix
    real(rk8), intent(in)    :: a(lbc:ubc, lbj:ubj)
    ! "b" diagonal column for tridiagonal matrix
    real(rk8), intent(in)    :: b(lbc:ubc, lbj:ubj)
    ! "c" right off diagonal tridiagonal matrix
    real(rk8), intent(in)    :: c(lbc:ubc, lbj:ubj)
    ! "r" forcing term of tridiagonal matrix
    real(rk8), intent(in)    :: r(lbc:ubc, lbj:ubj)
    real(rk8), intent(inout) :: u(lbc:ubc, lbj:ubj)    ! solution

    integer(ik4), pointer, contiguous :: ctype(:) ! column type

    integer(ik4) :: j, ci, fc        !indices
    real(rk8) :: gam(lbc:ubc,lbj:ubj)  !temporary
    real(rk8) :: bet(lbc:ubc)          !temporary

    ! Assign local pointers to derived subtypes components (column-level)

    ctype => clm3%g%l%c%itype

    ! Solve the matrix
    do concurrent (fc = 1: numf)
      ci = filter(fc)
      bet(ci) = b(ci,jtop(ci))
      !$acc loop seq 
      do j = lbj, ubj
        if ( (ctype(ci) == icol_sunwall .or. &
              ctype(ci) == icol_shadewall .or. &
              ctype(ci) == icol_roof) .and. j <= nlevurb ) then
          if ( j >= jtop(ci) ) then
            if (j == jtop(ci)) then
              u(ci,j) = r(ci,j) / bet(ci)
            else
              gam(ci,j) = c(ci,j-1) / bet(ci)
              bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
              u(ci,j) = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)
            end if
          end if
        else if ( ctype(ci) /= icol_sunwall .and. &
                  ctype(ci) /= icol_shadewall .and. &
                  ctype(ci) /= icol_roof ) then
          if (j >= jtop(ci)) then
            if (j == jtop(ci)) then
              u(ci,j) = r(ci,j) / bet(ci)
            else
              gam(ci,j) = c(ci,j-1) / bet(ci)
              bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
              u(ci,j) = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)
            end if
          end if
        end if
      end do
      !$acc loop seq
      do j = ubj-1, lbj, -1
        if ( (ctype(ci) == icol_sunwall .or. &
              ctype(ci) == icol_shadewall .or. &
              ctype(ci) == icol_roof) .and. j <= nlevurb-1 ) then
          if ( j >= jtop(ci) ) then
            u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
          end if
        else if ( ctype(ci) /= icol_sunwall .and. &
                  ctype(ci) /= icol_shadewall .and. &
                  ctype(ci) /= icol_roof ) then
          if ( j >= jtop(ci) ) then
            u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
          end if
        end if
      end do
    end do

  end subroutine Tridiagonal_gpu

#endif

  subroutine Tridiagonal_cpu (lbc, ubc, lbj, ubj, jtop, numf, filter, &
                          a, b, c, r, u)
    use mod_clm_type
    use mod_clm_varpar, only : nlevurb
    use mod_clm_varcon, only : icol_roof, icol_sunwall, icol_shadewall
    implicit none
    integer(ik4), intent(in)    :: lbc, ubc ! lbinning and ubing column indices
    integer(ik4), intent(in)    :: lbj, ubj ! lbinning and ubing level indices
    integer(ik4), intent(in)    :: jtop(lbc:ubc) ! top level for each column
    integer(ik4), intent(in)    :: numf          ! filter dimension
    integer(ik4), intent(in)    :: filter(1:numf)  ! filter
    ! "a" left off diagonal of tridiagonal matrix
    real(rk8), intent(in)    :: a(lbc:ubc, lbj:ubj)
    ! "b" diagonal column for tridiagonal matrix
    real(rk8), intent(in)    :: b(lbc:ubc, lbj:ubj)
    ! "c" right off diagonal tridiagonal matrix
    real(rk8), intent(in)    :: c(lbc:ubc, lbj:ubj)
    ! "r" forcing term of tridiagonal matrix
    real(rk8), intent(in)    :: r(lbc:ubc, lbj:ubj)
    real(rk8), intent(inout) :: u(lbc:ubc, lbj:ubj)    ! solution

    integer(ik4), pointer, contiguous :: ctype(:) ! column type

    integer(ik4) :: j, ci, fc        !indices
    real(rk8) :: gam(lbc:ubc,lbj:ubj)  !temporary
    real(rk8) :: bet(lbc:ubc)          !temporary

    ! Assign local pointers to derived subtypes components (column-level)

    ctype => clm3%g%l%c%itype

    ! Solve the matrix
    

    do fc = 1, numf
      ci = filter(fc)
      bet(ci) = b(ci,jtop(ci))
    end do

    do j = lbj, ubj
      do fc = 1, numf
        ci = filter(fc)
        if ( (ctype(ci) == icol_sunwall .or. &
              ctype(ci) == icol_shadewall .or. &
              ctype(ci) == icol_roof) .and. j <= nlevurb ) then
          if ( j >= jtop(ci) ) then
            if (j == jtop(ci)) then
              u(ci,j) = r(ci,j) / bet(ci)
            else
              gam(ci,j) = c(ci,j-1) / bet(ci)
              bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
              u(ci,j) = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)
            end if
          end if
        else if ( ctype(ci) /= icol_sunwall .and. &
                  ctype(ci) /= icol_shadewall .and. &
                  ctype(ci) /= icol_roof ) then
          if (j >= jtop(ci)) then
            if (j == jtop(ci)) then
              u(ci,j) = r(ci,j) / bet(ci)
            else
              gam(ci,j) = c(ci,j-1) / bet(ci)
              bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
              u(ci,j) = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)
            end if
          end if
        end if
      end do
    end do

    do j = ubj-1, lbj, -1
      do fc = 1, numf
        ci = filter(fc)
        if ( (ctype(ci) == icol_sunwall .or. &
              ctype(ci) == icol_shadewall .or. &
              ctype(ci) == icol_roof) .and. j <= nlevurb-1 ) then
          if ( j >= jtop(ci) ) then
            u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
          end if
        else if ( ctype(ci) /= icol_sunwall .and. &
                  ctype(ci) /= icol_shadewall .and. &
                  ctype(ci) /= icol_roof ) then
          if ( j >= jtop(ci) ) then
            u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
          end if
        end if
      end do
    end do
    
  end subroutine Tridiagonal_cpu

end module mod_clm_tridiagonal
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

module mod_clm_regcm
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_regcm_types
  use mod_clm_initialize
  use mod_clm_atmlnd , only : clm_a2l , clm_l2a

  private

  public :: initclm45

  contains

  subroutine initclm45(cl,lm,lms)
    implicit none
    type (masked_comm) , intent(in) :: cl
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: i , j , n
    call initialize1(cl)
    call initialize2( )
    ! Compute simple LAND emissivity for RegCM radiation
    if ( iemiss == 1 ) then
      do n = 1 , nnsg
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              if ( lm%iveg1(n,j,i) == 8 ) then
                lms%emisv(n,j,i) = 0.76D0
              else if ( lm%iveg1(n,j,i) == 11 ) then
                lms%emisv(n,j,i) = 0.85D0
              else if ( lm%iveg1(n,j,i) == 12 ) then
                lms%emisv(n,j,i) = 0.97D0
              else
                ! We should take this from CLM !
              end if
            end if
          end do
        end do
      end do
    else
      do n = 1 , nnsg
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              lms%emisv(n,j,i) = 0.9995D0
            end if
          end do
        end do
      end do
    end if
  end subroutine initclm45

  subroutine runclm45(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms

    ! Fill clm_a2l

    ! Run CLM
    ! call clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)

    ! Get back data from clm_l2a

  end subroutine runclm45

  subroutine albedoclm45(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    ! Just get albedoes from clm_l2a
  end subroutine albedoclm45

end module mod_clm_regcm

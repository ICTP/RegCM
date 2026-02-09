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

module mod_che_cumtran
!
! Tracer convective transport
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams
  use mod_che_common

  implicit none (type, external)

  private

  public :: init_cumtran, cumtran

  logical, dimension(:,:), pointer, contiguous :: dotran
  real(rkx), dimension(:,:,:,:), pointer, contiguous :: chiten0

  interface cumtran
    module procedure cumtran1
    module procedure cumtran2
  end interface cumtran

  contains

  subroutine init_cumtran
    implicit none (type, external)
    integer(ik4) :: i, j
    call getmem(dotran,jci1,jci2,ici1,ici2,'cumtran:dotran')
    if ( ichdiag > 0 ) then
      call getmem(chiten0,jci1,jci2, &
                            ici1,ici2,1,kz,1,ntr,'che_common:chiten0')
    end if
    dotran(:,:) = .false.
    ! Emanuel anf Tiedtke do their transport internally
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cveg2d(j,i) == 14 .or. cveg2d(j,i) == 15 ) then
          if ( icup_ocn /= 4 .and. icup_ocn /= 5 ) then
            dotran(j,i) = .true.
          end if
        else
          if ( icup_lnd /= 4 .and. icup_lnd /= 5 ) then
            dotran(j,i) = .true.
          end if
        end if
      end do
    end do
  end subroutine init_cumtran

  subroutine cumtran1(mxc)
    implicit none (type, external)
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mxc
    real(rkx) :: chibar, deltas, cumfrc
    integer(ik4) :: i, j, k, kctop, n

    if ( ichdiag > 0 ) then
      do j = jci1, jci2
        do i = ici1, ici2
         chiten0(j,i,:,:) = mxc(j,i,:,:)
        end do
      end do
    end if

    do n = 1, ntr
      do i = ici1, ici2
        do j = jci1, jci2
          if ( .not. dotran(j,i) ) cycle
          if ( kcumtop(j,i) > 0 ) then
            deltas = d_zero
            chibar = d_zero
            kctop = max(kcumtop(j,i),4)
            do k = kctop, kz
              deltas = deltas + dsigma(k)
              chibar = chibar + mxc(j,i,k,n)*dsigma(k)
            end do
            do k = kctop, kz
              cumfrc = convcldfra(j,i,k)
              mxc(j,i,k,n) = mxc(j,i,k,n)*(d_one-cumfrc) + &
                           cumfrc*chibar/deltas
            end do
          end if
        end do
      end do
    end do
    ! here calculate a pseudo tendency.
    ! factor 2 is added since we are out of leap frog
    if ( ichdiag > 0 ) then
      do j = jci1, jci2
        do i = ici1, ici2
          cconvdiag(j,i,:,:) = cconvdiag(j,i,:,:) + &
            (mxc(j,i,:,:) - chiten0(j,i,:,:))/dt * d_two * cfdout
       end do
      end do
    end if
  end subroutine cumtran1

  subroutine cumtran2(amxc,bmxc)
    implicit none (type, external)
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: amxc, bmxc
    real(rkx) :: chiabar, chibbar, deltas, cumfrc
    integer(ik4) :: i, j, k, kctop, n

    if ( ichdiag > 0 ) then
      do j = jci1, jci2
        do i = ici1, ici2
         chiten0(j,i,:,:) = bmxc(j,i,:,:)
        end do
      end do
    end if

    do n = 1, ntr
      do i = ici1, ici2
        do j = jci1, jci2
          if ( .not. dotran(j,i) ) cycle
          if ( kcumtop(j,i) > 0 ) then
            deltas = d_zero
            chiabar = d_zero
            chibbar = d_zero
            kctop = max(kcumtop(j,i),4)
            do k = kctop, kz
              deltas = deltas + dsigma(k)
              chiabar = chiabar + amxc(j,i,k,n)*dsigma(k)
              chibbar = chibbar + bmxc(j,i,k,n)*dsigma(k)
            end do
            do k = kctop, kz
              cumfrc = convcldfra (j,i,k)
              amxc(j,i,k,n) = amxc(j,i,k,n)*(d_one-cumfrc) + &
                           cumfrc*chiabar/deltas
              bmxc(j,i,k,n) = bmxc(j,i,k,n)*(d_one-cumfrc) + &
                           cumfrc*chibbar/deltas
            end do
          end if
        end do
      end do
    end do
    ! here calculate a pseudo tendency.
    ! factor 2 is added since we are out of leap frog
    if ( ichdiag > 0 ) then
      do j = jci1, jci2
        do i = ici1, ici2
          cconvdiag(j,i,:,:) = cconvdiag(j,i,:,:) + &
            (bmxc(j,i,:,:) - chiten0(j,i,:,:))/dt * d_two * cfdout
       end do
      end do
    end if
  end subroutine cumtran2

end module mod_che_cumtran
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

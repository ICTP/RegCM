!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_cloud_tompkins

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  public :: tompkins_cldfrac

  contains
  !
  ! This subroutine computes the fractional cloud coverage
  ! which is used in radiation.
  !
  ! See A cloud scheme for data assimilation purposes
  ! ECMWF Technical Memorandum 410
  subroutine tompkins_cldfrac(rh,qc,p,ps,qcrit,fcc)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: rh , qc , p
    real(rkx) , pointer , dimension(:,:) , intent(in) :: ps , qcrit
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fcc
    integer(ik4) :: i , j , k
    real(rkx) :: rhrng , kappa , rhcrit , sig

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

#ifdef STDPAR
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
      local(rhrng,kappa,rhcrit,sig)
#else
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
#endif
          if ( qc(j,i,k) > qcrit(j,i) ) then
            rhrng = min(max(rh(j,i,k),0.001_rkx),0.999_rkx)
            sig = p(j,i,k)/ps(j,i)
            kappa = max(0.0_rkx,0.9_rkx*(sig-0.2_rkx)**0.2_rkx)
            ! Adjusted relative humidity threshold
            rhcrit = 0.70_rkx + 0.30 * sig * (1.0_rkx-sig) * &
                    (1.85_rkx + 0.95_rkx*(sig-0.5_rkx))
            if ( rhrng > rhcrit ) then
              fcc(j,i,k) = 1.0_rkx - sqrt((1.0_rkx-rhrng) / &
                     (1.0_rkx - rhcrit - kappa*(rhrng-rhcrit)))
            else
              fcc(j,i,k) = d_zero
            end if
          else
            fcc(j,i,k) = d_zero
          end if
#ifndef STDPAR
        end do
      end do
#endif
    end do

  end subroutine tompkins_cldfrac

end module mod_cloud_tompkins
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

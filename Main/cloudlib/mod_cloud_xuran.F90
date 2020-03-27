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

module mod_cloud_xuran

  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  public :: xuran_cldfrac

  contains
  !
  ! This subroutine computes the fractional cloud fraction
  ! using the semi-empirical formula of Xu and Randall (1996, JAS)
  !
  subroutine xuran_cldfrac(p,qc,qv,qs,rh,fcc)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: p , rh
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qc , qv , qs
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fcc
    integer(ik4) :: i , j , k
    real(rkx) :: botm , rm , qcld , rhrng

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( rh(j,i,k) >= rhmax ) then  ! full cloud cover
            fcc(j,i,k) = hicld
          else if ( rh(j,i,k) <= rhmin ) then
            fcc(j,i,k) = d_zero
          else
            qcld = qc(j,i,k)
            rhrng = max(rhmin,min(rhmax,rh(j,i,k)))
            if ( rhrng < rhmax ) then
              botm = exp(0.49_rkx*log((rhmax-rhrng)*qs(j,i,k)))
              rm = exp(0.25_rkx*log(rhrng))
              if ( 100._rkx*(qcld/botm) > 25.0_rkx ) then
                fcc(j,i,k) = rm
              else
                fcc(j,i,k) = rm*(d_one-exp(-100.0_rkx*(qcld/botm)))
              end if
              fcc(j,i,k) = max(fcc(j,i,k),hicld)
            else
              fcc(j,i,k) = hicld
            end if
          end if
        end do
      end do
    end do
    !
    ! Correction:
    !   Ivan Guettler, 14.10.2010.
    ! Based on: Vavrus, S. and Waliser D., 2008,
    ! An Improved Parameterization for Simulating Arctic Cloud Amount
    ! in the CCSM3 Climate Model, J. Climate
    !
    if ( larcticcorr ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! clouds below 750hPa, extremely cold conditions,
            !  when no cld microphy
            if ( p(j,i,k) >= 75000.0_rkx .and. qv(j,i,k) <= 0.003_rkx ) then
              fcc(j,i,k) = fcc(j,i,k) * &
                    max(0.15_rkx,min(d_one,qv(j,i,k)/0.003_rkx))
            end if
          end do
        end do
      end do
    end if

  end subroutine xuran_cldfrac

end module mod_cloud_xuran
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

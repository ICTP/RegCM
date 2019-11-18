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

module mod_che_pollen

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_che_common
  use mod_che_species
  use mod_che_indices

  implicit none

  private
  !
  ! Parameter usefull for wet and dry deposition of carbon aerosol
  ! densities in kg/m3
  !
  real(rkx) , public , parameter :: rhopollen = 1200.0_rkx
  !
  ! effctive dimaters ( and not radius!)  in micrometer
  ! ( should they be defined intercatively in the future ? )
  !
  real(rkx) , public , parameter :: reffpollen = 20._rkx
  !
  ! solubility of carbon aer for rain out param of giorgi and chameides
  !
  real(rkx) , parameter :: solpollen = 0.05_rkx
  !
  ! flowering factor, a raffiner en fonction calendrier floraison
  !
  real(rkx) , parameter :: ce = 1.e-4_rkx
  real(rkx) , parameter :: htc = d_one ! cover height
  real(rkx) , parameter :: uconv = d_zero

  public :: solpollen , pollen_emission

  contains

    subroutine pollen_emission(ustar, wind10, rh10, prec, convprec)
      implicit none
      real(rkx) , dimension(luc,jci1:jci2,ici1:ici2) , intent(in) :: ustar
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: wind10 , rh10
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: prec , convprec
      real(rkx) , dimension(jci1:jci2,ici1:ici2) :: emispol
      integer(ik4) :: i , j
      real (rkx) :: emispot , precip , fh , fw , fr
      !
      ! calculate the actual pollen flux corrected for meteo
      ! receive emission potential in grain/m2/hr
      !
      emispol(:,:) = d_zero

      do i = ici1 , ici2
        do j = jci1 , jci2

          ! in particle/m2 + derniere correction
          emispot = chemsrc(j,i,ipollen) * 24._rkx
          if ( emispot < 1.e-20_rkx ) cycle
          ! in kg/m2
          emispol(j,i) = emispot * mathpi / 6._rkx * &
                  (reffpollen * 1.e-6_rkx)**3 * rhopollen

          if (  rh10(j,i)*100.0_rkx < 50._rkx ) then
            fh = d_one
          else if (  rh10(j,i)*100.0_rkx > 80._rkx ) then
            fh = d_zero
          else
            fh = (80._rkx - rh10(j,i)*100.0_rkx )/ (80._rkx - 50._rkx)
          end if

          precip = (prec(j,i) + convprec(j,i) ) * 3600._rkx
          if (  precip < 1.0e-5_rkx ) then
            fr = d_one
          else if ( precip > 0.5_rkx ) then
            fr = d_zero
          else
            fr = (0.5_rkx - precip)/0.5_rkx
          end if

          ! Sofiev et al., 2006

          fw = 0.5_rkx + 1.0_rkx * ( 1._rkx - &
                     exp( -(wind10(j,i) + uconv) / 5._rkx ))

          emispol(j,i) = emispol(j,i) * ustar(1,j,i)/htc * ce * fh * fw * fr

        end do
      end do

      if ( ichdrdepo /= 2 ) then
        if ( idynamic == 3 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              chiten(j,i,kz,ipollen) = chiten(j,i,kz,ipollen) + &
                        emispol(j,i)/(cdzq(j,i,kz)*crhob3d(j,i,kz))
              ! diagnostic for source, cumul
              cemtrac(j,i,ipollen) = cemtrac(j,i,ipollen) + emispol(j,i)*cfdout
            end do
          end do
        else
          do i = ici1 , ici2
            do j = jci1 , jci2
              chiten(j,i,kz,ipollen) = chiten(j,i,kz,ipollen) + &
                        emispol(j,i)/(cdzq(j,i,kz)*crhob3d(j,i,kz))*cpsb(j,i)
              ! diagnostic for source, cumul
              cemtrac(j,i,ipollen) = cemtrac(j,i,ipollen) + emispol(j,i)*cfdout
            end do
          end do
        end if
      else if ( ichdrdepo == 2 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            !then emission is injected in the PBL scheme
            chifxuw(j,i,ipollen) = chifxuw(j,i,ipollen) + emispol(j,i)
            ! diagnostic for source, cumul
            cemtrac(j,i,ipollen) = cemtrac(j,i,ipollen) + emispol(j,i)*cfdout
          end do
        end do
      end if
    end subroutine pollen_emission

end module mod_che_pollen
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

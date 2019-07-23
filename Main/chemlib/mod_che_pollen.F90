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

  ! Parameter usefull for wet and dry deposition of carbon aerosol
  ! densities in kg/m3

  real(rkx) , public , parameter :: rhopollen   = 1200.0_rkx


  ! effctive dimaters ( and not radius!)  in micrometer
  ! ( should they be defined intercatively in the future ? )
  real(rkx) , public , parameter :: reffpollen   = 20._rkx


  !
  ! solubility of carbon aer for rain out param of giorgi and chameides
  !
  real(rkx) , parameter :: solpollen = 0.05_rkx

  public :: solpollen, pollen_emission

  contains

    subroutine pollen_emission(i, ustar, wind10, rh10, prec, convprec )
      implicit none
      integer, intent(in) :: i
      real(rkx) , dimension(jci1:jci2) , intent(in) :: ustar , wind10
      real(rkx) , dimension(jci1:jci2) , intent(in) :: rh10 , prec , convprec
      real(rkx) , dimension(jci1:jci2) :: precip , emispol
      integer(ik4) :: j
      real (rkx) :: emispot, fh,fw,fr,uconv,htc,ce

! calculate the actual pollen flux corrected for meteo
! receive emission potential in grain/m2/hr
!
      htc = d_one ! cover height
      uconv = d_zero
      precip = (prec + convprec ) * 3600._rkx
      emispol = d_zero
      ! flowering factor, a raffiner en fonction calendrier floraison
      ce = 1.e-4_rkx

      do j = jci1 , jci2

        ! in particle/m2 + derniere correction
        emispot = chemsrc(j,i,ipollen) * 24._rkx
        if (emispot < 1.e-20_rkx) cycle
        ! in kg/m2
        emispol(j) = emispot * mathpi / 6._rkx * &
                (reffpollen * 1.e-6_rkx)**3 * rhopollen

        if (  rh10(j)*100.0_rkx < 50._rkx ) then
          fh =d_one
        else if (  rh10(j)*100.0_rkx > 80._rkx ) then
          fh=d_zero
        else
          fh = (80._rkx - rh10(j)*100.0_rkx )/ (80._rkx - 50._rkx)
        end if

        if (  precip(j) < 1.0e-5_rkx ) then
          fr = d_one
        else if ( precip(j) > 0.5_rkx ) then
          fr = d_zero
        else
          fr = (0.5_rkx - precip(j))/0.5_rkx
        end if

        ! Sofiev et al., 2006

        fw = 0.5_rkx + 1.0_rkx * ( 1._rkx - &
                     exp( -(wind10(j) + uconv) / 5._rkx ))

        emispol(j) = emispol(j) * ustar(j)/ htc  * ce * fh * fw * fr

      end do


      if ( ichdrdepo /= 2 ) then
        do j = jci1 , jci2
          chiten(j,i,kz,ipollen) = chiten(j,i,kz,ipollen) + &
          emispol(j)*egrav/(dsigma(kz)*1.0e3_rkx)
          ! diagnostic for source, cumul
          cemtrac(j,i,ipollen) = cemtrac(j,i,ipollen) + emispol(j)*cfdout
        end do
      else if ( ichdrdepo == 2 ) then
        do j = jci1 , jci2
          !then emission is injected in the PBL scheme
          chifxuw(j,i,ipollen) = chifxuw(j,i,ipollen) + emispol(j)
          ! diagnostic for source, cumul
          cemtrac(j,i,ipollen) = cemtrac(j,i,ipollen) + emispol(j)*cfdout
        end do
      end if

    end subroutine pollen_emission
!
end module mod_che_pollen
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

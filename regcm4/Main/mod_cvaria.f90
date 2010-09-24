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

      module mod_cvaria
!
! Storage for the prognostic variables at tau+1,
!     decoupled variables, diagnostic variables and
!     working spaces needed in the model.
!
      use mod_dynparam
      use mod_main , only : atmstate , allocate_atmstate

      implicit none
!
      real(8) , allocatable , dimension(:,:,:) :: diffq , difft ,   &
                               & difuu , difuv , omega , td , xkc
      real(8) , allocatable , dimension(:,:) :: psc , pten
      real(8) , allocatable , dimension(:,:,:) :: phi
      real(8) , allocatable , dimension(:,:) :: psd
      real(8) , allocatable , dimension(:,:,:) :: qdot
      real(8) , allocatable , dimension(:,:) :: qvcs
!
      real(8) , allocatable , dimension(:,:,:,:) :: chi
      real(8) , allocatable , dimension(:,:,:,:) :: chic , chiten

      type(atmstate) , public :: atmx , atmc , atmt

      contains 

        subroutine allocate_mod_cvaria(lmpi)
          implicit none
          logical , intent(in) :: lmpi

          call allocate_atmstate(atmx,lmpi,0,1)
          call allocate_atmstate(atmc,lmpi,0,0)
          call allocate_atmstate(atmt,lmpi,0,0)

          if (lmpi) then
            allocate(diffq(iy,kz,jxp))
            allocate(difft(iy,kz,jxp))
            allocate(difuu(iy,kz,jxp))
            allocate(difuv(iy,kz,jxp))
            allocate(omega(iy,kz,jxp))
            allocate(td(iy,kz,jxp))
            allocate(xkc(iy,kz,jxp))
            allocate(psc(iy,jxp))
            allocate(pten(iy,jxp))
            allocate(phi(iy,kz,0:jxp))
            allocate(psd(iy,0:jxp+1))
            allocate(qdot(iy,kzp1,0:jxp+1))
            allocate(chi(iy,kz,0:jxp+1,ntr))
            allocate(chic(iy,kz,jxp,ntr))
            allocate(chiten(iy,kz,jxp,ntr))
          else
            allocate(diffq(iy,kz,jx))
            allocate(difft(iy,kz,jx))
            allocate(difuu(iy,kz,jx))
            allocate(difuv(iy,kz,jx))
            allocate(omega(iy,kz,jx))
            allocate(td(iy,kz,jx))
            allocate(xkc(iy,kz,jx))
            allocate(psc(iy,jx))
            allocate(pten(iy,jx))
            allocate(phi(iy,kz,jx))
            allocate(psd(iy,jx))
            allocate(chi(iy,kz,jx,ntr))
            allocate(chic(iy,kz,jx,ntr))
            allocate(chiten(iy,kz,jx,ntr))
            allocate(qdot(iy,kzp1,jx))
          end if
          allocate(qvcs(iy,kz))
!
          diffq = 0.0D0
          difft = 0.0D0
          difuu = 0.0D0
          difuv = 0.0D0
          omega = 0.0D0
          td = 0.0D0
          xkc = 0.0D0
          psc = 0.0D0
          pten = 0.0D0
          phi = 0.0D0
          psd = 0.0D0
          chi = 0.0D0
          chic = 0.0D0
          chiten = 0.0D0
          qdot = 0.0D0
          qvcs = 0.0D0
        end  subroutine allocate_mod_cvaria
!
      end module mod_cvaria

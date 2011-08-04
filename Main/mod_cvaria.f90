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
! Storage for the prognostic variables at tau+1, decoupled variables,
! diagnostic variables and working spaces needed in the model.
!
  use mod_runparams
  use mod_main , only : atmstate , allocate_atmstate
  use mod_memutil

  implicit none
!
  real(8) , pointer , dimension(:,:,:) :: diffq , difft , difuu , difuv
  real(8) , pointer , dimension(:,:) :: psc , pten , psd
  real(8) , pointer , dimension(:,:,:) :: phi , qdot , omega
!
  type(atmstate) , public :: atmx , atmc , aten , holtten

  contains 

    subroutine allocate_mod_cvaria
      implicit none

      call allocate_atmstate(atmx,.true.,0,1)
      call allocate_atmstate(atmc,.true.,0,0)
      call allocate_atmstate(aten,.true.,0,0)

      call getmem3d(diffq,1,iy,1,kz,1,jxp,'cvaria:diffq')
      call getmem3d(difft,1,iy,1,kz,1,jxp,'cvaria:difft')
      call getmem3d(difuu,1,iy,1,kz,1,jxp,'cvaria:difuu')
      call getmem3d(difuv,1,iy,1,kz,1,jxp,'cvaria:difuv')
      call getmem3d(omega,1,iy,1,kz,1,jxp,'cvaria:omega')
      call getmem2d(psc,1,iy,1,jxp,'cvaria:psc')
      call getmem2d(pten,1,iy,1,jxp,'cvaria:pten')
      call getmem3d(phi,1,iy,1,kz,0,jxp,'cvaria:phi')
      call getmem2d(psd,1,iy,0,jxp+1,'cvaria:psd')
      call getmem3d(qdot,1,iy,1,kzp1,0,jxp+1,'cvaria:qdot')
      if ( ibltyp == 99 ) then
        call allocate_atmstate(holtten,.true.,0,0)
      end if
    end  subroutine allocate_mod_cvaria
!
end module mod_cvaria

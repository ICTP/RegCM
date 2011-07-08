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

module mod_rad

  use mod_constants
  use mod_dynparam
  use mod_memutil

  implicit none

  real(8) , pointer , dimension(:,:) :: cldfra , cldlwc
  real(8) , pointer , dimension(:,:,:) :: heatrt
  real(8) , pointer , dimension(:,:,:) :: o3prof

  contains 

  subroutine allocate_mod_rad
    implicit none
    call getmem2d(cldfra,1,iym1,1,kz,'mod_rad:cldfra')
    call getmem2d(cldlwc,1,iym1,1,kz,'mod_rad:cldlwc')
    call getmem3d(heatrt,1,iym1,1,kz,1,jxp,'mod_rad:heatrt')
    call getmem3d(o3prof,1,iym1,1,kzp1,1,jxp,'mod_rad:o3prof')
  end subroutine  allocate_mod_rad

end module mod_rad

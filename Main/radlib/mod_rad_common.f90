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

module mod_rad_common

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_memutil

  public

  ! absnxt  - Nearest layer absorptivities
  ! abstot  - Non-adjacent layer absorptivites
  ! emstot  - Total emissivity

  ! Those need to be saved in output file
  real(rk8) , pointer , dimension(:,:,:) :: o3prof
  real(rk8) , pointer , dimension(:,:,:,:)  :: gasabsnxt
  real(rk8) , pointer , dimension(:,:,:,:)  :: gasabstot
  real(rk8) , pointer , dimension(:,:,:) :: gasemstot
  real(rk8) , pointer , dimension(:,:,:,:) :: taucldsp

  logical :: doabsems , dolw , dosw
  integer(ik4) :: ichso4 , ichbc , ichoc

  real(rk8) :: chfrovrradfr ! chfrq/rafrq

  contains 

  subroutine allocate_mod_rad_common
    implicit none
    call getmem3d(o3prof,jci1,jci2,ici1,ici2,1,kzp1,'rad:o3prof')
    if ( irrtm == 0 ) then
      call getmem4d(gasabsnxt,jci1,jci2,ici1,ici2,1,kz,1,4,'rad:gasabsnxt')
      call getmem4d(gasabstot,jci1,jci2,ici1,ici2,1,kzp1,1,kzp1,'rad:gasabstot')
      call getmem3d(gasemstot,jci1,jci2,ici1,ici2,1,kzp1,'rad:gasemstot')
    end if
    if ( ichem == 1 ) then
      call getmem4d(taucldsp,jci1,jci2,ici1,ici2,0,kz,1,nspi,'rad:taucldsp')
    end if
  end subroutine  allocate_mod_rad_common

end module mod_rad_common

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
  use mod_dynparam , only : nspi

  implicit none

  public

  ! absnxt  - Nearest layer absorptivities
  ! abstot  - Non-adjacent layer absorptivites
  ! emstot  - Total emissivity

  ! Those need to be saved in output file
  real(rkx) , pointer , dimension(:,:,:) :: o3prof
  real(rkx) , pointer , dimension(:,:,:,:)  :: gasabsnxt
  real(rkx) , pointer , dimension(:,:,:,:)  :: gasabstot
  real(rkx) , pointer , dimension(:,:,:) :: gasemstot
  real(rkx) , pointer , dimension(:,:,:,:) :: taucldsp

  real(rkx) , dimension(nspi) :: wavmax , wavmin
  real(rkx) , dimension(14) :: wavnm1 , wavnm2

  logical , save :: doabsems , dolw , dosw
  integer(ik4) , save :: ichso4 , ichbc , ichoc
  integer(ik4) , save :: kclimh , kth , ktf , ksf , kclimf
  real(rkx) :: chfrovrradfr ! chfrq/rafrq

  data wavmin/0.200_rkx , 0.245_rkx , 0.265_rkx , 0.275_rkx , 0.285_rkx ,&
              0.295_rkx , 0.305_rkx , 0.350_rkx , 0.640_rkx , 0.700_rkx ,&
              0.701_rkx , 0.701_rkx , 0.701_rkx , 0.701_rkx , 0.702_rkx ,&
              0.702_rkx , 2.630_rkx , 4.160_rkx , 4.160_rkx/

  data wavmax/0.245_rkx , 0.265_rkx , 0.275_rkx , 0.285_rkx , 0.295_rkx ,&
              0.305_rkx , 0.350_rkx , 0.640_rkx , 0.700_rkx , 5.000_rkx ,&
              5.000_rkx , 5.000_rkx , 5.000_rkx , 5.000_rkx , 5.000_rkx ,&
              5.000_rkx , 2.860_rkx , 4.550_rkx , 4.550_rkx/

  data wavnm1/ 2600._rkx, 3250._rkx, 4000._rkx, 4650._rkx, 5150._rkx, &
               6150._rkx, 7700._rkx, 8050._rkx,12850._rkx,16000._rkx, &
              22650._rkx,29000._rkx,38000._rkx,  820._rkx/

  data wavnm2/ 3250._rkx, 4000._rkx, 4650._rkx, 5150._rkx, 6150._rkx, &
               7700._rkx, 8050._rkx,12850._rkx,16000._rkx,22650._rkx, &
              29000._rkx,38000._rkx,50000._rkx, 2600._rkx/

end module mod_rad_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

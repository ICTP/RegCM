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

      module mod_param1

      use mod_regcm_param

      implicit none
 
      integer :: ibdyfrq , ifrabe , klake , nbatst , nslice
!
      real(8) :: dt , dt0 , dt2 , dtbat , dtlake , dtmin
      real(8) :: dx , dx2 , dx4 , dx8 , dx16 , dxsq
      real(8) :: abatm , abemh
      real(8) :: c200 , c201 , c203 
      real(8) :: fnudge , gnudge
      real(8) :: xkhmax , xkhz
      real(8) , dimension(nsplit) :: dtau

      end module mod_param1

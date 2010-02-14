!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_param1

      use mod_regcm_param

      implicit none
!
! COMMON /I8PARAM1/
!
      integer :: jyear , jyear0 , jyearr , ktau , ktaur , ntime
!
! COMMON /IPARAM1/
!
      integer :: ibdyfrq , ifrabe , klake , nbatst , nslice
!
! COMMON /PARAM1/
!
      real(8) :: abatm , abemh , c200 , c201 , c203 , dt , dt0 , dt2 ,  &
               & dtbat , dtlake , dtmin , dx , dx16 , dx2 , dx4 , dx8 , &
               & dxsq , fnudge , gnudge , xkhmax , xkhz , xtime
      real(8) , dimension(nsplit) :: dtau

      end module mod_param1

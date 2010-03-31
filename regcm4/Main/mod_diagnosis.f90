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

      module mod_diagnosis

      use mod_regcm_param

      implicit none
!
      real(8) :: tdadv , tdini , tqadv , tqeva , tqini , tqrai
!
      real(8) , dimension(ntr) :: tchiad , tchie , tchitb
      real(8) , dimension(ntr,2) :: tremcvc , tremdrd , tremlsc ,       &
                                  & trxsaq1 , trxsaq2 , trxsg , ttrace

      end module mod_diagnosis

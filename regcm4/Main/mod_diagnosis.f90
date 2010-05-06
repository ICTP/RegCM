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

      use mod_dynparam

      implicit none
!
      real(8) :: tdadv , tdini , tqadv , tqeva , tqini , tqrai
!
      real(8) , allocatable , dimension(:) :: tchiad , tchie , tchitb
      real(8) , allocatable , dimension(:,:) :: tremcvc , tremdrd ,     &
               & tremlsc , trxsaq1 , trxsaq2 , trxsg , ttrace
      contains

      subroutine allocate_mod_diagnosis
      implicit none
        allocate(tchiad(ntr))
        allocate(tchie(ntr))
        allocate(tchitb(ntr))
        allocate(tremcvc(ntr,2))
        allocate(tremdrd(ntr,2))
        allocate(tremlsc(ntr,2))
        allocate(trxsaq1(ntr,2))
        allocate(trxsaq2(ntr,2))
        allocate(trxsg(ntr,2))
        allocate(ttrace(ntr,2))
      end subroutine allocate_mod_diagnosis

      end module mod_diagnosis

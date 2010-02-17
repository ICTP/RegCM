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

      module mod_iunits
      implicit none
!
      character(17) :: oldsav
!
      integer :: iutbat , iutbc , iutchem , iutchsrc , iutdat , iutin , &
               & iutin1 , iutlak , iutopt , iutrad , iutrs , iutsav ,   &
               & iutsub , mindisp
!
      integer :: iin , iout , lcount , numpts
!
      integer :: nrcout
      integer :: nrcbat
      integer :: nrcsub
      integer :: nrcchem
      integer :: nrcrad

      end module mod_iunits

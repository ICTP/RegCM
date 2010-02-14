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

      module mod_bmparam

      implicit none

      integer , parameter :: itb = 100
      integer , parameter :: jtb = 150

      real(8) :: pl , rdp , rdq , rdth , rdthe , thl
      real(8) , dimension(itb,jtb) :: ptbl
      real(8) , dimension(jtb) :: qs0 , sqs , sthe , the0
      real(8) , dimension(jtb,itb) :: ttbl
!
      real(8) , parameter :: aliq = 613.3D0
      real(8) , parameter :: bliq = 17.502D0
      real(8) , parameter :: cliq = 4780.8D0
      real(8) , parameter :: dliq = 32.19D0
      real(8) , parameter :: aice = 613.2D0
      real(8) , parameter :: bice = 22.452D0
      real(8) , parameter :: cice1 = 6133.0D0
      real(8) , parameter :: dice = 0.61D0
      real(8) , parameter :: xls0 = 2.905D6
      real(8) , parameter :: xls1 = 259.532D0

      end module mod_bmparam

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
 
      subroutine ljchrn(chrv,nlshf)
!
!****     routine left justifies a character variable cnull filling on
!****     the right.  chrv is the character variable.  nlshf tells
!****     how many characters the variable was shifted left
!
      implicit none
!
! Dummy arguments
!
      character(*) :: chrv
      integer :: nlshf
      intent (inout) nlshf
!
! Local variables
!
      character(1) :: blank , cnull
      integer :: lnchrv
!
      data blank/' '/
!
      cnull = char(0)
      lnchrv = len(chrv)
      nlshf = 0
      do while ( chrv(1:1).eq.blank .or. chrv(1:1).eq.cnull )
        call lshfch(chrv,lnchrv,1,1)
        nlshf = nlshf + 1
        if ( nlshf.eq.lnchrv ) exit
      end do
      end subroutine ljchrn

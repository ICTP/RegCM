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

      subroutine comp(fields,bvoc)
      implicit none
!
! Dummy arguments
!
      integer :: fields
      logical :: bvoc
      intent (in) bvoc
      intent (out) fields
!
! Local variables
!
      integer :: numcompounds
 
      if ( bvoc ) then
 50     continue
        print * , ' '
        print * , ' '
        print * , '********************************************'
        print * , 'Creating biogenic emissions files'
        print * , 'ENTER NUMBER OF SPECIES'
        read (*,*) numcompounds
        fields = 14 + numcompounds
        if ( fields>50 ) then
          stop 999
        else if ( fields<14 ) then
          go to 50
        else
        end if
      else
        fields = 14
      end if
 
      print * , 'producing ' , fields , ' files'
      end subroutine comp

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
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine printsp(array,iy,jx,kz,iinc,jinc,kref,aname)
!
      implicit none
!
! Dummy arguments
!
      character(08) :: aname
      integer :: iinc , iy , jinc , jx , kref , kz
      real(8) , dimension(iy,jx,kz) :: array
      intent (in) aname , array , iinc , iy , jinc , jx , kref , kz
!
! Local variables
!
      character(18) :: cform
      integer :: i , ineg , j , jcheck , jnum
!
      ineg = -iinc
      jnum = (jx-2)/jinc + 1
!
      jcheck = 11*jnum
      if ( jcheck.gt.127 ) print * , 'too many values printed.'
!
      write (cform,99001) jnum
!
      print * , '  field = ' , aname
      write (6,cform) (i,(array(i,j,kref),j=1,jx-1,jinc),i=iy-1,1,ineg)
99001 format ('(1x,i2,2x,',i2,'g11.3)')
!
      end subroutine printsp

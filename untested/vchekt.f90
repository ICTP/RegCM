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
 
      subroutine vchekt(tbarh,sigmah,sigmaf,xkappa,pt,pd,nk,numerr)
      implicit none
!
! Dummy arguments
!
      integer :: nk , numerr
      real(8) :: pd , pt , xkappa
      real(8) , dimension(nk+1) :: sigmaf
      real(8) , dimension(nk) :: sigmah , tbarh
      intent (in) nk , pd , pt , sigmaf , sigmah , tbarh , xkappa
      intent (inout) numerr
!
! Local variables
!
      real(8) :: ds1 , ds2 , g1 , g2 , tb
      integer :: k
      logical :: lstab
!
!  check if tbar is stable.  this is not the actual stability condition
!  consistent with the model finite differences.
!
      lstab = .true.
      do k = 1 , nk - 1
        ds1 = sigmaf(k+1) - sigmaf(k)
        ds2 = sigmaf(k+2) - sigmaf(k+1)
        tb = (ds1*tbarh(k)+ds2*tbarh(k+1))/(ds1+ds2)
        g1 = xkappa*tb/(sigmaf(k+1)+pt/pd)
        g2 = (tbarh(k+1)-tbarh(k))/(sigmah(k+1)-sigmah(k))
        if ( g1-g2.lt.0. ) lstab = .false.
      end do
      if ( .not.lstab ) then
        numerr = numerr + 1
        print 99001
99001   format ('0 indication that tbarh statically unstable')
      end if
!
      end subroutine vchekt

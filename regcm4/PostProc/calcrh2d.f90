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

      subroutine calcrh2d(f2d,nx,ny,nfld,nta,nqa,npsrf,nrh,vmisdat)
 
      use mod_constants , only : svp1 , svp2 , svp3 , ep2
      implicit none
!
! Dummy arguments
!
      integer :: nfld , npsrf , nqa , nrh , nta , nx , ny
      real(4) :: vmisdat
      real(4) , dimension(nx,ny,nfld) :: f2d
      intent (in) nfld , npsrf , nqa , nrh , nta , nx , ny , vmisdat
      intent (inout) f2d
!
! Local variables
!
      integer :: i , j
      real(4) :: pres , qa , qs , satvp , ta
!
      do j = 1 , ny
        do i = 1 , nx
          pres = f2d(i,j,npsrf)
          if ( pres>0. ) then
            ta = f2d(i,j,nta)
            qa = f2d(i,j,nqa)
            if ( ta>273.15 ) then
              satvp = svp1*exp(svp2*(ta-273.15)/(ta-svp3))
                                                         ! SAT'N VAP PRES
            else
              satvp = svp1*exp(22.514-6.15E3/ta)
            end if
            qs = ep2*satvp/(pres-satvp)              ! SAT. MIXING RATIO
            f2d(i,j,nrh) = qa/qs
          else
            f2d(i,j,nrh) = vmisdat
          end if
        end do
      end do
 
      end subroutine calcrh2d

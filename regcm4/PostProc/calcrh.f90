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

      subroutine calcrh(fld2d,fld3d,nx,ny,kx,nfld2d,nfld3d,sigh,pt,nta, &
                      & nqva,npsa,nrh,ntd,nth,nx1,ny1)
 
      use mod_constants , only : svp1 , svp2 , svp3 , ep2 , rovcp
      implicit none
!
! Dummy arguments
!
      integer :: nfld2d , nfld3d , npsa , nqva , nrh , nta , ntd , nth ,&
               & nx , nx1 , ny , ny1 , kx
      real(4) :: pt
      real(4) , dimension(nx,ny,nfld2d) :: fld2d
      real(4) , dimension(nx,ny,kx,nfld3d) :: fld3d
      real(4) , dimension(kx) :: sigh
      intent (in) fld2d , nfld2d , nfld3d , npsa , nqva , nrh , nta ,   &
                & ntd , nth , nx , nx1 , ny , ny1 , kx , pt , sigh
      intent (inout) fld3d
!
! Local variables
!
      real(4) :: dpd , pres , q , qs , satvp , t , tmp , x
      integer :: i , j , k
!
      do k = 1 , kx
        do j = 1 , ny1
          do i = 1 , nx1
            pres = (fld2d(i,j,npsa)-pt*10.)*sigh(k) + pt*10.
                                                         ! PRES AT LEVEL K
            t = fld3d(i,j,k,nta)
            q = fld3d(i,j,k,nqva)
            if ( t>273.15 ) then
              satvp = svp1*exp(svp2*(t-273.15)/(t-svp3))
                                                     ! SATURATION VAP PRESS.
            else
              satvp = svp1*exp(22.514-6.15E3/t)
            end if
            qs = ep2*satvp/(pres-satvp)               ! SAT. MIXING RATIO
            fld3d(i,j,k,nrh) = q/qs
            x = 1. - fld3d(i,j,k,nrh)
            tmp = t - 273.15
            dpd = (14.55+0.144*tmp)*x + 2*((2.5+0.007*tmp)*x)           &
                & **3 + (15.9+0.117*tmp)*(x**14)
            fld3d(i,j,k,ntd) = t - dpd  ! DEW POINT TEMP
            fld3d(i,j,k,nth) = t*(1000./pres)**rovcp    ! POTENTIAL TEMP
 
          end do
        end do
      end do
      end subroutine calcrh

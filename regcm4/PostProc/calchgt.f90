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

      subroutine calchgt(fld2d,fld3d,nx,ny,kz,nfld2d,nfld3d,zs,sigf,    &
                       & sigh,pt,nta,nqva,npsa,nhgt,nx1,ny1)
 
      use mod_constants , only : ep1 , rovg
      implicit none
!
! Dummy arguments
!
      integer :: nfld2d , nfld3d , nhgt , npsa , nqva , nta , nx , nx1 ,&
               & ny , ny1 , kz
      real(4) :: pt
      real(4) , dimension(nx,ny,nfld2d) :: fld2d
      real(4) , dimension(nx,ny,kz,nfld3d) :: fld3d
      real(4) , dimension(kz+1) :: sigf
      real(4) , dimension(kz) :: sigh
      real(4) , dimension(nx,ny) :: zs
      intent (in) fld2d , nfld2d , nfld3d , nhgt , npsa , nqva , nta ,  &
                & nx , nx1 , ny , ny1 , kz , pt , sigf , sigh , zs
      intent (inout) fld3d
!
! Local variables
!
      real(4) , dimension(kz) :: dsig
      real(4) :: pf , ps , q , q1 , q2 , t , t1 , t2 , tv , tv1 , tv2
      integer :: i , j , k , k1 , k2
!
      do k = 1 , kz
        dsig(k) = sigf(k) - sigf(k+1)
      end do
 
      do i = 1 , nx1
        do j = 1 , ny1
          ps = fld2d(i,j,npsa)/10.
          t = fld3d(i,j,1,nta)
          q = fld3d(i,j,1,nqva)
          pf = pt/(ps-pt)
          tv = t*(1.0+ep1*q)
          fld3d(i,j,kz,nhgt) = zs(i,j)                                  &
                             & + tv*rovg*log((1.+pf)/(sigh(kz)+pf))
          do k1 = kz - 1 , 1 , -1
            k2 = k1 + 1
            t1 = fld3d(i,j,k1,nta)
            t2 = fld3d(i,j,k2,nta)
            q1 = fld3d(i,j,k1,nqva)
            q2 = fld3d(i,j,k2,nqva)
            tv1 = t1*(1.0+ep1*q1)
            tv2 = t2*(1.0+ep1*q2)
            tv = (tv1*dsig(k1)+tv2*dsig(k2))/(dsig(k1)+dsig(k2))
            fld3d(i,j,k1,nhgt) = fld3d(i,j,k2,nhgt)                     &
                               & + tv*rovg*log((sigh(k2)+pf)            &
                               & /(sigh(k1)+pf))
          end do
        end do
      end do
 
      end subroutine calchgt

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
 
      subroutine slice3d
 
      use mod_regcm_param
      use mod_param2
      use mod_param3
      use mod_main
      use mod_mainchem
      use mod_pbldim
      use mod_pmoist
      use mod_slice
      use mod_constants , only : rgas , rovcp , rovg , ep2 , svp1 ,     &
                               & svp2 , svp3
      implicit none
!
! Local variables
!
      real(8) :: cell , pl , pres , psrf , rh , satvp , thcon , tv
      integer :: i , idx , idxp1 , j , jdx , jdxp1 , k , kk , n
 
#ifdef MPP1
      do j = 1 , jendx
#else
      do j = 1 , jxm1
#endif
        do k = 1 , kx
          do i = 1 , ixm1
            tb3d(i,k,j) = tb(i,k,j)/psb(i,j)
            qvb3d(i,k,j) = qvb(i,k,j)/psb(i,j)
            qcb3d(i,k,j) = qcb(i,k,j)/psb(i,j)
            if ( ichem.eq.1 ) then
              do n = 1 , ntr
                chib3d(i,k,j,n) = chib(i,k,j,n)/psb(i,j)
              end do
            end if
          end do
        end do
      end do
 
#ifdef MPP1
      do j = jbegin , jendx
        jdx = j
        if ( myid.eq.0 ) jdx = max0(j,2)
        jdxp1 = j + 1
        if ( myid.eq.nproc-1 ) jdxp1 = min0(j+1,jendx)
#else
      do j = 2 , jxm1
        jdx = max0(j,2)
        jdxp1 = min0(j+1,jxm1)
#endif
        do k = 1 , kx
          do i = 1 , ixm1
            idx = max0(i,2)
            idxp1 = min0(i+1,ixm1)
            ubx3d(i,k,j) = 0.25*(ub(idx,k,jdx)+ub(idxp1,k,jdx)+ub(idx,k,&
                         & jdxp1)+ub(idxp1,k,jdxp1))/psb(i,j)
            vbx3d(i,k,j) = 0.25*(vb(idx,k,jdx)+vb(idxp1,k,jdx)+vb(idx,k,&
                         & jdxp1)+vb(idxp1,k,jdxp1))/psb(i,j)
          end do
        end do
      end do
 
#ifdef MPP1
      do j = 1 , jendl
#else
      do j = 1 , jx
#endif
        do k = 1 , kx
          do i = 1 , ix
            ubd3d(i,k,j) = ub(i,k,j)/pdotb(i,j)
            vbd3d(i,k,j) = vb(i,k,j)/pdotb(i,j)
          end do
        end do
      end do
 
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
        do k = 1 , kx
          do i = 2 , ixm1
            pl = a(k)*psb(i,j) + ptop
            thcon = ((psb(i,j)+ptop)/pl)**rovcp
            pb3d(i,k,j) = pl
            thx3d(i,k,j) = tb3d(i,k,j)*thcon
          end do
        end do
      end do
 
!-----compute the height at full (za) and half (zq) sigma levels:
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , ixm1
          zq(i,kxp1) = 0.
        end do
        do kk = 1 , kx
          k = kxp1 - kk
          do i = 2 , ixm1
            cell = ptop/psb(i,j)
            zq(i,k) = rovg*tb3d(i,k,j)                                  &
                    & *dlog((sigma(k+1)+cell)/(sigma(k)+cell))          &
                    & + zq(i,k+1)
          end do
        end do
!
        do k = 1 , kx
          do i = 2 , ixm1
            za(i,k,j) = 0.5*(zq(i,k)+zq(i,k+1))
            dzq(i,k,j) = zq(i,k) - zq(i,k+1)
          end do
        end do
 
!-----Calculate the relative humidity and air density
        do i = 2 , ixm1
          psrf = (psb(i,j)+ptop)*1000.
          tv = tb3d(i,kx,j)
          rhox2d(i,j) = psrf/(rgas*tv)
        end do
        do k = 1 , kx
          do i = 2 , ixm2
            pres = (a(k)*psb(i,j)+ptop)*1000.
            rhob3d(i,k,j) = pres/(rgas*tb3d(i,k,j)) !air density
            if ( tb3d(i,k,j).gt.273.15 ) then
              satvp = svp1*1.E3*dexp(svp2*(tb3d(i,k,j)-273.15)          &
                    & /(tb3d(i,k,j)-svp3))
            else
              satvp = .611*1.E3*dexp(22.514-6.15E3/tb3d(i,k,j))
            end if
            qsb3d(i,k,j) = ep2*satvp/(pres-satvp)
            rh = 0.
            if ( qsb3d(i,k,j).gt.0. ) rh = qvb3d(i,k,j)/qsb3d(i,k,j)
            rhb3d(i,k,j) = rh
          end do
        end do
      end do
 
      end subroutine slice3d

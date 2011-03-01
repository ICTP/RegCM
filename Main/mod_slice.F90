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

      module mod_slice
!
! Prepare and fill 3D spaces for calculations
!
      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_main
      use mod_mainchem
      use mod_pbldim
      use mod_pmoist
!
      private
!
      public :: chib3d , pb3d , qsb3d , rhb3d , rhob3d , ubx3d , vbx3d
      public :: qcb3d , qvb3d , tb3d , ubd3d , vbd3d
!
      public :: allocate_mod_slice , slice
!
      real(8) ,allocatable, dimension(:,:,:,:) :: chib3d
      real(8) ,allocatable, dimension(:,:,:) :: pb3d , qsb3d , rhb3d ,  &
                                      & rhob3d , ubx3d , vbx3d
      real(8) ,allocatable, dimension(:,:,:) :: qcb3d , qvb3d , tb3d ,  &
                                      & ubd3d , vbd3d
!
      contains 
!
      subroutine allocate_mod_slice
        implicit none   
#ifdef MPP1
        if (ichem == 1 ) then
          allocate(chib3d(iy,kz,-1:jxp+2,ntr))      
        end if
        allocate(pb3d(iy,kz,jxp))
        allocate(qsb3d(iy,kz,jxp))
        allocate(rhb3d(iy,kz,jxp))
        allocate(rhob3d(iy,kz,jxp))
        allocate(ubx3d(iy,kz,jxp))
        allocate(vbx3d(iy,kz,jxp))
        allocate(qcb3d(iy,kz,-1:jxp+2))
        allocate(qvb3d(iy,kz,-1:jxp+2))
        allocate(tb3d(iy,kz,-1:jxp+2))
        allocate(ubd3d(iy,kz,-1:jxp+2))
        allocate(vbd3d(iy,kz,-1:jxp+2))
#else
        if ( ichem == 1 ) then
          allocate(chib3d(iy,kz,jx,ntr))      
        end if
        allocate(pb3d(iy,kz,jx))
        allocate(qsb3d(iy,kz,jx))
        allocate(rhb3d(iy,kz,jx))
        allocate(rhob3d(iy,kz,jx))
        allocate(ubx3d(iy,kz,jx))
        allocate(vbx3d(iy,kz,jx))
        allocate(qcb3d(iy,kz,jx))
        allocate(qvb3d(iy,kz,jx))
        allocate(tb3d(iy,kz,jx))
        allocate(ubd3d(iy,kz,jx))
        allocate(vbd3d(iy,kz,jx))
#endif
        if ( ichem == 1 ) then
          chib3d = 0.0D0
        end if
        pb3d = 0.0D0
        qsb3d = 0.0D0
        rhb3d = 0.0D0
        rhob3d = 0.0D0
        ubx3d = 0.0D0
        vbx3d = 0.0D0
        qcb3d = 0.0D0
        qvb3d = 0.0D0
        tb3d = 0.0D0
        ubd3d = 0.0D0
        vbd3d = 0.0D0
!
      end subroutine allocate_mod_slice
!
      subroutine slice
 
      implicit none
!
      real(8) :: cell , pl , pres , psrf , rh , satvp , thcon , tv
      integer :: i , idx , idxp1 , j , jdx , jdxp1 , k , kk , n
 
#ifdef MPP1
      do j = 1 , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
#endif
        do k = 1 , kz
          do i = 1 , iym1
            tb3d(i,k,j) = atm2%t(i,k,j)/sps2%ps(i,j)
            qvb3d(i,k,j) = atm2%qv(i,k,j)/sps2%ps(i,j)
            qcb3d(i,k,j) = atm2%qc(i,k,j)/sps2%ps(i,j)
            if ( ichem == 1 ) then
              do n = 1 , ntr
                chib3d(i,k,j,n) = chib(i,k,j,n)/sps2%ps(i,j)
              end do
            end if
          end do
        end do
      end do
#ifdef BAND
#ifdef MPP1
      do j = jbegin , jendx
        jdx = j
        jdxp1 = j + 1
#else
      do j = 1 , jx
        jdx = j
        jdxp1 = j+1
        if(jdxp1 == jx+1) jdxp1 = 1
#endif
#else
#ifdef MPP1
      do j = jbegin , jendx
        jdx = j
        if ( myid == 0 ) jdx = max0(j,2)
        jdxp1 = j + 1
        if ( myid == nproc-1 ) jdxp1 = min0(j+1,jendx)
#else
      do j = 2 , jxm1
        jdx = max0(j,2)
        jdxp1 = min0(j+1,jxm1)
#endif
#endif
        do k = 1 , kz
          do i = 1 , iym1
            idx = max0(i,2)
            idxp1 = min0(i+1,iym1)
            ubx3d(i,k,j) = 0.25D0* & 
                (atm2%u(idx,k,jdx)+atm2%u(idxp1,k,jdx)+ &
                 atm2%u(idx,k,jdxp1)+atm2%u(idxp1,k,jdxp1))/sps2%ps(i,j)
            vbx3d(i,k,j) = 0.25D0* &
                (atm2%v(idx,k,jdx)+atm2%v(idxp1,k,jdx)+ &
                 atm2%v(idx,k,jdxp1)+atm2%v(idxp1,k,jdxp1))/sps2%ps(i,j)
          end do
        end do
      end do
 
#ifdef MPP1
      do j = 1 , jendl
#else
      do j = 1 , jx
#endif
        do k = 1 , kz
          do i = 1 , iy
            ubd3d(i,k,j) = atm2%u(i,k,j)/sps2%pdot(i,j)
            vbd3d(i,k,j) = atm2%v(i,k,j)/sps2%pdot(i,j)
          end do
        end do
      end do
 
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        do k = 1 , kz
          do i = 2 , iym1
            pl = a(k)*sps2%ps(i,j) + r8pt
            thcon = ((sps2%ps(i,j)+r8pt)/pl)**rovcp
            pb3d(i,k,j) = pl
            thx3d(i,k,j) = tb3d(i,k,j)*thcon
          end do
        end do
      end do
 
!-----compute the height at full (za) and half (zq) sigma levels:
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        do i = 2 , iym1
          zq(i,kzp1) = 0.0D0
        end do
        do kk = 1 , kz
          k = kzp1 - kk
          do i = 2 , iym1
            cell = r8pt/sps2%ps(i,j)
            zq(i,k) = rovg*tb3d(i,k,j)                                  &
                    & *dlog((sigma(k+1)+cell)/(sigma(k)+cell))          &
                    & + zq(i,k+1)
          end do
        end do
!
        do k = 1 , kz
          do i = 2 , iym1
            za(i,k,j) = 0.5D0*(zq(i,k)+zq(i,k+1))
            dzq(i,k,j) = zq(i,k) - zq(i,k+1)
          end do
        end do
 
!-----Calculate the relative humidity and air density

        do i = 2 , iym1
          psrf = (sps2%ps(i,j)+r8pt)*1000.0D0
          tv = tb3d(i,kz,j)
          rhox2d(i,j) = psrf/(rgas*tv)
        end do
        do k = 1 , kz
          do i = 2 , iym2
            pres = (a(k)*sps2%ps(i,j)+r8pt)*1000.0D0
            rhob3d(i,k,j) = pres/(rgas*tb3d(i,k,j)) !air density
            if ( tb3d(i,k,j) > tzero ) then
              satvp = svp1*1.D3*dexp(svp2*(tb3d(i,k,j)-tzero)           &
                    & /(tb3d(i,k,j)-svp3))
            else
              satvp = svp4*1.D3*dexp(svp5-svp6/tb3d(i,k,j))
            end if
            qsb3d(i,k,j) = ep2*satvp/(pres-satvp)
            rh = 0.0D0
            if ( qsb3d(i,k,j) > 0.0D0 ) rh = qvb3d(i,k,j)/qsb3d(i,k,j)
            rhb3d(i,k,j) = rh
          end do
        end do
      end do
 
      end subroutine slice
!
      end module mod_slice

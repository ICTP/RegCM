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
! Fill 3D spaces for calculations
!
  use mod_runparams
  use mod_atm_interface
  use mod_che_common
  use mod_pbldim
!
  private
!
  public :: mkslice
!
  contains 
!
  subroutine mkslice
 
  implicit none
!
  real(8) :: cell , pl , pres , psrf , rh , satvp , thcon , tv
  integer :: i , idx , idxp1 , j , jdx , jdxp1 , k , kk , n
 
  do j = 1 , jendx
    do k = 1 , kz
      do i = 1 , iym1
        atms%tb3d(i,k,j) = atm2%t(i,k,j)/sps2%ps(i,j)
        atms%qvb3d(i,k,j) = atm2%qv(i,k,j)/sps2%ps(i,j)
        atms%qcb3d(i,k,j) = atm2%qc(i,k,j)/sps2%ps(i,j)
        if ( ichem == 1 ) then
          do n = 1 , ntr
            atms%chib3d(i,k,j,n) = chib(i,k,j,n)/sps2%ps(i,j)
          end do
        end if
      end do
    end do
  end do
  do j = jbegin , jendx
    jdx = j
    jdxp1 = j + 1
#ifndef BAND
    if ( myid == 0 ) jdx = max0(j,2)
    if ( myid == nproc-1 ) jdxp1 = min0(j+1,jendx)
#endif
    do k = 1 , kz
      do i = 1 , iym1
        idx = max0(i,2)
        idxp1 = min0(i+1,iym1)
        atms%ubx3d(i,k,j) = d_rfour* & 
            (atm2%u(idx,k,jdx)+atm2%u(idxp1,k,jdx)+ &
             atm2%u(idx,k,jdxp1)+atm2%u(idxp1,k,jdxp1))/sps2%ps(i,j)
        atms%vbx3d(i,k,j) = d_rfour* &
            (atm2%v(idx,k,jdx)+atm2%v(idxp1,k,jdx)+ &
             atm2%v(idx,k,jdxp1)+atm2%v(idxp1,k,jdxp1))/sps2%ps(i,j)
      end do
    end do
  end do
 
  do j = 1 , jendl
    do k = 1 , kz
      do i = 1 , iy
        atms%ubd3d(i,k,j) = atm2%u(i,k,j)/sps2%pdot(i,j)
        atms%vbd3d(i,k,j) = atm2%v(i,k,j)/sps2%pdot(i,j)
      end do
    end do
  end do
 
  do j = jbegin , jendx
    do k = 1 , kz
      do i = 2 , iym1
        pl = a(k)*sps2%ps(i,j) + ptop
        thcon = ((sps2%ps(i,j)+ptop)/pl)**rovcp
        atms%pb3d(i,k,j) = pl
        atms%thx3d(i,k,j) = atms%tb3d(i,k,j)*thcon
      end do
    end do
  end do
 
!-----compute the height at full (za) and half (zq) sigma levels:
  do j = jbegin , jendx
    do i = 2 , iym1
      zq(i,kzp1) = d_zero
    end do
    do kk = 1 , kz
      k = kzp1 - kk
      do i = 2 , iym1
        cell = ptop/sps2%ps(i,j)
        zq(i,k) = zq(i,k+1) + rovg*atms%tb3d(i,k,j) *  &
                  dlog((sigma(k+1)+cell)/(sigma(k)+cell))
      end do
    end do
!
    do k = 1 , kz
      do i = 2 , iym1
        za(i,k,j) = d_half*(zq(i,k)+zq(i,k+1))
        dzq(i,k,j) = zq(i,k) - zq(i,k+1)
      end do
    end do
 
!-----Calculate the relative humidity and air density

    do i = 2 , iym1
      psrf = (sps2%ps(i,j)+ptop)*d_1000
      tv = atms%tb3d(i,kz,j)
      rhox2d(i,j) = psrf/(rgas*tv)
    end do
    do k = 1 , kz
      do i = 2 , iym2
        pres = (a(k)*sps2%ps(i,j)+ptop)*d_1000
        atms%rhob3d(i,k,j) = pres/(rgas*atms%tb3d(i,k,j)) !air density
        if ( atms%tb3d(i,k,j) > tzero ) then
          satvp = svp1*d_1000*dexp(svp2*(atms%tb3d(i,k,j)-tzero)           &
                & /(atms%tb3d(i,k,j)-svp3))
        else
          satvp = svp4*d_1000*dexp(svp5-svp6/atms%tb3d(i,k,j))
        end if
        atms%qsb3d(i,k,j) = ep2*satvp/(pres-satvp)
        rh = d_zero
        if ( atms%qsb3d(i,k,j) > d_zero ) then
          rh = atms%qvb3d(i,k,j)/atms%qsb3d(i,k,j)
        end if
        atms%rhb3d(i,k,j) = rh
      end do
    end do
  end do
 
  end subroutine mkslice
!
end module mod_slice

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
  use mod_che_interface
  use mod_pbl_interface
!
  private
!
  public :: mkslice
!
  contains 
!
  subroutine mkslice(jstart,jend,istart,iend)
  implicit none
  integer , intent(in) :: jstart , jend , istart , iend
!
  real(8) :: cell , pl , pres , psrf , rh , satvp , thcon , tv
  integer :: i , idx , idxp1 , j , jdx , jdxp1 , k , kk , n
 
  do k = 1 , kz
    do i = istart , iend
      do j = jstart , jend
        atms%tb3d(j,i,k) = atm2%t(i,k,j)/sps2%ps(j,i)
        atms%qvb3d(j,i,k) = atm2%qv(i,k,j)/sps2%ps(j,i)
        atms%qcb3d(j,i,k) = atm2%qc(i,k,j)/sps2%ps(j,i)
        if ( ichem == 1 ) then
          do n = 1 , ntr
            atms%chib3d(j,i,k,n) = chib(i,k,j,n)/sps2%ps(j,i)
          end do
        end if
      end do
    end do
  end do
  do k = 1 , kz
    do i = istart , iend
      do j = jbegin , jend
        jdx = j
        jdxp1 = j + 1
#ifndef BAND
        if ( myid == 0 ) jdx = max0(jdx,jbegin)
        if ( myid == nproc-1 ) jdxp1 = min0(jdxp1,jend)
#endif
        idx = max0(i,istart)
        idxp1 = min0(i+1,iend)
        atms%ubx3d(j,i,k) = d_rfour* & 
            (atm2%u(idx,k,jdx)+atm2%u(idxp1,k,jdx)+ &
             atm2%u(idx,k,jdxp1)+atm2%u(idxp1,k,jdxp1))/sps2%ps(j,i)
        atms%vbx3d(j,i,k) = d_rfour* &
            (atm2%v(idx,k,jdx)+atm2%v(idxp1,k,jdx)+ &
             atm2%v(idx,k,jdxp1)+atm2%v(idxp1,k,jdxp1))/sps2%ps(j,i)
      end do
    end do
  end do
 
  do k = 1 , kz
    do i = istart , iend+1
      do j = jstart , jend
        atms%ubd3d(j,i,k) = atm2%u(i,k,j)/sps2%pdot(j,i)
        atms%vbd3d(j,i,k) = atm2%v(i,k,j)/sps2%pdot(j,i)
      end do
    end do
  end do
 
  do k = 1 , kz
    do i = istart , iend
      do j = jstart , jend
        pl = a(k)*sps2%ps(j,i) + ptop
        thcon = ((sps2%ps(j,i)+ptop)/pl)**rovcp
        atms%pb3d(j,i,k) = pl
        atms%thx3d(j,i,k) = atms%tb3d(j,i,k)*thcon
      end do
    end do
  end do
 
!-----compute the height at full (za) and half (zq) sigma levels:
  do i = istart , iend
    do j = jbegin , jend
      zq(j,i,kzp1) = d_zero
    end do
  end do
  do kk = 1 , kz
    k = kzp1 - kk
    do i = istart , iend
      do j = jbegin , jend
        cell = ptop/sps2%ps(j,i)
        zq(j,i,k) = zq(j,i,k+1) + rovg*atms%tb3d(j,i,k) *  &
                  dlog((sigma(k+1)+cell)/(sigma(k)+cell))
      end do
    end do
  end do
!
  do k = 1 , kz
    do j = jbegin , jend
      do i = istart , iend
        za(j,i,k) = d_half*(zq(j,i,k)+zq(j,i,k+1))
        dzq(j,i,k) = zq(j,i,k) - zq(j,i,k+1)
      end do
    end do
  end do
 
!-----Calculate the relative humidity and air density

  do i = istart , iend
    do j = jbegin , jend
      psrf = (sps2%ps(j,i)+ptop)*d_1000
      tv = atms%tb3d(j,i,kz)
      rhox2d(j,i) = psrf/(rgas*tv)
    end do
  end do
  do k = 1 , kz
    do i = istart , iend
      do j = jstart , jend
        pres = (a(k)*sps2%ps(j,i)+ptop)*d_1000
        atms%rhob3d(j,i,k) = pres/(rgas*atms%tb3d(j,i,k)) !air density
        if ( atms%tb3d(j,i,k) > tzero ) then
          satvp = svp1*d_1000*dexp(svp2*(atms%tb3d(j,i,k)-tzero)           &
                & /(atms%tb3d(j,i,k)-svp3))
        else
          satvp = svp4*d_1000*dexp(svp5-svp6/atms%tb3d(j,i,k))
        end if
        atms%qsb3d(j,i,k) = ep2*satvp/(pres-satvp)
        rh = d_zero
        if ( atms%qsb3d(j,i,k) > d_zero ) then
          rh = atms%qvb3d(j,i,k)/atms%qsb3d(j,i,k)
        end if
        atms%rhb3d(j,i,k) = rh
      end do
    end do
  end do
 
  if ( ibltyp == 2 .or. ibltyp == 99 ) then
    do k = 1 , kzp1
      do i = istart , iend
        do j = jstart , jend
          atms%tkeb3d(j,i,k) = atm2%tke(i,k,j)
        end do
      end do
    end do
  end if
  end subroutine mkslice
!
end module mod_slice

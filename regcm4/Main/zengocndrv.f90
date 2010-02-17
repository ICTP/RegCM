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
!
      subroutine zengocndrv(j , k , ng , istart , iend)
!
      use mod_interfaces
      use mod_param1 , only : dtmin
      use mod_param2 , only : kbats
      use mod_param3 , only : ptop
      use mod_main , only : tgb , psb , zpbl
      use mod_pbldim , only : rhox2d , za
      use mod_slice , only : tb3d , ubx3d , vbx3d , qvb3d
      use mod_bats , only : tgb1d , ocld2d , tgb2d , tg1d , sent1d ,    &
                   & evpr1d , drag1d , u10m1d , v10m1d , t2m_1d , q2m_1d
      use mod_date , only : jyear , jyearr , ntime , ktau , ktaur
      use mod_constants , only : wlhv , tmelt

      implicit none

!
! Dummy arguments
!
      integer , intent (in) :: j , ng , istart , iend , k
!
! Local variables
!
      real(kind=8) :: dqh , dth , facttq , lh , psurf , q995 , qs , sh ,&
               & t995 , tau , tsurf , ustar , uv10 , uv995 , z995 , zi ,&
               & zo
      integer :: i , n
!
      do i = istart , iend
        do n = 1 , ng
          if ( ocld2d(n,i,j).lt.0.5 ) then
            uv995 = sqrt(ubx3d(i,k,j)**2+vbx3d(i,k,j)**2)
            tsurf = tgb(i,j) - tmelt
            t995 = tb3d(i,k,j) - tmelt
            q995 = qvb3d(i,k,j)/(1.+qvb3d(i,k,j))
            z995 = za(i,k,j)
            zi = zpbl(i,j)
            psurf = (psb(i,j)+ptop)*10.
            call zengocn(uv995,tsurf,t995,q995,z995,zi,psurf,qs,        &
                       & uv10,tau,lh,sh,dth,dqh,ustar,zo)
            tg1d(n,i) = tgb(i,j)
            tgb1d(n,i) = tgb(i,j)
            sent1d(n,i) = sh
            evpr1d(n,i) = lh/wlhv
!           Back out Drag Coefficient
            drag1d(n,i) = ustar**2*rhox2d(i,j)/uv995
            facttq = dlog(z995/2.)/dlog(z995/zo)
            u10m1d(n,i) = ubx3d(i,k,j)*uv10/uv995
            v10m1d(n,i) = vbx3d(i,k,j)*uv10/uv995
            t2m_1d(n,i) = t995 + tmelt - dth*facttq
!
            if ( mod(ntime+nint(dtmin*60.),kbats).eq.0 .or.             &
               & (jyear.eq.jyearr .and. ktau.eq.ktaur) ) then
              facttq = dlog(z995/2.)/dlog(z995/zo)
              q2m_1d(n,i) = q995 - dqh*facttq
              tgb2d(n,i,j) = tgb(i,j)
            end if
          end if
        end do
      end do
!
      end subroutine zengocndrv

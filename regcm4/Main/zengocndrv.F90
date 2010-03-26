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
      subroutine zengocndrv(j , ng , istart , iend , k)
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
      use mod_constants , only : wlhv , tzero

#ifdef CLM
      use mod_regcm_param , only : myid , jxp
      use clm_varsur , only : landmask
#endif
#ifdef DCSST
      use mod_param1 , only : dtbat
      use mod_bats , only : deltas , tdeltas , firstcall , dtskin ,     &
                 &          fsw2d , flw2d
      use mod_constants , only : sigm , rhoh2o , cpw0 , emsw , gti ,    &
                 &               vonkar
#endif
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
#ifdef CLM
      integer :: jj
#endif
#ifdef DCSST
!     Implement Zeng and Beljaars, GRL , 2005, ZB2005
!     Account for SST diurnal evoluation warm layer/ skin temperature
!     scheme
!     real(8) :: lwds , lwus
      real(8) :: rs , rd , td , tdelta , delta
      real(8) :: q , ustarw , fd , l , phidl , aa , bb , cc , lamb
      real(8) :: dtstend , dts , fs , tskin , dtsst
      real(8) , parameter :: a1 = 0.28D+00
      real(8) , parameter :: a2 = 0.27D+00
      real(8) , parameter :: a3 = 0.45D+00
      real(8) , parameter :: b1 = 71.5D+00
      real(8) , parameter :: b2 = 2.8D+00
      real(8) , parameter :: b3 = 0.07D+00
      real(8) , parameter :: alphaw = 0.207D-06
      real(8) , parameter :: nuw = 1.004D-06
      real(8) , parameter :: kw = 0.60
      real(8) , parameter :: nu = 0.3
      real(8) , parameter :: d = 3 ! reference depth for bulk SST
#endif
!
#ifdef CLM
      jj = (jxp*myid) + j
#endif
      do i = istart , iend
        do n = 1 , ng
#ifdef CLM
          if ( ocld2d(n,i,j).lt.0.5 .or. landmask(jj,i).eq.3 ) then
#else
          if ( ocld2d(n,i,j).lt.0.5 ) then
#endif
            uv995 = sqrt(ubx3d(i,k,j)**2+vbx3d(i,k,j)**2)
            tsurf = tgb(i,j) - tzero
            t995 = tb3d(i,k,j) - tzero
            q995 = qvb3d(i,k,j)/(1.+qvb3d(i,k,j))
            z995 = za(i,k,j)
            zi = zpbl(i,j)
            psurf = (psb(i,j)+ptop)*10.
            call zengocn(uv995,tsurf,t995,q995,z995,zi,psurf,qs,        &
                       & uv10,tau,lh,sh,dth,dqh,ustar,zo)
#ifdef DCSST
!           time step considered for the integration pof pronostic skin
!           temperature , equal to BATS time step
            dtsst = dtbat
!           handle the first call of the scheme
            if ( .not.firstcall(i,j) ) then
              deltas(i,j) = 0.001
              tdeltas(i,j) = tgb(i,j) - 0.001
              firstcall(i,j) = .true.
              td = tdeltas(i,j)
            end if
!           deep impact of aod on sst
!           if ( sum(aerext(i,:,j)).le.1 ) then
!             td = ts1(i,j) - sum(aerext(i,:,j))*0.8
!           else if ( sum(aerext(i,:,j)).gt.1 ) then
!             td = ts1(i,j)- 1.*0.8
!           end if
!
!           rs is the net surface sw flux (sw energy absorbed)
            rs = fsw2d(i,j)
!           rd is sw flux at 3m
            rd = rs*(a1*dexp(-d*b1) + a2*dexp(-d*b2) + a3*dexp(-d*b3))
!           ustar water (with air density ==1)
            ustarw = 0.5*ustar*(rhox2d(i,j)/rhoh2o)**0.5
!           lwds =  flwd2d(i,j)
!           lwus =  emsw*sigm*(tsurf+273.16)**4
!           q is the skin cooling term inckude net lw flux from
!           the radiative scheme
!           q = -(lh+sh+(lwus-lwds))
            q = -(lh+sh+flw2d(i,j))
!           fraction of solar radiation abosrbed in the sublayer
            fs = 0.065+11.*delta-(6.6e-5/delta)*(1-dexp(-delta/8.e-4))
!           dts= temperature difference between bulk level and skin level
!                determined from previous time step (via tdelta and td)
            dts = tdelta-td
!           m.o lenght calculation
            if ( dts.gt.0 ) then
              fd = (nu*gti*alphaw/(5*d))**0.5*                         &
                  &  rhoh2o*cpw0*ustarw**2*dts**0.5
            else
              fd = gti*alphaw*(q+rs-rd)
            end if
            l = rhoh2o*cpw0*ustarw**3/(vonkar*fd)
!           calulation of phidl (stability function)
            if ( (d/l).ge.0 ) then
              phidl = 1+5.*(d/l)
            else
              phidl = (1-16.*(d/l))**-0.5
            end if
!           prognostic evolution of dts
!           we can split the tendencies ddts/dt = a - b * dts
!           with a and b are ultimately function of dts through q
            aa = (q + rs - rd) / (d * cpw0 * rhoh2o * nu/(nu+1))
            bb = (nu+1) * vonkar * ustarw / (d*phidl)
!           exponential solution
            dtstend = aa - dts*(1-dexp(-bb*dtsst))/dtsst
!           update dts
            dts = dts + dtstend * dtsst
!           update tdelta
            tdelta = dts + td
!           update delta thickness  and cool skin tempearture
            aa = -16.*gti*alphaw*rhoh2o*cpw0*nuw**3./                   &
                &     (ustarw**4. *kw**2.)
            bb =  aa *(q+rs*fs)
            if ( bb.gt.0 ) then
!             case of cool skin layer correction
              cc= bb**(3./4.)
              lamb=6.*( (1.+(aa*(q+rs*fs))**0.75)**-0.333)
              delta = lamb*nuw/ustarw
              tskin= delta/(rhoh2o*cpw0*kw)*(q+rs*fs) + tdelta
            else
!             no cool skin layer in this case, tskin = warm layer
!             temperature
              tskin=tdelta
            end if
!           save the temperature difference and skin layer thickness
!           for next time step
            deltas(i,j) = delta
            tdeltas(i,j) = tdelta
            dtskin(i,j) = tskin-td
!           now feedback tskin in surface variable
            tgb(i,j) = tskin
#endif
            tg1d(n,i) = tgb(i,j)
            tgb1d(n,i) = tgb(i,j)
            sent1d(n,i) = sh
            evpr1d(n,i) = lh/wlhv
!           Back out Drag Coefficient
            drag1d(n,i) = ustar**2*rhox2d(i,j)/uv995
            facttq = dlog(z995/2.)/dlog(z995/zo)
            u10m1d(n,i) = ubx3d(i,k,j)*uv10/uv995
            v10m1d(n,i) = vbx3d(i,k,j)*uv10/uv995
            t2m_1d(n,i) = t995 + tzero - dth*facttq
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

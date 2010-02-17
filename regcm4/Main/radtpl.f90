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
 
      subroutine radtpl(tnm,ts,qnm,pnm,plh2o,tplnka,s2c,s2t,w,tplnke,   &
                      & tint,tint4,tlayr,tlayr4,pmln,piln)
!-----------------------------------------------------------------------
!
! Compute temperatures and path lengths for longwave radiation
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      L. Buja, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      use mod_regcm_param
      use mod_constants , only : ep2 , rga , sslp , rgsslp
      implicit none
!
!     Input arguments
!
! tnm    - Model level temperatures
! ts     - Surface skin temperature
! qnm    - Model level specific humidity
! pnm    - Pressure at model interfaces (dynes/cm2)
! plh2o  - Pressure weighted h2o path
!
!     Output arguments
!
! tplnka - Level temperature from interface temperatures
! s2c    - H2o continuum path length
! s2t    - H2o tmp and prs wghtd path length
! w      - H2o path length
! tplnke - Equal to tplnka
! tint   - Layer interface temperature
! tint4  - Tint to the 4th power
! tlayr  - K-1 level temperature
! tlayr4 - Tlayr to the 4th power
! pmln   - Ln(pmidm1)
! piln   - Ln(pintm1)
!
!
! Dummy arguments
!
      real(8) , dimension(ixm1,kxp1) :: piln , plh2o , pnm , s2c ,    &
           & s2t , tint , tint4 , tlayr , tlayr4 , tplnka , w
      real(8) , dimension(ixm1,kx) :: pmln , qnm , tnm
      real(8) , dimension(ixm1) :: tplnke , ts
      intent (in) piln , plh2o , pmln , pnm , qnm , tnm , ts
      intent (out) tint4 , tplnke
      intent (inout) s2c , s2t , tint , tlayr , tlayr4 , tplnka , w
!
!---------------------------Local variables-----------------------------
!
! i      - Longitude index
! k      - Level index
! r296   - Inverse stand temp for h2o continuum
! repsil - Inver ratio mol weight h2o to dry air
! dy     - Thickness of layer for tmp interp
! dpnm   - Pressure thickness of layer
! dpnmsq - Prs squared difference across layer
! rtnm   - Inverse level temperature
!
!-----------------------------------------------------------------------
!
! Local variables
!
      real(8) :: dpnm , dpnmsq , dy , r296 , repsil , rtnm
      integer :: i , k
!
      r296 = 1./296.D0
      repsil = 1./ep2
!
!     Set the top and bottom intermediate level temperatures,
!     top level planck temperature and top layer temp**4.
!
!     Tint is lower interface temperature
!     (not available for bottom layer, so use ground temperature)
!
      do i = 1 , ixm1
        tint(i,kxp1) = ts(i)
        tint4(i,kxp1) = tint(i,kx + 1)**4
        tplnka(i,1) = tnm(i,1)
        tint(i,1) = tplnka(i,1)
        tlayr4(i,1) = tplnka(i,1)**4
        tint4(i,1) = tlayr4(i,1)
      end do
!
!     Intermediate level temperatures are computed using temperature
!     at the full level below less dy*delta t,between the full level
!
      do k = 2 , kx
        do i = 1 , ixm1
          dy = (piln(i,k)-pmln(i,k))/(pmln(i,k-1)-pmln(i,k))
          tint(i,k) = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
          tint4(i,k) = tint(i,k)**4
        end do
      end do
!
!     Now set the layer temp=full level temperatures and establish a
!     planck temperature for absorption (tplnka) which is the average
!     the intermediate level temperatures.  Note that tplnka is not
!     equal to the full level temperatures.
!
      do k = 2 , kxp1
        do i = 1 , ixm1
          tlayr(i,k) = tnm(i,k-1)
          tlayr4(i,k) = tlayr(i,k)**4
          tplnka(i,k) = .5*(tint(i,k)+tint(i,k-1))
        end do
      end do
!
!     Calculate tplank for emissivity calculation.
!     Assume isothermal tplnke i.e. all levels=ttop.
!
      do i = 1 , ixm1
        tplnke(i) = tplnka(i,1)
        tlayr(i,1) = tint(i,1)
      end do
!
!     Now compute h2o path fields:
!
      do i = 1 , ixm1
        s2t(i,1) = plh2o(i,1)*tnm(i,1)
!       ccm3.2
!       w(i,1)   = (plh2o(i,1)*2.) / pnm(i,1)
!       s2c(i,1) = plh2o(i,1) * qnm(i,1) * repsil
 
!       ccm3.6.6
        w(i,1) = sslp*(plh2o(i,1)*2.)/pnm(i,1)
        rtnm = 1./tnm(i,1)
        s2c(i,1) = plh2o(i,1)*exp(1800.*(rtnm-r296))*qnm(i,1)*repsil
      end do
      do k = 1 , kx
        do i = 1 , ixm1
          dpnm = pnm(i,k+1) - pnm(i,k)
          dpnmsq = pnm(i,k+1)**2 - pnm(i,k)**2
          rtnm = 1./tnm(i,k)
          s2t(i,k+1) = s2t(i,k) + rgsslp*dpnmsq*qnm(i,k)*tnm(i,k)
          w(i,k+1) = w(i,k) + rga*qnm(i,k)*dpnm
          s2c(i,k+1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)                &
                     & *dexp(1800.*(rtnm-r296))*qnm(i,k)*repsil
        end do
      end do
!
      end subroutine radtpl

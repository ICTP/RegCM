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
 
      subroutine radinp(pmid,pint,h2ommr,cld,o3vmr,pmidrd,pintrd,plco2, &
                      & plh2o,tclrsf,eccf,o3mmr)

!-----------------------------------------------------------------------
!
! Set latitude and time dependent arrays for input to solar
! and longwave radiation.
!
! Convert model pressures to cgs, compute path length arrays needed for the
! longwave radiation, and compute ozone mixing ratio, needed for the solar
! radiation.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      use mod_regcm_param
      use mod_comtim
      use mod_constants , only : gtigts , twopi , dayspy , rgsslp
      use mod_crdcae , only : co2vmr
#ifdef CLM
      use mod_clm
#endif
      implicit none
!
! Dummy arguments
!
      real(8) :: eccf
      real(8) , dimension(ixm1,kxp1) :: cld , pint , pintrd , plco2 ,   &
           & plh2o , tclrsf
      real(8) , dimension(ixm1,kx) :: h2ommr , o3mmr , o3vmr , pmid ,   &
           & pmidrd
      intent (in) cld , h2ommr , o3vmr , pint , pmid
      intent (out) eccf , o3mmr , plco2 , pmidrd
      intent (inout) pintrd , plh2o , tclrsf
!
! Local variables
!
#ifndef CLM
      real(8) :: theta
#endif
      real(8) :: amco2 , amd , amo , cpwpl , p0 , vmmr
      integer :: i , k
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! pmid    - Pressure at model mid-levels (pascals)
! pint    - Pressure at model interfaces (pascals)
! h2ommr  - H2o mass mixing ratio
! cld     - Fractional cloud cover
! o3vmr   - ozone volume mixing ratio
!
!     Output arguments
!
! pmidrd  - Pressure at mid-levels (dynes/cm*2)
! pintrd  - Pressure at interfaces (dynes/cm*2)
! plco2   - Vert. pth lngth of co2 (prs-weighted)
! plh2o   - Vert. pth lngth h2o vap.(prs-weighted)
! tclrsf  - Product of clr-sky fractions from top of atmosphere to level.
! eccf    - Earth-sun distance factor
! o3mmr   - Ozone mass mixing ratio
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude loop index
! k       - Vertical loop index
! theta   - Earth orbit seasonal angle in radians
! p0      - Standard pressure (dynes/cm**2)
! amd     - Effective molecular weight of dry air (g/mol)
! amo     - Molecular weight of ozone (g/mol)
! amco2   - Molecular weight of co2   (g/mol)
! cpwpl   - Const in co2 mixing ratio to path length conversn
! vmmr    - Ozone volume mixing ratio
!
      data p0/1.01325E6/
      data amd/28.9644/
      data amo/48.0000/
      data amco2/44.0000/
!
!-----------------------------------------------------------------------
!
!     Compute solar distance factor and cosine solar zenith angle usi
!     day value where a round day (such as 213.0) refers to 0z at
!     Greenwich longitude.
!
!     Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
!     Processes in Meterology and Climatology, Elsevier Scientific
!     Publishing Company, New York  p. 57, p. 62,63.
!
!     Compute eccentricity factor (sun-earth distance factor)
!
#ifdef CLM
      eccf  = r2ceccf
#else
      theta = twopi*calday/dayspy
      eccf = 1.000110 + .034221*dcos(theta) + .001280*dsin(theta)       &
           & + .000719*dcos(2.*theta) + .000077*dsin(2.*theta)
#endif
!
!     Convert pressure from pascals to dynes/cm2
!
      do k = 1 , kx
        do i = 1 , ixm1
          pmidrd(i,k) = pmid(i,k)*10.0
          pintrd(i,k) = pint(i,k)*10.0
        end do
      end do
      do i = 1 , ixm1
        pintrd(i,kxp1) = pint(i,kx + 1)*10.0
      end do
!
!     Compute path quantities used in the longwave radiation:
!
      vmmr = amco2/amd
      cpwpl = vmmr*0.5/(gtigts*p0)
      do i = 1 , ixm1
        plh2o(i,1) = rgsslp*h2ommr(i,1)*pintrd(i,1)*pintrd(i,1)
        plco2(i,1) = co2vmr*cpwpl*pintrd(i,1)*pintrd(i,1)
        tclrsf(i,1) = 1.
      end do
      do k = 1 , kx
        do i = 1 , ixm1
          plh2o(i,k+1) = plh2o(i,k)                                     &
                       & + rgsslp*(pintrd(i,k+1)**2-pintrd(i,k)**2)     &
                       & *h2ommr(i,k)
          plco2(i,k+1) = co2vmr*cpwpl*pintrd(i,k+1)**2
          tclrsf(i,k+1) = tclrsf(i,k)*(1.-cld(i,k+1))
        end do
      end do
!
!     Convert ozone volume mixing ratio to mass mixing ratio:
!
      vmmr = amo/amd
      do k = 1 , kx
        do i = 1 , ixm1
          o3mmr(i,k) = vmmr*o3vmr(i,k)
        end do
      end do
!
      end subroutine radinp

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
 
      subroutine aermix(pint , rh , j , istart , iend , nx , nk , ntrac)
 
!-----------------------------------------------------------------------
! Set global mean tropospheric aerosol
!
! Specify aerosol mixing ratio and compute relative humidity for later
! adjustment of aerosol optical properties. Aerosol mass mixing ratio
! is specified so that the column visible aerosol optical depth is a
! specified global number (tauvis). This means that the actual mixing
! ratio depends on pressure thickness of the lowest three atmospheric
! layers near the surface.
!
! Optical properties and relative humidity parameterization are from:
!
! J.T. Kiehl and B.P. Briegleb  "The Relative Roles of Sulfate Aerosols
! and Greenhouse Gases in Climate Forcing"  Science  260  pp311-314
! 16 April 1993
!
! Visible (vis) here means 0.5-0.7 micro-meters
! Forward scattering fraction is taken as asymmetry parameter squared
!
!---------------------------Code history--------------------------------
!
! Original version:  B. Briegleb  March 1995
! Standarized:       L. Buja,     Feb 1996
! Reviewed:          B. Briegleb, Mar 1996
!
!-----------------------------------------------------------------------
 
      use mod_param2 , only : ichem
      use mod_slice
      use mod_main
      use mod_mainchem
      use mod_aerosol , only : aermmb , aermmr
      use mod_constants , only : gtigts
      implicit none
!
! Dummy arguments
!
      integer :: j , istart , iend , nx , nk , ntrac
!     Radiation level interface pressures (dynes/cm2)
      real(8) , dimension(nx,nk+1) :: pint
!     Radiation level relative humidity (fraction)
      real(8) , dimension(nx,nk) :: rh
!
      intent (in) j , pint , istart , iend , nk , ntrac
      intent (out) rh
!
!---------------------------Local variables-----------------------------
!
! i      - longitude index
! k      - level index
! mxaerl - max nmbr aerosol levels counting up from surface
! tauvis - visible optical depth
! kaervs - visible extinction coefficiant of aerosol (m2/g)
! omgvis - visible omega0
! gvis   - visible forward scattering asymmetry parameter
!
!-----------------------------------------------------------------------
 
      real(8) :: gvis , kaervs , omgvis , rhfac , tauvis
      integer :: i , itr , k , mxaerl
!
      data kaervs/5.3012D0/         ! multiplication factor for kaer
      data omgvis/0.999999D0/
      data gvis/0.694889D0/
      data rhfac/1.6718D0/          ! EES added for efficiency
!
!--------------------------------------------------------------------------
!
      mxaerl = 4
!
!fil  tauvis = 0.01D0
      tauvis = 0.04D0
!
!     Set relative humidity and factor; then aerosol amount:
!
      do k = 1 , nk
        do i = istart , iend
 
!added    July 13, 2000: needed for aerosols in radiation
          rh(i,k) = dmin1(rhb3d(i,k,j),0.99D0)
!EES:     do not change to 1.00:  wscoef(3,10) in radcsw = .9924 and is
!         divided by RH.  rh is limited to .99 to avoid dividing by zero
!added
 
!
!         Define background aerosol
!         Find constant aerosol mass mixing ratio for specified levels
!         in the column, converting units where appropriate
!         for the moment no more used
!
          if ( k.ge.nk + 1 - mxaerl ) then
            aermmb(i,k) = gtigts*tauvis/(1.0D4*kaervs*rhfac*(1.-omgvis* &
                        & gvis*gvis)                                    &
                        & *(pint(i,kxp1)-pint(i,kx + 1 - mxaerl)))
          else
            aermmb(i,k) = 0.0D0
          end if
 
!         if(ichem .eq. 1 .and. idirect .eq. 1) then
          if ( ichem.eq.1 ) then
            do itr = 1 , ntrac
!             aermmr(i,k,itr)= dmax1( chia(i,k,j,itr)/psa(i,j)
!             $                               ,aermmb(i,k) )
              aermmr(i,k,itr) = chia(i,k,j,itr)/psa(i,j)
            end do
          else if ( ehso4 ) then
            do itr = 1 , ntrac
!             aermmr(i,k,itr) = aermmb(i,k) + aermm(i,k,j)
!Dec.11       aermmr(i,k,itr) = aermm(i,k,j)
              aermmr(i,k,itr) = 0.0D0
!Dec.11_
            end do
          else
            do itr = 1 , ntrac
!             aermmr(i,k,itr)= 0.0D0
              aermmr(i,k,itr) = aermmb(i,k)
            end do
          end if
        end do
      end do
!
      end subroutine aermix

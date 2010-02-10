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
 
      subroutine cldefr(ioro,t,rel,rei,fice,ps,pmid)
!-----------------------------------------------------------------------
!
! Compute cloud drop size
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Kiehl, January 1993
!
!-----------------------------------------------------------------------
!
      use mod_regcm_param
      use mod_parrad
      implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! ioro   - nint(oro(i))
! t      - Temperature
! ps     - surface pressure
! pmid   - midpoint pressures
!
!     Output arguments
!
! rel    - liquid effective drop size (microns)
! rei    - ice effective drop size (microns)
! fice   - fractional ice content within cloud
! pirnge - nrmlzd pres range for ice particle changes
! picemn - normalized pressure below which rei=reimax
! rirnge - range of ice radii (reimax - 10 microns)
! reimax - maximum ice effective radius
! pnrml  - normalized pressure
! weight - coef. for determining rei as fn of P/PS
!
!
! Dummy arguments
!
      real(8) , dimension(plond,plev) :: fice , pmid , rei , rel , t
      integer , dimension(plond) :: ioro
      real(8) , dimension(plond) :: ps
      intent (in) ioro , pmid , ps , t
      intent (out) fice , rei , rel
!
!---------------------------Local workspace-----------------------------
!
! i, k   - longitude, level indices
! rliq   - temporary liquid drop size
!
!-----------------------------------------------------------------------
!
! Local variables
!
      integer :: i , k
      real(8) :: picemn , pirnge , pnrml , reimax , rirnge , rliq ,     &
               & weight
!
      do k = 1 , plev
        do i = 1 , plon
!
!         Define liquid drop size
!
          if ( ioro(i).ne.1 ) then
!
!           Effective liquid radius over ocean and sea ice
!
            rliq = 10.0
          else
!
!           Effective liquid radius over land
!
            rliq = 5.0 + 5.0*dmin1(1.D0,dmax1(0.D0,(263.16-t(i,k))*0.05)&
                 & )
          end if
!
          rel(i,k) = rliq
!fil
!         test radius = 10.0
!         rel(i,k) = 10.0
!fil
!+        rei(i,k) = 30.0
!
!         Determine rei as function of normalized pressure
!
          reimax = 30.0
          rirnge = 20.0
          pirnge = 0.4
          picemn = 0.4
!
          pnrml = pmid(i,k)/ps(i)
          weight = dmax1(dmin1((pnrml-picemn)/pirnge,1.D0),0.D0)
          rei(i,k) = reimax - rirnge*weight
!
!         Define fractional amount of cloud that is ice
!
!         if warmer than -10 degrees C then water phase
!
          if ( t(i,k).gt.263.16 ) fice(i,k) = 0.0
!
!         if colder than -10 degrees C but warmer than -30 C mixed phase
!
          if ( t(i,k).le.263.16 .and. t(i,k).ge.243.16 ) fice(i,k)      &
             & = (263.16-t(i,k))/20.0
!
!         if colder than -30 degrees C then ice phase
!
          if ( t(i,k).lt.243.16 ) fice(i,k) = 1.0
!
!         Turn off ice radiative properties by setting fice = 0.0
!
!fil      no-ice test
!         fice(i,k) = 0.0
!
        end do
      end do
!
      end subroutine cldefr

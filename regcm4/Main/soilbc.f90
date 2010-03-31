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
 
      subroutine soilbc

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   this subrout overwrites many of the soil constants
!   as a function of location(jlon,jlat)
!
      use mod_regcm_param
      use mod_bats
      implicit none
!
! Local variables
!
      real(8) :: ck , dmax , dmin , dmnor , phi0 , tweak1
      integer :: itex , n , i
!
!     ================================================================
!     new soils data as a fn of texture make porosity, soil suction,
!     hydraul conduc, wilting frac variables rather than consts
!
!     for explanation of params set here see block data
!
!     relfc is the ratio of field capacity to saturated water content,
!     defined so the rate of gravitational drainage at field
!     capacity is assumed to be 2 mm/day (baver et al., 1972)
!     ===============================================================
!
 
      do i = 2 , iym1
        do n = 1 , nnsg
 
          if ( lveg(n,i).ne.0 ) then
 
!           **********            lveg is set in subr. interf
            freza(lveg(n,i)) = 0.15*deprv(lveg(n,i))
            frezu(lveg(n,i)) = 0.15*depuv(lveg(n,i))
            itex = iexsol(lveg(n,i))
            texrat(n,i) = skrat(itex)
            porsl(n,i) = xmopor(itex)
            xkmx(n,i) = xmohyd(itex)
            bsw(n,i) = bee(itex)
            bfc(n,i) = 5.8 - bsw(n,i)*(0.8+0.12*(bsw(n,i)-4.)*          &
                      & dlog10(1.E2*xkmx(n,i)))
 
            phi0 = xmosuc(itex)
            dmax = bsw(n,i)*phi0*xkmx(n,i)/porsl(n,i)
            dmin = 1.E-3
            dmnor = 1550.*dmin/dmax
            tweak1 = (bsw(n,i)*(bsw(n,i)-6.)+10.3)                      &
                   & /(bsw(n,i)*bsw(n,i)+40.*bsw(n,i))
            ck = (1.+dmnor)*tweak1*0.23/0.02356
            evmx0(n,i) = 1.02*dmax*ck/dsqrt(depuv(lveg(n,i))*           &
                        & deprv(lveg(n,i)))
            gwmx0(n,i) = depuv(lveg(n,i))*porsl(n,i)
            gwmx1(n,i) = deprv(lveg(n,i))*porsl(n,i)
            gwmx2(n,i) = deptv(lveg(n,i))*porsl(n,i)
            wiltr(n,i) = xmowil(itex)
!           **********            force irrigated crop to be at field
!           capacity
            relfc(n,i) = xmofc(itex)
 
          end if
 
        end do
      end do
 
      end subroutine soilbc

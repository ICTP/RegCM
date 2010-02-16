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
      integer :: itex , n , np
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
 
      do np = np1 , npts
        do n = 1 , nnsg
 
          if ( lveg(n,np).ne.0 ) then
 
!           **********            lveg is set in subr. interf
            freza(lveg(n,np)) = 0.15*deprv(lveg(n,np))
            frezu(lveg(n,np)) = 0.15*depuv(lveg(n,np))
            itex = iexsol(lveg(n,np))
            texrat(n,np) = skrat(itex)
            porsl(n,np) = xmopor(itex)
            xkmx(n,np) = xmohyd(itex)
            bsw(n,np) = bee(itex)
            bfc(n,np) = 5.8 - bsw(n,np)                                 &
                      & *(0.8+0.12*(bsw(n,np)-4.)*dlog10(1.E2*xkmx(n,np)&
                      & ))
 
            phi0 = xmosuc(itex)
            dmax = bsw(n,np)*phi0*xkmx(n,np)/porsl(n,np)
            dmin = 1.E-3
            dmnor = 1550.*dmin/dmax
            tweak1 = (bsw(n,np)*(bsw(n,np)-6.)+10.3)                    &
                   & /(bsw(n,np)*bsw(n,np)+40.*bsw(n,np))
            ck = (1.+dmnor)*tweak1*0.23/0.02356
            evmx0(n,np) = 1.02*dmax*ck/dsqrt(depuv(lveg(n,np))*deprv(   &
                        & lveg(n,np)))
            gwmx0(n,np) = depuv(lveg(n,np))*porsl(n,np)
            gwmx1(n,np) = deprv(lveg(n,np))*porsl(n,np)
            gwmx2(n,np) = deptv(lveg(n,np))*porsl(n,np)
            wiltr(n,np) = xmowil(itex)
!           **********            force irrigated crop to be at field
!           capacity
            relfc(n,np) = xmofc(itex)
 
          end if
 
        end do
      end do
 
      end subroutine soilbc

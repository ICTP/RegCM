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
 
      subroutine deriv

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     derivatives of energy fluxes with respect to leaf temperature for
!     newton-raphson calculation of leaf temperature.
!     input: rs,ra,cdrd,rppq,efe.    output: qsatld,dcd.
!
!     approximate by derivatives of cdr and ef.  many weaker
!     dependences on leaf temperature are omitted, as convergence
!     rate is not affected.
!
      use mod_regcm_param
      use mod_bats , only : npts , ldoc1d , sigf , tlef1d , c , cdr ,   &
                    & ts1d , wtga , qsatl , wta0 , wtgaq , tg1d , cn1 , &
                    & wtg0 , df
      use mod_ictp01
      implicit none
!
! Local variables
!
      real(8) :: dne , hfl , xkb
      integer :: n , np
!
      do np = 1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
              dne = 1./(tlef1d(n,np)-b(n,np))
              qsatld(n,np) = qsatl(n,np)*a(n,np)*(c(67)-b(n,np))*dne**2
              xkb = cdrd(n,np)/cdr(n,np)
              hfl = df(n,np)                                            &
                  & *(wtga(n,np)*tlef1d(n,np)-wtg0(n,np)*tg1d(n,np)     &
                  & -wta0(n,np)*ts1d(n,np))
              dcd(n,np) = cn1(n,np)*rppq(n,np)*wtgaq(n,np)*qsatld(n,np) &
                        & + (1.-wtgaq(n,np))*efe(n,np)                  &
                        & *xkb + (1.-wtga(n,np))*hfl*xkb
              dcd(n,np) = dmax1(dcd(n,np),0.D0)
              dcd(n,np) = dmin1(dcd(n,np),500.D0)
            end if
          end if
        end do
      end do
 
      end subroutine deriv

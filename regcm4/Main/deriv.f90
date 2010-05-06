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
      use mod_dynparam
      use mod_bats , only : wtg0 , ldoc1d , sigf , tlef1d , df , cdr ,  &
                    & ts1d , wtga , qsatl , wta0 , wtgaq , tg1d , cn1
      use mod_constants , only : tzero
      use mod_ictp01
      implicit none
!
! Local variables
!
      real(8) :: dne , hfl , xkb
      integer :: n , i
!
      do i = 1 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              dne = 1./(tlef1d(n,i)-b(n,i))
              qsatld(n,i) = qsatl(n,i)*a(n,i)*(tzero-b(n,i))*dne**2
              xkb = cdrd(n,i)/cdr(n,i)
              hfl = df(n,i)*(wtga(n,i)*tlef1d(n,i)-wtg0(n,i)*tg1d(n,i)  &
                  & -wta0(n,i)*ts1d(n,i))
              dcd(n,i) = cn1(n,i)*rppq(n,i)*wtgaq(n,i)*qsatld(n,i)      &
                        & + (1.-wtgaq(n,i))*efe(n,i)                    &
                        & *xkb + (1.-wtga(n,i))*hfl*xkb
              dcd(n,i) = dmax1(dcd(n,i),0.D0)
              dcd(n,i) = dmin1(dcd(n,i),500.D0)
            end if
          end if
        end do
      end do
 
      end subroutine deriv

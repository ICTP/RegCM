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
 
      module mod_interfaces

      interface
        subroutine pcp(j , istart , iend , nk)
          integer :: j , istart , iend , nk
        end subroutine pcp
      end interface

      interface
        subroutine zengocndrv(j , ng , istart , iend , k)
          integer :: j , ng , istart , iend , k
        end subroutine zengocndrv
      end interface

      interface
        subroutine zengocn(u , ts , t , q , hgt , zi , ps , qs , u10 ,  &
                           tau , alh , ash , dth , dqh , ustar , zo)
          real(kind=8) :: hgt , q , t , u , zi , ts , ps
          real(kind=8) :: alh , ash , tau , u10
          real(kind=8) :: dqh , dth , qs , ustar , zo
        end subroutine zengocn
      end interface

      interface
        subroutine vecbats(j , k , istart , iend , ng)
          integer :: j , k , istart , iend , ng
        end subroutine vecbats
      end interface

      interface
        subroutine interf(ivers , j , k , istart , iend , ng)
          integer :: ivers , j , k , istart , iend , ng
        end subroutine interf
      end interface

      end module mod_interfaces

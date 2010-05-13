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

      module mod_culookup
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: jptlucu1 = 50000 , jptlucu2 = 370000
!
      logical :: lookupoverflow
      real(8) , dimension(jptlucu1:jptlucu2) :: tlucua , tlucuaw ,      &
           & tlucub , tlucuc
!
      real(8) :: cmfcmax , cmfcmin , cmfctop , cmfdeps , cprcon ,       &
               & entrdd , entrmid , entrpen , entrscv , rhcdd
      logical :: lamip2 , lmfdd , lmfdudv , lmfmid , lmfpen , lmfscv

      contains

      subroutine set_lookup_tables
 
!-----------------------------------------------------------------------
! *mo_convect_tables* - tables for convective adjustment code
!
!           d.salmond  cray (uk)   12/8/91
!
!
! table lookups replaced
!
!           A. Rhodin mpi 12/98
!
! When replacing the table lookups the following code has been used:
!
!   tlucua :  c2es*exp(merge(c3les,c3ies,lo)*(tt-tzero) &
!           /      (tt-merge(c4les,c4ies,lo)))
!
!   tlucub :  merge(c5alvcp,c5alscp,lo) / (tt-merge(c4les,c4ies,lo))**2
!
!   tlucuc :  merge(wlhvocp, wlhsocp, lo)
!
!   tlucuaw:  c2es*exp(c3les*(tt-tzero)*(1./(tt-c4les)))
!
!   with:     lo = tt > tzero
!
! compile with option -DNOLOOKUP in order to replace lookup tables
!-----------------------------------------------------------------------
 
      use mod_constants , only :  tzero , rgas , rwat , c3les , c3ies , &
         & c4les , c4ies , c2es , c5alvcp , c5alscp , wlhvocp , wlhsocp
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: clavm1 = -6096.9385 , clavm2 = 21.2409642 ,&
                           & clavm3 = -2.711193 , clavm4 = 1.673952 ,   &
                           & clavm5 = 2.433502 , ciavm1 = -6024.5282 ,  &
                           & ciavm2 = 29.32707 , ciavm3 = 1.0613868 ,   &
                           & ciavm4 = -1.3198825 , ciavm5 = -0.49382577
!
! Local variables
!
      integer :: it
      logical :: lo
      real(8) :: tt

!----------------------------
!     -- Initialise lookup tables
!     called from 'setphys'
!----------------------------
 
      tt = jptlucu1*0.001
      do it = jptlucu1 , jptlucu2
        lo = tt.gt.tzero
 
        if ( lamip2 ) then
          tlucua(it) = exp((merge(clavm1,ciavm1,lo)/tt+merge(clavm2,    &
                     & ciavm2,lo)+merge(clavm3,ciavm3,lo)               &
                     & *tt/100.+merge(clavm4,ciavm4,lo)*tt*tt/1.E5+     &
                     & merge(clavm5,ciavm5,lo)*log(tt)))*rgas/rwat
        else
          tlucua(it) = c2es*exp(merge(c3les,c3ies,lo)*(tt-tzero)        &
                     & *(1./(tt-merge(c4les,c4ies,lo))))
        end if
        tlucub(it) = merge(c5alvcp,c5alscp,lo)                          &
                   & *(1./(tt-merge(c4les,c4ies,lo)))**2
        tlucuc(it) = merge(wlhvocp,wlhsocp,lo)
 
        tlucuaw(it) = c2es*exp(c3les*(tt-tzero)*(1./(tt-c4les)))
 
        tt = tt + 0.001
      end do
 
      end subroutine set_lookup_tables
!
      end module mod_culookup

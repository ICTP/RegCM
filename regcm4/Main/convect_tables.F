C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE set_lookup_tables

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
!   tlucua :  c2es*EXP(MERGE(c3les,c3ies,lo)*(tt-tmelt) &
!           /      (tt-MERGE(c4les,c4ies,lo)))
!
!   tlucub :  MERGE(c5alvcp,c5alscp,lo) / (tt-MERGE(c4les,c4ies,lo))**2
!
!   tlucuc :  MERGE(alvdcp, alsdcp, lo)
!
!   tlucuaw:  c2es*EXP(c3les*(tt-tmelt)*(1./(tt-c4les)))
!
!   with:     lo = tt > tmelt
!
! compile with option -DNOLOOKUP in order to replace lookup tables
!-----------------------------------------------------------------------

      IMPLICIT NONE
#     include "include/culookup.cb"
#     include "include/tiedparam0.cb"
#     include "include/tied_constants.cb"

!----------------------------
! -- Initialise lookup tables
!      called from 'setphys'
!----------------------------

      real(kind=8), PARAMETER :: clavm1 = -6096.9385
      real(kind=8), PARAMETER :: clavm2 =    21.2409642
      real(kind=8), PARAMETER :: clavm3 =    -2.711193
      real(kind=8), PARAMETER :: clavm4 =     1.673952
      real(kind=8), PARAMETER :: clavm5 =     2.433502

      real(kind=8), PARAMETER :: ciavm1 = -6024.5282
      real(kind=8), PARAMETER :: ciavm2 =    29.32707
      real(kind=8), PARAMETER :: ciavm3 =     1.0613868
      real(kind=8), PARAMETER :: ciavm4 =    -1.3198825
      real(kind=8), PARAMETER :: ciavm5 =    -0.49382577

      real(kind=8)  :: tt
      INTEGER :: it  
      LOGICAL :: lo
    
      tt = jptlucu1 * 0.001
      DO it = jptlucu1, jptlucu2
        lo = tt > tmelt
      
        IF (lamip2) THEN
          tlucua(it) = EXP((MERGE(clavm1, ciavm1, lo)/tt
     &                  + MERGE(clavm2, ciavm2, lo)
     &                  + MERGE(clavm3, ciavm3, lo)*tt/100.
     &                  + MERGE(clavm4, ciavm4, lo)*tt*tt/1.e5
     &                  + MERGE(clavm5, ciavm5, lo)*LOG(tt)))*rdti/rvti
        ELSE
          tlucua(it) = c2es*EXP(MERGE(c3les,  c3ies,  lo)*(tt-tmelt)
     &             *  (1./(tt-MERGE(c4les,  c4ies,  lo))))
        end if
        tlucub(it)  =          MERGE(c5alvcp,c5alscp,lo)
     &            *  (1./(tt-MERGE(c4les,  c4ies,  lo)))**2
        tlucuc(it)  =          MERGE(alvdcp, alsdcp, lo)

        tlucuaw(it) = c2es*EXP(c3les*(tt-tmelt)*(1./(tt-c4les)))

        tt = tt + 0.001
      end do

      END SUBROUTINE set_lookup_tables

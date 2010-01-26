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
C
C This is for test AER11D1_048x052
C
C Domain dimension parameters

      INTEGER IX
      INTEGER ILX, ILXM
      parameter(IX     =   52)
      parameter(ILX    =  ix-1)
      parameter(ILXM   =  ix-2)

#ifdef MPP1
      INTEGER NPROC
      INTEGER MJX
      INTEGER JXP
      integer JXBB
#else
      INTEGER JX
      integer JLX,JLXM
#endif

#ifdef MPP1
      parameter(NPROC  =    2)
      parameter(MJX    =   48)
      parameter(JXP    = MJX/NPROC)
      parameter(JXBB   = mjx-1)
#else
      parameter(JX     =   48)
      parameter(JLX    =  jx-1)
      parameter(JLXM   =  jx-2)
#endif

      INTEGER KX
      INTEGER KXM, KXP1, KXP2
      parameter(KX     =    18)
      parameter(kxm    =  kx-1)
      parameter(kxp1   =  kx+1)
      parameter(kxp2   =  kx+2)

      INTEGER NSG
      INTEGER NNSG
      parameter(NSG    =     1)
      parameter(NNSG   =     1)

      CHARACTER*5 DATTYP
      LOGICAL     EHSO4 
      CHARACTER*4 LSMTYP
      CHARACTER*7 AERTYP
      parameter(DATTYP='ERAIN')
      parameter(EHSO4 =.false.)
      parameter(LSMTYP='BATS')
      parameter(AERTYP='AER00D0')

      INTEGER IBYTE
      parameter(IBYTE  =     4)

      INTEGER debug_level
      parameter(debug_level = 1)

C
C Set Buffer Zone Depth
C

      integer nspgx     ! (nspgx-1) represents the number of
C                       ! cross-point slices on the boundary affected by
C                       ! sponge or relaxation boundary conditions. (=6)

      integer nspgd     ! (nspgd-1) represents the number of
C                       ! dot-point slices on the boundary affected by
C                       ! sponge or relaxation boundary conditions. (=6)

C     parameter (nspgx=12,nspgd=12) ! # of slices in boundary conditions
      parameter (nspgx=12,nspgd=12) ! # of slices in boundary conditions
C
      integer nspgv     ! dimension of the equivalence arrays.
C                       !  (= (nspgd+nspgx)*8+8 )

      integer nspgp     ! dimension of the equivalence arrays for
C                       ! boundary values and the pressure tendency.
C                       !  (= nspgx*4)
      parameter (nspgv=(nspgd+nspgx)*8+8,nspgp=nspgx*4)               

      integer nsplit
      parameter( nsplit = 2)           ! # of split exp. modes
      integer lkpts
      parameter( lkpts = 10)           !# of lake points for lake model

C Bats parameters

      integer nbmax
      integer np1
      parameter (nbmax=ix-1,np1=2) 

C Tracer parameter

      integer ntr                      ! tracer number
      parameter( ntr = 10)
      integer nbin                     ! bins number ( dust)
      parameter(nbin = 4)

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

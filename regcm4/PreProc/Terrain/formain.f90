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

      subroutine formain(iy,jx,kz,nsg,ibyte,dattyp,ehso4,lsmtyp,aertyp, &
                       & nproc)
      implicit none
!
! Dummy arguments
!
      character(7) :: aertyp
      character(5) :: dattyp
      logical :: ehso4
      integer :: ibyte , iy , jx , kz , nproc , nsg
      character(4) :: lsmtyp
      intent (in) aertyp , dattyp , ehso4 , ibyte , iy , jx , kz ,      &
                & lsmtyp , nproc , nsg
!
      open (23,file='../../Main/regcm.param')
      write (23,'(a)') '      INTEGER IX'
      write (23,'(a)') '      INTEGER JX'
      write (23,'(a)') '      INTEGER KX'
      write (23,'(a)') '      INTEGER NSG'
      write (23,'(a)') '      INTEGER NNSG'
      write (23,'(a)') '      INTEGER IBYTE'
      write (23,'(a)') '      CHARACTER*5 DATTYP'
      write (23,'(a)') '      LOGICAL     EHSO4 '
      write (23,'(a)') '      CHARACTER*4 LSMTYP'
      write (23,'(a)') '      CHARACTER*7 AERTYP'
      write (23,'(a)') '      integer jlx,jlxm'
      write (23,99001) 'IX     =' , iy
      write (23,99001) 'JX     =' , jx
      write (23,99001) 'KX     =' , kz
      write (23,99001) 'NSG    =' , nsg
      write (23,99001) 'NNSG   =' , nsg*nsg
      write (23,99001) 'IBYTE  =' , ibyte
      write (23,99002) "DATTYP='" , dattyp
      if ( dattyp=='EH5OM' .and. ehso4 ) then
        write (23,'(a)') "      parameter(EHSO4 =.true. )"
      else
        write (23,'(a)') "      parameter(EHSO4 =.false.)"
      end if
      write (23,99003) "LSMTYP='" , lsmtyp
      write (23,99004) "AERTYP='" , aertyp
      write (23,'(a)') '      parameter(jlx=jx-1,jlxm=jx-2)'
      close (23)
!
      if ( nproc>0 ) then
        open (23,file='../../Main/regcm.param2')
        write (23,'(a)') '      INTEGER IX'
        write (23,'(a)') '      INTEGER NPROC'
        write (23,'(a)') '      INTEGER MJX'
        write (23,'(a)') '      INTEGER KX'
        write (23,'(a)') '      INTEGER NSG'
        write (23,'(a)') '      INTEGER NNSG'
        write (23,'(a)') '      INTEGER IBYTE'
        write (23,'(a)') '      INTEGER JXP'
        write (23,'(a)') '      CHARACTER*5 DATTYP'
        write (23,'(a)') '      LOGICAL     EHSO4 '
        write (23,'(a)') '      CHARACTER*4 LSMTYP'
        write (23,'(a)') '      CHARACTER*7 AERTYP'
        write (23,'(a)') '      integer jxbb'
        write (23,99001) 'IX     =' , iy
        write (23,99001) 'NPROC  =' , nproc
        write (23,99001) 'MJX    =' , jx
        if ( mod(jx,nproc)==0 ) then
          write (23,99005) 'JXP    = MJX/NPROC'
        else
          write (*,*) 'The present parallel code requirea that'
          write (*,*) 'MJX can be divided by NPROC'
          stop 'subroutine formain'
        end if
        write (23,99001) 'KX     =' , kz
        write (23,99001) 'NSG    =' , nsg
        write (23,99001) 'NNSG   =' , nsg*nsg
        write (23,99001) 'IBYTE  =' , ibyte
        write (23,99002) "DATTYP='" , dattyp
        if ( dattyp=='EH5OM' .and. ehso4 ) then
          write (23,'(a)') "      parameter(EHSO4 =.true. )"
        else
          write (23,'(a)') "      parameter(EHSO4 =.false.)"
        end if
        write (23,99003) "LSMTYP='" , lsmtyp
        write (23,99004) "AERTYP='" , aertyp
        write (23,'(a)') '      parameter(jxbb=mjx-1)'
        close (23)
      end if
99001 format ('      parameter(',a8,i6,')')
99002 format ('      parameter(',a8,a5,"')")
99003 format ('      parameter(',a8,a4,"')")
99004 format ('      parameter(',a8,a7,"')")
99005 format ('      parameter(',a18,")")
!
      end subroutine formain

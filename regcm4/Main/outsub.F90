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
 
      subroutine outsub
 
! ******    write bats fields to unit iutsub
 
      use mod_dynparam
      use mod_param1
      use mod_param2
      use mod_iunits
      use mod_date
      use mod_bats
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Local variables
!
#ifdef BAND
      real(4) , dimension(jxsg,iym2sg) :: v2b
#else
      real(4) , dimension(jxm2sg,iym2sg) :: v2b
#endif
      integer :: n

#ifdef MPP1
      if ( iotyp.eq.1 ) then
        do n = 1 , numsub
          nrcsub = nrcsub + 1
#ifdef BAND
          call reorder(fsub_io(1,1,1,n),v2b,jx,iym2,nsg)
#else
          call reorder(fsub_io(1,1,1,n),v2b,jxm2,iym2,nsg)
#endif
          write (iutsub,rec=nrcsub) v2b
        end do
      else
        write (iutsub) idatex
        do n = 1 , numsub
#ifdef BAND
          call reorder(fsub_io(1,1,1,n),v2b,jx,iym2,nsg)
#else
          call reorder(fsub_io(1,1,1,n),v2b,jxm2,iym2,nsg)
#endif
          write (iutsub) v2b
        end do
      end if
#else
      if ( iotyp.eq.1 ) then
        do n = 1 , numsub
          nrcsub = nrcsub + 1
#ifdef BAND
          call reorder(fsub(1,1,1,n),v2b,jx,iym2,nsg)
#else
          call reorder(fsub(1,1,1,n),v2b,jxm2,iym2,nsg)
#endif
          write (iutsub,rec=nrcsub) v2b
        end do
      else
        write (iutsub) idatex
        do n = 1 , numsub
#ifdef BAND
          call reorder(fsub(1,1,1,n),v2b,jx,iym2,nsg)
#else
          call reorder(fsub(1,1,1,n),v2b,jxm2,iym2,nsg)
#endif
          write (iutsub) v2b
        end do
      end if
#endif

      write (*,*) 'sub_BATS variables written at ' , idatex , xtime

      end subroutine outsub

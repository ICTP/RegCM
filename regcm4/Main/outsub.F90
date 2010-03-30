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
 
      subroutine outsub
 
! ******    write bats fields to unit iutsub
 
      use mod_regcm_param
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
      real(4) , dimension(jxm2sg,ixm2sg) :: v2b
      integer :: n

#ifdef MPP1
      if ( iotyp.eq.1 ) then
        do n = 1 , numsub
          nrcsub = nrcsub + 1
          call reorder(fsub_io(1,1,1,n),v2b,jxm2,ixm2,nsg)
          write (iutsub,rec=nrcsub) v2b
        end do
      else
        write (iutsub) idatex
        do n = 1 , numsub
          call reorder(fsub_io(1,1,1,n),v2b,jxm2,ixm2,nsg)
          write (iutsub) v2b
        end do
      end if
#else
      if ( iotyp.eq.1 ) then
        do n = 1 , numsub
          nrcsub = nrcsub + 1
          call reorder(fsub(1,1,1,n),v2b,jxm2,ixm2,nsg)
          write (iutsub,rec=nrcsub) v2b
        end do
      else
        write (iutsub) idatex
        do n = 1 , numsub
          call reorder(fsub(1,1,1,n),v2b,jxm2,ixm2,nsg)
          write (iutsub) v2b
        end do
      end if
#endif

      write (*,*) 'sub_BATS variables written at ' , idatex , xtime

      end subroutine outsub

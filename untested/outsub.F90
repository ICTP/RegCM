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
#ifdef MPP1
      integer :: n
      real(4) , dimension((mjx-2)*nsg,(ix-2)*nsg) :: v2b
#else
      real(4) , dimension((jx-2)*nsg,(ix-2)*nsg) :: v2b
#endif
!
!     ****** check if at desired output time for bats variables
      write (*,*) 'sub_BATS variables written at ' , idatex , xtime
#ifdef MPP1
      if ( iotyp.eq.1 ) then
        do n = 1 , numsub
          nrcsub = nrcsub + 1
          call reorder(fsub_io(1,1,1,n),v2b,mjx-2,ix-2,nsg)
          write (iutsub,rec=nrcsub) v2b
        end do
      else
        write (iutsub) idatex
        do n = 1 , numsub
          call reorder(fsub_io(1,1,1,n),v2b,mjx-2,ix-2,nsg)
          write (iutsub) v2b
        end do
      end if
#else
      if ( iotyp.eq.1 ) then
        nrcsub = nrcsub + 1
        call reorder(u10m_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(v10m_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(drag_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(tg_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(tlef_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(t2m_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(q2m_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(ssw_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(rsw_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(tpr_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(evpa_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(rnos_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(scv_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(sena_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(prcv_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
        nrcsub = nrcsub + 1
        call reorder(ps_s,v2b,jx-2,ix-2,nsg)
        write (iutsub,rec=nrcsub) v2b
      else if ( iotyp.eq.2 ) then
        write (iutsub) idatex
        call reorder(u10m_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(v10m_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(drag_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(tg_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(tlef_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(t2m_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(q2m_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(ssw_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(rsw_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(tpr_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(evpa_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(rnos_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(scv_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(sena_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(prcv_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
        call reorder(ps_s,v2b,jx-2,ix-2,nsg)
        write (iutsub) v2b
      else
      end if
#endif
 
      end subroutine outsub

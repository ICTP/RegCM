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

      subroutine outsub

c ******    write bats fields to unit iutsub

      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/date.cb"
#ifdef MPP1
      integer nrcbat,nrcsub
      COMMON /BATUNIT/ nrcbat,nrcsub
      integer numsub
      parameter( numsub=16 )
      real(kind=4)  fsub_io(NNSG,mjx-2,ix-2,numsub)
      COMMON /SUBOUTio/ fsub_io
#else
#     include "include/bats2.cb"
#endif
#     include "include/main.cb"

#ifdef MPP1
      real(kind=4) v2b((mjx-2)*NSG,(ix-2)*NSG)
      integer n
#else
      real(kind=4) v2b((jx-2)*NSG,(ix-2)*NSG)
#endif

c ****** check if at desired output time for bats variables
      write(*,*) 'sub_BATS variables written at ',idatex,xtime
#ifdef MPP1
      if (iotyp.eq.1) then
        do n=1,numsub
          nrcsub=nrcsub+1
          call reorder(fsub_io(1,1,1,n),v2b,mjx-2,ix-2,NSG)      
          write(iutsub,rec=nrcsub) v2b
        end do
      else
        write(iutsub) idatex
        do n=1,numsub
          call reorder(fsub_io(1,1,1,n),v2b,mjx-2,ix-2,NSG)      
          write(iutsub) v2b
        end do
      end if
#else
      if (iotyp.eq.1) then
        nrcsub=nrcsub+1
        call reorder(u10m_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(v10m_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(drag_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(tg_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(tlef_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(t2m_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(q2m_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(ssw_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(rsw_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(tpr_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(evpa_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(rnos_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(scv_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(sena_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(prcv_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
        nrcsub=nrcsub+1
        call reorder(ps_s,v2b,jx-2,ix-2,NSG)
        write(iutsub,rec=nrcsub) v2b
      else if (iotyp.eq.2) then
        write(iutsub) idatex
        call reorder(u10m_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(v10m_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(drag_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(tg_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(tlef_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(t2m_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(q2m_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(ssw_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(rsw_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(tpr_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(evpa_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(rnos_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(scv_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(sena_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(prcv_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
        call reorder(ps_s,v2b,jx-2,ix-2,NSG)
        write(iutsub) v2b
      end if
#endif

      return
      end subroutine outsub

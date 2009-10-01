      subroutine outsub

c ******    write bats fields to unit iutsub

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param2.cb'
      include 'Commons/param3.cb'
      include 'Commons/iunits.cb'
      include 'Commons/date.cb'
#ifdef MPP1
      integer nrcbat,nrcsub
      COMMON /BATUNIT/ nrcbat,nrcsub
      integer numsub
      parameter( numsub=16 )
      real*4  fsub_io(NNSG,mjx-2,ix-2,numsub)
      COMMON /SUBOUTio/ fsub_io
#else
      include 'Commons/bats.cb2'
#endif
#ifdef MPP1
      include 'Common2/main.cb'
#else
      include 'Commons/main.cb'
#endif

#ifdef MPP1
      real*4 v2b((mjx-2)*NSG,(ix-2)*NSG)
      integer n
#else
      real*4 v2b((jx-2)*NSG,(ix-2)*NSG)
#endif

c ****** check if at desired output time for bats variables
      write(*,*) 'sub_BATS variables written at ',idatex,xtime
#ifdef MPP1
      if (iotyp.eq.1) then
        do n=1,numsub
          nrcsub=nrcsub+1
          call reorder(fsub_io(1,1,1,n),v2b,mjx-2,ix-2,NSG)      
          write(iutsub,rec=nrcsub) v2b
        enddo
      else
        write(iutsub) idatex
        do n=1,numsub
          call reorder(fsub_io(1,1,1,n),v2b,mjx-2,ix-2,NSG)      
          write(iutsub) v2b
        enddo
      endif
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
      elseif (iotyp.eq.2) then
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
      end

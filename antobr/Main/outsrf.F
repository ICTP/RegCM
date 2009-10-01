      subroutine outsrf

c ******    write bats fields to unit iutbat

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param2.cb'
      include 'Commons/iunits.cb'
      include 'Commons/date.cb'

      integer i,j,n

      integer numbat
      parameter( numbat=21+6 )
      integer nrcbat,nrcsub
      COMMON /BATUNIT/ nrcbat,nrcsub
#ifdef MPP1
      real*4  fbat_io(mjx-2,ix-2,numbat)
      COMMON /BATOUTio/fbat_io
      real*4  fbat(jxp,ix-2,numbat)
#else
      real*4  fbat(jx-2,ix-2,numbat)
#endif
      COMMON /BATOUT/fbat

c ****** check if at desired output time for bats variables
      write(*,*) 'BATS variables written at ',idatex,xtime
      if (iotyp.eq.2) write(iutbat) idatex
      do n=1,numbat
        if (iotyp.eq.1) then
          nrcbat=nrcbat+1
#ifdef MPP1
       write(iutbat,rec=nrcbat)((fbat_io(j,i,n),j=1,mjx-2),i=1,ix-2)
#else
          write(iutbat,rec=nrcbat)((fbat(j,i,n),j=1,jx-2),i=1,ix-2)
#endif
        elseif (iotyp.eq.2) then
#ifdef MPP1
          write(iutbat) ((fbat_io(j,i,n),j=1,mjx-2),i=1,ix-2)
#else
          write(iutbat) ((fbat(j,i,n),j=1,jx-2),i=1,ix-2)
#endif
        end if
      enddo

      return
      end

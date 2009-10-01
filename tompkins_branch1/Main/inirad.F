c-----------------------------------------------------------------------
      subroutine inirad

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param3.cb'
#ifdef MPP1
      include 'Common2/rad.cb'
      include 'mpiregcm.h'
#else
      include 'Commons/rad.cb'
#endif
      integer i,j,k
c
c  compute ozone mixing ratio distribution
c
      if (jyear.eq.jyear0.and.ktau.eq.0) then
       do 10 k=1,kx
#ifdef MPP1
       do 10 j=1,jendl
#else
       do 10 j=1,jx-1
#endif
       do 10 i=1,ix-1
        heatrt(i,k,j)=0.
        o3prof(i,k,j)=0.
10     continue
#ifdef MPP1
       do 11 j=1,jendl
#else
       do 11 j=1,jx-1
#endif
       do 11 i=1,ix-1
        o3prof(i,kxp1,j)=0.
11     continue
       call o3data
#ifdef MPP1
       if (myid .eq. 0) then
#endif
       print *,'ozone profiles'
       do 34 k=1,kx+1
        write(6,1212)o3prof(3,k,3)
34     continue
1212   format (1x,7e12.4)
#ifdef MPP1
        endif
#endif
      endif

      return
      end

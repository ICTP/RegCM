ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tstep(extime,dtinc,deltmx)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c                                                                     c
c     this subroutine makes a refined start to the model, i.e. the    c
c     model is started with a series of small time steps.             c
c                                                                     c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8  extime,dtinc,deltmx
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      integer idtmax
      real*8  tscale,deltmn
c---------------------------------------------------------------------
c
      if (jyear.eq.jyear0.and.ktau.eq.0) then
       idtmax=1
       tscale=5.*dt
       deltmn=0.1*dt
       deltmx=dt
      endif
      if (idtmax .eq. 1) then
       dt=deltmx*(1.-dexp(-extime/tscale))+deltmn
       dtinc=dt
       dtmin=dt/60.
       idtmax=2
      endif
      dt2=2.*dt
      if (.not.(jyear.eq.jyear0.and.ktau.eq.0)) dt=dt2
c
      return
      end

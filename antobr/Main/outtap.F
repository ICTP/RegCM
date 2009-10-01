ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outtap
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine writes the model output to tape or disk for use c
c     in dataflow analyses.                                           c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
#ifdef MPP1
      include 'Common2/cvaria.cb'
      include 'Common2/bats.cb2'
      include 'Common2/io.cb'
      real*8  out1(mjx,ix)
      real*4  fout(mjx-2,ix-2)
#else
      include 'Commons/cvaria.cb'
      include 'Commons/bats.cb2'
      include 'Commons/main.cb'
      real*8  out1(jx,ix)
      real*4  fout(jx-2,ix-2)
#endif
      include 'Commons/date.cb'
      integer nrcout
      COMMON /GENREC/ nrcout
      real*8  mmpd
c
      integer i,j,k,n,nn
      real*8  xh
c
c---------------------------------------------------------------------
c-----output RegCM-domain variables:
c
      xh=xtime/60.
      if (iotyp.eq.2) write (iutdat) idatex
c
c ******  write one time of data to mm4 output file.

      do k=kx,1,-1
         do i=1,ix
#ifdef MPP1
           do j=1,mjx
             out1(j,i) = ua_io(i,k,j)
#else
           do j=1,jx
             out1(j,i) = ua(i,k,j)
#endif
           enddo
         enddo
         do i=1,ix-2
#ifdef MPP1
           do j=1,mjx-2
             fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)
     &                        +out1(j+1,i+2)+out1(j+2,i+2))
     &                        /psa_io(i+1,j+1)
#else
           do j=1,jx-2
             fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)
     &                        +out1(j+1,i+2)+out1(j+2,i+2))
     &                        /psa(i+1,j+1)
#endif
           enddo
         enddo 
         if (iotyp.eq.1) then
           nrcout=nrcout+1
           write(iutdat,rec=nrcout) fout
         elseif (iotyp.eq.2) then
           write(iutdat) fout
         end if
      enddo
      do k=kx,1,-1
         do i=1,ix
#ifdef MPP1
           do j=1,mjx
             out1(j,i) = va_io(i,k,j)
#else
           do j=1,jx
             out1(j,i) = va(i,k,j)
#endif
           enddo
         enddo
         do i=1,ix-2
#ifdef MPP1
           do j=1,mjx-2
             fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)
     &                        +out1(j+1,i+2)+out1(j+2,i+2))
     &                       /psa_io(i+1,j+1)
#else
           do j=1,jx-2
             fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)
     &                        +out1(j+1,i+2)+out1(j+2,i+2))
     &                       /psa(i+1,j+1)
#endif
           enddo
         enddo
         if (iotyp.eq.1) then
           nrcout=nrcout+1
           write(iutdat,rec=nrcout) fout
         elseif (iotyp.eq.2) then
           write(iutdat) fout
         end if
      enddo
      do k=kx,1,-1
        do i=1,ix-2 
#ifdef MPP1
          do j=1,mjx-2
            fout(j,i) = omega_io(i+1,k,j+1)
#else
          do j=1,jx-2
            fout(j,i) = omega(i+1,k,j+1)
#endif
          enddo
        enddo
        if (iotyp.eq.1) then
          nrcout=nrcout+1
          write(iutdat,rec=nrcout) fout
        elseif (iotyp.eq.2) then
          write(iutdat) fout
        end if
      enddo
      do k=kx,1,-1
        do i=1,ix-2 
#ifdef MPP1
          do j=1,mjx-2
            fout(j,i) = ta_io(i+1,k,j+1)/psa_io(i+1,j+1)
#else
          do j=1,jx-2
            fout(j,i) = ta(i+1,k,j+1)/psa(i+1,j+1)
#endif
          enddo
        enddo
        if (iotyp.eq.1) then
          nrcout=nrcout+1
          write(iutdat,rec=nrcout) fout
        elseif (iotyp.eq.2) then
          write(iutdat) fout
        end if
      enddo
      do k=kx,1,-1
        do i=1,ix-2 
#ifdef MPP1
          do j=1,mjx-2
            fout(j,i) = qva_io(i+1,k,j+1)/psa_io(i+1,j+1)
#else
          do j=1,jx-2
            fout(j,i) = qva(i+1,k,j+1)/psa(i+1,j+1)
#endif
          enddo
        enddo
        if (iotyp.eq.1) then
          nrcout=nrcout+1
          write(iutdat,rec=nrcout) fout
        elseif (iotyp.eq.2) then
          write(iutdat) fout
        end if
      enddo
      do k=kx,1,-1
        do i=1,ix-2 
#ifdef MPP1
          do j=1,mjx-2
            fout(j,i) = qca_io(i+1,k,j+1)/psa_io(i+1,j+1)
#else
          do j=1,jx-2
            fout(j,i) = qca(i+1,k,j+1)/psa(i+1,j+1)
#endif
          enddo
        enddo
        if (iotyp.eq.1) then
          nrcout=nrcout+1
          write(iutdat,rec=nrcout) fout
        elseif (iotyp.eq.2) then
          write(iutdat) fout
        end if
      enddo
      do i=1,ix-2 
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = (psa_io(i+1,j+1)+ptop)*10.
#else
        do j=1,jx-2
          fout(j,i) = (psa(i+1,j+1)+ptop)*10.
#endif
        enddo
      enddo
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      elseif (iotyp.eq.2) then
        write(iutdat) fout
      end if
      mmpd = 24./tapfrq
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = (rainc_io(i+1,j+1)+rainnc_io(i+1,j+1))*mmpd
#else
        do j=1,jx-2
          fout(j,i) = (rainc(i+1,j+1)+rainnc(i+1,j+1))*mmpd
#endif
        enddo 
      enddo 
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      elseif (iotyp.eq.2) then
        write(iutdat) fout
      end if

c ****** write out the following surface fields:
c                  1.  temp of lower soil layer (17)
c                  2.  total soil water in mm h2o (13)
c                  3.  accum infiltration (30)
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = tgb2d_io(1,i+1,j+1)
          do n=2,NNSG
            fout(j,i) = fout(j,i) + tgb2d_io(n,i+1,j+1)
          enddo
          fout(j,i) = fout(j,i)/float(NNSG)
#else
        do j=1,jx-2
          fout(j,i) = tgb2d(1,i+1,j+1)
          do n=2,NNSG
            fout(j,i) = fout(j,i) + tgb2d(n,i+1,j+1)
          enddo
          fout(j,i) = fout(j,i)/float(NNSG)
#endif
        enddo 
      enddo 
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      elseif (iotyp.eq.2) then
        write(iutdat) fout
      end if

      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = 0.0
          nn=0
          do n=1,NNSG
            if (ocld2d_io(n,i,j).ge.0.5) then
              fout(j,i) = fout(j,i)+swt2d_io(n,i+1,j+1)
              nn=nn+1
            endif
          enddo
#else
        do j=1,jx-2
          fout(j,i) = 0.0
          nn=0
          do n=1,NNSG
            if (ocld2d(n,i,j).ge.0.5) then
              fout(j,i) = fout(j,i)+swt2d(n,i+1,j+1)
              nn=nn+1
            endif
          enddo
#endif
          if(nn.ge.max0(NNSG/2,1)) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.e34
          endif
        enddo 
      enddo 
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      elseif (iotyp.eq.2) then
        write(iutdat) fout
      end if

      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = 0.0
          nn=0
          do n=1,NNSG
            if (ocld2d_io(n,i,j).ge.0.5) then
              fout(j,i) = fout(j,i)+rno2d_io(n,i+1,j+1)
              nn=nn+1
            endif
          enddo
#else
        do j=1,jx-2
          fout(j,i) = 0.0
          nn=0
          do n=1,NNSG
            if (ocld2d(n,i,j).ge.0.5) then
              fout(j,i) = fout(j,i)+rno2d(n,i+1,j+1)
              nn=nn+1
            endif
          enddo
#endif
          if(nn.ge.max0(NNSG/2,1)) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.e34
          endif
        enddo 
      enddo 
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      elseif (iotyp.eq.2) then
        write(iutdat) fout
      end if

c     changes for accum infiltration

      do i=1,ix-1
#ifdef MPP1
        do j=1,mjx-1
          do n=1,NNSG
            rno2d_io(n,i,j) = 0.
          enddo
          rainc_io(i,j) = 0.
          rainnc_io(i,j) = 0.
#else
        do j=1,jx-1
          do n=1,NNSG
            rno2d(n,i,j) = 0.
          enddo
          rainc(i,j) = 0.
          rainnc(i,j) = 0.
#endif
        enddo
      enddo

      print *,'OUT-history written date = ',ldatez+xtime/1440.

      return
      end

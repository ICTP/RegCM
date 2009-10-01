      subroutine slice(j)

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
#ifdef MPP1
      include 'Common2/main.cb'
      include 'Common2/cvaria.cb'
      include 'Common2/pbldim.cb'
      include 'Common2/pmoist.cb'
      include 'Common2/slice.cb'
      include 'Common2/bats.cb2'
#else
      include 'Commons/main.cb'
      include 'Commons/cvaria.cb'
      include 'Commons/pbldim.cb'
      include 'Commons/pmoist.cb'
      include 'Commons/slice.cb'
      include 'Commons/bats.cb2'
#endif
      include 'Commons/bats.cb'

      integer j,n,np
      real*8  amxtem,sfac

c ********* For albedov
c **********            This is taken from subroutine interf so that
c **********            radiation can be called in tend (not vecbats).
      npts=nbmax
      do np=np1,npts
      do n=1,NNSG
        ldoc1d(n,np)=ocld2d(n,np,j)
        sice1d(n,np)=sice2d(n,np,j)
        tgb1d(n,np)=tgb2d(n,np,j)
        ssw1d(n,np)=ssw2d(n,np,j)
        lveg(n,np) = nint(veg2d1(n,np,j))
        amxtem=dmax1(298.-tgb1d(n,np),0.d0)
        sfac=1.-dmax1(0.d0,1.-0.0016*amxtem**2)
        if(lveg(n,np).eq.0)then
         veg1d(n,np)=0.
        else
         veg1d(n,np)=vegc(lveg(n,np))-seasf(lveg(n,np))*sfac
        endif
        ts1d(n,np)=thx3d(np,kx,j)-6.5e-3/g*(ht1(n,np,j)-ht(np,j))
        scv1d(n,np)=scv2d(n,np,j)
        sag1d(n,np)=sag2d(n,np,j)
      enddo
      end do

      return
      end

      subroutine cldfrac(j)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     This subroutine computes the fractional cloud coverage and      c
c     liquid water content (in cloud value).  Both are use in         c
c     radiation.                                                      c
c                                                                     c
c     The fractional coverage of large scale clouds is a function of  c
c     relative humidity, using the relationship of sundqvist et       c
c     al., 1989.  The relative humidity at which clouds begin to      c
c     form is lower over land than ocean, due to the greater number   c
c     of cloud condensation nucleii.                                  c
c                                                                     c
c     The fracional coverage of convective clouds is passed in from   c
c     the convection scheme.                                          c
c                                                                     c
c     The large-scale and convective clouds are combined as follows:  c
c     1) If the convective cloud fraction > large scale fraction, the c
c     convective fraction and water content are used (this occurs     c
c     infrequently).                                                  c
c     2) Otherwise, the cloud fraction equals the large-scale         c
c     fraction AND the water content is a weighted average of both    c
c     types.                                                          c
c                                                                     c
c     Note: the incloud water content (g/m3) is passed to radiation   c
c                                                                     c
c     See Pal et al (200) for more info.                              c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer j
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
      include 'Common2/pmoist.cb'
      include 'Common2/rad.cb'
      include 'Common2/slice.cb'
      include 'Common2/bats.cb2'  ! so you have the land mask
#else
      include 'Commons/main.cb'
      include 'Commons/cvaria.cb'
      include 'Commons/pmoist.cb'
      include 'Commons/rad.cb'
      include 'Commons/slice.cb'
      include 'Commons/bats.cb2'  ! so you have the land mask
#endif
      integer i,k
      real*8  exlwc
      real*8  rh0adj  ! Adjusted relative humidity threshold


c--------------------------------------------------------------------
c    1.  Determine large-scale cloud fraction
c--------------------------------------------------------------------
      do k=1,kx
      do i=2,ilxm
        if (tb3d(i,k,j).gt.tc0) then
          rh0adj = rh0(i,j)
        else   ! high cloud (less subgrid variability)
          rh0adj = rhmax-(rhmax-rh0(i,j))/(1.0+0.15*(tc0-tb3d(i,k,j)))
        end if
        if (rhb3d(i,k,j).ge.rhmax) then        ! full cloud cover
          fcc(i,k,j) = 1.
        else if (rhb3d(i,k,j).le.rh0adj) then  ! no cloud cover
          fcc(i,k,j) = 0.
        else                                  ! partial cloud cover
          fcc(i,k,j) = 1.-sqrt(1.-(rhb3d(i,k,j)-rh0adj)/(rhmax-rh0adj))
          fcc(i,k,j) = dmin1(dmax1(fcc(i,k,j),0.01D0),0.99D0)
        endif                        !  rh0 threshold
      end do
      end do

c--------------------------------------------------------------------
c    2.  Combine large-scale and convective fraction and liquid water
c        to be passed into radiation.
c--------------------------------------------------------------------
      do k=1,kx
      do i=2,ilxm
        exlwc = qcb3d(i,k,j)*rhob3d(i,k,j)*1.e03  ! Cloud Water Volume
c units:  = kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
        cldlwc(i,k) = (cldfra(i,k)*cldlwc(i,k)*1.e3 + fcc(i,k,j)*exlwc)
     &              / dmax1(cldfra(i,k)+fcc(i,k,j),0.01d0)
        cldfra(i,k) = dmin1(dmax1(cldfra(i,k),fcc(i,k,j)),fcmax)
      end do
      end do

      return
      end

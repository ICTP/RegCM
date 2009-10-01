
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine trcplk(tint, tlayr, tplnke, emplnk, abplnk1, abplnk2)
c----------------------------------------------------------------------
c   Calculate Planck factors for absorptivity and emissivity of
c   CH4, N2O, CFC11 and CFC12
c
c-----------------------------------------------------------------------
c
c $Id: trcplk.F,v 1.2 1995/02/17 21:29:00 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
c
      include 'Commons/parrad.cb'
C
C------------------------------Commons----------------------------------
      include 'Commons/crdcon.cb'
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real*8  tint(plond,plevp),  ! interface temperatures
     $     tlayr(plond,plevp), ! k-1 level temperatures
     $     tplnke(plond)       ! Top Layer temperature
c
c output arguments
c
      real*8  emplnk(14,plond),        ! emissivity Planck factor
     $     abplnk1(14,plond,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plond,plevp)  ! nearest layer factor
c
c local workspace
c
      integer wvl                   ! wavelength index
      integer i,k
      real*8  f1(14),                  ! Planck function factor
     $     f2(14),                  !        "
     $     f3(14)                   !        "
c
      data f1 /5.85713d8,7.94950d8,1.47009d9,1.40031d9,1.34853d8,
     $         1.05158d9,3.35370d8,3.99601d8,5.35994d8,8.42955d8,
     $         4.63682d8,5.18944d8,8.83202d8,1.03279d9/
      data f2 /2.02493d11,3.04286d11,6.90698d11,6.47333d11,
     $         2.85744d10,4.41862d11,9.62780d10,1.21618d11,
     $         1.79905d11,3.29029d11,1.48294d11,1.72315d11,
     $         3.50140d11,4.31364d11/
      data f3 /1383.d0,1531.d0,1879.d0,1849.d0,848.d0,1681.d0,
     $         1148.d0,1217.d0,1343.d0,1561.d0,1279.d0,1328.d0,
     $         1586.d0,1671.d0/
c
c Calculate emissivity Planck factor
c
      do wvl = 1,14
         do i = 1,plon
            emplnk(wvl,i) = f1(wvl)/
     $                   (tplnke(i)**4.0*(dexp(f3(wvl)/tplnke(i))-1.0))
         end do
      end do
c
c Calculate absorptivity Planck factor for tint and tlayr temperatures
c
      do wvl = 1,14
         do k = 1, plevp
            do i = 1, plon
c non-nearlest layer function
               abplnk1(wvl,i,k) = (f2(wvl)*dexp(f3(wvl)/tint(i,k)))
     $              /(tint(i,k)**5.0*(dexp(f3(wvl)/tint(i,k))-1.0)**2.0)
c nearest layer function
               abplnk2(wvl,i,k) = (f2(wvl)*dexp(f3(wvl)/tlayr(i,k)))
     $            /(tlayr(i,k)**5.0*(dexp(f3(wvl)/tlayr(i,k))-1.0)**2.0)
            end do
         end do
      end do
      return
      end

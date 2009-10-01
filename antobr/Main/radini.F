
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine radini(gravx   ,cpairx  ,epsilox ,stebolx )
C-----------------------------------------------------------------------
C
C Initialize various constants for radiation scheme; note that
C the radiation scheme uses cgs units.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radini.F,v 1.2 1995/02/17 21:28:44 jhack Exp $
c $Author: jhack $
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
      include 'Commons/comozp.cb'
      include 'Commons/crdcae.cb'
      include 'Commons/crdcon.cb'
C_IPCC
      include 'Commons/tracer0.cb'
      include 'Commons/tracer.cb'
      include 'Commons/date.cb'
CIPCC_
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real*8  gravx,     ! Acceleration of gravity (MKS)
     $     cpairx,    ! Specific heat of dry air (MKS)
     $     epsilox,   ! Ratio of mol. wght of H2O to dry air
     $     stebolx    ! Stefan-Boltzmann's constant (MKS)
C
C---------------------------Local variables-----------------------------
C
      integer iband   ! H2O band index
      real*8  v0,        ! Volume of a gas at stp (m**3/kmol)
     $     p0,        ! Standard pressure (pascals)
     $     amd,       ! Effective molecular weight of dry air (kg/kmol)
     $     goz        ! Acceleration of gravity (m/s**2)
C
C-----------------------------------------------------------------------
C
C Set general radiation consts; convert to cgs units where appropriate:
C
      gravit  =  100.*gravx
      rga     =  1./gravit
      cpair   =  1.e4*cpairx
      epsilo  =  epsilox
      sslp    =  1.013250e6
      stebol  =  1.e3*stebolx
      rgsslp  =  0.5/(gravit*sslp)
cIPCC
c1991-1995
c     co2vmr  =  3.55e-4
c     ch40 = 0.55241 * 1.714e-6
c     n2o0 = 1.51913 * 0.311e-6
c1961-1965
c     co2vmr  =  3.10e-4
c     ch40 = 0.55241 * 1.414e-6
c     n2o0 = 1.51913 * 0.287e-6

c     cfc110 = 4.69548 * 0.280e-9
c     cfc120 = 4.14307 * 0.503e-9
c     co2mmr = 1.51913 * co2vmr
      if (lyear.ge.1750 .and. lyear.le.2100) then
        co2vmr = cgas(2,lyear)*1.e-6
        co2mmr = co2vmr*44.0/28.9644
        ch40   = cgas(3,lyear)*1.e-9*0.55241
        n2o0   = cgas(4,lyear)*1.e-9*1.51913
        cfc110 = cgas(5,lyear)*1.e-12*4.69548
        cfc120 = cgas(6,lyear)*1.e-12*4.14307
      else
        print*,'CONCENTRATION VALUES OUTSIDE OF DATE RANGE (1750-2100)'
        print*,'  Simulation date:  ',lyear
        stop 'SUBROUTINE RADINI'
      end if
c     print*,'IN RADINI (TOP)'
c     print*,'  co2vmr= ',co2vmr
c     print*,'  co2mmr= ',co2mmr
c     print*,'  ch40  = ',ch40
c     print*,'  n2o0  = ',n2o0
c     print*,'  cfc110= ',cfc110
c     print*,'  cfc120= ',cfc120
cIPCC_
cIPCC_

      dpfo3   =  2.5e-3
      dpfco2  =  5.0e-3
      dayspy  =  365.
      pie     =  4.*datan(1.d0)
C
C Coefficients for h2o emissivity and absorptivity.
C
      do iband=1,4
         c1(iband) = coefe(3,iband)/coefe(2,iband)
         c2(iband) = coefb(3,iband)/coefb(2,iband)
         c3(iband) = coefb(4,iband)/coefb(3,iband)
         c4(iband) = coefd(3,iband)/coefd(2,iband)
         c5(iband) = coefd(4,iband)/coefd(3,iband)
         c6(iband) = coefa(3,iband)/coefa(2,iband)
         c7(iband) = coefc(3,iband)/coefc(2,iband)
      end do
      c8   = coeff(3,1)/coeff(2,1)
      c9   = coeff(3,2)/coeff(2,2)
      c10  = coeff(4,1)/coeff(3,1)
      c11  = coeff(4,2)/coeff(3,2)
      c12  = coeff(5,1)/coeff(4,1)
      c13  = coeff(5,2)/coeff(4,2)
      c14  = coeff(6,1)/coeff(5,1)
      c15  = coeff(6,2)/coeff(5,2)
      c16  = coefj(3,1)/coefj(2,1)
      c17  = coefk(3,1)/coefk(2,1)
      c18  = coefi(3,1)/coefi(2,1)
      c19  = coefi(3,2)/coefi(2,2)
      c20  = coefi(4,1)/coefi(3,1)
      c21  = coefi(4,2)/coefi(3,2)
      c22  = coefi(5,1)/coefi(4,1)
      c23  = coefi(5,2)/coefi(4,2)
      c24  = coefi(6,1)/coefi(5,1)
      c25  = coefi(6,2)/coefi(5,2)
      c26  = coefj(3,2)/coefj(2,2)
      c27  = coefk(3,2)/coefk(2,2)
      c28  = .5
      c29  = .002053
      c30  = .1
      c31  = 3.0e-5
      cfa1 = .61
C
C Initialize further longwave constants referring to far wing
C correction; R&D refers to:
C
C            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
C            Emissivity and Absorptivity Formulation for Water Vapor
C            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
C
      fwcoef = .1           ! See eq(33) R&D
      fwc1   = .30          ! See eq(33) R&D
      fwc2   = 4.5          ! See eq(33) and eq(34) in R&D
      fc1    = 2.6          ! See eq(34) R&D
C
C Initialize ozone data.
C
      v0  = 22.4136         ! Volume of a gas at stp (m**3/kmol)
      p0  = 0.1*sslp        ! Standard pressure (pascals)
      amd = 28.9644         ! Molecular weight of dry air (kg/kmol)
      goz = gravx           ! Acceleration of gravity (m/s**2)
C
C Constants for ozone path integrals (multiplication by 100 for unit
C conversion to cgs from mks):
C
      cplos = v0/(amd*goz)       *100.0
      cplol = v0/(amd*goz*p0)*0.5*100.0
C
      return
      end

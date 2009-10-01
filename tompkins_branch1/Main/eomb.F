
*****************************FUNCTION EOMB******************************
*          COMPUTES AIR VAPOR PRESSURE AS A FUNCTION OF TEMP (in mb)   *
************************************************************************
      function eomb(x)
      implicit none
      real*8  eomb,x
      real*8  tr1
c
      tr1=1.0-(373.15/(x+273.15))
      eomb=1013.25*dexp(13.3185*tr1-1.976*tr1**2-0.6445*tr1**3
     &   -0.1299*tr1**4)
      return
      end

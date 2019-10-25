!
! The following formula gives the water vapour saturation pressure to
! sufficient accuracy between -100 C and 373 C
!
! W. Wagner and A. Pruss
!
! The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary
! Water Substance for General and Scientific Use,
! Journal of Physical and Chemical Reference Data, June 2002 ,Volume 31,
! Issue 2, pp.  387535
!
!DIR$ ATTRIBUTES FORCEINLINE :: pfesat
  pure elemental real(rkx) function pfesat(t) result(es) ! Pressure in Pa
   implicit none
   real(rkx) , intent(in) :: t     ! Temperature (K)
   real(rk8) :: xt , rt , theta , lnp , ps
   real(rk8) , parameter :: tc = 647.096_rk8
   real(rk8) , parameter :: tn = 273.16_rk8
   real(rk8) , parameter :: pc = 220640.0_rk8
   real(rk8) , parameter :: pn = 6.11657_rk8
   real(rk8) , parameter :: c1 = -7.85951783_rk8
   real(rk8) , parameter :: c2 =  1.84408259_rk8
   real(rk8) , parameter :: c3 = -11.7866497_rk8
   real(rk8) , parameter :: c4 =  22.6807411_rk8
   real(rk8) , parameter :: c5 = -15.9618719_rk8
   real(rk8) , parameter :: c6 =  1.80122502_rk8
   real(rk8) , parameter :: e1 = 1.0_rk8
   real(rk8) , parameter :: e2 = 1.5_rk8
   real(rk8) , parameter :: e3 = 3.0_rk8
   real(rk8) , parameter :: e4 = 3.5_rk8
   real(rk8) , parameter :: e5 = 4.0_rk8
   real(rk8) , parameter :: e6 = 7.5_rk8
   real(rk8) , parameter :: a0 = -13.928169_rk8
   real(rk8) , parameter :: a1 =  34.707823_rk8
   real(rk8) , parameter :: p1 = -1.50_rk8
   real(rk8) , parameter :: p2 = -1.25_rk8

   xt = real(t,rk8)
   if ( xt > 273.15_rk8 ) then
     rt = xt/tc
     theta = 1.0_rk8 - rt
     lnp = (c1*theta**e1 + c2*theta**e2 + &
            c3*theta**e3 + c4*theta**e4 + &
            c5*theta**e5 + c6*theta**e6)/rt
     ps = pc * exp(lnp)
   else
     theta = xt/tn
     lnp = a0*(1.0_rk8-theta**p1)+a1*(1.0_rk8-theta**p2)
     ps = pn * exp(lnp)
   end if
   es = real(100.0_rk8*ps,rkx)
 end function pfesat

 ! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

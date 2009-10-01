
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine aermix(pint  ,h2ommr  ,rh ,j)
C-----------------------------------------------------------------------
C Set global mean tropospheric aerosol
C
C Specify aerosol mixing ratio and compute relative humidity for later
C adjustment of aerosol optical properties. Aerosol mass mixing ratio
C is specified so that the column visible aerosol optical depth is a
C specified global number (tauvis). This means that the actual mixing
C ratio depends on pressure thickness of the lowest three atmospheric
C layers near the surface.
C
C Optical properties and relative humidity parameterization are from:
C
C J.T. Kiehl and B.P. Briegleb  "The Relative Roles of Sulfate Aerosols
C and Greenhouse Gases in Climate Forcing"  Science  260  pp311-314
C 16 April 1993
C
C Visible (vis) here means 0.5-0.7 micro-meters
C Forward scattering fraction is taken as asymmetry parameter squared
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb  March 1995
C Standarized:       L. Buja,     Feb 1996
C Reviewed:          B. Briegleb, Mar 1996
C
C-----------------------------------------------------------------------
c
c $Id: aermix.F,v 1.1 1998/04/01 07:20:42 ccm Exp $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
      integer j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/parrad.cb'

C------------------------------Commons----------------------------------
      include 'Commons/param2.cb'
      include 'Commons/crdcon.cb'
#ifdef MPP1
      include 'Common2/slice.cb'
      include 'Common2/trachem.cb'
      include 'Common2/main.cb'
      include 'Common2/mainchem.cb'
      real*8  aermm
      common /aerm/ aermm(ix-1,kx,jxp)
#else
      include 'Commons/slice.cb'
      include 'Commons/trachem.cb'
      include 'Commons/main.cb'
      include 'Commons/mainchem.cb'
      real*8  aermm
      common /aerm/ aermm(ix-1,kx,jx-1)
#endif
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real*8  pint(plond,plevrp) ! Radiation level interface pressures
C                                ! (dynes/cm2)
      real*8  h2ommr(plond,plev) ! Radiation level specific humidity (g/g)
C
C Output arguments
C
     
      real*8  rh(plon,plevr)     ! Radiation level relative humidity (fraction)

Cchem2
      
      real*8  aermmr        ! Radiation level aerosol mass mixing ratio
      real*8  aermmb        ! Background aerosol mmr

      common/aermr2/aermmr(plond,plevr,ntr),
     &             aermmb(plon,plevr)  
      

Cchem2_
c
C
C---------------------------Local variables-----------------------------
C
      integer i,     ! longitude index
     $        k,     ! level index
     $        mxaerl, ! max nmbr aerosol levels counting up from surface
     $        itr
C
      real*8  tauvis,  ! visible optical depth
     $        kaervs,  ! visible extinction coefficiant of aerosol (m2/g)
     $        omgvis,  ! visible omega0
     $        gvis,    ! visible forward scattering asymmetry parameter
     $        rhcnst   ! constant relative humidity factor
C
C Relative humidity factor
C
      real*8  rhfac,              ! multiplication factor for kaer
     $        a0,                 ! constant in rh mult factor
     $        a1,                 ! constant in rh mult factor
     $        a2,                 ! constant in rh mult factor
     $        a3                  ! constant in rh mult factor
c
      data    a0 / -9.2906106183    /
      data    a1 /  0.52570211505   /
      data    a2 / -0.0089285760691 /
      data    a3 /  5.0877212432e-05/
C
      data    kaervs / 5.3012 /
      data    omgvis / 0.999999 /
      data    gvis   / 0.694889 /
      data    rhcnst / .80 /
C
      data    rhfac / 1.6718 /  ! EES added for efficiency
C
C--------------------------------------------------------------------------
      mxaerl = 4
C
Cfil  tauvis = .01
      tauvis = .04

C
C Set relative humidity and factor; then aerosol amount:
C
      do k=1,plevr
      do i=1,plon
Cadded July 13, 2000:
        rh(i,k)=dmin1(rhb3d(i,k,j),.99d0) ! needed for aerosols in radiation
CEES: do not change to 1.00:  wscoef(3,10) in radcsw = .9924 and is
c     divided by RH.  rh is limited to .99 to avoid dividing by zero
Cadded_
C
c define background aerosol
C Find constant aerosol mass mixing ratio for specified levels
C in the column, converting units where appropriate
C for thr momrnt no more used
        if( k .ge. plevrp-mxaerl ) then
          aermmb(i,k) = gravit*tauvis /
     $          (1.e4*kaervs*rhfac*(1.-omgvis*gvis*gvis)*
     $          (pint(i,plevrp)-pint(i,plevrp-mxaerl)))
        else
          aermmb(i,k) = 0.0
        endif

c       if(ichem .eq. 1 .and. idirect .eq. 1) then 
        if(ichem .eq. 1 ) then 
          do itr=1,ntr
c           aermmr(i,k,itr)= dmax1( chia(i,k,j,itr)/psa(i,j) 
c    &                             ,aermmb(i,k) )
            aermmr(i,k,itr)= chia(i,k,j,itr)/psa(i,j) 
          end do
        else
          if(EHSO4) then
             do itr=1,ntr
c              aermmr(i,k,itr) = aermmb(i,k) + aermm(i,k,j)
!Dec.11        aermmr(i,k,itr) = aermm(i,k,j)
               aermmr(i,k,itr) = 0.0d0
!Dec.11_
             end do
          else
             do itr=1,ntr
c              aermmr(i,k,itr)= 0.0d0
               aermmr(i,k,itr) = aermmb(i,k)
             end do
          endif
        endif
      enddo
      enddo
C       
      return
      end

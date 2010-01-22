C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine aermix(pint  ,rh ,j)
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

#     include "include/parrad.cb"
#     include "include/param2.cb"
#     include "include/crdcon.cb"
#     include "include/slice.cb"
#     include "include/trachem.cb"
#     include "include/main.cb"
#     include "include/mainchem.cb"

#ifdef MPP1
      real(kind=8)  aermm
      common /aerm/ aermm(ix-1,kx,jxp)
#else
      real(kind=8)  aermm
      common /aerm/ aermm(ix-1,kx,jx-1)
#endif

C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(kind=8)  pint(plond,plevrp) ! Radiation level interface
c                                pressures (dynes/cm2)
C
C Output arguments
C
     
      real(kind=8)  rh(plon,plevr)  ! Radiation level relative 
c                                     humidity (fraction)

Cchem2
      
      real(kind=8)  aermmr        ! Radiation level aerosol mass
c                                   mixing ratio
      real(kind=8)  aermmb        ! Background aerosol mmr

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
      real(kind=8)  tauvis,  ! visible optical depth
     $        kaervs,  ! visible extinction coefficiant of
c                        aerosol (m2/g)
     $        omgvis,  ! visible omega0
     $        gvis     ! visible forward scattering asymmetry parameter
c
C     integer rhcnst    ! constant relative humidity factor
C
C Relative humidity factor
C
      real(kind=8)  rhfac         ! multiplication factor for kaer
c
c     real(kind=8) a0,            ! constant in rh mult factor
c    $        a1,                 ! constant in rh mult factor
c    $        a2,                 ! constant in rh mult factor
c    $        a3                  ! constant in rh mult factor
c
c     data    a0 / -9.2906106183    /
c     data    a1 /  0.52570211505   /
c     data    a2 / -0.0089285760691 /
c     data    a3 /  5.0877212432e-05/
C
      data    kaervs / 5.3012 /
      data    omgvis / 0.999999 /
      data    gvis   / 0.694889 /
c     data    rhcnst / .80 /
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
        rh(i,k)=dmin1(rhb3d(i,k,j),.99d0) ! needed for aerosols in
C                                           radiation
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

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

      subroutine radout(solin  ,sabtp  ,frsa   ,clrst  ,clrss  ,
     $                  qrs    ,firtp  ,frla   ,clrlt  ,clrls  ,
     $                  qrl, slwd  ,srfrad,
     $                  sols,soll,solsd,solld,
Cadd by bixq
     &     alb,albc,fsds,fsnirt,fsnrtc,fsnirtsq,
Cadd_
     $                  jslc, h2ommr, cld, clwp)
c
c copy radiation output quantities to model buffer
c
c change units of the radiative fluxes from cgs to mks
c
c compute the total radiative heat flux at the surface for
c the surface temperature computation
c
      implicit none
      integer jslc
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/parrad.cb"
#     include "include/bats.cb"
#     include "include/bats2.cb"
#     include "include/rad.cb"
#     include "include/date.cb"

      integer nrad2d,nrad3d
      parameter(nrad2d=21,nrad3d=5)
#ifdef MPP1
#     include "mpiregcm.inc"
      real(kind=4)  frad2d(jxp,ix-2,nrad2d),frad3d(jxp,ix-2,kx,nrad3d)
#else
      real(kind=4)  frad2d(jx-2,ix-2,nrad2d),frad3d(jx-2,ix-2,kx,nrad3d)
#endif
      COMMON /OUTRAD/ frad2d,frad3d
      integer nrcrad
      COMMON /RADREC/ nrcrad
C
C Radiation resolution and I/O parameters
C
c
c input/output arguments
c
      real(kind=8)  solin(plond),     ! instantaneous incident solar
     $     sabtp(plond),     ! total column absorbed solar flux
     $     frsa(plond),      ! surface absorbed solar flux
     $     clrst(plond),     ! clear sky total column abs solar flux
     $     clrss(plond),     ! clear sky surface absorbed solar flux
     $     qrs(plond,plev),  ! solar heating rate
     $     firtp(plond),     ! net up flux top of model (up-dwn flx)
     $     frla(plond),      ! longwave cooling of surface (up-dwn flx)
     $     clrlt(plond),     ! clr sky net up flx top of model (up-dwn f
     $     clrls(plond),     ! clr sky lw cooling of srf (up-dwn flx)
     $     qrl(plond,plev),  ! longwave cooling rate
     $     slwd(plond),      ! surface longwave down flux
     $     srfrad(plond),     ! surface radiative heat flux (frsa+slwd)
     $     h2ommr(plond,plev), ! ozone mixing ratio
     $     cld(plond,plev),   ! cloud fractional cover
     $     clwp(plond,plev),  ! cloud liquid water path
     $     soll(plond),         ! Downward solar rad onto surface (lw direct)
     $     solld(plond),        ! Downward solar rad onto surface (lw diffuse)
     $     sols(plond),         ! Downward solar rad onto surface (sw direct)
     $     solsd(plond)         ! Downward solar rad onto surface (sw diffuse)
c
CEES next 3 added, they are calculated in radcsw
      real(kind=8) fsnirt(plond)      ! Near-IR flux absorbed at toa
      real(kind=8) fsnrtc(plond)      ! Clear sky near-IR flux absorbed at toa
      real(kind=8) fsnirtsq(plond)    ! Near-IR flux absorbed at toa >= 0.7 microns
      real(kind=8) fsds(plond)        ! Flux Shortwave Downwelling Surface
      real(kind=8) alb(plond)
      real(kind=8) albc(plond)
c

c local variables
c
      integer       i    ! latitude index
      integer       k    ! level index
      integer n,nll
c
c      real(kind=8)  cgsmks     ! conversion factor for fluxes from cgs to mks
c
c      data cgsmks / 1.e-3 /
c
c compute total radiative heating flux for the surface,
c converting units from cgs to mks:
c
      do 10 i=1,plon
cKN         srfrad(i) = (frsa(i) + slwd(i)) * cgsmks
         srfrad(i) = frsa(i) + slwd(i)
10    continue
c
c convert units from cgs to mks in solar fluxes:
c
cKN      do 20 i=1,plon
cKN         solin(i) = solin(i) * cgsmks
cKN         sabtp(i) = sabtp(i) * cgsmks
cKN         frsa(i)  = frsa(i)  * cgsmks
cKN         clrst(i) = clrst(i) * cgsmks
cKN         clrss(i) = clrss(i) * cgsmks
cKN20    continue
c
c convert units from cgs to mks in longwave fluxes:
c
cKN      do 30 i=1,plon
cKN         firtp(i) = firtp(i) * cgsmks
cKN         frla(i)  = frla(i)  * cgsmks
cKN         clrlt(i) = clrlt(i) * cgsmks
cKN         clrls(i) = clrls(i) * cgsmks
cKN30    continue
c------
c------total heating rate in deg/s
c------
      do 19 nll=1,plev
      do 19 n=1,plond
       heatrt(n,nll,jslc)=qrs(n,nll)+qrl(n,nll)
19    continue
c------
c------surface absorbed solar flux in watts/m2
c------
      do 24 n=1,plond
       fsw2d(n,jslc)=frsa(n)
24    continue
c------
c------net up longwave flux at the surface
c------
      do 27 n=1,plond
       flw2d(n,jslc)=frla(n)
       flwd2d(n,jslc)=slwd(n)                           ! BATS Output
27    continue
c------
c------for coupling with bats
c------
c   for now assume sabveg (solar absorbed by vegetation) is equal
c   to frsa (solar absorbed by surface). possible problems are
c   over sparsely vegetated areas in which vegetation and ground
c   albedo are significantly different
      do 44 n=1,plond
       sabv2d(n,jslc)=sabveg(n)
       sol2d(n,jslc)=solis(n)
       sinc2d(n,jslc)=soll(n)+sols(n)+solsd(n)+solld(n)
       solvs2d(n,jslc)=solvs(n)
       solvd2d(n,jslc)=solvd(n)
c       sinc2d(n,jslc)=solin(n)
44    continue
c
      if (ifrad) then
        if (mod(ntime+nint(dtmin*60.),nradisp).eq.0
     &      .or. (jyear.eq.jyearr.and.ktau.eq.ktaur)) then
          do k=1,kx
          do i=2,ix-1
#ifdef MPP1
            frad3d(jslc,i-1,k,1)=h2ommr(i,k)  ! skip
            frad3d(jslc,i-1,k,2)=   cld(i,k)
            frad3d(jslc,i-1,k,3)=  clwp(i,k)
            frad3d(jslc,i-1,k,4)=   qrs(i,k)
            frad3d(jslc,i-1,k,5)=   qrl(i,k)
#else
            frad3d(jslc-1,i-1,k,1)=h2ommr(i,k)  ! skip
            frad3d(jslc-1,i-1,k,2)=   cld(i,k)
            frad3d(jslc-1,i-1,k,3)=  clwp(i,k)
            frad3d(jslc-1,i-1,k,4)=   qrs(i,k)
            frad3d(jslc-1,i-1,k,5)=   qrl(i,k)
#endif
          end do
          end do

          do i=2,ix-1
#ifdef MPP1
            frad2d(jslc,i-1,1)  =    frsa(i)  ! write
            frad2d(jslc,i-1,2)  =    frla(i)  ! write
            frad2d(jslc,i-1,3)  =   clrst(i)  ! write
            frad2d(jslc,i-1,4)  =   clrss(i)  ! write
            frad2d(jslc,i-1,5)  =   clrlt(i)  ! write
            frad2d(jslc,i-1,6)  =   clrls(i)  ! write
            frad2d(jslc,i-1,7)  =   solin(i)  ! write
            frad2d(jslc,i-1,8)  =   sabtp(i)  ! write
            frad2d(jslc,i-1,9)  =   firtp(i)  ! write
            frad2d(jslc,i-1,10) =     alb(i)  ! skip
            frad2d(jslc,i-1,11) =    albc(i)  ! skip
            frad2d(jslc,i-1,12) =    fsds(i)  ! skip
            frad2d(jslc,i-1,13) =  fsnirt(i)  ! skip
            frad2d(jslc,i-1,14) =  fsnrtc(i)  ! skip
            frad2d(jslc,i-1,15) =fsnirtsq(i)  ! skip
            frad2d(jslc,i-1,16) =    soll(i)  ! skip
            frad2d(jslc,i-1,17) =    sols(i)  ! skip
            frad2d(jslc,i-1,18) =   solsd(i)  ! skip
            frad2d(jslc,i-1,19) =   solld(i)  ! skip
            frad2d(jslc,i-1,20) =   solis(i)  ! skip
            frad2d(jslc,i-1,21) =  sabveg(i)  ! skip
#else
            frad2d(jslc-1,i-1,1)  =    frsa(i)  ! write
            frad2d(jslc-1,i-1,2)  =    frla(i)  ! write
            frad2d(jslc-1,i-1,3)  =   clrst(i)  ! write
            frad2d(jslc-1,i-1,4)  =   clrss(i)  ! write
            frad2d(jslc-1,i-1,5)  =   clrlt(i)  ! write
            frad2d(jslc-1,i-1,6)  =   clrls(i)  ! write
            frad2d(jslc-1,i-1,7)  =   solin(i)  ! write
            frad2d(jslc-1,i-1,8)  =   sabtp(i)  ! write
            frad2d(jslc-1,i-1,9)  =   firtp(i)  ! write
            frad2d(jslc-1,i-1,10) =     alb(i)  ! skip
            frad2d(jslc-1,i-1,11) =    albc(i)  ! skip
            frad2d(jslc-1,i-1,12) =    fsds(i)  ! skip
            frad2d(jslc-1,i-1,13) =  fsnirt(i)  ! skip
            frad2d(jslc-1,i-1,14) =  fsnrtc(i)  ! skip
            frad2d(jslc-1,i-1,15) =fsnirtsq(i)  ! skip
            frad2d(jslc-1,i-1,16) =    soll(i)  ! skip
            frad2d(jslc-1,i-1,17) =    sols(i)  ! skip
            frad2d(jslc-1,i-1,18) =   solsd(i)  ! skip
            frad2d(jslc-1,i-1,19) =   solld(i)  ! skip
            frad2d(jslc-1,i-1,20) =   solis(i)  ! skip
            frad2d(jslc-1,i-1,21) =  sabveg(i)  ! skip
#endif
          end do
        end if
      end if
      return
      end subroutine radout

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_outrad

      use mod_dynparam
      use mod_runparams
      use mod_bats
      use mod_rad
      use mod_date
#ifdef CLM
      use mod_clm
#endif

      private

      public :: frad2d , frad3d
      public :: allocate_mod_outrad , radout

      real(4) ,allocatable, dimension(:,:,:) :: frad2d
      real(4) ,allocatable, dimension(:,:,:,:) :: frad3d
#ifndef MPP1
      real(4) ,allocatable, dimension(:,:) :: radpsa
      public :: radpsa
#endif

      contains

      subroutine allocate_mod_outrad
        implicit none
#ifdef MPP1
        allocate(frad2d(jxp,iym2,nrad2d))
        allocate(frad3d(jxp,iym2,kz,nrad3d))
#else
#ifdef BAND
        allocate(frad2d(jx,iym2,nrad2d))
        allocate(frad3d(jx,iym2,kz,nrad3d))
        allocate(radpsa(jx,iym2))
#else
        allocate(frad2d(jxm2,iym2,nrad2d))
        allocate(frad3d(jxm2,iym2,kz,nrad3d))
        allocate(radpsa(jxm2,iym2))
#endif
        radpsa = 0.0D0
#endif
        frad2d = 0.0D0
        frad3d = 0.0D0
        end subroutine allocate_mod_outrad
! 
        subroutine radout(solin,sabtp,frsa,clrst,clrss,qrs,firtp,frla, &
                          & clrlt,clrls,qrl,slwd,srfrad,sols,soll,     &
                          & solsd,solld,alb,albc,fsds,fsnirt,fsnrtc,   &
                          & fsnirtsq,jslc,h2ommr,cld,clwp)
!
! copy radiation output quantities to model buffer
!
! change units of the radiative fluxes from cgs to mks
!
! compute the total radiative heat flux at the surface for
! the surface temperature computation
!
          implicit none
!
!     input/output arguments
!
! solin  - instantaneous incident solar
! sabtp  - total column absorbed solar flux
! frsa   - surface absorbed solar flux
! clrst  - clear sky total column abs solar flux
! clrss  - clear sky surface absorbed solar flux
! qrs    - solar heating rate
! firtp  - net up flux top of model (up-dwn flx)
! frla   - longwave cooling of surface (up-dwn flx)
! clrlt  - clr sky net up flx top of model (up-dwn f
! clrls  - clr sky lw cooling of srf (up-dwn flx)
! qrl    - longwave cooling rate
! slwd   - surface longwave down flux
! srfrad - surface radiative heat flux (frsa+slwd)
! h2ommr - ozone mixing ratio
! cld    - cloud fractional cover
! clwp   - cloud liquid water path
! soll   - Downward solar rad onto surface (lw direct)
! solld  - Downward solar rad onto surface (lw diffuse)
! sols   - Downward solar rad onto surface (sw direct)
! solsd  - Downward solar rad onto surface (sw diffuse)
!
!EES  next 3 added, they are calculated in radcsw
! fsnirt   - Near-IR flux absorbed at toa
! fsnrtc   - Clear sky near-IR flux absorbed at toa
! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
! fsds     - Flux Shortwave Downwelling Surface
!
 
!
! Dummy arguments
!
          integer :: jslc
          real(8) , dimension(iym1) :: alb , albc , clrls , clrlt , &
                      clrss , clrst , firtp , frla , frsa , fsds ,  &
                      fsnirt , fsnirtsq , fsnrtc , sabtp , slwd ,   &
                      solin , soll , solld , sols , solsd , srfrad
          real(8) , dimension(iym1,kz) :: cld , clwp , h2ommr , &
                                          qrl , qrs
          intent (in) alb , albc , cld , clrls , clrlt , clrss , clrst ,&
                    & clwp , firtp , frla , frsa , fsds , fsnirt ,      &
                    & fsnirtsq , fsnrtc , h2ommr , jslc , qrl , qrs ,   &
                    & sabtp , slwd , solin , soll , solld , sols , solsd
          intent (out) srfrad
!
! Local variables
!
          integer :: i , k , n , nll
!
!     compute total radiative heating flux for the surface,
!     converting units from cgs to mks:
!
          do i = 1 , iym1    ! level index
!KN     srfrad(i) = (frsa(i) + slwd(i)) * cgsmks
            srfrad(i) = frsa(i) + slwd(i)
          end do
!
!     convert units from cgs to mks in solar fluxes:
!
!KN   do 20 i=1,iym1
!KN   solin(i) = solin(i) * cgsmks
!KN   sabtp(i) = sabtp(i) * cgsmks
!KN   frsa(i)  = frsa(i)  * cgsmks
!KN   clrst(i) = clrst(i) * cgsmks
!KN   clrss(i) = clrss(i) * cgsmks
!KN20 continue
!
!     convert units from cgs to mks in longwave fluxes:
!
!KN   do 30 i=1,iym1
!KN   firtp(i) = firtp(i) * cgsmks
!KN   frla(i)  = frla(i)  * cgsmks
!KN   clrlt(i) = clrlt(i) * cgsmks
!KN   clrls(i) = clrls(i) * cgsmks
!KN30 continue
!------
!------total heating rate in deg/s
!------
          do nll = 1 , kz
            do n = 1 , iym1
              heatrt(n,nll,jslc) = qrs(n,nll) + qrl(n,nll)
            end do
          end do
!------
!------surface absorbed solar flux in watts/m2
!------
          do n = 1 , iym1
            fsw2d(n,jslc) = frsa(n)
          end do
!------
!------net up longwave flux at the surface
!------
          do n = 1 , iym1
            flw2d(n,jslc) = frla(n)
            flwd2d(n,jslc) = slwd(n)  ! BATS Output
          end do
!------
!------for coupling with bats
!------
!     for now assume sabveg (solar absorbed by vegetation) is equal
!     to frsa (solar absorbed by surface). possible problems are
!     over sparsely vegetated areas in which vegetation and ground
!     albedo are significantly different
          do n = 1 , iym1
            sabv2d(n,jslc) = sabveg(n)
            sol2d(n,jslc) = solis(n)
            sinc2d(n,jslc) = soll(n) + sols(n) + solsd(n) + solld(n)
            solvs2d(n,jslc) = solvs(n)
            solvd2d(n,jslc) = solvd(n)
!           sinc2d(n,jslc)=solin(n)
#ifdef CLM
            sols2d(n,jslc) = sols(n)
            soll2d(n,jslc) = soll(n)
            solsd2d(n,jslc) = solsd(n)
            solld2d(n,jslc) = solld(n)
#endif
          end do
!
          if ( ifrad ) then
            if ( mod(ntime+idnint(dtmin*60.0D0),nradisp).eq.0 .or.  &
                ( jyear.eq.jyear0 .and. ktau.eq.0 ) .or.         &
                ( ifrest .and. .not. done_restart) ) then
              do k = 1 , kz
                do i = 2 , iym1
#ifdef MPP1
                  frad3d(jslc,i-1,k,1) = h2ommr(i,k)
                  frad3d(jslc,i-1,k,2) = cld(i,k)
                  if (clwp(i,k) > 1D-30) then
                    frad3d(jslc,i-1,k,3) = clwp(i,k)
                  else
                    frad3d(jslc,i-1,k,3) = 1D-30
                  end if
                  frad3d(jslc,i-1,k,4) = qrs(i,k)
                  frad3d(jslc,i-1,k,5) = qrl(i,k)
#else
#ifdef BAND
                  frad3d(jslc,i-1,k,1) = h2ommr(i,k)
                  frad3d(jslc,i-1,k,2) = cld(i,k)
                  if (clwp(i,k) > 1D-30) then
                    frad3d(jslc,i-1,k,3) = clwp(i,k)
                  else
                    frad3d(jslc,i-1,k,3) = 1D-30
                  end if
                  frad3d(jslc,i-1,k,4) = qrs(i,k)
                  frad3d(jslc,i-1,k,5) = qrl(i,k)
#else
                  frad3d(jslc-1,i-1,k,1) = h2ommr(i,k)
                  frad3d(jslc-1,i-1,k,2) = cld(i,k)
                  if (clwp(i,k) > 1D-30) then
                    frad3d(jslc-1,i-1,k,3) = clwp(i,k)
                  else
                    frad3d(jslc-1,i-1,k,3) = 1D-30
                  end if
                  frad3d(jslc-1,i-1,k,4) = qrs(i,k)
                  frad3d(jslc-1,i-1,k,5) = qrl(i,k)
#endif
#endif
                end do
              end do
 
              do i = 2 , iym1
#ifdef MPP1
                frad2d(jslc,i-1,1) = frsa(i)      ! write
                frad2d(jslc,i-1,2) = frla(i)      ! write
                frad2d(jslc,i-1,3) = clrst(i)     ! write
                frad2d(jslc,i-1,4) = clrss(i)     ! write
                frad2d(jslc,i-1,5) = clrlt(i)     ! write
                frad2d(jslc,i-1,6) = clrls(i)     ! write
                frad2d(jslc,i-1,7) = solin(i)     ! write
                frad2d(jslc,i-1,8) = sabtp(i)     ! write
                frad2d(jslc,i-1,9) = firtp(i)     ! write
                frad2d(jslc,i-1,10) = alb(i)      ! skip
                frad2d(jslc,i-1,11) = albc(i)     ! skip
                frad2d(jslc,i-1,12) = fsds(i)     ! skip
                frad2d(jslc,i-1,13) = fsnirt(i)   ! skip
                frad2d(jslc,i-1,14) = fsnrtc(i)   ! skip
                frad2d(jslc,i-1,15) = fsnirtsq(i) ! skip
                if ( soll(i) .lt. 1D-30 ) then
                  frad2d(jslc,i-1,16) = 0.0D0
                else
                  frad2d(jslc,i-1,16) = soll(i)     ! skip
                end if
                if ( sols(i) .lt. 1D-30 ) then
                  frad2d(jslc,i-1,17) = 0.0D0
                else
                  frad2d(jslc,i-1,17) = sols(i)     ! skip
                end if
                if ( solsd(i) .lt. 1D-30 ) then
                  frad2d(jslc,i-1,18) = 0.0D0
                else
                  frad2d(jslc,i-1,18) = solsd(i)     ! skip
                end if
                if ( solld(i) .lt. 1D-30 ) then
                  frad2d(jslc,i-1,19) = 0.0D0
                else
                  frad2d(jslc,i-1,19) = solld(i)     ! skip
                end if
                frad2d(jslc,i-1,20) = solis(i)    ! skip
                frad2d(jslc,i-1,21) = sabveg(i)   ! skip
#else
#ifdef BAND
                frad2d(jslc,i-1,1) = frsa(i)      ! write
                frad2d(jslc,i-1,2) = frla(i)      ! write
                frad2d(jslc,i-1,3) = clrst(i)     ! write
                frad2d(jslc,i-1,4) = clrss(i)     ! write
                frad2d(jslc,i-1,5) = clrlt(i)     ! write
                frad2d(jslc,i-1,6) = clrls(i)     ! write
                frad2d(jslc,i-1,7) = solin(i)     ! write
                frad2d(jslc,i-1,8) = sabtp(i)     ! write
                frad2d(jslc,i-1,9) = firtp(i)     ! write
                frad2d(jslc,i-1,10) = alb(i)      ! skip
                frad2d(jslc,i-1,11) = albc(i)     ! skip
                frad2d(jslc,i-1,12) = fsds(i)     ! skip
                frad2d(jslc,i-1,13) = fsnirt(i)   ! skip
                frad2d(jslc,i-1,14) = fsnrtc(i)   ! skip
                frad2d(jslc,i-1,15) = fsnirtsq(i) ! skip
                if ( soll(i) .lt. 1D-30 ) then
                  frad2d(jslc,i-1,16) = 0.0D0
                else
                  frad2d(jslc,i-1,16) = soll(i)     ! skip
                end if
                if ( sols(i) .lt. 1D-30 ) then
                  frad2d(jslc,i-1,17) = 0.0D0
                else
                  frad2d(jslc,i-1,17) = sols(i)     ! skip
                end if
                if ( solsd(i) .lt. 1D-30 ) then
                  frad2d(jslc,i-1,18) = 0.0D0
                else
                  frad2d(jslc,i-1,18) = solsd(i)     ! skip
                end if
                if ( solld(i) .lt. 1D-30 ) then
                  frad2d(jslc,i-1,19) = 0.0D0
                else
                  frad2d(jslc,i-1,19) = solld(i)     ! skip
                end if
                frad2d(jslc,i-1,20) = solis(i)    ! skip
                frad2d(jslc,i-1,21) = sabveg(i)   ! skip
#else
                frad2d(jslc-1,i-1,1) = frsa(i)      ! write
                frad2d(jslc-1,i-1,2) = frla(i)      ! write
                frad2d(jslc-1,i-1,3) = clrst(i)     ! write
                frad2d(jslc-1,i-1,4) = clrss(i)     ! write
                frad2d(jslc-1,i-1,5) = clrlt(i)     ! write
                frad2d(jslc-1,i-1,6) = clrls(i)     ! write
                frad2d(jslc-1,i-1,7) = solin(i)     ! write
                frad2d(jslc-1,i-1,8) = sabtp(i)     ! write
                frad2d(jslc-1,i-1,9) = firtp(i)     ! write
                frad2d(jslc-1,i-1,10) = alb(i)      ! skip
                frad2d(jslc-1,i-1,11) = albc(i)     ! skip
                frad2d(jslc-1,i-1,12) = fsds(i)     ! skip
                frad2d(jslc-1,i-1,13) = fsnirt(i)   ! skip
                frad2d(jslc-1,i-1,14) = fsnrtc(i)   ! skip
                frad2d(jslc-1,i-1,15) = fsnirtsq(i) ! skip
                if ( soll(i) .lt. 1D-30 ) then
                  frad2d(jslc-1,i-1,16) = 0.0D0
                else
                  frad2d(jslc-1,i-1,16) = soll(i)     ! skip
                end if
                if ( sols(i) .lt. 1D-30 ) then
                  frad2d(jslc-1,i-1,17) = 0.0D0
                else
                  frad2d(jslc-1,i-1,17) = sols(i)     ! skip
                end if
                if ( solsd(i) .lt. 1D-30 ) then
                  frad2d(jslc-1,i-1,18) = 0.0D0
                else
                  frad2d(jslc-1,i-1,18) = solsd(i)     ! skip
                end if
                if ( solld(i) .lt. 1D-30 ) then
                  frad2d(jslc-1,i-1,19) = 0.0D0
                else
                  frad2d(jslc-1,i-1,19) = solld(i)     ! skip
                end if
                frad2d(jslc-1,i-1,20) = solis(i)    ! skip
                frad2d(jslc-1,i-1,21) = sabveg(i)   ! skip
#endif
#endif
              end do
            end if
          end if
        end subroutine radout

      end module mod_outrad

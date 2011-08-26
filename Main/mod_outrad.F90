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

  use mod_runparams
  use mod_lm_interface
  use mod_rad
#ifdef CLM
  use mod_clm
#endif

  private

  public :: frad2d , frad3d
  public :: allocate_mod_outrad , radout
  public :: nrad2d , nrad3d

  integer , parameter :: nrad2d = 21
  integer , parameter :: nrad3d = 5

  real(4) , pointer , dimension(:,:,:) :: frad2d
  real(4) , pointer , dimension(:,:,:,:) :: frad3d

  contains

  subroutine allocate_mod_outrad
    implicit none
    call getmem3d(frad2d,1,jxp,1,iym2,1,nrad2d,'mod_outrad:frad2d')
    call getmem4d(frad3d,1,jxp,1,iym2,1,kz,1,nrad3d,'mod_outrad:frad3d')
  end subroutine allocate_mod_outrad
! 
  subroutine radout(solin,sabtp,frsa,clrst,clrss,qrs,firtp,frla, &
                        clrlt,clrls,qrl,slwd,srfrad,sols,soll,     &
                        solsd,solld,alb,albc,fsds,fsnirt,fsnrtc,   &
                        fsnirtsq,j,h2ommr,cld,clwp)
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
    integer :: j
    real(8) , dimension(iym1) :: alb , albc , clrls , clrlt , &
                clrss , clrst , firtp , frla , frsa , fsds ,  &
                fsnirt , fsnirtsq , fsnrtc , sabtp , slwd ,   &
                solin , soll , solld , sols , solsd , srfrad
    real(8) , dimension(iym1,kz) :: cld , clwp , h2ommr , qrl , qrs
    intent (in) alb , albc , cld , clrls , clrlt , clrss , clrst ,&
                clwp , firtp , frla , frsa , fsds , fsnirt ,      &
                fsnirtsq , fsnrtc , h2ommr , j , qrl , qrs ,      &
                sabtp , slwd , solin , soll , solld , sols , solsd
    intent (out) srfrad
!
    integer :: i , k
!
!     compute total radiative heating flux for the surface
!
    do i = 1 , iym1    ! level index
      srfrad(i) = frsa(i) + slwd(i)
    end do
!------
!------total heating rate in deg/s
!------
    do k = 1 , kz
      do i = 1 , iym1
        heatrt(i,k,j) = qrs(i,k) + qrl(i,k)
      end do
    end do
!------
!------surface absorbed solar flux in watts/m2
!------
    do i = 1 , iym1
      fsw2d(i,j) = frsa(i)
    end do
!------
!------net up longwave flux at the surface
!------
    do i = 1 , iym1
      flw2d(i,j) = frla(i)
      flwd2d(i,j) = slwd(i)  ! BATS Output
    end do
!------
!------for coupling with bats
!------
!     for now assume sabveg (solar absorbed by vegetation) is equal
!     to frsa (solar absorbed by surface). possible problems are
!     over sparsely vegetated areas in which vegetation and ground
!     albedo are significantly different
    do i = 1 , iym1
      sabv2d(i,j) = sabveg(i)
      sol2d(i,j) = solis(i)
      sinc2d(i,j) = soll(i) + sols(i) + solsd(i) + solld(i)
      solvs2d(i,j) = solvs(i)
      solvd2d(i,j) = solvd(i)
#ifdef CLM
      sols2d(i,j) = sols(i)
      soll2d(i,j) = soll(i)
      solsd2d(i,j) = solsd(i)
      solld2d(i,j) = solld(i)
#endif
    end do
!
    if ( ifrad ) then
      if ( mod(nbdytime+ntsec,nradfrq) == 0 .or.  &
          ktau == 0 .or. doing_restart ) then
        do k = 1 , kz
          do i = 2 , iym1
            frad3d(j,i-1,k,1) = real(h2ommr(i,k))
            frad3d(j,i-1,k,2) = real(cld(i,k))
            if (clwp(i,k) > dlowval) then
              frad3d(j,i-1,k,3) = real(clwp(i,k))
            else
              frad3d(j,i-1,k,3) = slowval
            end if
            frad3d(j,i-1,k,4) = real(qrs(i,k))
            frad3d(j,i-1,k,5) = real(qrl(i,k))
          end do
        end do
 
        do i = 2 , iym1
          frad2d(j,i-1,1) = real(frsa(i))      ! write
          frad2d(j,i-1,2) = real(frla(i))      ! write
          frad2d(j,i-1,3) = real(clrst(i))     ! write
          frad2d(j,i-1,4) = real(clrss(i))     ! write
          frad2d(j,i-1,5) = real(clrlt(i))     ! write
          frad2d(j,i-1,6) = real(clrls(i))     ! write
          frad2d(j,i-1,7) = real(solin(i))     ! write
          frad2d(j,i-1,8) = real(sabtp(i))     ! write
          frad2d(j,i-1,9) = real(firtp(i))     ! write
          frad2d(j,i-1,10) = real(alb(i))      ! skip
          frad2d(j,i-1,11) = real(albc(i))     ! skip
          frad2d(j,i-1,12) = real(fsds(i))     ! skip
          frad2d(j,i-1,13) = real(fsnirt(i))   ! skip
          frad2d(j,i-1,14) = real(fsnrtc(i))   ! skip
          frad2d(j,i-1,15) = real(fsnirtsq(i)) ! skip
          if ( soll(i) < dlowval ) then
            frad2d(j,i-1,16) = 0.0
          else
            frad2d(j,i-1,16) = real(soll(i))   ! skip
          end if
          if ( sols(i) < dlowval ) then
            frad2d(j,i-1,17) = 0.0
          else
            frad2d(j,i-1,17) = real(sols(i))   ! skip
          end if
          if ( solsd(i) < dlowval ) then
            frad2d(j,i-1,18) = 0.0
          else
            frad2d(j,i-1,18) = real(solsd(i))  ! skip
          end if
          if ( solld(i) < dlowval ) then
            frad2d(j,i-1,19) = 0.0
          else
            frad2d(j,i-1,19) = real(solld(i))  ! skip
          end if
          frad2d(j,i-1,20) = real(solis(i))    ! skip
          frad2d(j,i-1,21) = real(sabveg(i))   ! skip
        end do
      end if
    end if
  end subroutine radout

end module mod_outrad

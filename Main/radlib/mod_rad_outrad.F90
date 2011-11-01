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

module mod_rad_outrad

  use mod_dynparam
  use mod_rad_common

  private

  public :: frad2d , frad3d
  public :: allocate_mod_rad_outrad , radout
  public :: nrad2d , nrad3d

  integer , parameter :: nrad2d = 21
  integer , parameter :: nrad3d = 5

  real(4) , pointer , dimension(:,:,:) :: frad2d
  real(4) , pointer , dimension(:,:,:,:) :: frad3d

  contains

  subroutine allocate_mod_rad_outrad
    implicit none
    call getmem3d(frad2d,1,jxp,1,iym2,1,nrad2d,'mod_outrad:frad2d')
    call getmem4d(frad3d,1,jxp,1,iym2,1,kz,1,nrad3d,'mod_outrad:frad3d')
  end subroutine allocate_mod_rad_outrad
! 
  subroutine radout(jstart,jend,i,lout,solin,sabtp,frsa,clrst,clrss, &
                    qrs,firtp,frla,clrlt,clrls,qrl,slwd,srfrad,sols, &
                    soll,solsd,solld,alb,albc,fsds,fsnirt,fsnrtc,    &
                    fsnirtsq,h2ommr,cld,clwp)
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
    integer , intent(in) :: jstart , jend , i
    logical , intent(in) :: lout ! Preapre data for outfile
    real(8) , pointer , dimension(:) :: alb , albc , clrls , clrlt ,  &
                clrss , clrst , firtp , frla , frsa , fsds , fsnirt , &
                fsnirtsq , fsnrtc , sabtp , slwd , solin , soll ,     &
                solld , sols , solsd , srfrad
    real(8) , pointer , dimension(:,:) :: cld , clwp , h2ommr , qrl , qrs
    intent (in) alb , albc , cld , clrls , clrlt , clrss , clrst ,&
                clwp , firtp , frla , frsa , fsds , fsnirt ,      &
                fsnirtsq , fsnrtc , h2ommr , qrl , qrs , sabtp ,  &
                slwd , solin , soll , solld , sols , solsd
    intent (out) srfrad
!
    integer :: j , k
    logical , save :: firstin
!
    data firstin /.true./
!
!     compute total radiative heating flux for the surface
!
    do j = jstart , jend
      srfrad(j) = frsa(j) + slwd(j)
    end do
!------
!------total heating rate in deg/s
!------
    do k = 1 , kz
      do j = jstart , jend
        heatrt(j,i,k) = qrs(j,k) + qrl(j,k)
      end do
    end do
!------
!------surface absorbed solar flux in watts/m2
!------
    do j = jstart , jend
      srfabswflx(i,j) = frsa(j)
    end do
!------
!------net up longwave flux at the surface
!------
    do j = jstart , jend
      srflwflxup(i,j) = frla(j)
      srflwflxdw(i,j) = slwd(j)  ! BATS Output
    end do
!------
!------for coupling with bats
!------
!     for now assume abveg (solar absorbed by vegetation) is equal
!     to frsa (solar absorbed by surface). possible problems are
!     over sparsely vegetated areas in which vegetation and ground
!     albedo are significantly different
    do j = jstart , jend
      abveg2d(i,j) = abveg(j,i)
      solar2d(i,j) = solar(j,i)
      totsol2d(i,j) = soll(j) + sols(j) + solsd(j) + solld(j)
      soldir2d(i,j) = sols(j)
      soldif2d(i,j) = solsd(j)
#ifdef CLM
      solswdir(i,j) = sols(j)
      sollwdir(i,j) = soll(j)
      solswdif(i,j) = solsd(j)
      sollwdif(i,j) = solld(j)
#endif
    end do
!
    if ( ifrad ) then
      if ( lout .or. firstin ) then
        do k = 1 , kz
          do j = jstart , jend
            frad3d(j,i-1,k,1) = real(h2ommr(j,k))
            frad3d(j,i-1,k,2) = real(cld(j,k))
            if (clwp(j,k) > dlowval) then
              frad3d(j,i-1,k,3) = real(clwp(j,k))
            else
              frad3d(j,i-1,k,3) = slowval
            end if
            frad3d(j,i-1,k,4) = real(qrs(j,k))
            frad3d(j,i-1,k,5) = real(qrl(j,k))
          end do
        end do
 
        do j = jstart , jend
          frad2d(j,i-1,1) = real(frsa(j))      ! write
          frad2d(j,i-1,2) = real(frla(j))      ! write
          frad2d(j,i-1,3) = real(clrst(j))     ! write
          frad2d(j,i-1,4) = real(clrss(j))     ! write
          frad2d(j,i-1,5) = real(clrlt(j))     ! write
          frad2d(j,i-1,6) = real(clrls(j))     ! write
          frad2d(j,i-1,7) = real(solin(j))     ! write
          frad2d(j,i-1,8) = real(sabtp(j))     ! write
          frad2d(j,i-1,9) = real(firtp(j))     ! write
          frad2d(j,i-1,10) = real(alb(j))      ! skip
          frad2d(j,i-1,11) = real(albc(j))     ! skip
          frad2d(j,i-1,12) = real(fsds(j))     ! skip
          frad2d(j,i-1,13) = real(fsnirt(j))   ! skip
          frad2d(j,i-1,14) = real(fsnrtc(j))   ! skip
          frad2d(j,i-1,15) = real(fsnirtsq(j)) ! skip
          if ( soll(j) < dlowval ) then
            frad2d(j,i-1,16) = 0.0
          else
            frad2d(j,i-1,16) = real(soll(j))   ! skip
          end if
          if ( sols(j) < dlowval ) then
            frad2d(j,i-1,17) = 0.0
          else
            frad2d(j,i-1,17) = real(sols(j))   ! skip
          end if
          if ( solsd(j) < dlowval ) then
            frad2d(j,i-1,18) = 0.0
          else
            frad2d(j,i-1,18) = real(solsd(j))  ! skip
          end if
          if ( solld(j) < dlowval ) then
            frad2d(j,i-1,19) = 0.0
          else
            frad2d(j,i-1,19) = real(solld(j))  ! skip
          end if
          frad2d(j,i-1,20) = real(solar(j,i))    ! skip
          frad2d(j,i-1,21) = real(abveg(j,i))   ! skip
        end do
      end if
    end if
    firstin = .false.
  end subroutine radout

end module mod_rad_outrad

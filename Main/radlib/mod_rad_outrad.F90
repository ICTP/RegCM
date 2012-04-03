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

  integer , parameter :: nrad2d = 24
  integer , parameter :: nrad3d = 5

  real(sp) , pointer , dimension(:,:,:) :: frad2d
  real(sp) , pointer , dimension(:,:,:,:) :: frad3d

  contains

  subroutine allocate_mod_rad_outrad
    implicit none
    call getmem3d(frad2d,jci1,jci2,ici1,ici2,1,nrad2d,'mod_outrad:frad2d')
    call getmem4d(frad3d,jci1,jci2,ici1,ici2,1,kz,1,nrad3d,'mod_outrad:frad3d')
  end subroutine allocate_mod_rad_outrad
!
  subroutine radout(jstart,jend,i,lout,solin,sabtp,frsa,            &
                    clrst,clrss,qrs,firtp,frla,clrlt,clrls,qrl,     &
                    slwd,sols,soll,solsd,solld,alb,albc,            &
                    fsds,fsnirt,fsnrtc,fsnirtsq,totcf,totcl,totci,  &
                    h2ommr,cld,clwp,aeradfo,aeradfos,aerlwfo,       &
                    aerlwfos,tauxar3d,tauasc3d,gtota3d)
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
    real(dp) , pointer , dimension(:) :: alb , albc , clrls , clrlt ,  &
                clrss , clrst , firtp , frla , frsa , fsds , fsnirt ,  &
                fsnirtsq , fsnrtc , sabtp , slwd , solin , soll ,      &
                solld , sols , solsd , totcf , totcl , totci
    real(dp) , pointer , dimension(:,:) :: cld , clwp , h2ommr , qrl , qrs
    real(dp) , pointer , dimension(:,:,:) :: tauxar3d , tauasc3d , gtota3d
    real(dp) , pointer , dimension(:) :: aeradfo , aeradfos, aerlwfo , aerlwfos
    intent (in) alb , albc , cld , clrls , clrlt , clrss , clrst ,&
                clwp , firtp , frla , frsa , fsds , fsnirt ,      &
                fsnirtsq , fsnrtc , h2ommr , qrl , qrs , sabtp ,  &
                slwd , solin , soll , solld , sols , solsd ,      &
                totcf , totcl , totci , aeradfo , aeradfos,       &
                aerlwfo , aerlwfos
!
    integer :: j , k
    real(dp) :: rntim
!
!   total heating rate in deg/s
!
    do k = 1 , kz
      do j = jstart , jend
        heatrt(j,i,k) = qrs(j,k) + qrl(j,k)
      end do
    end do
!
!   surface absorbed solar flux in watts/m2
!
    do j = jstart , jend
      srfabswflx(j,i) = frsa(j)
    end do
!
!   net up longwave flux at the surface
!
    do j = jstart , jend
      srflwflxup(j,i) = frla(j)
      srflwflxdw(j,i) = slwd(j)  ! BATS Output
    end do
!
!   for coupling with bats
!
!   for now assume abveg (solar absorbed by vegetation) is equal
!   to frsa (solar absorbed by surface). possible problems are
!   over sparsely vegetated areas in which vegetation and ground
!   albedo are significantly different
!
    do j = jstart , jend
      totsol(j,i) = soll(j) + sols(j) + solsd(j) + solld(j)
      soldir(j,i) = sols(j)
      soldif(j,i) = solsd(j)
#ifdef CLM
      solswdir(j,i) = sols(j)
      sollwdir(j,i) = soll(j)
      solswdif(j,i) = solsd(j)
      sollwdif(j,i) = solld(j)
#endif
    end do
    if ( lchem ) then
      do k = 1 , kz
        do j = jstart , jend
          aerext(j,i,k) = tauxar3d(j,k,8)
          aerssa(j,i,k) = tauasc3d(j,k,8)
          aerasp(j,i,k) = gtota3d(j,k,8)
        end do
      end do
      rntim = d_one/(d_1000*minph*chfrovrradfr)
      do j = jstart , jend
        aertarf(j,i)   = aertarf(j,i)   + aeradfo(j)  * rntim
        aersrrf(j,i)   = aersrrf(j,i)   + aeradfos(j) * rntim
        aertalwrf(j,i) = aertalwrf(j,i) + aerlwfo(j)  * rntim
        aersrlwrf(j,i) = aersrlwrf(j,i) + aerlwfos(j) * rntim
      end do
    end if
!
    if ( ifrad ) then
      if ( lout ) then
        do k = 1 , kz
          do j = jstart , jend
            frad3d(j,i,k,1) = real(cld(j,k))    ! write
            if (clwp(j,k) > dlowval) then
              frad3d(j,i,k,2) = real(clwp(j,k)) ! write
            else
              frad3d(j,i,k,2) = slowval
            end if
            frad3d(j,i,k,3) = real(qrs(j,k))    ! write
            frad3d(j,i,k,4) = real(qrl(j,k))    ! write
            frad3d(j,i,k,5) = real(h2ommr(j,k)) ! skip
          end do
        end do
   
        do j = jstart , jend
          frad2d(j,i,1) = real(frsa(j))      ! write
          frad2d(j,i,2) = real(frla(j))      ! write
          frad2d(j,i,3) = real(clrst(j))     ! write
          frad2d(j,i,4) = real(clrss(j))     ! write
          frad2d(j,i,5) = real(clrlt(j))     ! write
          frad2d(j,i,6) = real(clrls(j))     ! write
          frad2d(j,i,7) = real(solin(j))     ! write
          frad2d(j,i,8) = real(sabtp(j))     ! write
          frad2d(j,i,9) = real(totcf(j))     ! write
          frad2d(j,i,10) = real(totcl(j))    ! write
          frad2d(j,i,11) = real(totci(j))    ! write
          frad2d(j,i,12) = real(firtp(j))    ! write
          frad2d(j,i,13) = real(alb(j))      ! skip
          frad2d(j,i,14) = real(albc(j))     ! skip
          frad2d(j,i,15) = real(fsds(j))     ! skip
          frad2d(j,i,16) = real(fsnirt(j))   ! skip
          frad2d(j,i,17) = real(fsnrtc(j))   ! skip
          frad2d(j,i,18) = real(fsnirtsq(j)) ! skip
          if ( soll(j) < dlowval ) then
            frad2d(j,i,19) = 0.0
          else
            frad2d(j,i,19) = real(soll(j))   ! skip
          end if
          if ( sols(j) < dlowval ) then
            frad2d(j,i,20) = 0.0
          else
            frad2d(j,i,20) = real(sols(j))   ! skip
          end if
          if ( solsd(j) < dlowval ) then
            frad2d(j,i,21) = 0.0
          else
            frad2d(j,i,21) = real(solsd(j))  ! skip
          end if
          if ( solld(j) < dlowval ) then
            frad2d(j,i,22) = 0.0
          else
            frad2d(j,i,22) = real(solld(j))  ! skip
          end if
          frad2d(j,i,23) = real(solar(j,i))    ! skip
          frad2d(j,i,24) = real(abveg(j,i))   ! skip
        end do
      end if
    end if
  end subroutine radout

end module mod_rad_outrad

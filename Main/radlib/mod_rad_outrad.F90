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
  use mod_mpmessage
  use mod_rad_common

  private

  public :: frad2d , frad3d
  public :: allocate_mod_rad_outrad , radout
  public :: nrad2d , nrad3d

  integer , parameter :: nrad2d = 24
  integer , parameter :: nrad3d = 5
  integer :: npr

  real(sp) , pointer , dimension(:,:,:) :: frad2d
  real(sp) , pointer , dimension(:,:,:,:) :: frad3d

  contains

  subroutine allocate_mod_rad_outrad
    implicit none
    npr = (jci2-jci1+1)*(ici2-ici1+1)
    call getmem3d(frad2d,jci1,jci2,ici1,ici2,1,nrad2d,'mod_outrad:frad2d')
    call getmem4d(frad3d,jci1,jci2,ici1,ici2,1,kz,1,nrad3d,'mod_outrad:frad3d')
  end subroutine allocate_mod_rad_outrad
!
  subroutine radout(lout,solin,sabtp,frsa,clrst,clrss,qrs,firtp,         &
                    frla,clrlt,clrls,qrl,slwd,sols,soll,solsd,solld,alb, &
                    albc,fsds,fsnirt,fsnrtc,fsnirtsq,totcf,totcl,totci,  &
                    h2ommr,cld,clwp,abv,sol,aeradfo,aeradfos,aerlwfo,    &
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
    logical , intent(in) :: lout ! Preapre data for outfile
    real(dp) , pointer , dimension(:) :: alb , albc , clrls , clrlt ,  &
                clrss , clrst , firtp , frla , frsa , fsds , fsnirt ,  &
                fsnirtsq , fsnrtc , sabtp , slwd , solin , soll ,      &
                solld , sols , solsd , totcf , totcl , totci , abv , sol
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
    integer :: i , j , k , n
    real(dp) :: rntim
    !
    ! total heating rate in deg/s
    !
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          heatrt(j,i,k) = qrs(n,k) + qrl(n,k)
          n = n + 1
        end do
      end do
    end do
!
!   surface absorbed solar flux in watts/m2
!   net up longwave flux at the surface
!
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        srfabswflx(j,i) = frsa(n)
        abveg(j,i) = abv(n)
        solar(j,i) = sol(n)
        srflwflxup(j,i) = frla(n)
        srflwflxdw(j,i) = slwd(n)  ! BATS Output
        n = n + 1
      end do
    end do
!
!   for coupling with bats
!
!   for now assume abveg (solar absorbed by vegetation) is equal
!   to frsa (solar absorbed by surface). possible problems are
!   over sparsely vegetated areas in which vegetation and ground
!   albedo are significantly different
!
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        totsol(j,i) = soll(n) + sols(n) + solsd(n) + solld(n)
        soldir(j,i) = sols(n)
        soldif(j,i) = solsd(n)
#ifdef CLM
        solswdir(j,i) = sols(n)
        sollwdir(j,i) = soll(n)
        solswdif(j,i) = solsd(n)
        sollwdif(j,i) = solld(n)
#endif
        n = n + 1
      end do
    end do

    if ( lchem ) then
      rntim = d_one/(d_1000*minph*chfrovrradfr)
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            aerext(j,i,k) = tauxar3d(n,k,8)
            aerssa(j,i,k) = tauasc3d(n,k,8)
            aerasp(j,i,k) = gtota3d(n,k,8)
            n = n + 1
          end do
        end do
      end do
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          aertarf(j,i)   = aertarf(j,i)   + aeradfo(n)  * rntim
          aersrrf(j,i)   = aersrrf(j,i)   + aeradfos(n) * rntim
          aertalwrf(j,i) = aertalwrf(j,i) + aerlwfo(n)  * rntim
          aersrlwrf(j,i) = aersrlwrf(j,i) + aerlwfos(n) * rntim
          n = n + 1
        end do
      end do
    end if
!
    if ( ifrad ) then
      if ( lout ) then
        do k = 1 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              frad3d(j,i,k,1) = real(cld(n,k))    ! write
              if (clwp(j,k) > dlowval) then
                frad3d(j,i,k,2) = real(clwp(n,k)) ! write
              else
                frad3d(j,i,k,2) = slowval
              end if
              frad3d(j,i,k,3) = real(qrs(n,k))    ! write
              frad3d(j,i,k,4) = real(qrl(n,k))    ! write
              frad3d(j,i,k,5) = real(h2ommr(n,k)) ! skip
              n = n + 1
            end do
          end do
        end do

        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            frad2d(j,i,1) = real(frsa(n))      ! write
            frad2d(j,i,2) = real(frla(n))      ! write
            frad2d(j,i,3) = real(clrst(n))     ! write
            frad2d(j,i,4) = real(clrss(n))     ! write
            frad2d(j,i,5) = real(clrlt(n))     ! write
            frad2d(j,i,6) = real(clrls(n))     ! write
            frad2d(j,i,7) = real(solin(n))     ! write
            frad2d(j,i,8) = real(sabtp(n))     ! write
            frad2d(j,i,9) = real(totcf(n))     ! write
            frad2d(j,i,10) = real(totcl(n))    ! write
            frad2d(j,i,11) = real(totci(n))    ! write
            frad2d(j,i,12) = real(firtp(n))    ! write
            frad2d(j,i,13) = real(alb(n))      ! skip
            frad2d(j,i,14) = real(albc(n))     ! skip
            frad2d(j,i,15) = real(fsds(n))     ! skip
            frad2d(j,i,16) = real(fsnirt(n))   ! skip
            frad2d(j,i,17) = real(fsnrtc(n))   ! skip
            frad2d(j,i,18) = real(fsnirtsq(n)) ! skip
            if ( soll(n) < dlowval ) then
              frad2d(j,i,19) = 0.0
            else
              frad2d(j,i,19) = real(soll(n))   ! skip
            end if
            if ( sols(n) < dlowval ) then
              frad2d(j,i,20) = 0.0
            else
              frad2d(j,i,20) = real(sols(n))   ! skip
            end if
            if ( solsd(n) < dlowval ) then
              frad2d(j,i,21) = 0.0
            else
              frad2d(j,i,21) = real(solsd(n))  ! skip
            end if
            if ( solld(n) < dlowval ) then
              frad2d(j,i,22) = 0.0
            else
              frad2d(j,i,22) = real(solld(n))  ! skip
            end if
            frad2d(j,i,23) = real(solar(j,i))  ! skip
            frad2d(j,i,24) = real(abveg(j,i))  ! skip
            n = n + 1
          end do
        end do
      end if
    end if
  end subroutine radout

end module mod_rad_outrad

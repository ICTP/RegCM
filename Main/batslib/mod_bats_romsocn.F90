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

module mod_bats_romsocn
!
! ROMS ocean model
!
  use mod_realkinds
  use mod_dynparam
  use mod_service
  use mod_bats_common
  use mod_bats_zengocn
  use mod_mppparam
!
  private
!
  public :: romsocndrv
  public :: allocate_mod_bats_romsocn
!
  ! missing value
  real(dp), parameter :: MISSING_R8 = 1.0d20
  ! sst
  real(dp), public, pointer, dimension(:,:) :: sst2d
#ifdef ROMSICE
  ! minimum ice depth in mm: less that this is removed
  real(dp) , parameter :: iceminh = d_10
  ! reference hgt in mm for latent heat removal from ice
  real(dp) , parameter :: href = d_two * iceminh
  ! steepness factor of latent heat removal
  real(dp) , parameter :: steepf = 1.0D0  ! Tuning needed
  real(dp), public, pointer, dimension(:,:) :: hice2d
#endif
!
  contains

  subroutine allocate_mod_bats_romsocn()
    implicit none
    call getmem2d(sst2d,jci1,jci2,ici1,ici2,'roms:sst2d') 
    sst2d  = MISSING_R8
!
#ifdef ROMSICE
    call getmem2d(hice2d,jci1,jci2,ici1,ici2,'roms:hice2d') 
    hice2d = MISSING_R8
#endif
  end subroutine allocate_mod_bats_romsocn

  subroutine romsocndrv
    implicit none
!
    real(dp) :: uv995, tsurf, t995, q995, z995, zi, psurf
    real(dp) :: qs, uv10, tau, lh, sh, dth, dqh, ustar, zo
    real(dp) :: facttq
#ifdef ROMSICE
    real(dp) :: toth 
#endif
!
    integer i , j , n
    character (len=64) :: subroutine_name='romsocndrv'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    ! feedback from ocean -> atmospheric model
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          ! update only inland water/ocean grids
          if (iveg1(n,j,i) == 14 .or. iveg1(n,j,i) == 15) then
            ! update sea surface temperature
            if (sst2d(j,i) .lt. MISSING_R8) then  
              ! update ground temperature
              tgrd(n,j,i) = sst2d(j,i)
              tgbrd(n,j,i) = sst2d(j,i)
              ! calculate bulk fluxes using zengocean bulk flux alg.
              uv995 = dsqrt(uatm(j,i,kz)**d_two+vatm(j,i,kz)**d_two)
              tsurf = sst2d(j,i)-tzero
              t995 = tatm(j,i,kz)-tzero
              q995 = qxatm(j,i,kz,iqv)/(d_one+qxatm(j,i,kz,iqv))
              z995 = hgt(j,i,kz)
              zi = hpbl(j,i)
              psurf = (sfps(j,i)+ptop)*d_10
              call zengocn(uv995,tsurf,t995,q995,z995,zi,psurf,qs, &
                           uv10,tau,lh,sh,dth,dqh,ustar,zo)
              ! update surface variables
              sent(n,j,i) = sh
              evpr(n,j,i) = lh/wlhv
              drag(n,j,i) = ustar**d_two*rho(j,i)/uv995
              facttq = dlog(z995*d_half)/dlog(z995/zo)
              u10m(n,j,i) = uatm(j,i,kz)*uv10/uv995
              v10m(n,j,i) = vatm(j,i,kz)*uv10/uv995
              t2m(n,j,i) = t995 + tzero - dth*facttq
              q2m(n,j,i) = q995 - dqh*facttq 
            end if
#ifdef ROMSICE
            ! update land use, mask and fluxes if ice is formed
            if ((hice2d(j,i) .lt. MISSING_R8) .and. (hice2d(j,i) .gt. iceminh)) then
              ldmsk1(n,j,i) = 2
              lveg(n,j,i) = 12
              sfice(n,j,i) = hice2d(j,i) 
              ! reduce sensible heat flux for ice presence
              toth = sfice(n,j,i) + sncv(n,j,i)
              if ( toth > href ) then
                sent(n,j,i) = sent(n,j,i) * (href/toth)**steepf
              end if
            else
              ldmsk1(n,j,i) = 0
              lveg(n,j,i) = iveg1(n,j,i) 
              sfice(n,j,i) = d_zero
              sncv(n,j,i) = d_zero
              snag(n,j,i) = d_zero
            end if
#endif
          end if
        end do
      end do
    end do
!
    call time_end(subroutine_name,idindx)
  end subroutine romsocndrv
!
  subroutine print_matrix_r8(inp, iskip, jskip, pet, header)
    implicit none
!
    real(dp), intent(in) :: inp(:,:)
    integer, intent(in) ::  iskip, jskip, pet
    character(len=*), intent(in) :: header
!
    integer :: i, j, imin, imax, jmin, jmax
    character(100) :: fmt_123
!
    jmin = lbound(inp, dim=1)
    jmax = ubound(inp, dim=1)
    imin = lbound(inp, dim=2)
    imax = ubound(inp, dim=2)
!
    write(6, fmt="('PET(',I2,') - ',A)") pet, trim(header)
!
    write(fmt_123, fmt="('(/, 5X, ', I3, 'I10)')") (jmax-jmin)+1
    write(6, fmt=trim(fmt_123))  (j, j=jmin, jmax, jskip)
!  
    write(fmt_123, fmt="('(I5, ', I3, 'F10.2)')") jmax
    do i=imin, imax, iskip
      write(6, fmt=trim(fmt_123)) i, (inp(j,i),j=jmin, jmax, jskip)
    end do
!
    return
  end subroutine print_matrix_r8
!
end module mod_bats_romsocn

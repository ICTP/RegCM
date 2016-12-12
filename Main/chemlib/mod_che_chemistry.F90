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

module mod_che_chemistry

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams , only : iqv , calday
  use mod_che_common
  use mod_che_indices
  use mod_che_species
  use mod_cbmz_Global
  use mod_cbmz_Parameters
  use mod_cbmz_main1
  use mod_che_molwg

  implicit none

  private

  real(rkx) :: dtchsolv
  integer(ik8) :: kchsolv

  integer(ik4) , parameter :: kmin = 2
  real(rkx) , parameter :: kb = 1.380658e-19_rkx
  real(rkx) , parameter :: mwa = 28.97_rkx

  public :: chemistry , dtchsolv , kchsolv

  contains

    subroutine chemistry(j)
      implicit none
      integer(ik4) , intent(in) :: j
      real(rkx) :: cfactor , pfact
      real(rk8) :: change
      integer(ik4) :: i , k , ic , n

      time = dtchsolv

      ! Begining of i , k loop
      ! do not solve chemistry anyway for topmost layer
      do k = kmin , kz
        do i = ici1 , ici2
          altmid   = cpb3d(j,i,k)
          ! Skip stratosphere ( ? Should we ? )
          if ( altmid < cptrop(j,i) - d_100 ) cycle
          temp     = ctb3d(j,i,k)
          zenith   = acos(czen(j,i))*raddeg
          cfactor  = crhob3d(j,i,k) * 1.e-3_rkx * navgdr
          dens     = cfactor / amd
          C_M      = altmid*10.0_rkx/(kb*temp)
          deptha   = d_zero
          depthb   = d_zero
          altabove = d_zero
          altbelow = d_zero

          if ( ichjphcld == 1 ) then
            deptha   = sum(ctaucld(j,i,1:k-1,8))
            altabove = sum(cdzq(j,i,1:k-1)*ctaucld(j,i,1:k-1,8))
            depthb   = sum(ctaucld(j,i,k+1:kz,8))
            altbelow = sum(cdzq(j,i,k+1:kz)*ctaucld(j,i,k+1:kz,8))
            if ( deptha > d_zero ) altabove = altabove / deptha
            if ( depthb > d_zero ) altbelow = altbelow / depthb
          end if

          ! call the chemistry solver
          xr(:) = d_zero
          xrin(:) = d_zero
          xrout(:) = d_zero
          ! 1 : initialise xrin with the concentrations from
          !     previous chemsolv step
!  FAB: this fix a stability bug , but the solver might slower 
!  other option is to transport all the species.
          !if ( ktau > 0 ) then
          !  do ic = 1 , totsp
          !    xrin(ic) = real(chemall(j,i,k,ic),rk8)
          !  end do
          !end if
          ! 2 : update input concentrations for transported species only
          do n = 1 , ntr
            if ( trac%indcbmz(n) > 0 ) then
               xrin(trac%indcbmz(n)) = chib3d(j,i,k,n) * cfactor / trac%mw(n)
            end if
          end do
          ! update for water vapor
          xrin(ind_H2O)  = cqxb3d(j,i,k,iqv) * cfactor / amw

          call chemmain(real(calday,rk8),real(dtchsolv,rk8))

          ! save the concentrations of all species for next chemistry step
          do ic = 1 , totsp
            if ( abs(xrout(ic)) > 1.e-20_rkx ) then
              chemall(j,i,k,ic) = real(xrout(ic),rkx)
            else
              chemall(j,i,k,ic) = d_zero
            end if
          end do
          ! Store photolysis rates for diagnostic
          do ic = 1 , nphoto
            if ( c_jval(1,ic) > 1.e-20_rkx ) then
              jphoto(j,i,k,ic) = real(c_jval(1,ic),rkx)
            else
              jphoto(j,i,k,ic) = d_zero
            end if
          end do
          !
          ! Now calculate chemical tendencies
          ! mole.cm-3.s-1  to kg.kg-1.s-1.ps (consistency with chiten unit)
          pfact = cpsb(j,i) / cfactor / dtchsolv

          do n = 1 , ntr
            if ( trac%indcbmz(n) > 0 ) then
              change = xrout(trac%indcbmz(n)) - xrin(trac%indcbmz(n))
              if ( abs(change) > 1.e-20_rkx ) then
                chemten(j,i,k,n) = real(change,rkx) * pfact * trac%mw(n)
              else
                chemten(j,i,k,n) = d_zero
              end if
            end if
          end do
        end do ! end i , k loop
      end do
    end subroutine chemistry

end module mod_che_chemistry
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

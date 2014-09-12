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
  use mod_runparams , only : iqv
  use mod_che_common
  use mod_che_indices
  use mod_che_species
  use mod_cbmz_Global
  use mod_cbmz_Parameters
  use mod_cbmz_main1
  use mod_che_molwg

  implicit none

  private

  real(rk8) , parameter :: dtchsolv = 900.0D0
  real(rk8) , parameter :: kb = 1.380658D-19
  real(rk8) , parameter :: mwa = 28.97D0

  public :: chemistry , dtchsolv

  contains

    subroutine chemistry(j,lyear,lmonth,lday )
      implicit none
      integer(ik4) , intent(in) :: j
      integer(ik4) , intent(in) :: lyear , lmonth , lday
      ! pressk    is pressure at level k
      ! presskp1  is pressure at level k+1
      ! presskm1  is pressue  at level k-1
      ! heightk   is the approximated height at level k
      ! heightkp1 is the approximated height at level k+1
      ! heightkm1 is the approximated height at level k-1
      real(rk8) :: cfactor , pfact , pss
      real(rk8) :: pressk , presskp1 , presskm1
      real(rk8) :: heightk , heightkp1 , heightkm1
      real(rk8) , parameter :: scaleH = 7.6 !km
      integer(ik4) :: i , k , kbl , kab , ic , n

      time = dtchsolv
      !! idate = (lyear-1900)*10000+lmonth*100+lday
      
      ! Begining of i , k loop
      ! do not solve chemistry for stratosphere (k == 1)
      do k = 2 , kz
        do i = ici1 , ici2
          ! care here pressure4 is considered ???
          altmid = (cpsb(j,i)*hsigma(k)+ptop)

          temp   = ctb3d(j,i,k)
          zenith = dacos(czen(j,i))*raddeg
          dens   = crhob3d(j,i,k) * 1.D-03 * navgdr / 28.97D0
          C_M    = altmid*10.0/(kb*temp)

          deptha = d_zero
          depthb = d_zero
          altabove = d_zero
          altbelow = d_zero      

          if ( ichjphcld == 1 ) then
            if ( k == kz ) then
              pressk    = (cpsb(j,i)*hsigma(k)+ptop)
              heightk   = -1.0*scaleH*log(pressk/(cpsb(j,i)+ptop))
            end if  
            if ( k > 2 ) then
              if ( ctaucld(j,i,k-1,8) > 0.0 ) then
                do kab = k-1 , 2 , -1
                  deptha    = deptha + ctaucld(j,i,kab,8)
                  pressk    = (cpsb(j,i)*hsigma(kab)+ptop)
                  presskp1  = (cpsb(j,i)*hsigma(kab+1)+ptop)
                  heightkp1 = -1.0*scaleH*log(presskp1/(cpsb(j,i)+ptop))
                  presskm1  = (cpsb(j,i)*hsigma(kab-1)+ptop)
                  heightkm1 = -1.0*scaleH*log(presskm1/(cpsb(j,i)+ptop))
                  heightk   = 0.5*(heightkp1+heightkm1)
                  altabove  = altabove + heightk*ctaucld(j,i,kab,8)
                end do           
                altabove  = altabove/deptha
              end if
            end if
            if ( k > kz ) then
              if ( ctaucld(j,i,k+1,8) > 0.0 ) then
                do kbl = kz , k+1 , -1
                  depthb    = depthb + ctaucld(j,i,kbl,8)
                  pressk    = (cpsb(j,i)*hsigma(kbl)+ptop)
                  presskp1  = (cpsb(j,i)*hsigma(kbl+1)+ptop)
                  heightkp1 = -1.0*scaleH*log(presskp1/(cpsb(j,i)+ptop))
                  presskm1  = (cpsb(j,i)*hsigma(kbl-1)+ptop)
                  heightkm1 = -1.0*scaleH*log(presskm1/(cpsb(j,i)+ptop))
                  heightk   = 0.5*(heightkp1+heightkm1)
                  altbelow  = altbelow + heightk*ctaucld(j,i,kbl,8)
                end do 
                altbelow  = altbelow/depthb
              end if
            end if
          end if

          ! call the chemistry solver
          xr(:) = d_zero
          xrin(:) = d_zero
          ! 1 : initialise xrin with the concentrations from
          !     previous chemsolv step
          do ic = 1 , totsp
            xrin(ic) = chemall(j,i,k,ic) 
          end do
          ! 2 : update input concentrations for transported species only
          cfactor =  crhob3d(j,i,k) * 1.D-03 * navgdr

          do n = 1 , ntr
            if ( trac%indcbmz(n) > 0 ) then
              xrin(trac%indcbmz(n)) = chib3d(j,i,k,n)*cfactor/trac%mw(n)
            end if
          end do 
          ! update for water vapor  
          xrin(ind_H2O)  = cqxb3d(j,i,k,iqv)*cfactor / 18.D00

          call chemmain

          ! save the concentrations of all species for next chemistry step
          do ic = 1 , totsp
            chemall(j,i,k,ic) = xrout(ic)
          end do
          ! Store photolysis rates for diagnostic
          do ic = 1 , nphoto
            jphoto(j,i,k,ic) = c_jval(1,ic)
          end do
          !
          ! Now calculate chemical tendencies       
          ! mole.cm-3.s-1  to kg.kg-1.s-1.ps (consistency with chiten unit)
          cfactor =  crhob3d(j,i,k) * 1.D-03 * navgdr
          pfact = cpsb(j,i)/cfactor/dtchsolv

          do n = 1 , ntr
            if ( trac%indcbmz(n) > 0 ) then
              chemten(j,i,k,n) = (xrout(trac%indcbmz(n)) - &
                                  xrin(trac%indcbmz(n)))*pfact* trac%mw(n)
            end if
          end do
        end do ! end i , k loop
      end do
    end subroutine chemistry

end module mod_che_chemistry

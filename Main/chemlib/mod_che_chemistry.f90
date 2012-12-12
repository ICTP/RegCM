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
  use mod_cbmz_boxvars
  use mod_cbmz_chemmech
  use mod_cbmz_chemvars
  use mod_cbmz_molwg
  use mod_cbmz_main1
  private

  real(rk8) , parameter :: dtchsolv = 900.0D0
  real(rk8) , parameter :: kb = 1.380658D-19

  public :: chemistry , dtchsolv

  contains

    subroutine chemistry(j,lyear,lmonth,lday )

      implicit none

      integer(ik4) , intent(in) :: j
      integer(ik4) , intent(in) :: lyear , lmonth , lday

      real(rk8) :: cfactor , pfact
      integer(ik4) :: i , k , kbl , kab ,ic

      time = dtchsolv
      idate = (lyear-1900)*10000+lmonth*100+lday
      c_numitr = 20
      kmax = 1
      
      ! Begining of i , k loop
      ! do not solve chemistry for stratosphere (k == 1)
      do k = 2 , kz
        do i = ici1 , ici2
          ! care here pressure4 is considered ???
          altmid(1) = (cpsb(j,i)*hsigma(k)+ptop)
          temp(1) = ctb3d(j,i,k)
          zenith = dacos(czen(j,i))*raddeg
          dens(1) = crhob3d(j,i,k) * 1.D-03 * navgdr / 28.97D0
          deptha = d_zero
          depthb = d_zero
          altabove = d_zero
          altbelow = d_zero      
          if ( k ==  1 ) then
            ! (add the half layer ctaucld, should be no cloud in this layer ) 
            deptha =  ctaucld(j,i,k,8) *d_half
            depthb =  ctaucld(j,i,k,8) *d_half 
            ! altabove = cdzq(j,i,k) / 2 !altitude or pressure ?? 
            ! altabove, altbelow are altitude above an below weighted
            ! by cloud optical depth 
            ! here altitude is taken in kpa to be consistent with altmid 
            ! WOULD not it BE BETTER TO CONSIDER ALTITUDE in M  ?
            altabove = dsigma(k)* cpsb(j,i) * d_half *  deptha  
            altbelow = dsigma(k)* cpsb(j,i) * d_half *  depthb  
            do kbl = k+1 , kz 
              depthb = depthb + ctaucld(j,i,kbl,8)
              altbelow = altbelow + dsigma(kbl)*cpsb(j,i)*ctaucld(j,i,kbl,8)
            end do
            ! FAB TEST consider the visible taucld
          else if ( k == kz ) then 
            depthb =  ctaucld(j,i,k,8) *d_half
            deptha =  ctaucld(j,i,k,8) *d_half
            altabove = dsigma(k)* cpsb(j,i) * d_half *  deptha  
            altbelow = dsigma(k)* cpsb(j,i) * d_half *  depthb  
            do kab = 1 , k-1 
              deptha = deptha + ctaucld(j,i,kab,8)
              altabove = altabove + dsigma(kab)*cpsb(j,i)*ctaucld(j,i,kab,8)
            end do
          else
            depthb =  ctaucld(j,i,k,8) *d_half
            deptha =  ctaucld(j,i,k,8) *d_half
            altabove = dsigma(k)* cpsb(j,i) * d_half *  deptha  
            altbelow = dsigma(k)* cpsb(j,i) * d_half *  depthb  
            do kbl = k+1 , kz 
              depthb = depthb + ctaucld(j,i,kbl,8)
              altbelow = altbelow + dsigma(kbl)*cpsb(j,i)*ctaucld(j,i,kbl,8)
            end do
            do kab = 1 , k-1 
              deptha = deptha + ctaucld(j,i,kab,8)
              altabove = altabove + dsigma(kab)*cpsb(j,i)*ctaucld(j,i,kab,8)
            end do
          endif
          ! normalise the weighted altitude above and bleow cloud        
          if ( depthb > d_zero ) altbelow = altbelow / depthb 
          if ( deptha > d_zero ) altabove = altabove / deptha      

!         call the chemistry solver         
          xr(1,:) = d_zero
          xrin(1,:) = d_zero
!         1 : initialise xrin with the concentrations from 
!         previous chemsolv step
          do ic = 1 , totsp
            xrin(1,ic) = chemall(j,i,k,ic) 
          end do
!         2  : update input concentrations for transported species only  
          cfactor =  crhob3d(j,i,k) * 1.D-03 * navgdr
          xrin(1,ind_H2O)  = cqxb3d(j,i,k,iqv)*cfactor / 18.D00
          xrin(1,ind_O3)   = chib3d(j,i,k,io3)*cfactor/W_O3
          xrin(1,ind_NO2)  = chib3d(j,i,k,ino2)*cfactor /W_NO2
          xrin(1,ind_NO)   = chib3d(j,i,k,ino)*cfactor/W_NO
          xrin(1,ind_CO)   = chib3d(j,i,k,ico)*cfactor/W_CO
          xrin(1,ind_H2O2) = chib3d(j,i,k,ih2o2)*cfactor/W_H2O2
          xrin(1,ind_HNO3) = chib3d(j,i,k,ihno3)*cfactor/W_HNO3
          xrin(1,ind_N2O5) = chib3d(j,i,k,in2o5)*cfactor/W_N2O5
          xrin(1,ind_SO2)  = chib3d(j,i,k,iso2)*cfactor/W_SO2
          xrin(1,ind_SULF) = chib3d(j,i,k,iso4)*cfactor/W_SULF
          xrin(1,ind_HCHO) = chib3d(j,i,k,ihcho)*cfactor/W_HCHO
          xrin(1,ind_ALD2) = chib3d(j,i,k,iald2)*cfactor/W_ALD2
          xrin(1,ind_ISOP) = chib3d(j,i,k,iisop)*cfactor/W_ISOP
          xrin(1,ind_C2H6) = chib3d(j,i,k,ic2h6)*cfactor/W_C2H6
          xrin(1,ind_PAR)  = chib3d(j,i,k,ipar)*cfactor/W_C3H8
          xrin(1,ind_ETHE) = chib3d(j,i,k,iethe)*cfactor/W_ETHENE
          xrin(1,ind_PRPE) = chib3d(j,i,k,iolt)*cfactor/W_OLT
          xrin(1,ind_BUTE) = chib3d(j,i,k,ioli)*cfactor/W_OLI
          xrin(1,ind_TOLU) = chib3d(j,i,k,itolue)*cfactor/W_TOLU
          xrin(1,ind_XYLE) = chib3d(j,i,k,ixyl)*cfactor/W_XYLE
          xrin(1,ind_PAN)  = chib3d(j,i,k,ipan)*cfactor/W_PAN
          xrin(1,ind_CH4)  = chib3d(j,i,k,ich4)*cfactor/W_CH4
          xrin(1,ind_MOH)  = chib3d(j,i,k,imoh)*cfactor/W_MOH
          xrin(1,ind_ACET) = chib3d(j,i,k,iacet)*cfactor/W_ACET
          xrin(1,ind_RCOOH) = chib3d(j,i,k,ircooh)*cfactor/W_RCOOH

          xrin(1,ind_mgly) = chib3d(j,i,k,imgly)*cfactor/w_MGLY
          xrin(1,ind_cres) = chib3d(j,i,k,icres)*cfactor/W_CRES
          xrin(1,ind_open) = chib3d(j,i,k,iopen)*cfactor/W_OPEN
          xrin(1,ind_isoprd) = chib3d(j,i,k,iisoprd)*cfactor/W_ISOPRD
          xrin(1,ind_onit) = chib3d(j,i,k,ionit)*cfactor/W_ONIT
          xrin(1,ind_hcooh ) = chib3d(j,i,k,ihcooh)*cfactor/W_HCOOH
          xrin(1,ind_ch3ooh) = chib3d(j,i,k,ich3ooh)*cfactor/W_CH3OOH
          xrin(1,ind_ethooh) = chib3d(j,i,k,iethooh)*cfactor/W_ETHOOH
          xrin(1,ind_rooh) = chib3d(j,i,k,irooh)*cfactor/W_ROOH
          xrin(1,ind_hono) = chib3d(j,i,k,ihono)*cfactor/W_HONO
          xrin(1,ind_hno4) = chib3d(j,i,k,ihno4)*cfactor/W_HNO4
          xrin(1,ind_xo2 ) = chib3d(j,i,k,ixo2)*cfactor/W_XO2

!         solver work with xr 
          xr(:,:) = xrin(:,:)

          call chemmain

          xrout(:,:) = xr(:,:)
!         save the concentrations of all species for next chemistry step
          do ic = 1 , totsp
            chemall(j,i,k,ic) = xrout(1,ic)
          end do
          ! Store photolysis rates for diagnostic
          do ic = 1 , nphoto
            jphoto(j,i,k,ic) = c_jval(1,ic)
          end do
          !
          ! Now calculate chemical tendencies       
          !
          ! mole.cm-3.s-1  to kg.kg-1.s-1.ps (consistency with chiten unit)
          cfactor =  crhob3d(j,i,k) * 1.D-03 * navgdr
          pfact = cpsb(j,i)/cfactor/dtchsolv
          chemten(j,i,k,io3)    = &
            (xrout(1,ind_O3)   - xrin(1,ind_O3))  *pfact*W_O3
          chemten(j,i,k,ino2)   = &
            (xrout(1,ind_NO2)  - xrin(1,ind_NO2)) *pfact*W_NO2
          chemten(j,i,k,ino)    = &
            (xrout(1,ind_NO)   - xrin(1,ind_NO))  *pfact*W_NO
          chemten(j,i,k,ico)    = &
            (xrout(1,ind_CO)   - xrin(1,ind_CO))  *pfact*W_CO
          chemten(j,i,k,ih2o2)  = &
            (xrout(1,ind_H2O2) - xrin(1,ind_H2O2))*pfact*W_H2O2
          chemten(j,i,k,ihno3)  = &
            (xrout(1,ind_HNO3) - xrin(1,ind_HNO3))*pfact*W_HNO3
          chemten(j,i,k,in2o5)  = &
            (xrout(1,ind_N2O5) - xrin(1,ind_N2O5))*pfact*W_N2O5
          chemten(j,i,k,iso2)   = &
            (xrout(1,ind_SO2)  - xrin(1,ind_SO2)) *pfact*W_SO2
          chemten(j,i,k,iso4)   = &
            (xrout(1,ind_SULF) - xrin(1,ind_SULF))*pfact*W_SULF
          chemten(j,i,k,ihcho)  = &
            (xrout(1,ind_HCHO) - xrin(1,ind_HCHO))*pfact*W_HCHO
          chemten(j,i,k,iald2)  = &
            (xrout(1,ind_ALD2) - xrin(1,ind_ALD2))*pfact*W_ALD2
          chemten(j,i,k,iisop)  = &
            (xrout(1,ind_ISOP) - xrin(1,ind_ISOP))*pfact*W_ISOP
          chemten(j,i,k,ic2h6)  = &
            (xrout(1,ind_C2H6) - xrin(1,ind_C2H6))*pfact*W_C2H6
          chemten(j,i,k,ipar)   = &
            (xrout(1,ind_PAR)  - xrin(1,ind_PAR)) *pfact*W_C3H8
          chemten(j,i,k,itolue) = &
            (xrout(1,ind_TOLU) - xrin(1,ind_TOLU))*pfact*W_TOLU
          chemten(j,i,k,ixyl)   = &
            (xrout(1,ind_XYLE) - xrin(1,ind_XYLE))*pfact*W_XYLE
          chemten(j,i,k,iethe)  = &
            (xrout(1,ind_ETHE) - xrin(1,ind_ETHE))*pfact*W_ETHENE
          chemten(j,i,k,ipan)   = &
            (xrout(1,ind_PAN)  - xrin(1,ind_PAN)) *pfact*W_PAN
          chemten(j,i,k,ich4)   = &
            (xrout(1,ind_CH4)  - xrin(1,ind_CH4)) *pfact*W_CH4
          chemten(j,i,k,iolt)   = &
            (xrout(1,ind_PRPE) - xrin(1,ind_PRPE))*pfact*W_OLT
          chemten(j,i,k,ioli)   = &
            (xrout(1,ind_BUTE) - xrin(1,ind_BUTE))*pfact*W_OLI
          chemten(j,i,k,imoh)   = &
            (xrout(1,ind_MOH)  - xrin(1,ind_MOH)) *pfact*W_MOH
          chemten(j,i,k,iacet)  = &
            (xrout(1,ind_ACET) - xrin(1,ind_ACET))*pfact*W_ACET
          chemten(j,i,k,ircooh)  = &
            (xrout(1,ind_RCOOH) - xrin(1,ind_RCOOH))*pfact*W_RCOOH
          chemten(j,i,k,imgly)  = &
            (xrout(1,ind_mgly) - xrin(1,ind_mgly))*pfact*W_mgly
           chemten(j,i,k,icres)  = &
            (xrout(1,ind_cres) - xrin(1,ind_cres))*pfact*W_cres
           chemten(j,i,k,iopen)  = &
            (xrout(1,ind_open) - xrin(1,ind_open))*pfact*W_open
           chemten(j,i,k,iisoprd)  = &
            (xrout(1,ind_isoprd) - xrin(1,ind_isoprd))*pfact*W_isoprd
           chemten(j,i,k,ionit)  = &
            (xrout(1,ind_onit) - xrin(1,ind_onit))*pfact*W_onit
           chemten(j,i,k,ihcooh)  = &
            (xrout(1,ind_hcooh) - xrin(1,ind_hcooh))*pfact*W_hcooh
           chemten(j,i,k,ich3ooh)  = &
            (xrout(1,ind_ch3ooh) - xrin(1,ind_ch3ooh))*pfact*W_ch3ooh
           chemten(j,i,k,iethooh)  = &
            (xrout(1,ind_ethooh) - xrin(1,ind_ethooh))*pfact*W_ethooh
           chemten(j,i,k,irooh)  = &
            (xrout(1,ind_rooh) - xrin(1,ind_rooh))*pfact*W_rooh
           chemten(j,i,k,ihono)  = &
            (xrout(1,ind_hono) - xrin(1,ind_hono))*pfact*W_hono
           chemten(j,i,k,ihno4)  = &
            (xrout(1,ind_hno4) - xrin(1,ind_hno4))*pfact*W_hno4
           chemten(j,i,k,ixo2)  = &
            (xrout(1,ind_xo2) - xrin(1,ind_xo2))*pfact*W_xo2


        end do ! end i , k loop
      end do

      if (ichdiag > 0 )chemdiag(j,:,:,:) = chemdiag(j,:,:,:) + &
              chemten(j,:,:,:) * dble(dtchsolv) / (3600D0 * dble(chemfrq)) 



    end subroutine chemistry
!
end module mod_che_chemistry

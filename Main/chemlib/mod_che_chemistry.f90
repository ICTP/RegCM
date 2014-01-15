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

!pressk    is pressure at level k
!presskp1  is pressure at level k+1
!presskm1  is pressue  at level k-1
!heightk   is the approximated height at level k
!heightkp1   is the approximated height at level k+1
!heightkm1   is the approximated height at level k-1
      real(rk8) :: cfactor , pfact,pss
      real(rk8) :: pressk,presskp1,presskm1
      real(rk8) :: heightk,heightkp1,heightkm1
      real(rk8), parameter  ::scaleH=7.6 !km
      integer(ik4) :: i , k , kbl , kab ,ic

      time = dtchsolv
!!     idate = (lyear-1900)*10000+lmonth*100+lday
      
      ! Begining of i , k loop
      ! do not solve chemistry for stratosphere (k == 1)
      do k = 2 , kz
        do i = ici1 , ici2
          ! care here pressure4 is considered ???
          altmid = (cpsb(j,i)*hsigma(k)+ptop)

          if (k .eq. kz) then
          pressk    = (cpsb(j,i)*hsigma(k)+ptop)
          heightk   = -1.0*scaleH*log(pressk/(cpsb(j,i)+ptop))
          end if 
!          if (k .lt. kz .and. k .ge. 2) then
!          pressk    = (cpsb(j,i)*hsigma(k)+ptop)
!          heightk   = -1.0*scaleH*log(pressk/(cpsb(j,i)+ptop))
!          presskp1  = (cpsb(j,i)*hsigma(k+1)+ptop)
!          heightkp1   = -1.0*scaleH*log(presskp1/(cpsb(j,i)+ptop))
!          presskm1    = (cpsb(j,i)*hsigma(k-1)+ptop)
!          heightkm1   = -1.0*scaleH*log(presskm1/(cpsb(j,i)+ptop))
!          end if

          temp      = ctb3d(j,i,k)
          zenith    = dacos(czen(j,i))*raddeg
          dens   = crhob3d(j,i,k) * 1.D-03 * navgdr / 28.97D0
          C_M       = altmid*10.0/(kb*temp)

          deptha = d_zero
          depthb = d_zero
          altabove = d_zero
          altbelow = d_zero      
!!          if ( k ==  2 ) then
            ! (add the half layer ctaucld, should be no cloud in this layer ) 
!!            deptha =  ctaucld(j,i,k,8) *d_half
!!            depthb =  ctaucld(j,i,k,8) *d_half 
            ! altabove = cdzq(j,i,k) / 2 !altitude or pressure ?? 
            ! altabove, altbelow are altitude above an below weighted
            ! by cloud optical depth 
            ! here altitude is taken in kpa to be consistent with altmid 
            ! WOULD not it BE BETTER TO CONSIDER ALTITUDE in M  ?
!!            altabove = dsigma(k)* cpsb(j,i) * d_half *  deptha  
!!            altbelow = dsigma(k)* cpsb(j,i) * d_half *  depthb  
!!            do kbl = k+1 , kz 
!!              depthb = depthb + ctaucld(j,i,kbl,8)
!!              altbelow = altbelow + dsigma(kbl)*cpsb(j,i)*ctaucld(j,i,kbl,8)
!!            end do
            ! FAB TEST consider the visible taucld
!!          else if ( k == kz ) then 
!!            depthb =  ctaucld(j,i,k,8) *d_half
!!            deptha =  ctaucld(j,i,k,8) *d_half
!!            altabove = dsigma(k)* cpsb(j,i) * d_half *  deptha  
!!            altbelow = dsigma(k)* cpsb(j,i) * d_half *  depthb  
            if (k .gt. 2 ) then
               if (ctaucld(j,i,k-1,8) .gt. 0.0) then
            do kab = k-1, 2,-1
              deptha    = deptha + ctaucld(j,i,kab,8)
              pressk    = (cpsb(j,i)*hsigma(kab)+ptop)
              presskp1  = (cpsb(j,i)*hsigma(kab+1)+ptop)
              heightkp1 = -1.0*scaleH*log(presskp1/(cpsb(j,i)+ptop))
              presskm1  = (cpsb(j,i)*hsigma(kab-1)+ptop)
              heightkm1 = -1.0*scaleH*log(presskm1/(cpsb(j,i)+ptop))
              heightk   = 0.5*(heightkp1+heightkm1)
              altabove  = altabove + heightk*ctaucld(j,i,kab,8)
            end do
            deptha    = deptha 
            altabove  = altabove/deptha
               end if
            end if
            if (k .lt. kz ) then
               if (ctaucld(j,i,k+1,8) .gt. 0.0) then
            do kbl = kz,k+1,-1
            depthb    = depthb + ctaucld(j,i,kbl,8)
            pressk    = (cpsb(j,i)*hsigma(kbl)+ptop)
            presskp1  = (cpsb(j,i)*hsigma(kbl+1)+ptop)
            heightkp1 = -1.0*scaleH*log(presskp1/(cpsb(j,i)+ptop))
            presskm1  = (cpsb(j,i)*hsigma(kbl-1)+ptop)
            heightkm1 = -1.0*scaleH*log(presskm1/(cpsb(j,i)+ptop))
            heightk   = 0.5*(heightkp1+heightkm1)
            altbelow  = altbelow + heightk*ctaucld(j,i,kbl,8)
            end do 
            depthb    = depthb
            altbelow  = altbelow/depthb
              end if
            end if
!!          else
!!            depthb =  ctaucld(j,i,k,8) *d_half
!!            deptha =  ctaucld(j,i,k,8) *d_half
!!            altabove = dsigma(k)* cpsb(j,i) * d_half *  deptha  
!!            altbelow = dsigma(k)* cpsb(j,i) * d_half *  depthb  
!!            do kbl = k+1 , kz 
!!              depthb = depthb + ctaucld(j,i,kbl,8)
!!              altbelow = altbelow + dsigma(kbl)*cpsb(j,i)*ctaucld(j,i,kbl,8)
!!            end do
!!            do kab = 1 , k-1 
!!              deptha = deptha + ctaucld(j,i,kab,8)
!!              altabove = altabove + dsigma(kab)*cpsb(j,i)*ctaucld(j,i,kab,8)
!!            end do
!!          endif
          ! normalise the weighted altitude above and bleow cloud        
!!          if ( depthb > d_zero ) altbelow = altbelow / depthb 
!!          if ( deptha > d_zero ) altabove = altabove / deptha      
!       if(myid==3 .and. j==35 .and.i==10)then
!       if(i==10)then

!       write(*,*)k,depthb,altbelow,ctaucld(j,i,k,8)
!       end if

!         call the chemistry solver         
          xr(:) = d_zero
          xrin(:) = d_zero
!         1 : initialise xrin with the concentrations from 
!         previous chemsolv step
          do ic = 1 , totsp
!            xrin(ic) = chemall(j,i,k,ic) 
          end do
!         2  : update input concentrations for transported species only  
          cfactor =  crhob3d(j,i,k) * 1.D-03 * navgdr
          xrin(ind_H2O)  = cqxb3d(j,i,k,iqv)*cfactor / 18.D00

          xrin(ind_NO)   = chib3d(j,i,k,ino)*cfactor/W_NO
          xrin(ind_NO2)  = chib3d(j,i,k,ino2)*cfactor /W_NO2
          xrin(ind_N2O5) = chib3d(j,i,k,in2o5)*cfactor/W_N2O5
          xrin(ind_HNO2) = chib3d(j,i,k,ihono)*cfactor/W_HONO
          xrin(ind_HNO3) = chib3d(j,i,k,ihno3)*cfactor/W_HNO3
          xrin(ind_HNO4) = chib3d(j,i,k,ihno4)*cfactor/W_HNO4
          xrin(ind_O3)   = chib3d(j,i,k,io3)*cfactor/W_O3
          xrin(ind_OH)   = chib3d(j,i,k,ioh)*cfactor/W_OH
          xrin(ind_HO2)  = chib3d(j,i,k,iho2)*cfactor/W_HO2
          xrin(ind_H2O2) = chib3d(j,i,k,ih2o2)*cfactor/W_H2O2
          xrin(ind_CO)   = chib3d(j,i,k,ico)*cfactor/W_CO
          xrin(ind_SO2)  = chib3d(j,i,k,iso2)*cfactor/W_SO2
          xrin(ind_H2SO4)= chib3d(j,i,k,ih2so4)*cfactor/W_SULF

          xrin(ind_CH4)  = chib3d(j,i,k,ich4)*cfactor/W_CH4
          xrin(ind_C2H6) = chib3d(j,i,k,ic2h6)*cfactor/W_C2H6
          xrin(ind_PAR)  = chib3d(j,i,k,ipar)*cfactor/W_C3H8
          xrin(ind_CH3OH)= chib3d(j,i,k,imoh)*cfactor/W_MOH
          xrin(ind_HCHO) = chib3d(j,i,k,ihcho)*cfactor/W_HCHO
          xrin(ind_ALD2) = chib3d(j,i,k,iald2)*cfactor/W_ALD2
          xrin(ind_AONE) = chib3d(j,i,k,iacet)*cfactor/W_AONE

          xrin(ind_ETH) = chib3d(j,i,k,iethe)*cfactor/W_ETHENE
          xrin(ind_OLET) = chib3d(j,i,k,iolt)*cfactor/W_OLT
          xrin(ind_OLEI) = chib3d(j,i,k,ioli)*cfactor/W_OLI

          xrin(ind_TOL) = chib3d(j,i,k,itolue)*cfactor/W_TOLU
          xrin(ind_XYL) = chib3d(j,i,k,ixyl)*cfactor/W_XYLE

          xrin(ind_ISOP) = chib3d(j,i,k,iisop)*cfactor/W_ISOP
          xrin(ind_ONIT) = chib3d(j,i,k,ionit)*cfactor/W_ONIT
          xrin(ind_PAN)  = chib3d(j,i,k,ipan)*cfactor/W_PAN

          xrin(ind_HCOOH) = chib3d(j,i,k,ihcooh)*cfactor/W_HCOOH
          xrin(ind_RCOOH) = chib3d(j,i,k,ircooh)*cfactor/W_RCOOH
          xrin(ind_CH3OOH)= chib3d(j,i,k,ich3ooh)*cfactor/W_CH3OOH
          xrin(ind_ETHOOH)= chib3d(j,i,k,iethooh)*cfactor/W_ETHOOH
          xrin(ind_ROOH)  = chib3d(j,i,k,irooh)*cfactor/W_ROOH

          xrin(ind_RO2)  = chib3d(j,i,k,iro2)*cfactor/W_RO2
          xrin(ind_XO2 ) = chib3d(j,i,k,ixo2)*cfactor/W_XO2
          xrin(ind_DMS ) = chib3d(j,i,k,idms)*cfactor/W_DMS


!         solver work with xr 
!          xr(:) = xrin(:)
!we do not need xr(:,:) any more

          call chemmain

!          xrout(:) = xr(:)
!         save the concentrations of all species for next chemistry step
          do ic = 1 , totsp
!            chemall(j,i,k,ic) = xrout(ic)
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

          chemten(j,i,k,ino)    = (xrout(ind_NO)   - xrin(ind_NO))  *pfact*W_NO
          chemten(j,i,k,ino2)   = (xrout(ind_NO2)  - xrin(ind_NO2)) *pfact*W_NO2
          chemten(j,i,k,in2o5)  = (xrout(ind_N2O5) - xrin(ind_N2O5))*pfact*W_N2O5
          chemten(j,i,k,ihono)  = (xrout(ind_HNO2) - xrin(ind_HNO2))*pfact*W_hono
          chemten(j,i,k,ihno3)  = (xrout(ind_HNO3) - xrin(ind_HNO3))*pfact*W_HNO3
          chemten(j,i,k,ihno4)  = (xrout(ind_hno4) - xrin(ind_hno4))*pfact*W_hno4



          chemten(j,i,k,io3)    = (xrout(ind_O3)   - xrin(ind_O3))   *pfact*W_O3
          chemten(j,i,k,ioh)    = (xrout(ind_OH)   - xrin(ind_OH))   *pfact*W_OH
          chemten(j,i,k,iho2)   = (xrout(ind_HO2)  - xrin(ind_HO2))  *pfact*W_HO2
          chemten(j,i,k,ih2o2)  = (xrout(ind_H2O2) - xrin(ind_H2O2)) *pfact*W_H2O2
          chemten(j,i,k,ico)    = (xrout(ind_CO)   - xrin(ind_CO))   *pfact*W_CO
          chemten(j,i,k,iso2)   = (xrout(ind_SO2)  - xrin(ind_SO2))  *pfact*W_SO2
          chemten(j,i,k,ih2so4) = (xrout(ind_H2SO4)- xrin(ind_H2SO4))*pfact*W_H2SO4


          chemten(j,i,k,ich4)   = (xrout(ind_CH4)  - xrin(ind_CH4))  *pfact*W_CH4
          chemten(j,i,k,ic2h6)  = (xrout(ind_C2H6) - xrin(ind_C2H6)) *pfact*W_C2H6
          chemten(j,i,k,ipar)   = (xrout(ind_PAR)  - xrin(ind_PAR))  *pfact*W_PAR
          chemten(j,i,k,imoh)   = (xrout(ind_CH3OH)- xrin(ind_CH3OH))*pfact*W_MOH
          chemten(j,i,k,ihcho)  = (xrout(ind_HCHO) - xrin(ind_HCHO)) *pfact*W_HCHO
          chemten(j,i,k,iald2)  = (xrout(ind_ALD2) - xrin(ind_ALD2)) *pfact*W_ALD2

          chemten(j,i,k,iacet)  = (xrout(ind_AONE) - xrin(ind_AONE)) *pfact*W_AONE
          chemten(j,i,k,iethe)  = (xrout(ind_ETH)  - xrin(ind_ETH))  *pfact*W_ETHENE
          chemten(j,i,k,iolt)   = (xrout(ind_OLET) - xrin(ind_OLET)) *pfact*W_OLT
          chemten(j,i,k,ioli)   = (xrout(ind_OLEI) - xrin(ind_OLEI)) *pfact*W_OLI

          chemten(j,i,k,itolue) = (xrout(ind_TOL) - xrin(ind_TOL))*pfact*W_TOLU
          chemten(j,i,k,ixyl)   = (xrout(ind_XYL) - xrin(ind_XYL))*pfact*W_XYLE


          chemten(j,i,k,iisop)  = (xrout(ind_ISOP) - xrin(ind_ISOP))*pfact*W_ISOP
          chemten(j,i,k,ionit)  = (xrout(ind_ONIT) - xrin(ind_ONIT))*pfact*W_ONIT
          chemten(j,i,k,ipan)   = (xrout(ind_PAN)  - xrin(ind_PAN)) *pfact*W_PAN


          chemten(j,i,k,ihcooh) = (xrout(ind_HCOOH) - xrin(ind_HCOOH)) *pfact*W_HCOOH
          chemten(j,i,k,ircooh) = (xrout(ind_RCOOH) - xrin(ind_RCOOH)) *pfact*W_RCOOH
          chemten(j,i,k,ich3ooh)= (xrout(ind_CH3OOH)- xrin(ind_CH3OOH))*pfact*W_CH3OOH
          chemten(j,i,k,iethooh)= (xrout(ind_ETHOOH)- xrin(ind_ETHOOH))*pfact*W_ETHOOH
          chemten(j,i,k,irooh)  = (xrout(ind_ROOH)  - xrin(ind_ROOH))  *pfact*W_ROOH



          chemten(j,i,k,iro2)   = (xrout(ind_RO2) - xrin(ind_RO2))*pfact*W_RO2
          chemten(j,i,k,ixo2)   = (xrout(ind_XO2) - xrin(ind_XO2))*pfact*W_XO2
          chemten(j,i,k,idms)   = (xrout(ind_DMS) - xrin(ind_DMS))*pfact*W_DMS


           pss=c_jval(1,4)*xrout(ind_no2)/(1.9e-14*xrout(ind_no)*xrout(ind_o3))
       if(myid==3 .and. j==35 .and. k == kz .and. i==10)then
!       write(*,*)tochar(idatex),xrout(1,ind_NO2)/2.5e10,xrout(1,ind_O3)/2.5e10,xrout(1,ind_OH),xrout(1,ind_HO2)
!       write(*,*)xrout(1,ind_NO2)/2.5e10,xrout(1,ind_O3)/2.5e10,xrout(1,ind_OH),xrout(1,ind_HO2)
!       write(*,555)ktau,xrin(1,ind_NO2)/2.5e10,xrout(1,ind_O3)/2.5e10,xrout(1,ind_OH),xrout(1,ind_HO2)
!       write(*,*)xrin(1,ind_NO2)/2.5e10,c_xcin(1,ind_NO2)/2.5e10,xrin(1,ind_O3)/2.5e10,c_xcin(1,ind_O3)/2.5e10
!       write(*,*)xrout(1,ind_hono)/2.5e10,c_xcout(1,ind_hono)/2.5e10,xrout(1,ind_co)/2.5e10,c_xcout(1,ind_co)/2.5e10
       write(*,*)xrout(ind_oh)/2.4e10,xrout(ind_no2)/2.4e10,xrout(ind_xo2)/2.4e10,pss,c_jval(1,4)
       end if

        end do ! end i , k loop
      end do
!      write(*,*)maxval(ctaucld(j,:,:,:)),maxval(ctaucld(j,:,:,8))
!      write(*,*)ici1,ici2,jci1,jci2

    end subroutine chemistry
!
end module mod_che_chemistry

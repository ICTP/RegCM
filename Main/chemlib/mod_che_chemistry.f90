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
  
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_che_common
  use mod_che_indices
  use mod_che_species
  use mod_cbmz_boxvars
  use mod_cbmz_chemmech
  use mod_cbmz_chemvars
  use mod_cbmz_molwg
  use mod_cbmz_main1
  private

  real(dp) , parameter :: dtchsolv = 900.0D0
! 

  public :: chemistry , dtchsolv

  integer , parameter :: jvO2 = 1
  integer , parameter :: jvO3a = 2
  integer , parameter :: jvO3b = 3
  integer , parameter :: jvNO2 = 4
  integer , parameter :: jvNO3a = 5
  integer , parameter :: jvNO3b = 6
  integer , parameter :: jvN2O5a = 7
  integer , parameter :: jvN2O5b = 8
  integer , parameter :: jvN2O = 9
  integer , parameter :: jvHO2 = 10
  integer , parameter :: jvH2O2 = 11
  integer , parameter :: jvHNO2 = 12
  integer , parameter :: jvHNO3 = 13
  integer , parameter :: jvHNO4 = 14
  integer , parameter :: jvCH2Oa = 15
  integer , parameter :: jvCH2Ob = 16
  integer , parameter :: jvCH3CHOa = 17
  integer , parameter :: jvCH3CHOb = 18
  integer , parameter :: jvCH3CHOc = 19
  integer , parameter :: jvC2H5CHO = 20
  integer , parameter :: jvCHOCHO = 21
  integer , parameter :: jvCH3COCHO = 22
  integer , parameter :: jvCH3COCH3 = 23
  integer , parameter :: jvCH3OOH = 24
  integer , parameter :: jvCH3ONO2 = 25
  integer , parameter :: jvPAN = 26

  real(dp) , parameter :: kb = 1.380658D-19

  contains

    subroutine chemistry(j,secofday,lyear,lmonth,lday )

      implicit none

      integer , intent(in) :: j
      integer , intent(in) :: lyear , lmonth , lday
      real(dp) , intent(in) :: secofday 
      real(dp) , dimension(ici1:ici2,1:kz,1:56) :: jphoto
      real(dp) :: cfactor , pfact
      integer :: i , k , kbl , kab ,ic

      time = dtchsolv
      idate = (lyear-1900)*10000+lmonth*100+lday
      xhour =   secofday/3600.0D0   ! abt added for time of day
      c_numitr = 20
      kmax = 1
      
      ! initialize jphoto to zero
      ! initialize jphoto to zero
      jphoto(:,:,:) = d_zero

      ! Begining of i , k loop
      ! do not solve chemistry for stratosphere (k == 1)
      do k = 2 , kz
        do i = ici1 , ici2
          ! care here pressure4 is considered ???
          altmid(1) = (cpsb(j,i)*hlev(k)+chptop)
          temp(1) = ctb3d(j,i,k)
          zenith = dacos(czen(j,i))*raddeg
          dens(1) = crhob3d(j,i,k) * 1.D-03 * navgdr / 28.97D0
          deptha = d_zero
          depthb = d_zero
          altabove= d_zero
          altbelow= d_zero
          if ( 1 == 2 ) then             
            if ( k ==  1 ) then
              ! (add the half layer ctaucld, should be no cloud in this layer ) 
              deptha =  ctaucld(j,i,k,8) *d_half
              depthb =  ctaucld(j,i,k,8) *d_half 
!             altabove = cdzq(j,i,k) / 2 !altitude or pressure ?? 
              ! altabove, altbelow are altitude above an below weighted
              ! by cloud optical depth 
              ! here altitude is taken in kpa to be consistent with altmid 
              ! WOULD not it BE BETTER TO CONSIDER ALTITUDE in M  ?
              altabove = cdsigma(k)* cpsb(j,i) * d_half *  deptha  
              altbelow = cdsigma(k)* cpsb(j,i) * d_half *  depthb  
              do kbl = k+1 , kz 
                depthb = depthb + ctaucld(j,i,kbl,8)
                altbelow = altbelow + cdsigma(kbl)*cpsb(j,i)*ctaucld(j,i,kbl,8)
              end do
              ! FAB TEST consider the visible taucld
            else if ( k == kz ) then 
              depthb =  ctaucld(j,i,k,8) *d_half
              deptha =  ctaucld(j,i,k,8) *d_half
              altabove = cdsigma(k)* cpsb(j,i) * d_half *  deptha  
              altbelow = cdsigma(k)* cpsb(j,i) * d_half *  depthb  
              do kab = 1 , k-1 
                deptha = deptha + ctaucld(j,i,kab,8)
                altabove = altabove + cdsigma(kab)*cpsb(j,i)*ctaucld(j,i,kab,8)
              end do
            else
              depthb =  ctaucld(j,i,k,8) *d_half
              deptha =  ctaucld(j,i,k,8) *d_half
              altabove = cdsigma(k)* cpsb(j,i) * d_half *  deptha  
              altbelow = cdsigma(k)* cpsb(j,i) * d_half *  depthb  
              do kbl = k+1 , kz 
                depthb = depthb + ctaucld(j,i,kbl,8)
                altbelow = altbelow + cdsigma(kbl)*cpsb(j,i)*ctaucld(j,i,kbl,8)
              end do
              do kab = 1 , k-1 
                deptha = deptha + ctaucld(j,i,kab,8)
                altabove = altabove + cdsigma(kab)*cpsb(j,i)*ctaucld(j,i,kab,8)
              end do
            endif
            ! normalise the weighted altitude above and bleow cloud        
            if (depthb > d_zero)  altbelow = altbelow / depthb 
            if (deptha >d_zero )  altabove = altabove / deptha      
          end if
          do ic = 1 , totsp
            xr(1,ic) = d_zero
          end do
          do ic = 1 , totsp
            xr(1,ic) = chemall(j,i,k,ic) 
          end do
          cfactor =  crhob3d(j,i,k) * 1.D-03 * navgdr
          xrin(1,ind_H2O)  = cqvb3d(j,i,k)*cfactor / 18.D00
          xrin(1,ind_O3)   = chib3d(j,i,k,io3)*cfactor/W_O3
          xrin(1,ind_NO2)  = chib3d(j,i,k,ino2)*cfactor /W_NO2
          xrin(1,ind_NO)   = chib3d(j,i,k,ino)*cfactor/W_NO
          xrin(1,ind_CO)   = chib3d(j,i,k,ico)*cfactor/W_CO
          xrin(1,ind_H2O2) = chib3d(j,i,k,ih2o2)*cfactor/W_H2O2
          xrin(1,ind_HNO3) = chib3d(j,i,k,ihno3)*cfactor/W_HNO3
          xrin(1,ind_N2O5) = chib3d(j,i,k,in2o5)*cfactor/W_N2O5
          xrin(1,ind_SO2)  = chib3d(j,i,k,iso2)*cfactor/W_SO2
          xrin(1,ind_SULF) = chib3d(j,i,k,iso4)*cfactor/W_SULF
          xrin(1,ind_DMS)  = chib3d(j,i,k,idms)*cfactor/W_DMS
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
          xrin(1,ind_rcooh) = chib3d(j,i,k,ircooh)*cfactor/W_RCOOH
          xr(:,:) = xrin(:,:)

          call chemmain

          xrout(:,:) = xr(:,:)

          do ic = 1 , totsp
            chemall(j,i,k,ic) = xrout(1,ic)
          end do
          ! Store photolysis rates
          do ic = 1 , 56
            jphoto(i,k,ic) = c_jval(1,ic)
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
          chemten(j,i,k,idms)   = &
            (xrout(1,ind_DMS)  - xrin(1,ind_DMS)) *pfact*W_DMS
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
        end do ! end i , k loop
      end do

      chemdiag(j,:,:,:) = chemdiag(j,:,:,:) + &
              chemten(j,:,:,:) * dble(dtchsolv) / (3600D0 * dble(chfrq)) 

    end subroutine chemistry
!
end module mod_che_chemistry

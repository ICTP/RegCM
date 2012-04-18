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


  real(dp) , parameter :: dtchsolv=900.E00
! 

  public :: gas_phase, dtchsolv

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

  real(dp) , parameter :: kb = 1.380658E-19

  contains

    subroutine chemistry(jj,chemin,chemox,taa,psaa,zena,idatein,tod)

      implicit none

      integer , intent(in) :: jj
      integer , intent(in) :: idatein
      real(dp) , intent(in) :: tod ! abt added for time of day
      real(dp) , dimension(2:iym2,1:kz,totsp) , intent(in) :: chemin
      real(dp) , dimension(2:iym2,1:kz,totsp) , intent(out) :: chemox
      real(dp) , dimension(2:iym2,1:kz) , intent(in) :: taa , psaa
      real(dp) , dimension(2:iym2) , intent(in) :: zena
      ! LOCAL VARIABLES
      real(dp) , dimension(2:iym2,1:kz) :: cfactor
      real(dp) , dimension(2:iym2,1:kz,1:56) :: jphoto
      !ah  variable for photolysis
      real(dp) , dimension(1:kz) :: taucab , taucbl
      real(dp) , dimension(1:iy,1:kz,1:jxp) :: taucld2
      real(dp) , dimension(1:iy,0:kz+1) :: psaa2 , taa2
      real(dp) , dimension(1:kz) :: hcbl , hcab
      real(dp) :: levav
      integer :: i , k , kbl , kab , ll,ic

      time = dtchsolv
      idate = idatein
      xhour = tod        !abt added for time of day
      c_numitr = 20

      ! initialize jphoto to zero
      jphoto(:,:,:) = d_zero

      ! Reorder from top-down to bottom-up
      do k = 1 , kz
        do i = 2 , iym2
          ll = (kz+1)-k
          taucld2(i,ll,jj) = taucld(i,k,jj)
          psaa2(i,ll) = psaa(i,k)
          taa2(i,ll) = taa(i,k)
        end do
      end do

      do k = kz , 1 , -1
        do i = 2 , iym2
          ll = (kz+1)-k
          altmid(1) = psaa2(i,k)
          temp(1) = taa2(i,k)
          zenith = zena(i) 
          cfactor(i,k) = psaa2(i,k)*d_10/(kb*taa2(i,k))
          dens(1)   = cfactor(i,k)
          taucab(ll) = d_zero
          hcab(ll) = d_zero
          taucbl(ll) = d_zero
          hcbl(ll) = d_zero
          if ( ll < kz ) then
            do kab = ll+1 , kz
              taucab(ll) = taucab(ll)+taucld2(i,kab,jj)
              levav = d_half*(psaa2(i,kab)+psaa2(i,kab+1))
              hcab(ll) = hcab(ll) + taucld2(i,kab,jj)*levav
            end do
          end if
          taucab(ll) = taucab(ll) +taucld2(i,ll,jj)*d_half
          levav = d_half*(psaa2(i,ll)+psaa2(i,ll+1))
          hcab(ll) = hcab(ll) + taucld2(i,ll,jj)*d_half*levav
          if ( dabs(taucab(ll)) < dlowval ) then
            hcab(ll) = d_zero
          else
            hcab(ll) = hcab(ll)/taucab(ll)
          end if
          deptha = taucab(ll)
          altabove = hcab(ll)
          if ( ll > 1 ) then
            do kbl = 1 , ll-1
              taucbl(ll) = taucbl(ll)+taucld2(i,kbl,jj)
              levav = d_half*(psaa2(i,kbl)+psaa2(i,kbl-1))
              hcbl(ll) = hcbl(ll) + taucld2(i,kbl,jj)*levav
            end do
          end if
          taucbl(ll) = taucbl(ll) + taucld2(i,ll,jj)*d_half
          levav = d_half*(psaa2(i,ll)+psaa2(i,ll-1))
          hcbl(ll) = hcbl(ll) + taucld2(i,ll,jj)*d_half*levav
          if (dabs(taucbl(ll)) < dlowval ) then
            hcbl(ll) = d_zero
          else
            hcbl(ll) = hcbl(ll)/taucbl(ll)
          end if
          depthb = taucbl(ll)
          altbelow = hcbl(ll)
          kmax = 1

          do ic = 1 , totsp
            xr(1,ic) = d_zero
          end do

          do ic = 1 , totsp
            xr(1,ic) = chemall(jj,i,k,ic) 
          end do

          xh2o           = chemin(i,k,ind_H2O)
          xr(1,ind_H2O)  = xh2o
          xr(1,ind_H2)   = 0.1E13
          xr(1,ind_O3)   = chemin(i,k,ind_O3) !0.094E+13
          xr(1,ind_NO2)  = chemin(i,k,ind_NO2) !0.300E+09
          xr(1,ind_NO)   = chemin(i,k,ind_NO) !0.300E+09
          xr(1,ind_CO)   = chemin(i,k,ind_CO) !0.100E+13
          xr(1,ind_H2O2) = chemin(i,k,ind_H2O2) !0.200E+11
          xr(1,ind_HNO3) = chemin(i,k,ind_HNO3) !0.200E+10
          xr(1,ind_N2O5) = chemin(i,k,ind_N2O5) !0.100E+08
          xr(1,ind_SO2)  = chemin(i,k,ind_SO2)  !0.200E+11
          xr(1,ind_SULF) = chemin(i,k,ind_SULF) !0.200E+11
          xr(1,ind_DMS)  = chemin(i,k,ind_DMS) !0.200E+10
          xr(1,ind_HCHO) = chemin(i,k,ind_HCHO) !0.200E+10
          xr(1,ind_ALD2) = chemin(i,k,ind_ALD2) !0.200E+10
          xr(1,ind_ISOP) = chemin(i,k,ind_ISOP) !0.500E+10
          xr(1,ind_C2H6) = chemin(i,k,ind_C2H6) !0.500E+10
          xr(1,ind_PAR)  = chemin(i,k,ind_PAR) !0.200E+08
          xr(1,ind_ACET) = chemin(i,k,ind_ACET) !0.200E+10
          xr(1,ind_MOH)  = chemin(i,k,ind_MOH) !0.200E+10
          xr(1,ind_PRPE) = chemin(i,k,ind_PRPE) !0.200E+09
          xr(1,ind_BUTE) = chemin(i,k,ind_BUTE) !0.200E+07
          xr(1,ind_TOLU) = chemin(i,k,ind_TOLU) !0.200E+07
          xr(1,ind_XYLE) = chemin(i,k,ind_XYLE) !0.000E+10
          xr(1,ind_CH4)  = chemin(i,k,ind_CH4)
          xr(1,ind_PAN)  = chemin(i,k,ind_PAN) !0.750E+10
          xr(1,ind_ETHE) = chemin(i,k,ind_ETHE) !10.200E+09

 
!fab deb
!!$          xr(1,ind_O3)   =  0.094E+13
!!$          xr(1,ind_NO2)  = 0.300E+09
!!$          xr(1,ind_NO)   = 0.300E+09
!!$          xr(1,ind_CO)   = 0.100E+13
!!$          xr(1,ind_H2O2) = 0.200E+11
!!$          xr(1,ind_HNO3) = 0.200E+10
!!$          xr(1,ind_N2O5) = 0.100E+08
!!$          xr(1,ind_SO2)  = 0.200E+11
!!$          xr(1,ind_SULF) = 0.200E+11
!!$          xr(1,ind_DMS)  = 0.200E+10
!!$          xr(1,ind_HCHO) = 0.200E+10
!!$          xr(1,ind_ALD2) = 0.200E+10
!!$          xr(1,ind_ISOP) = 0.500E+10
!!$          xr(1,ind_C2H6) = 0.500E+10
!!$          xr(1,ind_PAR)  = 0.200E+08
!!$          xr(1,ind_ACET) = 0.200E+10
!!$          xr(1,ind_MOH)  = 0.200E+10
!!$          xr(1,ind_PRPE) = 0.200E+09
!!$          xr(1,ind_BUTE) = 0.200E+07
!!$          xr(1,ind_TOLU) = 0.200E+07
!!$          xr(1,ind_XYLE) = 0.000E+10
!!$       
!!$          xr(1,ind_PAN)  =  0.750E+10
!!$          xr(1,ind_ETHE) = 10.200E+09
!!$
!!$
!!$         xr(1,:) = d_zero

          call chemmain

 
          

          do ic = 1 , totsp
            chemall(jj,i,k,ic) = xr(1,ic)
          end do
          ! Store photolysis rates
          do ic = 1 , 56
            jphoto(i,k,ic) = c_jval(1,ic)
          end do

          chemox(i,k,ind_O3)   = xr(1,ind_O3)
          chemox(i,k,ind_NO2)  = xr(1,ind_NO2)
          chemox(i,k,ind_NO)   = xr(1,ind_NO)
          chemox(i,k,ind_CO)   = xr(1,ind_CO)
          chemox(i,k,ind_H2O2) = xr(1,ind_H2O2)
          chemox(i,k,ind_HNO3) = xr(1,ind_HNO3)
          chemox(i,k,ind_N2O5) = xr(1,ind_N2O5)
          chemox(i,k,ind_SO2)  = xr(1,ind_SO2)
          chemox(i,k,ind_SULF) = xr(1,ind_SULF)
          chemox(i,k,ind_DMS)  = xr(1,ind_DMS)
          chemox(i,k,ind_HCHO) = xr(1,ind_HCHO)
          chemox(i,k,ind_ALD2) = xr(1,ind_ALD2)
          chemox(i,k,ind_ISOP) = xr(1,ind_ISOP)
          chemox(i,k,ind_C2H6) = xr(1,ind_C2H6)
          chemox(i,k,ind_PAR)  = xr(1,ind_PAR)
          chemox(i,k,ind_TOLU) = xr(1,ind_TOLU)
          chemox(i,k,ind_XYLE) = xr(1,ind_XYLE)
          chemox(i,k,ind_ETHE) = xr(1,ind_ETHE)
          chemox(i,k,ind_PAN)  = xr(1,ind_PAN)
          chemox(i,k,ind_CH4)  = xr(1,ind_CH4)
          chemox(i,k,ind_PRPE) = xr(1,ind_PRPE)
          chemox(i,k,ind_BUTE) = xr(1,ind_BUTE)
          chemox(i,k,ind_MOH)  = xr(1,ind_MOH)
          chemox(i,k,ind_ACET) = xr(1,ind_ACET)
        end do
      end do
    end subroutine chemistry
!
!-------------------------------------------------------------------
!
    subroutine gas_phase(j,secofday,lyear,lmonth,lday)
      implicit none
      integer , intent(in) :: j
      integer , intent(in) :: lyear , lmonth , lday
      real(dp) , intent(in) :: secofday

      real(dp) , dimension(ici1:ici2,1:kz) :: taa , psaa
      real(dp) , dimension(jci1:jci2,ici1:ici2,1:kz,totsp) :: chemin , chemox
      real(dp) :: cfactor , pfact
      real(dp) , dimension(ici1:ici2) :: zena
      real(dp) :: tod
      integer :: idatein
      integer :: i , k

      chemin(j,:,:,:) = d_zero
      chemox(j,:,:,:) = d_zero
     
      do k = 1 , kz
        do i = 2 , iym2
          taa(i,k) = ctb3d(j,i,k)
          psaa(i,k)= (cpsb(j,i)*hlev(k)+chptop)
        end do
      end do 

      zena(:) = dacos(czen(j,:)*degrad)

      do k = 1 , kz
        do i = 2 , iym2

!work with chib3d arrays which are in kg.kg-1 
!convert into molec.cm-3 

          cfactor =  crhob3d(j,i,k) * 1.D-03 * navgdr
          chemin(j,i,k,ind_H2O)  = cqvb3d(j,i,k)*cfactor / 18.D00

          chemin(j,i,k,ind_O3)   = chib3d(j,i,k,io3)*cfactor/W_O3
          chemin(j,i,k,ind_NO2)  = chib3d(j,i,k,ino2)*cfactor /W_NO2
          chemin(j,i,k,ind_NO)   = chib3d(j,i,k,ino)*cfactor/W_NO
          chemin(j,i,k,ind_CO)   = chib3d(j,i,k,ico)*cfactor/W_CO
          chemin(j,i,k,ind_H2O2) = chib3d(j,i,k,ih2o2)*cfactor/W_H2O2
          chemin(j,i,k,ind_HNO3) = chib3d(j,i,k,ihno3)*cfactor/W_HNO3
          chemin(j,i,k,ind_N2O5) = chib3d(j,i,k,in2o5)*cfactor/W_N2O5
          chemin(j,i,k,ind_SO2)  = chib3d(j,i,k,iso2)*cfactor/W_SO2
          chemin(j,i,k,ind_SULF) = chib3d(j,i,k,iso4)*cfactor/W_SULF
          chemin(j,i,k,ind_DMS)  = chib3d(j,i,k,idms)*cfactor/W_DMS
          chemin(j,i,k,ind_HCHO) = chib3d(j,i,k,ihcho)*cfactor/W_HCHO
          chemin(j,i,k,ind_ALD2) = chib3d(j,i,k,iald2)*cfactor/W_ALD2
          chemin(j,i,k,ind_ISOP) = chib3d(j,i,k,iisop)*cfactor/W_ISOP
          chemin(j,i,k,ind_C2H6) = chib3d(j,i,k,ic2h6)*cfactor/W_C2H6
          chemin(j,i,k,ind_PAR)  = chib3d(j,i,k,ipar)*cfactor/W_C3H8
          chemin(j,i,k,ind_ETHE) = chib3d(j,i,k,iethe)*cfactor/W_ETHENE
          chemin(j,i,k,ind_PRPE) = chib3d(j,i,k,iolt)*cfactor/W_OLT
          chemin(j,i,k,ind_BUTE) = chib3d(j,i,k,ioli)*cfactor/W_OLI
          chemin(j,i,k,ind_TOLU) = chib3d(j,i,k,itolue)*cfactor/W_TOLU
          chemin(j,i,k,ind_XYLE) = chib3d(j,i,k,ixyl)*cfactor/W_XYLE
          chemin(j,i,k,ind_PAN)  = chib3d(j,i,k,ipan)*cfactor/W_PAN
          chemin(j,i,k,ind_CH4)  = chib3d(j,i,k,ich4)*cfactor/W_CH4
          chemin(j,i,k,ind_MOH)  = chib3d(j,i,k,imoh)*cfactor/W_MOH
          chemin(j,i,k,ind_ACET) = chib3d(j,i,k,iacet)*cfactor/W_ACET
        end do
      end do

      tod = secofday/3600.0D0
     
      idatein = (lyear-1900)*10000+lmonth*100+lday

      call chemistry(j,chemin(j,:,:,:),chemox(j,:,:,:), &
                     taa,psaa,zena,idatein,tod)

      ! FAB :  Now save the chemistry tendency 
      ! be carefull should be multiplied by surface pressure
      ! ( consistency with chib unit)

      do k = 1 , kz
        do i = ici1 , ici2
          ! convection factor to get the tendency from 
          ! mole.cm-3.s-1  to kg.kg-1.s-1.ps (consistency with chiten unit)
          cfactor =  crhob3d(j,i,k) * 1.D-03 * navgdr
          pfact = cpsb(j,i)/cfactor/dtchsolv
          chemten(j,i,k,io3)    = &
             (chemox(j,i,k,ind_O3)  - chemin(j,i,k,ind_O3))  *pfact*W_O3
          chemten(j,i,k,ino2)   = &
             (chemox(j,i,k,ind_NO2) - chemin(j,i,k,ind_NO2)) *pfact*W_NO2
          chemten(j,i,k,ino)    = &
            (chemox(j,i,k,ind_NO)   - chemin(j,i,k,ind_O3))  *pfact*W_NO
          chemten(j,i,k,ico)    = &
            (chemox(j,i,k,ind_CO)   - chemin(j,i,k,ind_CO))  *pfact*W_CO
          chemten(j,i,k,ih2o2)  = &
            (chemox(j,i,k,ind_H2O2) - chemin(j,i,k,ind_H2O2))*pfact*W_H2O2
          chemten(j,i,k,ihno3)  = &
            (chemox(j,i,k,ind_HNO3) - chemin(j,i,k,ind_HNO3))*pfact*W_HNO3
          chemten(j,i,k,in2o5)  = &
            (chemox(j,i,k,ind_N2O5) - chemin(j,i,k,ind_N2O5))*pfact*W_N2O5
          chemten(j,i,k,iso2)   = &
            (chemox(j,i,k,ind_SO2)  - chemin(j,i,k,ind_SO2)) *pfact*W_SO2
          chemten(j,i,k,iso4)   = &
            (chemox(j,i,k,ind_SULF) - chemin(j,i,k,ind_SULF))*pfact*W_SULF
          chemten(j,i,k,idms)   = &
            (chemox(j,i,k,ind_DMS)  - chemin(j,i,k,ind_DMS)) *pfact*W_DMS
          chemten(j,i,k,ihcho)  = &
            (chemox(j,i,k,ind_HCHO) - chemin(j,i,k,ind_HCHO))*pfact*W_HCHO
          chemten(j,i,k,iald2)  = &
            (chemox(j,i,k,ind_ALD2) - chemin(j,i,k,ind_ALD2))*pfact*W_ALD2
          chemten(j,i,k,iisop)  = &
            (chemox(j,i,k,ind_ISOP) - chemin(j,i,k,ind_ISOP))*pfact*W_ISOP
          chemten(j,i,k,ic2h6)  = &
            (chemox(j,i,k,ind_C2H6) - chemin(j,i,k,ind_C2H6))*pfact*W_C2H6
          chemten(j,i,k,ipar)   = &
            (chemox(j,i,k,ind_PAR)  - chemin(j,i,k,ind_PAR)) *pfact*W_C3H8
          chemten(j,i,k,itolue) = &
            (chemox(j,i,k,ind_TOLU) - chemin(j,i,k,ind_TOLU))*pfact*W_TOLU
          chemten(j,i,k,ixyl)   = &
            (chemox(j,i,k,ind_XYLE) - chemin(j,i,k,ind_XYLE))*pfact*W_XYLE
          chemten(j,i,k,iethe)  = &
            (chemox(j,i,k,ind_ETHE) - chemin(j,i,k,ind_ETHE))*pfact*W_ETHENE
          chemten(j,i,k,ipan)   = &
            (chemox(j,i,k,ind_PAN)  - chemin(j,i,k,ind_PAN)) *pfact*W_PAN
          chemten(j,i,k,ich4)   = &
            (chemox(j,i,k,ind_CH4)  - chemin(j,i,k,ind_CH4)) *pfact*W_CH4
          chemten(j,i,k,iolt)   = &
            (chemox(j,i,k,ind_PRPE) - chemin(j,i,k,ind_PRPE))*pfact*W_OLT
          chemten(j,i,k,ioli)   = &
            (chemox(j,i,k,ind_BUTE) - chemin(j,i,k,ind_BUTE))*pfact*W_OLI
          chemten(j,i,k,imoh)   = &
            (chemox(j,i,k,ind_MOH)  - chemin(j,i,k,ind_MOH)) *pfact*W_MOH
          chemten(j,i,k,iacet)  = &
            (chemox(j,i,k,ind_ACET) - chemin(j,i,k,ind_ACET))*pfact*W_ACET
       end do
    end do

end subroutine gas_phase

end module mod_che_chemistry

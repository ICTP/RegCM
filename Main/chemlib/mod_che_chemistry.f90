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

  private

  public :: gas_phase

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

    subroutine chemistry(jj,chemin,chemox,taa,psaa,zena, &
                         ktau,idatein,tod)

      implicit none

      integer , intent(in) :: jj
      integer(8) , intent(in) :: ktau
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
      integer :: i , k , kbl , kab , ll

      time = 900.0D0
      idate = idatein
      ktaubx = ktau
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
            xr(1,ic) = chemall(i,k,jj,ic) 
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

          call chemmain

          do ic = 1 , totsp
            chemall(i,k,jj,ic) = xr(1,ic)
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
!c -------------------------------------------------------------------
    subroutine gas_phase(j,ktau,secofday,lyear,lmonth,lday)

      implicit none

      integer , intent(in) :: j
      integer , intent(in) :: lyear , lmonth , lday
      integer(8) :: ktau
      real(dp) , intent(in) :: secofday

      real(dp) , dimension(2:iym2,1:kz) :: taa , psaa , rhoa
      real(dp) , dimension(2:iym2,1:kz,1:jxp,totsp) :: chemin , chemox
      real(dp) , dimension(1:iym2,1:kz,56) :: jphoto
      real(dp) , dimension(iy,jxp) :: psdot
      real(dp) :: airmw , cfactor , ccfactor , pfact , ro3 , rh2o2 , kb
      real(dp) , dimension(2:iym2) :: zena
      real(dp) , dimension(iym1,ntr) :: vdep
      real(dp) :: fact1 , fact2 , srctemp , tod
      integer :: idatein , tstart
      integer :: itr , i , k , nn

      chemin(:,:,j,:) = d_zero
      chemox(:,:,j,:) = d_zero

      do k = 1 , kz
        do i = 2 , iym2
          taa(i,k) = ctb3d(j,i,k)
          psaa(i,k)= (cpsb(i,j)*hlev(k)+chptop)
        end do
      end do 

      do k = 1 , kz
        do i = 2 , iym2
          cfactor = psaa(i,k)*d_10/(kb*taa(i,k))
          cfactor = cfactor/cpsb(i,j)
          zena(i) = dacos(czen(j,i)*degrad)
          chemin(i,k,j,ind_H2O)  = cqvb3d(j,i,k)*cfactor
          chemin(i,k,j,ind_O3)   = chib(i,k,j,io3)*cfactor*amd/W_O3
          chemin(i,k,j,ind_NO2)  = chib(i,k,j,ino2)*cfactor*amd/W_NO2
          chemin(i,k,j,ind_NO)   = chib(i,k,j,ino)*cfactor*amd/W_NO
          chemin(i,k,j,ind_CO)   = chib(i,k,j,ico)*cfactor*amd/W_CO
          chemin(i,k,j,ind_H2O2) = chib(i,k,j,ih2o2)*cfactor*amd/W_H2O2
          chemin(i,k,j,ind_HNO3) = chib(i,k,j,ihno3)*cfactor*amd/W_HNO3
          chemin(i,k,j,ind_N2O5) = chib(i,k,j,in2o5)*cfactor*amd/W_N2O5
          chemin(i,k,j,ind_SO2)  = chib(i,k,j,iso2)*cfactor*amd/W_SO2
          chemin(i,k,j,ind_SULF) = chib(i,k,j,iso4)*cfactor*amd/W_SULF
          chemin(i,k,j,ind_DMS)  = chib(i,k,j,idms)*cfactor*amd/W_DMS
          chemin(i,k,j,ind_HCHO) = chib(i,k,j,ihcho)*cfactor*amd/W_HCHO
          chemin(i,k,j,ind_ALD2) = chib(i,k,j,iald2)*cfactor*amd/W_ALD2
          chemin(i,k,j,ind_ISOP) = chib(i,k,j,iisop)*cfactor*amd/W_ISOP
          chemin(i,k,j,ind_C2H6) = chib(i,k,j,ic2h6)*cfactor*amd/W_C2H6
          chemin(i,k,j,ind_PAR)  = chib(i,k,j,ipar)*cfactor*amd/W_C3H8
          chemin(i,k,j,ind_ETHE) = chib(i,k,j,iethe)*cfactor*amd/W_ETHENE
          chemin(i,k,j,ind_PRPE) = chib(i,k,j,iolt)*cfactor*amd/W_OLT
          chemin(i,k,j,ind_BUTE) = chib(i,k,j,ioli)*cfactor*amd/W_OLI
          chemin(i,k,j,ind_TOLU) = chib(i,k,j,itolue)*cfactor*amd/W_TOLU
          chemin(i,k,j,ind_XYLE) = chib(i,k,j,ixyl)*cfactor*amd/W_XYLE
          chemin(i,k,j,ind_PAN)  = chib(i,k,j,ipan)*cfactor*amd/W_PAN
          chemin(i,k,j,ind_CH4)  = chib(i,k,j,ich4)*cfactor*amd/W_CH4
          chemin(i,k,j,ind_MOH)  = chib(i,k,j,imoh)*cfactor*amd/W_MOH
          chemin(i,k,j,ind_ACET) = chib(i,k,j,iacet)*cfactor*amd/W_ACET
        end do
      end do

      tod = secofday/3600.0D0
      idatein = (lyear-1900)*10000+lmonth*100+lday

      call chemistry(j,chemin(:,:,j,:),chemox(:,:,j,:),      &
                     taa,psaa,zena,ktau,idatein,tod)

      do k = 1 , kz
        do i = 2 , iym2
          cfactor = psaa(i,k)*d_10/(kb*taa(i,k))
          pfact = cpsb(i,j)/(cfactor*amd)
          ! pfact to convert chemox to kg/kg
          chia(i,k,j,io3)    = chemox(i,k,j,ind_O3)*pfact*W_O3
          chia(i,k,j,ino2)   = chemox(i,k,j,ind_NO2)*pfact*W_NO2
          chia(i,k,j,ino)    = chemox(i,k,j,ind_NO)*pfact*W_NO
          chia(i,k,j,ico)    = chemox(i,k,j,ind_CO)*pfact*W_CO
          chia(i,k,j,ih2o2)  = chemox(i,k,j,ind_H2O2)*pfact*W_H2O2
          chia(i,k,j,ihno3)  = chemox(i,k,j,ind_HNO3)*pfact*W_HNO3
          chia(i,k,j,in2o5)  = chemox(i,k,j,ind_N2O5)*pfact*W_N2O5
          chia(i,k,j,iso2)   = chemox(i,k,j,ind_SO2)*pfact*W_SO2
          chia(i,k,j,iso4)   = chemox(i,k,j,ind_SULF)*pfact*W_SULF
          chia(i,k,j,idms)   = chemox(i,k,j,ind_DMS)*pfact*W_DMS
          chia(i,k,j,ihcho)  = chemox(i,k,j,ind_HCHO)*pfact*W_HCHO
          chia(i,k,j,iald2)  = chemox(i,k,j,ind_ALD2)*pfact*W_ALD2
          chia(i,k,j,iisop)  = chemox(i,k,j,ind_ISOP)*pfact*W_ISOP
          chia(i,k,j,ic2h6)  = chemox(i,k,j,ind_C2H6)*pfact*W_C2H6
          chia(i,k,j,ipar)   = chemox(i,k,j,ind_PAR)*pfact*W_C3H8
          chia(i,k,j,itolue) = chemox(i,k,j,ind_TOLU)*pfact*W_TOLU
          chia(i,k,j,ixyl)   = chemox(i,k,j,ind_XYLE)*pfact*W_XYLE
          chia(i,k,j,iethe)  = chemox(i,k,j,ind_ETHE)*pfact*W_ETHENE
          chia(i,k,j,ipan)   = chemox(i,k,j,ind_PAN)*pfact*W_PAN
          chia(i,k,j,ich4)   = chemox(i,k,j,ind_CH4)*pfact*W_CH4
          chia(i,k,j,iolt)   = chemox(i,k,j,ind_PRPE)*pfact*W_OLT
          chia(i,k,j,ioli)   = chemox(i,k,j,ind_BUTE)*pfact*W_OLI
          chia(i,k,j,imoh)   = chemox(i,k,j,ind_MOH)*pfact*W_MOH
          chia(i,k,j,iacet)  = chemox(i,k,j,ind_ACET)*pfact*W_ACET
        end do
      end do
    end subroutine gas_phase

end module mod_che_chemistry

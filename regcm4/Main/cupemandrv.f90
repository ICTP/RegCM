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
 
      subroutine cupemandrv(j)
 
! **********************************************
! **** Driver for Emanuel Convection Scheme ****
! **********************************************
 
      use mod_regcm_param
      use mod_param1 , only : dt , nbatst
      use mod_param3 , only : ptop , sigma
      use mod_main
      use mod_pmoist
      use mod_cvaria
      use mod_slice
      use mod_rad
      use mod_bats , only : pptc
      use mod_date , only : jyear , jyear0 , ktau
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ntra = 0 , nl = kz - 1
!
! Dummy arguments
!
      integer :: j
      intent (in) j
!
! Local variables
!
      real(8) :: akclth , aprdiv , cbmf , dtime , fppt , qprime ,       &
               & tprime , uconv , wd
      real(8) , dimension(kz) :: fq , ft , fu , fv , pcup , qcup ,      &
                               & qscup , tcup , ucup , vcup
      real(8) , dimension(kz,1) :: ftra , tra
      integer :: i , iconj , iflag , k , kbase , kclth , kk , ktop
      real(8) , dimension(kzp1) :: phcup
!
      dtime = dt
      uconv = 0.5*dt
      aprdiv = 1.0/float(nbatst)
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) aprdiv = 1.
      iconj = 0
      do i = 2 , iym2
        do k = 1 , kz
          kk = kzp1 - k
          cldlwc(i,k) = 0.          ! Zero out cloud water content
          cldfra(i,k) = 0.          ! Zero out cloud fraction coverage
          tcup(k) = tb3d(i,kk,j)                          ! [k]
          qcup(k) = qvb3d(i,kk,j)/(1.+qvb3d(i,kk,j))      ! [kg/kg]
          qscup(k) = qsb3d(i,kk,j)/(1.+qsb3d(i,kk,j))     ! [kg/kg]
          ucup(k) = ubx3d(i,kk,j)                         ! [m/s]
          vcup(k) = vbx3d(i,kk,j)                         ! [m/s]
          pcup(k) = pb3d(i,kk,j)*10.                      ! [hPa]
        end do
        do k = 1 , kzp1
          kk = kzp1 - k + 1
          phcup(k) = (sigma(kk)*psb(i,j)+ptop)*10.        ! [hPa]
        end do
        cbmf = cbmf2d(i,j)                                ! [(kg/m**2)/s]
 
        call cupeman(tcup,qcup,qscup,ucup,vcup,tra,pcup,phcup,kz,kzp1,  &
                   & nl,ntra,dtime,iflag,ft,fq,fu,fv,ftra,fppt,wd,      &
                   & tprime,qprime,cbmf,kbase,ktop)
 
        cbmf2d(i,j) = cbmf
 
!       iflag=0: No moist convection; atmosphere stable or surface
!       temperature < 250K or surface humidity is negative.
!       iflag=1: Moist convection occurs.
!       iflag=2: No moist convection: lifted condensation level above
!       200 mb. iflag=3: No moist convection: cloud base higher than
!       level kzm2. iflag=4: Moist convection occurs, but CFL condition
!       on the subsidence warming is violated. (Does not terminate
!       scheme.)
        if ( iflag.eq.1 .or. iflag.eq.4 ) then
                                              ! If moist convection
 
!         if (iflag.eq.4) then                ! If CFL violation
!         print*,'EMAN CFL VIOLATION: ',i,j,cbmf
!         end if
 
!         **** Tendencies
          do k = 1 , kz
            kk = kzp1 - k
            tten(i,kk,j) = ft(k)*psb(i,j) + tten(i,kk,j)
            qvten(i,kk,j) = fq(k)/(1.-fq(k))*psb(i,j) + qvten(i,kk,j)
!           There is a bit of an inconsistency here...  The wind
!           tendencies from convection are on cross points, but the
!           model wants them on dot points.
            uten(i,kk,j) = fu(k)*psb(i,j) + uten(i,kk,j)
            vten(i,kk,j) = fv(k)*psb(i,j) + vten(i,kk,j)
          end do
 
!         **** Cloud fraction and cloud water
          kclth = ktop - kbase + 1
          akclth = 1./float(kclth)
          do k = kbase , ktop
            kk = kzp1 - k
            cldlwc(i,kk) = cllwcv
            cldfra(i,kk) = 1. - (1.-clfrcv)**akclth
          end do
 
!         **** Precipitation
          if ( fppt.gt.0. ) then
            rainc(i,j) = rainc(i,j) + fppt*uconv ! mm
            pptc(i,j) = pptc(i,j) + fppt*aprdiv  ! mm/s
            iconj = iconj + 1
          end if
 
        end if
 
      end do
 
      icon(j) = iconj
 
      end subroutine cupemandrv
